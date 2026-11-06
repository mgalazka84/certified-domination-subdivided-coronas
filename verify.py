"""Reproduce the checks for certified domination in subdivided coronas.

The direct subset search and the integer models use the full constructed graph
and the defining constraints. They do not assume the paper's reduction.
"""
from __future__ import annotations
import csv
import itertools as it
import json
import platform
import random
import time
from collections import deque
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.optimize._highspy import _core as highs_core
from scipy.sparse import lil_matrix


def unlabelled_graphs(n):
    pairs = list(it.combinations(range(n), 2))
    index = {edge: i for i, edge in enumerate(pairs)}
    permutations = [[index[tuple(sorted((p[u], p[v])))] for u, v in pairs]
                    for p in it.permutations(range(n))]
    seen = set()
    for mask in range(1 << len(pairs)):
        if mask in seen:
            continue
        bits = [i for i in range(len(pairs)) if mask >> i & 1]
        seen.update(sum(1 << p[i] for i in bits) for p in permutations)
        yield [pairs[i] for i in bits]


def construct(n, edges, r):
    adjacency = [set() for _ in range(n)]
    supports, leaves, middle = [], [], []
    def vertex():
        adjacency.append(set())
        return len(adjacency) - 1
    def link(u, v):
        adjacency[u].add(v)
        adjacency[v].add(u)
    for v in range(n):
        sv, lv = [], []
        for _ in range(r[v]):
            s, leaf = vertex(), vertex()
            link(v, s)
            link(s, leaf)
            sv.append(s)
            lv.append(leaf)
        supports.append(sv)
        leaves.append(lv)
    for u, v in edges:
        w = vertex()
        link(u, w)
        link(v, w)
        middle.append(w)
    return adjacency, supports, leaves, middle


def is_certified(adj, D):
    return (all(v in D or adj[v] & D for v in range(len(adj)))
            and all(len(adj[v] - D) != 1 for v in D))


def exact_model(adj, certified=True):
    """Binary variables x encode D, y select the >=2 alternative."""
    N = len(adj)
    width, height = (2*N, 4*N) if certified else (N, N)
    A = lil_matrix((height, width), dtype=float)
    lo, hi = np.full(height, -np.inf), np.zeros(height)
    for v, neighbours in enumerate(adj):
        row = 4*v if certified else v
        A[row, v] = 1
        for u in neighbours:
            A[row, u] = 1
        lo[row], hi[row] = 1, np.inf
        if not certified:
            continue
        for u in neighbours:
            A[row+1, u] = 1
        A[row+1, N+v], hi[row+1] = 2, len(neighbours)
        A[row+2, v] = len(neighbours)
        for u in neighbours:
            A[row+2, u] = -1
        A[row+2, N+v] = -len(neighbours)
        A[row+3, N+v], A[row+3, v] = 1, -1
    cost = np.r_[np.ones(N), np.zeros(width-N)]
    result = milp(cost, integrality=np.ones(width), bounds=Bounds(0, 1),
                  constraints=LinearConstraint(A.tocsc(), lo, hi),
                  options={'mip_rel_gap': 0.0})
    if not result.success:
        raise RuntimeError(result.message)
    D = {i for i in range(N) if result.x[i] > 0.5}
    assert all(v in D or adj[v] & D for v in range(N))
    if certified:
        assert is_certified(adj, D)
    assert abs(result.fun-len(D)) < 1e-6
    return len(D)


def direct_all_optima(adj):
    """Examine full vertex subsets in increasing size, with no forced vertices."""
    N = len(adj)
    ordinary = None
    masks = [sum(1 << u for u in neighbours) for neighbours in adj]
    full = (1 << N) - 1
    for k in range(N+1):
        solutions = []
        for vertices in it.combinations(range(N), k):
            D = sum(1 << v for v in vertices)
            outside = full ^ D
            unseen = outside
            dominated = True
            while unseen:
                bit = unseen & -unseen
                v = bit.bit_length()-1
                if not masks[v] & D:
                    dominated = False
                    break
                unseen -= bit
            if not dominated:
                continue
            if ordinary is None:
                ordinary = k
            if all((masks[v] & outside).bit_count() != 1 for v in vertices):
                solutions.append(D)
        if solutions:
            return ordinary, k, set(solutions)
    raise AssertionError('The full vertex set is always feasible.')


def base_formula(n, edges, r):
    b = [x+1 for x in r]
    best, minimizers = float('inf'), []
    cover = n
    for X in range(1 << n):
        residual = sum(not ((X >> u & 1) or (X >> v & 1)) for u, v in edges)
        value = sum(b[v] for v in range(n) if X >> v & 1) + residual
        if value < best:
            best, minimizers = value, [X]
        elif value == best:
            minimizers.append(X)
        if residual == 0:
            cover = min(cover, X.bit_count())
    return sum(r)+best, minimizers, sum(r)+cover


def canonical(X, supports, leaves, middle, edges):
    D = {s for group in supports for s in group}
    for v, group in enumerate(leaves):
        if X >> v & 1:
            D.add(v)
            D.update(group)
    D.update(w for w, (u, v) in zip(middle, edges)
             if not ((X >> u & 1) or (X >> v & 1)))
    return D


def bipartition(n, edges):
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    color = [-1]*n
    for root in range(n):
        if color[root] != -1:
            continue
        color[root] = 0
        queue = deque([root])
        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if color[v] == -1:
                    color[v] = 1-color[u]
                    queue.append(v)
                elif color[v] == color[u]:
                    return None
    return color


def bounded_edges(n, edges, b, color):
    """Integral augmenting paths; return matching and minimum-cut witnesses."""
    s, t = n, n+1
    neighbours = [set() for _ in range(n+2)]
    residual = {}
    def arc(u, v, c):
        neighbours[u].add(v)
        neighbours[v].add(u)
        residual[u, v], residual[v, u] = c, 0
    for v in range(n):
        if color[v] == 0:
            arc(s, v, b[v])
        else:
            arc(v, t, b[v])
    oriented = [(u, v) if color[u] == 0 else (v, u) for u, v in edges]
    for u, v in oriented:
        arc(u, v, 1)
    value = 0
    while True:
        parent = {s: None}
        queue = deque([s])
        while queue and t not in parent:
            u = queue.popleft()
            for v in sorted(neighbours[u]):
                if v not in parent and residual[u, v] > 0:
                    parent[v] = u
                    queue.append(v)
        if t not in parent:
            break
        increment, v = float('inf'), t
        while v != s:
            u = parent[v]
            increment = min(increment, residual[u, v])
            v = u
        v = t
        while v != s:
            u = parent[v]
            residual[u, v] -= increment
            residual[v, u] += increment
            v = u
        value += increment
    reachable = set(parent)
    X = sum(1 << v for v in range(n)
            if (color[v] == 0 and v not in reachable)
            or (color[v] == 1 and v in reachable))
    F = [e for e, (u, v) in zip(edges, oriented) if residual[u, v] == 0]
    degree = [0]*n
    for u, v in F:
        degree[u] += 1
        degree[v] += 1
    assert len(F) == value and all(degree[v] <= b[v] for v in range(n))
    return value, X, F


def verify_case(stage, n, edges, r, direct=False):
    adj, supports, leaves, middle = construct(n, edges, r)
    predicted, optimal_X, ordinary = base_formula(n, edges, r)
    actual = exact_model(adj)
    actual_ordinary = exact_model(adj, certified=False)
    assert actual == predicted, (n, edges, r, actual, predicted)
    assert ordinary == actual_ordinary
    predicted_sets = set()
    for X in optimal_X:
        D = canonical(X, supports, leaves, middle, edges)
        assert len(D) == actual and is_certified(adj, D)
        predicted_sets.add(sum(1 << v for v in D))
    assert len(predicted_sets) == len(optimal_X)
    if direct:
        g, c, all_sets = direct_all_optima(adj)
        assert (g, c) == (ordinary, actual)
        assert all_sets == predicted_sets, (n, edges, r, 'all-optima mismatch')
    color = bipartition(n, edges)
    flow_value, certificate = None, None
    if color is not None:
        flow_value, X, certificate = bounded_edges(n, edges, [z+1 for z in r], color)
        assert sum(r)+flow_value == actual and X in optimal_X
    deg = [0]*n
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1
    if all(deg[v] <= r[v]+1 for v in range(n)):
        assert actual == sum(r)+len(edges)
        tight = sum(1 << v for v in range(n) if deg[v] == r[v]+1)
        independents = {X for X in range(1 << n)
                        if X & ~tight == 0 and
                        all(not (X >> u & 1 and X >> v & 1) for u, v in edges)}
        assert independents == set(optimal_X)
    else:
        assert actual < sum(r)+len(edges)
    return {'stage': stage, 'n': n, 'm': len(edges), 'edges': edges, 'r': list(r),
            'constructed_order': len(adj), 'gamma': actual_ordinary,
            'gamma_cer': actual, 'optimal_original_vertex_masks': optimal_X,
            'bipartite': color is not None, 'matching_value': flow_value,
            'matching_edges': certificate, 'all_optima_directly_checked': direct}


def main():
    started = time.time()
    results, stage_counts, graph_counts = [], {}, {}
    for n in range(1, 7):
        graphs = list(unlabelled_graphs(n))
        graph_counts[n] = len(graphs)
        for edges in graphs:
            results.append(verify_case('all_graphs_one_leaf', n, edges, [1]*n, direct=n<=4))
        print(f'all graphs: n={n}, representatives={len(graphs)}', flush=True)
    assert list(graph_counts.values()) == [1, 2, 4, 11, 34, 156]
    for n in range(1, 5):
        for edges in unlabelled_graphs(n):
            for r in it.product((1, 2, 3), repeat=n):
                direct = n <= 3 and max(r) <= 2
                results.append(verify_case('nonuniform_leaf_counts', n, edges, r, direct=direct))
        print(f'nonuniform leaves: n={n} finished', flush=True)
    rng = random.Random(20261106)
    for i in range(60):
        n = rng.randint(7, 12)
        if i < 30:
            split = rng.randint(1, n-1)
            pairs = [(u,v) for u in range(split) for v in range(split,n)]
        else:
            pairs = list(it.combinations(range(n), 2))
        p = rng.uniform(0.1, 0.85)
        edges = [e for e in pairs if rng.random() < p]
        r = [rng.randint(1,4) for _ in range(n)]
        results.append(verify_case('seeded_larger_cases', n, edges, r))
    fib = [0,1]
    for _ in range(20):
        fib.append(fib[-1]+fib[-2])
    enumeration_checks = []
    for n in range(2, 14):
        edges = [(v,v+1) for v in range(n-1)]
        value, xs, _ = base_formula(n, edges, [1]*n)
        assert value == 2*n-1 and len(xs) == fib[n]
        enumeration_checks.append({'family':'path','n':n,'value':value,'count':len(xs)})
    for n in range(3, 14):
        for family, edges, value_expected, count_expected in [
            ('cycle',[(v,(v+1)%n) for v in range(n)],2*n,fib[n-1]+fib[n+1]),
            ('complete',list(it.combinations(range(n),2)),3*n-3,
             n*(n-1)//2+n*(n-1)*(n-2)//6)]:
            value, xs, _ = base_formula(n, edges, [1]*n)
            assert value == value_expected and len(xs) == count_expected
            enumeration_checks.append({'family':family,'n':n,'value':value,'count':len(xs)})
    for row in results:
        stage_counts[row['stage']] = stage_counts.get(row['stage'],0)+1
    summary = {'status':'all passed','seed':20261106,'graph_counts':graph_counts,
               'case_counts':stage_counts,'total_cases':len(results),
               'milp_solves':2*len(results),
               'bipartite_cases':sum(x['bipartite'] for x in results),
               'direct_all_optima_cases':sum(x['all_optima_directly_checked'] for x in results),
               'enumeration_formula_checks':len(enumeration_checks),
               'versions':{'python':platform.python_version(),'numpy':np.__version__,
                           'scipy':scipy.__version__,
                           'highs':'.'.join(str(getattr(highs_core, 'HIGHS_VERSION_'+p))
                                            for p in ('MAJOR','MINOR','PATCH'))},
               'elapsed_seconds':round(time.time()-started,3)}
    destination = Path(__file__).resolve().parent/'results'
    destination.mkdir(exist_ok=True)
    (destination/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (destination/'cases.json').write_text(json.dumps(results,indent=2)+'\n')
    (destination/'enumeration.json').write_text(json.dumps(enumeration_checks,indent=2)+'\n')
    with (destination/'cases.csv').open('w',newline='') as stream:
        fields=['stage','n','m','r','constructed_order','gamma','gamma_cer',
                'bipartite','matching_value','all_optima_directly_checked']
        writer=csv.DictWriter(stream,fieldnames=fields,extrasaction='ignore')
        writer.writeheader()
        writer.writerows(results)
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == '__main__':
    main()
