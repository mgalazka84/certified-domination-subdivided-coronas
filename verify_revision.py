"""Additional reproducible checks for the revised DML manuscript.

Hardness statements are proved in the manuscript, not inferred from these tests.
The gadget count is checked on base subsets; full-graph MILP independently
checks its optimum value. Constrained min-cut listing is compared with exhaustive
base-subset minimization and with the defining certification constraints.
"""
from __future__ import annotations
import itertools as it
import json
from collections import deque
from pathlib import Path
from verify import (unlabelled_graphs, construct, is_certified, exact_model,
                    base_formula, canonical, bipartition, bounded_edges,
                    direct_all_optima)


def degrees(n, edges):
    d = [0]*n
    for u, v in edges:
        d[u] += 1
        d[v] += 1
    return d


def independent_masks(n, edges):
    return [X for X in range(1 << n)
            if all(not (X >> u & 1 and X >> v & 1) for u, v in edges)]


def constrained_cut(n, edges, b, colors, prefix):
    """Max flow with forced source/sink assignments; prefix True means source."""
    source, sink = n, n+1
    cap = [[0]*(n+2) for _ in range(n+2)]
    for v in range(n):
        if colors[v] == 0:
            cap[source][v] += b[v]
        else:
            cap[v][sink] += b[v]
    for u, v in edges:
        if colors[u] == 1:
            u, v = v, u
        cap[u][v] += 1
    for v, side in enumerate(prefix):
        if side:
            cap[source][v] += len(edges)+1
        else:
            cap[v][sink] += len(edges)+1
    flow = 0
    while True:
        parent = [-1]*(n+2)
        parent[source] = source
        queue = deque([source])
        while queue and parent[sink] == -1:
            u = queue.popleft()
            for v, capacity in enumerate(cap[u]):
                if capacity and parent[v] == -1:
                    parent[v] = u
                    queue.append(v)
        if parent[sink] == -1:
            return flow
        v, amount = sink, float('inf')
        while v != source:
            u = parent[v]
            amount = min(amount, cap[u][v])
            v = u
        v = sink
        while v != source:
            u = parent[v]
            cap[u][v] -= amount
            cap[v][u] += amount
            v = u
        flow += amount


def enumerate_min_cuts(n, edges, r):
    colors = bipartition(n, edges)
    assert colors is not None
    b = [v+1 for v in r]
    optimum = constrained_cut(n, edges, b, colors, ())
    calls, output = 1, []
    def visit(prefix):
        nonlocal calls
        if len(prefix) == n:
            X = sum(1 << v for v in range(n)
                    if (colors[v] == 0 and not prefix[v])
                    or (colors[v] == 1 and prefix[v]))
            output.append(X)
            return
        for side in (False, True):
            child = prefix+(side,)
            calls += 1
            if constrained_cut(n, edges, b, colors, child) == optimum:
                visit(child)
    visit(())
    return optimum, output, calls


def is_forest(n, edges):
    parent = list(range(n))
    def root(v):
        while parent[v] != v:
            v = parent[v]
        return v
    for u, v in edges:
        a, b = root(u), root(v)
        if a == b:
            return False
        parent[a] = b
    return True


def main():
    output = {'status': 'all passed', 'hardness_identity': [],
              'counting_gadget': [], 'forest_identity': [],
              'minimum_cut_listing': [], 'star_values': [], 'path_one': {}}
    for n in range(1, 7):
        for edges in unlabelled_graphs(n):
            deg = degrees(n, edges)
            if min(deg) >= 3:
                r = [d-2 for d in deg]
                alpha = max(X.bit_count() for X in independent_masks(n, edges))
                expected = sum(r)+len(edges)-alpha
                value, xs, _ = base_formula(n, edges, r)
                adj, *_ = construct(n, edges, r)
                actual = exact_model(adj)
                assert actual == value == expected
                for X in range(1 << n):
                    f = sum(r[v]+1 for v in range(n) if X >> v & 1)
                    f += sum(not (X >> u & 1 or X >> v & 1) for u, v in edges)
                    internal = sum(bool(X >> u & 1 and X >> v & 1) for u,v in edges)
                    assert f == len(edges)-X.bit_count()+internal
                output['hardness_identity'].append(dict(n=n, edges=edges, r=r,
                    alpha=alpha, gamma_cer=actual, cubic=all(d == 3 for d in deg)))
            if is_forest(n, edges):
                # Exhaustive edge subsets give an independent path-cover computation.
                best = 0
                for mask in range(1 << len(edges)):
                    subset = [e for j,e in enumerate(edges) if mask >> j & 1]
                    if max(degrees(n, subset)) <= 2:
                        best = max(best, len(subset))
                pc = n-best
                adj, *_ = construct(n, edges, [1]*n)
                actual = exact_model(adj)
                assert actual == 2*n-pc
                output['forest_identity'].append(dict(n=n, edges=edges,
                                                       path_cover=pc, gamma_cer=actual))
            colors = bipartition(n, edges)
            if colors is not None and n <= 5:
                assignments = it.product((1,2), repeat=n) if n <= 4 else [[1]*n]
                for r in assignments:
                    optimum, xs, calls = enumerate_min_cuts(n, edges, r)
                    expected, expected_xs, _ = base_formula(n, edges, r)
                    assert optimum+sum(r) == expected
                    assert len(xs) == len(set(xs)) and set(xs) == set(expected_xs)
                    adj, supports, leaves, middle = construct(n, edges, r)
                    for X in xs:
                        D = canonical(X, supports, leaves, middle, edges)
                        assert len(D) == expected and is_certified(adj, D)
                    output['minimum_cut_listing'].append(dict(n=n, edges=edges,
                        r=list(r), min_cut=optimum, minimizing_masks=xs, flow_calls=calls))
            if colors is not None and n <= 4:
                # H is Q plus two new leaves per vertex; original vertices stay 0..n-1.
                h_edges = list(edges)+[(v, n+2*v+i) for v in range(n) for i in (0,1)]
                r = [d+1 for d in deg]+[1]*(2*n)
                expected_masks = set(independent_masks(n, edges))
                value, xs, _ = base_formula(3*n, h_edges, r)
                assert set(xs) == expected_masks
                assert bipartition(3*n, h_edges) is not None
                adj, supports, leaves, middle = construct(3*n, h_edges, r)
                assert len(adj) == 11*n+5*len(edges)
                assert exact_model(adj) == value == sum(r)+len(h_edges)
                for X in xs:
                    D = canonical(X, supports, leaves, middle, h_edges)
                    assert is_certified(adj,D) and len(D) == value
                direct = n == 1
                if direct:
                    _, actual, ds = direct_all_optima(adj)
                    predicted = {sum(1 << v for v in canonical(X,supports,leaves,middle,h_edges))
                                 for X in xs}
                    assert actual == value and ds == predicted
                output['counting_gadget'].append(dict(q_order=n, q_edges=edges,
                    h_order=3*n, h_edges=h_edges, r=r, constructed_order=len(adj),
                    independent_masks=sorted(expected_masks), gamma_cer=value,
                    minimum_set_count=len(xs), full_subset_enumeration=direct))
    for q in range(1, 14):
        edges = [(0,v) for v in range(1,q+1)]
        value, xs, _ = base_formula(q+1,edges,[1]*(q+1))
        assert value == q+1+min(q,2)
        output['star_values'].append(dict(q=q,gamma_cer=value,count=len(xs)))
    adj, *_ = construct(1,[],[1])
    ordinary, value, ds = direct_all_optima(adj)
    assert value == len(ds) == 1
    output['path_one'] = dict(gamma_cer=value,count=len(ds),ordinary=ordinary)
    counts = {k:len(v) for k,v in output.items() if isinstance(v,list)}
    counts['cubic_bases'] = sum(v['cubic'] for v in output['hardness_identity'])
    counts['path_one'] = 1
    output['counts'] = counts
    path = Path(__file__).resolve().parent/'results'/'revision_checks.json'
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'status':output['status'],'counts':counts},indent=2))


if __name__ == '__main__':
    main()
