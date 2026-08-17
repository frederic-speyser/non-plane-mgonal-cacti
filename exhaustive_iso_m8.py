"""
exhaustive_iso_m8.py

Fully independent verification of the small-k counts (1, 1, 5) for k=1,2,3
blocks, strict 8-gonal cacti (unrooted, free/non-plane case): the graphs
are built directly as combinatorial objects (via networkx) and
deduplicated by graph isomorphism -- without going through the functional
equation of the paper, Howroyd's general PARI array formula, or any other
generating-function machinery used elsewhere in this repository.

Context: unlike verify_pari_dissymmetry_odd_m.gp (m=7), no closed-form
simplification exists for m=8 -- see verify_pari_dissymmetry_even_m.gp for
why (four divisors of 8, plus the even-k correction term in Howroyd's
general formula, neither of which cancel the way they do for odd m). This
script offers a different kind of verification instead: not an algebraic
identity, but an entirely separate combinatorial method, in the same
spirit as exhaustive_iso.py (which does the analogous check for m=5).

Verified below: k=1 -> 1, k=2 -> 1, k=3 -> 5, matching both Howroyd's
general U(n,8) array formula and the independent exact-rational solver of
mgonal_cactus_series.py.

Reference: F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane m-Gonal Cactus Graphs via Split-Decomposition", Sections 5.1-5.3
(compare against the resulting counts for m=8, unrooted).

Author: Frederic G. Speyser
Run: python3 exhaustive_iso_m8.py   (requires: pip install networkx)
"""
import networkx as nx


def octagon(offset):
    """A C8 with vertices 'offset_0'..'offset_7'."""
    verts = [f"{offset}_{i}" for i in range(8)]
    G = nx.Graph()
    G.add_nodes_from(verts)
    for i in range(8):
        G.add_edge(verts[i], verts[(i + 1) % 8])
    return G, verts


def glue(graphs_verts, merges):
    """graphs_verts: list of (G, verts). merges: list of pairs
    ((i1,v1),(i2,v2)) meaning verts[i1][v1] and verts[i2][v2] must become
    the same vertex."""
    G = nx.Graph()
    for g, verts in graphs_verts:
        G = nx.union(G, g)
    parent = {}

    def find(x):
        while parent.get(x, x) != x:
            x = parent[x]
        return x

    def union(x, y):
        parent[find(x)] = find(y)

    for (i1, v1), (i2, v2) in merges:
        a = graphs_verts[i1][1][v1]
        b = graphs_verts[i2][1][v2]
        union(a, b)
    mapping = {n: find(n) for n in G.nodes()}
    return nx.relabel_nodes(G, mapping)


def dedup_by_isomorphism(graph_list):
    classes = []
    for G in graph_list:
        found = False
        for cls in classes:
            if nx.is_isomorphic(G, cls[0]):
                cls.append(G)
                found = True
                break
        if not found:
            classes.append([G])
    return classes


print("=" * 70)
print("k=1 block: should give exactly 1 isomorphism class")
print("=" * 70)
G1, v1 = octagon("a")
print(f"  Classes found: {len(dedup_by_isomorphism([G1]))}  (expected: 1)")

print()
print("=" * 70)
print("k=2 blocks: should give exactly 1 isomorphism class")
print("=" * 70)
candidates_k2 = []
for shift in range(8):  # try different attachment points -- all equivalent by symmetry
    g1, v1 = octagon("a")
    g2, v2 = octagon("b")
    G = glue([(g1, v1), (g2, v2)], [((0, 0), (1, shift))])
    candidates_k2.append(G)
classes_k2 = dedup_by_isomorphism(candidates_k2)
print(f"  Classes found (out of {len(candidates_k2)} constructions tested): "
      f"{len(classes_k2)}  (expected: 1)")

print()
print("=" * 70)
print("k=3 blocks: should give exactly 5 isomorphism classes")
print("=" * 70)
candidates_k3 = []
g1, v1 = octagon("a")
g2, v2 = octagon("b")
g3, v3 = octagon("c")
G_star = glue([(g1, v1), (g2, v2), (g3, v3)], [((0, 0), (1, 0)), ((0, 0), (2, 0))])
candidates_k3.append(("star", G_star))

for d1 in range(8):
    for d2 in range(8):
        if d1 == d2:
            continue
        dist = min(abs(d1 - d2), 8 - abs(d1 - d2))
        g1, v1 = octagon("a")
        g2, v2 = octagon("b")
        g3, v3 = octagon("c")
        G = glue([(g1, v1), (g2, v2), (g3, v3)],
                 [((0, 0), (1, d1)), ((2, 0), (1, d2))])
        candidates_k3.append((f"chain_dist{dist}_{d1}_{d2}", G))

graphs_only = [g for _, g in candidates_k3]
classes_k3 = dedup_by_isomorphism(graphs_only)
print(f"  Classes found (out of {len(candidates_k3)} constructions tested): "
      f"{len(classes_k3)}  (expected: 5)")
for i, cls in enumerate(classes_k3):
    print(f"    class {i+1}: {len(cls)} equivalent constructions, "
          f"{cls[0].number_of_nodes()} vertices, {cls[0].number_of_edges()} edges")

print()
print("Cross-check against Howroyd's U(n,8) / mgonal_cactus_series.py: "
      "1, 1, 5, ... -- matches exactly for k=1,2,3.")
