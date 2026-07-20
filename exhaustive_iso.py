"""
exhaustive_iso.py

Fully independent verification of the small-k counts (1, 1, 3) for k=1,2,3
blocks (strict 5-gonal cacti, free/non-plane case): the graphs are built
directly as combinatorial objects (via networkx) and deduplicated by graph
isomorphism -- without going through the functional equation of the paper,
or any other code used elsewhere in this repository.

Reference: F. G. Speyser, "Enumeration of Non-Plane m-Gonal Cactus Graphs
via Split-Decomposition", Sections 5.1-5.3 (compare against the resulting
counts s_n for m=5).

Author: Frederic G. Speyser
Run: python3 exhaustive_iso.py   (requires: pip install networkx)
"""
import networkx as nx
from itertools import combinations


def pentagon(offset):
    """A C5 with vertices 'offset_0'..'offset_4'."""
    verts = [f"{offset}_{i}" for i in range(5)]
    G = nx.Graph()
    G.add_nodes_from(verts)
    for i in range(5):
        G.add_edge(verts[i], verts[(i + 1) % 5])
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
G1, v1 = pentagon("a")
print(f"  Classes found: {len(dedup_by_isomorphism([G1]))}  (expected: 1)")

print()
print("=" * 70)
print("k=2 blocks: should give exactly 1 isomorphism class")
print("=" * 70)
candidates_k2 = []
for shift in range(5):  # try different attachment points -- all equivalent by symmetry
    g1, v1 = pentagon("a")
    g2, v2 = pentagon("b")
    G = glue([(g1, v1), (g2, v2)], [((0, 0), (1, shift))])
    candidates_k2.append(G)
classes_k2 = dedup_by_isomorphism(candidates_k2)
print(f"  Classes found (out of {len(candidates_k2)} constructions tested): "
      f"{len(classes_k2)}  (expected: 1)")

print()
print("=" * 70)
print("k=3 blocks: should give exactly 3 isomorphism classes")
print("=" * 70)
candidates_k3 = []
# (i) star: all 3 blocks share a SINGLE cut vertex
g1, v1 = pentagon("a")
g2, v2 = pentagon("b")
g3, v3 = pentagon("c")
G_star = glue([(g1, v1), (g2, v2), (g3, v3)], [((0, 0), (1, 0)), ((0, 0), (2, 0))])
candidates_k3.append(("star", G_star))

# (ii) chain: the two cut vertices of the middle block, at every relative offset
for d1 in range(5):
    for d2 in range(5):
        if d1 == d2:
            continue
        dist = min(abs(d1 - d2), 5 - abs(d1 - d2))
        g1, v1 = pentagon("a")
        g2, v2 = pentagon("b")
        g3, v3 = pentagon("c")
        G = glue([(g1, v1), (g2, v2), (g3, v3)],
                 [((0, 0), (1, d1)), ((2, 0), (1, d2))])
        candidates_k3.append((f"chain_dist{dist}_{d1}_{d2}", G))

graphs_only = [g for _, g in candidates_k3]
classes_k3 = dedup_by_isomorphism(graphs_only)
print(f"  Classes found (out of {len(candidates_k3)} constructions tested): "
      f"{len(classes_k3)}  (expected: 3)")
for i, cls in enumerate(classes_k3):
    print(f"    class {i+1}: {len(cls)} equivalent constructions, "
          f"{cls[0].number_of_nodes()} vertices, {cls[0].number_of_edges()} edges")
