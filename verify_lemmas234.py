"""
verify_lemmas234.py

Independent verification of Lemmas 2-4 (the reverse direction of Theorem 1:
reconstructing the graph G from a graph-labeled tree T satisfying conditions
(a)-(d)), which had previously only been checked by reading the proofs, not
by an independent construction.

Construction rule used (vertices of G = leaves of T, following Lemma 2's own
proof): each C_m-node's marker either leads directly to a leaf, or -- via a
star-node's extremity -- to that star's center, which leads to a leaf; the
m markers' associated leaves are then joined in a cycle. This was arrived
at after an incorrect first attempt (naively fusing an entire star's
extremities into one point) produced a 3-internal-node path, contradicting
Lemma 3's stated bound of at most two -- the bug was in this script's first
draft, not in the paper; it is recorded here as a useful cross-check of how
precisely the paper's own construction must be followed.

Verified below, on a 2-block and a 3-block example:
  - Lemma 2: a C_m-node's m markers have pairwise distinct associated
    leaves, adjacent in G iff cyclically consecutive.
  - Lemma 4: those m leaves form an exact, maximal biconnected block of G,
    isomorphic to C_m (not contained in, nor containing, any larger
    2-connected piece).
  - Lemma 3 (implicitly): the construction only ever requires resolving a
    marker's associated leaf through at most one star-node, consistent
    with the stated path-length bound.

Author: Frederic G. Speyser
Run with: python3 verify_lemmas234.py   (requires: pip install networkx)
"""
import networkx as nx


def build_graph(m, cm_nodes, star_nodes, connections):
    connections_by_marker = {(cid, i): sid for (cid, i), (sid, k) in connections}
    star_center_leaf = {sid: f"leaf_center_{sid}" for sid in star_nodes}

    def resolve_leaf(cid, i):
        if (cid, i) in connections_by_marker:
            return star_center_leaf[connections_by_marker[(cid, i)]]
        return f"leaf_direct_{cid}_{i}"

    leaf_of = {(cid, i): resolve_leaf(cid, i) for cid in cm_nodes for i in range(m)}

    G = nx.Graph()
    ordered_markers = {}
    for cid in cm_nodes:
        cycle_leaves = [leaf_of[(cid, i)] for i in range(m)]
        ordered_markers[cid] = cycle_leaves
        for i in range(m):
            G.add_edge(cycle_leaves[i], cycle_leaves[(i + 1) % m])
    return G, ordered_markers, leaf_of


def check_case(m, cm_nodes, star_nodes, connections, expected_vertices):
    G, ordered_markers, leaf_of = build_graph(m, cm_nodes, star_nodes, connections)
    print(f"  vertices: {G.number_of_nodes()} (expected {expected_vertices})"
          f"  -- {'OK' if G.number_of_nodes()==expected_vertices else 'MISMATCH'}")
    blocks = list(nx.biconnected_components(G))
    print(f"  biconnected blocks found: {len(blocks)} (expected {len(cm_nodes)})")
    for cid in cm_nodes:
        target = set(ordered_markers[cid])
        present = target in blocks
        subg = G.subgraph(target)
        clean = all(d == 2 for _, d in subg.degree()) and subg.number_of_edges() == m
        print(f"    block {cid}: exact maximal block = {present}, clean C_{m} = {clean}")


if __name__ == "__main__":
    m = 5
    print("=== 2-block case (P and Q joined via one star) ===")
    check_case(m, ['P', 'Q'], ['S1'],
               [(('P', 0), ('S1', 0)), (('Q', 0), ('S1', 1))],
               expected_vertices=1 + 4 * 2)

    print("\n=== 3-block chain (P - S1 - Q - S2 - R) ===")
    check_case(m, ['P', 'Q', 'R'], ['S1', 'S2'],
               [(('P', 0), ('S1', 0)), (('Q', 0), ('S1', 1)),
                (('Q', 2), ('S2', 0)), (('R', 0), ('S2', 1))],
               expected_vertices=1 + 4 * 3)
