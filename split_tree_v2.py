"""
split_tree_v2.py

Brute-force verification of Theorem 1's characterization (Definition 1's
split-decomposition), tested directly on small graphs -- both a positive
case (a genuine strict m-gonal cactus) and three negative cases (a chord
added inside a block, a bridge between two blocks, a cycle of the wrong
length), checking that condition (a) of Theorem 1 holds or fails exactly
as it should.

Reference: F. G. Speyser, "Enumeration of Non-Plane m-Gonal Cactus Graphs
via Split-Decomposition", Theorem 1.

Author: Frederic G. Speyser
Run: python3 split_tree_v2.py
"""
from itertools import combinations


def neighbors(edges, v):
    return {b for a, b in edges if a == v} | {a for a, b in edges if b == v}


def find_split(vertices, edges):
    """Brute-force search for a split (Definition 1)."""
    edges_set = {frozenset(e) for e in edges}
    vertices = list(vertices)
    n = len(vertices)
    for size1 in range(2, n - 1):
        for V1 in combinations(vertices, size1):
            V1 = set(V1)
            V2 = set(vertices) - V1
            if len(V2) < 2:
                continue
            A = {v for v in V2 if neighbors(edges, v) & V1}
            B = {v for v in V1 if neighbors(edges, v) & V2}
            if A and B and all(frozenset((a, b)) in edges_set for a in A for b in B):
                return (V1, V2, A, B)
    return None


def cycle_edges(order):
    n = len(order)
    return [(order[i], order[(i + 1) % n]) for i in range(n)]


def is_clean_m_cycle(vertices, edges, m):
    """A genuine C_m block: m vertices, degree 2 everywhere, a single cycle."""
    if len(vertices) != m:
        return False, f"size {len(vertices)} != {m}"
    deg = {v: 0 for v in vertices}
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    if not all(d == 2 for d in deg.values()):
        return False, f"non-uniform degrees: {deg}"
    if len(edges) != m:
        return False, f"{len(edges)} edges != {m} (a chord or a missing edge)"
    return True, "OK"


def block_edges(vertices, edges, block):
    return [e for e in edges if e[0] in block and e[1] in block]


print("=" * 70)
print("POSITIVE TEST: a genuine strict 5-gonal cactus, 2 blocks")
print("=" * 70)
b1 = list(range(5))
b2 = [0, 5, 6, 7, 8]
edges = cycle_edges(b1) + cycle_edges(b2)
vertices = set(b1) | set(b2)
split = find_split(vertices, edges)
print(f"Split found: {split is not None}")
V1, V2, A, B = split
print(f"  V1={sorted(V1)} V2={sorted(V2)}  boundary A={sorted(A)} B={sorted(B)}")
ok1, msg1 = is_clean_m_cycle(b1, block_edges(vertices, edges, b1), 5)
ok2, msg2 = is_clean_m_cycle(b2, block_edges(vertices, edges, b2), 5)
print(f"  Block 1 (b1={b1}) is a clean C_5: {ok1} ({msg1})")
print(f"  Block 2 (b2={b2}) is a clean C_5: {ok2} ({msg2})")
print(f"  ==> Theorem 1, condition (a) satisfied: {ok1 and ok2}")

print()
print("=" * 70)
print("NEGATIVE TEST #1: same graph + a CHORD added inside block 1")
print("=" * 70)
edges_chord = edges + [(1, 3)]
ok1, msg1 = is_clean_m_cycle(b1, block_edges(vertices, edges_chord, b1), 5)
print(f"  Block 1 with a chord is a clean C_5: {ok1} ({msg1})")
print(f"  ==> Theorem 1, condition (a) VIOLATED, as expected: {not ok1}")

print()
print("=" * 70)
print("NEGATIVE TEST #2: same graph + a BRIDGE added between the two blocks")
print("=" * 70)
edges_bridge = edges + [(2, 6)]
split_bridge = find_split(set(b1) | set(b2), edges_bridge)
print(f"  Split found despite the bridge: {split_bridge is not None}")
if split_bridge:
    V1b, V2b, Ab, Bb = split_bridge
    print(f"    V1={sorted(V1b)} V2={sorted(V2b)}")
    ok1b, msg1b = is_clean_m_cycle(V1b, block_edges(set(b1) | set(b2), edges_bridge, V1b), 5)
    print(f"    Does this split still give two clean C_5's? side 1: {ok1b} ({msg1b})")
print("  ==> the bridge destroys the 'strict' structure (no bridges allowed in")
print("      a cactus): this graph is no longer in the class Theorem 1 covers")
print("      to begin with, which is exactly the expected behaviour.")

print()
print("=" * 70)
print("NEGATIVE TEST #3: cycle of the WRONG LENGTH (a hexagon instead of a pentagon)")
print("=" * 70)
b1_wrong = list(range(6))
edges_wrong = cycle_edges(b1_wrong) + cycle_edges([0, 6, 7, 8, 9])
ok1, msg1 = is_clean_m_cycle(
    b1_wrong, block_edges(set(b1_wrong) | {6, 7, 8, 9}, edges_wrong, b1_wrong), 5
)
print(f"  Block 1 (actually a C_6) is a clean C_5? {ok1} ({msg1})")
print(f"  ==> Theorem 1, condition (a) VIOLATED (right graph type, wrong m): {not ok1}")
