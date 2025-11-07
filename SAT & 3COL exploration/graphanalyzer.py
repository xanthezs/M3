#!/usr/bin/env python3
"""
graph_degree_detail.py
----------------------
For every *.txt graph passed on the command line (or every *.txt file in the
current directory if none are provided) this script prints:

* total vertices / edges
* exact counts for vertices with degree 0, 1, 2, …, 30
* a single bucket for degrees 31-99
* a list of every vertex whose degree ≥100 with its exact degree

Input format per file (same as earlier examples):
    first line:  N  M
    next M lines:  u  v   (undirected edge)

Usage
-----
    python3 graph_degree_detail.py          # scans all *.txt in cwd
    python3 graph_degree_detail.py G1.txt … # explicit files
"""

import sys
import glob
from collections import defaultdict

MAX_EXACT = 30             # print individual counts for 0..30

def analyse(path: str) -> None:
    with open(path) as f:
        N, M_claimed = map(int, f.readline().split())
        deg = [0]*(N+1)

        # Use a set to track unique edges (avoid double-counting duplicates)
        unique_edges = set()
        lines_read = 0

        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) != 2:
                continue

            try:
                u, v = map(int, parts)
            except ValueError:
                print(f"  Warning: Skipping invalid line: {line}")
                continue

            lines_read += 1

            # Validate vertices are in range
            if not (1 <= u <= N and 1 <= v <= N):
                print(f"  Warning: Edge ({u}, {v}) has vertex outside range 1-{N}")
                continue

            # Skip self-loops
            if u == v:
                print(f"  Warning: Skipping self-loop ({u}, {v})")
                continue

            # Normalize edge representation (smaller vertex first)
            edge = (min(u, v), max(u, v))

            # Only count each unique edge once
            if edge not in unique_edges:
                unique_edges.add(edge)
                deg[u] += 1
                deg[v] += 1

    # Calculate actual edge count
    M_actual = len(unique_edges)

    exact = defaultdict(int)   # degree → count (for 0..30)
    over = 0                   # 31..99
    hi_vertices = []           # (id, degree) for ≥100

    for v in range(1, N+1):
        d = deg[v]
        if d <= MAX_EXACT:
            exact[d] += 1
        elif d < 100:
            over += 1
        else:
            hi_vertices.append((v, d))

    # ────────── report ──────────
    print(f"\nFile: {path}")
    print(f"  vertices: {N:,}, edges: {M_actual:,}", end="")
    if M_actual != M_claimed:
        print(f" (file header claimed {M_claimed:,})")
    else:
        print()

    if lines_read != M_claimed:
        print(f"  Warning: Expected {M_claimed:,} edge lines, read {lines_read:,}")

    # Verify edge count using degree sum
    degree_sum = sum(deg[1:])
    if degree_sum != 2 * M_actual:
        print(f"  ERROR: Degree sum {degree_sum} != 2 * edges {2 * M_actual}")

    for d in range(0, MAX_EXACT+1):
        count = exact.get(d, 0)
        if count > 0:  # Only print non-zero counts
            print(f"  deg={d:2}: {count:,}")
    print(f"  31-99 : {over:,}")
    print(f"  >=100 : {len(hi_vertices):,}")

    if hi_vertices:
        print("    vertices with degree ≥100 (id → degree):")
        for v, d in sorted(hi_vertices, key=lambda x: (-x[1], x[0])):
            print(f"      {v:>10,} → {d:,}")

def main():
    files = sys.argv[1:] or glob.glob("*.txt")
    if not files:
        print("No .txt graph files found.")
        return
    for path in files:
        try:
            analyse(path)
        except Exception as e:
            print(f"\nError processing {path}: {e}")

if __name__ == "__main__":
    main()
