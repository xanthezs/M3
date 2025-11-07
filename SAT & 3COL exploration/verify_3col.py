#!/usr/bin/env python3
"""
3-Coloring Verification Script

Simple verification that checks if adjacent vertices have different colors.

Usage: python verify_coloring.py graph_file coloring_file
graph_file - File containing the graph edges
coloring_file - File containing the vertex coloring (format: vertex color)

Prints verification result to console.
"""

import sys
import argparse


def read_graph_file(filename):
    edges = {}
    vertices_seen = set()
    edges_read = 0
    header_skipped = False
    first_few_edges = []

    with open(filename, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith('#') or line.startswith('c '):
                continue

            parts = line.split()

            # Better header detection - skip first line if it looks like metadata
            if line_num == 1 and len(parts) == 2:
                try:
                    num1, num2 = int(parts[0]), int(parts[1])
                    # More robust header detection
                    if num1 > 10 and num2 >= num1:  # Likely "vertices edges" format
                        header_skipped = True
                        print(f"Header detected and skipped: {num1} vertices, {num2} edges expected")
                        continue
                except ValueError:
                    pass

            # Handle edge lines - need exactly 2 parts
            if len(parts) >= 2:
                # Convert to integers for consistent comparison
                try:
                    vertex = int(parts[0])
                    nbr = int(parts[1])

                    # Skip self-loops
                    if vertex == nbr:
                        continue

                    # Convert back to strings for dictionary keys (consistent with coloring)
                    vertex_str = str(vertex)
                    nbr_str = str(nbr)

                    # Store first few edges for debugging
                    if len(first_few_edges) < 5:
                        first_few_edges.append(f"{vertex_str}-{nbr_str}")

                    # Add edge in both directions for undirected graph
                    if vertex_str not in edges:
                        edges[vertex_str] = []
                    edges[vertex_str].append(nbr_str)

                    if nbr_str not in edges:
                        edges[nbr_str] = []
                    edges[nbr_str].append(vertex_str)

                    vertices_seen.add(vertex_str)
                    vertices_seen.add(nbr_str)
                    edges_read += 1

                except ValueError:
                    continue

    print(f"Graph parsing complete:")
    print(f"  Header skipped: {header_skipped}")
    print(f"  Vertices: {len(vertices_seen)}")
    print(f"  Edges: {edges_read}")
    print(f"  First few edges: {first_few_edges}")

    return edges


def read_coloring_file(filename):
    coloring = {}
    first_few_colors = []

    with open(filename, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) >= 2:
                try:
                    # Convert to integer then back to string for consistency with graph
                    vertex = str(int(parts[0]))
                    color = parts[1]

                    if vertex in coloring:
                        return None

                    # Store first few for debugging
                    if len(first_few_colors) < 5:
                        first_few_colors.append(f"{vertex}:{color}")

                    coloring[vertex] = color

                except ValueError:
                    continue

    print(f"Coloring parsing complete:")
    print(f"  Vertices colored: {len(coloring)}")
    print(f"  First few assignments: {first_few_colors}")

    return coloring


def verify_3_coloring(graph_file, coloring_file):
    edges = read_graph_file(graph_file)
    coloring = read_coloring_file(coloring_file)

    if coloring is None:
        return False

    print("\nStarting verification...")

    # Check first few vertices for debugging
    vertices_checked = 0

    for vertex, nbrlist in edges.items():
        if vertex not in coloring:
            print(f"Error: Vertex {vertex} in graph but not in coloring")
            return False

        vertex_color = coloring[vertex]

        # Show debug info for first few vertices
        if vertices_checked < 3:
            print(f"Checking vertex {vertex} (color {vertex_color}) with {len(nbrlist)} neighbors: {nbrlist[:5]}...")

        for nbr in nbrlist:
            if nbr not in coloring:
                print(f"Error: Neighbor {nbr} of {vertex} not in coloring")
                return False

            if coloring[nbr] == vertex_color:
                print(f"❌ VERIFICATION FAILED")
                print(f"Adjacent vertices {vertex} and {nbr} both have color {vertex_color}")
                return False

        vertices_checked += 1

        # Stop debug output after first few vertices
        if vertices_checked == 3:
            print("(Continuing verification without debug output...)")

    print("✅ VERIFICATION SUCCESSFUL")
    return True


def main():
    parser = argparse.ArgumentParser(description="Verify a 3-coloring solution")
    parser.add_argument("graph_file", help="File containing the graph edges")
    parser.add_argument("coloring_file", help="File containing the vertex coloring")

    args = parser.parse_args()

    print("3-Coloring Verification Tool")
    print("=" * 40)

    try:
        result = verify_3_coloring(args.graph_file, args.coloring_file)
        return 0 if result else 1
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        return 2
    except Exception as e:
        print(f"Error: {e}")
        return 2


if __name__ == "__main__":
    sys.exit(main())
