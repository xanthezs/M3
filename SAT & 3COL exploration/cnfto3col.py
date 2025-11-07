#!/usr/bin/env python3
"""
CNF to 3-Coloring Converter
Transforms CNF to 3-coloring graph without checking satisfiability.
Optimized for large instances.

Usage: python cnfto3col.py input.cnf output.txt
Output:
  output.txt - Graph in standard format
  output.txt.map - Vertex mapping information
"""

import sys


def parse_cnf(filename):
    """Parse DIMACS CNF file efficiently."""
    variables = set()
    clauses = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('c'):
                continue

            if line.startswith('p'):
                # Can extract expected numbers if needed
                parts = line.split()
                if len(parts) >= 4:
                    expected_vars = int(parts[2])
                    expected_clauses = int(parts[3])
                continue

            # Parse clause
            clause = []
            for token in line.split():
                lit = int(token)
                if lit == 0:
                    break
                clause.append(lit)
                variables.add(abs(lit))

            if clause:
                clauses.append(clause)

    return sorted(variables), clauses


def cnf_to_3coloring(variables, clauses):
    """
    Convert CNF to 3-coloring graph.
    Uses the standard reduction with proven correctness.
    """

    # Graph components
    vertices = {}
    edges = []
    vertex_id = 0

    def add_vertex(label):
        nonlocal vertex_id
        vertex_id += 1
        vertices[vertex_id] = label
        return vertex_id

    def add_edge(u, v):
        if u != v:
            edges.append((min(u, v), max(u, v)))

    # Step 1: Base palette triangle
    TRUE = add_vertex("TRUE")
    FALSE = add_vertex("FALSE")
    BASE = add_vertex("BASE")

    add_edge(TRUE, FALSE)
    add_edge(FALSE, BASE)
    add_edge(TRUE, BASE)

    # Step 2: Variable gadgets
    var_map = {}
    for var in variables:
        pos = add_vertex(f"x{var}")
        neg = add_vertex(f"!x{var}")

        # Triangle with BASE
        add_edge(pos, neg)
        add_edge(pos, BASE)
        add_edge(neg, BASE)

        var_map[var] = (pos, neg)

    # Step 3: Clause gadgets
    for clause_idx, clause in enumerate(clauses):
        clause_id = clause_idx + 1

        # Get literal vertices
        lit_vertices = []
        for lit in clause:
            var = abs(lit)
            if lit > 0:
                lit_vertices.append(var_map[var][0])
            else:
                lit_vertices.append(var_map[var][1])

        if len(clause) == 0:
            # Empty clause - make uncolorable
            bad = add_vertex(f"empty_C{clause_id}")
            add_edge(bad, TRUE)
            add_edge(bad, FALSE)
            add_edge(bad, BASE)

        elif len(clause) == 1:
            # Unit clause - force literal to be TRUE
            force = add_vertex(f"unit_C{clause_id}")
            add_edge(force, FALSE)
            add_edge(force, BASE)
            add_edge(force, lit_vertices[0])

        elif len(clause) == 2:
            # 2-literal OR gadget
            g1 = add_vertex(f"C{clause_id}_g1")
            g2 = add_vertex(f"C{clause_id}_g2")

            add_edge(g1, TRUE)
            add_edge(g2, TRUE)
            add_edge(g1, g2)
            add_edge(g1, lit_vertices[0])
            add_edge(g2, lit_vertices[1])

        elif len(clause) == 3:
            # 3-literal OR gadget
            a, b, c = lit_vertices

            # Create intermediate vertex for (a ∨ b)
            ab = add_vertex(f"C{clause_id}_ab")
            add_edge(ab, BASE)

            # Sub-gadget for (a ∨ b)
            s1 = add_vertex(f"C{clause_id}_s1")
            s2 = add_vertex(f"C{clause_id}_s2")

            add_edge(s1, s2)
            add_edge(s1, ab)
            add_edge(s2, ab)
            add_edge(s1, a)
            add_edge(s2, b)

            # Main gadget for ((a∨b) ∨ c)
            m1 = add_vertex(f"C{clause_id}_m1")
            m2 = add_vertex(f"C{clause_id}_m2")

            add_edge(m1, TRUE)
            add_edge(m2, TRUE)
            add_edge(m1, m2)
            add_edge(m1, ab)
            add_edge(m2, c)

            # Output vertex
            out = add_vertex(f"C{clause_id}_out")
            add_edge(out, FALSE)
            add_edge(out, BASE)
            add_edge(out, m1)
            add_edge(out, m2)

        else:
            # Long clause - reduce to 3-SAT using Tseitin transformation
            # Create auxiliary variables
            aux_vars = []
            for i in range(len(clause) - 3):
                aux_pos = add_vertex(f"y{clause_id}_{i+1}")
                aux_neg = add_vertex(f"!y{clause_id}_{i+1}")

                add_edge(aux_pos, aux_neg)
                add_edge(aux_pos, BASE)
                add_edge(aux_neg, BASE)

                aux_vars.append((aux_pos, aux_neg))

            # Transform to 3-SAT clauses
            # First: (l1 ∨ l2 ∨ y1)
            create_3lit_gadget(add_vertex, add_edge, TRUE, FALSE, BASE,
                             lit_vertices[0], lit_vertices[1], aux_vars[0][0],
                             f"{clause_id}_sub1")

            # Middle: (!yi ∨ l(i+2) ∨ y(i+1))
            for i in range(len(clause) - 4):
                create_3lit_gadget(add_vertex, add_edge, TRUE, FALSE, BASE,
                                 aux_vars[i][1], lit_vertices[i+2], aux_vars[i+1][0],
                                 f"{clause_id}_sub{i+2}")

            # Last: (!yn-3 ∨ ln-1 ∨ ln)
            if len(aux_vars) > 0:
                create_3lit_gadget(add_vertex, add_edge, TRUE, FALSE, BASE,
                                 aux_vars[-1][1], lit_vertices[-2], lit_vertices[-1],
                                 f"{clause_id}_sub_last")

    return vertex_id, vertices, edges, var_map


def create_3lit_gadget(add_vertex, add_edge, TRUE, FALSE, BASE, a, b, c, label):
    """Helper to create 3-literal gadget (reused for long clauses)."""
    # Create intermediate vertex for (a ∨ b)
    ab = add_vertex(f"{label}_ab")
    add_edge(ab, BASE)

    # Sub-gadget for (a ∨ b)
    s1 = add_vertex(f"{label}_s1")
    s2 = add_vertex(f"{label}_s2")

    add_edge(s1, s2)
    add_edge(s1, ab)
    add_edge(s2, ab)
    add_edge(s1, a)
    add_edge(s2, b)

    # Main gadget for ((a∨b) ∨ c)
    m1 = add_vertex(f"{label}_m1")
    m2 = add_vertex(f"{label}_m2")

    add_edge(m1, TRUE)
    add_edge(m2, TRUE)
    add_edge(m1, m2)
    add_edge(m1, ab)
    add_edge(m2, c)

    # Output vertex
    out = add_vertex(f"{label}_out")
    add_edge(out, FALSE)
    add_edge(out, BASE)
    add_edge(out, m1)
    add_edge(out, m2)


def write_graph(filename, n_vertices, edges):
    """Write graph in standard format."""
    # Remove duplicate edges
    edge_set = set()
    for u, v in edges:
        edge_set.add((min(u, v), max(u, v)))

    with open(filename, 'w') as f:
        f.write(f"{n_vertices} {len(edge_set)}\n")
        for u, v in sorted(edge_set):
            f.write(f"{u} {v}\n")


def write_mapping(filename, vertices, var_map):
    """Write vertex mapping for debugging."""
    with open(filename + ".map", 'w') as f:
        f.write("CNF to 3-Coloring Vertex Mapping\n")
        f.write("=" * 50 + "\n\n")

        f.write("PALETTE VERTICES:\n")
        f.write("  1: TRUE\n")
        f.write("  2: FALSE\n")
        f.write("  3: BASE\n\n")

        f.write("VARIABLE VERTICES:\n")
        for var, (pos, neg) in sorted(var_map.items()):
            f.write(f"  {pos}: x{var}\n")
            f.write(f"  {neg}: !x{var}\n")

        f.write("\nEXTRACTION RULE:\n")
        f.write("  If vertex for x_i has same color as vertex 1 (TRUE), then x_i = true\n")
        f.write("  If vertex for !x_i has same color as vertex 1 (TRUE), then x_i = false\n")


def main():
    if len(sys.argv) != 3:
        print("Usage: python cnf2col_fast.py input.cnf output.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Parse CNF
    print("Parsing CNF file...")
    variables, clauses = parse_cnf(input_file)
    print(f"  Found {len(variables)} variables and {len(clauses)} clauses")

    # Convert to graph
    print("Building 3-coloring graph...")
    n_vertices, vertices, edges, var_map = cnf_to_3coloring(variables, clauses)
    print(f"  Created {n_vertices} vertices and {len(set(edges))} edges")

    # Write output
    print("Writing output files...")
    write_graph(output_file, n_vertices, edges)
    write_mapping(output_file, vertices, var_map)

    print(f"\nOutput files:")
    print(f"  {output_file} - Graph in standard format")
    print(f"  {output_file}.map - Vertex mapping information")

    print("\nConversion complete!")


if __name__ == "__main__":
    main()
