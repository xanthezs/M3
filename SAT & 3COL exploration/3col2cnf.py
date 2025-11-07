#!/usr/bin/env python3
"""
3-Coloring to SAT Assignment Converter with File Output
Extracts SAT assignment from 3-coloring and saves results to files.

Usage: python 3col2cnf.py cnf_file graph_file coloring_file
Graph should be DIMACS format in .txt
Coloring should be in .txt.coloring format
"""

import sys
import os
from datetime import datetime


def parse_cnf_file(filename):
    """Parse a DIMACS CNF format file."""
    clauses = []
    num_vars = 0
    num_clauses = 0

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip comments and empty lines
            if not line or line.startswith('c'):
                continue

            # Parse problem line
            if line.startswith('p'):
                parts = line.split()
                num_vars = int(parts[2])
                num_clauses = int(parts[3])
                continue

            # Parse clause
            literals = []
            for token in line.split():
                lit = int(token)
                if lit == 0:  # End of clause
                    break
                literals.append(lit)

            if literals:
                clauses.append(literals)

    return num_vars, num_clauses, clauses


def parse_graph_file(filename):
    """Parse a graph file to understand the structure."""
    edges = []
    with open(filename, 'r') as f:
        first_line = f.readline().strip().split()
        n_vertices = int(first_line[0])
        n_edges = int(first_line[1])

        for line in f:
            u, v = map(int, line.strip().split())
            edges.append((u, v))

    return n_vertices, n_edges, edges


def parse_coloring_file(filename):
    """Parse a coloring file."""
    coloring = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split()
            vertex = int(parts[0])
            color = int(parts[1])
            coloring[vertex] = color

    return coloring


def detect_palette_mapping(coloring, num_vars):
    """
    Try to detect which palette vertex is TRUE/FALSE/BASE by checking constraints.
    Returns (true_vertex, false_vertex, base_vertex) or None if cannot determine.
    """
    # The three palette vertices form a triangle (vertices 1, 2, 3)
    palette_colors = {
        1: coloring.get(1),
        2: coloring.get(2),
        3: coloring.get(3)
    }

    # Check that palette has three different colors
    if len(set(palette_colors.values())) != 3:
        return None

    # For each variable, check which palette vertex corresponds to BASE
    # Variables form triangles with BASE: xi - !xi - BASE
    base_candidates = {1: 0, 2: 0, 3: 0}

    for var in range(1, min(num_vars + 1, 10)):  # Check first few variables
        pos_v = 3 + 2*(var-1) + 1
        neg_v = 3 + 2*(var-1) + 2

        pos_color = coloring.get(pos_v)
        neg_color = coloring.get(neg_v)

        if pos_color is None or neg_color is None:
            continue

        # Check which palette vertex could be BASE
        # BASE vertex should have different color from both xi and !xi
        for pv in [1, 2, 3]:
            if palette_colors[pv] != pos_color and palette_colors[pv] != neg_color:
                base_candidates[pv] += 1

    # The vertex with most votes is likely BASE
    base_vertex = max(base_candidates, key=base_candidates.get)

    # The other two are TRUE and FALSE
    tf_vertices = [v for v in [1, 2, 3] if v != base_vertex]

    # Try to determine which is TRUE and which is FALSE by checking unit clauses
    # or by convention (lower number = TRUE)
    true_vertex = tf_vertices[0]
    false_vertex = tf_vertices[1]

    return true_vertex, false_vertex, base_vertex


def extract_sat_assignment(coloring, num_vars, verbose=True, palette_override=None):
    """
    Extract SAT assignment from 3-coloring.

    palette_override: Optional tuple (true_vertex, false_vertex, base_vertex) to override detection
    """
    color_names = {0: "RED", 1: "GREEN", 2: "BLUE"}

    # Detect or use provided palette mapping
    if palette_override:
        true_vertex, false_vertex, base_vertex = palette_override
        if verbose:
            print("\n=== Using Manual Palette Override ===")
    else:
        detected = detect_palette_mapping(coloring, num_vars)
        if detected:
            true_vertex, false_vertex, base_vertex = detected
            if verbose:
                print("\n=== Auto-Detected Palette Mapping ===")
        else:
            # Fall back to standard assumption
            true_vertex, false_vertex, base_vertex = 1, 2, 3
            if verbose:
                print("\n=== Using Standard Palette Mapping (fallback) ===")

    # Get the colors
    true_color = coloring.get(true_vertex)
    false_color = coloring.get(false_vertex)
    base_color = coloring.get(base_vertex)

    if verbose:
        print(f"TRUE vertex: {true_vertex} (color {color_names.get(true_color, true_color)})")
        print(f"FALSE vertex: {false_vertex} (color {color_names.get(false_color, false_color)})")
        print(f"BASE vertex: {base_vertex} (color {color_names.get(base_color, base_color)})")

    # Store detailed mapping info
    mapping_info = {
        'palette': {
            'TRUE': {'vertex': true_vertex, 'color': true_color, 'color_name': color_names.get(true_color, str(true_color))},
            'FALSE': {'vertex': false_vertex, 'color': false_color, 'color_name': color_names.get(false_color, str(false_color))},
            'BASE': {'vertex': base_vertex, 'color': base_color, 'color_name': color_names.get(base_color, str(base_color))}
        },
        'variables': {}
    }

    # Verify palette forms a valid triangle (all different colors)
    if len({true_color, false_color, base_color}) != 3:
        print("ERROR: Palette vertices don't have three different colors!")
        return None, mapping_info

    assignment = {}

    if verbose:
        print("\n=== Variable Assignments ===")

    # Track statistics for debugging
    true_count = 0
    false_count = 0

    for var in range(1, num_vars + 1):
        # Calculate vertex numbers for this variable
        pos_vertex = 3 + 2*(var-1) + 1  # x_var
        neg_vertex = 3 + 2*(var-1) + 2  # !x_var

        pos_color = coloring.get(pos_vertex)
        neg_color = coloring.get(neg_vertex)

        if pos_color is None or neg_color is None:
            print(f"ERROR: Variable {var} vertices not found in coloring!")
            return None, mapping_info

        # Store variable mapping info
        mapping_info['variables'][var] = {
            'positive': {'vertex': pos_vertex, 'color': pos_color, 'color_name': color_names.get(pos_color, str(pos_color))},
            'negative': {'vertex': neg_vertex, 'color': neg_color, 'color_name': color_names.get(neg_color, str(neg_color))}
        }

        # Verify variable gadget constraints
        if pos_color == neg_color:
            print(f"ERROR: x{var} and !x{var} have the same color!")
            return None, mapping_info

        if pos_color == base_color or neg_color == base_color:
            print(f"ERROR: Variable {var} vertex has BASE color!")
            return None, mapping_info

        # Determine assignment based on which vertex has TRUE color
        if pos_color == true_color and neg_color == false_color:
            assignment[var] = True
            true_count += 1
            if verbose and var <= 5:  # Only show first few
                print(f"Variable {var}: x{var}=TRUE (vertex {pos_vertex} has TRUE color)")
        elif neg_color == true_color and pos_color == false_color:
            assignment[var] = False
            false_count += 1
            if verbose and var <= 5:  # Only show first few
                print(f"Variable {var}: x{var}=FALSE (!x{var} at vertex {neg_vertex} has TRUE color)")
        else:
            print(f"ERROR: Variable {var} has invalid color assignment!")
            print(f"  x{var} (vertex {pos_vertex}): color {pos_color}")
            print(f"  !x{var} (vertex {neg_vertex}): color {neg_color}")
            print(f"  Expected one to be TRUE color ({true_color}) and other FALSE color ({false_color})")
            return None, mapping_info

    if verbose and num_vars > 5:
        print(f"... (showing first 5 variables, {num_vars} total)")

    # Show statistics for debugging
    if verbose:
        print(f"\nAssignment statistics: {true_count} TRUE, {false_count} FALSE")
        if true_count == 1 and false_count == num_vars - 1:
            true_var = [v for v, val in assignment.items() if val][0]
            print(f"NOTE: Only x{true_var} is TRUE, all others are FALSE")
        elif false_count == 1 and true_count == num_vars - 1:
            false_var = [v for v, val in assignment.items() if not val][0]
            print(f"NOTE: Only x{false_var} is FALSE, all others are TRUE")
            print("This might indicate the TRUE/FALSE vertices are swapped!")

    return assignment, mapping_info


def verify_sat_assignment(clauses, assignment, verbose=True):
    """Verify that the assignment satisfies all CNF clauses."""
    if verbose:
        print("\n=== Verifying SAT Assignment ===")

    all_satisfied = True
    unsatisfied_clauses = []
    clause_results = []

    for i, clause in enumerate(clauses):
        clause_satisfied = False
        satisfied_literals = []

        for literal in clause:
            var = abs(literal)
            if literal > 0:
                # Positive literal
                if assignment.get(var, False):
                    clause_satisfied = True
                    satisfied_literals.append(f"x{var}")
            else:
                # Negative literal
                if not assignment.get(var, True):
                    clause_satisfied = True
                    satisfied_literals.append(f"!x{var}")

        clause_results.append({
            'clause_num': i + 1,
            'clause': clause,
            'satisfied': clause_satisfied,
            'satisfied_by': satisfied_literals
        })

        if verbose:
            status = "✓ SATISFIED" if clause_satisfied else "✗ NOT SATISFIED"
            print(f"Clause {i+1}: {clause} - {status}")
            if clause_satisfied and satisfied_literals:
                print(f"  Satisfied by: {', '.join(satisfied_literals)}")

        if not clause_satisfied:
            all_satisfied = False
            unsatisfied_clauses.append(i+1)

    if not all_satisfied and verbose:
        print(f"\nERROR: Clauses {unsatisfied_clauses} are not satisfied!")

    return all_satisfied, clause_results


def save_assignment_file(filename, assignment, mapping_info, satisfies_cnf, clause_results, cnf_file):
    """Save the extracted assignment and verification results to a file."""
    with open(filename, 'w') as f:
        f.write("# 3-Coloring to SAT Assignment Results\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Source CNF: {cnf_file}\n")
        f.write("#" + "="*60 + "\n\n")

        # Write palette information
        f.write("# PALETTE MAPPING\n")
        for name, info in mapping_info['palette'].items():
            f.write(f"# {name}: vertex {info['vertex']} = color {info['color']} ({info['color_name']})\n")
        f.write("\n")

        # Write variable mapping
        f.write("# VARIABLE VERTEX MAPPING\n")
        for var, info in mapping_info['variables'].items():
            f.write(f"# Variable {var}:\n")
            f.write(f"#   x{var} (vertex {info['positive']['vertex']}): color {info['positive']['color']} ({info['positive']['color_name']})\n")
            f.write(f"#   !x{var} (vertex {info['negative']['vertex']}): color {info['negative']['color']} ({info['negative']['color_name']})\n")
        f.write("\n")

        # Write assignment
        f.write("# EXTRACTED SAT ASSIGNMENT\n")
        for var in sorted(assignment.keys()):
            f.write(f"x{var} = {assignment[var]}\n")
        f.write("\n")

        # Write DIMACS format
        f.write("# DIMACS FORMAT\n")
        if satisfies_cnf:
            f.write("s SATISFIABLE\n")
        else:
            f.write("s UNSATISFIABLE\n")

        dimacs_vars = []
        for var in sorted(assignment.keys()):
            if assignment[var]:
                dimacs_vars.append(str(var))
            else:
                dimacs_vars.append(f"-{var}")
        f.write(f"v {' '.join(dimacs_vars)} 0\n\n")

        # Write clause verification results
        f.write("# CLAUSE VERIFICATION\n")
        for result in clause_results:
            status = "SATISFIED" if result['satisfied'] else "NOT SATISFIED"
            f.write(f"# Clause {result['clause_num']}: {result['clause']} - {status}\n")
            if result['satisfied'] and result['satisfied_by']:
                f.write(f"#   Satisfied by: {', '.join(result['satisfied_by'])}\n")
        f.write("\n")

        # Write final verdict
        f.write("#" + "="*60 + "\n")
        f.write(f"# FINAL VERDICT: {'SATISFIES CNF' if satisfies_cnf else 'DOES NOT SATISFY CNF'}\n")
        f.write("#" + "="*60 + "\n")


def main():
    if len(sys.argv) < 4:
        print("Usage: python col2sat_with_output.py cnf_file graph_file coloring_file [options]")
        print("Options:")
        print("  --quiet                  Suppress verbose output")
        print("  --swap-tf                Swap TRUE/FALSE vertices (use 2=TRUE, 1=FALSE)")
        print("  --palette T,F,B          Manual palette override (e.g., --palette 2,1,3)")
        print("Example: python col2sat_with_output.py REP3.cnf REP3.txt REP3.txt.coloring")
        sys.exit(1)

    cnf_file = sys.argv[1]
    graph_file = sys.argv[2]
    coloring_file = sys.argv[3]

    # Parse options
    verbose = "--quiet" not in sys.argv
    palette_override = None

    if "--swap-tf" in sys.argv:
        palette_override = (2, 1, 3)  # Swap TRUE and FALSE
        if verbose:
            print("Using --swap-tf: TRUE=vertex 2, FALSE=vertex 1")

    for i, arg in enumerate(sys.argv):
        if arg == "--palette" and i + 1 < len(sys.argv):
            try:
                parts = sys.argv[i + 1].split(',')
                palette_override = (int(parts[0]), int(parts[1]), int(parts[2]))
                if verbose:
                    print(f"Using manual palette: TRUE={parts[0]}, FALSE={parts[1]}, BASE={parts[2]}")
            except:
                print("ERROR: Invalid palette format. Use --palette T,F,B (e.g., --palette 2,1,3)")
                sys.exit(1)

    # Generate output filename
    base_name = os.path.splitext(coloring_file)[0]
    assignment_file = base_name + ".assignment"

    if verbose:
        print("\n=== 3-Coloring to SAT Assignment Converter ===")
        print(f"CNF file: {cnf_file}")
        print(f"Graph file: {graph_file}")
        print(f"Coloring file: {coloring_file}")
        print(f"Output file: {assignment_file}")

    # Parse CNF to know the structure
    num_vars, num_clauses, clauses = parse_cnf_file(cnf_file)
    if verbose:
        print(f"\nCNF: {num_vars} variables, {num_clauses} clauses")

    # Parse graph (optional, for verification)
    n_vertices, n_edges, edges = parse_graph_file(graph_file)
    if verbose:
        print(f"Graph: {n_vertices} vertices, {n_edges} edges")

    # Parse coloring
    coloring = parse_coloring_file(coloring_file)
    if verbose:
        print(f"Coloring: {len(coloring)} vertices colored")

    # Extract assignment
    assignment, mapping_info = extract_sat_assignment(coloring, num_vars, verbose, palette_override)

    if assignment is None:
        print("\nERROR: Failed to extract valid SAT assignment!")
        print("\nTroubleshooting suggestions:")
        print("1. Try --swap-tf option if all variables seem inverted")
        print("2. Check the debug output to see which vertex colors are assigned")
        print("3. Verify the CNF-to-3COL reduction used the standard vertex numbering")
        sys.exit(1)

    # Verify assignment
    satisfies, clause_results = verify_sat_assignment(clauses, assignment, verbose)

    # Save results to file
    save_assignment_file(assignment_file, assignment, mapping_info, satisfies, clause_results, cnf_file)

    # Output results
    print("\n" + "="*70)
    print("                           FINAL RESULT")
    print("="*70)

    if satisfies:
        print("✓ SUCCESS: The 3-coloring corresponds to a SATISFYING assignment!")
        print(f"✓ The assignment SATISFIES the CNF formula in {cnf_file}")

        # Output assignment summary
        true_vars = [v for v, val in assignment.items() if val]
        false_vars = [v for v, val in assignment.items() if not val]

        print(f"\nAssignment summary:")
        print(f"  TRUE variables ({len(true_vars)}): ", end="")
        if len(true_vars) <= 10:
            print(f"x{', x'.join(map(str, sorted(true_vars)))}")
        else:
            print(f"x{', x'.join(map(str, sorted(true_vars)[:5]))}, ... ({len(true_vars)} total)")

        print(f"  FALSE variables ({len(false_vars)}): ", end="")
        if len(false_vars) <= 10:
            print(f"x{', x'.join(map(str, sorted(false_vars)))}")
        else:
            print(f"x{', x'.join(map(str, sorted(false_vars)[:5]))}, ... ({len(false_vars)} total)")
    else:
        print("✗ FAILURE: The extracted assignment does NOT satisfy the CNF!")
        print(f"✗ The assignment DOES NOT SATISFY the CNF formula in {cnf_file}")

        # Show which clauses failed
        failed_clauses = [r['clause_num'] for r in clause_results if not r['satisfied']]
        print(f"\nUnsatisfied clauses: {failed_clauses}")

        # Show assignment summary for debugging
        true_vars = [v for v, val in assignment.items() if val]
        false_vars = [v for v, val in assignment.items() if not val]

        print(f"\nExtracted assignment summary:")
        print(f"  TRUE: {len(true_vars)} variables")
        print(f"  FALSE: {len(false_vars)} variables")

        if len(true_vars) == 1:
            print(f"  Only x{true_vars[0]} is TRUE - if this should be FALSE, try --swap-tf")
        elif len(false_vars) == 1:
            print(f"  Only x{false_vars[0]} is FALSE - if this should be TRUE, try --swap-tf")

        print("\nPossible issues:")
        print("  1. TRUE/FALSE vertices might be swapped - try --swap-tf option")
        print("  2. The graph doesn't correspond to the CNF file")
        print("  3. There's a bug in the CNF-to-3COL reduction")

    print("\n" + "="*70)
    print(f"Results saved to: {assignment_file}")
    print("="*70)

    # Exit with appropriate code
    sys.exit(0 if satisfies else 1)


if __name__ == "__main__":
    main()
