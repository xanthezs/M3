"""
CNF Evaluator
Verifies if a CNF assignment is correct

Loads CNF in DIMACS format and variable assignments in either
readable format (variable: value) or DIMACS solution format.

Usage: python evaluate_cnf.py cnf_file assignment_file
Where:
  cnf_file - Path to DIMACS .cnf file
  assignment_file - Path to variable assignment file (formatted with variable: value)

Outputs evaluation results and detailed analysis of any failed clauses.
"""


import itertools
import os

def load_cnf_dimacs(filename):
    """Load CNF clauses from a DIMACS .cnf file."""
    if not os.path.exists(filename):
        raise FileNotFoundError(f"CNF file not found: {filename}")

    clauses = []
    num_vars = 0
    num_clauses = 0

    with open(filename, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            if line.startswith('c'):
                continue  # Skip comments

            if line.startswith('p cnf'):
                # Parse problem line: p cnf <num_vars> <num_clauses>
                parts = line.split()
                if len(parts) >= 4:
                    num_vars = int(parts[2])
                    num_clauses = int(parts[3])
                    print(f"📋 DIMACS header: {num_vars} variables, {num_clauses} clauses expected")
                continue

            # Parse clause
            try:
                literals = [int(lit) for lit in line.split() if lit != '0']
                if literals:  # Only add non-empty clauses
                    clauses.append(literals)
            except ValueError as e:
                print(f"⚠️  Warning: Malformed clause at line {line_num}: {line}")
                continue

    print(f"📁 Loaded {len(clauses)} clauses from {filename}")
    if num_clauses > 0 and len(clauses) != num_clauses:
        print(f"⚠️  Warning: Expected {num_clauses} clauses, but loaded {len(clauses)}")

    return clauses, num_vars

def load_assignment_readable_format(filename):
    """
    Load variable assignment from readable format file.
    Handles formats:
    - 'variable: value' (colon-separated, from exact_sha256_assignments.txt)
    - 'variable value' (space-separated)
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Assignment file not found: {filename}")

    assignment = {}
    with open(filename, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            try:
                # Handle colon format: "variable: value"
                if ':' in line:
                    parts = line.split(':')
                    if len(parts) >= 2:
                        var = int(parts[0].strip())
                        val = int(parts[1].strip())
                # Handle space format: "variable value"
                else:
                    parts = line.split()
                    if len(parts) >= 2:
                        var = int(parts[0])
                        val = int(parts[1])
                    else:
                        continue

                if val not in [0, 1]:
                    print(f"⚠️  Warning: Invalid value {val} for variable {var} at line {line_num}")
                    continue
                assignment[var] = val

            except ValueError as e:
                print(f"⚠️  Warning: Could not parse line {line_num}: {line}")
                continue

    print(f"📁 Loaded assignments for {len(assignment)} variables from {filename}")
    return assignment

def load_assignment_dimacs_solution(filename):
    """
    Load variable assignment from DIMACS solution format.
    Format: sequence of literals ending with 0, where positive = true, negative = false
    Example: "1 -2 3 -4 5 0" means var1=true, var2=false, var3=true, var4=false, var5=true
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Solution file not found: {filename}")

    assignment = {}

    with open(filename, 'r') as f:
        content = f.read().strip()

        # Parse all literals from the file
        literals = []
        for token in content.split():
            try:
                lit = int(token)
                if lit == 0:
                    break  # End of solution
                literals.append(lit)
            except ValueError:
                continue

        # Convert literals to assignments
        for lit in literals:
            var = abs(lit)
            val = 1 if lit > 0 else 0
            assignment[var] = val

    print(f"📁 Loaded DIMACS solution with {len(assignment)} variable assignments from {filename}")
    return assignment

def load_assignment(filename):
    """
    Auto-detect format and load variable assignment.
    Tries DIMACS solution format first, then readable format.
    """
    try:
        # First, try to detect if it's a DIMACS solution format
        with open(filename, 'r') as f:
            first_line = f.readline().strip()

            # If first non-comment line contains only numbers and possibly 0, it's likely DIMACS solution
            if first_line and not first_line.startswith('#'):
                tokens = first_line.split()
                if all(token.lstrip('-').isdigit() for token in tokens):
                    print(f"🔍 Detected DIMACS solution format")
                    return load_assignment_dimacs_solution(filename)

        # Otherwise, try readable format
        print(f"🔍 Detected readable assignment format")
        return load_assignment_readable_format(filename)

    except Exception as e:
        print(f"❌ Error loading assignment: {e}")
        raise

def evaluate_clause(clause, assignment):
    """Return True if the clause is satisfied under assignment."""
    satisfied_literals = []
    unsatisfied_literals = []
    missing_vars = []

    for lit in clause:
        var = abs(lit)
        if var not in assignment:
            missing_vars.append(var)
            continue

        val = assignment[var]
        # Positive literal: satisfied if variable is true
        # Negative literal: satisfied if variable is false
        if (lit > 0 and val == 1) or (lit < 0 and val == 0):
            satisfied_literals.append(lit)
        else:
            unsatisfied_literals.append(lit)

    if missing_vars:
        return False, f"Missing variables: {missing_vars}"

    if satisfied_literals:
        return True, f"Satisfied by literal(s): {satisfied_literals}"
    else:
        return False, f"No literals satisfied. Unsatisfied: {unsatisfied_literals}"

def evaluate_cnf(clauses, assignment, verbose=False, max_failures=10):
    """Evaluate CNF formula under assignment."""
    failed_clauses = []
    satisfied_count = 0

    for i, clause in enumerate(clauses):
        satisfied, reason = evaluate_clause(clause, assignment)
        if satisfied:
            satisfied_count += 1
        else:
            failed_clauses.append((i, clause, reason))
            if verbose and len(failed_clauses) <= max_failures:
                print(f"❌ Clause #{i}: {clause} - {reason}")

    success_rate = (satisfied_count / len(clauses)) * 100 if clauses else 0
    print(f"📊 Evaluation results: {satisfied_count}/{len(clauses)} clauses satisfied ({success_rate:.1f}%)")

    if failed_clauses:
        print(f"❌ {len(failed_clauses)} clauses failed out of {len(clauses)} total")
        if len(failed_clauses) > max_failures:
            print(f"   (showing first {max_failures} failures)")
        return False, failed_clauses

    return True, []

def analyze_assignment_coverage(clauses, assignment):
    """Analyze how well the assignment covers the CNF variables."""
    all_vars = set(abs(lit) for clause in clauses for lit in clause)
    assigned_vars = set(assignment.keys())
    missing_vars = all_vars - assigned_vars
    extra_vars = assigned_vars - all_vars

    print(f"\n📊 Assignment Coverage Analysis:")
    print(f"   Variables in CNF: {len(all_vars)}")
    print(f"   Variables assigned: {len(assigned_vars)}")
    print(f"   Missing assignments: {len(missing_vars)}")
    print(f"   Extra assignments: {len(extra_vars)}")

    if all_vars:
        var_range = f"{min(all_vars)}-{max(all_vars)}"
        print(f"   CNF variable range: {var_range}")

    if assigned_vars:
        assigned_range = f"{min(assigned_vars)}-{max(assigned_vars)}"
        print(f"   Assignment range: {assigned_range}")

    if missing_vars:
        print(f"   First 10 missing: {sorted(missing_vars)[:10]}")
    if extra_vars:
        print(f"   First 10 extra: {sorted(extra_vars)[:10]}")

    return len(missing_vars) == 0

def analyze_failed_clauses(failed_clauses, assignment, max_analysis=5):
    """Provide detailed analysis of failed clauses."""
    if not failed_clauses:
        return

    print(f"\n🔍 Detailed Analysis of Failed Clauses:")

    for i, (clause_idx, clause, reason) in enumerate(failed_clauses[:max_analysis]):
        print(f"\n❌ Failed Clause #{clause_idx}: {clause}")
        print(f"   Reason: {reason}")

        # Show variable assignments for this clause
        print(f"   Variable assignments:")
        for lit in clause:
            var = abs(lit)
            if var in assignment:
                val = assignment[var]
                expected = "1" if lit > 0 else "0"
                actual = str(val)
                status = "✅" if actual == expected else "❌"
                print(f"     var{var} = {val} (literal {lit} expects {expected}) {status}")
            else:
                print(f"     var{var} = MISSING (literal {lit})")

def main():
    # File paths - update these to match your files
    cnf_file = "thunder.cnf"  # Your CNF file
    assignment_file = "thunder_msb.txt"  # Readable format
    # assignment_file = "exact_sha256_solution.txt"  # DIMACS solution format

    try:
        # Load files
        print("🔍 Loading CNF file...")
        clauses, expected_vars = load_cnf_dimacs(cnf_file)

        print("🔍 Loading assignment file...")
        assignment = load_assignment(assignment_file)

        # Analyze coverage
        complete_assignment = analyze_assignment_coverage(clauses, assignment)

        if not complete_assignment:
            print("\n⚠️  Warning: Assignment is incomplete. Some clauses may fail due to missing variables.")

        # Evaluate CNF
        print(f"\n🧪 Evaluating {len(clauses)} clauses...")
        satisfied, failed_clauses = evaluate_cnf(clauses, assignment, verbose=True, max_failures=5)

        # Analyze failures if any
        if failed_clauses:
            analyze_failed_clauses(failed_clauses, assignment, max_analysis=3)

        # Final result
        print(f"\n{'='*50}")
        if satisfied:
            print("✅ SUCCESS: CNF is satisfied by the assignment!")
            print("   Your CNF encoding and variable assignment are working correctly.")
            print("   This means:")
            print("   - CNF generation logic is correct")
            print("   - Assignment generation logic is correct")
            print("   - Variable numbering is consistent")
        else:
            print("❌ FAILURE: CNF is NOT satisfied by the assignment.")
            print(f"   {len(failed_clauses)} clause(s) failed.")
            print("   This suggests:")
            print("   - Bug in CNF generation logic, OR")
            print("   - Bug in assignment generation logic, OR")
            print("   - Mismatch between CNF and assignment variable numbering")

            # Provide debugging suggestions
            print(f"\n🔧 Debugging suggestions:")
            if complete_assignment:
                print("   - Assignment covers all variables, so likely a logic error")
                print("   - Check the first few failed clauses for patterns")
                print("   - Verify SHA-256 computation in assignment generator")
            else:
                print("   - Assignment is incomplete, fix missing variables first")
                print("   - Ensure CNF and assignment use same variable numbering")

    except FileNotFoundError as e:
        print(f"❌ File error: {e}")
        print(f"   Make sure these files exist:")
        print(f"   - CNF file: {cnf_file}")
        print(f"   - Assignment file: {assignment_file}")
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
