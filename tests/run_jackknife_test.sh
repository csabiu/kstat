#!/bin/bash
#------------------------------------------------------------
# Jackknife Test Script
#
# This script:
# 1. Compiles the Fortran test program
# 2. Generates test data using the Fortran program
# 3. Runs the shell script jk.sh on the same data
# 4. Compares the results
#
# Usage: ./run_jackknife_test.sh [Ngal] [Nran] [Ndiv] [Ndim] [seed]
#        Default: 100 200 3 3 12345
#------------------------------------------------------------

set -e

# Parameters with defaults
NGAL=${1:-100}
NRAN=${2:-200}
NDIV=${3:-3}
NDIM=${4:-3}
SEED=${5:-12345}

# Directories
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SRC_DIR="$SCRIPT_DIR/../src"
BIN_DIR="$SCRIPT_DIR/../bin"
TEST_DIR="$SCRIPT_DIR"

echo "================================================"
echo "Jackknife Module Comparison Test"
echo "================================================"
echo "Parameters:"
echo "  Ngal = $NGAL"
echo "  Nran = $NRAN"
echo "  Ndiv = $NDIV"
echo "  Ndim = $NDIM"
echo "  Seed = $SEED"
echo "  Expected regions = $((NDIV**NDIM))"
echo "================================================"

# Change to test directory
cd "$TEST_DIR"

# Step 1: Compile the Fortran test program
echo ""
echo "Step 1: Compiling Fortran test program..."
gfortran -O2 -o test_jackknife "$SRC_DIR/jackknife.F90" "$TEST_DIR/test_jackknife.F90"
echo "  Compiled successfully."

# Step 2: Run the Fortran test program
echo ""
echo "Step 2: Running Fortran test program..."
./test_jackknife $NGAL $NRAN $NDIV $NDIM $SEED

# Step 3: Run the shell script
echo ""
echo "Step 3: Running shell script jk.sh..."
bash "$BIN_DIR/jk.sh" test_gal.dat test_ran.dat $NDIV $NDIM

# Step 4: Compare results
echo ""
echo "Step 4: Comparing results..."
echo ""

# Function to compare two files and count differences
compare_files() {
    local fortran_file="$1"
    local shell_file="$2"
    local label="$3"

    if [ ! -f "$fortran_file" ]; then
        echo "ERROR: Fortran output file '$fortran_file' not found!"
        return 1
    fi
    if [ ! -f "$shell_file" ]; then
        echo "ERROR: Shell output file '$shell_file' not found!"
        return 1
    fi

    # Extract just the region IDs for comparison
    # Fortran format: x y [z] [weight] region_id
    # Shell format: x y [z] [weight] region_id

    local fortran_regions=$(mktemp)
    local shell_regions=$(mktemp)

    if [ "$NDIM" == "2" ]; then
        awk '{print $NF}' "$fortran_file" > "$fortran_regions"
        awk '{print $NF}' "$shell_file" > "$shell_regions"
    else
        awk '{print $NF}' "$fortran_file" > "$fortran_regions"
        awk '{print $NF}' "$shell_file" > "$shell_regions"
    fi

    # Count lines
    local fortran_lines=$(wc -l < "$fortran_regions")
    local shell_lines=$(wc -l < "$shell_regions")

    echo "  $label:"
    echo "    Fortran output: $fortran_lines lines"
    echo "    Shell output:   $shell_lines lines"

    if [ "$fortran_lines" != "$shell_lines" ]; then
        echo "    WARNING: Line count mismatch!"
        rm -f "$fortran_regions" "$shell_regions"
        return 1
    fi

    # The shell script sorts data differently, so direct comparison won't work.
    # Instead, we compare the distribution of points per region.

    local fortran_dist=$(mktemp)
    local shell_dist=$(mktemp)

    sort "$fortran_regions" | uniq -c | sort -k2 -n > "$fortran_dist"
    sort "$shell_regions" | uniq -c | sort -k2 -n > "$shell_dist"

    echo ""
    echo "    Region distribution comparison:"
    echo "    Fortran              | Shell Script"
    echo "    Count  Region        | Count  Region"
    echo "    ----   ------        | ----   ------"

    paste "$fortran_dist" "$shell_dist" | while read line; do
        echo "    $line"
    done

    # Check if distributions match
    if diff -q "$fortran_dist" "$shell_dist" > /dev/null 2>&1; then
        echo ""
        echo "    RESULT: Distributions MATCH!"
        rm -f "$fortran_regions" "$shell_regions" "$fortran_dist" "$shell_dist"
        return 0
    else
        echo ""
        echo "    RESULT: Distributions DIFFER (see above)"
        rm -f "$fortran_regions" "$shell_regions" "$fortran_dist" "$shell_dist"
        return 1
    fi
}

# Compare galaxy files
GAL_RESULT=0
compare_files "test_fortran_gal.jk" "test_gal.dat.jk" "Galaxy comparison" || GAL_RESULT=1

echo ""

# Compare random files
RAN_RESULT=0
compare_files "test_fortran_ran.jk" "test_ran.dat.jk" "Random comparison" || RAN_RESULT=1

echo ""
echo "================================================"
echo "Test Summary"
echo "================================================"

if [ "$GAL_RESULT" == "0" ] && [ "$RAN_RESULT" == "0" ]; then
    echo "ALL TESTS PASSED!"
    echo ""
    echo "The Fortran jackknife implementation produces the same"
    echo "region distribution as the jk.sh shell script."
    EXIT_CODE=0
else
    echo "SOME TESTS FAILED!"
    echo ""
    echo "Please review the differences above."
    EXIT_CODE=1
fi

echo ""
echo "Test files preserved for manual inspection:"
echo "  test_gal.dat          - Input galaxy coordinates"
echo "  test_ran.dat          - Input random coordinates"
echo "  test_fortran_gal.jk   - Fortran JK output (galaxies)"
echo "  test_fortran_ran.jk   - Fortran JK output (randoms)"
echo "  test_gal.dat.jk       - Shell JK output (galaxies)"
echo "  test_ran.dat.jk       - Shell JK output (randoms)"
echo "================================================"

exit $EXIT_CODE
