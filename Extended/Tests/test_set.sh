#!/bin/sh
# Tests of the Extended model (see Tests/README.md). Called by 'make test'
# (Makefile: test: sh Tests/test_set.sh).
# The sources are not modified: everything is compiled in temporary directories.
#   1. Unit tests (Fortran): numerical_utils and the physical functions of global.f90
#   2. Regression of the terminal output of the BD+43 case
#   3. Physical checks (energy budgets, prescriptions, One-Zone, OpenMP, ...)
# Extra arguments are passed to Tests/test_physics.py (e.g. --only BD43 HVS).

# Works from any directory: moves to Extended/ (the parent of Tests/)
cd "$(dirname "$0")/.." || exit 1

failed=""

echo "=== 1. Unit tests ==="
unit_dir=$(mktemp -d "${TMPDIR:-/tmp}/code-bs-unit-XXXXXX")
cp constants.f90 common_parameters.f90 numerical_utils.f90 global.f90 Tests/test_numerical_utilities.f90 \
   Tests/test_global_functions.f90 "$unit_dir"/
cp Tests/cases/BD43/system_parameters.f90 "$unit_dir"/
FFLAGS="-fopenmp -ffree-line-length-none -fbackslash -O2"
if (cd "$unit_dir" && \
    gfortran $FFLAGS -c constants.f90 numerical_utils.f90 system_parameters.f90 common_parameters.f90 global.f90 && \
    gfortran $FFLAGS test_numerical_utilities.f90 constants.o numerical_utils.o -o test_numerical_utilities && \
    gfortran $FFLAGS test_global_functions.f90 constants.o numerical_utils.o system_parameters.o common_parameters.o global.o \
             -o test_global_functions && \
    ./test_numerical_utilities && ./test_global_functions); then
    :
else
    failed="$failed unit"
fi
rm -rf "$unit_dir"

echo
echo "=== 2. Terminal-output regression (BD+43) ==="
python3 Tests/test_terminal_output.py || failed="$failed regression"

echo
echo "=== 3. Physical checks ==="
python3 Tests/test_physics.py "$@" || failed="$failed physics"

echo
if [ -n "$failed" ]; then
    echo "FAILED:$failed"
    exit 1
fi
echo "All the tests passed."
