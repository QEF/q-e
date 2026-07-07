#!/bin/sh
#
# run-tests.sh -- test harness for the helpdoc "--strict" validation mode
#                 and for parsing/rendering of the modernized .def
#                 constructs (dimensionality, units, default kind, case).
#
# It checks that:
#   * a fully-migrated .def file (good.def) passes --strict;
#   * good.def also builds in non-strict (legacy) mode;
#   * each bad_*.def file is REJECTED by --strict;
#   * each bad_*.def file still BUILDS in non-strict mode (lenient).
#
# Exit status is 0 iff every expectation holds.
#

set -u

testdir=$(cd "$(dirname "$0")" && pwd)
devtools=$(cd "$testdir/.." && pwd)
helpdoc="$devtools/helpdoc"

# use the test-local vocabulary stub
HELPDOC_VOCABULARIES="$testdir/vocabularies.tcl"
export HELPDOC_VOCABULARIES

workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT
cp "$devtools/input_xx.xsl" "$workdir/"

fail=0
pass=0

run_helpdoc() {
    # run_helpdoc <strict|lenient> <deffile>
    mode=$1
    def=$2
    cp "$testdir/$def" "$workdir/"
    if [ "$mode" = strict ]; then
        ( cd "$workdir" && "$helpdoc" -strict "$def" ) > "$workdir/log" 2>&1
    else
        ( cd "$workdir" && "$helpdoc" "$def" ) > "$workdir/log" 2>&1
    fi
    return $?
}

expect_pass() {
    desc=$1; mode=$2; def=$3
    if run_helpdoc "$mode" "$def"; then
        echo "PASS: $desc"
        pass=$((pass + 1))
    else
        echo "FAIL: $desc (helpdoc returned non-zero, expected success)"
        sed 's/^/    | /' "$workdir/log"
        fail=$((fail + 1))
    fi
}

expect_fail() {
    desc=$1; mode=$2; def=$3
    if run_helpdoc "$mode" "$def"; then
        echo "FAIL: $desc (helpdoc returned success, expected rejection)"
        fail=$((fail + 1))
    else
        echo "PASS: $desc"
        pass=$((pass + 1))
    fi
}

echo "=== helpdoc --strict test harness ==="

# the migrated file must pass strict, and also build leniently
expect_pass "good.def passes --strict"            strict  good.def
expect_pass "good.def builds in non-strict mode"  lenient good.def

# a childless bare variable re-reference (e.g. "var nwf" with no body, used
# in an alternate card syntax branch) must NOT trip the mandatory
# dimensionality/units checks under --strict
expect_pass "good_bare_ref.def passes --strict (childless bare re-reference)" \
            strict  good_bare_ref.def
expect_pass "good_bare_ref.def builds in non-strict mode" \
            lenient good_bare_ref.def

# each malformed file must be rejected by --strict ...
for def in bad_missing_dim.def bad_missing_units.def bad_vocab.def \
           bad_index_range.def bad_default_mixed.def bad_default_nofallback.def \
           bad_units_base.def bad_units_exponent.def bad_enum_default.def \
           bad_dim_base.def bad_units_nofallback.def \
           bad_dim_on_nonreal.def bad_required_default.def
do
    expect_fail "$def rejected by --strict"        strict  "$def"
    # ... but must still build leniently (migration-friendliness)
    expect_pass "$def still builds in non-strict mode" lenient "$def"
done

echo "=== $pass passed, $fail failed ==="
[ "$fail" -eq 0 ]
