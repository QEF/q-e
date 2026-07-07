#!/bin/bash
#
# Report every hardcoded QE version string next to include/qe_version.h,
# the single source of truth, so drift is visible before a release.
# See dev-tools/version-update-checklist.md for the full rationale and
# which locations are supposed to be manual vs automatic.
#
# NOT checked here: EPW/src/epw.f90. Its version_number is maintained
# independently by the EPW group and is intentionally decoupled from the
# QE suite version -- it is out of scope by design, not an oversight.
#
# Usage: dev-tools/check_versions.sh
# Exit status: 0 if everything matches include/qe_version.h, 1 otherwise.
#
# Besides the human-readable report on stdout, this writes a small
# machine-readable status file (default: .version-check-status, override
# with $VERSION_CHECK_STATUS_FILE) with one PASS/FAIL line per section --
# BUILD_STATUS, DOCS_STATUS, NOTES_STATUS -- so CI can run this script
# once and have separate jobs report on each section without re-running
# the checks.

set -euo pipefail
cd "$(dirname "$0")/.."

REF_RAW=$(grep -oP "version_number\s*=\s*'\K[^']+" include/qe_version.h)
# Strip any pre-release marker (-dev, -rc1, _rc2, rc2, -beta, whatever
# literal happens to be in use -- with or without a '-'/'_' separator)
# to get the bare X.Y[.Z] base.
REF_BASE=$(sed -E 's/[-_]?[a-zA-Z].*//' <<< "$REF_RAW")

status=0
status_build=0
status_docs=0
status_notes=0
section=

# Compares only the bare X.Y[.Z] base -- any pre-release marker (on either
# side) is ignored, since most of these locations can't carry one at all
# (CMake's project(VERSION ...) is numeric-only) or shouldn't by
# convention (docs/release-notes describe a released version, never an
# RC candidate).
check() {
  local label=$1 value=$2
  local value_base
  value_base=$(sed -E 's/[-_]?[a-zA-Z].*//' <<< "$value")
  if [ "$value_base" = "$REF_BASE" ]; then
    printf '  [OK]      %-45s %s\n' "$label" "$value"
  else
    printf '  [DIFFERS] %-45s %-12s (expected %s)\n' "$label" "$value" "$REF_BASE"
    status=1
    case "$section" in
      build) status_build=1 ;;
      docs)  status_docs=1 ;;
      notes) status_notes=1 ;;
    esac
  fi
}

echo "Reference (include/qe_version.h): $REF_RAW   [base: $REF_BASE]   (pre-release markers are ignored in all comparisons)"
echo

section=build
echo "Build system (must match to keep package metadata correct):"
cmake_v=$(grep -A3 '^project(qe' CMakeLists.txt | grep -oP 'VERSION\s+\K[0-9.]+' | head -1 || echo '?')
check "CMakeLists.txt (project VERSION)" "$cmake_v"

ac_v=$(grep -oP 'AC_INIT\(ESPRESSO,\s*\K[^,\s]+' install/configure.ac || echo '?')
check "install/configure.ac (AC_INIT)" "$ac_v"

if [ -f configure ]; then
  cfg_v=$(grep -oP "PACKAGE_VERSION='\K[^']+" configure || echo '?')
  check "configure (generated -- regenerate, don't hand-edit)" "$cfg_v"
fi

section=docs
echo
echo "Documentation (\\def\\version{...}, intentionally hand-maintained --"
echo "reported for visibility before a release, not auto-fixed):"
while IFS=: read -r file _line content; do
  file=${file#./}
  value=$(sed -E 's/.*\\def\\version\{([^}]*)\}.*/\1/' <<< "$content")
  check "$file" "$value"
done < <(grep -rn '\\def\\version{' --include='*.tex' . 2>/dev/null \
           | grep -vE '/(test-suite|external|build[^/]*|GUI/QE-modes)/')

section=notes
echo
echo "Release notes (checks only that the top section's version label matches"
echo "include/qe_version.h -- does NOT check that its content is complete):"
notes_v=$(grep -oP '[0-9]+\.[0-9]+(\.[0-9]+)?' Doc/release-notes | head -1 || echo '?')
check "Doc/release-notes (top entry label)" "$notes_v"
echo
section=

STATUS_FILE=${VERSION_CHECK_STATUS_FILE:-.version-check-status}
{
  echo "BUILD_STATUS=$status_build"
  echo "DOCS_STATUS=$status_docs"
  echo "NOTES_STATUS=$status_notes"
} > "$STATUS_FILE"

if [ "$status" -eq 0 ]; then
  echo "All checked locations match include/qe_version.h."
else
  echo "Some locations differ from include/qe_version.h -- see dev-tools/version-update-checklist.md."
fi

exit "$status"
