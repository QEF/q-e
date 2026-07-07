# helpdoc --strict test harness

Tests for the modernized `.def` toolchain: the `dimensionality` / `units`
elements, the restructured `default` (with `kind` attribute and `case`
children), and the `helpdoc --strict` validation pass.

## Running

    ./run-tests.sh

Exit status is `0` iff every expectation holds.

## What it checks

* `good.def` -- a fully-migrated input description -- passes `--strict`
  and also builds in non-strict (legacy) mode.
* each `bad_*.def` file is rejected by `--strict`, one per rule:
  * `bad_missing_dim.def`        -- REAL var without `<dimensionality>`
  * `bad_dim_on_nonreal.def`     -- non-REAL var carrying a `<dimensionality>`
  * `bad_required_default.def`   -- `status {required}` var that also has a `<default>`
  * `bad_missing_units.def`      -- dimensionful var without `<units>`
  * `bad_vocab.def`              -- dimensionality term outside the vocabulary
  * `bad_index_range.def`        -- keyed-clist index out of `-start/-end` bounds
  * `bad_default_mixed.def`      -- `<default>` malformed (text body + cases)
  * `bad_default_nofallback.def` -- conditional `<default>` without a fallback `case`
  * `bad_units_base.def`         -- `<units>` base token outside the vocabulary
  * `bad_units_exponent.def`     -- `<units>` malformed exponent
* every `bad_*.def` file still *builds* in non-strict mode.

## Files

* `vocabularies.tcl` -- a test stub of the closed vocabularies, pointed to
  via the `HELPDOC_VOCABULARIES` environment variable.
