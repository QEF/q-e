# Version update checklist

Detailed companion to item 4 ("update the release number in all the
documentation that contains references to it") and item 9 ("update the
version number in `include/qe_version.h`") of `release-checklist.md`.
This enumerates every place a QE version string is *set* — either the
single source of truth, a place that reads it automatically, or a place
that hardcodes its own copy and must be edited by hand.

Run **`dev-tools/check_versions.sh`** before tagging a release: it reads
`include/qe_version.h` and reports every location in sections listed
below that currently disagrees with it, without changing anything. Use
it to see, at a glance, which of these checkboxes still need attention.

## 1. Single source of truth

- [ ] **`include/qe_version.h`** — `version_number = '7.6-dev'` (Fortran
  `CHARACTER` constant). This is the master value. Everything in
  section 2 below reads it automatically at build/doc-generation time;
  everything in section 3 does **not** and must be edited separately.

## 2. Derived automatically from `qe_version.h` — no manual edit needed

These all `awk`/`#include` the value out of `qe_version.h`, so once step 1
is done they stay in sync by themselves. Listed so you can sanity-check
them after a build, not because they need hand-editing:

- `MODULE global_version` (`Modules/version.f90`) — `#include`s
  `qe_version.h` directly; every executable's `USE global_version, ONLY:
  version_number` gets the value from here.
- Runtime "Program `<code>` v.X.Y starts" banner —
  `Modules/environment.f90:74` (`code_version = TRIM(code)//" v."//TRIM(version_number)`),
  printed by `opening_message`. Drives the header for **every** QE
  executable (pw.x, ph.x, cp.x, pp.x, etc.) in one place.
- XML/UPF/output metadata that embeds the version: `Modules/qexsd.f90`
  (XML schema `version=` attribute), `upflib/write_upf_new.f90` (UPF
  `VERSION` attribute), `PP/src/write_hamiltonians.f90`,
  `atomic/src/export_upf.f90`, `TDDFPT/src/lr_magnons_main.f90`,
  `PHonon/PH/ph_restart.f90`, `external/qe-gipaw/src/output_magres.f90`.
- Doc-generation `Makefile`s that compute `VERSION := $(shell awk ... qe_version.h)`
  and feed it to `dev-tools/helpdoc -version $(VERSION)` or pandoc
  `--metadata title`: `HP/Doc/Makefile`, `KCW/Doc/Makefile`,
  `QEHeat/Doc/Makefile`, `TDDFPT/Doc/Makefile`, `atomic/Doc/Makefile`,
  `PWCOND/Doc/Makefile`, `CPV/Doc/Makefile`, `PW/Doc/Makefile`,
  `PP/Doc/Makefile`, `NEB/Doc/Makefile`, `PHonon/Doc/Makefile`.
- `GUI/QE-modes/Doc/Makefile` and `GUI/PWgui/Makefile` — generate a
  `version.tex` / `VERSION` file from `qe_version.h` (this is the pattern
  the hand-edited `.tex` files below should ideally be switched to).
- `git-rev.h` (`build*/git-rev.h`) — generated at build time from `git
  describe`; not release-checklist material.

## 3. Hardcoded independently — MUST be edited by hand every release

These do **not** read `qe_version.h` and were found already out of sync
with each other (mix of `7.3.1`, `7.4`, `7.5`) and with the current
`7.6-dev`, confirming they get missed. `dev-tools/check_versions.sh`
checks all of them against `include/qe_version.h` on every run.

- [ ] **`CMakeLists.txt:17`** — `project(qe VERSION 7.5)`. Feeds
  `PACKAGE_VERSION` used in `qeConfigVersion.cmake`. *Possible fix, not
  yet implemented — needs a decision first:* CMake can extract the
  version at configure time with `file(STRINGS ...)`/regex before the
  `project()` call, e.g. parse `include/qe_version.h` and pass the
  numeric part in as `project(qe VERSION ${QE_VERSION_NUM} ...)`. Caveat:
  CMake's `VERSION` field is numeric-only (up to 4 dot-separated
  components), so a `-dev` suffix would need to be stripped before
  passing it in — worth confirming that's acceptable before wiring it up,
  since it touches the build system for every user.
- [ ] **`install/configure.ac:9`** — `AC_INIT(ESPRESSO, 7.5, , espresso)`.
  *Possible fix, not yet implemented — same caveat as above:* autoconf
  supports pulling the version in via `m4_esyscmd_s([...])` run against
  `include/qe_version.h` at `autoreconf` time, so `configure.ac` itself
  never hardcodes a version. Would need re-running `autoreconf` and a
  build-system sanity check before adopting.
- [ ] **`configure`** (generated, `PACKAGE_VERSION='7.5'`) — do not edit
  directly, regenerate with `autoconf`/`autoreconf` from `configure.ac`
  per release-checklist item 8, then verify the version string came out
  right.
- [ ] **`PW/Doc/user_guide.tex:2`** — `\def\version{7.5}`
- [ ] **`PP/Doc/user_guide.tex:2`** — `\def\version{7.4}`
- [ ] **`PHonon/Doc/user_guide.tex:2`** — `\def\version{7.4}`
- [ ] **`NEB/Doc/user_guide.tex:2`** — `\def\version{7.4}`
- [ ] **`CPV/Doc/tp.tex:2`** — `\def\version{7.4}`
- [ ] **`Doc/user_guide.tex:2`** — `\def\version{7.5.0}`
- [ ] **`Doc/brillouin_zones.tex:2`** — `\def\version{7.3.1}`
- [ ] **`Doc/Hubbard_input.tex:3`** — `\def\version{7.3.1}`
- [ ] **`Doc/release-notes`** — the topmost section header carries the
  in-progress version label (currently `7.6`, matching `qe_version.h`'s
  `7.6-dev` base) and should be updated whenever a new version starts.
  `dev-tools/check_versions.sh` checks that this label matches
  `qe_version.h`, but only the label — it cannot check whether the
  content underneath is actually complete (release-checklist item 6);
  that still needs a human read-through.

## 4. Stale/independent version literals — verify, don't blindly bump

Not part of the normal release bump, but flagged because they store a
version number disconnected from `qe_version.h` and could rot silently:

- **`EPW/src/epw.f90:88`** — `version_number = '6.1'` overwrites the
  global `version_number` right before `environment_start`/clocks are
  initialized, so `epw.x`'s startup banner and any header derived from
  it always print `EPW v.6.1` regardless of the actual QE release. This
  is **intentional**: EPW maintains its own version number, independent
  of the QE suite version, and any change to it is supervised by the EPW
  group. Do not "fix" this as part of a QE release bump — leave it out of
  scope for this checklist.
- **`PW/src/environ_pw_module.f90:95`** and
  **`CPV/src/environ_cp_module.f90:96`** — `IF (version_number == '6.3')`
  compatibility check against the external Environ library's expected QE
  API version, not the running QE version. Not a docs/header string, but
  worth a sanity check if Environ's supported version range has moved on.

## Suggested follow-up (optional, not required for a release)

Only `GUI/QE-modes/Doc/Makefile` generates its `\version` from
`qe_version.h` instead of hardcoding it. Auto-generating `\version` the
same way in the other `.tex` files (section 3) was considered and
**explicitly declined** by the maintainer: several of those documents may
be genuinely out of date in content, and silently auto-updating just the
version number could misleadingly imply the whole document is current.
The team wants to keep full manual responsibility for bumping these, so
that a version bump prompts an actual content review rather than being a
no-op. Do not re-propose this without new instruction — use
`dev-tools/check_versions.sh` for visibility instead.

The `CMakeLists.txt` / `install/configure.ac` auto-linking described in
section 3 is a separate, still-open question (build-system files, not
docs) — raise it if/when a release bump is being prepared.
