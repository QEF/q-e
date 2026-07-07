#
# vocabularies.tcl -- TEST stub of the helpdoc closed vocabularies, loaded by
# the --strict test harness via the HELPDOC_VOCABULARIES env var. Populates
# ::helpdoc::vocab with dim_base, units_base, kind, computed, dimensionless.
#

namespace eval ::helpdoc {
    variable vocab

    set vocab(dim_base) {
        energy length time mass charge temperature
    }

    set vocab(dimensionless) {
        dimensionless
    }

    set vocab(units_base) {
        bohr angstrom Ry eV Hartree bohr_magneton atomic_mass
        electron_mass e states alat
    }

    set vocab(kind) {
        literal ref expr computed conditional
    }

    set vocab(computed) {
        from_pseudo from_input automatic
    }
}
