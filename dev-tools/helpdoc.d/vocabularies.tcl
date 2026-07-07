# ----------------------------------------------------------------------------
# vocabularies.tcl  --  closed vocabularies for the helpdoc -strict pass
#
# Sourced by helpdoc.d/validate.tcl; its location may be overridden with the
# HELPDOC_VOCABULARIES environment variable. Defines the ::helpdoc::vocab
# array with the keys: dim_base, dimensionless, units_base, kind, computed.
# ----------------------------------------------------------------------------

namespace eval ::helpdoc {

    variable vocab

    # base dimensions: a <dimensionality> value is a product-of-powers of
    # these tokens (e.g. "energy length^-1" for force).
    set vocab(dim_base) {
        energy length time mass charge temperature
    }

    # the literal empty-product token: <dimensionality> for which no <units>
    # field is required.
    set vocab(dimensionless) {
        dimensionless
    }

    # base unit tokens: a <units> value is a product-of-powers of these.
    # Only tokens actually bound to a variable in the current .def files are
    # listed; units appearing only in prose are excluded.
    set vocab(units_base) {
        Ry Hartree eV
        bohr angstrom
        electron_mass amu
        e
        K
        kbar GPa
        THz
        cm-1
        s
        degrees
        states
        kcal mol
        alat tpiba
    }

    # allowed values of the optional "kind" attribute on a <default> element
    set vocab(kind) {
        literal ref expr computed conditional
    }

    # allowed sentinel tokens for the body of a <default kind="computed">
    set vocab(computed) {
        from_pseudopotential from_xml from_environment internal
    }
}
