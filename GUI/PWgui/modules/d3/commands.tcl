# ------------------------------------------------------------------------
#
# ------------------------------------------------------------------------
proc ::pwscf::d3ReadFilter {moduleObj channel} {
    return [readFilter::default $moduleObj $channel {logical amass}]
}