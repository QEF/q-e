proc ::pwscf::ppReadFilter {moduleObj channel} {
    return [readFilter::default $moduleObj $channel logical]
}
