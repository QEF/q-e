# Tcl package index file, version 1.0

set version [gets [open [file join $dir VERSION] r]]
package ifneeded Guib $version [list source [file join $dir init.tcl]]
