# Tcl package index file, version 1.0

set channel [open [file join $dir VERSION] r]
set version [gets $channel]
close $channel
package ifneeded Guib $version [list source [file join $dir init.tcl]]
