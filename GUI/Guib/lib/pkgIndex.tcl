
#****h* /TclTkLib
#
# DESCRIPTION

# This is Tone Kokalj's Tcl/Tk utility library. It consits of two
# parts: tcl-part (file: tclUtils.tcl, package: tclu) and tk-part
# (file: tkUtils.tcl, package: tku).
#
#****

set version 0.9
package ifneeded tclu $version [list source [file join $dir tclUtils.tcl]]
package ifneeded tku  $version [list source [file join $dir tkUtils.tcl]]
