#!/bin/sh
# next line restarts wish \
exec wish

# ------------------------------------------------------------------------
#  INITIALIZATION
# ------------------------------------------------------------------------

# check if GUIB environmental variable is defined

if { [info exists env(GUIB)] } {
    # instruct the Tcl where to search for GUIB package
    lappend auto_path [file join $env(GUIB)]
} else {
    # we assume that GUIB package is on some "standard" path and Tcl
    # will find it.
}

# load a Guib package
package require Guib 0.1.1

# withdraw the "." toplevel window, and bind the <Destroy> event
# to ::guib::exitApp
wm withdraw .
bind . <Destroy> ::guib::exitApp


# ------------------------------------------------------------------------
#  GUI construction
# ------------------------------------------------------------------------

# construct the GUI object
set gui [::guib::GUI \#auto -title "My 1st GUI" -appname MyGUI]

# Add modules. The syntax is:
# ------------
# $gui addModule module $moduleID $moduleLabel $moduleFile
#
$gui addModule module inp1 "My Input-1" myInput-1.def
$gui addModule module inp2 "My Input-2" myInput-2.def

# Add help. The syntax is:
# ---------
# $gui addHelp help $helpID $helpLabel $helpFile
#
$gui addHelp help help1 "Help for Input-1" myInput-1.html
$gui addHelp help help2 "Help for Input-2" myInput-2.html


#
# some extra configuration of the GUI
#
$gui extra {
    #
    # add a logo to toolbar panel
    #
    image create photo myLogo -format gif -file myLogo.gif

    set tb   [component toolbar]
    set logo [$tb add label logo -image myLogo]
    pack configure $logo -side right
}