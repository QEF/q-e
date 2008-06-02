#
# Initialization stuff for PWgui
#

package require Guib 0.5

namespace eval ::pwscf {
    variable pwscf
    variable settings

    set pwscf(PWD) [pwd]
}

# define here all pwscf's namespaces ...
namespace eval ::pwscf::edit      {
    variable edit
}
namespace eval ::pwscf::menustate {}
namespace eval ::pwscf::view      {}
namespace eval ::pwscf::run       {
    variable run
    set run(mode) nonblocking ; # possibilities: nonblocking || background
}

# load settings file ...
source $env(PWGUI)/pwgui.settings
if { [file exists $env(HOME)/.pwgui/pwgui.settings] } {
    # overwritte default settings by user-settings
    source $env(HOME)/.pwgui/pwgui.settings
}

lappend auto_path [file join $env(PWGUI) src]
