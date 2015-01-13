#!/bin/sh
# ----------------------------------------------------------------------
#  PROGRAM: pwgui_reformat
#  PURPOSE: reformat the input for nice lookout of a particular PWSCF module
# ----------------------------------------------------------------------
#  Anton Kokalj
#  Jozef Stefan Institute, Ljubljana, Slovenia
#  INFM DEMOCRITOS National Simulation Center, Trieste, Italy
#  Email: Tone.Kokalj@ijs.si
# ======================================================================
#  Copyright (c) 2003--2004 Anton Kokalj
# ======================================================================
#
#
# This file is distributed under the terms of the GNU General Public
# License. See the file `COPYING' in the root directory of the present
# distribution, or http://www.gnu.org/copyleft/gpl.txt .
#
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT.  IN NO EVENT SHALL ANTON KOKALJ BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

if { [info exists env(PWGUI)] } {
    set guib_dir [glob -nocomplain -directory [file join $env(PWGUI) lib] Guib-*]
    if { $guib_dir != "" } {
	set env(GUIB) $guib_dir
    }
    if { [info exists env(GUIB)] } {
	lappend auto_path $env(GUIB)
    }
    lappend auto_path [file join $env(PWGUI) src]
} else {
    puts stderr "   "
    puts stderr "   Please define the PWGUI enviromental variable !!!"
    puts stderr "   PWGUI should point to the package root directory."
    puts stderr "   "
    exit
}

#
# Usage: pwgui_reformat PWSCF_module input_file
#
if { $argc != 2 } {

    # look into the $PWGUI/modules for allowed PWSCF_modules

    cd [file join $env(PWGUI) modules]
    set modules ""
    foreach item [glob -nocomplain *] {
	if { [file isdirectory $item] && $item != "CVS" } {
	    append modules "          $item\n"
	}
    }
    puts stderr "   "
    puts stderr "   Usage: pwgui_reformat PWSCF_module input_file"
    puts stderr "   "
    puts stderr "   Valid PWSCF's modules (PWSCF_module) are:"
    puts stderr "   "
    puts stderr $modules
    exit
} else {
    set module [lindex $argv 0]
    set module [file join $::env(PWGUI) modules $module $module.tcl]
    set input  [lindex $argv 1]
    if { [file pathtype $input] != "absolute" } {
        set input [file join [pwd] $input]
    }
}
    
package require Guib 0.5
wm withdraw .
option readfile [file join $::guib::library guib.theme] startupFile

# DEBUGGING
set ::tclu::DEBUG 0
set ::tclu::DEBUG_FILE 0

#
# do the reformatting ....
#
namespace eval ::pwscf {}
namespace eval ::guib  {
    #set moduleObj [simpleTplwGUI $::module]
    set moduleObj [embedGUI $::module .]
    $moduleObj openFile $::input
    $moduleObj print 1
}

exit 0
