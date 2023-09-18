# ----------------------------------------------------------------------
#  PROGRAM: PWgui
#  PURPOSE: A GUI input builder for PWscf
# ----------------------------------------------------------------------
#  Anton Kokalj
#  Jozef Stefan Institute, Ljubljana, Slovenia
#  Email: tone.kokalj@ijs.si
# ======================================================================
#  Copyright (c) 2003--2023 Anton Kokalj
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


# PWgui version

set fid [open $env(PWGUI)/VERSION]
set ver [read -nonewline $fid]
close $fid


# make a list of which QE programs are supported

foreach item [glob -nocomplain [file join $env(PWGUI) modules *]] {
    if { [file isdirectory $item] } {
        set prog [file tail $item]
        if { $prog == "atomic" } {
            append progs "          [format {%-9s  or  %-9s} ld1.x atomic]\n"
        } else {
            append progs "          [format {%-9s  or  %-9s} $prog.x $prog]\n"
        }
        lappend prog_l [file rootname $prog]
    }
}

# usage

set usage [subst -nocommands {
    Usage:  pwgui [-h] |  [-r program inputFile] [-o program inputFile]

    Description of options:

    -h
    -?
    --help
          Print this usage message.

    -r  program  inputFile
    --reformat  program  inputFile
          Read the "inputFile" of the QE's "program" (e.g. pw.x, pp.x, ...),
          check its syntax for correctness, and print it nicely formatted.

    -o  program  inputFile
    --open  program  inputFile
          Open the "inputFile" of the QE's "program" (e.g. pw.x, pp.x, ...)
          in PWgui.

    Valid values for "program" are:

$progs
}]

set o_usage {
    Incorrect usage of $opt option, must be:
    
    pwgui $opt program inputFile
}


# parse command-line options

set reformat 0
set open     0

if { $argc > 0 } {
    set opt [lindex $argv 0]
    switch -glob -- $opt {
        -h - --help - -\\? {
            puts $usage
            exit 0
        }
        -r* - --ref* {
            # -r | --reformat
            if { $argc != 3 } { puts [subst $o_usage]; exit 1 }
            set reformat 1
        }
        -o* - --open {
            # -o | --open
            if { $argc != 3 } { puts [subst $o_usage]; exit 1 }
            set open 1
        }
        default {
            puts "\n    Incorrect option \"$opt\", aborting...\n$usage"
            exit 1
        }
    }

    set module    [file rootname [lindex $argv 1]]
    set inputFile [lindex $argv 2]
    if { ! [file exists $inputFile] } {
        puts "
   File \"$inputFile\" does not exist.\n"  
        exit 1
    }
    if { $module == "ld1" || $module == "atomic" } {
        set module      atomic
        set moduleIdent ld
        set moduleName  LD1.X
    } else {
        set moduleIdent $module
        set moduleName  [string toupper $module].X
    }
    set moduleFile [file join $env(PWGUI) modules $module $module.tcl]

    # check for valid module    
    
    if { [lsearch $prog_l $module] < 0 } {
        puts "[subst $o_usage]
    where program must be one of:
    
       [join $prog_l {, }]\n"
        exit 1
    }
}


if { [info exists env(PWGUI)] } {    
    set guib_dir [glob -nocomplain -directory [file join $env(PWGUI) lib] Guib-*]
    if { $guib_dir != "" } {
	set env(GUIB) $guib_dir
    } else {
        set guib_dir [file normalize [file join $env(PWGUI) .. Guib]]
	if { [file isdirectory $guib_dir] } {
            # we arrive here, if we are using PWgui inside the QE
            set env(GUIB) $guib_dir
        }
    }

    if { [info exists env(GUIB)] } {
	lappend auto_path $env(GUIB)
    } else {
        puts "
 Guib engine was not found. 
 You may consider to defined GUIB enviromental variable that points to Guib engine.

 Aborting.
"
        exit 1
    }

    if { ! $reformat } {
        puts "
 ==================================================
  PWgui version :   $ver
 --------------------------------------------------

 PWGUI       : $env(PWGUI)
 GUIB engine : $env(GUIB)\n"
    }
    
} else {
    puts stderr " "
    puts stderr " Please define the PWGUI enviromental variable !!!"
    puts stderr " PWGUI should point to the PWgui root directory."
    puts stderr " "
    exit 1
}


#
# all initialization stuff should go into init.tcl
#
source [file join $env(PWGUI) init.tcl]

wm withdraw .

if { ! $reformat } {
    #
    # open a GUI main window
    #
    bind . <Destroy> ::guib::exitApp
    source [file join $env(PWGUI) src pwscf.itcl]

    if { $open } {
        # open the specified input file
        $gui openInput $moduleIdent $moduleName $moduleFile $inputFile
    }
} else {
    #
    # only reformat the specified input file
    #
    set ::tclu::DEBUG 0
    set ::tclu::DEBUG_FILE 0
    option readfile [file join $::guib::library guib.theme] startupFile

    namespace eval ::guib  {
        set moduleObj [embedGUI $::moduleFile .]
        $moduleObj openFile $::inputFile
        $moduleObj print 1
    }
    exit 0
}
