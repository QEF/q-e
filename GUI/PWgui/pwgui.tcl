# ----------------------------------------------------------------------
#  PROGRAM: PWgui
#  PURPOSE: tries to be a GUI for the PWscf
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
    puts " PWGUI       : $env(PWGUI)"
    set guib_dir [glob -nocomplain -directory [file join $env(PWGUI) lib] Guib-*]
    if { $guib_dir != "" } {
	set env(GUIB) $guib_dir
    } else {
	# we arrive here, if we are using SVN version of code
	if { [file isdirectory [file join $env(PWGUI) .. Guib]] } {
	    puts "   "
	    puts "   It seems you are using SVN version of PWgui/Quantum-Espresso."
	    puts "   "
	    puts "   For the SVN version you need to do the following:"
	    puts "      * cd into GUI/PWgui directory, and"
	    puts "      * execute: make svninit"
	    puts "   "
	    exit
	}
    }

    if { [info exists env(GUIB)] } {
	lappend auto_path $env(GUIB)
        puts " GUIB engine : $env(GUIB)\n"
    }
} else {
    puts stderr "   "
    puts stderr "   Please define the PWGUI enviromental variable !!!"
    puts stderr "   PWGUI should point to the package root directory."
    puts stderr "   "
    exit
}

#
# all initialization stuff should go into init.tcl
#
source [file join $env(PWGUI) init.tcl]

#
# now go-ahead: launch the application
#
wm withdraw .
bind . <Destroy> ::guib::exitApp
source [file join $env(PWGUI) src pwscf.itcl]
