#
# $RCSfile: guib.tcl,v $ --
#
#      This is the main application file, i.e., might be called as
#      "wish guib.itcl".
#
# Copyright (c) 2003  Anton Kokalj   Email: tone.kokalj@ijs.si
#
#
# This file is distributed under the terms of the GNU General Public
# License. See the file `COPYING' in the root directory of the present
# distribution, or http://www.gnu.org/copyleft/gpl.txt .
#
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# ANTON KOKALJ BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
# AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#
# $Id: guib.tcl,v 1.3 2008-02-15 16:50:18 kokalj Exp $ 
#

if { [info exists env(GUIB)] } {
    lappend auto_path [file join $env(GUIB)]
} else {
	puts stderr "   "
	puts stderr "   Please define the GUIB enviromental variable !!!"
	puts stderr "   GUIB should point to the package root directory."
	puts stderr "   "
	exit
}

package require Guib
wm withdraw .

option readfile [file join $::guib::library guib.theme] startupFile

if { $argc < 1 } {
    puts stderr ""
    puts stderr "Usage:  guib  module-file"
    puts stderr ""
    exit
}

set moduleFile [lindex $argv 0]

set ::tclu::DEBUG      0
set ::tclu::DEBUG_FILE 0

namespace eval ::guib {    
    simpleTplwGUI $::moduleFile
}
