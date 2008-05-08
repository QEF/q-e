#
# $RCSfile: guib-keywords-def.tcl,v $ --
#
#      This file contains ...
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
# $Id: guib-keywords-def.tcl,v 1.5 2008-05-08 18:44:36 kokalj Exp $ 
#

#------------------------------------------------------------------------
#****v* guib.itcl/options
#  NAME
#    options -- this is the guib::options variables holding options of GUIB keywords !!!
#
#  DESCRIPTION
#    The allowed options of the GUIB keyowrds are stored in
# options(keyword) variable. For example, for the "var" keyword, the
# allowed options are stored in options(var) variable. The structure of
# the options(keyword) variables is accoring to the cmdline::getoptions 
# requirements.
# 
#  SOURCE


#
# BEWARE: if options has an argument we add ".arg" extension to it
#

#------------------------------------------------------------------------
# DEFINITION of standard "module" options.
#------------------------------------------------------------------------
set options(module) {
    {title.arg     GUI {Human readable name for the module}}
    {script.arg   {}   {Scirpt of the module}}
    {varscope.arg {}   {Scope for variables}}
}    



#------------------------------------------------------------------------
# DEFINITION of options for nameObj objects (line, namelist, group)
#------------------------------------------------------------------------
set options(page) {
    { name.arg      {}    {The name of the line} }
}
set options(line) {
    { name.arg      {}       {The name of the line} }
    { decor.arg     prefixed {Style of decoration-label: prefixed, normal or none} }
}
set options(namelist) {
    { name.arg      {}       {The name of the namelist} }
    { decor.arg     prefixed {Style of decoration-label: prefixed, normal or none} }
}
set options(group) {
    { name.arg      {}    {The name of the namelist} }
    { decor.arg     none  {Style of decoration-label: prefixed, normal or none} }
}



#------------------------------------------------------------------------
# DEFINITION of standard "var" options. Applies to the 
# following keywords: var, auxilvar, dimension
#------------------------------------------------------------------------
set standard_var_options {
    { variable.arg  {}    {The name of the variable} }
    { text.arg      {}    {Displayed top-text for the variable} }  
    { label.arg     {}	  {Displayed label for the variable} }
    { value.arg     {}    {Allowed values for the variable} }
    { textvalue.arg {}	  {Allowed values for the variable that will be diplayed in a widget like radiobutton} }
    { default.arg   {}    {Default value for the variable} }
    { fmt.arg       {}	  {Format speficication for the variable} }
    { validate.arg  {}    {The validation type for the value (i.e. numeric, alphabetic, integer, hexidecimal, real, and alphanumeric)} }
    { infmt.arg     {}	  {Input format speficication for the variable} }
    { outfmt.arg    {}    {Output format speficication for the variable} }
    { widget.arg    entry {Widget to be displayed for this variable} }
    { helptext.arg  {}    {help text for this var} }
    { helpfmt.arg   txt   {Format of the help.} } 
}   

set options(var)      $standard_var_options
set options(auxilvar) $standard_var_options

set options(scriptvar) {
    { value.arg     {}    {Allowed values for the variable} }
    { textvalue.arg {}	  {Allowed values for the variable that will be diplayed in a widget like radiobutton} }
}

#------------------------------------------------------------------------
# DEFINITION of additional dimension options
#------------------------------------------------------------------------
set special_dimension_options {
    {start.arg  0    {starting element of the dimension} }
    {end.arg    0    {last element of the dimension} }
    {pack.arg   top  {packing-side of the dimension widget} }
}
set options(dimension) [concat $standard_var_options $special_dimension_options]

#------------------------------------------------------------------------
# DEFINITION of table options
#------------------------------------------------------------------------
set options(table) {
    {caption.arg  {}    {The caption of the table}}
    {head.arg     {}    {The head of the table}}
    {variable.arg {}    {The name of the variable holding the table data}}
    {validate.arg {}    {Specification of the validation for the tables columns}}
    {cols.arg     1     {Number of columns}}
    {rows.arg     1     {Number of rows}}
    {widgets.arg  entry {List of widgets/per-columns to use for table cells}}
    {onvalues.arg  1    {List of on-values for the table cell widget like checkbutton}}
    {offvalues.arg 0    {List of off-values for the table cell widget like checkbutton}}
    {fmt.arg      {}    {Format speficication for the variable} }
    {infmt.arg    {}    {Input format speficication for the variable} }
    {outfmt.arg   {}    {Output format speficication for the variable} }
    {helptext.arg {}    {help text for this var.}}
    {helpfmt.arg  txt   {Format of the help.}}
    {optionalcols.arg -1 {Marks which columns are optional: if optionalcols > 0 then columns >= optionalcols are optional}}
}
 
#------------------------------------------------------------------------
# DEFINITION of help options
#------------------------------------------------------------------------
#{variable.arg  {}   {name of a variable for which this help is meant} }
set options(help) {    
    {helpfmt.arg   txt2html "Format of the help, must be one of html, txt2html, txt, or helpdoc"}
    {vartype.arg   {}       "Type of the variable"}
    {helptext.arg  {}       "Help-text"}
}

#------------------------------------------------------------------------
# definitions of separator options
#------------------------------------------------------------------------
set options(separator) {
    {label.org  {}  {Text for the label of separator widget}}
}


#------------------------------------------------------------------------
# definitions of text options
#------------------------------------------------------------------------
set options(text) {
    {caption.arg   {}    {The caption of the text.}}
    {label.arg     {}    {Text for the label of separator widget.}}
    {readvar.arg   {}    {Variable holding the content of the text upon the file-read. Requires the readfilter routine.}}
    {helptext.arg  {}    {help text for this var.}}
    {helpfmt.arg   txt   {Format of the help.}}
}

#****
#------------------------------------------------------------------------