# ------------------------------------------------------------------------
#  initialize the ::guib namespace and set the library & version variables
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
#  load Tcl/Tk + [Incr Tcl]/[Incr Tk]/[Incr Widgets]
# ------------------------------------------------------------------------

package require Tk       
#package require Itcl
package require Itk      
package require Iwidgets 

# We need to import all of the itcl functions into the global
# namespace.
namespace import -force itcl::*     


namespace eval ::guib {
    namespace export module

    variable library [file dirname [info script]]
    variable version [gets [open [file join $library VERSION] r]]

    variable options    
    variable module
    variable library

    # associative-dimension "settings" is used for user-spefic seetings

    variable settings

    #------------------------------------------------------------------------
    # CONVENTION for settings() array elements:
    #------------------------------------------------------------------------
    # The elements of the settings() array have the following names:
    # settings(WHAT.details). For example, all elements that are
    # dealing with filename, have the following names:
    # settings(FILENAME.*)
    #------------------------------------------------------------------------            

    # regular expresion for the for the endlist-string
    set settings(NAMELIST.end_regexp) {^ *&end|^ */}

    # string to write for the end-of-namelist
    set settings(NAMELIST.end_string) { /}

    # case-sensitivity of namelists variable names (i.e., _nocase like
    # the "-nocase" tcl-option)    
    set settings(NAMELIST.varname_nocase) 1

    # format for printing namelist's variables names
    set settings(NAMELIST.varname_format) {   %15s}

    # support for undefined namelist variables (0|1)
    set settings(NAMELIST.variable_support_undefined) 0

    # do we quotre the string [charagetr*(*)] type of variables (0|1)
    set settings(NAMELIST.quote_strings) 1
    set settings(NAMELIST.quote_char)    '

    # case-sensitivity of input (0|1)
    set settings(INPUT.nocase) 0
    
    #---
    # For case-insensitive input ( $settings(INPUT.nocase) == 1 ):
    
    #   print the keywords as "upper|lower|unchanged"
    set settings(INPUT.nocase_preference_keyword)  upper; # not USED yet
    #   print the variable values as "upper|lower|unchanged"
    set settings(INPUT.nocase_preference_varvalue) unchanged; # not USED yet

    #---

    # whether to add a trailing slash to dirnames
    set settings(DIRNAME.trailing_slash)  0

    # whether to use only file-tail for filenames
    set settings(FILENAME.only_tail)  0


    # ------------------------------------------------------------------------
    #  IMAGES
    # ------------------------------------------------------------------------
    
    # create toolbar images ...    
    
    image create photo filenew -format gif \
	-file [file join $env(GUIB) images filenew2.gif]
    
    image create photo fileopen -format gif \
	-file [file join $env(GUIB) images fileopen2.gif]  
    
    image create photo filesave -format gif \
	-file [file join $env(GUIB) images filesave2.gif]  
    
    image create photo filesaveas -format gif \
	-file [file join $env(GUIB) images filesaveas2.gif]
    
    image create photo fileclose  -format gif \
	-file [file join $env(GUIB) images fileclose2.gif]
    
    image create photo exitApp -format gif \
	-file [file join $env(GUIB) images exit2.gif]     
}


# ------------------------------------------------------------------------
#  load TclLib:
#               we only need the cmdline package of TclLib
# ------------------------------------------------------------------------
lappend auto_path [file join $guib::library external lib]
package require cmdline 1.2


#------------------------------------------------------------------------
# load TCLU & TKU library
#------------------------------------------------------------------------
lappend auto_path [file join $guib::library lib]
package require tclu
package require tku

#------------------------------------------------------------------------
# load GUIB
#------------------------------------------------------------------------
lappend auto_path [file join $guib::library src]


# ------------------------------------------------------------------------
#  provide package GUIB
# ------------------------------------------------------------------------
package provide Guib $guib::version


namespace eval ::guib {

    # we should source the file guib-keywords-def.tcl in order to load
    # the definiton of GUIB keywords options
    source [file join $library src guib-keywords-def.tcl]

    # we should source the file widgets.itcl in order to load the
    # definiton of ::guib::widgets namespace variables !!!
    source [file join $library src widgets.itcl]

    source [file join $library src keywidgets.itcl]

    #source [file join $library src keywordObj.itcl]
    #source [file join $library src moduleObj.itcl]
}

