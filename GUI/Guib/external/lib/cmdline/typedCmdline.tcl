# typedCmdline.tcl --
#
#    This package provides a utility for parsing typed command
#    line arguments that may be processed by various applications.
#
# Copyright (c) 2000 by Ross Palmer Mohn.
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
# 
# RCS: @(#) $Id: typedCmdline.tcl,v 1.1 2004-02-18 11:29:17 kokalj Exp $

namespace eval cmdline {
    namespace export typedGetopt typedGetoptions typedUsage

    # variable cmdline::charclasses --
    #
    #    Create regexp list of allowable character classes
    #    from "string is" error message.
    #
    # Results:
    #    String of character class names separated by "|" characters.

    variable charclasses
    catch {string is . .} charclasses
    regexp -- {must be (.+)$} $charclasses dummy charclasses
    regsub -all -- {, (or )?} $charclasses {|} charclasses

}

# cmdline::typedGetopt --
#
#	The cmdline::typedGetopt works in a fashion like the standard
#	C based getopt function.  Given an option string and a
#	pointer to a list of args this command will process the
#	first argument and return info on how to procede. In addition,
#	you may specify a type for the argument to each option.
#
# Arguments:
#	argvVar		Name of the argv list that you want to process.
#			If options are found, the arg list is modified
#			and the processed arguments are removed from the
#			start of the list.
#
#	optstring	A list of command options that the application
#			will accept.  If the option ends in ".xxx", where
#			xxx is any valid character class to the tcl
#			command "string is", then typedGetopt routine will
#			use the next argument as a typed argument to the
#			option. The argument must match the specified
#			character classes (e.g. integer, double, boolean,
#			xdigit, etc.). Alternatively, you may specify
#			".arg" for an untyped argument.
#
#	optVar		Upon success, the variable pointed to by optVar
#			contains the option that was found (without the
#			leading '-' and without the .xxx extension).  If
#			typedGetopt fails the variable is set to the empty
#			string. SOMETIMES! Different for each -value!
#
#	argVar		Upon success, the variable pointed to by argVar
#			contains the argument for the specified option.
#			If typedGetopt fails, the variable is filled with
#			an error message.
#
# Argument type syntax:
#	Option that takes no argument.
#		foo
#
#	Option that takes a typeless argument.
#		foo.arg
#
#	Option that takes a typed argument. Allowable types are all
#	valid character classes to the tcl command "string is".
#	Currently must be one of alnum, alpha, ascii, control,
#	boolean, digit, double, false, graph, integer, lower, print,
#	punct, space, true, upper, wordchar, or xdigit.
#		foo.double
#
#	Option that takes an argument from a list.
#		foo.(bar|blat)
#
# Argument quantifier syntax:
#	Option that takes an optional argument.
#		foo.arg?
#
#	Option that takes a list of arguments terminated by "--".
#		foo.arg+
#
#	Option that takes an optional list of arguments terminated by "--".
#		foo.arg*
#
#	Argument quantifiers work on all argument types, so, for
#	example, the following is a valid option specification.
#		foo.(bar|blat|blah)?
#
# Argument syntax miscellany:
#	Options may be specified on the command line using a unique,
#	shortened version of the option name. Given that program foo
#	has an option list of {bar.alpha blah.arg blat.double},
#	"foo -b fob" returns an error, but "foo -ba fob"
#	successfully returns {bar fob}
#
# Results:
#	The typedGetopt function returns one of the following:
#	 1	a valid option was found
#	 0	no more options found to process
#	-1	invalid option
#	-2	missing argument to a valid option
#	-3	argument to a valid option does not match type
#
# Known Bugs:
#	When using options which include special glob characters,
#	you must use the exact option. Abbreviating it can cause
#	an error in the "cmdline::prefixSearch" procedure.

proc cmdline::typedGetopt {argvVar optstring optVar argVar} {
    variable charclasses

    upvar $argvVar argsList

    upvar $optVar retvar
    upvar $argVar optarg

    # default settings for a normal return
    set optarg ""
    set retvar ""
    set retval 0

    # check if we're past the end of the args list
    if {[llength $argsList] != 0} {

        # if we got -- or an option that doesn't begin with -, return (skipping
        # the --).  otherwise process the option arg.
        switch -glob -- [set arg [lindex $argsList 0]] {
            "--" {
                set argsList [lrange $argsList 1 end]
            }

            "-*" {
                # Create list of options without their argument extentions

                set optstr ""
                foreach str $optstring {
                    lappend optstr [file rootname $str]
                }

                set _opt [string range $arg 1 end]

                set i [cmdline::prefixSearch $optstr [file rootname $_opt]]
                if {$i != -1} {
                    set opt [lindex $optstring $i]

                    set quantifier "none"
                    if {[regexp -- {\.[^.]+([?+*])$} $opt dummy quantifier]} {
                        set opt [string range $opt 0 end-1]
                    }

                    if {[string first . $opt] == -1} {
                        set retval 1
                        set retvar $opt
                        set argsList [lrange $argsList 1 end]

                    } elseif {[regexp -- "\\.(arg|$charclasses)\$" $opt dummy charclass]
                            || [regexp -- {\.\(([^)]+)\)} $opt dummy charclass]} {
				if {[string equal arg $charclass]} {
                            set type arg
			} elseif {[regexp -- "^($charclasses)\$" $charclass]} {
                            set type class
                        } else {
                            set type oneof
                        }

                        set argsList [lrange $argsList 1 end]
                        set opt [file rootname $opt]

                        while {1} {
                            if {[llength $argsList] == 0
                                    || [string equal "--" [lindex $argsList 0]]} {
                                if {[string equal "--" [lindex $argsList 0]]} {
                                    set argsList [lrange $argsList 1 end]
                                }

                                set oneof ""
                                if {$type == "arg"} {
                                    set charclass an
                                } elseif {$type == "oneof"} {
                                    set oneof ", one of $charclass"
                                    set charclass an
                                }
    
                                if {$quantifier == "?"} {
                                    set retval 1
                                    set retvar $opt
                                    set optarg ""
                                } elseif {$quantifier == "+"} {
                                    set retvar $opt
                                    if {[llength $optarg] < 1} {
                                        set retval -2
                                        set optarg "Option requires at least one $charclass argument$oneof -- $opt"
                                    } else {
                                        set retval 1
                                    }
                                } elseif {$quantifier == "*"} {
                                    set retval 1
                                    set retvar $opt
                                } else {
                                    set optarg "Option requires $charclass argument$oneof -- $opt"
                                    set retvar $opt
                                    set retval -2
                                }
                                set quantifier ""
                            } elseif {($type == "arg")
                                    || (($type == "oneof")
                                    && [string first "|[lindex $argsList 0]|" "|$charclass|"] != -1)
                                    || (($type == "class")
                                    && [string is $charclass [lindex $argsList 0]])} {
                                set retval 1
                                set retvar $opt
                                lappend optarg [lindex $argsList 0]
                                set argsList [lrange $argsList 1 end]
                            } else {
                                set oneof ""
                                if {$type == "arg"} {
                                    set charclass an
                                } elseif {$type == "oneof"} {
                                    set oneof ", one of $charclass"
                                    set charclass an
                                }
                                set optarg "Option requires $charclass argument$oneof -- $opt"
                                set retvar $opt
                                set retval -3
    
                                if {$quantifier == "?"} {
                                    set retval 1
                                    set optarg ""
                                }
                                set quantifier ""
                            }
                             if {![regexp -- {[+*]} $quantifier]} {
                                break;
                            }
                        }
                    } else {
                        error "Illegal option type specification:\
                                must be one of $charclasses"
                    }
                } else {
                    set optarg "Illegal option -- $_opt"
                    set retvar $_opt
                    set retval -1
                }
            }
	    default {
		# Skip ahead
	    }
        }
    }

    return $retval
}

# cmdline::typedGetoptions --
#
#	Process a set of command line options, filling in defaults
#	for those not specified. This also generates an error message
#	that lists the allowed options if an incorrect option is
#	specified.
#
# Arguments:
#	arglistVar	The name of the argument list, typically argv
#	optlist		A list-of-lists where each element specifies an option
#			in the form:
#
#				option default comment
#
#			Options formatting is as described for the optstring
#			argument of typedGetopt. Default is for optionally
#			specifying a default value. Comment is for optionally
#			specifying a comment for the usage display. The
#			options "-help" and "-?" are automatically included
#			in optlist.
#
# Argument syntax miscellany:
#	Options formatting and syntax is as described in typedGetopt.
#	There are two additional suffixes that may be applied when
#	passing options to typedGetoptions.
#
#	You may add ".multi" as a suffix to any option. For options
#	that take an argument, this means that the option may be used
#	more than once on the command line and that each additional
#	argument will be appended to a list, which is then returned
#	to the application.
#		foo.double.multi
#
#	If a non-argument option is specified as ".multi", it is
#	toggled on and off for each time it is used on the command
#	line.
#		foo.multi
#
#	If an option specification does not contain the ".multi"
#	suffix, it is not an error to use an option more than once.
#	In this case, the behavior for options with arguments is that
#	the last argument is the one that will be returned. For
#	options that do not take arguments, using them more than once
#	has no additional effect.
#
#	Options may also be hidden from the usage display by
#	appending the suffix ".secret" to any option specification.
#	Please note that the ".secret" suffix must be the last suffix,
#	after any argument type specification and ".multi" suffix.
#		foo.xdigit.multi.secret
#
# Results
#	Name value pairs suitable for using with array set.

proc cmdline::typedGetoptions {arglistVar optlist {usage options:}} {
    variable charclasses

    upvar 1 $arglistVar argv

    set opts {? help}
    foreach opt $optlist {
        set name [lindex $opt 0]
        if {[regsub -- {\.secret$} $name {} name] == 1} {
            # Remove this extension before passing to typedGetopt.
        }
        if {[regsub -- {\.multi$} $name {} name] == 1} {
            # Remove this extension before passing to typedGetopt.

            regsub -- {\..*$} $name {} temp
            set multi($temp) 1
        }
        lappend opts $name
        if {[regsub -- "\\.(arg|$charclasses|\\(.+).?\$" $name {} name] == 1} {
            # Set defaults for those that take values.
            # Booleans are set just by being present, or not

            set dflt [lindex $opt 1]
            if {$dflt != {}} {
                set defaults($name) $dflt
            }
        }
    }
    set argc [llength $argv]
    while {[set err [cmdline::typedGetopt argv $opts opt arg]]} {
        if {$err == 1} {
            if {[info exist result($opt)]
                    && [info exist multi($opt)]} {
                # Toggle boolean options or append new arguments

                if {$arg == ""} {
                    unset result($opt)
                } else {
                    set result($opt) "$result($opt) $arg"
                }
            } else {
                set result($opt) "$arg"
            }
        } elseif {($err == -1) || ($err == -3)} {
            error [cmdline::typedUsage $optlist $usage]
        } elseif {$err == -2 && ![info exist defaults($opt)]} {
            error [cmdline::typedUsage $optlist $usage]
        }
    }
    if {[info exist result(?)] || [info exists result(help)]} {
        error [cmdline::typedUsage $optlist $usage]
    }
    foreach {opt dflt} [array get defaults] {
        if {![info exist result($opt)]} {
            set result($opt) $dflt
        }
    }
    return [array get result]
}

# cmdline::typedUsage --
#
#	Generate an error message that lists the allowed flags,
#	type of argument taken (if any), default value (if any),
#	and an optional description.
#
# Arguments:
#	optlist		As for cmdline::typedGetoptions
#
# Results
#	A formatted usage message

proc cmdline::typedUsage {optlist {usage {options:}}} {
    variable charclasses

    set str "[cmdline::getArgv0] $usage\n"
    foreach opt [concat $optlist \
            {{help "Print this message"} {? "Print this message"}}] {
        set name [lindex $opt 0]
        if {[regsub -- {\.secret$} $name {} name] == 1} {
            # Hidden option

        } else {
            if {[regsub -- {\.multi$} $name {} name] == 1} {
                # Display something about multiple options
            }

            if {[regexp -- "\\.(arg|$charclasses)\$" $name dummy charclass]
                    || [regexp -- {\.\(([^)]+)\)} $opt dummy charclass]} {
                   regsub -- "\\..+\$" $name {} name
                set comment [lindex $opt 2]
                set default "<[lindex $opt 1]>"
                if {$default == "<>"} {
                    set default ""
                }
                append str [format " %-20s %s %s\n" "-$name $charclass" \
                        $comment $default]
            } else {
                set comment [lindex $opt 1]
		append str [format " %-20s %s\n" "-$name" $comment]
            }
        }
    }
    return $str
}

# cmdline::prefixSearch --
#
#	Search a Tcl list for a pattern; searches first for an exact match,
#	and if that fails, for a unique prefix that matches the pattern 
#	(ie, first "lsearch -exact", then "lsearch -glob $pattern*"
#
# Arguments:
#	list		list of words
#	pattern		word to search for
#
# Results:
#	Index of found word is returned. If no exact match or
#	unique short version is found then -1 is returned.

proc cmdline::prefixSearch {list pattern} {
    # Check for an exact match

    if {[set pos [::lsearch -exact $list $pattern]] > -1} {
        return $pos
    }

    # Check for a unique short version

    set slist [lsort $list]
    if {[set pos [::lsearch -glob $slist $pattern*]] > -1} {
        # What if there is nothting for the check variable?

        set check [lindex $slist [expr {$pos + 1}]]
        if {[string first $pattern $check] != 0} {
            return [::lsearch -exact $list [lindex $slist $pos]]
        }
    }
    return -1
}
