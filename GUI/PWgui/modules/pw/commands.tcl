
# ------------------------------------------------------------------------
#  ::pwscf::pwSelectPseudoDir --
# ------------------------------------------------------------------------

proc ::pwscf::pwSelectPseudoDir {moduleObj} {
    variable pwscf
    global env

    set _dir [$moduleObj varvalue pseudo_dir]
    if { [file isdirectory $_dir] } {
	set dir $_dir
    } elseif { [info exists pwscf($moduleObj,LASTDIR,pseudo_dir)] } {
	set dir $pwscf($moduleObj,LASTDIR,pseudo_dir)
    } else {
	set dir [file join $env(HOME) pw pseudo]
	if { ! [file exists $dir] } {
	    set dir pwscf(PWD)
	}
    }

    set dir [tk_chooseDirectory -initialdir $dir \
		 -title "Chose directory" -mustexist 0]
    if { $dir == "" } {
	return
    }
    set pwscf($moduleObj,LASTDIR,pseudo_dir) $dir    
    
    # add a trailing slash "/"
    set dir [string trimright [file join $dir _] _]
    # add a quotes
    set dir '$dir'
    $moduleObj varset pseudo_dir -value $dir

}


# ------------------------------------------------------------------------
#  ::pwscf::pwSelectPPFile_
#
#  This "master" routine is used by ::pwscf::pwSelectPseudopotential &
#  ::pwscf::pwSelectMOLfile because both PP and MOL files are read
#  from pseudo_dir
#  ------------------------------------------------------------------------

proc ::pwscf::pwSelectPPFile_ {moduleObj variable ir ic label} {
     variable pwscf
    global env
        
    set _dir [string trim [$moduleObj varvalue pseudo_dir] "'"]
    if { [file isdirectory $_dir] } {
	set dir $_dir
    } elseif { [info exists pwscf($moduleObj,LASTDIR,pseudopotential)] } {
	set dir $pwscf($moduleObj,LASTDIR,pseudopotential)
    } else {
	set dir [file join $env(HOME) pw pseudo]
	if { ! [file isdirectory $dir] } {
	    set dir $pwscf(PWD)
	}    
    }
        
    set file [tk_getOpenFile \
		  -initialdir $dir \
		  -title      $label]    
    if { $file == "" } {
	return
    }
    set pwscf($moduleObj,LASTDIR,pseudopotential) [file dirname $file]
    
    $moduleObj varset ${variable}($ir,$ic) -value [file tail $file]    
}
proc ::pwscf::pwSelectPseudopotential {moduleObj variable ir ic} {
    ::pwscf::pwSelectPPFile_ $moduleObj $variable $ir $ic "Select a Pseudopotential File" 
}
proc ::pwscf::pwSelectMOLfile {moduleObj variable ir ic} {
    # MOL files are read from psudo_dir, hence use the ::pwscf::pwSelectPPFile_ routine
    ::pwscf::pwSelectPPFile_ $moduleObj $variable $ir $ic "Select a Pseudopotential File" 
}


# ------------------------------------------------------------------------
#  ::pwscf::pwLoadAtomCoord --
# 
# This function loads atomic coordinates from file. It supposes the
# PWSCF syntax, i.e., searches for ATOMIC_POSTIONS keyword. This means
# that atomic coordinates can be extracted from PWSCF's input file.
# ------------------------------------------------------------------------
proc ::pwscf::pwLoadAtomCoor {moduleObj} {
    set file [pwLoadAtomCoor:_init $moduleObj]
    if { $file != "" } {
	pwLoadAtomCoor:_read $moduleObj $file atomic_coordinates
    }
}
proc ::pwscf::pwLoadAtomCoor:_init {moduleObj} {
    variable pwscf

    if { [info exists pwscf($moduleObj,LASTDIR,atomic_coor)] } {
	set dir $pwscf($moduleObj,LASTDIR,atomic_coor)
    } else {
	set dir pwscf(PWD)	
    }
    
    # query the filename
    
    set file [tk_getOpenFile \
		  -initialdir $dir \
		  -title "Load Atomic Coordinates"]
    if { $file == "" } {
	uplevel 1 return
    }
    set pwscf($moduleObj,LASTDIR,atomic_coor) [file dirname $file]
    return $file
}
proc ::pwscf::pwLoadAtomCoor:_read {moduleObj file coorVar} {
    variable pwscf

    set _readCoor 0
    set IA(0)     0
    set ntyp_list {}
    set _UNIT     ""    
    set nimage 0

    # read the file

    set channel [open $file r]
        
    while {1} {
	set res [gets $channel _line]
	if { $res == -1 && $_line == "" } {
	    # end of file occurred
	    break
	}    

	set _len [llength $_line]

	# skip empty lines
	if { $_len == 0 } {
	    continue
	}
	
	if { [string match "ATOMIC_POSITIONS*" $_line] } {
	    set _line [readFilter::purifyCardLine $_line]
	    set _UNIT [lindex $_line 1]
	    set _readCoor 1
	    continue
	}	
	if { $_readCoor } {
	    if { $_len == 4 || $_len == 7 } {
		# read coordinates
		incr IA($nimage)
		set ia $IA($nimage)
		set len($ia,$nimage) $_len
		for {set i 1} {$i <= $_len} {incr i} {
		    set Atoms($ia,$i,$nimage) [lindex $_line [expr $i - 1]]
		}
		if { [lsearch -exact $ntyp_list $Atoms($ia,1,$nimage)] == -1 } {
		    lappend ntyp_list $Atoms($ia,1,$nimage)
		}
	    } elseif { $_len == 1 } {
		# might be first_image/intermediate_image/last_image string (NEB || SMD)
		if { [string match *first_image* $_line] } {
		    # do nothing so far ...
		    ;
		} elseif { [string match *intermediate_image* $_line] } {
		    incr nimage
		    set IA($nimage) 0
		} elseif { [string match *last_image* $_line] } {
		    incr nimage
		    set IA($nimage) 0
		} else {
		    # no, it is not first_image/intermediate_image/last_image string 
		    break
		}
	    } else {
		# record does not match; probably end-of-ATOMIC_POSITIONS
		break
	    }
	}
    }
    #/reading done


    # assign the "atomic-position" variables
    
    set ia $IA(0)
    for {set i 1} {$i <= $nimage} {incr i} {
	if { $IA($nimage) > $ia } { set ia $IA($nimage) }
    }
    $moduleObj varset nat -value $ia
    #$moduleObj varset ntyp -value $NTYP

    if { $_UNIT != "" } {
	$moduleObj varset ATOMIC_POSITIONS_flags -value [$moduleObj valueToTextvalue ATOMIC_POSITIONS_flags $_UNIT]
    }

    # load atomic-labels from ntyp_list if the "atomic_species" is not yet defined
    
    set ntyp [llength $ntyp_list]
    $moduleObj varset ntyp -value $ntyp
    for {set i 1} {$i <= $ntyp} {incr i} {
	set empty($i) 1
    }
    foreach type $ntyp_list {
	set new_type 1
	for {set i 1} {$i <= $ntyp} {incr i} {
	    if { $type == [$moduleObj varvalue atomic_species($i,1)] } {
		set new_type  0
		set empty($i) 0
	    }
	}
	if { $new_type } {
	    for {set i 1} {$i <= $ntyp} {incr i} {
		if { $empty($i) } {
		    $moduleObj varset atomic_species($i,1) -value $type
		    set empty($i) 0
		    break
		}
	    }
	}
    }

    # load the "atomic_coordinates" or "atomic_coordinates_inter_*" or "atomic_coordinates_last" table 

    if { $coorVar == "atomic_coordinates"      } { set start_index 0 }
    if { $coorVar == "atomic_coordinates_last" } { set start_index $nimage }
    if { [string match "atomic_coordinates_*_inter" $coorVar] } { 
	set start_index [regsub {[a-zA-Z_]*} $coorVar {}]
	if { $start_index > $nimage } { 
	    set start_index $nimage 
	}
    }
	    
    for {set ii $start_index} {$ii <= $nimage} {incr ii} { 
	# loop over images ...

	for {set ia 1} {$ia <= $IA($ii)} {incr ia} {
	    # loop over atoms ...

	    if { ! [info exists len($ia,$ii)] } { 
		break 
	    }

	    for {set i 1} {$i <= $len($ia,$ii)} {incr i} {
		# loop over fields ...

		if { ! [info exists Atoms($ia,$i,$ii)] } { 
		    break 
		}
		
		if { $coorVar == "atomic_coordinates" } {
		    # in this case load all the coordinates

		    if { $ii == 0 } {
			$moduleObj varset "atomic_coordinates($ia,$i)"  -value $Atoms($ia,$i,$ii)
		    } elseif { $ii == $nimage } {
			$moduleObj varset "atomic_coordinates_last($ia,$i)"  -value $Atoms($ia,$i,$ii)
		    } else {			
			$moduleObj varset "atomic_coordinates_${ii}_inter($ia,$i)"  -value $Atoms($ia,$i,$ii)
		    }
		} elseif { [string match "atomic_coordinates_*_inter" $coorVar] } {
		    # in this case load only the coordinates from start_index --> nimage
		    
		    if { $ii < $start_index || $start_index == $nimage } {
			set ind $ii
			# pathological case: ind == 0 , but smallest possible is 1
			if { $ind == 0 } { set ind 1 }
			$moduleObj varset "atomic_coordinates_${ind}_inter($ia,$i)"  -value $Atoms($ia,$i,$ii)
		    } else {
			$moduleObj varset "atomic_coordinates_last($ia,$i)"  -value $Atoms($ia,$i,$ii)
		    }
		} elseif { $coorVar == "atomic_coordinates_last" } {
		    # load only the last coordinates

		    $moduleObj varset "atomic_coordinates_last($ia,$i)"  -value $Atoms($ia,$i,$ii)
		}
	    }
	}
    }
}


# ------------------------------------------------------------------------
#  ::pwscf::pwLoadKPoints --
# ------------------------------------------------------------------------

proc ::pwscf::pwLoadKPoints {moduleObj} {
    variable pwscf

    if { [info exists pwscf($moduleObj,LASTDIR,k_points)] } {
	set dir $pwscf($moduleObj,LASTDIR,k_points)
    } else {
	set dir pwscf(PWD)	
    }

    #
    # query filename
    #
    set file [tk_getOpenFile -initialdir [pwd] -title "Load K-Points"]
    if { $file == "" } {
	return
    }
    set pwscf($moduleObj,LASTDIR,k_points) [file dirname $file]

    #
    # read the file
    #
    set channel [open $file r]
    # find the K_POINTS card
    while {1} {
	set _line [_getsNonEmptyLine $channel]
	
	if { [string match "K_POINTS*" $_line] } {
	    set _line [readFilter::purifyCardLine $_line]
	    set _UNIT [lindex $_line 1]
	    # assing the K_POINTS_flags variable
	    $moduleObj varset K_POINTS_flags -value [$moduleObj valueToTextvalue K_POINTS_flags $_UNIT]
	    break
	}
    }	    
    # read NKS
    set NKS [_getsNonEmptyLine $channel]
    if { [string is integer $NKS] } {
	$moduleObj varset nks -value $NKS
    } else {
	# TODO: raise an error
	return
    }
    # read K-POINTS
    for {set ia 1} {$ia <= $NKS} {incr ia} {
	set _line [_getsNonEmptyLine $channel]
	if { [llength $_line] != 4 } {
	    # TODO: raise an error
	}
	for {set i 1} {$i <= 4} {incr i} {
	    $moduleObj varset kpoints($ia,$i) -value [lindex $_line [expr $i - 1]]
	}
    }
}


# ------------------------------------------------------------------------
#  ::pwscf::pwLoadAddKPoints --
# ------------------------------------------------------------------------

proc ::pwscf::pwLoadAddKPoints {moduleObj} {
    variable pwscf

    if { [info exists pwscf($moduleObj,LASTDIR,add_k_points)] } {
	set dir $pwscf($moduleObj,LASTDIR,add_k_points)
    } else {
	set dir pwscf(PWD)	
    }

    #
    # query filename
    #
    set file [tk_getOpenFile -initialdir [pwd] -title "Load Additional K-Points"]
    if { $file == "" } {
	return
    }
    set pwscf($moduleObj,LASTDIR,k_points) [file dirname $file]

    #
    # read the file
    #
    set channel [open $file r]
    # find the K_POINTS card
    while {1} {
	set _line [_getsNonEmptyLine $channel]

	if { [string match "ADDITIONAL_K_POINTS*" $_line] } {
	    set _line [readFilter::purifyCardLine $_line]
	    set _UNIT [lindex $_line 1]
	    # assing the ADDITIONAL_K_POINTS_flags variable
	    $moduleObj varset ADDITIONAL_K_POINTS_flags -value [$moduleObj valueToTextvalue ADDITIONAL_K_POINTS_flags $_UNIT]
	    break
	}
    }	    
    # read NKS
    set NKS [_getsNonEmptyLine $channel]
    if { [string is integer $NKS] } {
	$moduleObj varset nks_add -value $NKS
    } else {
	# TODO: raise an error
	return
    }
    # read Add. K-POINTS
    for {set ia 1} {$ia <= $NKS} {incr ia} {
	set _line [_getsNonEmptyLine $channel]
	if { [llength $_line] != 4 } {
	    # TODO: raise an error
	}
	for {set i 1} {$i <= 4} {incr i} {
	    $moduleObj varset add_kpoints($ia,$i) -value [lindex $_line [expr $i - 1]]
	}
    }
}


# ------------------------------------------------------------------------
#  ::pwscf::pwLoadAtomicForces --
# ------------------------------------------------------------------------

proc ::pwscf::pwLoadAtomicForces {moduleObj} {
    variable pwscf

    if { [info exists pwscf($moduleObj,LASTDIR,atomic_forces)] } {
	set dir $pwscf($moduleObj,LASTDIR,atomic_forces)
    } else {
	set dir pwscf(PWD)	
    }

    #
    # query filename
    #
    set file [tk_getOpenFile -initialdir [pwd] -title "Load Atomic Forces"]
    if { $file == "" } {
	return
    }
    set pwscf($moduleObj,LASTDIR,atomic_forces) [file dirname $file]
    
    #
    # read the file
    #
    set channel [open $file r]

    # find the ATOMIC_FORCES card

    set ifor 0
    set _readForces 0
    
    while {1} {
	
	set _line [_getsNonEmptyLine $channel]
	
	if { [string match "ATOMIC_FORCES" $_line] } {
	    set _readForces 1
	    continue
	} 
	
	if { $_readForces } {
	    
	    incr ifor
	    
	    if {[llength $_line] == 4 } {
		$moduleObj varset atomic_forces($ifor,1)  -value [lindex $_line 0]
		$moduleObj varset atomic_forces($ifor,2)  -value [lindex $_line 1]
		$moduleObj varset atomic_forces($ifor,3)  -value [lindex $_line 2]
		$moduleObj varset atomic_forces($ifor,4)  -value [lindex $_line 3]
	    } else {
		# TODO: raise an error
	    }
	}
    }
}


# ------------------------------------------------------------------------
#  ::pwscf::pwLoadAtomicVelocities --
# ------------------------------------------------------------------------

proc ::pwscf::pwLoadAtomicVelocities {moduleObj} {
    variable pwscf

    if { [info exists pwscf($moduleObj,LASTDIR,atomic_velocities)] } {
	set dir $pwscf($moduleObj,LASTDIR,atomic_velocities)
    } else {
	set dir pwscf(PWD)	
    }

    #
    # query filename
    #
    set file [tk_getOpenFile -initialdir [pwd] -title "Load Atomic Velocities"]
    if { $file == "" } {
	return
    }
    set pwscf($moduleObj,LASTDIR,atomic_velocities) [file dirname $file]
    
    #
    # read the file
    #
    set channel [open $file r]

    # find the ATOMIC_VELOCITIES card

    set ifor 0
    set _readVelocities 0
    
    while {1} {
	
	set _line [_getsNonEmptyLine $channel]
	
	if { [string match "ATOMIC_VELOCITIES" $_line] } {
	    set _readVelocities 1
	    continue
	} 
	
	if { $_readVelocities } {
	    
	    incr ifor
	    
	    if {[llength $_line] == 4 } {
		$moduleObj varset atomic_velocities($ifor,1)  -value [lindex $_line 0]
		$moduleObj varset atomic_velocities($ifor,2)  -value [lindex $_line 1]
		$moduleObj varset atomic_velocities($ifor,3)  -value [lindex $_line 2]
		$moduleObj varset atomic_velocities($ifor,4)  -value [lindex $_line 3]
	    } else {
		# TODO: raise an error
	    }
	}
    }
}

# ------------------------------------------------------------------------
#  ::pwscf::pwReadFilter --
#  
# TODO: check is &SYSTEM exists, if not the input file is not pw.x input
#       (this is done by checking the existence of SYSTEM_namelist_content
#       variable !!!
# ------------------------------------------------------------------------

proc ::pwscf::pwReadFilter {moduleObj channel} {
    variable pwscf

    # clear the head & tail
    set pwscf($moduleObj,inputHeadContent) {}
    #set pwscf($moduleObj,inputTailContent) {}
    set pwscf($moduleObj,OCCUPATIONS)      {}
    set pwscf($moduleObj,HUBBARD)          {}
    set hubbard_line                       {}
    
    #
    # Check if input file is a pw.x formatted ...
    #    

    # Any properly formated pw.x input should have the following namelists:
    #    CONTROL, SYSTEM, and ELECTRONS 
    
    set namelists {CONTROL SYSTEM ELECTRONS}
    
    # Any properly formated pw.x input should have the following cards:
    #    ATOMIC_SPECIES, ATOMIC_POSITIONS, and K_POINTS
    
    set cards {ATOMIC_SPECIES ATOMIC_POSITIONS K_POINTS}

    set status [::pwscf::readFilter::findNamelistsAndCards $moduleObj $channel pw.x $namelists $cards errMsg]
    if { $status == 0 } {
	$moduleObj readFileWrongFormat pw.x $errMsg
	return $channel
    }


    #
    # check if lattice is specified by celldm() or A,B,C,...
    #
    set system 0
    while { ! [eof $channel] } {
	gets $channel _line
	if { [string match -nocase "*&SYSTEM*" $_line] } {
	    set system 1
	    continue
	}
	if { $system } {
	    if { [regexp -nocase -- $::guib::settings(NAMELIST.end_regexp) $_line] } {
		break
	    }
	    append SYSTEM_namelist_content ${_line}\n
	}
    }

    foreach record [split $SYSTEM_namelist_content ,\n] {
	set var [string trim [lindex [split $record =] 0] { \t}]
        if { $var != "" } {
            if { [string match -nocase celldm* $var] } {
                $moduleObj varset how_lattice -value celldm
            } elseif  { [string match -nocase A $var] } {
                $moduleObj varset how_lattice -value abc
            }
        }
    }

    seek $channel 0 start

    #
    # check if there are non-empty lines before starting &CONTROL
    # namelist
    # 
    while {1} {
	set res [gets $channel _line]
	if { $res == -1 && $_line == "" } {
	    # end of file occurred
	    ::tclu::ERROR "end of file occurred, while reading PW.X input"
	    return -code return
	}	
	if { [string match -nocase "*&CONTROL*" $_line] } {
	    # we found "&CONTROL" line
	    set Line(1)  $_line
	    break
	} else {
	    append pwscf($moduleObj,inputHeadContent) "$_line\n"
	}
    }
    

    # Re-order the cards in the following order:
    # --------------------------------------------------

    #   CELL_PARAMETERS
    #   ATOMIC_SPECIES
    #   ATOMIC_POSITIONS
    #   K_POINTS
    #   ADDITIONAL_K_POINTS
    #   CONSTRAINTS
    #   HUBBARD
    #   OCCUPATIONS
    #   ATOMIC_VELOCITIES
    #   ATOMIC_FORCES

    # The content of OCCUPATIONS & HUBBARD cards is managed by the "text"
    # keyword, hence we have to store the content of OCCUPATIONS & HUBBARD

    set what  {}
    set ind   1
    set _read 1
    
    while {1} {
	if { $_read } {
	    set res [gets $channel _line]
	    if { $res == -1 && $_line == "" } {
		# end of file occurred
		break
	    }    
	}
	if { [string match "CELL_PARAMETERS*" $_line] } {	    
	    set what CELL_PARAMETERS	    
	    set _line [readFilter::purifyCardLine $_line]	    
	} elseif { [string match "ATOMIC_SPECIES*" $_line] } {
	    set what ATOMIC_SPECIES	    
	    set _line [readFilter::purifyCardLine $_line]
	} elseif { [string match "ATOMIC_POSITIONS*" $_line] } {
	    set what ATOMIC_POSITIONS	    
	    set _line [readFilter::purifyCardLine $_line]
	} elseif { [string match "K_POINTS*" $_line] } {
	    set what K_POINTS
	    set _line [readFilter::purifyCardLine $_line]
	} elseif { [string match "ADDITIONAL_K_POINTS*" $_line] } {
	    set what ADDITIONAL_K_POINTS
	    set _line [readFilter::purifyCardLine $_line]
	} elseif { [string match "CONSTRAINTS*" $_line] } {
	    set what CONSTRAINTS
	    set _line [readFilter::purifyCardLine $_line]
	} elseif { [string match "HUBBARD*" $_line] } {
	    set what HUBBARD
	    set _line [readFilter::purifyCardLine $_line]
	} elseif { [string match "OCCUPATIONS*" $_line] } {
	    set what OCCUPATIONS
	    set _line [readFilter::purifyCardLine $_line]
	} elseif { [string match "ATOMIC_VELOCITIES*" $_line] } {
	    set what ATOMIC_VELOCITIES
	    set _line [readFilter::purifyCardLine $_line]
	} elseif { [string match "ATOMIC_FORCES*" $_line] } {
	    set what ATOMIC_FORCES
	    set _line [readFilter::purifyCardLine $_line]
	}

	if { $what == {} } {	    
		
	    #---------------------------------------------
	    # VARIABLES: handle multiple flags options
	    #---------------------------------------------
	    foreach {var optList} {
		smearing {
		    {'gaussian' 'gauss'}
		    {'methfessel-paxton' 'm-p' 'mp'}
		    {'marzari-vanderbilt' 'cold' 'm-v' 'mv'}
		    {'fermi-dirac' 'f-d' 'fd'}
		}
		assume_isolated {
		    {'makov-payne' 'm-p' 'mp'}
		    {'martyna-tuckerman' 'm-t' 'mt'}
		}
		vdw_corr {
		    {'grimme-d2' 'Grimme-D2' 'DFT-D'  'dft-d'}
		    {'grimme-d3' 'Grimme-D3' 'DFT-D3' 'dft-d3'}
		    {'TS'        'ts'        'ts-vdw' 'ts-vdW' 'tkatchenko-scheffler'}
                    {'MBD'       'mbd'       'many-body-dispersion' 'mbd_vdw'}
		    {'XDM'       'xdm' }
		}
                diagonalization {
                    'david'
                    'cg'
                    'ppcg'
                    {'paro' 'ParO'}
                    {'rmm-davidson' 'rmm-paro'}
                }                    
	    } {		
		set _line [readFilter::replaceVarFlag $_line $var $optList]
	    }
    
	    # logical VARIABLES: use only .true. and .false.
	    #-----------------------------------------------
	    set _line [readFilter::logicalFlag $_line]

	    incr ind
	    set Line($ind) $_line
	} else {
	    # fortranreal --> real translation
	    regsub -all -nocase {([0-9]|[0-9].)(d)([+-]?[0-9]+)} $_line {\1e\3} _transline

            # the OCCUPATIONS & HUBBARD are treated specially (see below)            
            if { [string match "HUBBARD*" $_line] } {
                set hubbard_line $_line
            } elseif { ! [string match "OCCUPATIONS*" $_line] } {
		append $what "$_transline\n"
	    }
	}
	set _read 1
    }
    
    # close the old channel
    close $channel

    # open a new channel (i.e. temporary file)     
    set tmpfile    [::tclu::tempFile name pw_input]
    set newChannel [open $tmpfile w+]

    #
    # write the file:
    # ---------------
    # write the namelists
    for {set i 1} {$i <= $ind} {incr i} {
	puts $newChannel $Line($i)
    }
    # write the CELL_PARAMETERS
    if { [info exists CELL_PARAMETERS] } {
	puts $newChannel $CELL_PARAMETERS
    }
    # write the ATOMIC_SPECIES
    if { [info exists ATOMIC_SPECIES] } {
	puts $newChannel $ATOMIC_SPECIES
    }
    # write the ATOMIC_POSITIONS
    if { [info exists ATOMIC_POSITIONS] } {
	puts $newChannel $ATOMIC_POSITIONS
    }    
    # write the K_POINTS
    if { [info exists K_POINTS] } {
	puts $newChannel $K_POINTS
    }
    # write the ADDITIONAL_K_POINTS
    if { [info exists ADDITIONAL_K_POINTS] } {
        $moduleObj varset specify_add_kpoints -value .true.
	puts $newChannel $ADDITIONAL_K_POINTS
    }
    # write the CONSTRAINTS
    if { [info exists CONSTRAINTS] } {
        $moduleObj varset constraints_enable -value Yes
	puts $newChannel $CONSTRAINTS
    }
    # write the ATOMIC_VELOCITIES
    if { [info exists ATOMIC_VELOCITIES] } {
	puts $newChannel $ATOMIC_VELOCITIES
    }
    # write the ATOMIC_FORCES
    if { [info exists ATOMIC_FORCES] } {
	$moduleObj varset specify_atomic_forces -value .true.
	puts $newChannel $ATOMIC_FORCES
    }

    # store the HUBBARD record
    if { [info exists HUBBARD] } {
	puts $newChannel $hubbard_line\n
        $moduleObj varset hubbard_enable -value Yes
	set pwscf($moduleObj,HUBBARD) [string trim $HUBBARD \n]
    }

    # store the OCCUPATIONS record
    if { [info exists OCCUPATIONS] } {
	puts $newChannel "OCCUPATIONS\n"
	set pwscf($moduleObj,OCCUPATIONS) [string trim $OCCUPATIONS \n]
    }
    
    flush $newChannel

    # rewind the newChannel
    seek $newChannel 0 start

    return $newChannel
}


# ------------------------------------------------------------------------
#  ::pwscf::pwWriteFilter --
# ------------------------------------------------------------------------

proc ::pwscf::pwWriteFilter {moduleObj outputContent} {
    variable pwscf
    # HEAD
    if { [info exists pwscf($moduleObj,inputHeadContent)] } {
	append output $pwscf($moduleObj,inputHeadContent)
    }
    # BODY
    append output $outputContent
    # TAIL
    #if { [info exists pwscf($moduleObj,inputTailContent)] } {
    #	append output $pwscf($moduleObj,inputTailContent)
    #}
    return $output    
}


# ------------------------------------------------------------------------
# Reads next non-empty line from channel. If EOF occurs, then it
# returns return-code.
# ------------------------------------------------------------------------
proc ::pwscf::_getsNonEmptyLine {channel} {
    # while loop for skipping empty lines !!!
    while {1} {
	set res [gets $channel _line]
	if { $res == -1 && $_line == "" } {
	    # end of file occurred
	    return -code return
	}    
	if { [llength $_line] != 0 } {
	    return $_line
	}
    }
}

