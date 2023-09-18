
# ------------------------------------------------------------------------
#  Page: CONTROL
# ------------------------------------------------------------------------

tracevar calculation w {

    set nat    [varvalue nat]
    set calc   [varvalue calculation]    

    set ion_dynamics [varvalue ion_dynamics]
    set fcp_dynamics [varvalue fcp_dynamics]
    set widget       [getWidgetFromVarident ion_dynamics]
    
    set all {ions cell fcp vc_md constraints_card}
    
    set disable {}
    set enable  {}

    switch -exact -- $calc {
	'scf' - 
	'nscf' -
        'bands' {
	    set disable $all
	    varset ion_dynamics -value {}
	}
        
	'relax' {
	    set enable  {ions fcp constraints_card}
	    set disable {cell vc_md for_fcp_dynamics}
	    
	    widget ion_dynamics enable
	    widgetconfigure ion_dynamics -textvalues {
		"BFGS quasi-newton method for structural optimization  <bfgs>"
		"damped dynamics (quick-min Verlet) for structural optimization  <damp>"
                "FIRE minimization algorithm  <fire>"
	    }
	    if { ! [regexp bfgs|damp|fire $ion_dynamics] } {
		varset ion_dynamics -value {}
	    }

            widgetconfigure fcp_dynamics -textvalues {
                "BFGS quasi-newton algorithm  <bfgs>"
                "Newton-Raphson algorithm with DIIS  <newton>"
                "damped (quick-min Verlet) dynamics  <damp>"
                "line-minimization algorithm  <lm>"                
            }
	    if { ! [regexp bfgs|newton|damp|lm $fcp_dynamics] } {
		varset fcp_dynamics -value {}
	    }
	}
        
	'vc-relax' {
	    set enable  {ions cell vc_md constraints_card}
            set disable {fcp}

	    widget ion_dynamics enable
	    widgetconfigure ion_dynamics -textvalues {
		"BFGS quasi-newton method for structural optimization  <bfgs>"
		"Beeman algorithm for variable cell damped dynamics  <damp>"
	    }
	    if { ! [regexp bfgs|damp $ion_dynamics] } {
		varset ion_dynamics -value {}
	    }

	    widget cell_dynamics enable
	    widgetconfigure cell_dynamics -textvalues {
		"None  <none>"
		"Damped (Beeman) dynamics of the Parrinello-Raman extended lagrangian  <damp-pr>"
		"Damped (Beeman) dynamics of the new Wentzcovitch extended lagrangian  <damp-w>"
		"BFGS quasi-newton algorithm (ion_dynamics must be 'bfgs' too)  <bfgs>"
	    }
	}
        
	'md' {
	    set enable  {ions fcp constraints_card}
	    set disable {cell vc_md}	
	    
	    widget ion_dynamics enable
	    widgetconfigure ion_dynamics -textvalues {
		"Verlet algorithm for molecular dynamics  <verlet>"
		"over-damped Langevin dynamics  <langevin>"
		"over-damped Langevin with Smart Monte Carlo  <langevin-smc>"
	    }	    
	    if { ! [regexp verlet|langevin $ion_dynamics] } {
		varset ion_dynamics -value {}
	    }

            widgetconfigure fcp_dynamics -textvalues {
                "Velocity-Verlet algorithm  <velocity-verlet>"
                "Verlet algorithm  <verlet>"                
            }
            if { ! [regexp verlet $fcp_dynamics] } {
		varset fcp_dynamics -value {}
	    }
	}
        
	'vc-md' {
	    set enable  {ions cell vc_md constraints_card}
            set disable {fcp}

	    widget ion_dynamics enable
	    widgetconfigure ion_dynamics -textvalues {
		"Beeman algorithm for variable cell MD  <beeman>"
	    }

	    if { ! [regexp beeman $ion_dynamics] } {
		varset ion_dynamics -value {}
	    }

	    #widgetconfigure ion_temperature -textvalues {
	    #    "velocity rescaling via tempw&tolp  <rescaling>"
	    #    "not controlled  <not_controlled>"
	    #}
            #
	    #if { ! [regexp rescaling [varvalue ion_temperature]] } {
	    #    varset ion_temperature -value {}
	    #}

	    widget cell_dynamics enable
	    widgetconfigure cell_dynamics -textvalues {
		"None  <none>"
		"(Beeman) dynamics of the Parrinello-Raman extended lagrangian  <pr>"
		"(Beeman) dynamics of the new Wentzcovitch extended lagrangian  <w>"
	    }
	}
    }

    foreach group $enable {
        if { $group != "fcp" } {
            groupwidget $group enable
        } else {
            if { [vartextvalue lfcp] == "Yes" } {
                groupwidget $group enable
            }
        }
    }
    foreach group $disable {
	groupwidget $group disable
    }
    
    # force to update the state of widgets by resetting corresponding variables

    varset ion_dynamics       -value [varvalue ion_dynamics]
    #varset ion_temperature    -value [varvalue ion_temperature]
    varset cell_dynamics      -value [varvalue cell_dynamics]
    varset constraints_enable -value [varvalue constraints_enable]

    widgetconfigure atomic_coordinates -caption "Enter atomic coordinates:"    	
}

tracevar lfcp w {
    if { [vartextvalue lfcp] == "Yes" } {
	groupwidget fcp enable
        if { [varvalue calculation] != "'md'" } {
            groupwidget for_fcp_dynamics disable
        }
    } else {
	groupwidget fcp disable
    }
}

tracevar trism w {
    if { [vartextvalue trism] == "Yes" } {
	groupwidget rism enable
        groupwidget solvents_card enable
    } else {
       	groupwidget rism disable
        groupwidget solvents_card disable
    } 
}
  
tracevar gate w {
    if { [vartextvalue gate] == "Yes" } {
	groupwidget gate_group enable
    } else {
	groupwidget gate_group disable
    }
}

tracevar twochem w {
    if { [vartextvalue twochem] == "Yes" } {
	groupwidget twochem_group enable
    } else {
	groupwidget twochem_group disable
    }
}


# ------------------------------------------------------------------------
#  Page: SYSTEM
# ------------------------------------------------------------------------

tracevar how_lattice w {
    switch -exact -- [varvalue how_lattice] {
	celldm { 
	    widget celldm enable 
	    groupwidget abc    disable 
	    groupwidget cosABC disable 
	}
	abc {
	    widget celldm disable 
	    groupwidget abc    enable
	    groupwidget cosABC enable
	}
    }
    varset ibrav -value [varvalue ibrav]
}	    

tracevar ibrav w {
    switch -exact -- [vartextvalue ibrav] {
	"Free lattice"             -
	"Cubic P (sc)"             -
	"Cubic F (fcc)"            -
	"Cubic I (bcc)"           { set uses {1 A} }

	"Hexagonal and Trigonal P" -
	"Tetragonal P (st)"        -       
	"Tetragonal I (bct)"      { set uses {1 3   A C} }
	"Trigonal R, 3fold axis c" -
	"Trigonal R, 3fold axis <111>" { set uses {1 4   A cosAB} }

	"Orthorhombic P"                  -
	"Orthorhombic base-centered(bco)" - 
	"Orthorhombic face-centered"      -
	"Orthorhombic body-centered"     { set uses {1 2 3   A B C} }

	"Monoclinic P, unique axis c" -
	"Monoclinic P, unique axis b" -
	"Monoclinic base-centered" { set uses {1 2 3 4   A B C cosAB} }
	"Triclinic P" { set uses { 1 2 3 4 5 6   A B C cosAB cosAC cosBC } }
	default {
	    set uses {}
	}
    }

    switch -exact -- [varvalue how_lattice] {
	celldm {
	    for {set i 1} {$i <= 6} {incr i} {
		if { [lsearch -exact $uses $i] < 0 } {		
		    # the lattice doesn't need the $i-th parameter
		    widget celldm($i) disable 
		} else {
		    widget celldm($i) enable 
		}
	    }
	}
	abc {
	    foreach par {A B C  cosAB cosAC cosBC} {
		if { [lsearch -exact $uses $par] < 0 } {
		    # the lattice doesn't need the $par parameter
		    widget $par disable
		} else {
		    widget $par enable
		}
	    }
	}
    }

    if { [vartextvalue ibrav] == "Free lattice" } {
	groupwidget cards__CELL_PARAMETERS enable
    } else {
	groupwidget cards__CELL_PARAMETERS disable
    }
}

tracevar nat w {
    set nat [varvalue nat]
    widgetconfigure atomic_coordinates  -rows $nat
    widgetconfigure atomic_forces       -rows $nat
    widgetconfigure atomic_velocities   -rows $nat    
    varset specify_atomic_forces -value [varvalue specify_atomic_forces]
    varset ion_velocities -value [varvalue ion_velocities]
}

tracevar ntyp w {
    set ntyp [varvalue ntyp]
    widgetconfigure atomic_species -rows $ntyp

    widgetconfigure starting_charge -end $ntyp
    widgetconfigure starting_magnetization -end $ntyp; # nspin-dependent
    varset nspin -value [varvalue nspin]

    widgetconfigure angle1 -end $ntyp; # noncolin
    widgetconfigure angle2 -end $ntyp
    varset noncolin -value [varvalue noncolin]
    
    widgetconfigure Hubbard_occ   -rows $ntyp
    widgetconfigure Hubbard_alpha -end $ntyp
    widgetconfigure Hubbard_beta  -end $ntyp

    widgetconfigure london_c6 -end $ntyp
    widgetconfigure london_rvdw -end $ntyp
    
    # hack
    varset london -value [varvalue london]
}


tracevar nspin w {
    if { [vartextvalue nspin] == "Yes" } {
	groupwidget spin_polarization enable
	widgetconfigure starting_magnetization -end [varvalue ntyp]

	varset noncolin -value ""
	groupwidget noncolin_group disable

	foreach var {fixed_magnetization B_field} {
	    widgetconfigure $var -start 3 -end 3
	}    
    } else {
	if { [vartextvalue noncolin] != "Yes" } {
	    groupwidget spin_polarization disable
	}
	groupwidget noncolin_group disable
    }

    # constrained/fixed magnetization
    if { [vartextvalue nspin] == "Yes" || [vartextvalue noncolin] == "Yes" } {
	groupwidget constrained_magnetization_group enable
    } else {
	groupwidget constrained_magnetization_group disable
    }
}


tracevar noncolin w {
    if { [vartextvalue noncolin] == "Yes" } {
	groupwidget spin_polarization enable
	widgetconfigure starting_magnetization -end [varvalue ntyp]

	varset nspin -value ""
	groupwidget noncolin_group enable
	
	foreach var {fixed_magnetization B_field} {
	    widgetconfigure $var -start 1 -end 3
	}
    } else {
	if { [vartextvalue nspin] != "Yes" } {
	    groupwidget spin_polarization disable
	}
	groupwidget noncolin_group disable
    }

    # constrained/fixed magnetization
    if { [vartextvalue nspin] == "Yes" || [vartextvalue noncolin] == "Yes" } {
	groupwidget constrained_magnetization_group enable
    } else {
	groupwidget constrained_magnetization_group disable
    }
}

tracevar tefield w {
    switch -- [vartextvalue tefield] {
	Yes     { groupwidget tefield_group enable  }
	default { groupwidget tefield_group disable }
    }
}

tracevar lelfield w {
    switch -- [vartextvalue lelfield] {
	Yes     { foreach w {nberrycyc efield efield_cart} {widget $w enable}; groupwidget elfield_group enable }
	default { foreach w {nberrycyc efield efield_cart} {widget $w disable};groupwidget elfield_group disable  }
    }
}

tracevar occupations w {
    if { [varvalue occupations] == "'from_input'" } {
	groupwidget occupations_card enable
    } else {
	groupwidget occupations_card disable
    }
}


tracevar assume_isolated w {
    if { [varvalue assume_isolated] == "'esm'" } {
	groupwidget ESM enable
    } else {
	groupwidget ESM disable
    }
}

tracevar london w {
    if { [vartextvalue london] == "Yes" } {
	groupwidget dftdG enable
	widgetconfigure london_c6   -end [varvalue ntyp]
	widgetconfigure london_rvdw -end [varvalue ntyp]
	
	varset xdm -value .false.
	groupwidget xdmG  disable 
	groupwidget tsG   disable	
    } else  {
	groupwidget dftdG disable 
    }
}

tracevar xdm w {
    if { [vartextvalue xdm] == "Yes" } { 
	groupwidget xdmG enable
	
	varset london -value .false.
	groupwidget dftdG disable 
	groupwidget tsG   disable
    } else  {
	groupwidget xdmG disable 
    }
}

tracevar vdw_corr w {
    groupwidget dftdG  disable
    groupwidget dftd3G disable
    groupwidget xdmG   disable 
    groupwidget tsG    disable
    
    if { [varvalue vdw_corr] == "'grimme-d2'" } {
	groupwidget dftdG enable 
    } elseif { [varvalue vdw_corr] == "'grimme-d3'" } {
	groupwidget dftd3G enable 
    } elseif { [varvalue vdw_corr] == "'XDM'" } {
	groupwidget xdmG  enable
    } elseif { [varvalue vdw_corr] == "'TS'" } {
	groupwidget tsG   enable
    } elseif { [varvalue vdw_corr] == "'MBD'" } {
        widget ts_vdw_isolated enable
    }
}

tracevar lgcscf w {
    if { [vartextvalue lgcscf] == "Yes" } {
	groupwidget gcscf_group enable
    } else {
	groupwidget gcscf_group disable
    }
}

# ------------------------------------------------------------------------
#  Page: ELECTRONS
# ------------------------------------------------------------------------

tracevar adaptive_thr w {
    if { [vartextvalue adaptive_thr] == "Yes" } { 
	groupwidget adaptive_thr_setup enable 
    } else  {
	groupwidget adaptive_thr_setup disable
    }
}


tracevar diagonalization w {
    switch -glob -- [varvalue diagonalization] {
	'david*' {
	    widget diago_david_ndim enable
	    widget diago_cg_maxiter disable
	    widget diago_ppcg_maxiter disable
            widget diago_rmm_ndim disable
            widget diago_rmm_conv disable
	}
	'cg' {
	    widget diago_cg_maxiter enable
	    widget diago_ppcg_maxiter disable
	    widget diago_david_ndim disable
            widget diago_rmm_ndim disable
            widget diago_rmm_conv disable
	}
        'ppcg' {
            widget diago_ppcg_maxiter enable
            widget diago_cg_maxiter disable
	    widget diago_david_ndim disable
            widget diago_rmm_ndim disable
            widget diago_rmm_conv disable
        }
        'rmm*' {
            widget diago_david_ndim enable
            widget diago_rmm_ndim enable
            widget diago_rmm_conv enable
            widget diago_cg_maxiter disable
            widget diago_ppcg_maxiter disable
        }
	default {
            # 'paro'
	    widget diago_cg_maxiter disable
	    widget diago_david_ndim disable
	}
    }
}


# ------------------------------------------------------------------------
#  Page: IONS
# ------------------------------------------------------------------------

tracevar ion_dynamics w {
    set calc [varvalue calculation]
    set iond [varvalue ion_dynamics]

    # MD
    switch -exact -- $calc {
	'scf' - 'nscf' - 'bands' {
	    groupwidget md disable
	}
	'relax' - 'vc-relax' {	    
	    # check !!!
	    switch -exact -- $iond {
		'damp' - 'verlet' - 'langevin' - 'langevin-smc' - 'beeman' {
		    groupwidget md enable
		} 
		default {
		    groupwidget md disable
		}
	    }
	}
	'md' - 'vc-md' {	    
	    # check !!!
	    groupwidget md enable
	}
    }

    # BFGS
    if { $iond == "'bfgs'" && ($calc == "'relax'" || $calc == "'vc-relax'") } {
	groupwidget bfgs enable
    } else {
	groupwidget bfgs disable
    }

    # FIRE
    if { $iond == "'fire'" && $calc == "'relax'" } {
	groupwidget fire enable
    } else {
	groupwidget fire disable
    }
}

tracevar n_fe_step w {
    widgetconfigure fe_step -end [varvalue n_fe_step]
}

# ------------------------------------------------------------------------
# Page: CELL_PARAMETERS, ATOMIC_SPECIES, ATOMIC_POSITIONS
# ------------------------------------------------------------------------


#
# ------------------------------------------------------------------------
#  Page: K-POINT DATA
# ------------------------------------------------------------------------

tracevar K_POINTS_flags w {
    switch -exact -- [varvalue K_POINTS_flags] {
	tpiba -	crystal - crystal_b - tpiba_b - {} {
	    groupwidget nks_line   enable
	    groupwidget kmesh_line disable
	    widget kpoints enable
	    widgetconfigure kpoints -cols 4 -rows [varvalue nks]    
	}

	automatic {
	    groupwidget nks_line   disable
	    groupwidget kmesh_line enable
	    widget kpoints disable
	    widgetconfigure kpoints -cols 4 -rows 0
	}

	gamma {
	    groupwidget nks_line   disable
	    groupwidget kmesh_line disable
	    widget kpoints disable
	    widgetconfigure kpoints -cols 4 -rows 0
	}

	crystal_c - tpiba_c - {} {
	    varset nks -value 3
	    groupwidget nks_line   enable
	    groupwidget kmesh_line disable
	    widget kpoints enable
	    widgetconfigure kpoints -cols 4 -rows [varvalue nks]  
	}
    }
}

tracevar nks w {
    widgetconfigure kpoints -rows [varvalue nks]    
}
tracevar nks_add w {
    widgetconfigure add_kpoints -rows [varvalue nks_add]    
}


# ------------------------------------------------------------------------
# Page: Other cards
# ------------------------------------------------------------------------

tracevar hubbard_enable w {
    if { [varvalue hubbard_enable] == "Yes" } {
	groupwidget hubbard_group enable
    } else {
	groupwidget hubbard_group disable
    }
}

tracevar constraints_enable w {

    set calc [varvalue calculation]

    if { [regexp relax|md $calc] } {
	
	widget constraints_enable enable

	if { [varvalue constraints_enable] } {
	    groupwidget constraints_card enable
	} else {
	    groupwidget constraints_card disable
	}	    	
    } else {
	widget constraints_enable disable
	groupwidget constraints_card disable		
    }
}

tracevar nconstr w {
    set nc [varvalue nconstr]    
    widgetconfigure constraints_table -rows $nc
}

tracevar specify_atomic_forces w {
    if { [vartextvalue specify_atomic_forces] == "Yes" } {
	groupwidget atomic_forces_specs enable
    } else {
	groupwidget atomic_forces_specs disable
    }
}
tracevar specify_add_kpoints w {
    if { [vartextvalue specify_add_kpoints] == "Yes" } {
	groupwidget add_kpoints_specs enable
    } else {
	groupwidget add_kpoints_specs disable
    }
}
tracevar ion_velocities w {
    if { [varvalue ion_velocities] == "'from_input'" } {
	groupwidget atomic_velocities_group enable
    } else {
	groupwidget atomic_velocities_group disable
    }
}

tracevar laue_both_hands w {
    if { [vartextvalue laue_both_hands] == "Yes" } {
        widget laue_one_hand_table forget
        widget laue_both_hands_table create
        widget laue_both_hands_table enable
    } else {
        widget laue_one_hand_table create
        widget laue_one_hand_table enable
        widget laue_both_hands_table forget
    }
}

tracevar nsolv w {
    set nsolv [varvalue nsolv]
    if { $nsolv == "" || ! [string is integer $nsolv] || $nsolv < 1 } {
        set nsolv 1
    }
    widgetconfigure laue_one_hand_table -rows $nsolv
    widgetconfigure laue_both_hands_table -rows $nsolv
}

# ------------------------------------------------------------------------
# POST-PROCESSING: assign default values for "traced" variables, ...
# ------------------------------------------------------------------------
postprocess {    
    varset calculation     -value 'scf'
    varset lfcp            -value {}
    varset trism           -value {}
    varset laue_both_hands -value {}
    varset twochem         -value {}
    varset gate            -value {}
    varset ibrav           -value {}
    varset how_lattice     -value celldm
    varset nspin           -value {}
    varset tefield         -value {}
    varset occupations     -value {}
    varset assume_isolated -value {}
    varset vdw_corr        -value {}
    varset london          -value {}
    varset lgcscf          -value {}
    varset xdm             -value {}
    varset adaptive_thr    -value {}
    varset diagonalization -value {}
    varset ion_dynamics    -value {}
    varset ion_velocities  -value {}
    varset K_POINTS_flags  -value automatic
    varset CELL_PARAMETERS_flags -value {}
    varset hubbard_enable -value No
    varset HUBBARD_flags  -value ortho-atomic
    # unused variables
    groupwidget unused_1 disable
    #groupwidget vdw_obsolete disable
    
    varset specify_atomic_forces -value .false.
    varset specify_add_kpoints   -value .false.
}
