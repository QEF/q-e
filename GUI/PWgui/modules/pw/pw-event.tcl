
# ------------------------------------------------------------------------
#  Page: CONTROL
# ------------------------------------------------------------------------

tracevar calculation w {

    set calc   [varvalue calculation]
    set widget [getWidgetFromVarident ion_dynamics]
     switch -exact -- $calc {
	 'scf' - 
	 'nscf' {
	     groupwidget ions   disable
	     groupwidget cell   disable
	     groupwidget phonon disable
	     groupwidget vc_md  disable
	 }
	 'phonon' {
	     groupwidget ions   disable
	     groupwidget cell   disable
	     groupwidget phonon enable
	     groupwidget vc_md  disable
	 }
	 'relax' {
	     groupwidget ions enable
	     widget ion_dynamics enable
	     widgetconfigure ion_dynamics -textvalues {
		 "new BFGS algorithm, based on the trust radius procedure <bfgs>"
		 "old BFGS quasi-newton method for structural relaxation  <old-bfgs>"
		 "BFGS with the CONSTRAINT coded into \"constraint.f90\"  <constrained-bfgs>"
		 "damped (Beeman) dynamics for structural relaxation  <damp>"
	     }
	     groupwidget neb    disable
	     groupwidget cell   disable
	     groupwidget phonon disable
	     groupwidget vc_md  disable
	 }
	 'vc-relax' {
	     groupwidget ions enable
	     widget ion_dynamics enable
	     widgetconfigure ion_dynamics -textvalues {
		 "damped (Beeman) dynamics for structural relaxation  <damp>"
	     }
	     groupwidget neb    disable
	     groupwidget cell   enable
	     groupwidget phonon disable	  
	     groupwidget vc_md  enable
	 }
	 'md' {
	     groupwidget ions enable
	     widget ion_dynamics enable
	     widgetconfigure ion_dynamics -textvalues {
		 "Verlet algorithm for Molecular dynamics  <verlet>"
		 "Verlet-MD with the CONSTRAINT coded into \"constraint.f90\"  <constrained-verlet>"
		 "Beeman algorithm for MD  <beeman>"
	     }
	     groupwidget neb    disable
	     groupwidget cell   disable
	     groupwidget phonon disable
	     groupwidget vc_md  disable
	 }
	 'vc-md' {
	     groupwidget ions enable
	     widget ion_dynamics enable
	     widgetconfigure ion_dynamics -textvalues {
		 "Beeman algorithm for MD  <beeman>"
	     }
	     groupwidget neb    disable
	     groupwidget cell   enable
	     groupwidget phonon disable	     
	     groupwidget vc_md  enable
	 }
	 'neb' {
	     groupwidget ions     enable
	     widget ion_dynamics  disable
	     groupwidget neb      enable
	     groupwidget cell     disable
	     groupwidget phonon   disable
	     groupwidget vc_md    disable
	 }
     }

    # force the update of the state of new_bfgs group widgets
    varset ion_dynamics -value [varvalue ion_dynamics]

    # force the update of the state of the CLIMBING_IMAGES card
    varset CI_scheme -value [varvalue CI_scheme]

    # take care of NEB coordinates
    if { $calc == "'neb'" } {
	keywordconfigure first_image enable
	keywordconfigure last_image  enable
	widget atomic_coordinates2   create

	widgetconfigure atomic_coordinates -caption "Enter atomic coordinates for 1st image:"
    } else {
	keywordconfigure first_image disable
	keywordconfigure last_image  disable
	widget atomic_coordinates2   forget

	widgetconfigure atomic_coordinates -caption "Enter atomic coordinates:"
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
}	    

tracevar ibrav w {
    switch -exact -- [vartextvalue ibrav] {
	"Free lattice"             -
	"Cubic P (sc)"             -
	"Cubic F (fcc)"            -
	"Cubic I (bcc)"           { set uses 1 }
	
	"Hexagonal and Trigonal P" -
	"Tetragonal P (st)"        -       
	"Tetragonal I (bct)"      { set uses {1 3} }
	"Trigonal R"              { set uses {1 4} }
	
	"Orthorhombic P"                  -
	"Orthorhombic base-centered(bco)" - 
	"Orthorhombic face-centered"      -
	"Orthorhombic body-centered"     { set uses {1 2 3} }
	
	"Monoclinic P"              -
	"Monoclinic base-centered" { set uses {1 2 3 4} }
	"Triclinic P" { set uses { 1 2 3 4 5 6 } }
	default {
	    set uses {}
	}
    }
    for {set i 1} {$i <= 6} {incr i} {
	if { [lsearch -exact $uses $i] < 0 } {		
	    
	    # the lattice doesn't need the $i-th parameter
	    
	    widget celldm($i) disable 
	} else {
	    widget celldm($i) enable 
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
    widgetconfigure atomic_coordinates2 -rows $nat
}
tracevar ntyp w {
    set ntyp [varvalue ntyp]
    widgetconfigure atomic_species -rows $ntyp
    widgetconfigure Hubbard_U     -end $ntyp
    widgetconfigure Hubbard_alpha -end $ntyp
}

tracevar nspin w {
    if { [vartextvalue nspin] == "Yes" } {
	widget starting_magnetization enable
	widgetconfigure starting_magnetization -end [varvalue nat]
    } else {
	widget starting_magnetization disable
    }
}

tracevar tefield w {
    switch -- [vartextvalue tefield] {
	Yes     { groupwidget tefield_group enable  }
	default { groupwidget tefield_group disable }
    }
}
	    
tracevar lda_plus_u w {
    switch -- [vartextvalue lda_plus_u] {
	Yes     { widget mixing_fixed_ns enable;  groupwidget hubbard enable }
	default { widget mixing_fixed_ns disable; groupwidget hubbard disable }
    }
}

tracevar occupations w {
    if { [varvalue occupations] == "'from_input'" } {
	groupwidget occupations_card enable
    } else {
	groupwidget occupations_card disable
    }
}


# ------------------------------------------------------------------------
#  Page: ELECTRONS
# ------------------------------------------------------------------------

tracevar diagonalization w {
    switch -glob -- [varvalue diagonalization] {
	'david*' {
	    widget diago_cg_maxiter disable
	    widget diago_david_ndim enable
	    groupwidget diis        disable
	}
	'diis*' {
	    widget diago_cg_maxiter disable
	    widget diago_david_ndim disable
	    groupwidget diis        enable
	}
	'cg' {
	    widget diago_cg_maxiter enable
	    widget diago_david_ndim disable
	    groupwidget diis        disable
	}
	default {
	    widget diago_cg_maxiter disable
	    widget diago_david_ndim disable
	    groupwidget diis        disable
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
    	'scf' - 'nscf' - 'phonon' - 'neb' {
    	    groupwidget md disable
    	}
    	'relax' - 'vc-relax' - 'md' - 'vc-md' {	    
    	    switch -exact -- $iond {
    	    	'damp' - 'verlet' - 'constrained-verlet' - 'beeman' {
    	    	    groupwidget md enable
    	    	} 
    	    	default {
    	    	    groupwidget md disable
    	    	}
    	    }
    	}
    }

    # NEW-BFGS
    if { $iond == "'bfgs'" && $calc == "'relax'" } {
	groupwidget new_bfgs enable
    } else {
	groupwidget new_bfgs disable
    }
}

tracevar CI_scheme w {
    if { [varvalue calculation] == "'neb'" && [varvalue CI_scheme] == "'manual'" } {
	groupwidget climbing_images enable
    } else {
	groupwidget climbing_images disable
    }
}
	

# ------------------------------------------------------------------------
#  Page: K-POINT DATA
# ------------------------------------------------------------------------

tracevar kpoint_type w {
    switch -exact -- [varvalue kpoint_type] {
    	tpiba -	crystal {
    	    #widget nks enable
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
    }
}

tracevar nks w {
    widgetconfigure kpoints -rows [varvalue nks]    
}


# ------------------------------------------------------------------------
# Page: Other cards
# ------------------------------------------------------------------------

# incoming in next PWscf version
#tracevar nconstr w {
#    widgetconfigure constraints_table -rows [varvalue nconstr]
#}


# ------------------------------------------------------------------------
# POST-PROCESSING: assign default values for "traced" variables
# ------------------------------------------------------------------------
postprocess {
    varset calculation     -value 'scf'
    varset how_lattice     -value celldm
    varset ibrav           -value {}
    varset nspin           -value {}
    varset tefield         -value {}
    varset lda_plus_u      -value {}
    varset occupations     -value {}
    varset diagonalization -value {}
    varset CI_scheme       -value {}
    varset kpoint_type     -value tpiba
    varset lattice_type    -value cubic
    varset ion_dynamics    -value {}
}
