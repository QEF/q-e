
# ------------------------------------------------------------------------
#  Page: CONTROL
# ------------------------------------------------------------------------

tracevar calculation w {

    set nat    [varvalue nat]
    set calc   [varvalue calculation]
    set widget [getWidgetFromVarident ion_dynamics]
    
    groupwidget ions   disable
    groupwidget cell   disable
    groupwidget phonon disable
    groupwidget vc_md  disable
    groupwidget path   disable
    groupwidget neb    disable
    #groupwidget smd    disable
    groupwidget metadyn disable

    switch -exact -- $calc {
	'scf' - 
	'nscf' {
	    groupwidget ions   disable
	    groupwidget cell   disable
	    groupwidget phonon disable
	    groupwidget vc_md  disable
	    groupwidget path   disable
	    groupwidget neb    disable
            #groupwidget smd    disable
	    groupwidget metadyn disable
	}
	'phonon' {
	    groupwidget ions   disable
	    groupwidget cell   disable
	    groupwidget phonon enable
	    groupwidget vc_md  disable
	    groupwidget path   disable
	    groupwidget neb    disable
            #groupwidget smd    disable
	    groupwidget metadyn disable
	}
	'relax' {
	    groupwidget ions   enable
	    groupwidget cell   disable
	    groupwidget phonon disable
	    groupwidget vc_md  disable
	    groupwidget path   disable
	    groupwidget neb    disable
            #groupwidget smd    disable
	    groupwidget metadyn disable	    
	    widget ion_dynamics enable
	    widgetconfigure ion_dynamics -textvalues {
		"BFGS quasi-newton method for structural optimization <bfgs>"
		"damped dynamics (quick-min Verlet) for structural optimization <damp>"
		"damped dynamics (quick-min Verlet) for structural optimization with the CONSTRAINT <constrained-damp>"
	    }
	}
	'vc-relax' {
	    groupwidget ions   enable
	    groupwidget cell   enable
	    groupwidget phonon disable
	    groupwidget vc_md  enable
	    groupwidget path   disable
	    groupwidget neb    disable
            #groupwidget smd    disable
	    groupwidget metadyn disable	    
	    widget ion_dynamics enable
	    widgetconfigure ion_dynamics -textvalues {
		"Beeman algorithm for variable cell damped dynamics <damp>"
	    }
	}
	'md' {
	    groupwidget ions   enable
	    groupwidget cell   disable
	    groupwidget phonon disable
	    groupwidget vc_md  disable
	    groupwidget path   disable
	    groupwidget neb    disable
            #groupwidget smd    disable
	    groupwidget metadyn disable	    
	    widget ion_dynamics enable
	    widgetconfigure ion_dynamics -textvalues {
		"Verlet algorithm for molecular dynamics <verlet>"
		"Verlet algorithm for constrained molecular dynamics <constrained-verlet>"
	    }
	}
	'vc-md' {
	    groupwidget ions   enable
	    groupwidget cell   enable
	    groupwidget phonon disable
	    groupwidget vc_md  enable
	    groupwidget path   disable
	    groupwidget neb    disable
            #groupwidget smd    disable
	    groupwidget metadyn disable	    
	    widget ion_dynamics enable
	    widgetconfigure ion_dynamics -textvalues {
		"Beeman algorithm for variable cell MD <beeman>"
	    }
	}
	'neb' {
	    groupwidget ions   enable
	    groupwidget cell   disable
	    groupwidget phonon disable
	    groupwidget vc_md  disable
	    groupwidget path   enable
	    groupwidget neb    enable
            #groupwidget smd    disable
	    groupwidget metadyn disable	    
	}
	'smd' {
	    groupwidget ions   enable
	    groupwidget cell   disable
	    groupwidget phonon disable
	    groupwidget vc_md  disable
	    groupwidget path   enable
	    groupwidget neb    disable
            #groupwidget smd    enable
	    groupwidget metadyn disable	    
	}
	'metadyn' {
	    groupwidget ions   enable
	    groupwidget cell   disable
	    groupwidget phonon disable
	    groupwidget vc_md  disable
	    groupwidget path   disable
	    groupwidget neb    disable
            #groupwidget smd    disable
	    groupwidget metadyn enable
	}
    }

    # force the update of the state of bfgs group widgets
    varset ion_dynamics -value [varvalue ion_dynamics]

    # force the update of the state of the CLIMBING_IMAGES card
    varset CI_scheme -value [varvalue CI_scheme]

    # take care of NEB || SMD coordinates

    #update
    set ni [varvalue  path_inter_nimages]; if { $ni == "" } { set ni 0 }
    
    if { $calc == "'neb'" || $calc == "'smd'" } {    	
    	widget path_inter_nimages  enable
    	widgetconfigure atomic_coordinates -caption "Enter atomic coordinates for the FIRST image:"
    	
    	keywordconfigure first_image  enable
    	keywordconfigure last_image   enable
    	
	for {set i 1} {$i <= $ni} {incr i} {
    	    keywordconfigure intermediate_image_$i  enable
    	    widget           atomic_coordinates_${i}_inter  create
	    widgetconfigure  atomic_coordinates_${i}_inter -rows $nat
    	}
    	widget atomic_coordinates_last  create
    	
    } else {
	widget path_inter_nimages  disable
    	widgetconfigure atomic_coordinates -caption "Enter atomic coordinates:"
    	
    	keywordconfigure first_image  disable
    	keywordconfigure last_image   disable
		
	for {set i 1} {$i <= $ni} {incr i} {
	    keywordconfigure intermediate_image_$i  disable
	    widget   atomic_coordinates_${i}_inter  forget
	}
	widget atomic_coordinates_last  forget	
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
    widgetconfigure atomic_coordinates_last -rows $nat

    set ni [varvalue path_inter_nimages]
    if { $ni == "" } { set ni 0 }
    for {set i 1} {$i <= $ni} {incr i} {
	widgetconfigure atomic_coordinates_${i}_inter -rows $nat
    }
}
tracevar ntyp w {
    set ntyp [varvalue ntyp]
    widgetconfigure atomic_species -rows $ntyp
    widgetconfigure starting_magnetization -end $ntyp
    widgetconfigure angle1 -end $ntyp
    widgetconfigure angle2 -end $ntyp
    widgetconfigure Hubbard_U     -end $ntyp
    widgetconfigure Hubbard_alpha -end $ntyp
}

#tracevar nspin w {
#    if { [vartextvalue nspin] == "Yes" || [vartextvalue nspin] == "Yes noncollinear"} {
#	widget starting_magnetization enable
#	widgetconfigure starting_magnetization -end [varvalue ntyp]
#        if { [vartextvalue nspin] == "Yes" } {
#	    groupwidget noncolin_group disable
#	    #widget angle1 disable
#	    #widget angle2 disable
#        } else {
#	    groupwidget noncolin_group enable
#	    #widget angle1 enable
#	    #widget angle2 enable
#	    widgetconfigure angle1 -end [varvalue ntyp]
#	    widgetconfigure angle2 -end [varvalue ntyp]
#        }
#    } else {
#	widget starting_magnetization disable
#	groupwidget noncolin_group disable
#	#widget angle1 disable
#	#widget angle2 disable
#    }
#}


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
	#widget angle1 disable
	#widget angle2 disable
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
	#widget angle1 disable
	#widget angle2 disable
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
	'scf' - 'nscf' - 'phonon' - 'neb' - 'smd' {
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

    # BFGS
    if { $iond == "'bfgs'" && $calc == "'relax'" } {
	groupwidget bfgs enable
    } else {
	groupwidget bfgs disable
    }

    # CONSTRAINTS
    if { ( $iond == "'constrained-damp'" && $calc == "'relax'" ) \
	     || ( $iond == "'constrained-verlet'" && $calc == "'md'" ) } {
	groupwidget constraints_card enable
    } else {
	groupwidget constraints_card disable
    }	     
}

tracevar CI_scheme w {
    if { [varvalue calculation] == "'neb'" && [varvalue CI_scheme] == "'manual'" } {
	groupwidget climbing_images enable
    } else {
	groupwidget climbing_images disable
    }
}

tracevar n_fe_step w {
    widgetconfigure fe_step -end [varvalue n_fe_step]
}

# ------------------------------------------------------------------------
# Page: CELL_PARAMETERS, ATOMIC_SPECIES, ATOMIC_POSITIONS
# ------------------------------------------------------------------------

tracevar path_inter_nimages w {
    # Note: this is a bit complicated ...
    

    set nat     [varvalue nat]
    set ni      [varvalue path_inter_nimages] 
    set ni_old  [varvalue old_path_inter_nimages]
    if { $nat == ""    } { set nat 1 }
    if { $ni == ""     } { set ni 0 }
    if { $ni_old == "" } { set ni_old 0 }
    
    if { $ni_old > $ni } {
	# delete tables ...

	for {set i $ni_old} { $i > $ni } {incr i -1} {
	    keywordconfigure intermediate_image_$i  disable
	    widget atomic_coordinates_${i}_inter    forget
	}
    } elseif { $ni_old < $ni } {
	# create tables ...
	
	widget atomic_coordinates_last  forget ; # this forces the right pack-order

	for {set i 1} {$i <= $ni} {incr i} {
	    keywordconfigure intermediate_image_$i  enable
	    widget           atomic_coordinates_${i}_inter  create
	    widgetconfigure  atomic_coordinates_${i}_inter -rows $nat
	}

    	widget atomic_coordinates_last  create
    }
	

    # remember current value of path_inter_nimages
    varset old_path_inter_nimages -value $ni
}


#
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
tracevar nconstr w {
    set nc [varvalue nconstr]

    widgetconfigure constraints_table -rows $nc

    #for {set i 1} {$i <= $nc} {incr i} {
    #	tableset constraints_table $i 1 -value 1
    #}
}


# ------------------------------------------------------------------------
# POST-PROCESSING: assign default values for "traced" variables, ...
# ------------------------------------------------------------------------
postprocess {
    # BEWARE: it is assumed that 50 intermediate images is the
    # largest allowed number (this is dirty)
    varset old_path_inter_nimages -value 50
    varset path_inter_nimages     -value 0
    
    varset calculation    -value 'scf'
    varset how_lattice     -value celldm
    varset ibrav           -value {}
    varset nspin           -value {}
    varset tefield         -value {}
    varset lda_plus_u      -value {}
    varset occupations     -value {}
    varset diagonalization -value {}
    varset CI_scheme       -value {}
    varset kpoint_type     -value automatic
    varset lattice_type    -value cubic
    varset ion_dynamics    -value {}

    # so far the only constraint-type is "1"
    tableset constraints_table 1 1 -value 1

    # unused variables
    groupwidget unused_1 disable
}
