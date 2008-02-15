
# ------------------------------------------------------------------------
#  Page: CONTROL
# ------------------------------------------------------------------------

tracevar calculation w {

    set nat    [varvalue nat]
    set calc   [varvalue calculation]
    set widget [getWidgetFromVarident ion_dynamics]
    
    set all {ions cell phonon vc_md path neb metadyn constraints_card collective_vars_card}
    
    #foreach group $all {
    #	groupwidget $group disable
    #}
	
    set disable {}
    set enable  {}

    widgetconfigure ion_temperature -textvalues {
	"velocity rescaling via tempw&tolp  <rescaling>"
	"velocity rescaling via tempw&nraise  <rescale-v>"
	"velocity rescaling via delta_t  <rescale-T>"
	"reduce ionic temperature via delta_t&nraise  <reduce-T>"
	"\"soft\" Berendsen velocity rescaling via tempw&nraise  <berendsen>"
	"use Andersen thermostat via tempw&nraise  <andersen>"	
	"not controlled  <not_controlled>"
    }

    switch -exact -- $calc {
	'scf' - 
	'nscf' {
	    set disable $all
	}
	'phonon' {
	    set disable $all
	}
	'relax' {
	    set enable  {ions constraints_card}
	    set disable {cell phonon vc_md path neb metadyn collective_vars_card}
	    
	    widget ion_dynamics enable
	    widgetconfigure ion_dynamics -textvalues {
		"BFGS quasi-newton method for structural optimization  <bfgs>"
		"damped dynamics (quick-min Verlet) for structural optimization  <damp>"
	    }
	}
	'vc-relax' {
	    set enable  {ions cell vc_md constraints_card}
	    set disable {phonon path neb metadyn collective_vars_card}

	    widget ion_dynamics enable
	    widgetconfigure ion_dynamics -textvalues {
		"BFGS quasi-newton method for structural optimization  <bfgs>"
		"Beeman algorithm for variable cell damped dynamics  <damp>"
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
	    set enable {ions constraints_card}
	    set disable {cell phonon vc_md path neb metadyn collective_vars_card}	
	    
	    widget ion_dynamics enable
	    widgetconfigure ion_dynamics -textvalues {
		"Verlet algorithm for molecular dynamics  <verlet>"
		"over-damped Langevin dynamics  <langevin>"
	    }
	}
	'vc-md' {
	    set enable {ions cell vc_md constraints_card}
	    set disable {phonon path neb metadyn collective_vars_card}

	    widget ion_dynamics enable
	    widgetconfigure ion_dynamics -textvalues {
		"Beeman algorithm for variable cell MD  <beeman>"
	    }

	    widgetconfigure ion_temperature -textvalues {
		"velocity rescaling via tempw&tolp  <rescaling>"
		"not controlled  <not_controlled>"
	    }

	    widget cell_dynamics enable
	    widgetconfigure cell_dynamics -textvalues {
		"None  <none>"
		"(Beeman) dynamics of the Parrinello-Raman extended lagrangian  <pr>"
		"(Beeman) dynamics of the new Wentzcovitch extended lagrangian  <w>"
	    }
	}
	'neb' {
	    set enable  {ions path neb constraints_card}
	    set disable {cell phonon vc_md metadyn collective_vars_card}

	    widget opt_scheme enable
	    widgetconfigure opt_scheme -textvalues {
		"optimization algorithm based on molecular dynamics  <quick-min>"
		"second Broyden method  <broyden>"
		"steepest descent  <sd>"
	    }
	}
	'smd' {
	    set enable  {ions path constraints_card collective_vars_card}
	    set disable {cell phonon vc_md neb metadyn}

	    widget opt_scheme enable
	    widgetconfigure opt_scheme -textvalues {
		"optimization algorithm based on molecular dynamics  <quick-min>"
		"second Broyden method  <broyden>"
		"steepest descent  <sd>"
		"finite temperature langevin dynamics  <langevin>"
	    }
	}
	'metadyn' {
	    set enable  {ions metadyn neb constraints_card collective_vars_card}
	    set disable {cell phonon vc_md path}
	}
    }

    foreach group $enable {
	groupwidget $group enable
    }
    foreach group $disable {
	groupwidget $group disable
    }
    
    # force to update the state of widgets by resetting corresponding variables

    varset ion_dynamics           -value [varvalue ion_dynamics]
    varset opt_scheme             -value [varvalue opt_scheme]
    varset CI_scheme              -value [varvalue CI_scheme]
    varset constraints_enable     -value [varvalue constraints_enable]
    varset collective_vars_enable -value [varvalue collective_vars_enable]


    # take care of NEB || SMD coordinates

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
	    widget ortho_para       enable
	    #groupwidget diis        disable
	}
	'cg' {
	    widget diago_cg_maxiter enable
	    widget diago_david_ndim disable
	    widget ortho_para       disable
	    #groupwidget diis        disable
	}
	default {
	    widget diago_cg_maxiter disable
	    widget diago_david_ndim disable
	    widget ortho_para       disable
	    #groupwidget diis        disable
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
	'relax' - 'vc-relax' - 'md' - 'vc-md' - 'metadyn' {	    
	    switch -exact -- $iond {
		'damp' - 'verlet' - 'langevin' - 'beeman' {
		    groupwidget md enable
		} 
		default {
		    groupwidget md disable
		}
	    }
	}
    }

    # BFGS
    if { $iond == "'bfgs'" && ($calc == "'relax'" || $calc == "'vc-relax'") } {
	groupwidget bfgs enable
    } else {
	groupwidget bfgs disable
    }

    ## CONSTRAINTS
    #if { ( $iond == "'constrained-damp'" && $calc == "'relax'" ) \
    #	     || ( $iond == "'constrained-verlet'" && $calc == "'md'" ) } {
    #	groupwidget constraints_card enable
    #} else {
    #	groupwidget constraints_card disable
    #}	     
}

tracevar opt_scheme w {
    if { [regexp smd [varvalue calculation]] && [regexp langevin [varvalue opt_scheme]] } {
	widget temp_req enable
    } else {
	widget temp_req disable
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
	tpiba -	crystal - {} {
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

tracevar constraints_enable w {

    set calc [varvalue calculation]

    if { [regexp relax|md|metadyn $calc] } {
	
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

    #set nc_old  [varvalue old_nconstr]
    #if { $nc_old == "" } { set nc_old 0 }
    #if { $nc == "" } { set nc 0 }
    #if { $nc < 0 } { set nc 0 }
    #if { $nc_old < 0 } { set nc_old 0 }
    #
    #if { $nc_old > $nc } {
    #	for {set i $nc_old} { $i > $nc } {incr i -1} {
    #	    puts stderr "*** i= $i nc=$nc old=$nc_old forget"
    #	    widget constraint_type.$i forget
    #	    widget constraint.$i forget
    #	    widget constr.${i}_1 forget
    #	    widget constr.${i}_2 forget
    #	    widget constr.${i}_3 forget
    #	    widget constr.${i}_4 forget
    #	    widget constr_target_$i forget
    #	}
    #} elseif { $nc_old < $nc } {
    #	for {set i $nc} {$i > $nc_old} {incr i -1} {
    #	    puts stderr "*** i= $i create"
    #	    widget constraint_type.$i create
    #	    widget constraint.$i create
    #	    widget constr.${i}_1 create
    #	    widget constr.${i}_2 create
    #	    widget constr.${i}_3 create
    #	    widget constr.${i}_4 create
    #	    widget constr_target_$i create
    #	}
    #}
    #
    #varset old_nconstr -value $nc
}

tracevar collective_vars_enable w {
    set calc [varvalue calculation]
    
    if { [regexp smd|metadyn $calc] } {
	
	widget collective_vars_enable enable

	if { $calc == "'metadyn'" } {
	    varset collective_vars_enable -value Yes
	}
	if { [varvalue collective_vars_enable] } {
	    groupwidget collective_vars_card enable
	} else {
	    groupwidget collective_vars_card disable
	}
	     
    } else {
	
	widget collective_vars_enable disable
	groupwidget collective_vars_card disable	
	
    }
}

tracevar ncolvars w {
    set nc [varvalue ncolvars]    
    widgetconfigure colvars_table -rows $nc
}

# ------------------------------------------------------------------------
# POST-PROCESSING: assign default values for "traced" variables, ...
# ------------------------------------------------------------------------
postprocess {
    # BEWARE: it is assumed that 50 intermediate images is the
    # largest allowed number (this is dirty)
    varset old_path_inter_nimages -value 50
    varset path_inter_nimages     -value 0

    #varset old_nconstr -value 50
    #varset nconstr     -value 0
    
    varset calculation     -value 'scf'
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
    #tableset constraints_table 1 1 -value 1

    # unused variables
    groupwidget unused_1 disable
}
