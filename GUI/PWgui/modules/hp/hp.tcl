source commands.tcl

module HP\#auto -title "PWSCF GUI: module HP.x" -script {
    
    readfilter  ::pwscf::hpReadFilter

    #
    # Namelist: INPUTHP
    #
    
    namelist inputhp -name "INPUTHP" {

        separator -label "--- Control specs ---"
        
        var prefix {
            -label    "Prefix of outdir files saved by program PW.X (prefix):" 
            -fmt      %S -validate string
        }

        var outdir {
            -label    "Outdir directory where PW.X files resides (outdir):"
            -widget   entrydirselectquote
            -fmt      %S -validate string
        }

        var iverbosity {
            -label "Verbosity (iverbosity):"
            -textvalue { "1 = low" "2 = medium" "3 = high" "4 = highest" }
            -value { 1 2 3 4 }
            -validate int
            -widget optionmenu
        }

        var max_seconds {
            -label    "Maximum CPU time \[in seconds\] (max_seconds):"
            -validate posint
        }

        separator -label "--- Q-points ---"

        group q_points {
            packwidgets left
            var nq1 -label "nq1:"  -widget spinint  -validate posint
            var nq2 -label "nq2:"  -widget spinint  -validate posint
            var nq3 -label "nq3:"  -widget spinint  -validate posint
        }    

        var skip_equivalence_q {
            -label "Skip the equivalence analysis of q points (skip_equivalence_q):"
            -widget    radiobox
            -textvalue { Yes No }	      
            -value     { .true. .false. }
        }


        group qp -decor normal {
            auxilvar only_q {
                -label "Computes only the q points from start_q to last_q:"
                -widget    radiobox
                -textvalue { Yes No }	      
                -value     { .true. .false. }
            }
            group start_last_q {
                packwidgets left
                var start_q {
                    -label "First q-point (start_q):"
                    -validate int
                }            
                var last_q {
                    -label "Last q-point (last_q):"
                    -validate int
                }
            }
            var sum_pertq {
                -label "Collect pieces of the response occupation matrices (sum_pertq):"
                -widget    radiobox
                -textvalue { Yes No }	      
                -value     { .true. .false. }
            }        
        }

        separator -label "--- Perturbations ---"

        var determine_num_pert_only {
            -label "Determine the number of perturbations (determine_num_pert_only):"
            -widget    radiobox
            -textvalue { Yes No }
            -value     { .true. .false. }
        }
        
        var find_atpert {
            -label "Method for searching which atoms must be perturbed (find_atpert):"
            -textvalue {
                "1 = find atoms by analyzing unperturbed occupations"
                "2 = find atoms from different Hubbard atomic types"
                "3 = find atoms by symmetry"
            }
            -value { 1 2 3 }                    
            -widget optionmenu
        }

        var docc_thr {
            -label "Threshold for a comparison of unperturbed occupations (docc_thr):"
            -validate fortranreal
        }

        separator -label "--- Atom types specs ---"

        group ntyp_vars {
            auxilvar ntyp {
                -label "Number of type of atoms:"
                -validate posint
                -widget spinint
                -default 1
            }
            
            dimension skip_type {
                -label "Skip i-th atom-type in linear-response calculation (skip_type):"
                -widget    radiobox
                -textvalue { Yes No }	      
                -value     { .true. .false. }
                -start 1 -end 1
            }
            
            dimension equiv_type {
                -label "Make atom-type i equivalent to atom-type j (equiv_type):"
                -validate posint
                -start 1 -end 1
            }
        
            dimension perturb_only_atom {
                -label "Perturb only i-th atom types (perturb_only_atom):"
                -widget    radiobox
                -textvalue { .true. .false. }
                -value     { .true. .false. }
                -start 1 -end 1
            }
        }

        separator -label "--- Miscellaneous ---"
        
        var compute_hp {
            -label "Collect pieces of the chi0 and chi matrices (compute_hp):"
            -widget    radiobox
            -textvalue { Yes No }	      
            -value     { .true. .false. }
        }
        
        var conv_thr_chi {
            -label "Convergence threshold for the response function chi (conv_thr_chi):"
            -validate fortranreal
        }

        var thresh_init {
            -label "Threshold for the 1st iteration of the linear system (thresh_init):"
            -validate fortranreal
        }

        var ethr_nscf {
            -label "Threshold for the convergence of eigenvalues in NSCF (ethr_nscf):"
            -validate fortranreal
        }

        var niter_max {
            -label "Maximum number of iterations for linear-response calc. (niter_max):"
            -validate int
        }

        var alpha_mix {
            -variable alpha_mix(1) 
            -label "Mixing parameter (alpha_mix(1)):"
            -validate fortranreal
        }
        
        var nmix {
            -label "Number of iterations to use in Broyden mixing (nmix):"
            -validate int
        }

        var num_neigh {
            -label "Number of neighbors of every Hubbard atom to consider for V (num_neigh):"
            -validate int
        }

        var lmin {
            -label "Minimum orbital quantum number of the Hubbard atoms (lmin):"
            -validate int
        }
        
        var rmax {
            -label "Maximum neighbor distance (in Bohr) between two atoms (rmax):"
            -validate fortranreal
        }
    }

    # ----------------------------------------------------------------------
    # take care of specialties
    # ----------------------------------------------------------------------
    source hp-event.tcl
    
    # ------------------------------------------------------------------------
    # source the HELP file
    # ------------------------------------------------------------------------
    source hp-help.tcl
}
