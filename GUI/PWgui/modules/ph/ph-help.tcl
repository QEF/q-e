help title_ph -vartype character -helpfmt txt2html -helptext {
            One line of comment describing the job.
}

help amass -vartype real -helpfmt txt2html -helptext {
    Atomic mass [amu] of each atomic type.
}

help outdir -vartype character -helpfmt txt2html -helptext {     
            scratch directory where data files written by pw.x reside
<p> ( default = './' )
}

help recover -vartype logical -helpfmt -helptext {   
            if .true. restart from an interrupted run 
<p> ( default = .false. )
}

help prefix -vartype character -helpfmt txt2html -helptext {
            prepended to input/output filenames
            must be the same used in the calculation
            of unperturbed system with pw.x
<p> ( default = 'pwscf' )
}

help niter_ph -vartype integer -helpfmt txt2html -helptext {
            maximum number of iterations in a scf step             
<p> ( default = 50 )
}

help tr2_ph -vartype real -helpfmt txt2html -helptext {
            threshold for selfconsistency                     
<p> ( default = 1e-10 )
}

help alpha_mix1 -vartype real -helpfmt txt2html -helptext {
            mixing factor (for each iteration) for   
            updating the scf potential:
            vnew(in) = alpha_mix*vold(out) + (1-alpha_mix)*vold(in)
<p> ( default = 0.7 )
}

help nmix_ph -vartype integer -helpfmt txt2html -helptext {
            number of iterations used in potential mixing          
<p> ( default = 4 )
}

help iverbosity -vartype integer -helpfmt txt2html -helptext {
              iverbosity=0 short output, iverbosity=1 verbose output 0
}

help max_seconds  -vartype integer -helpfmt txt2html -helptext {
              maximum allowed run time before the job stops smoothly
}

help fildyn   -vartype character -helpfmt txt2html -helptext {
              file where the dynamical matrix is written           
<p> ( default = 'matdyn' )
}

help fildrho  -vartype character -helpfmt txt2html -helptext {
              file where the charge density variations is written     
<p> ( default = ' ' )
}

help filelph  -vartype character -helpfmt txt2html -helptext {
              file where electron-phonon matrix elements are written  
<p> ( default = ' ' )
}


help fildvscf -vartype character -helpfmt txt2html -helptext {
              file where the the potential variation is written       
              (for later use in electron-phonon calculation)
<p> ( default = ' ' )
}

help maxirr   -vartype integer -helpfmt txt2html -helptext {
            maximum number of irreducible representation to be      
            computed in a single run
            maxirr=0 --> compute all representations
}

help epsil   -vartype logical -helpfmt txt2html -helptext {
            if .true. in a q=0 calculation for a non metal the     
            macroscopic dielectric constant of the system is 
            computed. 
<p> ( default =  .false. )
}

help trans   -vartype logical -helpfmt txt2html -helptext {
            if .true. the phonons are computed                      
            If trans.and.epsil effective charges are calculated            
<p> ( default =  .true. )
}

help elph    -vartype logical -helpfmt txt2html -helptext {
            if .true. electron-phonon lambda coeffs are computed    
            For metals only, requires gaussian smearing.
            If elph.and.trans, the lambdas are calculated in the same
            run, using the same k-point grid for phonons and lambdas
            If elph.and..not.trans, the lambdas are calculated using
            previously saved DeltaVscf in fildvscf, previously saved
            dynamical matrix, and the present punch file. This allows
            the use of a different (larger) k-point grid.
<p> ( default =  .false. )
}

help zue     -vartype logical -helpfmt txt2html -helptext {
            if .true. in a q=0 calculation for a non metal the 
            effective charges are computed from the phonon
            density responses. Note that if trans.and.epsil
            effective charges are calculated using a different
            algorithm. The results should be the same within
            numerical noise.
}

set _xq {
            the phonon wavevector; must be equal to the one used
            in the non-selfconsistent calculation.
}
foreach var {xq1 xq2 xq3} {
    help $var -vartype real -helpfmt txt2html -helptext $_xq
}

help lraman -vartype logical -helpfmt txt2html -helptext {
            if .true. calculate nonresonant Raman coefficients      
            using second-order response as in:
            M. Lazzeri and F. Mauri, Phys. Rev. Lett. 90, 036401 (2003)
<p> ( default = .false. )

            Optional variables for Raman:

    eth_rps threshold for calculation of  Pc R |psi>                1.0d-9
    eth_ns  threshold for non-scf wavefunction calculation          1.0e-12
    dek     delta_xk used for wavefunction derivation wrt k         1.0e-3
}

help eth_rps -vartype real -helpfmt txt2html -helptext {
    threshold for calculation of  Pc R |psi>, used in Raman calculation
    Default value: 1.0d-9
}
help eth_ns -vartype real -helpfmt txt2html -helptext {
    threshold for non-selfconsistent calculation of wavefunctions
    at k+dk, used in Raman calculation
    Default value: 1.0d-12
}
help dek -vartype real -helpfmt txt2html -helptext {
    dk used for finite-difference derivation of wavefunctions: dpsi_k/dk, 
    for Raman calculation
    Default value: 1.0e-3
}


help fpol -vartype real -helpfmt txt2html -helptext {
    if .true. calculate dynamic polarizabilities            
    ( experimantal stage, see example33 for calculation
      of methane )
<p> ( default = .false. )
}



help nrapp -vartype integer -helpfmt txt2html -helptext {                             
    The representations to do.
    (nrapp and nat_todo are not compatible)
<p> ( default = 0 )
}

help nat_todo    -vartype integer -helpfmt txt2html -helptext {  
    Number of atom to be displaced.
    (nrapp and nat_todo are not compatible)
<p> ( default = 0 )
}

help modenum  -vartype integer -helpfmt txt2html -helptext { 
    Used for single mode calculation.
    If not set here, will be read from file .save           
    (if modenum.ne.0 => nat_todo = 0; nrapp = 1; list(1) = 
     modenum)
<p> ( default = -1 )
}


