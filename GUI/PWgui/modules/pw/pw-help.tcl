# ------------------------------------------------------------------------
#
# This File is the HELP file for the PWSCF.
#
# ------------------------------------------------------------------------

::tclu::DEBUG help: reading help


# =============================================================================
# NAMELIST &CONTROL
# =============================================================================

help calculation -vartype character -helpfmt html -helptext {
    A string describing the task to be performed.
    <p>
    <b>Allowed values:</b> 'scf', 'nscf', 'phonon', 'relax', 
                           'md', 'vc-relax', 'vc-md', 'neb', 
                           'smd' (vc = variable-cell). 
    <p>
    <b>Default:</b> 'scf'
}

help title -vartype character -helpfmt txt2html -helptext {
             reprinted on output. Default: ' '
}

help verbosity    -vartype character -helpfmt txt2html -helptext { 
             Allowed values: 'high' | 'default' | 'low' | 'minimal'
}

help restart_mode -vartype character -helpfmt txt2html -helptext { 
             'from_scratch'  : from scratch ( default )
	                       NEB and SMD only: the starting path is obtained
                               with a linear interpolation between the images
                               specified in the ATOMIC_POSITIONS card.
			       Note that in the linear interpolation
                               periodic boundary conditions ARE NON USED.
             'restart'       : from previous interrupted run
}

help nstep        -vartype integer -helpfmt txt2html -helptext {                          
             number of ionic+electronic steps
             default:  1 if calculation = 'scf', 'nscf'
                       0 if calculation = 'neb', 'smd'
                      50 for the other cases
}

help iprint       -vartype integer -helpfmt txt2html -helptext { 
             band energies are written every iprint iterations
             default: write only at convergence
}

help tstress      -vartype logical -helpfmt txt2html -helptext { 
             calculate stress. Set to .TRUE. if calculation = 'vc-md'
}

help tprnfor      -vartype logical -helpfmt txt2html -helptext { 
             print forces. Set to .TRUE. if calculation = 'relax', 'md', 'vc-md'
}
            
help dt           -vartype real -helpfmt txt2html -helptext { 
             time step for molecular dynamics, in Rydberg atomic units
             (1 a.u.=4.8378 * 10^-17 s : beware, CP and FPMD codes use
              Hartree atomic units, half that much!!!)
}

help outdir       -vartype character -helpfmt txt2html -helptext { 
             input, temporary, output files are found in this directory
             Default: current directory ('./')
}
            
help prefix       -vartype character -helpfmt txt2html -helptext { 
             prepended to input/output filenames: 
             prefix.wfc, prefix.rho, etc. Default: 'pwscf'
}

help max_seconds  -vartype integer -helpfmt txt2html -helptext { 
             jobs stops after max_seconds CPU time
}

help etot_conv_thr     -vartype real -helpfmt txt2html -helptext { 
             convergence threshold on total energy (a.u)
             for ionic minimization. Default: 1.D-4
             See also forc_conv_thr - both criteria must be satisfied
}

help forc_conv_thr     -vartype real -helpfmt txt2html -helptext { 
             convergence threshold on forces (a.u)
             for ionic minimization. Default: 1.D-3
             Ssee also etot_conv_thr - both criteria must be satisfied
}

help disk_io      -vartype character -helpfmt txt2html -helptext { 
             'high', 'default', 'low', 'minimal'
}

help pseudo_dir   -vartype character -helpfmt txt2html -helptext {
             directory containing pseudopotential files
             Default: '$HOME/pw/pseudo/'
}
	
help tefield      -vartype logical -helpfmt txt2html -helptext {  
             If .TRUE. a sawlike potential is added to the 
               bare ionic potential.  Default: .FALSE.
}	



# =============================================================================
# NAMELIST &SYSTEM
# =============================================================================

set ibravCelldm {
The specification of the <b>celldm</b> dimension depends on the value
of <b>ibrav</b> variable.

<pre>
     ibrav        structure                       celldm(2)-celldm(6)

       0          "free", see cards 2bis above      not used
       1          cubic P (sc)                      not used
       2          cubic F (fcc)                     not used   
       3          cubic I (bcc)                     not used

       4          Hexagonal and Trigonal P        celldm(3)=c/a
       5          Trigonal R                      celldm(4)=cos(aalpha)
       6          Tetragonal P (st)               celldm(3)=c/a
       7          Tetragonal I (bct)              celldm(3)=c/a
       8          Orthorhombic P                  celldm(2)=b/a,
                                                  celldm(3)=c/a

       9          Orthorhombic base-centered(bco) celldm(2)=b/a,
	                                          celldm(3)=c/a

      10          Orthorhombic face-centered      celldm(2)=b/a,
	                                          celldm(3)=c/a

      11          Orthorhombic body-centered      celldm(2)=b/a,
                                                  celldm(3)=c/a

      12          Monoclinic P                    celldm(2)=b/a,
                                                  celldm(3)=c/a,
                                                  celldm(4)=cos(ab)

      13          Monoclinic base-centered        celldm(2)=b/a,
                                                  celldm(3)=c/a,
                                                  celldm(4)=cos(ab)

      14          Triclinic P                     celldm(2)= b/a,
                                                  celldm(3)= c/a,
                                                  celldm(4)= cos(bc),
                                                  celldm(5)= cos(ac),
                                                  celldm(6)= cos(ab)
</pre>

   The special axis is the z-axis, one basal-plane vector is along x, 
   and the other basal-plane vector is at angle beta for monoclinic 
   (beta is not actually used), at 120 degrees for trigonal and hexagonal(p)
   groups, and at 90 degrees for remaining groups, excepted fcc, bcc, 
   tetragonal(i), for which the crystallographic vectors are as follows:

<h2>FCC bravais lattice:</h2>

   a1=(a/2)(-1,0,1), a2=(a/2)(0,1,1), a3=(a/2)(-1,1,0).

<h2>BCC bravais lattice:</h2>

   a1=(a/2)(1,1,1), a2=(a/2)(-1,1,1), a3=(a/2)(-1,-1,1).

<h2>TETRAGONAL (I) bravais lattices:</h2>

   a1=(a/2,a/2,c/2), a2=(a/2,-a/2,c/2), a3=(-a/2,-a/2,c/2).

<h2>TRIGONAL(R) groups:</h2>

   For these groups, the z-axis is chosen as the 3-fold axis, but the
   crystallographic vectors form a three-fold star around the z-axis,
   and the primitive cell is a simple rhombohedron. if c is the cosine
   of the angle between any pair of crystallographic vectors, and if
   tx=sqrt((1-c)/2), ty=sqrt((1-c)/6), tz=sqrt((1+2c)/3), the crystal-
   lographic vectors are:

<p>
         a1=a(0,2ty,tz),  a2=a(tx,-ty,tz),  a3=a(-tx,-ty,tz).

<h2>BCO base centered orthorombic:</h2>

   a1=(a/2,b/2,0), a2=(-a/2,b/2,0), a3=(0,0,c)
}

help ibrav  -vartype integer -helpfmt html -helptext $ibravCelldm
help celldm -vartype real    -helpfmt html -helptext $ibravCelldm

foreach var {A B C cosAB cosAC cosBC} {
    help $var -vartype real -helpfmt txt2html -helptext { 
	     Traditional crystallographic constants (A,B,C in ANGSTROM)
             specify either these or celldm but not both.
    }
}
 
help nat          -vartype integer -helpfmt txt2html -helptext { 
             number of atoms in the unit cell - must be specified
}


help ntyp         -vartype integer -helpfmt txt2html -helptext { 
             number of types of atoms in the unit cell - must be specified
}

help nbnd -vartype integer -helpfmt html -helptext {
	Number of electronic states (bands) to be calculated.
	<p>
	<b>Default:</b>
        <ul>
        <li>for an insulator, nbnd = "number of valence bands"
                      (nbnd=nelec/2, see below for nelec)</li>
        <li>for a metal, 20% more (minimum 4 more) for a metal</li>
	</ul>
        Note that in spin-polarized calculations the number of
        k-point, not the number of bands per k-point, is doubled
}

help nelec        -vartype real -helpfmt txt2html -helptext { 
             number of electron in the unit cell
             (may be noninteger if you wish)
             Default: the same as ionic charge (neutral cell)
             A compensating jellium background is inserted
             to remove divergencies if the cell is not neutral
}

help ecutwfc      -vartype real -helpfmt txt2html -helptext { 
             kinetic energy cutoff (Ry) for wavefunctions
             (must be specified)
}

help ecutrho      -vartype real -helpfmt txt2html -helptext { 
             kinetic energy cutoff (Ry) for charge density and potential
             Default: 4*ecutwfc . May be larger (for ultrasoft PP)
             or somewhat smaller (but not much smaller) than this.
}

foreach var {nr1 nr2 nr3} {
    help $var  -vartype integer -helpfmt txt2html -helptext { 
             three-dimensional FFT mesh (hard grid) for charge
             density (and scf potential). If not specified
             the grid is calculated based on the cutoff for
             charge density (see also "ecutrho")
    }
}

foreach var {nr1s nr2s nr3s} {
    help $var  -vartype integer -helpfmt txt2html -helptext {
             three-dimensional mesh for wavefunction FFT
             and for the smooth part of charge density
             (smooth grid). Coincides with nr1,nr2,nr3 if
             ecutrho=4*ecutwfc (default)
    }
}

help nosym        -vartype logical -helpfmt txt2html -helptext { 
             if (.TRUE.) symmetry is not used. Note that a k-point grid
             provided in input is used "as is"; an automatically generated
             k-point grid will contain only points in the irreducible BZ 
             of the lattice.  Use with care in low-symmetry large cells 
             if you cannot afford a k-point grid with the correct symmetry.
             Default: .FALSE.
}

help starting_magnetization -vartype real -helpfmt txt2html -helptext { 
             starting spin polarization (values between -1 and 1)
             on atomic type 'i' in a lsda calculation. Breaks the
             symmetry and provides a starting point for self-consistency.
             The default value is zero, BUT a value MUST be specified for
             AT LEAST one atomic type in spin polarized calculations.
             If zero starting magnetization is specified, zero final
             magnetization will be obtained.
}

help occupations  -vartype character -helpfmt txt2html -helptext { 
             'smearing':    gaussian smearing for metals
                            requires a value for degauss
             'tetrahedra' : for metals and DOS calculation
                            (see PRB49, 16223 (1994))
                            Requires uniform grid of k-points,
                            automatically generated (see below)
             'fixed' :      for insulators with a gap
             'from_input':  The occupation are read from input file.
                            Presently works only with one k-point 
                            (LSDA allowed). 
}

help degauss      -vartype real -helpfmt txt2html -helptext { 
             value of the gaussian spreading (Ry) for brillouin-zone
             integration in metals. Default: 0.D0
}


help smearing     -vartype character -helpfmt txt2html -helptext { 
             'gaussian', 'gauss':  
                  ordinary Gaussian spreading (Default)
             'methfessel-paxton', 'm-p', 'mp':
                  Methfessel-Paxton first-order spreading
                  (see PRB 40, 3616 (1989)).
             'marzari-vanderbilt', 'cold', 'm-v', 'mv':
                  Marzari-Vanderbilt cold smearing
                  (see PRL 82, 3296 (1999))
             'fermi-dirac', 'f-d', 'fd':
                  smearing with Fermi-Dirac function
}

foreach var {nelup neldw} {
    help $var      -vartype real -helpfmt txt2html -helptext { 
             Number of spin-up and spin-down electrons, respectively.
             The sum must yield nelec !!! NOT YET USED !!!
    }
}

help nspin        -vartype integer -helpfmt txt2html -helptext { 
             nspin=1 :  non-polarized calculation (default)
             nspin=2 : spin-polarized calculation
}

help ecfixed      -vartype real -helpfmt txt2html -helptext { 
             parameters for modified functional to be used in
             variable-cell molecular dynamics (or in stress calculation)
 	     Deafult:      40.0
}
help qcutz        -vartype real -helpfmt txt2html -helptext {
             parameters for modified functional to be used in
             variable-cell molecular dynamics (or in stress calculation)
             Deafult:       0.0
}
help q2sigma      -vartype real -helpfmt txt2html -helptext {
             parameters for modified functional to be used in
             variable-cell molecular dynamics (or in stress calculation)
             Deafult:       0.1
}

help xc_type      -vartype character -helpfmt txt2html -helptext { 
             Exchange-correlation functional.
             Presently unused: XC functional is read from PP files
}

set _hubbard {
             parameters for LDA+U calculations
             If lda_plus_u = .TRUE. you must specify, for species I,
             the parameters U and (optionally) alpha of the Hubbard 
             model (both in eV). See:
             Anisimov, Zaanen, and Andersen, PRB 44, 943 (1991);
             Anisimov et al., PRB 48, 16929 (1993);
             Liechtenstein, Anisimov, and Zaanen, PRB 52, R5467 (1994); 
             Cococcioni and de Gironcoli (to be published)
}
help lda_plus_u     -vartype logical -helpfmt txt2html -helptext "$_hubbard\n          (default: .false.)"
help Hubbard_U      -vartype real -helpfmt txt2html -helptext "$_hubbard\n             (default: 0.0 for all species)"
help Hubbard_alpha  -vartype real -helpfmt txt2html -helptext "$_hubbard\n             (default: 0.0 for all species)"
help U_projector_type -vartype character -helpfmt txt2html -helptext {
    Only active when lda_plus_U is .true., specifies the type
    of projector on localized orbital to be used in the LDA+U
    scheme.

    Currently available choices:
    'atomic': use atomic wfc's (as they are) to build the projector 
    'ortho-atomic': use Lowdin orthogonalized atomic wfc's
    'file': use the information from file "prefix".atwfc that must
        have been generated previously, for instance by pmw.x
        (see PP/poormanwannier.f90 for details)

    NB: forces and stress currently implemented only for the
    'atomic' choice.
}

help edir         -vartype integer -helpfmt txt2html -helptext { 
             1, 2 or 3. Used only if tefield is .TRUE. The direction of the
             electric field is parallel to the bg(.,edir) reciprocal 
             lattice vector (So the potential is constant in planes 
             defined by the mesh points).
}

help emaxpos      -vartype real -helpfmt txt2html -helptext { 
             Position of the maximum of the sawlike potential within the 
             unit cell (0<emaxpos<1). Default: 0.5D0
}

help eopreg       -vartype real -helpfmt txt2html -helptext { 
             Part of the unit cell where the sawlike potential decreases.
             ( 0 < eopreg < 1 ). Default: 0.1D0
}

help eamp         -vartype real -helpfmt txt2html -helptext { 
             Amplitude of the electric field (in a.u.) ( 1 a.u. = 51.44 10^11 V/m )
             Default: 0.001 a.u.
}



# =============================================================================
# NAMELIST &ELECTRONS
# =============================================================================

help electron_maxstep  -vartype integer -helpfmt txt2html -helptext { 
             maximum number of iterations in a scf step
             (default: 50)
}

help conv_thr     -vartype real -helpfmt txt2html -helptext { 
             Convergence threshold for selfconsistency: 
             estimated energy error < conv_thr ( Default 1.D-6 )
}

help mixing_mode  -vartype character -helpfmt txt2html -helptext { 
             'plain'    : charge density Broyden mixing (default)
             'TF'       : as above, with simple Thomas-Fermi screening
                          (for highly homogeneous systems)
             'local-TF' : as above, with local-density-dependent TF screening
                          (for highly inhomogeneous systems)
             'potential': (obsolete) potential mixing
}

help mixing_beta  -vartype real -helpfmt txt2html -helptext { 
             mixing factor for self-consistency (default: 0.7D0)
}

help mixing_ndim  -vartype integer -helpfmt txt2html -helptext { 
             number of iterations used in mixing scheme (default: 8)
}

help mixing_fixed_ns   -vartype integer -helpfmt txt2html -helptext { 
             For LDA+U : number of iterations with fixed ns (ns is the
             atomic density appearing in the Hubbard term) (default: 0)
}

help diagonalization   -vartype character -helpfmt txt2html -helptext { 
            'david': Davidson iterative diagonalization with overlap matrix
                     (default)
            'diis' : DIIS-like diagonalization
            'cg'   : conjugate-gradient-like band-by-band diagonalization
}

help diago_thr_init -vartype real -helpfmt txt2html -helptext {
               Convergence threshold for the firts iterative diagonalization.
               The threshold (ethr) is automatically updated along the
               self consistency loop.
               <p>Default: 1.D-2
}

help diago_cg_maxiter  -vartype integer -helpfmt txt2html -helptext { 
             For conjugate gradient diagonalization:
             max number of iterations
}

help diago_david_ndim  -vartype integer -helpfmt txt2html -helptext { 
             For Davidson diagonalization: dimension of workspace 
             (number of wavefunction packets, at least 2 needed). 
             A larger value may yield a faster algorithm but uses 
             more memory.<p>Default: 4
}

help diago_diis_ndim   -vartype integer -helpfmt txt2html -helptext { 
             For DIIS: dimension of the reduced space.<p>Default: 3
}


help startingpot  -vartype character -helpfmt txt2html -helptext { 
             'atomic': starting potential from atomic charge superposition
                       (default for scf,*relax,*md)
             'file'  : start from existing "prefix".pot file
                       (default and only possibility for nscf and phonon)
}

help startingwfc  -vartype character -helpfmt txt2html -helptext { 
             'atomic': start from superposition of atomic orbitals (default)
                       If not enough atomic orbitals are available,
                       fill with random numbers the remaining wfcs
             'random': start from random wfcs
             'file'  : start from a wavefunction file
}



# =============================================================================
# NAMELIST &IONS
# =============================================================================

help ion_dynamics      -vartype character -helpfmt txt2html -helptext { 
               specify the type of ionic dynamics. 
               For different type of calculation different possibilities are 
               allowed and different default values apply:
               
	       CASE ( calculation = 'relax' )
                 'bfgs' :   (default)   a new BFGS quasi-newton algorithm, based
                                        on the trust radius procedure, is used 
                                        for structural relaxation (experimental)
                 'old-bfgs' :           use the old BFGS quasi-newton method for
                                        structural relaxation
		 'damp' :               use damped (quick-min velocity Verlet) 
                                        dynamics for structural relaxation
                 'constrained-damp' :   use damped (quick-min velocity Verlet) 
                                        dynamics for structural relaxation with 
                                        the constraint specified in the 
                                        CONSTRAINTS CARD
               CASE ( calculation = 'md' )
                 'verlet' : (default)   use velocity Verlet algorithm to 
                                        integrate Newton's equation
                 'constrained-verlet' : use velocity Verlet algorithm to do 
                                        molecular dynamics with the constraint
                                        specified in the CONSTRAINTS CARD 
               CASE ( calculation = 'vc-relax' )
                 'damp' :   (default)   use damped (Beeman) dynamics for 
		                        structural relaxation
               CASE ( calculation = 'vc-md' )
                 'beeman' : (default)   use Beeman algorithm to integrate 
                                        Newton's equation
}
             
help ion_temperature   -vartype character -helpfmt txt2html -helptext { 
             'nose'          : Nose' thermostat, not implemented
             'rescaling'     : velocity rescaling (sort of implemented)
             'not_controlled': default
}

help tempw        -vartype real -helpfmt txt2html -helptext { 
             starting temperature in MD runs
}

help tolp         -vartype real -helpfmt txt2html -helptext { 
             tolerance for velocity rescaling. Velocities are
             not rescaled if the ratio of the run-averaged and 
             target temperature differs from unit less than tolp
             (default: 1.D-3)
}

help upscale      -vartype real -helpfmt txt2html -helptext { 
               max reduction factor for conv_thr during structural optimization
               conv_thr is automatically reduced when the relaxation 
               approaches convergence so that forces are still accurate,
               but conv_thr will not be reduced to less that 
               conv_thr / upscale. (default: 10.D0)
}

help potential_extrapolation     -vartype character -helpfmt txt2html -helptext { 
               used to extrapolate the potential and the wave-functions 
               from preceding ionic step(s)
               
               'none':   no extrapolation
               'atomic': extrapolate the potential as if it was a sum of
                         atomic-like orbitals (default for calculation='relax')
               'wfc':    extrapolate the potential as above
                         extrapolate wave-functions with first-order formula
                         (default for calcualtion = 'md', 'neb', 'smd')
               'wfc2':   as above, with second order formula
}


help lbfgs_ndim    -vartype integer -helpfmt txt2html -helptext {
               Number of old forces and displacements vectors used in the
	       linear scaling BFGS algorithm. When lbfgs_ndim = 1 the complete
	       inverse Hessian is stored (suggested for small/medium-size 
	       systems).<p>
	       On large systems (some hundreds of atoms) a good performance can 
	       be achieved with only 4 or 6 old vectors 
               (bfgs only)
<p> ( default = 1 )
}

help trust_radius_max -vartype  real -helpfmt txt2html -helptext {
               maximum ionic displacement in the structural relaxation 
	       (bfgs only)
<p> ( default = 0.5D0 BOHR ) 
}

help trust_radius_min -vartype  real -helpfmt txt2html -helptext {
	       minimum ionic displacement in the structural relaxation
               BFGS is reset when trust_radius < trust_radius_min
	       (bfgs only)
<p> ( default = 1.D-5 BOHR )
}

help trust_radius_ini -vartype real -helpfmt txt2html -helptext {
               initial ionic displacement in the structural relaxation
               (bfgs only)
<p> ( default = 0.5D0 BOHR )
}

help trust_radius_end -vartype real -helpfmt txt2html -helptext {
               BFGS is stopped when trust_radius < trust_radius_end
	       trust_radius_end is not intended to be used as a criterium
	       for convergence (bfgs only)
<p> ( default = 1.D-7 BOHR )
}


set _w12 {
	       parameters used in line search based on the Wolfe conditions
	       (bfgs only)
<p> ( default : w_1 = 1.D-5, w_2 = 0.2D0 )

}     
foreach var {w_1 w_2} {
    help $var -vartype real -helpfmt txt2html -helptext $_w12
}

help num_of_images -vartype integer -helpfmt txt2html -helptext {
               number of points used to discrtize the path
<p> ( default = 0 )
}

help CI_scheme    -vartype   character -helpfmt txt2html -helptext {
               specify the type of Climbing Image scheme
               "no-CI"      : climbing image is not used
               "highest-TS" : original CI scheme. The image highest in energy 
	                      does not feel the effect of springs and is 
			      allowed to climb along the path
               "manual"     : images that have to climb are manually selected. 
	                      See also CLIMBING_IMAGES card 
<p> ( default = "no-CI" )
}

help first_last_opt  -vartype  logical -helpfmt txt2html -helptext {
               also the first and the last configurations are optimized
               "on the fly" 
	       (these images do not feel the effect of the springs)
<p> ( default = .FALSE. )
}

help minimization_scheme  -vartype    character -helpfmt txt2html -helptext {
               specify the type of optimization scheme      
               "sd"         : steepest descent
	       "quick-min"  : a minimization algorithm based on
	                      molecular dynamics (suggested)
               "damped-dyn" : damped molecular dynamics. See also the 
	                      keyword damp
               "mol-dyn"    : constant temperature molecular dynamics. See 
	                      also the keyword temp_req.
	                      Note that, in order to perform such molecular 
			      dynamics, spring forces are NOT projected 
			      along the path. 
<p> ( default = "quick-min" )
}

help damp        -vartype    real -helpfmt txt2html -helptext {
               Damping coefficent. Ignored when "minimization_scheme"
               is different from "damped-verlet"
<p> ( default = 1.D0 )
}

help temp_req      -vartype  real -helpfmt txt2html -helptext {
               temperature associated to the elastic band. Each image has its 
	       own thermostat. The temperature in the output is the average 
	       temperature of the elastic band computed before the 
	       thermalization
               ignored when "minimization_scheme" is different from
               "mol-dyn"
<p> ( default = 0.D0 Kelvin )
}

help ds          -vartype    real -helpfmt txt2html -helptext {
               optimization step length ( Hartree atomic units )
<p> ( default = 1.5D0 )
}
	
set _k {              
               set them to use a Variable Elastic Constants scheme 
	       elastic constants are in the range [ k_min, k_max ] 
	       this is useful to rise the resolution around the saddle point
<p> ( default = 0.1D0 Hartree atomic units )
}
foreach var {k_max k_min} {
    help $var -vartype real -helpfmt txt2html -helptext $_k
}

help path_thr      -vartype   real -helpfmt txt2html -helptext {
               the simulation stops when the error ( the norm of the force 
	       orthogonal to the path in eV/A ) is less than path_thr.
<p> ( default = 0.05D0 eV / Angstrom )
}

help reset_vel     -vartype logical -helpfmt txt2html -helptext {
                used to reset quick-min velocities at restart time
               (sort of clean-up of the history)
<p> ( default = .FALSE. )
}

# =============================================================================
# NAMELIST &CELL
# =============================================================================

help cell_dynamics     -vartype character -helpfmt txt2html -helptext { 
               specify the type of dynamics for the cell. 
               For different type of calculation different possibilities 
               are allowed and different default values apply:

               CASE ( calculation = 'vc-relax' )
                 'none':    default 
                 'sd':      steepest descent ( not implemented )
                 'damp-pr': damped (Beeman) dynamics of the Parrinello-Raman 
                            extended lagrangian
                 'damp-w':  damped (Beeman) dynamics of the new Wentzcovitch
                            extended lagrangian
               CASE ( calculation = 'vc-md' )
                 'none': default 
                 'pr':      (Beeman) molecular dynamics of the Parrinello-Raman 
                            extended lagrangian
                 'w':       (Beeman) molecular dynamics of the new Wentzcovitch
                            extended lagrangian
}

help press        -vartype real -helpfmt txt2html -helptext { 
             target pressure [KBar] in a variable-cell md simulation
             (default: 0.D0)
}

help wmass        -vartype real -helpfmt txt2html -helptext { 
             ficticious cell mass for variable-cell md simulations
}

help cell_factor  -vartype real -helpfmt txt2html -helptext { 
             used in the construction of the pseudopotential tables. 
             It should exceed the maximum linear contraction of the
             cell during a simulation. Default: 1.2D0
}


# =============================================================================
# NAMELIST &PHONON
# =============================================================================

help modenum      -vartype integer -helpfmt txt2html -helptext { 
             for single-mode phonon calculation (default:0)
}

help xqq       -vartype real -helpfmt txt2html -helptext { 
             q-point (units 2pi/a) for phonon calculation
}


# ==============================================================================
# CARD: CELL_PARAMETERS
# ==============================================================================

set _lattice {
CELL_PARAMETERS { cubic | hexagonal }
  optional card, needed only if ibrav = 0 is specified
  cubic    : assume cubic symmetry or a subset (default)  
  hexagonal: assume hexagonal symmetry or a subset
  Next cards:
    a(1,1) a(2,1) a(3,1)
    a(1,2) a(2,2) a(3,2)
    a(1,3) a(2,3) a(3,3)

  a(:,1) = crystal axis 1    alat units if celldm(1) was specified
      2                 2    a.u. if celldm(1)=0
      3                 3
}
help lattice -vartype {real, real, real} -helpfmt txt2html -helptext $_lattice
help lattice_type -vartype character -helpfmt txt2html -helptext $_lattice


# ==============================================================================
# CARD: ATOMIC_SPECIES
# ==============================================================================

help atomic_species -vartype {character, real, character} -helpfmt txt2html -helptext {
Syntax:

ATOMIC_SPECIES
 X(1)     Mass_X(1)     PseudoPot_X(ntyp)
 X(2)     Mass_X(2)     PseudoPot_X(ntyp)
 ...
 X(ntyp)  Mass_X(ntyp)  PseudoPot_X(ntyp)


Description:
 X           Character: label of the atom
 Mass_X      Real     : mass of the atomic species
                        not used if calculation='scf', 'nscf', 'phonon'
 PseudoPot_X Character: file containing PP for this species
}


# ==============================================================================
# CARD: ATOMIC_POSITIONS
# ==============================================================================

set _pos {
    Atomic positions can be loaded by pressing the 
    "Load atomic coordinates from file ..." button. The coordinates 
    should be written in PW.X syntax. This means that the "ATOMIC_POSITIONS"
    card will be parsed. Below is the description of its syntax.
    
    Syntax:
    
    ATOMIC_POSITIONS { alat | bohr | crystal | angstrom }
    if ( calculation = 'neb' .OR. 'smd' )
      first_image
      X 0.0  0.0  0.0  {if_pos(1) if_pos(2) if_pos(3)}
      Y 0.5  0.0  0.0
      Z O.0  0.2  0.2
      last_image
      X 0.0  0.0  0.0
      Y 0.7  0.0  0.0
      Z O.0  0.5  0.2 
    else
      X 0.0  0.0  0.0  {if_pos(1) if_pos(2) if_pos(3)}
      Y 0.5  0.0  0.0
      Z O.0  0.2  0.2
    
    Description:
    
       alat    : atomic positions are in units of alat (default)
       bohr    : atomic positions are in a.u.
       crystal : atomic positions are in crystal coordinates (see below)
       angstrom: atomic positions are in ANGSTROMS    
}
help atomic_coordinates \
    -vartype {character, real, real, real, \[ingeter, ingteger, integer\]} \
    -helpfmt txt2html -helptext $_pos
help atomic_coordinates_last \
    -vartype {character, real, real, real} \
    -helpfmt txt2html -helptext $_pos

# BEWARE: it is assumed that 50 intermediate images is the
# largest allowed number (this is dirty)
for {set i 1} {$i <= 50} {incr i} {
    help atomic_coordinates_${i}_inter \
	-vartype {character, real, real, real} \
	-helpfmt txt2html -helptext $_pos
}



# ==============================================================================
# CARD: K_POINTS 
# ==============================================================================

set _kpoints0 {
    K-points coordinates can be loaded by pressing the 
    "Load K-point coordinates from file ..." button. The coordinates 
    should be written in PW.X syntax. This means that the "K_POINTS"
    card will be parsed. Below is the description of its syntax.
}
set _kpoints {
Syntax:

K_POINTS { tpiba | automatic | crystal | gamma }
if (gamma)
   nothing to read
elseif (automatic) 
   nk1  nk2  nk3  sk1  sk2  sk3
else
   nks
   xk_x(1)    xk_y(1)    xk_z(1)     wk(1)
   xk_x(2)    xk_y(3)    xk_z(4)     wk(5)
   ...
   xk_x(nks)  xk_y(nks)  xk_z(nks)   wk(nks)
 

Description:

K_POINTS { tpiba | automatic | crystal | gamma }

   gamma    : use k = 0 ( do not read anything after this card )
              Note that a set of subroutines optimized for clculations at 
              the gamma point are used so that both memory and cpu requirements
              are reduced
   automatic: automatically generated uniform grid of k-points
              next card:
   nk1, nk2, nk3, sk1, sk2, sk3 :
              generates ( nk1, nk2, nk3 ) mesh with ( sk1, sk2, sk3 ) offset
              nk1, nk2, nk3 as in Monkhorst-Pack grids
              sk1, sk2, sk3 must be 0 ( no offset ) or 1 ( grid displaced 
              by half a grid step in the corresponding direction )
              The mesh with offset may not work with tetrahedra.
   crystal  : read k-points in crystal coordinates
   tpiba    : read k-points in 2pi/a units ( default )
              next card:
   nks
              number of supplied special points
   do i=1,nks
   ! xk_x, xk_y, xk_z, wk      : k-point coordinates and weight
   !           special points in the irreducible Brillouin Zone
   !           of the lattice (with all symmetries) and weights
   !           If the symmetry is lower than the full symmetry 
   !           of the lattice, additional points with appropriate
   !           weights are generated
   enddo 
}
help kpoint_type -vartype character -helpfmt txt2html -helptext $_kpoints
foreach var {nks nk1 nk2 nk3 sk1 sk2 sk3} {
    help $var -vartype integer -helpfmt txt2html -helptext $_kpoints
}
help kpoints -vartype {real, real, real, real} -helpfmt txt2html -helptext ${_kpoints0}\n$_kpoints


# ==============================================================================
# CARD: CLIMBING_IMAGES
# ==============================================================================

help climbing_images_var -vartype {list of comma-separated integers} -helpfmt txt2html -helptext {
Syntax:

CLIMBING_IMAGES
   index1, index2, ..., indexN

Description:
  Optional card, needed only if calculation = 'neb' and CI_scheme = 'manual',    
  where index1, index2, ..., indexN are the indices of the images to which 
  apply the Climbing Image procedure. If more than an image is specified they
  must be separated by a comma
}

# ==============================================================================
# CARD: CONSTRAINTS
# ==============================================================================

set constraints_help {

   Ionic Constraints

 Syntax:

    CONSTRAINTS
      nconstr   constr_tol
      constr_type(.)   constr(1,.)   constr(2,.)

 Where:

      nconstr (INTEGER)         INTEGER, number of constraints
      constr_tol                REAL,    tolerance for keeping the constraints 
                                         satisfied
      constr_type(.)            INTEGER, type of constrain
      constr(1,.) constr(2,.)   INTEGER, atoms indices object of the constraint.
                                        
                          I.E.: 1 ia1 ia2 "1" is the constrain type 
                                (fixed distance) "ia1 ia2" are the 
                                indices of the atoms (as they appear 
                                in the 'ATOMIC_POSITION' CARD) whose 
                                distance has to be kept constant

}

foreach var { nconstr constr_tol constraints_table } { 
help $var -helpfmt txt2html -helptext $constraints_help 
}
