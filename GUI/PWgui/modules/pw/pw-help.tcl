
#
# Help-file automatically created by helpdoc utility
#
#    !!! DO NOT EDIT: CHANGES WILL BE LOST !!!
#
	

# ------------------------------------------------------------------------
help calculation -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>calculation</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'scf'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
a string describing the task to be performed:
   'scf',
   'nscf',
   'bands',
   'relax',
   'md',
   'vc-relax',
   'vc-md'

   (vc = variable-cell).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help title -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>title</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
reprinted on output.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help verbosity -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>verbosity</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'low'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Currently two verbosity levels are implemented:
  'high' and 'low'. 'debug' and 'medium' have the same
  effect as 'high'; 'default' and 'minimal', as 'low'
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help restart_mode -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>restart_mode</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'from_scratch'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
'from_scratch'  : from scratch. This is the normal way
                  to perform a PWscf calculation
'restart'       : from previous interrupted run. Use this
                  switch only if you want to continue an
                  interrupted calculation, not to start a
                  new one, or to perform non-scf calculations.
                  Works only if the calculation was cleanly
                  stopped using variable "max_seconds", or
                  by user request with an "exit file" (i.e.:
                  create a file "prefix".EXIT, in directory
                  "outdir"; see variables "prefix", "outdir")
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help wf_collect -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>wf_collect</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
This flag controls the way wavefunctions are stored to disk :

.TRUE.  collect wavefunctions from all processors, store them
        into the output data directory "outdir"/"prefix".save,
        one wavefunction per k-point in subdirs K000001/,
        K000001/, etc.. Use this if you want wavefunctions
        to be readable on a different number of processors.

.FALSE. do not collect wavefunctions, leave them in temporary
        local files (one per processor). The resulting format
        will be readable only by jobs running on the same
        number of processors and pools. Requires less I/O
        than the previous case.

Note that this flag has no effect on reading, only on writing.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nstep -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nstep</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em>
1  if calculation = 'scf', 'nscf', 'bands';
50 for the other cases
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
number of ionic + electronic steps
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help iprint -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>iprint</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> write only at convergence
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
band energies are written every "iprint" iterations
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tstress -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tstress</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
calculate stress. It is set to .TRUE. automatically if
calculation='vc-md' or 'vc-relax'
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tprnfor -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tprnfor</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
calculate forces. It is set to .TRUE. automatically if
calculation='relax','md','vc-md'
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help dt -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>dt</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 20.D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
time step for molecular dynamics, in Rydberg atomic units
(1 a.u.=4.8378 * 10^-17 s : beware, the CP code uses
 Hartree atomic units, half that much!!!)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help outdir -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>outdir</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em>
value of the ESPRESSO_TMPDIR environment variable if set;
current directory ('./') otherwise
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
input, temporary, output files are found in this directory,
see also "wfcdir"
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help wfcdir -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>wfcdir</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> same as "outdir"
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
this directory specifies where to store files generated by
each processor (*.wfc{N}, *.igk{N}, etc.). Useful for
machines without a parallel file system: set "wfcdir" to
a local file system, while "outdir" should be a parallel
or networkfile system, visible to all processors. Beware:
in order to restart from interrupted runs, or to perform
further calculations using the produced data files, you
may need to copy files to "outdir". Works only for pw.x.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help prefix -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>prefix</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'pwscf'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
prepended to input/output filenames:
prefix.wfc, prefix.rho, etc.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lkpoint_dir -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lkpoint_dir</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .true.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .false. a subdirectory for each k_point is not opened
in the "prefix".save directory; Kohn-Sham eigenvalues are
stored instead in a single file for all k-points. Currently
doesn't work together with "wf_collect"
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help max_seconds -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>max_seconds</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.D+7, or 150 days, i.e. no time limit
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
jobs stops after "max_seconds" CPU time. Use this option
in conjunction with option "restart_mode" if you need to
split a job too long to complete into shorter jobs that
fit into your batch queues.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help etot_conv_thr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>etot_conv_thr</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.0D-4
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
convergence threshold on total energy (a.u) for ionic
minimization: the convergence criterion is satisfied
when the total energy changes less than "etot_conv_thr"
between two consecutive scf steps. Note that "etot_conv_thr"
is extensive, like the total energy.
See also "forc_conv_thr" - both criteria must be satisfied
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help forc_conv_thr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>forc_conv_thr</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.0D-3
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
convergence threshold on forces (a.u) for ionic minimization:
the convergence criterion is satisfied when all components of
all forces are smaller than "forc_conv_thr".
See also "etot_conv_thr" - both criteria must be satisfied
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help disk_io -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>disk_io</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> see below
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Specifies the amount of disk I/O activity
'high':   save all data to disk at each SCF step

'medium': save wavefunctions at each SCF step unless
          there is a single k-point per process (in which
          case the behavior is the same as 'low')

'low' :   store wfc in memory, save only at the end

'none':   do not save anything, not even at the end
          ('scf', 'nscf', 'bands' calculations; some data
           may be written anyway for other calculations)

Default is 'low' for the scf case, 'medium' otherwise.
Note that the needed RAM increases as disk I/O decreases!
It is no longer needed to specify 'high' in order to be able
to restart from an interrupted calculation (see "restart_mode")
but you cannot restart in disk_io='none'
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help pseudo_dir -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>pseudo_dir</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em>
value of the $ESPRESSO_PSEUDO environment variable if set;
'$HOME/espresso/pseudo/' otherwise
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
directory containing pseudopotential files
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tefield -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tefield</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE. a saw-like potential simulating an electric field
is added to the bare ionic potential. See variables "edir",
"eamp", "emaxpos", "eopreg" for the form and size of
the added potential.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help dipfield -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>dipfield</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE. and tefield=.TRUE. a dipole correction is also
added to the bare ionic potential - implements the recipe
of L. Bengtsson, PRB 59, 12301 (1999). See variables "edir",
"emaxpos", "eopreg" for the form of the correction. Must
be used ONLY in a slab geometry, for surface calculations,
with the discontinuity FALLING IN THE EMPTY SPACE.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lelfield -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lelfield</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE. a homogeneous finite electric field described
through the modern theory of the polarization is applied.
This is different from "tefield=.true." !
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nberrycyc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nberrycyc</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
In the case of a finite electric field  ( lelfield == .TRUE. )
it defines the number of iterations for converging the
wavefunctions in the electric field Hamiltonian, for each
external iteration on the charge density
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lorbm -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lorbm</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE. perform orbital magnetization calculation.
If finite electric field is applied (lelfield=.true.)
only Kubo terms are computed
[for details see New J. Phys. 12, 053032 (2010)].
The type of calculation is 'nscf' and should be performed
on an automatically generated uniform grid of k points.
Works ONLY with norm-conserving pseudopotentials.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lberry -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lberry</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE. perform a Berry phase calculation
See the header of PW/src/bp_c_phase.f90 for documentation
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help gdir -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>gdir</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
For Berry phase calculation: direction of the k-point
strings in reciprocal space. Allowed values: 1, 2, 3
1=first, 2=second, 3=third reciprocal lattice vector
For calculations with finite electric fields
(lelfield==.true.) "gdir" is the direction of the field
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nppstr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nppstr</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
For Berry phase calculation: number of k-points to be
calculated along each symmetry-reduced string
The same for calculation with finite electric fields
(lelfield=.true.)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ibrav -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ibrav</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Status: </em> REQUIRED
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
  Bravais-lattice index. If ibrav /= 0, specify EITHER
  [ celldm(1)-celldm(6) ] OR [ A,B,C,cosAB,cosAC,cosBC ]
  but NOT both. The lattice parameter "alat" is set to
  alat = celldm(1) (in a.u.) or alat = A (in Angstrom);
  see below for the other parameters.
  For ibrav=0 specify the lattice vectors in CELL_PARAMETER,
  optionally the lattice parameter alat = celldm(1) (in a.u.)
  or = A (in Angstrom), or else it is taken from CELL_PARAMETERS

ibrav      structure                   celldm(2)-celldm(6)
                                     or: b,c,cosab,cosac,cosbc
  0          free
      crystal axis provided in input: see card CELL_PARAMETERS

  1          cubic P (sc)
      v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,1)

  2          cubic F (fcc)
      v1 = (a/2)(-1,0,1),  v2 = (a/2)(0,1,1), v3 = (a/2)(-1,1,0)

  3          cubic I (bcc)
      v1 = (a/2)(1,1,1),  v2 = (a/2)(-1,1,1),  v3 = (a/2)(-1,-1,1)

  4          Hexagonal and Trigonal P        celldm(3)=c/a
      v1 = a(1,0,0),  v2 = a(-1/2,sqrt(3)/2,0),  v3 = a(0,0,c/a)

  5          Trigonal R, 3fold axis c        celldm(4)=cos(alpha)
      The crystallographic vectors form a three-fold star around
      the z-axis, the primitive cell is a simple rhombohedron:
      v1 = a(tx,-ty,tz),   v2 = a(0,2ty,tz),   v3 = a(-tx,-ty,tz)
      where c=cos(alpha) is the cosine of the angle alpha between
      any pair of crystallographic vectors, tx, ty, tz are:
        tx=sqrt((1-c)/2), ty=sqrt((1-c)/6), tz=sqrt((1+2c)/3)
 -5          Trigonal R, 3fold axis &lt;111&gt;    celldm(4)=cos(alpha)
      The crystallographic vectors form a three-fold star around
      &lt;111&gt;. Defining a' = a/sqrt(3) :
      v1 = a' (u,v,v),   v2 = a' (v,u,v),   v3 = a' (v,v,u)
      where u and v are defined as
        u = tz - 2*sqrt(2)*ty,  v = tz + sqrt(2)*ty
      and tx, ty, tz as for case ibrav=5
      Note: if you prefer x,y,z as axis in the cubic limit,
            set  u = tz + 2*sqrt(2)*ty,  v = tz - sqrt(2)*ty
            See also the note in flib/latgen.f90

  6          Tetragonal P (st)               celldm(3)=c/a
      v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,c/a)

  7          Tetragonal I (bct)              celldm(3)=c/a
      v1=(a/2)(1,-1,c/a),  v2=(a/2)(1,1,c/a),  v3=(a/2)(-1,-1,c/a)

  8          Orthorhombic P                  celldm(2)=b/a
                                             celldm(3)=c/a
      v1 = (a,0,0),  v2 = (0,b,0), v3 = (0,0,c)

  9          Orthorhombic base-centered(bco) celldm(2)=b/a
                                             celldm(3)=c/a
      v1 = (a/2, b/2,0),  v2 = (-a/2,b/2,0),  v3 = (0,0,c)
 -9          as 9, alternate description
      v1 = (a/2,-b/2,0),  v2 = (a/2,-b/2,0),  v3 = (0,0,c)

 10          Orthorhombic face-centered      celldm(2)=b/a
                                             celldm(3)=c/a
      v1 = (a/2,0,c/2),  v2 = (a/2,b/2,0),  v3 = (0,b/2,c/2)

 11          Orthorhombic body-centered      celldm(2)=b/a
                                             celldm(3)=c/a
      v1=(a/2,b/2,c/2),  v2=(-a/2,b/2,c/2),  v3=(-a/2,-b/2,c/2)

 12          Monoclinic P, unique axis c     celldm(2)=b/a
                                             celldm(3)=c/a,
                                             celldm(4)=cos(ab)
      v1=(a,0,0), v2=(b*cos(gamma),b*sin(gamma),0),  v3 = (0,0,c)
      where gamma is the angle between axis a and b.
-12          Monoclinic P, unique axis b     celldm(2)=b/a
                                             celldm(3)=c/a,
                                             celldm(5)=cos(ac)
      v1 = (a,0,0), v2 = (0,b,0), v3 = (c*cos(beta),0,c*sin(beta))
      where beta is the angle between axis a and c

 13          Monoclinic base-centered        celldm(2)=b/a
                                             celldm(3)=c/a,
                                             celldm(4)=cos(ab)
      v1 = (  a/2,         0,                -c/2),
      v2 = (b*cos(gamma), b*sin(gamma), 0),
      v3 = (  a/2,         0,                  c/2),
      where gamma is the angle between axis a and b

 14          Triclinic                       celldm(2)= b/a,
                                             celldm(3)= c/a,
                                             celldm(4)= cos(bc),
                                             celldm(5)= cos(ac),
                                             celldm(6)= cos(ab)
      v1 = (a, 0, 0),
      v2 = (b*cos(gamma), b*sin(gamma), 0)
      v3 = (c*cos(beta),  c*(cos(alpha)-cos(beta)cos(gamma))/sin(gamma),
           c*sqrt( 1 + 2*cos(alpha)cos(beta)cos(gamma)
                     - cos(alpha)^2-cos(beta)^2-cos(gamma)^2 )/sin(gamma) )
  where alpha is the angle between axis b and c
         beta is the angle between axis a and c
        gamma is the angle between axis a and b
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help celldm -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>celldm(i), i=1,6</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>See: </em> ibrav
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Crystallographic constants - see the "ibrav" variable.
Specify either these OR A,B,C,cosAB,cosBC,cosAC NOT both.
Only needed values (depending on "ibrav") must be specified
alat = celldm(1) is the lattice parameter "a" (in BOHR)
If ibrav=0, only celldm(1) is used if present;
cell vectors are read from card CELL_PARAMETERS
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {A B C cosAB cosAC cosBC} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>A, B, C, cosAB, cosAC, cosBC</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Traditional crystallographic constants: a,b,c in ANGSTROM
  cosAB = cosine of the angle between axis a and b (gamma)
  cosAC = cosine of the angle between axis a and c (beta)
  cosBC = cosine of the angle between axis b and c (alpha)
The axis are chosen according to the value of "ibrav".
Specify either these OR "celldm" but NOT both.
Only needed values (depending on "ibrav") must be specified
The lattice parameter alat = A (in ANGSTROM )
If ibrav = 0, only A is used if present;
cell vectors are read from card CELL_PARAMETERS
            </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
help nat -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nat</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Status: </em> REQUIRED
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
number of atoms in the unit cell
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ntyp -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ntyp</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Status: </em> REQUIRED
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
number of types of atoms in the unit cell
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nbnd -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nbnd</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em>
for an insulator, nbnd = number of valence bands
(nbnd = # of electrons /2);
for a metal, 20% more (minimum 4 more)
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
number of electronic states (bands) to be calculated.
Note that in spin-polarized calculations the number of
k-point, not the number of bands per k-point, is doubled
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tot_charge -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tot_charge</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
total charge of the system. Useful for simulations with charged cells.
By default the unit cell is assumed to be neutral (tot_charge=0).
tot_charge=+1 means one electron missing from the system,
tot_charge=-1 means one additional electron, and so on.

In a periodic calculation a compensating jellium background is
inserted to remove divergences if the cell is not neutral.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tot_magnetization -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tot_magnetization</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> -1 [unspecified]
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
total majority spin charge - minority spin charge.
Used to impose a specific total electronic magnetization.
If unspecified then tot_magnetization variable is ignored and
the amount of electronic magnetization is determined during
the self-consistent cycle.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help starting_magnetization -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>starting_magnetization(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
starting spin polarization on atomic type 'i' in a spin
polarized calculation. Values range between -1 (all spins
down for the valence electrons of atom type 'i') to 1
(all spins up). Breaks the symmetry and provides a starting
point for self-consistency. The default value is zero, BUT a
value MUST be specified for AT LEAST one atomic type in spin
polarized calculations, unless you constrain the magnetization
(see "tot_magnetization" and "constrained_magnetization").
Note that if you start from zero initial magnetization, you
will invariably end up in a nonmagnetic (zero magnetization)
state. If you want to start from an antiferromagnetic state,
you may need to define two different atomic species
corresponding to sublattices of the same atomic type.
starting_magnetization is ignored if you are performing a
non-scf calculation, if you are restarting from a previous
run, or restarting from an interrupted run.
If you fix the magnetization with "tot_magnetization",
you should not specify starting_magnetization.
In the spin-orbit case starting with zero
starting_magnetization on all atoms imposes time reversal
symmetry. The magnetization is never calculated and
kept zero (the internal variable domag is .FALSE.).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ecutwfc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ecutwfc</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Status: </em> REQUIRED
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
kinetic energy cutoff (Ry) for wavefunctions
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ecutrho -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ecutrho</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 4 * ecutwfc
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
kinetic energy cutoff (Ry) for charge density and potential
For norm-conserving pseudopotential you should stick to the
default value, you can reduce it by a little but it will
introduce noise especially on forces and stress.
If there are ultrasoft PP, a larger value than the default is
often desirable (ecutrho = 8 to 12 times ecutwfc, typically).
PAW datasets can often be used at 4*ecutwfc, but it depends
on the shape of augmentation charge: testing is mandatory.
The use of gradient-corrected functional, especially in cells
with vacuum, or for pseudopotential without non-linear core
correction, usually requires an higher values of ecutrho
to be accurately converged.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ecutfock -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ecutfock</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> ecutrho
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
kinetic energy cutoff (Ry) for the exact exchange operator in
EXX type calculations. By default this is the same as ecutrho
but in some EXX calculations significant speed-up can be found
by reducing ecutfock, at the expense of some loss in accuracy.
Currently only implemented for the optimized gamma point only
calculations.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {nr1 nr2 nr3} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>nr1, nr2, nr3</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
three-dimensional FFT mesh (hard grid) for charge
density (and scf potential). If not specified
the grid is calculated based on the cutoff for
charge density (see also "ecutrho")
Note: you must specify all three dimensions for this setting to
be used.
         </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
grouphelp {nr1s nr2s nr3s} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>nr1s, nr2s, nr3s</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
three-dimensional mesh for wavefunction FFT and for the smooth
part of charge density ( smooth grid ).
Coincides with nr1, nr2, nr3 if ecutrho = 4 * ecutwfc ( default )
Note: you must specify all three dimensions for this setting to
be used.
         </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
help nosym -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nosym</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if (.TRUE.) symmetry is not used. Note that
- if the k-point grid is provided in input, it is used "as is"
  and symmetry-inequivalent k-points are not generated;
- if the k-point grid is automatically generated, it will
  contain only points in the irreducible BZ for the bravais
  lattice, irrespective of the actual crystal symmetry.
A careful usage of this option can be advantageous
- in low-symmetry large cells, if you cannot afford a k-point
  grid with the correct symmetry
- in MD simulations
- in calculations for isolated atoms
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nosym_evc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nosym_evc</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if(.TRUE.) symmetry is not used but the k-points are
forced to have the symmetry of the Bravais lattice;
an automatically generated k-point grid will contain
all the k-points of the grid and the points rotated by
the symmetries of the Bravais lattice which are not in the
original grid. If available, time reversal is
used to reduce the k-points (and the q =&gt; -q symmetry
is used in the phonon code). To disable also this symmetry set
noinv=.TRUE..
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help noinv -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>noinv</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if (.TRUE.) disable the usage of k =&gt; -k symmetry
(time reversal) in k-point generation
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help no_t_rev -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>no_t_rev</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if (.TRUE.) disable the usage of magnetic symmetry operations
that consist in a rotation + time reversal.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help force_symmorphic -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>force_symmorphic</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if (.TRUE.) force the symmetry group to be symmorphic by disabling
symmetry operations having an associated fractionary translation
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help use_all_frac -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>use_all_frac</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if (.TRUE.) do not discard symmetry operations with an
associated fractionary translation that does not send the
real-space FFT grid into itself. These operations are
incompatible with real-space symmetrization but not with the
new G-space symmetrization. BEWARE: do not use for phonons!
The phonon code still uses real-space symmetrization.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help occupations -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>occupations</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
'smearing':     gaussian smearing for metals
                see variables 'smearing' and 'degauss'

'tetrahedra' :  especially suited for calculation of DOS
                (see P.E. Bloechl, PRB49, 16223 (1994))
                Requires uniform grid of k-points,
                automatically generated (see below)
                Not suitable (because not variational) for
                force/optimization/dynamics calculations

'fixed' :       for insulators with a gap

'from_input' :  The occupation are read from input file,
                card OCCUPATIONS. Option valid only for a
                single k-point, requires "nbnd" to be set
                in input. Occupations should be consistent
                with the value of "tot_charge".
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help one_atom_occupations -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>one_atom_occupations</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
This flag is used for isolated atoms (nat=1) together with
occupations='from_input'. If it is .TRUE., the wavefunctions
are ordered as the atomic starting wavefunctions, independently
from their eigenvalue. The occupations indicate which atomic
states are filled.
The order of the states is written inside the UPF
pseudopotential file.
In the scalar relativistic case:
S -&gt; l=0, m=0
P -&gt; l=1, z, x, y
D -&gt; l=2, r^2-3z^2, xz, yz, xy, x^2-y^2
In the noncollinear magnetic case (with or without spin-orbit),
each group of states is doubled. For instance:
P -&gt; l=1, z, x, y for spin up, l=1, z, x, y for spin down.
Up and down is relative to the direction of the starting
magnetization.
In the case with spin-orbit and time-reversal
(starting_magnetization=0.0) the atomic wavefunctions are
radial functions multiplied by spin-angle functions.
For instance:
P -&gt; l=1, j=1/2, m_j=-1/2,1/2. l=1, j=3/2,
     m_j=-3/2, -1/2, 1/2, 3/2.
In the magnetic case with spin-orbit the atomic wavefunctions
can be forced to be spin-angle functions by setting
starting_spin_angle to .TRUE..
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help starting_spin_angle -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>starting_spin_angle</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
In the spin-orbit case when domag=.TRUE., by default,
the starting wavefunctions are initialized as in scalar
relativistic noncollinear case without spin-orbit.
By setting starting_spin_angle=.TRUE. this behaviour can
be changed and the initial wavefunctions are radial
functions multiplied by spin-angle functions.
When domag=.FALSE. the initial wavefunctions are always
radial functions multiplied by spin-angle functions
independently from this flag.
When lspinorb is .FALSE. this flag is not used.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help degauss -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>degauss</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.D0 Ry
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
value of the gaussian spreading (Ry) for brillouin-zone
integration in metals.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help smearing -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>smearing</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'gaussian'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
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
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nspin -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nspin</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
nspin = 1 :  non-polarized calculation (default)

nspin = 2 :  spin-polarized calculation, LSDA
             (magnetization along z axis)

nspin = 4 :  spin-polarized calculation, noncollinear
             (magnetization in generic direction)
             DO NOT specify nspin in this case;
             specify "noncolin=.TRUE." instead
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help noncolin -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>noncolin</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true. the program will perform a noncollinear calculation.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ecfixed -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ecfixed</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.0
         </li>
<br><li> <em>See: </em> q2sigma
         </li>
<br>
</ul>      
      
}


# ------------------------------------------------------------------------
help qcutz -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>qcutz</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.0
         </li>
<br><li> <em>See: </em> q2sigma
         </li>
<br>
</ul>      
      
}


# ------------------------------------------------------------------------
help q2sigma -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>q2sigma</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
ecfixed, qcutz, q2sigma:  parameters for modified functional to be
used in variable-cell molecular dynamics (or in stress calculation).
"ecfixed" is the value (in Rydberg) of the constant-cutoff;
"qcutz" and "q2sigma" are the height and the width (in Rydberg)
of the energy step for reciprocal vectors whose square modulus
is greater than "ecfixed". In the kinetic energy, G^2 is
replaced by G^2 + qcutz * (1 + erf ( (G^2 - ecfixed)/q2sigma) )
See: M. Bernasconi et al, J. Phys. Chem. Solids 56, 501 (1995)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help input_dft -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>input_dft</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> read from pseudopotential files
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Exchange-correlation functional: eg 'PBE', 'BLYP' etc
See Modules/functionals.f90 for allowed values.
Overrides the value read from pseudopotential files.
Use with care and if you know what you are doing!
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help exx_fraction -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>exx_fraction</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> it depends on the specified functional
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Fraction of EXX for hybrid functional calculations. In the case of
input_dft='PBE0', the default value is 0.25, while for input_dft='B3LYP'
the exx_fraction default value is 0.20.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help screening_parameter -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>screening_parameter</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.106
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
screening_parameter for HSE like hybrid functionals.
See J. Chem. Phys. 118, 8207 (2003)
and J. Chem. Phys. 124, 219906 (2006) for more informations.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help exxdiv_treatment -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>exxdiv_treatment</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> gygi-baldereschi
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Specific for EXX. It selects the kind of approach to be used
for treating the Coulomb potential divergencies at small q vectors.

gygi-baldereschi : appropriate for cubic and quasi-cubic supercells
vcut_spherical : appropriate for cubic and quasi-cubic supercells
vcut_ws : appropriate for strongly anisotropic supercells, see also
          ecutvcut.
none : sets Coulomb potential at G,q=0 to 0.0 (required for GAU-PBE)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help x_gamma_extrapolation -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>x_gamma_extrapolation</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .true.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Specific for EXX. If true, extrapolate the G=0 term of the
potential (see README in examples/EXX_example for more)
Set this to .false. for GAU-PBE.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ecutvcut -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ecutvcut</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.0 Ry
         </li>
<br><li> <em>See: </em> exxdiv_treatment
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Reciprocal space cutoff for correcting
Coulomb potential divergencies at small q vectors.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {nqx1 nqx2 nqx3} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>nqx1, nqx2, nqx3</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
three-dimensional mesh for q (k1-k2) sampling of
the Fock operator (EXX). Can be smaller than
the number of k-points.

Currently this defaults to the size of the k-point mesh used.
 In QE =&lt; 5.0.2 it defaulted to nqx1=nqx2=nqx3=1.
         </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
help lda_plus_u -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lda_plus_u</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Status: </em>
DFT+U (formerly known as LDA+U) currently works only for
a few selected elements. Modify flib/set_hubbard_l.f90 and
PW/src/tabd.f90 if you plan to use DFT+U with an element that
is not configured there.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Specify lda_plus_u = .TRUE. to enable DFT+U calculations
See: Anisimov, Zaanen, and Andersen, PRB 44, 943 (1991);
     Anisimov et al., PRB 48, 16929 (1993);
     Liechtenstein, Anisimov, and Zaanen, PRB 52, R5467 (1994).
You must specify, for each species with a U term, the value of
U and (optionally) alpha, J of the Hubbard model (all in eV):
see lda_plus_u_kind, Hubbard_U, Hubbard_alpha, Hubbard_J
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lda_plus_u_kind -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lda_plus_u_kind</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Specifies the type of DFT+U calculation:
                  0   simplified version of Cococcioni and de Gironcoli,
                      PRB 71, 035105 (2005), using Hubbard_U
                  1   rotationally invariant scheme of Liechtenstein et al.,
                      using Hubbard_U and Hubbard_J
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help Hubbard_U -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>Hubbard_U(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.D0 for all species
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Hubbard_U(i): U parameter (eV) for species i, DFT+U calculation
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help Hubbard_J0 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>Hubbard_J0(i), i=1,ntype</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.D0 for all species
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Hubbard_J0(i): J0 parameter (eV) for species i, DFT+U+J calculation,
see PRB 84, 115108 (2011) for details.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help Hubbard_alpha -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>Hubbard_alpha(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.D0 for all species
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Hubbard_alpha(i) is the perturbation (on atom i, in eV)
used to compute U with the linear-response method of
Cococcioni and de Gironcoli, PRB 71, 35105 (2005)
(only for lda_plus_u_kind=0)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help Hubbard_beta -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>Hubbard_beta(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.D0 for all species
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Hubbard_beta(i) is the perturbation (on atom i, in eV)
used to compute J0 with the linear-response method of
Cococcioni and de Gironcoli, PRB 71, 35105 (2005)
(only for lda_plus_u_kind=0). See also
PRB 84, 115108 (2011).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help U_projection_type -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>U_projection_type</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'atomic'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Only active when lda_plus_U is .true., specifies the type
of projector on localized orbital to be used in the DFT+U
scheme.

Currently available choices:
'atomic': use atomic wfc's (as they are) to build the projector

'ortho-atomic': use Lowdin orthogonalized atomic wfc's

'norm-atomic':  Lowdin normalization of atomic wfc. Keep in mind:
                atomic wfc are not orthogonalized in this case.
                This is a "quick and dirty" trick to be used when
                atomic wfc from the pseudopotential are not
                normalized (and thus produce occupation whose
                value exceeds unity). If orthogonalized wfc are
                not needed always try 'atomic' first.

'file':         use the information from file "prefix".atwfc that must
                have been generated previously, for instance by pmw.x
                (see PP/src/poormanwannier.f90 for details).

'pseudo':       use the pseudopotential projectors. The charge density
                outside the atomic core radii is excluded.
                N.B.: for atoms with +U, a pseudopotential with the
                all-electron atomic wavefunctions is required (i.e.,
                as generated by ld1.x with lsave_wfc flag).

NB: forces and stress currently implemented only for the
'atomic' and 'pseudo' choice.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help edir -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>edir</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The direction of the electric field or dipole correction is
parallel to the bg(:,edir) reciprocal lattice vector, so the
potential is constant in planes defined by FFT grid points;
edir = 1, 2 or 3. Used only if tefield is .TRUE.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help emaxpos -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>emaxpos</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.5D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Position of the maximum of the saw-like potential along crystal
axis "edir", within the  unit cell (see below), 0 &lt; emaxpos &lt; 1
Used only if tefield is .TRUE.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help eopreg -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>eopreg</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.1D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Zone in the unit cell where the saw-like potential decreases.
( see below, 0 &lt; eopreg &lt; 1 ). Used only if tefield is .TRUE.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help eamp -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>eamp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.001 a.u.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Amplitude of the electric field, in ***Hartree*** a.u.;
1 a.u. = 51.4220632*10^10 V/m). Used only if tefield=.TRUE.
The saw-like potential increases with slope "eamp" in the
region from (emaxpos+eopreg-1) to (emaxpos), then decreases
to 0 until (emaxpos+eopreg), in units of the crystal
vector "edir". Important: the change of slope of this
potential must be located in the empty region, or else
unphysical forces will result.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help angle1 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>angle1(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The angle expressed in degrees between the initial
magnetization and the z-axis. For noncollinear calculations
only; index i runs over the atom types.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help angle2 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>angle2(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The angle expressed in degrees between the projection
of the initial magnetization on x-y plane and the x-axis.
For noncollinear calculations only.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help constrained_magnetization -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>constrained_magnetization</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'none'
         </li>
<br><li> <em>See: </em> lambda, fixed_magnetization
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Used to perform constrained calculations in magnetic systems.
Currently available choices:

'none':
         no constraint

'total':
         total magnetization is constrained by
         adding a penalty functional to the total energy:

         LAMBDA * SUM_{i} ( magnetization(i) - fixed_magnetization(i) )**2

         where the sum over i runs over the three components of
         the magnetization. Lambda is a real number (see below).
         Noncolinear case only. Use "tot_magnetization" for LSDA

'atomic':
         atomic magnetization are constrained to the defined
         starting magnetization adding a penalty:

         LAMBDA * SUM_{i,itype} ( magnetic_moment(i,itype) - mcons(i,itype) )**2

         where i runs over the cartesian components (or just z
         in the collinear case) and itype over the types (1-ntype).
         mcons(:,:) array is defined from starting_magnetization,
         (and angle1, angle2 in the non-collinear case). lambda is
         a real number

'total direction':
          the angle theta of the total magnetization
          with the z axis (theta = fixed_magnetization(3))
          is constrained:

          LAMBDA * ( arccos(magnetization(3)/mag_tot) - theta )**2

          where mag_tot is the modulus of the total magnetization.

'atomic direction':
          not all the components of the atomic
          magnetic moment are constrained but only the cosine
          of angle1, and the penalty functional is:

          LAMBDA * SUM_{itype} ( mag_mom(3,itype)/mag_mom_tot - cos(angle1(ityp)) )**2

N.B.: symmetrization may prevent to reach the desired orientation
      of the magnetization. Try not to start with very highly symmetric
      configurations or use the nosym flag (only as a last remedy)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fixed_magnetization -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>fixed_magnetization(i), i=1,3</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.d0
         </li>
<br><li> <em>See: </em> constrained_magnetization
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
total magnetization vector (x,y,z components) to be kept
fixed when constrained_magnetization='total'
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lambda -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lambda</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.d0
         </li>
<br><li> <em>See: </em> constrained_magnetization
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
parameter used for constrained_magnetization calculations
N.B.: if the scf calculation does not converge, try to reduce lambda
      to obtain convergence, then restart the run with a larger lambda
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help report -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>report</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
It is the number of iterations after which the program
write all the atomic magnetic moments.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lspinorb -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lspinorb</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .TRUE. the noncollinear code can use a pseudopotential with
spin-orbit.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help assume_isolated -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>assume_isolated</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'none'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Used to perform calculation assuming the system to be
isolated (a molecule or a cluster in a 3D supercell).

Currently available choices:

'none' (default): regular periodic calculation w/o any correction.

'makov-payne', 'm-p', 'mp' : the Makov-Payne correction to the
         total energy is computed. An estimate of the vacuum
         level is also calculated so that eigenvalues can be
         properly aligned. ONLY FOR CUBIC SYSTEMS (ibrav=1,2,3)
         Theory:
         G.Makov, and M.C.Payne,
         "Periodic boundary conditions in ab initio
         calculations" , Phys.Rev.B 51, 4014 (1995)

'martyna-tuckerman', 'm-t', 'mt' : Martyna-Tuckerman correction
         to both total energy and scf potential. Adapted from:
         G.J. Martyna, and M.E. Tuckerman,
         "A reciprocal space based method for treating long
         range interactions in ab-initio and force-field-based
         calculation in clusters", J.Chem.Phys. 110, 2810 (1999)

'esm' :  Effective Screening Medium Method.
         For polarized or charged slab calculation, embeds
         the simulation cell within an effective semi-
         infinite medium in the perpendicular direction
         (along z). Embedding regions can be vacuum or
         semi-infinite metal electrodes (use 'esm_bc' to
         choose boundary conditions). If between two
         electrodes, an optional electric field
         ('esm_efield') may be applied. Method described in
         M. Otani and O. Sugino, "First-principles
         calculations of charged surfaces and interfaces:
         A plane-wave nonrepeated slab approach," PRB 73,
         115407 (2006).
         NB: Requires cell with a_3 lattice vector along z,
         normal to the xy plane, with the slab centered
         around z=0. Also requires symmetry checking to be
         disabled along z, either by setting 'nosym' = .TRUE.
         or by very slight displacement (i.e., 5e-4 a.u.)
         of the slab along z.
         See 'esm_bc', 'esm_efield', 'esm_w', 'esm_nfit'.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help esm_bc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>esm_bc</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'pbc'
         </li>
<br><li> <em>See: </em> assume_isolated
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If assume_isolated = 'esm', determines the boundary
conditions used for either side of the slab.

Currently available choices:

'pbc' (default): regular periodic calculation (no ESM).

'bc1' : Vacuum-slab-vacuum (open boundary conditions)

'bc2' : Metal-slab-metal (dual electrode configuration).
        See also 'esm_efield'.

'bc3' : Vacuum-slab-metal
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help esm_w -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>esm_w</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.d0
         </li>
<br><li> <em>See: </em> assume_isolated
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If assume_isolated = 'esm', determines the position offset
[in a.u.] of the start of the effective screening region,
measured relative to the cell edge. (ESM region begins at
z = +/- [L_z/2 + esm_w] ).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help esm_efield -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>esm_efield</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.d0
         </li>
<br><li> <em>See: </em> assume_isolated, esm_bc
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If assume_isolated = 'esm' and esm_bc = 'bc2', gives the
magnitude of the electric field [Ry/a.u.] to be applied
between semi-infinite ESM electrodes.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help esm_nfit -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>esm_nfit</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 4
         </li>
<br><li> <em>See: </em> assume_isolated
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If assume_isolated = 'esm', gives the number of z-grid points
for the polynomial fit along the cell edge.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help vdw_corr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>vdw_corr</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'none'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Type of Van der Waals correction. Allowed values:

   'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d': semiempirical Grimme's DFT-D2.
    Optional variables: "london_s6", "london_rcut"
    S. Grimme, J. Comp. Chem. 27, 1787 (2006),
    V. Barone et al., J. Comp. Chem. 30, 934 (2009).

    'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler': Tkatchenko-Scheffler
     dispersion corrections with first-principle derived C6 coefficients
     (implemented in CP only). Optional variables: "ts_vdw_econv_thr", "ts_vdw_isolated"
     See A. Tkatchenko and M. Scheffler, Phys. Rev. Lett. 102, 073005 (2009)

    'XDM', 'xdm': Exchange-hole dipole-moment model. Optional variables: "xdm_a1", "xdm_a2"
     A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007)
         A. Otero de la Roza, E. R. Johnson, J. Chem. Phys. 136, 174109 (2012)

Note that non-local functionals (eg vdw-DF) are NOT specified here but in "input_dft"
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help london -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>london</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
OBSOLESCENT, same as vdw_corr='DFT-D'
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help london_s6 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>london_s6</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.75
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
global scaling parameter for DFT-D. Default is good for PBE.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help london_rcut -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>london_rcut</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 200
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
cutoff radius (a.u.) for dispersion interactions
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help xdm -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>xdm</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
OBSOLESCENT, same as vdw_corr='xdm'
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help xdm_a1 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>xdm_a1</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.6836
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Damping function parameter a1 (adimensional). This value should change
with the exchange-correlation functional. The default corresponds to
PW86PBE.
For other functionals, see:
   http://gatsby.ucmerced.edu/wiki/XDM_damping_function_parameters
   A. Otero de la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help xdm_a2 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>xdm_a2</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.5045
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Damping function parameter a2 (angstrom). This value should change
with the exchange-correlation functional. The default corresponds to
PW86PBE.
For other functionals, see:
   http://gatsby.ucmerced.edu/wiki/XDM_damping_function_parameters
   A. Otero de la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help space_group -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>space_group</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The number of the space group of the crystal, as given
                in the International Tables of Crystallography A (ITA).
                This allows to give in input only the inequivalent atomic
                positions. The positions of all the symmetry equivalent atoms
                are calculated by the code. Used only when the atomic positions
                are of type crystal_sg.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help uniqueb -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>uniqueb</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Used only for monoclinic lattices. If .TRUE. the b
                 unique ibrav (-12 or -13) are used, and symmetry
                 equivalent positions are chosen assuming that the
                 two fold axis or the mirror normal is parallel to the
                 b axis. If .FALSE. it is parallel to the c axis.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help origin_choice -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>origin_choice</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Used only for space groups that in the ITA allow
                 the use of two different origins. origin_choice=1,
                 means the first origin, while origin_choice=2 is the
                 second origin.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rhombohedral -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rhombohedral</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .TRUE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Used only for rhombohedral space groups.
                 When .TRUE. the coordinates of the inequivalent atoms are
                 given with respect to the rhombohedral axes, when .FALSE.
                 the coordinates of the inequivalent atoms are given with
                 respect to the hexagonal axes. They are converted internally
                 to the rhombohedral axes and ibrav=5 is used in both cases.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help electron_maxstep -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>electron_maxstep</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 100
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
maximum number of iterations in a scf step
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help scf_must_converge -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>scf_must_converge</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .TRUE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .false. do not stop molecular dynamics or ionic relaxation
when electron_maxstep is reached. Use with care.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help conv_thr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>conv_thr</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.D-6
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Convergence threshold for selfconsistency:
   estimated energy error &lt; conv_thr
(note that conv_thr is extensive, like the total energy).
For non-self-consistent calculations, conv_thr is used
to set the default value of the threshold (ethr) for
iterative diagonalizazion: see diago_thr_init
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help adaptive_thr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>adaptive_thr</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE. this turns on the use of an adaptive conv_thr for
the inner scf loops when using EXX.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help conv_thr_init -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>conv_thr_init</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.D-3
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
When adaptive_thr = .TRUE. this is the convergence threshold
used for the first scf cycle.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help conv_thr_multi -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>conv_thr_multi</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.D-1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
When adaptive_thr = .TRUE. the convergence threshold for
each scf cycle is given by:
max( conv_thr, conv_thr_multi * dexx )
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help mixing_mode -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>mixing_mode</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'plain'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
'plain' :    charge density Broyden mixing

'TF' :       as above, with simple Thomas-Fermi screening
            (for highly homogeneous systems)

'local-TF':  as above, with local-density-dependent TF screening
             (for highly inhomogeneous systems)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help mixing_beta -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>mixing_beta</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.7D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
mixing factor for self-consistency
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help mixing_ndim -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>mixing_ndim</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 8
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
number of iterations used in mixing scheme.
If you are tight with memory, you may reduce it to 4 or so.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help mixing_fixed_ns -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>mixing_fixed_ns</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
For DFT+U : number of iterations with fixed ns ( ns is the
  atomic density appearing in the Hubbard term ).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help diagonalization -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>diagonalization</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'david'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
'david':  Davidson iterative diagonalization with overlap matrix
          (default). Fast, may in some rare cases fail.

'cg' :    conjugate-gradient-like band-by-band diagonalization
          Typically slower than 'david' but it uses less memory
          and is more robust (it seldom fails)

'cg-serial', 'david-serial': obsolete, use "-ndiag 1 instead"
          The subspace diagonalization in Davidson is performed
          by a fully distributed-memory parallel algorithm on
          4 or more processors, by default. The allocated memory
          scales down with the number of procs. Procs involved
          in diagonalization can be changed with command-line
          option "-ndiag N". On multicore CPUs it is often
          convenient to let just one core per CPU to work
          on linear algebra.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ortho_para -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ortho_para</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Status: </em> OBSOLETE: use command-line option " -ndiag XX" instead
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>

         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help diago_thr_init -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>diago_thr_init</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Convergence threshold (ethr) for iterative diagonalization
(the check is on eigenvalue convergence).
For scf calculations: default is 1.D-2 if starting from a
superposition of atomic orbitals; 1.D-5 if starting from a
charge density. During self consistency the threshold
is automatically reduced (but never below 1.D-13) when
approaching convergence.
For non-scf calculations: default is (conv_thr/N elec)/10.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help diago_cg_maxiter -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>diago_cg_maxiter</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
For conjugate gradient diagonalization:
max number of iterations
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help diago_david_ndim -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>diago_david_ndim</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 4
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
For Davidson diagonalization: dimension of workspace
(number of wavefunction packets, at least 2 needed).
A larger value may yield a somewhat faster algorithm
but uses more memory. The opposite holds for smaller values.
Try diago_david_ndim=2 if you are tight on memory or if
your job is large: the speed penalty is often negligible
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help diago_full_acc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>diago_full_acc</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE. all the empty states are diagonalized at the same level
of accuracy of the occupied ones. Otherwise the empty states are
diagonalized using a larger threshold (this should not affect
total energy, forces, and other ground-state properties).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help efield -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>efield</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Amplitude of the finite electric field (in Ry a.u.;
1 a.u. = 36.3609*10^10 V/m). Used only if lelfield=.TRUE.
and if k-points (K_POINTS card) are not automatic.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help efield_cart -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>efield_cart(i), i=1,3</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> (0.D0, 0.D0, 0.D0)
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Finite electric field (in Ry a.u.=36.3609*10^10 V/m) in
cartesian axis. Used only if lelfield=.TRUE. and if
k-points (K_POINTS card) are automatic.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help startingpot -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>startingpot</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
'atomic': starting potential from atomic charge superposition
          ( default for scf, *relax, *md )

'file'  : start from existing "charge-density.xml" file in the
          directory specified by variables "prefix" and "outdir"
          For nscf and bands calculation this is the default
          and the only sensible possibility.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help startingwfc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>startingwfc</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'atomic+random'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
'atomic': start from superposition of atomic orbitals
          If not enough atomic orbitals are available,
          fill with random numbers the remaining wfcs
          The scf typically starts better with this option,
          but in some high-symmetry cases one can "loose"
          valence states, ending up in the wrong ground state.

'atomic+random': as above, plus a superimposed "randomization"
          of atomic orbitals. Prevents the "loss" of states
          mentioned above.

'random': start from random wfcs. Slower start of scf but safe.
          It may also reduce memory usage in conjunction with
          diagonalization='cg'

'file':   start from an existing wavefunction file in the
          directory specified by variables "prefix" and "outdir"
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tqr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tqr</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true., use the real-space algorithm for augmentation
charges in ultrasoft pseudopotentials.
Must faster execution of ultrasoft-related calculations,
but numerically less accurate than the default algorithm.
Use with care and after testing!
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ion_dynamics -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ion_dynamics</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Specify the type of ionic dynamics.

For different type of calculation different possibilities are
allowed and different default values apply:

CASE ( calculation = 'relax' )
    'bfgs' :   (default)   use BFGS quasi-newton algorithm,
                           based on the trust radius procedure,
                           for structural relaxation
    'damp' :               use damped (quick-min Verlet)
                           dynamics for structural relaxation
                           Can be used for constrained
                           optimisation: see CONSTRAINTS card

CASE ( calculation = 'md' )
    'verlet' : (default)   use Verlet algorithm to integrate
                           Newton's equation. For constrained
                           dynamics, see CONSTRAINTS card
    'langevin'             ion dynamics is over-damped Langevin
    'langevin-smc'         over-damped Langevin with Smart Monte Carlo:
                           see R.J.Rossky, JCP, 69, 4628(1978)


CASE ( calculation = 'vc-relax' )
    'bfgs' :   (default)   use BFGS quasi-newton algorithm;
                           cell_dynamics must be 'bfgs' too
    'damp' :               use damped (Beeman) dynamics for
                           structural relaxation
CASE ( calculation = 'vc-md' )
    'beeman' : (default)   use Beeman algorithm to integrate
                           Newton's equation
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ion_positions -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ion_positions</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'default'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
'default '  : if restarting, use atomic positions read from the
              restart file; in all other cases, use atomic
              positions from standard input.

'from_input' : restart the simulation with atomic positions read
              from standard input, even if restarting.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help pot_extrapolation -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>pot_extrapolation</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'atomic'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
   Used to extrapolate the potential from preceding ionic steps.

   'none'        :  no extrapolation

   'atomic'      :  extrapolate the potential as if it was a sum of
                    atomic-like orbitals

   'first_order' :  extrapolate the potential with first-order
                    formula

   'second_order':  as above, with second order formula

Note: 'first_order' and 'second-order' extrapolation make sense
only for molecular dynamics calculations
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help wfc_extrapolation -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>wfc_extrapolation</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'none'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
    Used to extrapolate the wavefunctions from preceding ionic steps.

   'none'        :  no extrapolation

   'first_order' :  extrapolate the wave-functions with first-order
                    formula.

   'second_order':  as above, with second order formula.

Note: 'first_order' and 'second-order' extrapolation make sense
only for molecular dynamics calculations
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help remove_rigid_rot -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>remove_rigid_rot</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
This keyword is useful when simulating the dynamics and/or the
thermodynamics of an isolated system. If set to true the total
torque of the internal forces is set to zero by adding new forces
that compensate the spurious interaction with the periodic
images. This allows for the use of smaller supercells.

BEWARE: since the potential energy is no longer consistent with
the forces (it still contains the spurious interaction with the
repeated images), the total energy is not conserved anymore.
However the dynamical and thermodynamical properties should be
in closer agreement with those of an isolated system.
Also the final energy of a structural relaxation will be higher,
but the relaxation itself should be faster.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ion_temperature -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ion_temperature</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'not_controlled'
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
'rescaling'   control ionic temperature via velocity rescaling
              (first method) see parameters "tempw", "tolp", and
              "nraise" (for VC-MD only). This rescaling method
              is the only one currently implemented in VC-MD

'rescale-v'   control ionic temperature via velocity rescaling
              (second method) see parameters "tempw" and "nraise"

'rescale-T'   control ionic temperature via velocity rescaling
              (third method) see parameter "delta_t"

'reduce-T'    reduce ionic temperature every "nraise" steps
              by the (negative) value "delta_t"

'berendsen'   control ionic temperature using "soft" velocity
              rescaling - see parameters "tempw" and "nraise"

'andersen'    control ionic temperature using Andersen thermostat
              see parameters "tempw" and "nraise"

'initial'     initialize ion velocities to temperature "tempw"
              and leave uncontrolled further on

'not_controlled' (default) ionic temperature is not controlled
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tempw -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tempw</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 300.D0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Starting temperature (Kelvin) in MD runs
target temperature for most thermostats.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tolp -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tolp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 100.D0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Tolerance for velocity rescaling. Velocities are rescaled if
the run-averaged and target temperature differ more than tolp.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help delta_t -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>delta_t</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.D0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if ion_temperature='rescale-T':
       at each step the instantaneous temperature is multiplied
       by delta_t; this is done rescaling all the velocities.

if ion_temperature='reduce-T':
       every 'nraise' steps the instantaneous temperature is
       reduced by -delta_T (i.e. delta_t &lt; 0 is added to T)

The instantaneous temperature is calculated at the end of
every ionic move and BEFORE rescaling. This is the temperature
reported in the main output.

For delta_t &lt; 0, the actual average rate of heating or cooling
should be roughly C*delta_t/(nraise*dt) (C=1 for an
ideal gas, C=0.5 for a harmonic solid, theorem of energy
equipartition between all quadratic degrees of freedom).
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nraise -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nraise</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if ion_temperature='reduce-T':
       every 'nraise' steps the instantaneous temperature is
       reduced by -delta_T (.e. delta_t is added to the temperature)

if ion_temperature='rescale-v':
       every 'nraise' steps the average temperature, computed from
       the last nraise steps, is rescaled to tempw

if ion_temperature='rescaling' and calculation='vc-md':
       every 'nraise' steps the instantaneous temperature
       is rescaled to tempw

if ion_temperature='berendsen':
       the "rise time" parameter is given in units of the time step:
       tau = nraise*dt, so dt/tau = 1/nraise

if ion_temperature='andersen':
       the "collision frequency" parameter is given as nu=1/tau
       defined above, so nu*dt = 1/nraise
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help refold_pos -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>refold_pos</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
This keyword applies only in the case of molecular dynamics or
damped dynamics. If true the ions are refolded at each step into
the supercell.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help upscale -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>upscale</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 100.D0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Max reduction factor for conv_thr during structural optimization
conv_thr is automatically reduced when the relaxation
approaches convergence so that forces are still accurate,
but conv_thr will not be reduced to less that
conv_thr / upscale.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help bfgs_ndim -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>bfgs_ndim</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Number of old forces and displacements vectors used in the
PULAY mixing of the residual vectors obtained on the basis
of the inverse hessian matrix given by the BFGS algorithm.
When bfgs_ndim = 1, the standard quasi-Newton BFGS method is
used.
(bfgs only)
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help trust_radius_max -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>trust_radius_max</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.8D0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Maximum ionic displacement in the structural relaxation.
(bfgs only)
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help trust_radius_min -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>trust_radius_min</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.D-3
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Minimum ionic displacement in the structural relaxation
BFGS is reset when trust_radius &lt; trust_radius_min.
(bfgs only)
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help trust_radius_ini -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>trust_radius_ini</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.5D0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Initial ionic displacement in the structural relaxation.
(bfgs only)
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help w_1 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>w_1</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.01D0
            </li>
<br><li> <em>See: </em> w_2
            </li>
<br>
</ul>      
      
}


# ------------------------------------------------------------------------
help w_2 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>w_2</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.5D0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Parameters used in line search based on the Wolfe conditions.
(bfgs only)
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help cell_dynamics -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>cell_dynamics</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Specify the type of dynamics for the cell.
For different type of calculation different possibilities
are allowed and different default values apply:

CASE ( calculation = 'vc-relax' )
  'none':    no dynamics
  'sd':      steepest descent ( not implemented )
  'damp-pr': damped (Beeman) dynamics of the Parrinello-Rahman
             extended lagrangian
  'damp-w':  damped (Beeman) dynamics of the new Wentzcovitch
             extended lagrangian
  'bfgs':    BFGS quasi-newton algorithm (default)
             ion_dynamics must be 'bfgs' too
CASE ( calculation = 'vc-md' )
  'none':    no dynamics
  'pr':      (Beeman) molecular dynamics of the Parrinello-Rahman
             extended lagrangian
  'w':       (Beeman) molecular dynamics of the new Wentzcovitch
             extended lagrangian
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help press -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>press</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Target pressure [KBar] in a variable-cell md or relaxation run.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help wmass -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>wmass</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em>
0.75*Tot_Mass/pi**2 for Parrinello-Rahman MD;
0.75*Tot_Mass/pi**2/Omega**(2/3) for Wentzcovitch MD
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Fictitious cell mass [amu] for variable-cell simulations
(both 'vc-md' and 'vc-relax')
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help cell_factor -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>cell_factor</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.2D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Used in the construction of the pseudopotential tables.
It should exceed the maximum linear contraction of the
cell during a simulation.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help press_conv_thr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>press_conv_thr</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.5D0 Kbar
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Convergence threshold on the pressure for variable cell
relaxation ('vc-relax' : note that the other convergence
thresholds for ionic relaxation apply as well).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help cell_dofree -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>cell_dofree</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'all'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Select which of the cell parameters should be moved:

all     = all axis and angles are moved
x       = only the x component of axis 1 (v1_x) is moved
y       = only the y component of axis 2 (v2_y) is moved
z       = only the z component of axis 3 (v3_z) is moved
xy      = only v1_x and v2_y are moved
xz      = only v1_x and v3_z are moved
yz      = only v2_y and v3_z are moved
xyz     = only v1_x, v2_y, v3_z are moved
shape   = all axis and angles, keeping the volume fixed
volume  = the volume changes, keeping all angles fixed (i.e. only celldm(1) changes)
2Dxy    = only x and y components are allowed to change
2Dshape = as above, keeping the area in xy plane fixed

BEWARE: if axis are not orthogonal, some of these options do not
 work (symmetry is broken). If you are not happy with them,
 edit subroutine init_dofree in file Modules/cell_base.f90
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help atomic_species -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variable: </em><big><b>X</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
label of the atom. Acceptable syntax:
chemical symbol X (1 or 2 characters, case-insensitive)
or chemical symbol plus a number or a letter, as in
"Xn" (e.g. Fe1) or "X_*" or "X-*" (e.g. C1, C_h;
max total length cannot exceed 3 characters)
                  </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>Mass_X</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
mass of the atomic species [amu: mass of C = 12]
Used only when performing Molecular Dynamics run
or structural optimization runs using Damped MD.
Not actually used in all other cases (but stored
in data files, so phonon calculations will use
these values unless other values are provided)
                  </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>PseudoPot_X</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
File containing PP for this species.

The pseudopotential file is assumed to be in the new UPF format.
If it doesn't work, the pseudopotential format is determined by
the file name:

*.vdb or *.van     Vanderbilt US pseudopotential code
*.RRKJ3            Andrea Dal Corso's code (old format)
none of the above  old PWscf norm-conserving format
                  </pre></blockquote>
</ul>   
    
}


# ------------------------------------------------------------------------
help ATOMIC_POSITIONS_flags -helpfmt helpdoc -helptext {
      <h2>Description of ATOMIC_POSITIONS card's flags</h2><pre>
alat    : atomic positions are in cartesian coordinates, in
          units of the lattice parameter (either celldm(1)
          or A). If no option is specified, 'alat' is assumed;
          not specifying units is DEPRECATED and will no
          longer be allowed in the future

bohr    : atomic positions are in cartesian coordinate,
          in atomic units (i.e. Bohr radii)

angstrom: atomic positions are in cartesian coordinates,
          in Angstrom

crystal : atomic positions are in crystal coordinates, i.e.
          in relative coordinates of the primitive lattice
          vectors as defined either in card CELL_PARAMETERS
          or via the ibrav + celldm / a,b,c... variables

crystal_sg : atomic positions are in crystal coordinates, i.e.
          in relative coordinates of the primitive lattice.
          This option differs from the previous one because
          in this case only the symmetry inequivalent atoms
          are given. The variable space_group must indicate
          the space group number used to find the symmetry
          equivalent atoms. The other variables that control
          this option are uniqueb, origin_choice, and
          rhombohedral.
         </pre>
      
}


# ------------------------------------------------------------------------
help atomic_coordinates -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variable: </em><big><b>X</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> label of the atom as specified in ATOMIC_SPECIES
                        </pre></blockquote>
</ul><ul>
<li> <em>Variables: </em><big><b>x, y, z</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
atomic positions

NOTE: each atomic coordinate can also be specified as a simple algebraic expression.
      To be interpreted correctly expression must NOT contain any blank
      space and must NOT start with a "+" sign. The available expressions are:

        + (plus), - (minus), / (division), * (multiplication), ^ (power)

     All numerical constants included are considered as double-precision numbers;
     i.e. 1/2 is 0.5, not zero. Other functions, such as sin, sqrt or exp are
     not available, although sqrt can be replaced with ^(1/2).

    Example:
                 C  1/3   1/2*3^(-1/2)   0

    is equivalent to

            C  0.333333  0.288675  0.000000

    Please note that this feature is NOT supported by XCrysDen (which will
    display a wrong structure, or nothing at all).
    When atomic positions are of type crystal_sg coordinates can be given
    in the following three forms (Wyckoff positions):
    C  1a
    C  8g   x
    C  24m  x y
    The first form must be used when the Wyckoff letter determines uniquely
    all three coordinates, the second form when the Wyckoff letter
    and one coordinate are needed, the third form when one letter and two
    coordinates are needed.
    The forms:
    C 8g  x  x  x
    C 24m x  x  y
    are not allowed, but
    C x x x
    C x x y
    are correct.
                        </pre></blockquote>
</ul><ul>
<li> <em>Variables: </em><big><b>if_pos(1), if_pos(2), if_pos(3)</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
                           </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
component i of the force for this atom is multiplied by if_pos(i),
which must be either 0 or 1.  Used to keep selected atoms and/or
selected components fixed in MD dynamics or
structural optimization run.

With crystal_sg atomic coordinates the constraints are copied in all equivalent
atoms.
                           </pre></blockquote>
</ul>   
    
}


# ------------------------------------------------------------------------
help K_POINTS_flags -helpfmt helpdoc -helptext {
      <h2>Description of K_POINTS card's flags</h2><pre>
 tpiba    : read k-points in cartesian coordinates,
            in units of 2 pi/a (default)

 automatic: automatically generated uniform grid of k-points, i.e,
            generates ( nk1, nk2, nk3 ) grid with ( sk1, sk2, sk3 ) offset.
            nk1, nk2, nk3 as in Monkhorst-Pack grids
            k1, k2, k3 must be 0 ( no offset ) or 1 ( grid displaced
            by half a grid step in the corresponding direction )
            BEWARE: only grids having the full symmetry of the crystal
            work with tetrahedra. Some grids with offset may not work.

 crystal  : read k-points in crystal coordinates, i.e. in relative
            coordinates of the reciprocal lattice vectors

 gamma    : use k = 0 (no need to list k-point specifications after card)
            In this case wavefunctions can be chosen as real,
            and specialized subroutines optimized for calculations
            at the gamma point are used (memory and cpu requirements
            are reduced by approximately one half).

 tpiba_b  : Used for band-structure plots.
            k-points are in units of  2 pi/a.
            nks points specify nks-1 lines in reciprocal space.
            Every couple of points identifies the initial and
            final point of a line. pw.x generates N
            intermediate points of the line where N is the
            weight of the first point.

 crystal_b: as tpiba_b, but k-points are in crystal coordinates.

 tpiba_c  : Used for band-structure contour plots.
            k-points are in units of  2 pi/a. nks must be 3.
            3 k-points k_0, k_1, and k_2 specify a rectangle
            in reciprocal space of vertices k_0, k_1, k_2,
            k_1 + k_2 - k_0: k_0 + \alpha (k_1-k_0)+
            \beta (k_2-k_0) with 0&lt;\alpha,\beta &lt; 1.
            The code produces a uniform mesh n1 x n2
            k points in this rectangle. n1 and n2 are
            the weights of k_1 and k_2. The weight of k_0
            is not used.

crystal_c: as tpiba_c, but k-points are in crystal coordinates.
         </pre>
      
}


# ------------------------------------------------------------------------
help nks -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nks</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Number of supplied special k-points.
                     </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help kpoints -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>xk_x, xk_y, xk_z, wk</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Special k-points (xk_x/y/z) in the irreducible Brillouin Zone
(IBZ) of the lattice (with all symmetries) and weights (wk)
See the literature for lists of special points and
the corresponding weights.

If the symmetry is lower than the full symmetry
of the lattice, additional points with appropriate
weights are generated. Notice that such procedure
assumes that ONLY k-points in the IBZ are provided in input

In a non-scf calculation, weights do not affect the results.
If you just need eigenvalues and eigenvectors (for instance,
for a band-structure plot), weights can be set to any value
(for instance all equal to 1).
                        </pre></blockquote>
</ul>   
    
}


# ------------------------------------------------------------------------
grouphelp {nk1 nk2 nk3} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>nk1, nk2, nk3</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
These parameters specify the k-point grid
(nk1 x nk2 x nk3) as in Monkhorst-Pack grids.
                     </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
grouphelp {sk1 sk2 sk3} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>sk1, sk2, sk3</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The grid offsets;  sk1, sk2, sk3 must be
0 ( no offset ) or 1 ( grid displaced by
half a grid step in the corresponding direction ).
                     </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
help CELL_PARAMETERS_flags -helpfmt helpdoc -helptext {
      <h2>Description of CELL_PARAMETERS card's flags</h2><pre>
'bohr'/'angstrom': lattice vectors in bohr radii / angstrom.
   In this case the lattice parameter alat = sqrt(v1*v1).
'alat' / nothing specified: lattice vectors in units of the
lattice parameter (either celldm(1) or a). Not specifying
units is DEPRECATED and will not be allowed in the future.
If nothing specified and no lattice parameter specified,
'bohr' is assumed - DEPRECATED, will no longer be allowed
         </pre>
      
}


# ------------------------------------------------------------------------
help lattice -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>v1, v2, v3</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Crystal lattice vectors (in cartesian axis):
    v1(1)  v1(2)  v1(3)    ... 1st lattice vector
    v2(1)  v2(2)  v2(3)    ... 2nd lattice vector
    v3(1)  v3(2)  v3(3)    ... 3rd lattice vector
                  </pre></blockquote>
</ul>   
    
}


# ------------------------------------------------------------------------
help nconstr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nconstr</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Number of constraints.
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help constr_tol -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>constr_tol</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Tolerance for keeping the constraints satisfied.
                  </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help constraints_table -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variable: </em><big><b>constr_type</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Type of constrain :

'type_coord'      : constraint on global coordination-number, i.e. the
                    average number of atoms of type B surrounding the
                    atoms of type A. The coordination is defined by
                    using a Fermi-Dirac.
                    (four indexes must be specified).

'atom_coord'      : constraint on local coordination-number, i.e. the
                    average number of atoms of type A surrounding a
                    specific atom. The coordination is defined by
                    using a Fermi-Dirac.
                    (four indexes must be specified).

'distance'        : constraint on interatomic distance
                    (two atom indexes must be specified).

'planar_angle'    : constraint on planar angle
                    (three atom indexes must be specified).

'torsional_angle' : constraint on torsional angle
                    (four atom indexes must be specified).

'bennett_proj'    : constraint on the projection onto a given direction
                    of the vector defined by the position of one atom
                    minus the center of mass of the others.
                    G.Roma,J.P.Crocombette: J.Nucl.Mater.403,32(2010)
                  </pre></blockquote>
</ul><ul>
<li> <em>Variables: </em><big><b>constr(1), constr(2), constr(3), constr(4)</b></big>
</li>
<br><li> <em>Type: </em>
</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
                      These variables have different meanings
                      for different constraint types:

                     'type_coord' : constr(1) is the first index of the
                                    atomic type involved
                                    constr(2) is the second index of the
                                    atomic type involved
                                    constr(3) is the cut-off radius for
                                    estimating the coordination
                                    constr(4) is a smoothing parameter

                     'atom_coord' : constr(1) is the atom index of the
                                    atom with constrained coordination
                                    constr(2) is the index of the atomic
                                    type involved in the coordination
                                    constr(3) is the cut-off radius for
                                    estimating the coordination
                                    constr(4) is a smoothing parameter

                       'distance' : atoms indices object of the
                                    constraint, as they appear in
                                    the 'ATOMIC_POSITION' CARD

'planar_angle', 'torsional_angle' : atoms indices object of the
                                    constraint, as they appear in the
                                    'ATOMIC_POSITION' CARD (beware the
                                    order)

                   'bennett_proj' : constr(1) is the index of the atom
                                    whose position is constrained.
                                    constr(2:4) are the three coordinates
                                    of the vector that specifies the
                                    constraint direction.
                  </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>constr_target</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Target for the constrain ( angles are specified in degrees ).
This variable is optional.
                     </pre></blockquote>
</ul>   
    
}


# ------------------------------------------------------------------------
help occupations_table -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variable: </em><big><b>f_inp1</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Occupations of individual states (MAX 10 PER ROW).
For spin-polarized calculations, these are majority spin states.
                  </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>f_inp2</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Occupations of minority spin states (MAX 10 PER ROW)
To be specified only for spin-polarized calculations.
                     </pre></blockquote>
</ul>   
    
}


# ------------------------------------------------------------------------
help atomic_forces -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variable: </em><big><b>X</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> label of the atom as specified in ATOMIC_SPECIES
                  </pre></blockquote>
</ul><ul>
<li> <em>Variables: </em><big><b>fx, fy, fz</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
external force on atom X (cartesian components, Ry/a.u. units)
                  </pre></blockquote>
</ul>   
    
}

