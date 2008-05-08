
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
'scf', 'nscf', 'bands', 'phonon', 'relax', 'md',
'vc-relax', 'vc-md', 'neb', 'smd', 'metadyn'
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
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
'high' | 'default' | 'low' | 'minimal'
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
'from_scratch'  : from scratch
                  NEB and SMD only: the starting path is obtained
                  with a linear interpolation between the images
                  specified in the ATOMIC_POSITIONS card.
                  Note that in the linear interpolation
                  periodic boundary conditions ARE NON USED.

'restart'       : from previous interrupted run
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

.TRUE.  collect wavefunctions from all processors and store
        them into the output data directory outdir/prefix.save

.FALSE. do not collect wavefunctions, leave them in temporary
        local files (one per processor). The resulting format
        will be readable only by jobs running on the same
        number of processors and pools. Useful if you do not
        need the wavefunction or if you want to reduce the I/O
        or the disk occupancy.
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
0  if calculation = 'neb', 'smd';
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
band energies are written every iprint iterations
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tstress -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tstress</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
calculate stress. Set to .TRUE. if calculation='vc-md'
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
print forces. Set to .TRUE. if calculation='relax','md','vc-md'
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
(1 a.u.=4.8378 * 10^-17 s : beware, CP and FPMD codes use
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
see also 'wfcdir'
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help wfcdir -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>wfcdir</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> same as outdir
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
this directory specifies where to store files generated by
each processor (*.wfc{N}, *.igk{N}, etc.). The idea here is
to be able to separately store the largest files, while
the files necessary for restarting still go into 'outdir'
(for now only works for stand alone PW )
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
If .false. it does not open a subdirectory for each k_point
in the prefix.save directory.
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
jobs stops after max_seconds CPU time
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
when the total energy changes less than etot_conv_thr
between two consecutive scf steps.
See also forc_conv_thr - both criteria must be satisfied
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
convergence threshold on forces (a.u) for ionic
minimization: the convergence criterion is satisfied
when all components of all forces are smaller than
forc_conv_thr.
See also etot_conv_thr - both criteria must be satisfied
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help disk_io -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>disk_io</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'default'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Specifies the amount of disk I/O activity
'high':    save all data at each SCF step

'default': save wavefunctions at each SCF step unless
           there is a single k-point per process

'low' :    store wfc in memory, save only at the end

'none':    do not save wfc, not even at the end

If restarting from an interrupted calculation, the code
will try to figure out what is available on disk. The
more you write, the more complete the restart will be.
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
If .TRUE. a sawlike potential simulating an electric field
is added to the bare ionic potential. See variables
edir, eamp, emaxprog, eopreg for the form and size of
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
of L. Bengtsson, PRB 59, 12301 (1999). See variables edir,
emaxprog, eopreg for the form of the correction, that must
be used only in a slab geometry, for surface calculations,
with the discontinuity in the empty space.
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
See the header of PW/bp_c_phase.f90 for documentation
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
(lelfield==.true.), gdir is the direction of the field
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
(lelfield==.true.)
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
Bravais-lattice index:

  ibrav        structure                   celldm(2)-celldm(6)

    0          "free", see above                 not used
    1          cubic P (sc)                      not used
    2          cubic F (fcc)                     not used
    3          cubic I (bcc)                     not used
    4          Hexagonal and Trigonal P        celldm(3)=c/a
    5          Trigonal R                      celldm(4)=cos(alpha)
    6          Tetragonal P (st)               celldm(3)=c/a
    7          Tetragonal I (bct)              celldm(3)=c/a
    8          Orthorhombic P                  celldm(2)=b/a,celldm(3)=c/a
    9          Orthorhombic base-centered(bco) celldm(2)=b/a,celldm(3)=c/a
   10          Orthorhombic face-centered      celldm(2)=b/a,celldm(3)=c/a
   11          Orthorhombic body-centered      celldm(2)=b/a,celldm(3)=c/a
   12          Monoclinic P                    celldm(2)=b/a,celldm(3)=c/a,
                                               celldm(4)=cos(ab)
   13          Monoclinic base-centered        celldm(2)=b/a,celldm(3)=c/a,
                                               celldm(4)=cos(ab)
   14          Triclinic                       celldm(2)= b/a,
                                               celldm(3)= c/a,
                                               celldm(4)= cos(bc),
                                               celldm(5)= cos(ac),
                                               celldm(6)= cos(ab)

For P lattices: the special axis (c) is the z-axis, one basal-plane
vector (a) is along x, the other basal-plane vector (b) is at angle
gamma for monoclinic, at 120 degrees for trigonal and hexagonal
lattices, at 90 degrees for cubic, tetragonal, orthorhombic lattices

sc simple cubic
====================
v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,1)

fcc face centered cubic
====================
v1 = (a/2)(-1,0,1),  v2 = (a/2)(0,1,1), v3 = (a/2)(-1,1,0).

bcc body entered cubic
====================
v1 = (a/2)(1,1,1),  v2 = (a/2)(-1,1,1),  v3 = (a/2)(-1,-1,1).

simple hexagonal and trigonal(p)
====================
v1 = a(1,0,0),  v2 = a(-1/2,sqrt(3)/2,0),  v3 = a(0,0,c/a).

trigonal(r)
===================
for these groups, the z-axis is chosen as the 3-fold axis, but the
crystallographic vectors form a three-fold star around the z-axis,
and the primitive cell is a simple rhombohedron. The crystallographic
vectors are:
      v1 = a(tx,-ty,tz),   v2 = a(0,2ty,tz),   v3 = a(-tx,-ty,tz).
where c=cos(alpha) is the cosine of the angle alpha between any pair
of crystallographic vectors, tc, ty, tz are defined as
     tx=sqrt((1-c)/2), ty=sqrt((1-c)/6), tz=sqrt((1+2c)/3)

simple tetragonal (p)
====================
   v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,c/a)

body centered tetragonal (i)
================================
   v1 = (a/2)(1,-1,c/a),  v2 = (a/2)(1,1,c/a),  v3 = (a/2)(-1,-1,c/a).

simple orthorhombic (p)
=============================
   v1 = (a,0,0),  v2 = (0,b,0), v3 = (0,0,c)

bco base centered orthorhombic
=============================
   v1 = (a/2,b/2,0),  v2 = (-a/2,b/2,0),  v3 = (0,0,c)

face centered orthorhombic
=============================
   v1 = (a/2,0,c/2),  v2 = (a/2,b/2,0),  v3 = (0,b/2,c/2)

body centered orthorhombic
=============================
   v1 = (a/2,b/2,c/2),  v2 = (-a/2,b/2,c/2),  v3 = (-a/2,-b/2,c/2)

monoclinic (p)
=============================
   v1 = (a,0,0), v2= (b*cos(gamma), b*sin(gamma), 0),  v3 = (0, 0, c)
where gamma is the angle between axis a and b

base centered monoclinic
=============================
   v1 = (  a/2,         0,                -c/2),
   v2 = (b*cos(gamma), b*sin(gamma), 0),
   v3 = (  a/2,         0,                  c/2),
where gamma is the angle between axis a and b

triclinic
=============================
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
Crystallographic constants - see description of ibrav variable.

* alat = celldm(1) is the lattice parameter "a" (in BOHR)
* only needed celldm (depending on ibrav) must be specified
* if ibrav=0 only alat = celldm(1) is used (if present)
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
Traditional crystallographic constants (a,b,c in ANGSTROM),
cosab = cosine of the angle between axis a and b
specify either these OR celldm but NOT both.
If ibrav=0 only alat = a is used (if present)
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
(nbnd=nelec/2, see below for nelec);
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
help nelec -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nelec</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> the same as ionic charge (neutral cell)
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
number of electron in the unit cell
(may be noninteger if you wish)

A compensating jellium background is inserted
to remove divergences if the cell is not neutral
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tot_charge -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tot_charge</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
total system charge. Used only if nelec is unspecified,
otherwise it is ignored.
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
May be larger ( for ultrasoft PP ) or somewhat smaller
( but not much smaller ) than the default value. Note that
if you have norm-conserving PP only, setting it to a larger
value than the default is a waste of time.
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
if (.TRUE.) symmetry is not used. Note that a k-point grid
provided in input is used "as is"; an automatically generated
k-point grid will contain only points in the irreducible BZ
of the lattice.  Use with care in low-symmetry large cells
if you cannot afford a k-point grid with the correct symmetry.
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
if (.TRUE.) disable the usage of time reversal (q =&gt; -q)
symmetry in k-point generation
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
                requires a value for degauss

'tetrahedra' :  for metals and DOS calculation
                (see PRB49, 16223 (1994))
                Requires uniform grid of k-points,
                automatically generated (see below)
                Not suitable (because not variational) for
                force/optimization/dynamics calculations

'fixed' :       for insulators with a gap

'from_input' :  The occupation are read from input file.
                Presently works only with one k-point
                (LSDA allowed).
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
help starting_magnetization -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>starting_magnetization(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
starting spin polarization (values between -1 and 1)
on atomic type 'i' in a spin-polarized calculation.
Breaks the symmetry and provides a starting point for
self-consistency. The default value is zero, BUT a value
MUST be specified for AT LEAST one atomic type in spin
polarized calculations. Note that if start from zero
initial magnetization, you will get zero final magnetization
in any case. If you desire to start from an antiferromagnetic
state, you may need to define two different atomic species
corresponding to sublattices of the same atomic type.
If you fix the magnetization with "nelup/neldw" or with
"multiplicity" or with "tot_magnetization", you should
not specify starting_magnetization.
If you are restarting from a previous run, or from an
interrupted run, starting_magnetization is ignored.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {nelup neldw} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>nelup, neldw</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
number of spin-up and spin-down electrons, respectively
Note that this fixes the final value of the magnetization.
The sum must yield nelec that must also be specified
explicitly in this case. Not valid for spin-unpolarized
or noncollinear calculations, only for LSDA. Obsolescent:
use multiplicity or tot_magnetization instead.
         </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
help multiplicity -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>multiplicity</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0 [unspecified]
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
spin multiplicity (2s+1). 1 is singlet, 2 for doublet etc.
Note that this fixes the final value of the magnetization.
if unspecified or a non-zero value is specified in nelup/neldw
then multiplicity variable is ignored.
Do not specify both multiplicity and tot_magnetization.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tot_magnetization -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tot_magnetization</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> -1 [unspecified]
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
majority spin - minority spin (nelup - neldw).
if unspecified or a non-zero value is specified in nelup/neldw
then tot_magnetization variable is ignored.
Do not specify both multiplicity and tot_magnetization.
YES, there is redundancy! nelup/neldw are enough to specify
the spin state. However these variables are not very convenient
and will be eliminated from the input in future versions.
It is recommended to use either 'multiplicity' or equivalently
'tot_magnetization' to specify the spin state.
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
help xc_type -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>xc_type</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Status: </em> presently UNUSED: XC functional is read from PP files
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Exchange-correlation functional
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
<br><li> <em>See: </em> Hubbard_U
         </li>
<br>
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
<br><li> <em>See: </em> Hubbard_U
         </li>
<br>
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
<br><li> <em>Status: </em>
LDA+U works only for a few selected elements. Modify
PW/set_hubbard_l.f90 and PW/tabd.f90 if you plan to use LDA+U with an
element that is not configured there.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
lda_plus_u, Hubbard_alpha(i), Hubbard_U(i): parameters for LDA+U
calculations If lda_plus_u = .TRUE. you must specify, for species i,
the parameters U and (optionally) alpha of the Hubbard model (both in eV).
See: Anisimov, Zaanen, and Andersen, PRB 44, 943 (1991); Anisimov
et al., PRB 48, 16929 (1993); Liechtenstein, Anisimov, and Zaanen, PRB
52, R5467 (1994); Cococcioni and de Gironcoli, PRB 71, 035105 (2005).
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
of projector on localized orbital to be used in the LDA+U
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
                (see PP/poormanwannier.f90 for details)

NB: forces and stress currently implemented only for the
'atomic' choice.
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
Position of the maximum of the sawlike potential along crystal
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
Zone in the unit cell where the sawlike potential decreases.
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
Amplitude of the electric field (in a.u. = 51.44 10^10 V/m )
The sawlike potential increases with slope "eamp" in the
region from (emaxpos+eopreg-1) to (emaxpos), then decreases
to 0 until (emaxpos+eopreg), in units of the crystal
vector "edir". Used only if tefield is .TRUE.
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
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Used to perform constrained calculations in magnetic systems.
Currently available choices:

'none':
         no constraint

'total':
         total magnetization is constrained
         If nspin=4 (noncolin=.True.) constraint is imposed by
         adding a penalty functional to the total energy:

         LAMBDA * SUM_{i} ( magnetization(i) - fixed_magnetization(i) )**2

         where the sum over i runs over the three components of
         the magnetization. Lambda is a real number (see below).
         If nspin=2 constraint is imposed by defining two Fermi
         energies for spin up and down.
         Only fixed_magnetization(3) can be defined in this case.

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

          LAMBDA * ( magnetization(1) - magnetization(3)*tan(theta) )**2

'atomic direction':
          not all the components of the atomic
          magnetic moment are constrained but only the cosine
          of angle1, and the penalty functional is:

          LAMBDA * SUM_{itype} ( mag_mom(3,itype)/mag_mom_tot - cos(angle1(ityp)) )**2
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
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
value of the total magnetization to be maintained fixed when
constrained_magnetization='total'
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help B_field -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>B_field(i), i=1,3</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.d0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
A fixed magnetic field defined by the vector B_field is added
to the exchange and correlation magnetic field.
The three components of the magnetic field are given in Ry.
Only B_field(3) can be used if nspin=2.

In all calculations with a finite magnetic field,
we print the total energy WITHOUT the B dot M term.
In the calculations with the penalty functional we write
only the total energy, NOT the penalty functional.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lambda -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lambda</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
parameter used for constrained_magnetization calculations
NB: LAMBDA is reduced in the first iterations and is increased
    slowly up to the input value.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help report -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>report</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
It's the number of iterations after which the program
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
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .TRUE. the system is assumed to be isolated (a molecule or cluster
in a supercell) and the Makov-Payne correction to the total energy is
computed. An estimate of the vacuum level is also calculated so that
eigenvalues can be properly aligned.
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
number of iterations used in mixing scheme
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
For LDA+U : number of iterations with fixed ns ( ns is the
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
          and is more robust (it fails very seldom)

'david-serial': do not use parallel subspace diagonalization
          in Davidson algorithm (for testing purposes).
          The subspace diagonalization in Davidson is performed
          by a fully distributed-memory parallel algorithm on
          4 or more processors, by default. The allocated memory
          scales down with the number of procs. Procs involved
          in diagonalization can be changed with input parameter
          "ortho_para". On multicore CPUs often it is convenient
          to let only one core per CPU to work on linear algebra.
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
<br><li> <em>Status: </em> OBSOLESCENT: use command-line option " -ndiag XX" instead
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
meaningful for diagonalization='david' and parallel executables.
The number of processors to be used for the parallel subspace
diagonalization algorithm. With the default value (0) the code
tries to use as many processors as available. Note that the
algorithm uses a square number of processors (4, 9, 16, 25,...),
so the actual number of processors used will be the largest
square number less or equal to ortho_para (if set) or to the
total number of processors (if ortho_para is not set).
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
Convergence threshold for the first iterative diagonalization
(the check is on eigenvalue convergence).
For scf calculations, the default is 1.D-2 if starting from a
superposition of atomic orbitals; 1.D-5 if starting from a
charge density. During self consistency the threshold (ethr)
is automatically reduced when approaching convergence.
For non-scf calculations, this is the threshold used in the
iterative diagonalization. The default is conv_thr / nelec.
For 'phonon' calculations, diago_thr_init is ignored:
the threshold is always set to  conv_thr / nelec .
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
A larger value may yield a faster algorithm but uses
more memory
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
For finite electric field calculations (lelfield == .TRUE.),
it defines the intensity of the field in a.u.
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
          ( default for scf, *relax, *md, neb, smd )

'file'  : start from existing "prefix".pot file
          ( default, only possibility for nscf, bands, phonon )
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help startingwfc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>startingwfc</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'atomic'
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

'file':   start from a wavefunction file
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
If .true., use the (VERY EXPERIMENTAL) real-space algorithm
for augmentation charges in ultrasoft pseudopotentials.
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

For constrained dynamics or constrained optimisations add the
CONSTRAINTS card (when the card is present the SHAKE algorithm is
                  automatically used).

For different type of calculation different possibilities are
allowed and different default values apply:

CASE ( calculation = 'relax' )
    'bfgs' :   (default)   a new BFGS quasi-newton algorithm, based
                           on the trust radius procedure, is used
                           for structural relaxation (experimental)
    'damp' :               use damped (quick-min Verlet)
                           dynamics for structural relaxation

CASE ( calculation = 'md' )
    'verlet' : (default)   use Verlet algorithm to integrate
                           Newton's equation
    'langevin'             ion dynamics is over-damped Langevin

CASE ( calculation = 'vc-relax' )
    'damp' :   (default)   use damped (Beeman) dynamics for
                           structural relaxation

CASE ( calculation = 'vc-md' )
    'beeman' : (default)   use Beeman algorithm to integrate
                           Newton's equation
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help phase_space -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>phase_space</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'full'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
'full' :           the full phase-space is used for the ionic
                   dynamics.

'coarse-grained' : a coarse-grained phase-space, defined by a set
                   of constraints, is used for the ionic dynamics
                   (used for calculation of free-energy barriers)
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
                 formula - NOT IMPLEMENTED WITH USPP

'second_order':  as above, with second order formula
                 NOT IMPLEMENTED WITH USPP
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
images. This allowes for the use of smaller supercells.

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
              (first method) see parameters "tempw" and "tolp"
              This is the only method implemented in VC-MD

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
       reduced by -delta_T (.e. delta_t is added to the temperature)

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
<br><li> <em>Default: </em> 10.D0
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
help num_of_images -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>num_of_images</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Number of points used to discretize the path
(it must be larger than 3).
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help opt_scheme -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>opt_scheme</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'quick-min'
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Specify the type of optimization scheme:

'sd'         : steepest descent

'broyden'    : quasi-Newton Broyden's second method (suggested)

'quick-min'  : an optimisation algorithm based on the
               projected velovity Verlet scheme

'langevin'   : finite temperature langevin dynamics of the
               string (smd only). It is used to compute the
               average path and the free-energy profile.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help CI_scheme -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>CI_scheme</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'no-CI'
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Specify the type of Climbing Image scheme:

'no-CI'      : climbing image is not used

'auto'       : original CI scheme. The image highest in energy
               does not feel the effect of springs and is
               allowed to climb along the path

'manual'     : images that have to climb are manually selected.
               See also CLIMBING_IMAGES card
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help first_last_opt -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>first_last_opt</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Also the first and the last configurations are optimised
"on the fly" (these images do not feel the effect of the springs).
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help temp_req -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>temp_req</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.D0 Kelvin
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Temperature used for the langevin dynamics of the string.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ds -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ds</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.D0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Optimisation step length ( Hartree atomic units ).
If opt_scheme="broyden", ds is used as a guess for the diagonal
part of the Jacobian matrix.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {k_max k_min} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>k_max, k_min</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.1D0 Hartree atomic units
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Set them to use a Variable Elastic Constants scheme
elastic constants are in the range [ k_min, k_max ]
this is useful to rise the resolution around the saddle point.
            </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
help path_thr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>path_thr</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.05D0 eV / Angstrom
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The simulation stops when the error ( the norm of the force
orthogonal to the path in eV/A ) is less than path_thr.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help use_masses -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>use_masses</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If. TRUE. the optimisation of the path is performed using
mass-weighted coordinates.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help use_freezing -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>use_freezing</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If. TRUE. the images are optimised according to their error:
only those images with an error larger than half of the largest
are optimised. The other images are kept frozen.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fe_step -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>fe_step(i), i=1,ncolvar</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.04
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Meta-dynamics step length (in principle different for each
collective variable), defined using the same units used
to define the collective variables themselves.
The step also defines the spread of the Gaussian-like bias
potential.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help g_amplitude -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>g_amplitude</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.005 Hartree
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Amplitude of the gaussians used in meta-dynamics.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fe_nstep -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fe_nstep</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 100
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Maximum number of steps used to evaluate the potential of
mean force.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help sw_nstep -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>sw_nstep</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 10
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Number of steps used to switch to the new values of the
collective variables.
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
  'none':    default
  'sd':      steepest descent ( not implemented )
  'damp-pr': damped (Beeman) dynamics of the Parrinello-Rahman
             extended lagrangian
  'damp-w':  damped (Beeman) dynamics of the new Wentzcovitch
             extended lagrangian
  'bfgs':    new BFGS quasi-newton algorithm, based on the trust
             radius procedure, is used for structural relaxation
             (experimental), ion_dynamics must be 'bfgs' too
CASE ( calculation = 'vc-md' )
  'none': default
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

all     = all axis and angles are propagated
volume  = the cell is simply rescaled, without changing the shape
x       = only the x axis is moved
y       = only the y axis is moved
z       = only the z axis is moved
xy      = only the x and y axis are moved, angles are unchanged
xz      = only the x and z axis are moved, angles are unchanged
yz      = only the y and z axis are moved, angles are unchanged
xyz     = x, y and z axis are moved, angles are unchanged
xyt     = x1, x2, y2 (i.e. lower xy triangle of the 2 vectors)
xys     = x1, y1, x2, y2 (i.e. xy square of the 2 vectors)
xyzt    = x1, x2, y2, x3, y3, z3 (i.e. lower xyz triangle of
          the 3 vectors)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help modenum -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>modenum</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
For single-mode phonon calculation : modenum is the index of the
irreducible representation (irrep) into which the reducible
representation formed by the 3*nat atomic displacements are
decomposed in order to perform the phonon calculation.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help xqq -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>xqq(i), i=1,3</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
q-point (units 2pi/a) for phonon calculation.
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
<blockquote><pre> label of the atom
                  </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>Mass_X</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
mass of the atomic species [amu: mass of C = 12]
not used if calculation='scf','nscf', 'bands', 'phonon'
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
alat    : atomic positions are in cartesian coordinates,
          in units of the lattice parameter "a" (default)

bohr    : atomic positions are in cartesian coordinate,
          in atomic units (i.e. Bohr)

angstrom: atomic positions are in cartesian coordinates,
          in Angstrom

crystal : atomic positions are in crystal coordinates, i.e.
          in relative coordinates of the primitive lattice vectors (see below)
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
<blockquote><pre> atomic positions
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
selected components fixed in meta-dynamics, neb, smd, MD dynamics or
structural optimization run.
                           </pre></blockquote>
</ul>   
    
}


# ------------------------------------------------------------------------
help atomic_coordinates_last_image -helpfmt helpdoc -helptext {
       
    
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
           coordinates of the reciprocal lattive vectors

gamma    : use k = 0 (no need to list k-point specifications after card)
           In this case wavefunctions can be chosen as real,
           and specialized subroutines optimized for calculations
           at the gamma point are used (memory and cpu requirements
           are reduced by approximately one half).
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
of the lattice (with all symmetries) and weights (wk)
See the literature for lists of special points and
the corresponding weights.

If the symmetry is lower than the full symmetry
of the lattice, additional points with appropriate
weights are generated.

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
The grid offests;  sk1, sk2, sk3 must be
0 ( no offset ) or 1 ( grid displaced by
half a grid step in the corresponding direction ).
                     </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
help CELL_PARAMETERS_flags -helpfmt helpdoc -helptext {
      <h2>Description of CELL_PARAMETERS card's flags</h2><pre>
Flag "cubic" or "hexagonal" specify if you want to look for symmetries
derived from the cubic symmetry group (default) or from the hexagonal
symmetry group (assuming c axis as the z axis, a axis along the x axis).
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
Crystal lattice vectors:
    v1(1)  v1(2)  v1(3)    ... 1st lattice vector
    v2(1)  v2(2)  v2(3)    ... 2nd lattice vector
    v3(1)  v3(2)  v3(3)    ... 3rd lattice vector

In alat units if celldm(1) was specified or in a.u. otherwise.
                  </pre></blockquote>
</ul>   
    
}


# ------------------------------------------------------------------------
help climbing_images_list -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b> index1, index2, ... indexN
            </b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
index1, index2, ..., indexN are indices of the images to which the
Climbing-Image procedure apply. If more than one image is specified
they must be separated by a comma.
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
                    ( Ch.H. Bennett in Diffusion in Solids, Recent
                      Developments, Ed. by A.S. Nowick and J.J. Burton,
                      New York 1975 ).
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
help ncolvar -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ncolvar</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Number of collective variables.
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help colvar_tol -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>colvar_tol</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Tolerance used for SHAKE.
                  </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help collective_vars_table -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variable: </em><big><b>colvar_type</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>See: </em> constr_type
                  </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
See the definition of constr_type in the CONSTRAINTS card.
                  </pre></blockquote>
</ul><ul>
<li> <em>Variables: </em><big><b>colvar(1), colvar(2), colvar(3), colvar(4)</b></big>
</li>
<br><li> <em>Type: </em>
</li>
<br><li> <em>See: </em> constr(1)
                  </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
These variables have different meanings for
different collective variable types. See the
definition of constr in the CONSTRAINTS card.
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
Occupations of individual states.
For spin-polarized calculation, these are majority spin states.
                  </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>f_inp2</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Occupations of minority spin states for spin-polarized calculation;
specify only for spin-polarized calculation.
                     </pre></blockquote>
</ul>   
    
}

