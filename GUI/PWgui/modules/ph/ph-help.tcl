
#
# Help-file automatically created by helpdoc utility
#
#    !!! DO NOT EDIT: CHANGES WILL BE LOST !!!
#
	

# ------------------------------------------------------------------------
help title_ph -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>title_ph</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Title of the job, i.e., a line that is reprinted on output.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help amass -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>amass(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Atomic mass [amu] of each atomic type.
If not specified, masses are read from data file
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help outdir -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>outdir</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> './'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Scratch directory.
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
Prepended to input/output filenames;  must be the same
used in the calculation of unperturbed system.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help niter_ph -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>niter_ph</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 50
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Maximum number of iterations in a scf step.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tr2_ph -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tr2_ph</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1e-10
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Threshold for selfconsistency.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nmix_ph -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nmix_ph</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 4
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Number of iterations used in potential mixing.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help iverbosity -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>iverbosity</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
0 = short output
1 = verbose output
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help reduce_io -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>reduce_io</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Reduce I/O to the strict minimum.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help max_seconds -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>max_seconds</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.d7
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Maximum allowed run time before the job stops smoothly.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fildyn -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fildyn</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'matdyn'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> File where the dynamical matrix is written.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fildrho -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fildrho</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> File where the charge density responses are written.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fildvscf -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fildvscf</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
File where the the potential variation is written
(for later use in electron-phonon calculation).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help epsil -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>epsil</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. in a q=0 calculation for a non metal the
macroscopic dielectric constant of the system is
computed. Do not set epsil to .true. if you have a
metallic system or q/=0: the code will complain and stop
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lrpa -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lrpa</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. the dielectric constant is calculated at the
RPA level with DV_xc=0.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lnoloc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lnoloc</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. the dielectric constant is calculated without
local fields, i.e. by setting DV_H=0 and DV_xc=0.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help trans -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>trans</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .true.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true. the phonons are computed
if trans .and. epsil effective charges are calculated
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lraman -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lraman</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true. calculate nonresonant Raman coefficients
using second-order response as in:
M. Lazzeri and F. Mauri, Phys. Rev. Lett. 90, 036401 (2003)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help eth_rps -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>eth_rps</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.0d-9
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> threshold for calculation of  Pc R |psi&gt;
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help eth_ns -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>eth_ns</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.0e-12
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> threshold for non-scf wavefunction calculation
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help dek -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>dek</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.0e-3
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> delta_xk used for wavefunction derivation wrt k
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help recover -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>recover</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> if .true. restart from an interrupted run
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help elph -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>elph</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true. electron-phonon lambda coeffs are computed

For metals only, requires gaussian smearing.

If elph .and. trans, the lambdas are calculated in the same
run, using the same k-point grid for phonons and lambdas
If elph.and..not.trans, the lambdas are calculated using
previously saved DeltaVscf in fildvscf, previously saved
dynamical matrix, and the present punch file. This allows
the use of a different (larger) k-point grid.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help zue -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>zue</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true. in a q=0 calculation for a non metal the
effective charges are computed from the phonon
density responses. Note that if trans.and.epsil
effective charges are calculated using a different
algorithm. The results should be the same within
numerical noise.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help elop -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>elop</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if true calculate electro-optic tensor
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fpol -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fpol</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true. calculate dynamic polarizabilities
( experimantal stage, see example33 for calculation
  of methane )
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lnscf -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lnscf</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> If .TRUE. the run makes first a nscf calculation.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ldisp -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ldisp</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE. the run calculates phonons for a grid of
q-points specified by nq1, nq2, nq3  - for direct
calculation of the entire phonon dispersion.
The pw.x data file should not be produced using
"calculation='phonon'" in this case.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {nq1 nq2 nq3} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>nq1, nq2, nq3</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Parameters of the Monkhorst-Pack grid (no offset) used
when ldisp=.true. Same meaning as for nk1, nk2, nk3
in the input of pw.x.
         </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
grouphelp {iq1 iq2 iq3} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>iq1, iq2, iq3</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
These go together with nq1, nq2, nq3 and allow to choose
just one point out of the Monkhorst-Pack grid with ldisp=.true.
Note the the actual point chosen is something like
(iq1-1)/nq1, (iq2-1)/nq2, (iq3-1)/nq3 (so, check the
                                       output for what you get). Also make sure that PW left *.wfc
files behind (no 'phonon' is needed though).
         </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
help nrapp -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nrapp</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0, i.e. use all irreps
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Choose the subset of irreducible representations (irreps)
for which the linear response calculation is performed:
"nrapp" irreps, specified in input (see below) are used.

IMPORTANT:
   * nrapp must be &lt;= 3*nat
   * do not specify "nat_todo" together with "nrapp"
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help maxirr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>maxirr</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0, i.e.  use all irreps
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Perform calculations only up to the first "maxirr" irreps.

IMPORTANT:
   * maxirr must be &lt;= 3*nat
   * do not specify "nat_todo" or "nrapp" together with "maxirr"
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nat_todo -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nat_todo</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0, i.e. displace all atoms
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Choose the subset of atoms to be used in the linear response
calculation: "nat_todo" atoms, specified in input (see below)
are displaced.

IMPORTANT:
    * nat_todo &lt;= nat
    * do not specify "nrapp" together with "nat_todo"
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {xq1 xq2 xq3} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b> xq(1)  xq(2)  xq(3)
            </b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The phonon wavevector; must be equal to the one used
in the non-selfconsistent calculation (not read if
ldisp is true).
            </pre></blockquote>
</ul>  
    
}


# ------------------------------------------------------------------------
help irrep_list -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b> irrep(1) irrep(2) ... irrep(nrapp)
               </b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The list of indices of irreps used in the  calculation
if  "nrapp" is specified.
               </pre></blockquote>
</ul>      
    
}


# ------------------------------------------------------------------------
help nat_todo_list -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b> atom(1)  atom(2) ... atom(nat_todo)
               </b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Contains the list of indices of atoms used in the
calculation if "nat_todo" is specified.
               </pre></blockquote>
</ul>      
    
}

