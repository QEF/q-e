
#
# Help-file automatically created by helpdoc utility
#
#    !!! DO NOT EDIT: CHANGES WILL BE LOST !!!
#
	

# ------------------------------------------------------------------------
help title_line -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>title_line</b></big>
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
If not specified, masses are read from data file.
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
Prepended to input/output filenames; must be the same
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
<br><li> <em>Default: </em> 100
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
<br><li> <em>Default: </em> 1e-12
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Threshold for self-consistency.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help alpha_mix1 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>alpha_mix(niter)</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> alpha_mix(1)=0.7
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Mixing factor (for each iteration) for updating
the scf potential:

vnew(in) = alpha_mix*vold(out) + (1-alpha_mix)*vold(in)
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
metallic system or q/=0: the code will complain and stop.
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
If .true. the phonons are computed.
If trans .and. epsil are .true. effective charges are
calculated.
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
If .true. calculate non-resonant Raman coefficients
using second-order response as in:
M. Lazzeri and F. Mauri, Phys. Rev. Lett. 90, 036401 (2003).
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
<blockquote><pre> Threshold for calculation of  Pc R |psi&gt;.
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
<blockquote><pre> Threshold for non-scf wavefunction calculation.
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
<blockquote><pre> Delta_xk used for wavefunction derivation wrt k.
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
<blockquote><pre> If .true. restart from an interrupted run.
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
If .true. electron-phonon lambda coefficients are computed.

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
help zeu -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>zeu</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> zeu=epsil
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. in a q=0 calculation for a non metal the
effective charges are computed from the dielectric
response. This is the default algorithm. If epsil=.true.
and zeu=.false. only the dielectric tensor is calculated.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help zue -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>zue</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. in a q=0 calculation for a non metal the
effective charges are computed from the phonon
density responses. This is an alternative algorithm,
different from the default one (if trans .and. epsil )
The results should be the same within numerical noise.
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
If .true. calculate electro-optic tensor.
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
If .true. calculate dynamic polarizabilities
( experimental stage, see example33 for calculation
 of methane ).
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
If .true. the run calculates phonons for a grid of
q-points specified by nq1, nq2, nq3 - for direct
calculation of the entire phonon dispersion.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nogg -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nogg</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. disable the "gamma_gamma" trick used to speed
up calculations at q=0 (phonon wavevector) if the sum over
the Brillouin Zone includes k=0 only. The gamma_gamma
trick exploits symmetry and acoustic sum rule to reduce
the number of linear response calculations to the strict
minimum, as it is done in code phcg.x. This option MUST
BE USED if a run with ph.x is to be followed by a run
with d3.x for third-order terms calculation.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {nq1 nq2 nq3} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>nq1, nq2, nq3</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0,0,0
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
help start_irr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>start_irr</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
            </li>
<br><li> <em>See: </em> last_irr
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Perform calculations only from start_irr to last_irr
irreducible representations.

IMPORTANT:
   * start_irr must be &lt;= 3*nat
   * do not specify "nat_todo" together with
     "start_irr", "last_irr"
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help last_irr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>last_irr</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 3*nat
            </li>
<br><li> <em>See: </em> start_irr
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Perform calculations only from start_irr to last_irr
irreducible representations.

IMPORTANT:
   * start_irr must be &lt;= 3*nat
   * do not specify "nat_todo" together with
     "start_irr", "last_irr"
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
are displaced. Can be used to estimate modes for a molecule
adsorbed over a surface without performing a full fledged
calculation. Use with care, at your own risk,m and be aware
that this is an approximation and may not work.
IMPORTANT:
   * nat_todo &lt;= nat
   * if linear-response is calculated for a given atom, it
     should also be done for all symmetry-equivalent atoms,
     or else you will get incorrect results
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
Note that a single-mode calculation will not give you the
frequency of a single phonon mode: in general, the selected
"modenum" is not an eigenvector. What you get on output is
a column of the dynamical matrix.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help start_q -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>start_q</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
            </li>
<br><li> <em>See: </em> last_q
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Used only when ldisp=.true..
Computes only the q points from start_q to last_q.

IMPORTANT:
   * start_q must be &lt;= nqs (number of q points found)
   * do not specify "nat_todo" together with
     "start_q", "last_q"
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help last_q -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>last_q</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> number of q points
            </li>
<br><li> <em>See: </em> start_q
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Used only when ldisp=.true..
Computes only the q points from start_q to last_q.

IMPORTANT
   * last_q must be &lt;= nqs (number of q points)
   * do not specify "nat_todo" together with
     "start_q", "last_q"
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
The phonon wavevector, in units of 2pi/a0
(a0 = lattice parameter).
Not used if ldisp=.true.
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

