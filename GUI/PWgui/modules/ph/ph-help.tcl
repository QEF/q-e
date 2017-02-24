
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
<br><li> <em>Default: </em>
value of the <tt>ESPRESSO_TMPDIR</tt> environment variable if set;
<br> current directory ('./') otherwise
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Directory containing input, output, and scratch files;
must be the same as specified in the calculation of
the unperturbed system.
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
<br><li> <em>Default: </em> maxter=100
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Maximum number of iterations in a scf step. If you want
more than 100, edit variable "maxter" in PH/phcom.f90
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
help alpha_mix -helpfmt helpdoc -helptext {
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
help verbosity -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>verbosity</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'default'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre> Options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'debug'</b>, <b>'high'</b>, <b>'medium'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> verbose output
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'low'</b>, <b>'default'</b>, <b>'minimal'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> short output
            </pre></dd>
</dl>
</blockquote>
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
computed. Do not set "epsil" to .true. if you have a
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
If "trans" .and. "epsil" are .true. effective charges are
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
M. Lazzeri and F. Mauri, "PRL 90, 036401 (2003)".
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
help low_directory_check -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>low_directory_check</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. search in the phsave directory only the
                 quantities requested in input.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help only_init -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>only_init</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. only the bands and other initialization quantities are calculated.
(used for GRID parallelization)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help qplot -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>qplot</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> If .true. a list of q points is read from input.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help q2d -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>q2d</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. three q points and relative weights are
read from input. The three q points define the rectangle
q(:,1) + l (q(:,2)-q(:,1)) + m (q(:,3)-q(:,1)) where
0&lt; l,m &lt; 1. The weights are integer and those of points two
and three are the number of points in the two directions.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help q_in_band_form -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>q_in_band_form</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
This flag is used only when qplot is .true. and q2d is
.false.. When .true. each couple of q points q(:,i+1) and
q(:,i) define the line from q(:,i) to q(:,i+1) and nq
points are generated along that line. nq is the weigth of
q(:,i). When .false. only the list of q points given as
input is calculated. The weights are not used.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help electron_phonon -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>electron_phonon</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre>
Options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'simple'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Electron-phonon lambda coefficients are computed
for a given q and a grid of k-points specified by
the variables nk1, nk2, nk3, k1, k2, k3.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'interpolated'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Electron-phonon is calculated by interpolation
over the Brillouin Zone as in M. Wierzbowska, et
al. "arXiv:cond-mat/0504077"
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'lambda_tetra'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
The electron-phonon coefficient \lambda_{q \nu}
is calculated with the optimized tetrahedron method.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'gamma_tetra'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
The phonon linewidth \gamma_{q \nu} is calculated
from the electron-phonon interactions
using the optimized tetrahedron method.
            </pre></dd>
</dl>
<pre>
For metals only, requires gaussian smearing.

If "trans"=.true., the lambdas are calculated in the same
run, using the same k-point grid for phonons and lambdas.
If "trans"=.false., the lambdas are calculated using
previously saved DeltaVscf in "fildvscf", previously saved
dynamical matrix, and the present punch file. This allows
the use of a different (larger) k-point grid.
            </pre>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lshift_q -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lshift_q</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Use a wave-vector grid displaced by half a grid step
in each direction - meaningful only when ldisp is .true.
When this option is set, the q2r.x code cannot be used.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help zeu -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>zeu</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> zeu="epsil"
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. in a q=0 calculation for a non metal the
effective charges are computed from the dielectric
response. This is the default algorithm. If "epsil"=.true.
and "zeu"=.false. only the dielectric tensor is calculated.
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
different from the default one (if "trans" .and. "epsil" )
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
Requires "epsil"=.true. ( experimental stage:
see example09 for calculation of methane ).
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
q-points specified by "nq1", "nq2", "nq3" - for direct
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
help ldiag -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ldiag</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. forces the diagonalization of the dynamical
matrix also when only a part of the dynamical matrix
has been calculated. It is used together with "start_irr"
and "last_irr". If all modes corresponding to a
given irreducible representation have been calculated,
the phonon frequencies of that representation are
correct. The others are zero or wrong. Use with care.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lqdir -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lqdir</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. ph.x creates inside outdir a separate subdirectory
for each q vector. The flag is set to .true. when "ldisp"=.true.
and "fildvscf" /= ' ' or when an electron-phonon
calculation is performed. The induced potential is saved
separately for each q inside the subdirectories.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help search_sym -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>search_sym</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .true.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Set it to .false. if you want to disable the mode
symmetry analysis.
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
when "ldisp"=.true. Same meaning as for nk1, nk2, nk3
in the input of pw.x.
         </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
grouphelp {nk1 nk2 nk3 k1 k2 k3} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>nk1, nk2, nk3, k1, k2, k3</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0,0,0,0,0,0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
When these parameters are specified the phonon program
runs a pw non-self consistent calculation with a different
k-point grid thant that used for the charge density.
This occurs even in the Gamma case.
nk1,nk2,nk3 are the parameters of the Monkhorst-Pack grid
with offset determined by k1,k2,k3.
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
Perform calculations only from "start_irr" to "last_irr"
irreducible representations.

IMPORTANT:
   * "start_irr" must be &lt;= 3*nat
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
Perform calculations only from "start_irr" to "last_irr"
irreducible representations.

IMPORTANT:
   * "start_irr" must be &lt;= 3*nat
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
calculation. Use with care, at your own risk, and be aware
that this is an approximation and may not work.
IMPORTANT:
   * "nat_todo" &lt;= nat
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
Computes only the q points from "start_q" to "last_q".

IMPORTANT:
   * "start_q" must be &lt;= "nqs" (number of q points found)
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
Used only when "ldisp"=.true..
Computes only the q points from "start_q" to "last_q".

IMPORTANT
   * "last_q" must be &lt;= "nqs" (number of q points)
   * do not specify "nat_todo" together with
     "start_q", "last_q"
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {dvscf_star_open dvscf_star_dir dvscf_star_ext dvscf_star_basis dvscf_star_pat} -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>dvscf_star</b></big>
</li>
<br><li> <em>Type: </em>STRUCTURE</li>
<br><li> <em>Default: </em> disabled
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
It contains the following components:

<b>dvscf_star%open</b>  (logical, default: .false.)
<b>dvscf_star%dir</b>   (character, default: outdir//"Rotated_DVSCF" or the
                  ESPRESSO_FILDVSCF_DIR environment variable)
<b>dvscf_star%ext</b>   (character, default: "dvscf") the extension to use
                  for the name of the output files, see below
<b>dvscf_star%basis</b> (character, default: "cartesian") the basis on which
                  the rotated dvscf will be saved
<b>dvscf_star%pat</b>   (logical, default: false) save an optional file with the
                  displacement patterns and q vector for each dvscf file

IF dvscf_star%open is .true. use symmetry to compute and store the variation
of the self-consistent potential on every q* in the star of the present q.

The rotated dvscf will then be stored in directory dvscf_star%dir with name
prefix.dvscf_star%ext.q_name//"1". Where q_name is derived from the coordinates
of the q-point, expressed as fractions in crystalline coordinates
(notice that ph.x reads q-points in cartesian coordinates).
E.g. q_cryst= (0, 0.5, -0.25) -&gt; q_name = "0_1o2_-1o4"

The dvscf can be represented on a basis of cartesian 1-atom displacements
(dvscf_star%basis='cartesian') or on the basis of the modes at the rotated q-point
(dvscf_star%basis='modes'). Notice that the el-ph wannier code requires 'cartesian'.
Each dvscf file comes with a corresponding pattern file with an additional ".pat"
suffix; this file contains information about the basis and the q-point of the dvscf.

Note: rotating dvscf can require a large amount of RAM memory and can be i/o
      intensive; in its current implementation all the operations are done
      on a single processor.
Note2: this feature is currently untested with image parallelisation.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {drho_star_open drho_star_dir drho_star_ext drho_star_basis drho_star_pat} -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>drho_star</b></big>
</li>
<br><li> <em>Type: </em>STRUCTURE</li>
<br><li> <em>Default: </em> disabled
            </li>
<br><li> <em>See: </em> dvscf_star
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
It contains the following components:

<b>drho_star%open</b>  (logical, default: .false.)
<b>drho_star%dir</b>   (character, default: outdir//"Rotated_DRHO" or the
                 ESPRESSO_FILDRHO_DIR environment variable)
<b>drho_star%ext</b>   (character, default: "drho") the extension to use
                 for the name of the output files, see below
<b>drho_star%basis</b> (character, default: "modes") the basis on which
                 the rotated drho will be saved
<b>drho_star%pat</b>   (logical, default: true) save an optional file with the
                 displacement patterns and q vector for each drho file

Like "dvscf_star", but for the perturbation of the charge density.
Notice that the defaults are different.
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
Not used if "ldisp"=.true. or "qplot"=.true.
               </pre></blockquote>
</ul>  
    
}


# ------------------------------------------------------------------------
help nqs -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nqs</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Number of q points in the list. Used only if "qplot"=.true.
                     </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help qPoints -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>xq1, xq2, xq3</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
q-point coordinates; used only with "ldisp"=.true. and qplot=.true.
The phonon wavevector, in units of 2pi/a0 (a0 = lattice parameter).
The meaning of these q points and their weights nq depend on the
flags q2d and q_in_band_form. (NB: nq is integer)
                        </pre></blockquote>
</ul><ul>
<li> <em>Variable: </em><big><b>nq</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The weight of the q-point; the meaning of nq depends
on the flags q2d and q_in_band_form.
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

