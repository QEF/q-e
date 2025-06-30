
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
<blockquote>
<pre>
A string describing the task to be performed. Options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'scf'</b></tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'nscf'</b></tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'bands'</b></tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'relax'</b></tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'md'</b></tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'vc-relax'</b></tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'vc-md'</b></tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
            </pre></dd>
</dl>
<pre>
(vc = variable-cell).
            </pre>
</blockquote>
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
<blockquote>
<pre>
Currently two verbosity levels are implemented:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'high'</b></tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'low'</b></tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
            </pre></dd>
</dl>
<pre>
<b>'debug'</b> and <b>'medium'</b> have the same effect as <b>'high';</b>
<b>'default'</b> and <b>'minimal'</b> as <b>'low'</b>
            </pre>
</blockquote>
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
<blockquote>
<pre> Available options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'from_scratch'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
From scratch. This is the normal way to perform a PWscf calculation
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'restart'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
From previous interrupted run. Use this switch only if you want to
continue, using the same number of processors and parallelization,
an interrupted calculation. Do not use to start a new one, or to
perform a non-scf calculations.  Works only if the calculation was
cleanly stopped using variable "max_seconds", or by user request
with an "exit file" (i.e.: create a file "prefix".EXIT, in directory
"outdir"; see variables "prefix", "outdir"). The default for
"startingwfc" and "startingpot" is set to 'file'.
            </pre></dd>
</dl>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help wf_collect -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>wf_collect</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> OBSOLETE - NO LONGER IMPLEMENTED
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
1  if "calculation" == 'scf', 'nscf', 'bands';
50 for the other cases
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
number of molecular-dynamics or structural optimization steps
performed in this run. If set to 0, the code performs a quick
"dry run", stopping just after initialization. This is useful
to check for input correctness and to have the summary printed.
NOTE: in MD calculations, the code will perform "nstep" steps
even if restarting from a previously interrupted calculation.
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
When "calculation" == 'md' (molecular dynamics)
trajectory is written every <i>iprint</i> md steps.
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
"calculation" == 'vc-md' or 'vc-relax'
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
"calculation" == 'relax','md','vc-md'
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
This directory specifies where to store files generated by
each processor (*.wfc{N}, *.igk{N}, etc.). Useful for
machines without a parallel file system: set "wfcdir" to
a local file system, while "outdir" should be a parallel
or network file system, visible to all processors. Beware:
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
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
OBSOLETE - NO LONGER IMPLEMENTED
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
Jobs stops after "max_seconds" CPU time. Use this option
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
Convergence threshold on total energy (a.u) for ionic
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
Convergence threshold on forces (a.u) for ionic minimization:
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
<blockquote>
<pre>
Specifies the amount of disk I/O activity:
(only for binary files and xml data file in data directory;
other files printed at each molecular dynamics / structural
optimization step are not controlled by this option )
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'high'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
save charge to disk at each SCF step,
keep wavefunctions on disk (in "distributed" format),
save mixing data as well.
Do not use this option unless you have a good reason!
It is no longer needed to specify 'high' in order to be able
to restart from an interrupted calculation (see "restart_mode")
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'medium'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
save charge to disk at each SCF step,
keep wavefunctions on disk only if more than one k-point,
per process is present, otherwise keep them in memory;
save them to disk only at the end (in "portable" format)
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'low'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
save charge to disk at each SCF step,
keep wavefunctions in memory (for all k-points),
save them to disk only at the end (in "portable" format).
Reduces I/O but increases memory wrt the previous cases
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'nowf'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
save to disk only the xml data file and the charge density
at convergence, never save wavefunctions. Restarting from
an interrupted calculation is not possible with this option.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'minimal'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
save to disk only the xml data file at convergence
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'none'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
do not save anything to disk
            </pre></dd>
</dl>
<pre>
<b>Default</b> is <b>'low'</b> for the scf case, <b>'medium'</b> otherwise.
Note that the needed RAM increases as disk I/O decreases
            </pre>
</blockquote>
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
If .TRUE. and "tefield"==.TRUE. a dipole correction is also
added to the bare ionic potential - implements the recipe
of L. Bengtsson, "PRB 59, 12301 (1999)". See variables "edir",
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
This is different from "tefield" == .true. !
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
In the case of a finite electric field  ( "lelfield" == .TRUE. )
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
If <b>.TRUE.</b> perform orbital magnetization calculation.
If finite electric field is applied ("lelfield"==.true.) only Kubo terms are computed
[for details see New J. Phys. 12, 053032 (2010), "doi:10.1088/1367-2630/12/5/053032"].

The type of calculation is <b>'nscf'</b> and should be performed on an automatically
generated uniform grid of k points.

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
If .TRUE. perform a Berry phase calculation.
See the header of PW/src/bp_c_phase.f90 for documentation.
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
("lelfield"==.true.) "gdir" is the direction of the field.
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
calculated along each symmetry-reduced string.
The same for calculation with finite electric fields
("lelfield"==.true.).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help gate -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>gate</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>See: </em> zgate, relaxz, block, block_1, block_2, block_height
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
In the case of charged cells ("tot_charge" .ne. 0) setting gate = .TRUE.
represents the counter charge (i.e. -tot_charge) not by a homogeneous
background charge but with a charged plate, which is placed at "zgate"
(see below). Details of the gate potential can be found in
T. Brumme, M. Calandra, F. Mauri; "PRB 89, 245406 (2014)".
Note, that in systems which are not symmetric with respect to the plate,
one needs to enable the dipole correction! ("dipfield"=.true.).
Currently, symmetry can be used with gate=.true. but carefully check
that no symmetry is included which maps <i>z</i> to -<i>z</i> even if in principle one
could still use them for symmetric systems (i.e. no dipole correction).
For "nosym"=.false. verbosity is set to 'high'.
Note: this option was called "monopole" in v6.0 and 6.1 of pw.x
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help twochem -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>twochem</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>See: </em> nelec_cond, nbnd_cond, degauss_cond
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
IF .TRUE. , a two chemical potential calculation for the simulation of
photoexcited systems is performed, constraining a fraction of the
electrons in the conduction manifold.
See G. Marini, M. Calandra; "PRB 104, 144103 (2021)".
Note: requires "occupations" to be set to 'smearing'.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lfcp -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lfcp</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE. perform a constant bias potential (constant-mu) calculation
for a system with ESM method. See the header of PW/src/fcp_module.f90
for documentation. To perform the calculation, you must set a namelist FCP.

NB:
- The total energy displayed in output includes the potentiostat
  contribution (-mu*N).
- "calculation" must be 'relax' or 'md'.
- "assume_isolated" = 'esm' and "esm_bc" = 'bc2' or 'bc3' must be set
  in "SYSTEM" namelist.
- ESM-RISM is also supported ("assume_isolated" = 'esm' and "esm_bc" = 'bc1'
  and "trism" = .TRUE.).
- "ignore_wolfe" is always .TRUE., for BFGS.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help trism -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>trism</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE. perform a 3D-RISM-SCF calculation
[for details see H.Sato et al., JCP 112, 9463 (2000), "doi:10.1063/1.481564"].
The solvent's distributions are calculated by 3D-RISM,
though solute is treated as SCF. The charge density and
the atomic positions are optimized, simultaneously with
the solvents. To perform the calculation, you must set
a namelist "RISM" and a card "SOLVENTS".

If "assume_isolated" = 'esm' and "esm_bc" = 'bc1',
Laue-RISM is calculated instead of 3D-RISM
and coupled with ESM method (i.e. ESM-RISM).
[for details see S.Nishihara and M.Otani, "PRB 96, 115429 (2017)"].

The default of "mixing_beta" is 0.2
for both 3D-RISM and Laue-RISM.

For structural relaxation with BFGS,
"ignore_wolfe" is always .TRUE. .
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
  Bravais-lattice index. Optional only if space_group is set.
  If ibrav /= 0, specify EITHER [ "celldm"(1)-"celldm"(6) ]
  OR [ "A", "B", "C", "cosAB", "cosAC", "cosBC" ]
  but NOT both. The lattice parameter "alat" is set to
  alat = celldm(1) (in a.u.) or alat = A (in Angstrom);
  see below for the other parameters.
  For ibrav=0 specify the lattice vectors in "CELL_PARAMETERS",
  optionally the lattice parameter alat = celldm(1) (in a.u.)
  or = A (in Angstrom). If not specified, the lattice parameter is
  taken from "CELL_PARAMETERS"
  IMPORTANT NOTICE 1:
  with ibrav=0 lattice vectors must be given with a sufficiently large
  number of digits and with the correct symmetry, or else symmetry
  detection may fail and strange problems may arise in symmetrization.
  IMPORTANT NOTICE 2:
  do not use celldm(1) or A as a.u. to Ang conversion factor,
  use the true lattice parameters or nothing,
  specify units in "CELL_PARAMETERS" and "ATOMIC_POSITIONS"

ibrav      structure                   celldm(2)-celldm(6)
                                     or: b,c,cosbc,cosac,cosab
  0          free
      crystal axis provided in input: see card "CELL_PARAMETERS"

  1          cubic P (sc)
      v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,1)

  2          cubic F (fcc)
      v1 = (a/2)(-1,0,1),  v2 = (a/2)(0,1,1), v3 = (a/2)(-1,1,0)

  3          cubic I (bcc)
      v1 = (a/2)(1,1,1),  v2 = (a/2)(-1,1,1),  v3 = (a/2)(-1,-1,1)
 -3          cubic I (bcc), more symmetric axis:
      v1 = (a/2)(-1,1,1), v2 = (a/2)(1,-1,1),  v3 = (a/2)(1,1,-1)

  4          Hexagonal and Trigonal P        celldm(3)=c/a
      v1 = a(1,0,0),  v2 = a(-1/2,sqrt(3)/2,0),  v3 = a(0,0,c/a)

  5          Trigonal R, 3fold axis c        celldm(4)=cos(gamma)
      The crystallographic vectors form a three-fold star around
      the z-axis, the primitive cell is a simple rhombohedron:
      v1 = a(tx,-ty,tz),   v2 = a(0,2ty,tz),   v3 = a(-tx,-ty,tz)
      where c=cos(gamma) is the cosine of the angle gamma between
      any pair of crystallographic vectors, tx, ty, tz are:
        tx=sqrt((1-c)/2), ty=sqrt((1-c)/6), tz=sqrt((1+2c)/3)
 -5          Trigonal R, 3fold axis &lt;111&gt;    celldm(4)=cos(gamma)
      The crystallographic vectors form a three-fold star around
      &lt;111&gt;. Defining a' = a/sqrt(3) :
      v1 = a' (u,v,v),   v2 = a' (v,u,v),   v3 = a' (v,v,u)
      where u and v are defined as
        u = tz - 2*sqrt(2)*ty,  v = tz + sqrt(2)*ty
      and tx, ty, tz as for case ibrav=5
      Note: if you prefer x,y,z as axis in the cubic limit,
            set  u = tz + 2*sqrt(2)*ty,  v = tz - sqrt(2)*ty
            See also the note in Modules/latgen.f90

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
      v1 = (a/2,-b/2,0),  v2 = (a/2, b/2,0),  v3 = (0,0,c)
 91          Orthorhombic one-face base-centered A-type
                                             celldm(2)=b/a
                                             celldm(3)=c/a
      v1 = (a, 0, 0),  v2 = (0,b/2,-c/2),  v3 = (0,b/2,c/2)

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
             (unique axis c)                 celldm(3)=c/a,
                                             celldm(4)=cos(gamma)
      v1 = (  a/2,         0,          -c/2),
      v2 = (b*cos(gamma), b*sin(gamma), 0  ),
      v3 = (  a/2,         0,           c/2),
      where gamma=angle between axis a and b projected on xy plane

-13          Monoclinic base-centered        celldm(2)=b/a
             (unique axis b)                 celldm(3)=c/a,
                                             celldm(5)=cos(beta)
      v1 = (  a/2,       b/2,             0),
      v2 = ( -a/2,       b/2,             0),
      v3 = (c*cos(beta),   0,   c*sin(beta)),
      where beta=angle between axis a and c projected on xz plane
 IMPORTANT NOTICE: until QE v.6.4.1, axis for ibrav=-13 had a
 different definition: v1(old) =-v2(now), v2(old) = v1(now)

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
Specify either these OR "A","B","C","cosAB","cosBC","cosAC" NOT both.
Only needed values (depending on "ibrav") must be specified
alat = "celldm"(1) is the lattice parameter "a" (in BOHR)
If "ibrav"==0, only "celldm"(1) is used if present;
cell vectors are read from card "CELL_PARAMETERS"
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {A B C cosAB cosAC cosBC} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>A, B, C, cosAB, cosAC, cosBC</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>See: </em> ibrav
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Traditional crystallographic constants:

  a,b,c in ANGSTROM
  cosAB = cosine of the angle between axis a and b (gamma)
  cosAC = cosine of the angle between axis a and c (beta)
  cosBC = cosine of the angle between axis b and c (alpha)

The axis are chosen according to the value of "ibrav".
Specify either these OR "celldm" but NOT both.
Only needed values (depending on "ibrav") must be specified.

The lattice parameter alat = A (in ANGSTROM ).

If "ibrav" == 0, only A is used if present, and
cell vectors are read from card "CELL_PARAMETERS".
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
number of atoms in the unit cell (ALL atoms, except if
space_group is set, in which case, INEQUIVALENT atoms)
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
for an insulator, "nbnd" = number of valence bands
("nbnd" = # of electrons /2);
<br> for a metal, 20% more (minimum 4 more)
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Number of electronic states (bands) to be calculated.
Note that in spin-polarized calculations the number of
k-point, not the number of bands per k-point, is doubled
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nbnd_cond -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nbnd_cond</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em>
nbnd_cond = "nbnd" - # of electrons / 2 in the collinear case;
                     nbnd_cond = "nbnd" - # of electrons in the noncollinear case.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Number of electronic states in the conduction manifold
for a two chemical-potential calculation ("twochem"=.true.).
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
Total charge of the system. Useful for simulations with charged cells.
By default the unit cell is assumed to be neutral (tot_charge=0).
tot_charge=+1 means one electron missing from the system,
tot_charge=-1 means one additional electron, and so on.

In a periodic calculation a compensating jellium background is
inserted to remove divergences if the cell is not neutral.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help starting_charge -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>starting_charge(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
starting charge on atomic type 'i',
to create starting potential with "startingpot" = 'atomic'.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tot_magnetization -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tot_magnetization</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> -10000 [unspecified]
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Total majority spin charge - minority spin charge.
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
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Starting spin polarization on atomic type 'i' in a spin-polarized
(LSDA or non-collinear/spin-orbit) calculation.
The input values can have an absolute value greater than or equal to 1,
which will be interpreted as the site's magnetic moment.
Alternatively, the values can range between -1 and 1,
which will be interpreted as the site magnetization per valence electron.
For QE-v7.2 and older versions, only the second option is allowed.

If you expect a nonzero magnetization in your ground state,
you MUST either specify a nonzero value for at least one
atomic type, or constrain the magnetization using variable
"tot_magnetization" for LSDA, "constrained_magnetization"
for noncollinear/spin-orbit calculations. If you don't,
you will get a nonmagnetic (zero magnetization) state.
In order to perform LSDA calculations for an antiferromagnetic
state, define two different atomic species corresponding to
sublattices of the same atomic type.

<b>NOTE 1:</b> "starting_magnetization" is ignored in most BUT NOT ALL
cases in non-scf calculations: it is safe to keep the same
values for the scf and subsequent non-scf calculation.

<b>NOTE 2:</b> If you fix the magnetization with
"tot_magnetization", do not specify "starting_magnetization".

<b>NOTE 3:</b> In the noncollinear/spin-orbit case, starting with zero
starting_magnetization on all atoms imposes time reversal
symmetry. The magnetization is never calculated and is
set to zero (the internal variable domag is set to .FALSE.).
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
<br><li> <em>Default: </em> 4 * "ecutwfc"
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Kinetic energy cutoff (Ry) for charge density and potential
For norm-conserving pseudopotential you should stick to the
default value, you can reduce it by a little but it will
introduce noise especially on forces and stress.
If there are ultrasoft PP, a larger value than the default is
often desirable (ecutrho = 8 to 12 times "ecutwfc", typically).
PAW datasets can often be used at 4*"ecutwfc", but it depends
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
Kinetic energy cutoff (Ry) for the exact exchange operator in
EXX type calculations. By default this is the same as "ecutrho"
but in some EXX calculations, a significant speed-up can be obtained
by reducing ecutfock, at the expense of some loss in accuracy.
Must be .gt. "ecutwfc". Not implemented for stress calculation
and for US-PP and PAW pseudopotentials.
Use with care, especially in metals where it may give raise
to instabilities.
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
Three-dimensional FFT mesh (hard grid) for charge
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
Three-dimensional mesh for wavefunction FFT and for the smooth
part of charge density ( smooth grid ).
Coincides with "nr1", "nr2", "nr3" if "ecutrho" = 4 * ecutwfc ( default )
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
if (.TRUE.) symmetry is not used. Consequences:

- if a list of k points is provided in input, it is used
  "as is": symmetry-inequivalent k-points are not generated,
  and the charge density is not symmetrized;

- if a uniform (Monkhorst-Pack) k-point grid is provided in
  input, it is expanded to cover the entire Brillouin Zone,
  irrespective of the crystal symmetry.
  Time reversal symmetry is assumed so k and -k are considered
  as equivalent unless "noinv"=.true. is specified.

Do not use this option unless you know exactly what you want
and what you get. May be useful in the following cases:
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
if (.TRUE.) symmetry is not used, and k points are
forced to have the symmetry of the Bravais lattice;
an automatically generated Monkhorst-Pack grid will contain
all points of the grid over the entire Brillouin Zone,
plus the points rotated by the symmetries of the Bravais
lattice which were not in the original grid. The same
applies if a k-point list is provided in input instead
of a Monkhorst-Pack grid. Time reversal symmetry is assumed
so k and -k are equivalent unless "noinv"=.true. is specified.
This option differs from "nosym" because it forces k-points
in all cases to have the full symmetry of the Bravais lattice
(not all uniform grids have such property!)
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
if (.FALSE.) force real-space FFT grids to be commensurate with
fractionary translations of non-symmorphic symmetry operations,
if present (e.g.: if a fractional translation (0,0,c/4) exists,
the FFT dimension along the c axis must be multiple of 4).
if (.TRUE.) do not impose any constraints to FFT grids, even in
the presence of non-symmorphic symmetry operations.
BEWARE: use_all_frac=.TRUE. may lead to wrong results for
hybrid functionals and phonon calculations. Both cases use
symmetrization in real space that works for non-symmorphic
operations only if the real-space FFT grids are commensurate.
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
<blockquote>
<pre> Available options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'smearing'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
gaussian smearing for metals;
see variables "smearing" and "degauss"
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'tetrahedra'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Tetrahedron method, Bloechl's version:
P.E. Bloechl, "PRB 49, 16223 (1994)"
Requires uniform grid of k-points, to be
automatically generated (see card "K_POINTS").
Well suited for calculation of DOS,
less so (because not variational) for
force/optimization/dynamics calculations.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'tetrahedra_lin'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Original linear tetrahedron method.
To be used only as a reference;
the optimized tetrahedron method is more efficient.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'tetrahedra_opt'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Optimized tetrahedron method:
see M. Kawamura, "PRB 89, 094515 (2014)".
Can be used for phonon calculations as well.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'fixed'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
for insulators with a gap
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'from_input'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
The occupation are read from input file,
card "OCCUPATIONS". Option valid only for a
single k-point, requires "nbnd" to be set
in input. Occupations should be consistent
with the value of "tot_charge".
            </pre></dd>
</dl>
</blockquote>
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
This flag is used for isolated atoms ("nat"=1) together with
"occupations"='from_input'. If it is .TRUE., the wavefunctions
are ordered as the atomic starting wavefunctions, independently
from their eigenvalue. The occupations indicate which atomic
states are filled.

The order of the states is written inside the UPF pseudopotential file.
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
("starting_magnetization"=0.0) the atomic wavefunctions are
radial functions multiplied by spin-angle functions.
For instance:
P -&gt; l=1, j=1/2, m_j=-1/2,1/2. l=1, j=3/2,
     m_j=-3/2, -1/2, 1/2, 3/2.

In the magnetic case with spin-orbit the atomic wavefunctions
can be forced to be spin-angle functions by setting
"starting_spin_angle" to .TRUE..
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
In the spin-orbit case when "domag"=.TRUE., by default,
the starting wavefunctions are initialized as in scalar
relativistic noncollinear case without spin-orbit.

By setting "starting_spin_angle"=.TRUE. this behaviour can
be changed and the initial wavefunctions are radial
functions multiplied by spin-angle functions.

When "domag"=.FALSE. the initial wavefunctions are always
radial functions multiplied by spin-angle functions
independently from this flag.

When "lspinorb" is .FALSE. this flag is not used.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help degauss_cond -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>degauss_cond</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.D0 Ry
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
value of the gaussian spreading (Ry) for brillouin-zone
integration in the conduction manifold in a two-chemical
potential calculation ("twochem"=.true.).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nelec_cond -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nelec_cond</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Number of electrons placed in the conduction manifold in a two-chemical
potential calculation ("twochem"=.true.). Of the total # of
electrons nelec, nelec-nelec_cond will occupy the valence
manifold and nelec_cond will be constrained in the conduction manifold.
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
<blockquote>
<pre>
Available options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'gaussian'</b>, <b>'gauss'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
ordinary Gaussian spreading (Default)
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'methfessel-paxton'</b>, <b>'m-p'</b>, <b>'mp'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Methfessel-Paxton first-order spreading
(see "PRB 40, 3616 (1989)").
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'marzari-vanderbilt'</b>, <b>'cold'</b>, <b>'m-v'</b>, <b>'mv'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Marzari-Vanderbilt-DeVita-Payne cold smearing
(see "PRL 82, 3296 (1999)")
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'fermi-dirac'</b>, <b>'f-d'</b>, <b>'fd'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
smearing with Fermi-Dirac function
            </pre></dd>
</dl>
</blockquote>
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
             DO NOT specify "nspin" in this case;
             specify "noncolin"=.TRUE. instead
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help sic_gamma -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>sic_gamma</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Strength of the gammaDFT potential.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help pol_type -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>pol_type</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre>
Type of polaron in gammaDFT.
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'e'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> electron polaron
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'h'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> hole polaron
            </pre></dd>
</dl>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help sic_energy -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>sic_energy</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Enable the calculation of the total energy in gammaDFT. When .true.,
a preliminary calculation is performed to calculate the electron density
in the absence of the polaron. When .false., the total energy printed in
output should not be considered. For structural relaxations, it is
recommended to use .false. to avoid doubling the computational cost.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help sci_vb -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>sci_vb</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Valence band shift (in eV) through self-consistent
scissor operator. When performing gammaDFT calculations
of polarons, the polaron level is not shifted.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help sci_cb -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>sci_cb</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Conduction band band shift (in eV) through self-consistent
scissor operator. When performing gammaDFT calculations
of polarons, the polaron level is not shifted.
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
See: M. Bernasconi et al, J. Phys. Chem. Solids 56, 501 (1995),
"doi:10.1016/0022-3697(94)00228-2"
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
See Modules/funct.f90 for allowed values.
Overrides the value read from pseudopotential files.
Use with care and if you know what you are doing!
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ace -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ace</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> true
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Use Adaptively Compressed Exchange operator as in
Lin Lin, J. Chem. Theory Comput. 2016, 12, 2242--2249, "doi:10.1021/acs.jctc.6b00092"

Set to false to use standard Exchange (much slower)
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
"input_dft"='PBE0', the default value is 0.25, while for "input_dft"='B3LYP'
the "exx_fraction" default value is 0.20.
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
For more information, see:
J. Chem. Phys. 118, 8207 (2003), "doi:10.1063/1.1564060"
J. Chem. Phys. 124, 219906 (2006), "doi:10.1063/1.2204597"
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help exxdiv_treatment -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>exxdiv_treatment</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'gygi-baldereschi'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre>
Specific for EXX. It selects the kind of approach to be used
for treating the Coulomb potential divergencies at small q vectors.
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'gygi-baldereschi'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
appropriate for cubic and quasi-cubic supercells
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'vcut_spherical'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
appropriate for cubic and quasi-cubic supercells
(untested for non-orthogonal crystal axis)
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'vcut_ws'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
appropriate for strongly anisotropic supercells, see also "ecutvcut"
(untested for non-orthogonal crystal axis)
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'none'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
sets Coulomb potential at G,q=0 to 0.0 (required for GAU-PBE)
            </pre></dd>
</dl>
</blockquote>
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
Specific for EXX. If .true., extrapolate the G=0 term of the
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
Reciprocal space cutoff for correcting Coulomb potential
divergencies at small q vectors.
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
Three-dimensional mesh for q (k1-k2) sampling of
the Fock operator (EXX). Can be smaller than
the number of k-points.

Currently this defaults to the size of the k-point mesh used.
In QE =&lt; 5.0.2 it defaulted to nqx1=nqx2=nqx3=1.
         </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
help localization_thr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>localization_thr</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Overlap threshold over which the exchange integral over a pair of localized orbitals
is included in the evaluation of EXX operator. Any value greater than 0.0 triggers
the SCDM localization and the evaluation on EXX using the localized orbitals.
Very small value of the threshold should yield the same result as the default EXX
evaluation
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help Hubbard_occ -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>Hubbard_occ(ityp,i), (ityp,i) = (1,1) ... (ntyp,3)</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> read from pseudopotentials
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Hubbard occupations is the number of electrons in the
Hubbard manifold. By default they are initialized by
reading the occupations from pseudopotentials. If specified
from the input, then the values read from the pseudopotentials
will be overwritten.
The second index of the Hubbard_occ array corresponds to the
Hubbard manifold number. It is possible to specify up to
three Hubbard manifolds per Hubbard atom. However, if you want
to specify three manifolds then the second and the third manifolds
will be considered as one effective manifold (see Doc/Hubbard_input.pdf)
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
used to compute U (and V) with the linear-response method of
Cococcioni and de Gironcoli, "PRB 71, 035105 (2005)"
(only for DFT+U or DFT+U+V).

Note: Hubbard U and V can be computed using the HP code
which is based on density-functional perturbation theory,
and it gives exactly the same result as the method of
Cococcioni and de Gironcoli.
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
Cococcioni and de Gironcoli, "PRB 71, 035105 (2005)"
(only for DFT+U or DFT+U+V). See also
"PRB 84, 115108 (2011)".
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help dmft -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>dmft</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Status: </em>
Requires compilation with hdf5 support
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If true, nscf calculation will exit in restart mode, scf calculation
will restart from there if DMFT updates are provided as hdf5 archive.
Scf calculation should be used only with "electron_maxstep" = 1.
"K_POINTS" have to be identical and given explicitly with "nosym".
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help dmft_prefix -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>dmft_prefix</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> "prefix"
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
prepended to hdf5 archive: dmft_prefix.h5

DMFT update should be provided in group/dataset as:
- dft_misc_input/band_window with dimension [1, number of k-points, 2 (real + complex)]
- dft_update/delta_N with dimension [number of k-points, number of correlated orbitals,
number of correlated orbitals, 2 (real + complex)]
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ensemble_energies -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ensemble_energies</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If "ensemble_energies" = .true., an ensemble of xc energies
is calculated non-selfconsistently for perturbed
exchange-enhancement factors and LDA vs. PBE correlation
ratios after each converged electronic ground state
calculation.

Ensemble energies can be analyzed with the 'bee' utility
included with libbeef.

Requires linking against libbeef.
"input_dft" must be set to a BEEF-type functional
(e.g. input_dft = 'BEEF-vdW')
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
"edir" = 1, 2 or 3. Used only if "tefield" is .TRUE.
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
Used only if "tefield" is .TRUE.
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
( see below, 0 &lt; eopreg &lt; 1 ). Used only if "tefield" is .TRUE.
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
1 a.u. = 51.4220632*10^10 V/m. Used only if "tefield"==.TRUE.
The saw-like potential increases with slope "eamp" in the
region from ("emaxpos"+"eopreg"-1) to ("emaxpos"), then decreases
to 0 until ("emaxpos"+"eopreg"), in units of the crystal
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
help lforcet -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lforcet</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
When starting a non collinear calculation using an existing density
file from a collinear lsda calculation assumes previous density points in
<i>z</i> direction and rotates it in the direction described by "angle1" and
"angle2" variables for atomic type 1
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
<blockquote>
<pre>
Used to perform constrained calculations in magnetic systems.
Currently available choices:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'none'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
no constraint
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'total'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
total magnetization is constrained by
adding a penalty functional to the total energy:

LAMBDA * SUM_{i} ( magnetization(i) - fixed_magnetization(i) )**2

where the sum over i runs over the three components of
the magnetization. Lambda is a real number (see below).
Noncolinear case only. Use "tot_magnetization" for LSDA
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'atomic'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
atomic magnetization are constrained to the defined
starting magnetization adding a penalty:

LAMBDA * SUM_{i,itype} ( magnetic_moment(i,itype) - mcons(i,itype) )**2

where i runs over the cartesian components (or just z
in the collinear case) and itype over the types (1-ntype).
mcons(:,:) array is defined from starting_magnetization,
(also from angle1, angle2 in the noncollinear case).
lambda is a real number
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'total direction'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
the angle theta of the total magnetization
with the z axis (theta = fixed_magnetization(3))
is constrained:

LAMBDA * ( arccos(magnetization(3)/mag_tot) - theta )**2

where mag_tot is the modulus of the total magnetization.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'atomic direction'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
not all the components of the atomic
magnetic moment are constrained but only the cosine
of angle1, and the penalty functional is:

LAMBDA * SUM_{itype} ( mag_mom(3,itype)/mag_mom_tot - cos(angle1(ityp)) )**2
            </pre></dd>
</dl>
<pre>
N.B.: symmetrization may prevent to reach the desired orientation
of the magnetization. Try not to start with very highly symmetric
configurations or use the nosym flag (only as a last remedy)
            </pre>
</blockquote>
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
fixed when "constrained_magnetization"=='total'
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
<br><li> <em>Default: </em> -1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
determines when atomic magnetic moments are printed on output:
<b>report = 0</b>  never
<b>report =-1</b>  at the beginning of the scf and at convergence
<b>report = N</b>  as -1, plus every N scf iterations
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
<blockquote>
<pre>
Used to perform calculation assuming the system to be
isolated (a molecule or a cluster in a 3D supercell).

Currently available choices:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'none'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
(default): regular periodic calculation w/o any correction.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'makov-payne'</b>, <b>'m-p'</b>, <b>'mp'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
the Makov-Payne correction to the
total energy is computed. An estimate of the vacuum
level is also calculated so that eigenvalues can be
properly aligned. ONLY FOR CUBIC SYSTEMS ("ibrav"=1,2,3).
Theory: G.Makov, and M.C.Payne,
     "Periodic boundary conditions in ab initio
     calculations" , "PRB 51, 4014 (1995)".
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'martyna-tuckerman'</b>, <b>'m-t'</b>, <b>'mt'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Martyna-Tuckerman correction
to both total energy and scf potential. Adapted from:
G.J. Martyna, and M.E. Tuckerman,
"A reciprocal space based method for treating long
range interactions in ab-initio and force-field-based
calculation in clusters", J. Chem. Phys. 110, 2810 (1999),
"doi:10.1063/1.477923".
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'esm'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Effective Screening Medium Method.
For polarized or charged slab calculation, embeds
the simulation cell within an effective semi-
infinite medium in the perpendicular direction
(along z). Embedding regions can be vacuum or
semi-infinite metal electrodes (use "esm_bc" to
choose boundary conditions). If between two
electrodes, an optional electric field
("esm_efield") may be applied. Method described in
M. Otani and O. Sugino, "First-principles calculations
of charged surfaces and interfaces: A plane-wave
nonrepeated slab approach", "PRB 73, 115407 (2006)".

NB:
   - Two dimensional (xy plane) average charge density
     and electrostatic potentials are printed out to
     'prefix.esm1'.

   - Requires cell with a_3 lattice vector along z,
     normal to the xy plane, with the slab centered
     around z=0.

   - For bc2 with an electric field and bc3 boundary
     conditions, the inversion symmetry along z-direction
     is automatically eliminated.

   - In case of calculation='vc-relax', use
     "cell_dofree"='2Dxy' or other parameters so that
     c-vector along z-axis should not be moved.

See "esm_bc", "esm_efield", "esm_w", "esm_nfit".
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'2D'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Truncation of the Coulomb interaction in the z direction
for structures periodic in the x-y plane. Total energy,
forces and stresses are computed in a two-dimensional framework.
Linear-response calculations () done on top of a self-consistent
calculation with this flag will automatically be performed in
the 2D framework as well. Please refer to:
Sohier, T., Calandra, M., &amp; Mauri, F. (2017), "Density functional
perturbation theory for gated two-dimensional heterostructures:
Theoretical developments and application to flexural phonons in graphene",
"PRB, 96, 075448 (2017)".

NB:
   - The length of the unit-cell along the z direction should
     be larger than twice the thickness of the 2D material
     (including electrons). A reasonable estimate for a
     layer's thickness could be the interlayer distance in the
     corresponding layered bulk material. Otherwise,
     the atomic thickness + 10 bohr should be a safe estimate.
     There is also a lower limit of 20 bohr imposed by the cutoff
     radius used to read pseudopotentials (see read_pseudo.f90 in Modules).

   - As for ESM above, only in-plane stresses make sense and one
     should use "cell_dofree"= '2Dxy' in a <b>vc-relax</b> calculation.
            </pre></dd>
</dl>
</blockquote>
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
<blockquote>
<pre>
If "assume_isolated" = 'esm', determines the boundary
conditions used for either side of the slab.

Currently available choices:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'pbc'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> (default): regular periodic calculation (no ESM).
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'bc1'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> Vacuum-slab-vacuum (open boundary conditions).
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'bc2'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Metal-slab-metal (dual electrode configuration).
See also "esm_efield".
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'bc3'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> Vacuum-slab-metal
            </pre></dd>
</dl>
</blockquote>
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
If "assume_isolated" = 'esm', determines the position offset
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
<br><li> <em>See: </em> assume_isolated
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If "assume_isolated" = 'esm' and "esm_bc" = 'bc2', gives the
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
If "assume_isolated" = 'esm', gives the number of z-grid points
for the polynomial fit along the cell edge.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lgcscf -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lgcscf</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE. perform a constant bias potential (constant-mu) calculation
with Grand-Canonical SCF. (JCP 146, 114104 (2017), R.Sundararaman, et al.)

NB:
- The total energy displayed in output includes the potentiostat
  contribution (-mu*N).
- "assume_isolated" = 'esm' and "esm_bc" = 'bc2' or 'bc3' must be set
  in "SYSTEM" namelist.
- ESM-RISM is also supported ("assume_isolated" = 'esm' and "esm_bc" = 'bc1'
  and "trism" = .TRUE.).
- "mixing_mode" has to be 'TF' or 'local-TF', also its default is 'TF.'
- The default of "mixing_beta" is 0.1 with ESM-RISM, 0.2 without ESM-RISM.
- The default of "diago_thr_init" is 1.D-5.
- "diago_full_acc" is always .TRUE. .
- "diago_rmm_conv" is always .TRUE. .
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help gcscf_mu -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>gcscf_mu</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Status: </em> REQUIRED
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The target Fermi energy (eV) of GC-SCF. One can start
with appropriate total charge of the system by giving "tot_charge" .
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help gcscf_conv_thr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>gcscf_conv_thr</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.D-2
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Convergence threshold of Fermi energy (eV) for GC-SCF.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help gcscf_beta -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>gcscf_beta</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.05D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Mixing factor for GC-SCF.
Larger values are recommended,
if systems with small DOS on Fermi surface as graphite.
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
<br><li> <em>See: </em>
london_s6, london_rcut, london_c6, london_rvdw,
dftd3_version, dftd3_threebody, ts_vdw_econv_thr, ts_vdw_isolated, xdm_a1, xdm_a2
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre>
Type of Van der Waals correction. Allowed values:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'grimme-d2'</b>, <b>'Grimme-D2'</b>, <b>'DFT-D'</b>, <b>'dft-d'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Semiempirical Grimme's DFT-D2. Optional variables:
"london_s6", "london_rcut", "london_c6", "london_rvdw"
S. Grimme, J. Comp. Chem. 27, 1787 (2006), "doi:10.1002/jcc.20495"
V. Barone et al., J. Comp. Chem. 30, 934 (2009), "doi:10.1002/jcc.21112"
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'grimme-d3'</b>, <b>'Grimme-D3'</b>, <b>'DFT-D3'</b>, <b>'dft-d3'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Semiempirical Grimme's DFT-D3. Optional variables:
"dftd3_version", "dftd3_threebody"
S. Grimme et al, J. Chem. Phys 132, 154104 (2010), "doi:10.1063/1.3382344"
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'TS'</b>, <b>'ts'</b>, <b>'ts-vdw'</b>, <b>'ts-vdW'</b>, <b>'tkatchenko-scheffler'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Tkatchenko-Scheffler dispersion corrections with first-principle derived
C6 coefficients.
Optional variables: "ts_vdw_econv_thr", "ts_vdw_isolated"
See A. Tkatchenko and M. Scheffler, "PRL 102, 073005 (2009)".
J. Hermann et al., J. Chem. Phys. 159, 174802 (2023), "doi:10.1063/5.0170972"
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'MBD'</b>, <b>'mbd'</b>, <b>'many-body-dispersion'</b>, <b>'mbd_vdw'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Many-body dipersion (MBD) correction to long-range interactions.
Optional variables: "ts_vdw_isolated"
A. Ambrosetti, A. M. Reilly, R. A. DiStasio, A. Tkatchenko, J. Chem. Phys. 140
18A508 (2014).
J. Hermann et al., J. Chem. Phys. 159, 174802 (2023), "doi:10.1063/5.0170972"
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'XDM'</b>, <b>'xdm'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Exchange-hole dipole-moment model. Optional variables: "xdm_a1", "xdm_a2"
A. D. Becke et al., J. Chem. Phys. 127, 154108 (2007), "doi:10.1063/1.2795701"
A. Otero de la Roza et al., J. Chem. Phys. 136, 174109 (2012),
"doi:10.1063/1.4705760"
            </pre></dd>
</dl>
<pre> Note that non-local functionals (eg vdw-DF) are NOT specified here but in "input_dft"
            </pre>
</blockquote>
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
<br><li> <em>Status: </em>
OBSOLESCENT, same as "vdw_corr"='DFT-D'
         </li>
<br>
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
help london_c6 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>london_c6(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> standard Grimme-D2 values
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
atomic C6 coefficient of each atom type

( if not specified default values from S. Grimme, J. Comp. Chem. 27, 1787 (2006),
  "doi:10.1002/jcc.20495" are used; see file Modules/mm_dispersion.f90 )
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help london_rvdw -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>london_rvdw(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> standard Grimme-D2 values
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
atomic vdw radii of each atom type

( if not specified default values from S. Grimme, J. Comp. Chem. 27, 1787 (2006),
  "doi:10.1002/jcc.20495" are used; see file Modules/mm_dispersion.f90 )
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
help dftd3_version -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>dftd3_version</b></big>
</li>
<br><li> <em>Type: </em>integer</li>
<br><li> <em>Default: </em> 3
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre>
Version of Grimme implementation of Grimme-D3:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>dftd3_version = 2</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Original Grimme-D2 parametrization
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>dftd3_version = 3</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Grimme-D3 (zero damping)
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>dftd3_version = 4</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Grimme-D3 (BJ damping)
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>dftd3_version = 5</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Grimme-D3M (zero damping)
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>dftd3_version = 6</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Grimme-D3M (BJ damping)
            </pre></dd>
</dl>
<pre>
NOTE: not all functionals are parametrized.
            </pre>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help dftd3_threebody -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>dftd3_threebody</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> TRUE
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Turn three-body terms in Grimme-D3 on. If .false. two-body contributions
only are computed, using two-body parameters of Grimme-D3.
If dftd3_version=2, three-body contribution is always disabled.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ts_vdw_econv_thr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ts_vdw_econv_thr</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.D-6
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Optional: controls the convergence of the vdW energy (and forces). The default value
is a safe choice, likely too safe, but you do not gain much in increasing it
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ts_vdw_isolated -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ts_vdw_isolated</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Optional: set it to .TRUE. when computing the Tkatchenko-Scheffler vdW energy or the
Many-Body dispersion (MBD) energy for an isolated (non-periodic) system.
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
<br><li> <em>Status: </em>
OBSOLESCENT, same as "vdw_corr"='xdm'
         </li>
<br>
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
Damping function parameter a1 (adimensional). It is NOT necessary to give
a value if the functional is one of B86bPBE, PW86PBE, PBE, BLYP. For functionals
in this list, the coefficients are given in:
   "https://github.com/aoterodelaroza/postg/blob/master/xdm.param"
   or "https://erin-r-johnson.github.io/software/"
   A. Otero de la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013),
   "doi:10.1063/1.4705760"
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
Damping function parameter a2 (angstrom). It is NOT necessary to give
a value if the functional is one of B86bPBE, PW86PBE, PBE, BLYP. For functionals
in this list, the coefficients are given in:
   "https://github.com/aoterodelaroza/postg/blob/master/xdm.param"
   or "https://erin-r-johnson.github.io/software/"
   A. Otero de la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013),
   "doi:10.1063/1.4705760"
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
are of type crystal_sg. See also "uniqueb",
"origin_choice", "rhombohedral"
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
unique "ibrav" (-12 or -13) are used, and symmetry
equivalent positions are chosen assuming that the
twofold axis or the mirror normal is parallel to the
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
the use of two different origins. "origin_choice"=1,
means the first origin, while "origin_choice"=2 is the
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
to the rhombohedral axes and "ibrav"=5 is used in both cases.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help zgate -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>zgate</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.5
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
used only if "gate" = .TRUE.
Specifies the position of the charged plate which represents
the counter charge in doped systems ("tot_charge" .ne. 0).
In units of the unit cell length in <i>z</i> direction, "zgate" in ]0,1[
Details of the gate potential can be found in
T. Brumme, M. Calandra, F. Mauri; "PRB 89, 245406 (2014)".
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help relaxz -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>relaxz</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
used only if "gate" = .TRUE.
Allows the relaxation of the system towards the charged plate.
Use carefully and utilize either a layer of fixed atoms or a
potential barrier ("block"=.TRUE.) to avoid the atoms moving to
the position of the plate or the dipole of the dipole
correction ("dipfield"=.TRUE.).
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help block -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>block</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
used only if "gate" = .TRUE.
Adds a potential barrier to the total potential seen by the
electrons to mimic a dielectric in field effect configuration
and/or to avoid electrons spilling into the vacuum region for
electron doping. Potential barrier is from "block_1" to "block_2" and
has a height of block_height.
If "dipfield" = .TRUE. then "eopreg" is used for a smooth increase and
decrease of the potential barrier.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help block_1 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>block_1</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.45
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
used only if "gate" = .TRUE. and "block" = .TRUE.
lower beginning of the potential barrier, in units of the
unit cell size along <i>z,</i> "block_1" in ]0,1[
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help block_2 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>block_2</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.55
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
used only if "gate" = .TRUE. and "block" = .TRUE.
upper beginning of the potential barrier, in units of the
unit cell size along <i>z,</i> "block_2" in ]0,1[
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help block_height -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>block_height</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.1
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
used only if "gate" = .TRUE. and "block" = .TRUE.
Height of the potential barrier in Rydberg.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nextffield -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nextffield</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Number of activated external ionic force fields.
See Doc/ExternalForceFields.tex for further explanation and parameterizations
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
maximum number of iterations in a scf step. If exact exchange is active,
this will affect the inner loops.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help exx_maxstep -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>exx_maxstep</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 100
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
maximum number of outer iterations in a scf calculation with exact exchange.
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
iterative diagonalization: see "diago_thr_init"
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
If .TRUE. this turns on the use of an adaptive "conv_thr" for
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
When "adaptive_thr" = .TRUE. this is the convergence threshold
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
When "adaptive_thr" = .TRUE. the convergence threshold for
each scf cycle is given by:
max( "conv_thr", "conv_thr_multi" * dexx )
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
<blockquote>
<pre> Available options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'plain'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> charge density Broyden mixing
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'TF'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
as above, with simple Thomas-Fermi screening
(for highly homogeneous systems)
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'local-TF'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
as above, with local-density-dependent TF screening
(for highly inhomogeneous systems)
            </pre></dd>
</dl>
</blockquote>
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
<blockquote>
<pre> Available options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'david'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Davidson iterative diagonalization with overlap matrix
(default). Fast, may in some rare cases fail.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'cg'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Conjugate-gradient-like band-by-band diagonalization.
MUCH slower than 'david' but uses less memory and is
(a little bit) more robust.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'ppcg'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
PPCG iterative diagonalization (end support Dec 2024)
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'paro'</b>, <b>'ParO'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
ParO iterative diagonalization
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'rmm-davidson'</b>, <b>'rmm-paro'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
RMM-DIIS iterative diagonalization.
To stabilize the SCF loop
RMM-DIIS is alternated with calls to Davidson or
ParO  solvers depending on the string used.
Other variables that can be used to tune the behavior of
RMM-DIIS are:  "diago_rmm_ndim" and "diago_rmm_conv"
            </pre></dd>
</dl>
</blockquote>
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

For non-scf calculations: default is ("conv_thr"/N elec)/10.
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
For conjugate gradient diagonalization:  max number of iterations
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help diago_david_ndim -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>diago_david_ndim</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 2
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
For Davidson diagonalization: dimension of workspace
(number of wavefunction packets, at least 2 needed).
A larger value may yield a smaller number of iterations in
the algorithm but uses more memory and more CPU time in
subspace diagonalization (cdiaghg/rdiaghg). You may try
"diago_david_ndim"=4 if you are not tight on memory
and if the time spent in subspace diagonalization is small
compared to the time spent in h_psi
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help diago_rmm_ndim -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>diago_rmm_ndim</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 4
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
For RMM-DIIS diagonalization: dimension of workspace
(number of wavefunction packets, at least 2 needed).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help diago_rmm_conv -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>diago_rmm_conv</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE., RMM-DIIS is performed up to converge.
If .FALSE., RMM-DIIS is performed only once.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help diago_gs_nblock -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>diago_gs_nblock</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 16
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
For RMM-DIIS diagonalization:
blocking size of Gram-Schmidt orthogonalization
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
1 a.u. = 36.3609*10^10 V/m). Used only if "lelfield"==.TRUE.
and if k-points ("K_POINTS" card) are not automatic.
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
cartesian axis. Used only if "lelfield"==.TRUE. and if
k-points ("K_POINTS" card) are automatic.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help efield_phase -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>efield_phase</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'none'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre> Available options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'read'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
set the zero of the electronic polarization (with "lelfield"==.true..)
to the result of a previous calculation
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'write'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
write on disk data on electronic polarization to be read in another
calculation
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'none'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
none of the above points
            </pre></dd>
</dl>
</blockquote>
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
<blockquote>
<pre> Available options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'atomic'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
starting potential from atomic charge superposition
(default for scf, *relax, *md)
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'file'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
start from existing "charge-density.xml" file in the
directory specified by variables "prefix" and "outdir"
For nscf and bands calculation this is the default
and the only sensible possibility.
            </pre></dd>
</dl>
</blockquote>
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
<blockquote>
<pre> Available options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'atomic'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Start from superposition of atomic orbitals.
If not enough atomic orbitals are available,
fill with random numbers the remaining wfcs
The scf typically starts better with this option,
but in some high-symmetry cases one can "loose"
valence states, ending up in the wrong ground state.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'atomic+random'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
As above, plus a superimposed "randomization"
of atomic orbitals. Prevents the "loss" of states
mentioned above.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'random'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Start from random wfcs. Slower start of scf but safe.
It may also reduce memory usage in conjunction with
"diagonalization"='cg'.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'file'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Start from an existing wavefunction file in the
directory specified by variables "prefix" and "outdir".
            </pre></dd>
</dl>
</blockquote>
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
If .true., use a real-space algorithm for augmentation
charges of ultrasoft pseudopotentials and PAWsets.
Faster but numerically less accurate than the default
G-space algorithm. Use with care and after testing!
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help real_space -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>real_space</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true., exploit real-space localization to compute
matrix elements for nonlocal projectors. Faster and in
principle better scaling than the default G-space algorithm,
but numerically less accurate, may lead to some loss of
translational invariance. Use with care and after testing!
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
<blockquote>
<pre> Available options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'default'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
if restarting, use atomic positions read from the
restart file; in all other cases, use atomic
positions from standard input.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'from_input'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
read atomic positions from standard input, even if restarting.
            </pre></dd>
</dl>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ion_velocities -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ion_velocities</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'default'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre>
Initial ionic velocities. Available options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'default'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
start a new simulation from random thermalized
distribution of velocities if "tempw" is set,
with zero velocities otherwise; restart from
atomic velocities read from the restart file
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'from_input'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
start or continue the simulation with atomic
velocities read from standard input - see card
"ATOMIC_VELOCITIES"
            </pre></dd>
</dl>
</blockquote>
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
<blockquote>
<pre>
Specify the type of ionic dynamics.

For different type of calculation different possibilities are
allowed and different default values apply:

<b>CASE</b> ( "calculation" == 'relax' )
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'bfgs'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
<b>(default)</b>  use BFGS quasi-newton algorithm,
based on the trust radius procedure,
for structural relaxation
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'damp'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
use damped (quick-min Verlet)
dynamics for structural relaxation
Can be used for constrained
optimisation: see "CONSTRAINTS" card
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'fire'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
use the FIRE minimization algorithm employing the
semi-implicit Euler integration scheme
see:
   Bitzek et al.,"PRL, 97, 170201, (2006)", "doi: 10.1103/PhysRevLett.97.170201"
   Guenole et al.,CMS, 175, 109584, (2020), "doi: 10.1016/j.commatsci.2020.109584"

Can be used for constrained
optimisation: see "CONSTRAINTS" card
            </pre></dd>
</dl>
<pre>
<b>CASE</b> ( "calculation" == 'md' )
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'verlet'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
<b>(default)</b>  use Verlet algorithm to integrate
Newton's equation. For constrained
dynamics, see "CONSTRAINTS" card
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'langevin'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
ion dynamics is over-damped Langevin
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'langevin-smc'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
over-damped Langevin with Smart Monte Carlo:
see R.J. Rossky, JCP, 69, 4628 (1978), "doi:10.1063/1.436415"
            </pre></dd>
</dl>
<pre>
<b>CASE</b> ( "calculation" == 'vc-relax' )
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'bfgs'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
<b>(default)</b>  use BFGS quasi-newton algorithm;
"cell_dynamics" must be 'bfgs' too
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'damp'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
use damped (Beeman) dynamics for
structural relaxation
            </pre></dd>
</dl>
<pre>
<b>CASE</b> ( "calculation" == 'vc-md' )
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'beeman'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
<b>(default)</b>  use Beeman algorithm to integrate
Newton's equation
            </pre></dd>
</dl>
</blockquote>
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
<blockquote>
<pre>
Used to extrapolate the potential from preceding ionic steps.
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'none'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> no extrapolation
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'atomic'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
extrapolate the potential as if it was a sum of
atomic-like orbitals
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'first_order'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
extrapolate the potential with first-order
formula
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'second_order'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
as above, with second order formula
            </pre></dd>
</dl>
<pre>
Note: 'first_order' and 'second-order' extrapolation make sense
only for molecular dynamics calculations
            </pre>
</blockquote>
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
<blockquote>
<pre>
Used to extrapolate the wavefunctions from preceding ionic steps.
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'none'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> no extrapolation
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'first_order'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
extrapolate the wave-functions with first-order formula.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'second_order'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
as above, with second order formula.
            </pre></dd>
</dl>
<pre>
Note: <b>'first_order'</b> and <b>'second-order'</b> extrapolation make sense
only for molecular dynamics calculations
            </pre>
</blockquote>
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
<blockquote>
<pre> Available options are:
               </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'rescaling'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
control ionic temperature via velocity rescaling
(first method) see parameters "tempw", "tolp", and
"nraise" (for VC-MD only). This rescaling method
is the only one currently implemented in VC-MD
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'rescale-v'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
control ionic temperature via velocity rescaling
(second method) see parameters "tempw" and "nraise"
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'rescale-T'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
scale temperature of the thermostat every "nraise" steps
by "delta_t", starting from "tempw".
The temperature is controlled via velocitiy rescaling.
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'reduce-T'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
reduce temperature of the thermostat every "nraise" steps
by the (negative) value "delta_t", starting from "tempw".
If  "delta_t" is positive, the target temperature is augmented.
The temperature is controlled via velocitiy rescaling.
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'berendsen'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
control ionic temperature using "soft" velocity
rescaling - see parameters "tempw" and "nraise"
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'andersen'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
control ionic temperature using Andersen thermostat
see parameters "tempw" and "nraise"
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'svr'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
control ionic temperature using stochastic-velocity rescaling
(Donadio, Bussi, Parrinello, J. Chem. Phys. 126, 014101, 2007),
with parameters "tempw" and "nraise".
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'initial'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
initialize ion velocities to temperature "tempw"
and leave uncontrolled further on
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'not_controlled'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
(default) ionic temperature is not controlled
               </pre></dd>
</dl>
</blockquote>
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
if "ion_temperature" == 'rescale-T' :
       at each step the instantaneous temperature is multiplied
       by delta_t; this is done rescaling all the velocities.

if "ion_temperature" == 'reduce-T' :
       every 'nraise' steps the instantaneous temperature is
       reduced by -"delta_t" (i.e. "delta_t" &lt; 0 is added to T)

The instantaneous temperature is calculated at the end of
every ionic move and BEFORE rescaling. This is the temperature
reported in the main output.

For "delta_t" &lt; 0, the actual average rate of heating or cooling
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
if "ion_temperature" == 'reduce-T' :
       every "nraise" steps the instantaneous temperature is
       reduced by -"delta_t" (i.e. "delta_t" is added to the temperature)

if "ion_temperature" == 'rescale-v' :
       every "nraise" steps the average temperature, computed from
       the last "nraise" steps, is rescaled to "tempw"

if "ion_temperature" == 'rescaling' and "calculation" == 'vc-md' :
       every "nraise" steps the instantaneous temperature
       is rescaled to "tempw"

if "ion_temperature" == 'berendsen' :
       the "rise time" parameter is given in units of the time step:
       tau = nraise*dt, so dt/tau = 1/nraise

if "ion_temperature" == 'andersen' :
       the "collision frequency" parameter is given as nu=1/tau
       defined above, so nu*dt = 1/nraise

if "ion_temperature" == 'svr' :
       the "characteristic time" of the thermostat is set to
       tau = nraise*dt
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
Max reduction factor for "conv_thr" during structural optimization
"conv_thr" is automatically reduced when the relaxation
approaches convergence so that forces are still accurate,
but "conv_thr" will not be reduced to less that "conv_thr" / "upscale".
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
When "bfgs_ndim" = 1, the standard quasi-Newton BFGS method is
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
BFGS is reset when "trust_radius" &lt; "trust_radius_min".
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
help fire_alpha_init -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fire_alpha_init</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.2D0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Initial value of the alpha mixing factor in the FIRE minimization scheme;
recommended values are between 0.1 and 0.3
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fire_falpha -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fire_falpha</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.99D0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Scaling of the alpha mixing parameter for steps with P &gt; 0;
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fire_nmin -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fire_nmin</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 5
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Minimum number of steps with P &gt; 0 before increase of "dt"
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fire_f_inc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fire_f_inc</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.1D0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Factor for increasing "dt"
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fire_f_dec -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fire_f_dec</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.5D0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Factor for decreasing "dt"
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fire_dtmax -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fire_dtmax</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 10.D0
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Determines the maximum value of "dt" in the FIRE minimization;
dtmax = fire_dtmax*"dt"
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
<blockquote>
<pre>
Specify the type of dynamics for the cell.
For different type of calculation different possibilities
are allowed and different default values apply:

<b>CASE</b> ( "calculation" == 'vc-relax' )
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'none'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> no dynamics
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'sd'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> steepest descent ( not implemented )
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'damp-pr'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
damped (Beeman) dynamics of the Parrinello-Rahman extended lagrangian
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'damp-w'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
damped (Beeman) dynamics of the new Wentzcovitch extended lagrangian
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'bfgs'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
BFGS quasi-newton algorithm <b>(default)</b>
"ion_dynamics" must be <b>'bfgs'</b> too
            </pre></dd>
</dl>
<pre>
<b>CASE</b> ( "calculation" == 'vc-md' )
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'none'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> no dynamics
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'pr'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
(Beeman) molecular dynamics of the Parrinello-Rahman extended lagrangian
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'w'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
(Beeman) molecular dynamics of the new Wentzcovitch extended lagrangian
            </pre></dd>
</dl>
</blockquote>
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
<br><li> <em>Default: </em> 2.0 for variable-cell calculations, 1.0 otherwise
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
<blockquote>
<pre>
Select which of the cell parameters should be moved:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'all'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> all axis and angles are moved
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'ibrav'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
all axis and angles are moved,
               but the lattice remains consistent
               with the initial ibrav choice. You can use this option in combination
               with any other one by specifying "ibrav+option". Please note that some
               combinations do not make sense for some crystals and will guarantee that
               the relax will never converge. E.g. 'ibrav+2Dxy' is not a problem for
               hexagonal cells, but will never converge for cubic ones.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'a'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> the x component of axis 1 (v1_x) is fixed
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'b'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> the y component of axis 2 (v2_y) is fixed
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'c'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> the z component of axis 3 (v3_z) is fixed
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'fixa'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> axis 1 (v1_x,v1_y,v1_z) is fixed
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'fixb'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> axis 2 (v2_x,v2_y,v2_z) is fixed
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'fixc'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> axis 3 (v3_x,v3_y,v3_z) is fixed
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'x'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> only the x component of axis 1 (v1_x) is moved
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'y'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> only the y component of axis 2 (v2_y) is moved
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'z'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> only the z component of axis 3 (v3_z) is moved
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'xy'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> only v1_x and v2_y are moved
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'xz'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> only v1_x and v3_z are moved
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'yz'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> only v2_y and v3_z are moved
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'xyz'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> only v1_x, v2_y, v3_z are moved
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'shape'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> all axis and angles, keeping the volume fixed
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'volume'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> the volume changes, keeping all angles fixed (i.e. only celldm(1) changes)
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'2Dxy'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> only x and y components are allowed to change
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'2Dshape'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> as above, keeping the area in xy plane fixed
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'epitaxial_ab'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> fix axis 1 and 2 while allowing axis 3 to move
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'epitaxial_ac'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> fix axis 1 and 3 while allowing axis 2 to move
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'epitaxial_bc'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> fix axis 2 and 3 while allowing axis 1 to move
            </pre></dd>
</dl>
<pre>
BEWARE: if axis are not orthogonal, some of these options do not
        work (symmetry is broken). If you are not happy with them,
        edit subroutine init_dofree in file Modules/cell_base.f90
            </pre>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fcp_mu -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fcp_mu</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Status: </em> REQUIRED
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The target Fermi energy (eV). One can start
with appropriate total charge of the system by giving "tot_charge" .
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fcp_dynamics -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fcp_dynamics</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre>
Specify the type of dynamics for the Fictitious Charge Particle (FCP).

For different type of calculation different possibilities
are allowed and different default values apply:

<b>CASE</b> ( "calculation" == 'relax' )
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'bfgs'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
<b>(default)</b> BFGS quasi-newton algorithm, coupling with ions relaxation
"ion_dynamics" must be <b>'bfgs'</b> too
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'newton'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Newton-Raphson algorithm with DIIS
"ion_dynamics" must be <b>'damp'</b> too
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'damp'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
damped (quick-min Verlet) dynamics for FCP relaxation
"ion_dynamics" must be <b>'damp'</b> too
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'lm'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Line-Minimization algorithm for FCP relaxation
"ion_dynamics" must be <b>'damp'</b> too
            </pre></dd>
</dl>
<pre>
<b>CASE</b> ( "calculation" == 'md' )
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'velocity-verlet'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
<b>(default)</b> Velocity-Verlet algorithm to integrate Newton's equation.
"ion_dynamics" must be <b>'verlet'</b> too
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'verlet'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
<b>Verlet</b> algorithm to integrate Newton's equation.
"ion_dynamics" must be <b>'verlet'</b> too
            </pre></dd>
</dl>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fcp_conv_thr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fcp_conv_thr</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.D-2
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Convergence threshold on force (eV) for FCP relaxation.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fcp_ndiis -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fcp_ndiis</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 4
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Size of DIIS for FCP relaxation,
used only if "fcp_dynamics" = 'newton'.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fcp_mass -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fcp_mass</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em>
5.D+6 / (xy area) for ESM only;
5.D+4 / (xy area) for ESM-RISM
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Mass of the FCP.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fcp_velocity -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fcp_velocity</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> determined by "fcp_temperature"
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Initial velocity of the FCP.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fcp_temperature -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fcp_temperature</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> "ion_temperature"
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre> Available options are:
               </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'rescaling'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
control FCP's temperature via velocity rescaling
(first method) see parameters "fpc_tempw" and "fcp_tolp".
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'rescale-v'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
control FCP's temperature via velocity rescaling
(second method) see parameters "fcp_tempw" and "fcp_nraise"
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'rescale-T'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
control FCP's temperature via velocity rescaling
(third method) see parameter "fcp_delta_t"
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'reduce-T'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
reduce FCP's temperature every "fcp_nraise" steps
by the (negative) value "fcp_delta_t"
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'berendsen'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
control FCP's temperature using "soft" velocity
rescaling - see parameters "fcp_tempw" and "fcp_nraise"
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'andersen'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
control FCP's temperature using Andersen thermostat
see parameters "fcp_tempw" and "fcp_nraise"
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'initial'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
initialize FCP's velocities to temperature "fcp_tempw"
and leave uncontrolled further on
               </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'not_controlled'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
<b>(default)</b> FCP's temperature is not controlled
               </pre></dd>
</dl>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fcp_tempw -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fcp_tempw</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> "tempw"
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Starting temperature (Kelvin) in FCP dynamics runs
target temperature for most thermostats.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fcp_tolp -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fcp_tolp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> "tolp"
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
help fcp_delta_t -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fcp_delta_t</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> "delta_t"
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if "fcp_temperature" == 'rescale-T' :
       at each step the instantaneous temperature is multiplied
       by fcp_delta_t; this is done rescaling all the velocities.

if "fcp_temperature" == 'reduce-T' :
       every "fcp_nraise" steps the instantaneous temperature is
       reduced by -"fcp_delta_t" (i.e. "fcp_delta_t" &lt; 0 is added to T)

The instantaneous temperature is calculated at the end of
FCP's move and BEFORE rescaling. This is the temperature
reported in the main output.

For "fcp_delta_t" &lt; 0, the actual average rate of heating or cooling
should be roughly C*fcp_delta_t/(fcp_nraise*dt) (C=1 for an
ideal gas, C=0.5 for a harmonic solid, theorem of energy
equipartition between all quadratic degrees of freedom).
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fcp_nraise -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fcp_nraise</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> "nraise"
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if "fcp_temperature" == 'reduce-T' :
       every "fcp_nraise" steps the instantaneous temperature is
       reduced by -"fcp_delta_t" (i.e. "fcp_delta_t" is added to the temperature)

if "fcp_temperature" == 'rescale-v' :
       every "fcp_nraise" steps the average temperature, computed from
       the last "fcp_nraise" steps, is rescaled to "fcp_tempw"

if "fcp_temperature" == 'berendsen' :
       the "rise time" parameter is given in units of the time step:
       tau = fcp_nraise*dt, so dt/tau = 1/fcp_nraise

if "fcp_temperature" == 'andersen' :
       the "collision frequency" parameter is given as nu=1/tau
       defined above, so nu*dt = 1/fcp_nraise
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help freeze_all_atoms -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>freeze_all_atoms</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE., freeze all atoms
to perform relaxation or dynamics only with FCP.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nsolv -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nsolv</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Status: </em> REQUIRED
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The number of solvents (i.e. molecular species) in the unit cell
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help closure -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>closure</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'kh'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre>
Specify the type of closure equation:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'kh'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
The Kovalenko and Hirata's model.
[A.Kovalenko, F.Hirata, JCP 110, 10095 (1999), "doi:10.1063/1.478883"]
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'hnc'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
The HyperNetted-Chain model, which is
suitable only for solvents without charge.
[J.P.Hansen et al., Theory of simple liquids. Academic Press, London, 1990]
            </pre></dd>
</dl>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tempv -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tempv</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 300.D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Temperature (Kelvin) of solvents.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ecutsolv -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ecutsolv</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 4 * "ecutwfc"
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Kinetic energy cutoff (Ry) for solvent's correlation functions.
If a solute is an isolated system or slab, you may allowed to
use default value. For a frameworked or porous solute (e.g. Zeolite, MOF),
it is desirable to apply a larger value. Solvents confined in a framework
often have a high frequency.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help solute_lj -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>solute_lj(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'uff'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre>
Specify the Lennard-Jones potential of solute on atomic type 'i':
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'none'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
The Lennard-Jones potential is not specified here.
you must set "solute_epsilon" and "solute_sigma".
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'uff'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Universal Force Field.
[A.K.Rappe et al., JACS 144, 10024 (1992), "doi:10.1021/ja00051a040"]
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'clayff'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Clay's Force Field
[R.T.Cygan et al., JPC B 108, 1255 (2004), "doi:10.1021/jp0363287"]
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'opls-aa'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
OPLS-AA (generic parameters for QM/MM)
            </pre></dd>
</dl>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help solute_epsilon -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>solute_epsilon(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The Lennard-Jones potential of solute on atomic type 'i'.
Here, you can set the parameter 'epsilon' (kcal/mol).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help solute_sigma -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>solute_sigma(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The Lennard-Jones potential of solute on atomic type 'i'.
Here, you can set the parameter 'sigma' (Angstrom).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help starting1d -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>starting1d</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'zero'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Starting correlation functions of 1D-RISM from zero.
( default for scf, *relax, *md )
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'file'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Start from existing "1d-rism_csvv_r.xml" file in the
directory specified by variables "prefix" and "outdir".
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'fix'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Read from existing "1d-rism_csvv_r.xml" file in the
directory specified by variables "prefix" and "outdir",
and never calculate 1D-RISM.
For nscf and bands calculation this is the default.
            </pre></dd>
</dl>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help starting3d -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>starting3d</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'zero'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Starting correlation functions of 3D-RISM from zero.
( default for scf, *relax, *md )
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'file'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Start from existing "3d-rism_csuv_r.dat" file in the
directory specified by variables "prefix" and "outdir".
For nscf and bands calculation this is the default.
            </pre></dd>
</dl>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help smear1d -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>smear1d</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 2.D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Coulomb smearing radius (a.u.) for 1D-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help smear3d -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>smear3d</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 2.D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Coulomb smearing radius (a.u.) for 3D-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rism1d_maxstep -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rism1d_maxstep</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 50000
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Maximum number of iterations in a 1D-RISM step.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rism3d_maxstep -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rism3d_maxstep</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 5000
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Maximum number of iterations in a 3D-RISM step.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rism1d_conv_thr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rism1d_conv_thr</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.D-8
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Convergence threshold for 1D-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rism3d_conv_thr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rism3d_conv_thr</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em>
1.D-5 if "lgcscf" == .FALSE.;
5.D-6 if "lgcscf" == .TRUE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Convergence threshold for 3D-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help mdiis1d_size -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>mdiis1d_size</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 20
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Size of Modified DIIS (MDIIS) for 1D-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help mdiis3d_size -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>mdiis3d_size</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 10
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Size of Modified DIIS (MDIIS) for 3D-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help mdiis1d_step -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>mdiis1d_step</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.5D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Step of Modified DIIS (MDIIS) for 1D-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help mdiis3d_step -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>mdiis3d_step</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.8D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Step of Modified DIIS (MDIIS) for 3D-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rism1d_bond_width -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rism1d_bond_width</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Gaussian width of bonds to smear intra-molecular correlation for 1D-RISM.
If 3D-RISM calculation, default is 0.
If Laue-RISM calculation, default is 2 / SQRT("ecutwfc").
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rism1d_dielectric -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rism1d_dielectric</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> -1.0D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Dielectric constant for 1D-RISM.
If "rism1d_dielectric" &gt; 0, dielectrically consistent RISM (DRISM) is performed.

For details of DRISM, see:
J.S.Perkyns and B.M.Pettitt, CPL 1992, 190, 626, "doi:10.1016/0009-2614(92)85201-K"
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rism1d_molesize -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rism1d_molesize</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 2.0D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Size of solvent molecules (a.u.) for 1D-RISM.
This is used only if "rism1d_dielectric" &gt; 0.
If you have large molecules, you have to set ~ 20 a.u. .
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rism1d_nproc -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rism1d_nproc</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 128
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Number of processes to calculate 1D-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rism3d_conv_level -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rism3d_conv_level</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em>
0.1 if "laue_both_hands" == .FALSE. .AND. "lgcscf" == .FALSE.;
0.3 if "laue_both_hands" == .FALSE. .AND. "lgcscf" == .TRUE.;
0.5 if "laue_both_hands" == .TRUE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre>
Convergence level of 3D-RISM.
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>0.0</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Convergence level is 'low'.
Convergence threshold of 3D-RISM is greater than
"rism3d_conv_thr", when estimated energy error &gt;&gt; "conv_thr" .
The threshold becomes "rism3d_conv_thr", when
estimated energy error is enough small.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>0.0&lt;x&lt;1.0</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Convergence level is 'medium'.
Convergence threshold of 3D-RISM is intermediate value
between 'low' and 'high', where "rism3d_conv_level" is mixing rate.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>1.0</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Convergence level is 'high'.
Convergence threshold of 3D-RISM is always "rism3d_conv_thr" .
            </pre></dd>
</dl>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rism3d_planar_average -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rism3d_planar_average</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE., planar averages of solvent densities and potentials
are calculated and written to 'prefix.rism1'.
For 3D-RISM, default is .FALSE.
For Laue-RISM, default is .TRUE.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help laue_nfit -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>laue_nfit</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 4
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The number of z-grid points for the polynomial fit along the cell edge.
This is only for Laue-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help laue_expand_right -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>laue_expand_right</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> -1.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If positive value, set the ending position offset [in a.u.]
of the solvent region on right-hand side of the unit cell,
measured relative to the unit cell edge.
(the solvent region ends at z = + [L_z/2 + "laue_expand_right"].)
This is only for Laue-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help laue_expand_left -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>laue_expand_left</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> -1.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If positive value, set the ending position offset [in a.u.]
of the solvent region on left-hand side of the unit cell,
measured relative to the unit cell edge.
(the solvent region ends at z = - [L_z/2 + "laue_expand_left"].)
This is only for Laue-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help laue_starting_right -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>laue_starting_right</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Set the starting position [in a.u.] of the solvent region
on right-hand side of the unit cell. Then the solvent region is
defined as [ "laue_starting_right" , L_z/2 + "laue_expand_right" ],
where distribution functions are finite.
This is only for Laue-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help laue_starting_left -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>laue_starting_left</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Set the starting position [in a.u.] of the solvent region
on left-hand side of the unit cell. Then the solvent region is
defined as [ -L_z/2 - "laue_expand_left" , "laue_starting_left" ],
where distribution functions are finite.
This is only for Laue-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help laue_buffer_right -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>laue_buffer_right</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em>
 8.0 if "laue_expand_right" &gt; 0.0;
-1.0 if "laue_expand_right" &lt;= 0.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If positive value, set the buffering length [in a.u.]
of the solvent region on right-hand side of the unit cell.
Then correlation functions are defined inside of
[ "laue_starting_right" - "laue_buffer_right" , L_z/2 + "laue_expand_right" ].
This is only for Laue-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help laue_buffer_left -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>laue_buffer_left</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em>
 8.0 if "laue_expand_left" &gt; 0.0;
-1.0 if "laue_expand_left" &lt;= 0.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If positive value, set the buffering length [in a.u.]
of the solvent region on left-hand side of the unit cell.
Then correlation functions are defined inside of
[ -L_z/2 - "laue_expand_left" , "laue_starting_left" + "laue_buffer_left" ].
This is only for Laue-RISM.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help laue_both_hands -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>laue_both_hands</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE., you can set different densities
to the solvent regions of right-hand side and left-hand side.
See "SOLVENTS" card.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help laue_wall -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>laue_wall</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'auto'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre>
Set the repulsive wall with (1/r)^12 term of Lennard-Jones potential.
This is only for Laue-RISM.
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'none'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
The repulsive wall is not defined.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'auto'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
The repulsive wall is defined, whose edge position is set automatically.
One does not have to set "laue_wall_z" (the edge position).
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'manual'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
The repulsive wall is defined, whose edge position is set manually.
One have to set "laue_wall_z" (the edge position).
            </pre></dd>
</dl>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help laue_wall_z -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>laue_wall_z</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Set the edge position [in a.u.] of the repulsive wall.
If "laue_expand_right" &gt; 0.0, the repulsive wall is defined on [ -inf , "laue_wall_z" ].
If "laue_expand_left" &gt; 0.0, the repulsive wall is defined on [ "laue_wall_z" , inf ].
This is only for Laue-RISM and "laue_wall" == 'manual' .
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help laue_wall_rho -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>laue_wall_rho</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.01
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The density (1/bohr^3) of the repulsive wall.
This is only for Laue-RISM and "laue_wall" /= 'none' .
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help laue_wall_epsilon -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>laue_wall_epsilon</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The Lennard-Jones potential of the repulsive wall.
Here, you can set the parameter 'epsilon' (kcal/mol).
This is only for Laue-RISM and "laue_wall" /= 'none' .
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help laue_wall_sigma -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>laue_wall_sigma</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 4.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The Lennard-Jones potential of the repulsive wall.
Here, you can set the parameter 'sigma' (Angstrom).
This is only for Laue-RISM and "laue_wall" /= 'none' .
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help laue_wall_lj6 -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>laue_wall_lj6</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE., the attractive term -(1/r)^6 of Lennard-Jones potential is added.
This is only for Laue-RISM and "laue_wall" /= 'none' .
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
      <h2>Description of ATOMIC_POSITIONS card's flags</h2><li> <em>Description:</em>
</li><blockquote>
<pre>
Units for ATOMIC_POSITIONS:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>alat</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
atomic positions are in cartesian coordinates, in
units of the lattice parameter (either celldm(1)
or A). If no option is specified, 'alat' is assumed;
not specifying units is DEPRECATED and will no
longer be allowed in the future
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>bohr</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
atomic positions are in cartesian coordinate,
in atomic units (i.e. Bohr radii)
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>angstrom</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
atomic positions are in cartesian coordinates, in Angstrom
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>crystal</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
atomic positions are in crystal coordinates, i.e.
in relative coordinates of the primitive lattice
vectors as defined either in card "CELL_PARAMETERS"
or via the ibrav + celldm / a,b,c... variables
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>crystal_sg</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
atomic positions are in crystal coordinates, i.e.
in relative coordinates of the primitive lattice.
This option differs from the previous one because
in this case only the symmetry inequivalent atoms
are given. The variable "space_group" must indicate
the space group number used to find the symmetry
equivalent atoms. The other variables that control
this option are uniqueb, origin_choice, and
rhombohedral.
            </pre></dd>
</dl>
</blockquote>
      
}


# ------------------------------------------------------------------------
help atomic_coordinates -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variable: </em><big><b>X</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> label of the atom as specified in "ATOMIC_SPECIES"
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
      in the following four forms (Wyckoff positions):
         C  1a
         C  8g   x
         C  24m  x y
         C  48n  x y z
      The first form must be used when the Wyckoff letter determines uniquely
      all three coordinates, forms 2,3,4 when the Wyckoff letter and 1,2,3
      coordinates respectively are needed.

      The forms:
         C 8g  x  x  x
         C 24m x  x  y
      are not allowed, but
         C x x x
         C x x y
         C x y z
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
      <h2>Description of K_POINTS card's flags</h2><li> <em>Description:</em>
</li><blockquote>
<pre>
K_POINTS options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>tpiba</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
read k-points in cartesian coordinates,
in units of 2 pi/a (default)
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>automatic</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
automatically generated uniform grid of k-points, i.e,
generates ( nk1, nk2, nk3 ) grid with ( sk1, sk2, sk3 ) offset.
nk1, nk2, nk3 as in Monkhorst-Pack grids
k1, k2, k3 must be 0 ( no offset ) or 1 ( grid displaced
by half a grid step in the corresponding direction )
BEWARE: only grids having the full symmetry of the crystal
        work with tetrahedra. Some grids with offset may not work.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>crystal</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
read k-points in crystal coordinates, i.e. in relative
coordinates of the reciprocal lattice vectors
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>gamma</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
use k = 0 (no need to list k-point specifications after card)
In this case wavefunctions can be chosen as real,
and specialized subroutines optimized for calculations
at the gamma point are used (memory and cpu requirements
are reduced by approximately one half).
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>tpiba_b</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Used for band-structure plots.
See Doc/brillouin_zones.pdf for usage of BZ labels;
otherwise, k-points are in units of  2 pi/a.
nks points specify nks-1 lines in reciprocal space.
Every couple of points identifies the initial and
final point of a line. pw.x generates N intermediate
points of the line where N is the weight of the first point.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>crystal_b</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
As tpiba_b, but k-points are in crystal coordinates.
See Doc/brillouin_zones.pdf for usage of BZ labels.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>tpiba_c</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Used for band-structure contour plots.
k-points are in units of  2 <i>pi/a.</i> nks must be 3.
3 k-points k_0, k_1, and k_2 specify a rectangle
in reciprocal space of vertices k_0, k_1, k_2,
k_1 + k_2 - k_0: k_0 + \alpha (k_1-k_0)+
\beta (k_2-k_0) with 0 &lt;\alpha,\beta &lt; 1.
The code produces a uniform mesh n1 x n2
k points in this rectangle. n1 and n2 are
the weights of k_1 and k_2. The weight of k_0
is not used.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>crystal_c</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
As tpiba_c, but k-points are in crystal coordinates.
            </pre></dd>
</dl>
</blockquote>
      
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
    

    <ul>
<li> <em>Variables: </em><big><b>k_x, k_y, k_z, wk_</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
for the respective explanation, see the "xk_x", "xk_y", "xk_z", "wk"
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
help ADDITIONAL_K_POINTS_flags -helpfmt helpdoc -helptext {
      <h2>Description of ADDITIONAL_K_POINTS card's flags</h2><li> <em>Description:</em>
</li><blockquote><pre>
for the explanation of the K_POINTS' options, see "K_POINTS"
         </pre></blockquote>
      
}


# ------------------------------------------------------------------------
help nks_add -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nks_add</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Number of supplied "additional" k-points.
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help CELL_PARAMETERS_flags -helpfmt helpdoc -helptext {
      <h2>Description of CELL_PARAMETERS card's flags</h2><li> <em>Description:</em>
</li><blockquote><pre>
Unit for lattice vectors; options are:

<b>'bohr'</b> / <b>'angstrom':</b>
                     lattice vectors in bohr-radii / angstrom.
                     In this case the lattice parameter alat = sqrt(v1*v1).

<b>'alat'</b> / nothing specified:
                     lattice vectors in units of the lattice parameter (either
                     "celldm"(1) or "A"). Not specifying units is DEPRECATED
                     and will not be allowed in the future.

If neither unit nor lattice parameter are specified,
'bohr' is assumed - DEPRECATED, will no longer be allowed
         </pre></blockquote>
      
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
<li> <em>Variables: </em><big><b>constr(1), constr(2), constr(3), constr(4)</b></big>
</li>
<br><li> <em>Type: </em>
</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
These variables have different meanings for different constraint types:

<b>'type_coord'</b> :
               <i>constr(1)</i> is the first index of the atomic type involved
               <i>constr(2)</i> is the second index of the atomic type involved
               <i>constr(3)</i> is the cut-off radius for estimating the coordination
               <i>constr(4)</i> is a smoothing parameter

<b>'atom_coord'</b> :
               <i>constr(1)</i> is the atom index of the atom with constrained coordination
               <i>constr(2)</i> is the index of the atomic type involved in the coordination
               <i>constr(3)</i> is the cut-off radius for estimating the coordination
               <i>constr(4)</i> is a smoothing parameter

<b>'distance'</b> :
               atoms indices object of the constraint, as they appear in
               the "ATOMIC_POSITIONS" card

<b>'planar_angle',</b> <b>'torsional_angle'</b> :
               atoms indices object of the constraint, as they appear in the
               "ATOMIC_POSITIONS" card (beware the order)

<b>'bennett_proj'</b> :
               <i>constr(1)</i> is the index of the atom whose position is constrained.
               <i>constr(2:4)</i> are the three coordinates of the vector that specifies
               the constraint direction.
<b>'potential_wall'</b> :
               Formula is: External force = prefac * exponent * Exp(-exponent). Force is only applied
               on atoms within the cutoff.
               <i>constr(1)</i> is the prefactor
               <i>constr(2)</i> is the value in the exponent
               <i>constr(3)</i> is the cutoff (in a.u.)
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
help atomic_velocities -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variable: </em><big><b>V</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> label of the atom as specified in ATOMIC_SPECIES
                  </pre></blockquote>
</ul><ul>
<li> <em>Variables: </em><big><b>vx, vy, vz</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> atomic velocities along x y and z direction
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
<blockquote><pre> label of the atom as specified in "ATOMIC_SPECIES"
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


# ------------------------------------------------------------------------
help SOLVENTS_flags -helpfmt helpdoc -helptext {
      <h2>Description of SOLVENTS card's flags</h2><li> <em>Description:</em>
</li><blockquote>
<dl style="margin-left: 1.5em;">
<dt><tt><b>1/cell</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
solvent's densities are specified
as number of molecules in the unit cell.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>mol/L</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
solvent's densities are specified as molar concentrations.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>g/cm^3</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
solvent's densities are in gram per cm^3.
            </pre></dd>
</dl>
</blockquote>
      
}


# ------------------------------------------------------------------------
help HUBBARD_flags -helpfmt helpdoc -helptext {
      <h2>Description of HUBBARD card's flags</h2><li> <em>Description:</em>
</li><blockquote>
<pre>
<b>HUBBARD</b> options are:
            </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>atomic</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
use atomic orbitals (read from pseudopotential) to build the
Hubbard projectors
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>ortho-atomic</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
use Lowdin orthogonalized atomic orbitals. This option is
recommended to be used whenever possible instead of atomic
because it allows to avoid applying Hubbard corrections twice
in the orbital overlap regions.
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>norm-atomic</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
Lowdin normalization of atomic orbitals. Keep in mind:
atomic orbitals are not orthogonalized in this case.
This is a "quick and dirty" trick to be used when
atomic orbitals from the pseudopotential are not
normalized (and thus produce occupation whose
value exceeds unity).
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>wf</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
use Wannier functions to built Hubbard projectors.
The information about the Wannier functionas are read
from file "prefix".hub that must be generated using pmw.x
(see PP/src/poormanwannier.f90 for details).
Note: these are not maximally localized Wannier functions.
(see PP/examples/example05)
            </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>pseudo</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
use the pseudopotential projectors. The charge density
outside the atomic core radii is excluded.
N.B.: for atoms with +U, a pseudopotential with the
all-electron atomic orbitals are required (i.e.,
as generated by ld1.x with lsave_wfc flag).
            </pre></dd>
</dl>
<pre>
NB: forces and stress are currently implemented only for the
'atomic', 'ortho-atomic', and 'pseudo' Hubbard projectors.
            </pre>
<pre>
Check Doc/Hubbard_input.pdf to see how to specify Hubbard parameters
U, J0, J, B, E2, E3, V in the HUBBARD card.
            </pre>
</blockquote>
      
}

