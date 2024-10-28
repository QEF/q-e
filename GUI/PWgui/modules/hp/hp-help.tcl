
#
# Help-file automatically created by helpdoc utility
#
#    !!! DO NOT EDIT: CHANGES WILL BE LOST !!!
#
	

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
help iverbosity -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>iverbosity</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
= 1 : minimal output
= 2 : as above + symmetry matrices, final response
      matrices chi0 and chi1 and their inverse matrices,
      full U matrix
= 3 : as above + various detailed info about the NSCF
      calculation at k and k+q
= 4 : as above + response occupation matrices at every
      iteration and for every q point in the star
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
<blockquote><pre>
Maximum allowed run time before the job stops smoothly.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {nq1 nq2 nq3} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>nq1, nq2, nq3</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1,1,1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Parameters of the Monkhorst-Pack grid (no offset).
Same meaning as for nk1, nk2, nk3 in the input of pw.x.
         </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
help skip_equivalence_q -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>skip_equivalence_q</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. then the HP code will skip the equivalence
analysis of q points, and thus the full grid of q points
will be used. Otherwise the symmetry is used to determine
equivalent q points (star of q), and then perform
calculations only for inequivalent q points.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help determine_num_pert_only -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>determine_num_pert_only</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>See: </em> find_atpert
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. determines the number of perturbations
(i.e. which atoms will be perturbed) and exits smoothly
without performing any calculation. For DFT+U+V, it also
determines the indices of inter-site couples.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help determine_q_mesh_only -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>determine_q_mesh_only</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>See: </em> perturb_only_atom
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. determines the number of q points
for a given perturbed atom and exits smoothly.
This keyword can be used only if perturb_only_atom
is set to .true.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help find_atpert -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>find_atpert</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Method for searching of atoms which must be perturbed.
1 = Find how many inequivalent Hubbard atoms there are
    by analyzing unperturbed occupations.
2 = Find how many Hubbard atoms to perturb based on
    how many different Hubbard atomic types there are.
    Warning: atoms which have the same type but which
    are inequivalent by symmetry or which have different
    occupations will not be distinguished in this case
    (use option 1 or 3 instead).
3 = Find how many inequivalent Hubbard atoms
    there are using symmetry. Atoms which have the
    same type but are not equivalent by symmetry will
    be distinguished in this case.
4 = Perturb all Hubbard atoms (the most expensive option)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help docc_thr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>docc_thr</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 5.D-5
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Threshold for a comparison of unperturbed occupations
which is needed for the selection of atoms which must
be perturbed. Can be used only when "find_atpert" = 1.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help skip_type -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>skip_type(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> skip_type(i) = .false.
         </li>
<br><li> <em>See: </em> equiv_type
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
"skip_type"(i), where i runs over types of atoms.
If "skip_type"(i)=.true. then no linear-response
calculation will be performed for the i-th atomic type:
in this case "equiv_type"(i) must be specified, otherwise
the HP code will stop. This option is useful if the
system has atoms of the same type but opposite spin
pollarizations (anti-ferromagnetic case).
This keyword cannot be used when "find_atpert" = 1.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help equiv_type -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>equiv_type(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> equiv_type(i) = 0
         </li>
<br><li> <em>See: </em> skip_type
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
"equiv_type"(i), where i runs over types of atoms.
"equiv_type"(i)=j, will make type i equivalent to type j
(useful when nspin=2). Such a merging of types is done
only at the post-processing stage.
This keyword cannot be used when "find_atpert" = 1.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help perturb_only_atom -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>perturb_only_atom(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> perturb_only_atom(i) = .false.
         </li>
<br><li> <em>See: </em> compute_hp
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If "perturb_only_atom"(i)=.true. then only the i-th
atom will be perturbed and considered in the run.
This variable is useful when one wants to split
the whole calculation on parts.

<b>Note:</b> this variable has a higher priority than "skip_type".
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
<br><li> <em>See: </em> last_q, sum_pertq
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Computes only the q points from "start_q" to "last_q".

<b>IMPORTANT:</b> "start_q" must be smaller or equal to
the total number of q points found.
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
<br><li> <em>See: </em> start_q, sum_pertq
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Computes only the q points from "start_q" to "last_q".

<b>IMPORTANT:</b> "last_q" must be smaller or equal to
the total number of q points found.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help sum_pertq -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>sum_pertq</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>See: </em> start_q, last_q, perturb_only_atom
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If it is set to .true. then the HP code will collect
pieces of the response occupation matrices for all
q points. This variable should be used only when
"start_q", "last_q" and "perturb_only_atom" are used.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help compute_hp -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>compute_hp</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>See: </em> perturb_only_atom
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If it is set to .true. then the HP code will collect
pieces of the chi0 and chi matrices (which must have
been produced in previous runs) and then compute
Hubbard parameters. The HP code will look for files
tmp_dir/HP/prefix.chi.i.dat. Note that all files
prefix.chi.i.dat (where i runs over all perturbed
atoms) must be placed in one folder tmp_dir/HP/.
"compute_hp"=.true. must be used only when the
calculation was parallelized over perturbations.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help conv_thr_chi -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>conv_thr_chi</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.D-5
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Convergence threshold for the response function chi,
which is defined as a trace of the response
occupation matrix.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help thresh_init -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>thresh_init</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.D-14
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Initial threshold for the solution of the linear
system (first iteration). Needed to converge the
bare (non-interacting) response function chi0.
The specified value will be multiplied by the
number of electrons in the system.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ethr_nscf -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ethr_nscf</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.D-11
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Threshold for the convergence of eigenvalues during
the iterative diagonalization of the Hamiltonian in
the non-self-consistent-field (NSCF) calculation at
k and k+q points. Note, this quantity is NOT extensive.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help niter_max -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>niter_max</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 100
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Maximum number of iterations in the iterative
solution of the linear-response Kohn-Sham equations.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help alpha_mix -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>alpha_mix(i)</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> alpha_mix(1)=0.3
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Mixing parameter (for the i-th iteration) for updating
the response SCF potential using the modified Broyden
method. See: D.D. Johnson, "PRB 38, 12807 (1988)".
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nmix -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nmix</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 4
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Number of iterations used in potential mixing
using the modified Broyden method. See:
D.D. Johnson, "PRB 38, 12807 (1988)".
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help num_neigh -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>num_neigh</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 6
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Number of nearest neighbors of every Hubbard atom which
will be considered when writting Hubbard V parameters to
the file parameters.out, which can be used in the
subsequent DFT+U+V calculation. This keyword is used only
for DFT+U+V (post-processing stage).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lmin -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lmin</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 2
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Minimum value of the orbital quantum number of the Hubbard
atoms starting from which (and up to the maximum l in the
system) Hubbard V will be written to the file parameters.out.
"lmin" refers to the orbital quantum number of the atom
corresponding to the first site-index in Hubbard_V(:,:,:).
This keyword is used only for DFT+U+V and only
in the post-processing stage. Example: "lmin"=1 corresponds to
writing to file V between e.g. oxygen (with p states) and its
neighbors, and including V between transition metals (with d
states) and their neighbors. Instead, when "lmin"=2 only the
latter will be written to parameters.out.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help rmax -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>rmax</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 100.D0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Maximum distance (in Bohr) between two atoms to search
neighbors (used only at the postprocessing step for
DFT+U+V). This keyword is useful when there
are e.g. defects in the system.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help dist_thr -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>dist_thr</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 6.D-4
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Threshold (in Bohr) for comparing inter-atomic distances
when reconstructing the missing elements of the response
susceptibility in the post-processing step.
         </pre></blockquote>
</ul>      
      
}

