
#
# Help-file automatically created by helpdoc utility
#
#    !!! DO NOT EDIT: CHANGES WILL BE LOST !!!
#
	

# ------------------------------------------------------------------------
help string_method -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>string_method</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'neb'
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote>
<pre>
A string describing the task to be performed. Options are:
                  </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'neb'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> nudget-elastic-band
                  </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'smd'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> string-method-dynamics
                  </pre></dd>
</dl>
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
<pre> Options are:
                  </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'from_scratch'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> from scratch
                  </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'restart'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;"> from previous interrupted run
                  </pre></dd>
</dl>
</blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help nstep_path -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>nstep_path</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em>
1
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
number of ionic + electronic steps
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
<blockquote>
<pre>
Specify the type of optimization scheme:
                  </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'sd'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
steepest descent
                  </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'broyden'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
quasi-Newton Broyden's second method (suggested)
                  </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'broyden2'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
another variant of the quasi-Newton Broyden's
second method to be tested and compared with the
previous one.
                  </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'quick-min'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
an optimisation algorithm based on the
projected velocity Verlet scheme
                  </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'langevin'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
finite temperature langevin dynamics of the
string (smd only). It is used to compute the
average path and the free-energy profile.
                  </pre></dd>
</dl>
</blockquote>
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
<blockquote>
<pre>
Specify the type of Climbing Image scheme:
                  </pre>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'no-CI'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
climbing image is not used
                  </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'auto'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
original CI scheme. The image highest in energy
does not feel the effect of springs and is
allowed to climb along the path
                  </pre></dd>
</dl>
<dl style="margin-left: 1.5em;">
<dt><tt><b>'manual'</b> :</tt></dt>
<dd><pre style="margin-top: 0em; margin-bottom: -1em;">
images that have to climb are manually selected.
See also "CLIMBING_IMAGES" card
                  </pre></dd>
</dl>
</blockquote>
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
Also the first and the last configurations are optimized
"on the fly" (these images do not feel the effect of the springs).
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help minimum_image -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>minimum_image</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Assume a "minimum image criterion" to build the path. If an atom
moves by more than half the length of a crystal axis between one
image and the next in the input (before interpolation),
an appropriate periodic replica of that atom is chosen.
Useful to avoid jumps in the initial reaction path.
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
If "opt_scheme"=="broyden", ds is used as a guess for the
diagonal part of the Jacobian matrix.
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
mass-weighted coordinates. Useful together with quick-min
optimization scheme, if some bonds are much stiffer than
others. By assigning a larger (fictitious) mass to atoms
with stiff bonds, one may use a longer time step "ds"
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
help lfcpopt -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lfcpopt</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .FALSE.
               </li>
<br><li> <em>See: </em> fcp_mu
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .TRUE. perform a constant bias potential (constant-mu)
calculation with ESM method (assume_isolated = 'esm' and
esm_bc = 'bc2' or 'bc3' must be set in SYSTEM namelist).
"fcp_mu" gives the target Fermi energy.
See the header of PW/src/fcp.f90 for documentation
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fcp_mu -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fcp_mu</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.d0
               </li>
<br><li> <em>See: </em> lfcpopt
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If "lfcpopt" == .TRUE., gives the target Fermi energy [Ry].
One can specify the total charge of the system for the first
and last image by giving "fcp_tot_charge_first" and "fcp_tot_charge_last"
so that the Fermi energy of these systems will be the target value,
otherwise "first_last_opt" should be .TRUE.
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fcp_tot_charge_first -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fcp_tot_charge_first</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.d0
               </li>
<br><li> <em>See: </em> lfcpopt
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Total charge of the system ('tot_charge') for the first image.
Initial 'tot_charge' for intermediate images will be given by
linear interpolation of "fcp_tot_charge_first" and "fcp_tot_charge_last"
               </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fcp_tot_charge_last -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fcp_tot_charge_last</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.d0
               </li>
<br><li> <em>See: </em> lfcpopt
               </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Total charge of the system ('tot_charge') for the last image.
Initial 'tot_charge' for intermediate images will be given by
linear interpolation of "fcp_tot_charge_first" and "fcp_tot_charge_last"
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

