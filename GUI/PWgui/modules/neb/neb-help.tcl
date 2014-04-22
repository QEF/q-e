
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
<blockquote><pre>
a string describing the task to be performed:
   'neb',
   'smd'
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

'restart'       : from previous interrupted run
         </pre></blockquote>
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
<blockquote><pre>
Specify the type of optimization scheme:

'sd'         : steepest descent

'broyden'    : quasi-Newton Broyden's second method (suggested)

'broyden2'   : another variant of the quasi-Newton Broyden's
               second method to be tested and compared with the
               previous one.

'quick-min'  : an optimisation algorithm based on the
               projected velocity Verlet scheme

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
If opt_scheme="broyden", ds is used as a guess for the
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

