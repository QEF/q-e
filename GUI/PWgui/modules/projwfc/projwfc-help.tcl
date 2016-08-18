
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
Prefix of input file produced by pw.x
(wavefunctions are needed).
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
directory containing the input data, i.e. the same as in pw.x
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ngauss -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ngauss</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Type of gaussian broadening:
    0 ... Simple Gaussian (default)
    1 ... Methfessel-Paxton of order 1
   -1 ... Marzari-Vanderbilt "cold smearing"
  -99 ... Fermi-Dirac function
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help degauss -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>degauss</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 0.0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> gaussian broadening, Ry (not eV!)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {Emin Emax} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>Emin, Emax</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> (band extrema)
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> min &amp; max energy (eV) for DOS plot
         </pre></blockquote>
</ul>
    
}


# ------------------------------------------------------------------------
help DeltaE -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>DeltaE</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> energy grid step (eV)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lsym -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lsym</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .true.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true. the projections are symmetrized
if .false. partial density of state cannot be computed
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help pawproj -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>pawproj</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true. use PAW projectors and all-electron PAW basis
functions to calculate weight factors for the partial
densities of states. Following Bloechl, "PRB 50, 17953 (1994)",
Eq. (4 &amp; 6), the weight factors thus approximate the real
charge within the augmentation sphere of each atom.
Only for PAW, not implemented in the noncolinear case.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help filpdos -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>filpdos</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> (value of "prefix" variable)
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> prefix for output files containing PDOS(E)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help filproj -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>filproj</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> (standard output)
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
file containing the projections
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lwrite_overlaps -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lwrite_overlaps</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true., the overlap matrix of the atomic orbitals
prior to orthogonalization is written to the atomic_proj datafile.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lbinary_data -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lbinary_data</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true., the atomic_proj datafile is written in binary fmt.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help kresolveddos -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>kresolveddos</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true. the k-resolved DOS is computed: not summed over
all k-points but written as a function of the k-point index.
In this case all k-point weights are set to unity
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help tdosinboxes -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>tdosinboxes</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true. compute the local DOS integrated in volumes
volumes are defined as boxes with edges parallel
to the unit cell, containing the points of the
(charge density) FFT grid included within
"irmin" and "irmax", in the three dimensions:
from "irmin"(j,n) to "irmax"(j,n) for j=1,2,3
(n=1,"n_proj_boxes")
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help n_proj_boxes -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>n_proj_boxes</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
number of boxes where the local DOS is computed
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help irmin -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>irmin(3,n_proj_boxes)</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 1 for each box
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
first point of the given box

BEWARE: irmin is a 2D array of the form: irmin(3,"n_proj_boxes")
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help irmax -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>irmax(3,n_proj_boxes)</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0 for each box
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
last point of the given box;
( 0 stands for the last point in the FFT grid )

BEWARE: irmax is a 2D array of the form: irmax(3,"n_proj_boxes")
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help plotboxes -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>plotboxes</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if .true., the boxes are written in output as
as xsf files with 3D datagrids, valued 1.0
inside the box volume and 0 outside
(visualize them as isosurfaces with isovalue 0.5)
         </pre></blockquote>
</ul>      
      
}

