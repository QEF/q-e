
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
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
prefix of files saved by program pw.x
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
directory containing the input data,
i.e. the same as in pw.x
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help filband -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>filband</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
file "filband" contains the bands
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help spin_component -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>spin_component</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
In the lsda case select:

   1 = spin-up
   2 = spin-down
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lsigma -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>lsigma(i), i=1,3</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If true writes a file filband.i with the expectation
values of the spin operator on the spinor wave-functions.
(only in the noncollinear case).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lp -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lp</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. matrix elements of the momentum operator p between
conduction and valence bands are written to file "filp"
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help filp -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>filp</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'p_avg.dat'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If "lp" is set to .true., file name for matrix elements of p
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help lsym -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>lsym</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. the bands are classified according to the
irreducible representations of the small group of k. A
file "filband".rap with the same format of "filband"
is written.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help no_overlap -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>no_overlap</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. writes the eigenvalues in the output file
without changing their order.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help plot_2d -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>plot_2d</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. writes the eigenvalues in the output file
in a 2D format readable by gnuplot. Band ordering is not
changed. Each band is written in a different file called
filband.# with the format:
xk, yk, energy
xk, yk, energy
..  ..  ..
energies are written in eV and xk in units 2\pi/a.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {firstk lastk} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>firstk, lastk</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
if lsym=.true. makes the symmetry analysis only for k
points between firstk to lastk
         </pre></blockquote>
</ul>
    
}

