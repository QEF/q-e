
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
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
temporary directory where pw.x files resides
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

