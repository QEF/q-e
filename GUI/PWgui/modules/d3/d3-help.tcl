
#
# Help-file automatically created by helpdoc utility
#
#    !!! DO NOT EDIT: CHANGES WILL BE LOST !!!
#
	

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
<blockquote><pre>
The file containing the variation of the charge
density at the q point under consideration, this
file is produced by phonon.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fild0rho -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fild0rho</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> ' '
            </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The file containing the variation of the charge
density at q=0, this file is produced by phonon.
            </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help amass -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>amass(i), i=1,ntyp</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
ionic masses [atomic mass units]
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
<blockquote><pre> prefix for file names
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
Directory containing input, output, and scratch files;
must be the same as specified in the calculation of
the unperturbed system and for phonon calculation.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fildyn -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fildyn</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> 'd3dyn'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
The file where the derivative of the dynamical
matrix will be written
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help ethr_ph -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>ethr_ph</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> 1.0d-5
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Threshold for iterative diagonalization
(accuracy in ryd of the calculated eigenvalues).
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help q0mode_todo -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variables: </em><big><b>q0mode_todo(i), i=1,3*nat</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Status: </em> q0mode_todo is statically allocated to dimension 300
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
This array contains the list of the q=0 modes that
will be computed. If q0mode_todo(1).eq.0 the
program will compute every q=0 mode.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help wraux -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>wraux</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If .true. the program will write different terms
of the matrix on different files.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help recv -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>recv</b></big>
</li>
<br><li> <em>Type: </em>LOGICAL</li>
<br><li> <em>Default: </em> .false.
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> Specify .true. for a recover run.
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help istop -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>istop</b></big>
</li>
<br><li> <em>Type: </em>INTEGER</li>
<br><li> <em>Default: </em> 0
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
If this number is set different from zero the
program will stop after the specified routine
and will write the partial result in the recover
file.
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
<blockquote><pre> type of printing ( 0 few, 1 all )
         </pre></blockquote>
</ul>      
      
}

