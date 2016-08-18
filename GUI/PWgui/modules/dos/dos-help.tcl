
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
prefix of input file produced by pw.x
(wavefunctions are not needed)
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
<br><li> <em>Status: </em> optional
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
Type of gaussian broadening:

    =  0  Simple Gaussian (default)

    =  1  Methfessel-Paxton of order 1

    = -1  Marzari-Vanderbilt "cold smearing"

    =-99  Fermi-Dirac function
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help degauss -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>degauss</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
gaussian broadening, Ry (not eV!)
(see below)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
grouphelp {Emin Emax} -helpfmt helpdoc -helptext {
    <ul>
<li> <em>Variables: </em><big><b>Emin, Emax</b></big>
</li>
<br><li> <em>Type: </em>REAL</li>
<br><li> <em>Default: </em> band extrema
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
min, max energy (eV) for DOS plot. If unspecified, the
lower and/or upper band value, plus/minus 3 times the
value of the gaussian smearing if present, will be used.
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
<blockquote><pre>
energy grid step (eV)
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help fildos -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>fildos</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> '"prefix".dos'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre>
output file containing DOS(E)
         </pre></blockquote>
</ul>      
      
}

