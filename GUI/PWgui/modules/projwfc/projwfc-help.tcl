
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
<br><li> <em>Default: </em> './'
         </li>
<br><li> <em>Description:</em>
</li>
<blockquote><pre> directory containing the input file
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
if true the projections are symmetrized
         </pre></blockquote>
</ul>      
      
}


# ------------------------------------------------------------------------
help filpdos -helpfmt helpdoc -helptext {
      <ul>
<li> <em>Variable: </em><big><b>filpdos</b></big>
</li>
<br><li> <em>Type: </em>CHARACTER</li>
<br><li> <em>Default: </em> (value of prefix variable)
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

