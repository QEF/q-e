help prefix  -vartype character -helpfmt txt2html -helptext {    
     prefix of input files saved by program pwscf 
<p> ( default = 'pwscf' )
}


help outdir  -vartype character -helpfmt txt2html -helptext {    
     temporary directory where files resides 
<p> ( default = './' )
}

help ngauss -vartype integer -helpfmt txt2html -helptext {
     type of gaussian broadening (optional)    ( default = 0 )
            =  0  Simple Gaussian (default)
            =  1  Methfessel-Paxton of order 1
            = -1  Marzari-Vanderbilt "cold smearing"
            = 99  Fermi-Dirac function
}

help degauss -vartype real -helpfmt txt2html -helptext {
     gaussian broadening in Ry (not eV!)
<p> ( default = 0.0 )
}

set _e {
     projected DOS is plotted starting at E=Emin  
     (eV, DEFAULT: lowest eigenvalue)... 
     ...up to E=Emax (eV, DEFAULT: highest eigenvalue)... 
     ...in steps of DeltaE (eV, DEFAULT: 0.01) 
}

foreach var {Emin Emax DeltaE} {
    help $var -vartype real -helpfmt txt2html -helptext $_e
}

help filpdos -vartype character -helpfmt txt2html -helptext {  
     prefix for output files containing PDOS(E)
<p> ( default = using the value of <i>prefix</i> )
}
