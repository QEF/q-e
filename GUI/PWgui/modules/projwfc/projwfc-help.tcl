help prefix  -vartype character -helpfmt txt2html -helptext {    
     prefix of input files saved by program pwscf 
<p> ( default = 'pwscf' )
}


help outdir  -vartype character -helpfmt txt2html -helptext {    
     temporary directory where files resides 
<p> ( default = './' )
}

help io_choice  -vartype character -helpfmt txt2html -helptext { 
            'standard' : write projections to standard output 

            'file'     : write projected DOS, one file per atomic 
                         wavefunction (those that have been read 
                         from the pseudopotential file) 

            'both'     : do both of the above things (DEFAULT)
<p> ( default = 'both' )
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

help smoothing -vartype character -helpfmt txt2html -helptext {  
     gaussian broadening (eV, DEFAULT: DeltaE)
}

