help fildrho -vartype character -helpfmt txt2html -helptext {
              The file containing the variation of the charge     
              density at the q point under consideration, this
              file is produced by phonon.
<p> ( default = ' ' )
}
	       
help fild0rho -vartype character -helpfmt txt2html -helptext {
              the file containing the variation of the charge 
              density at q=0, this file is produced by phonon
<p> ( default = ' ' )
}

help amass -vartype real -helpfmt txt2html -helptext {
              ionic masses [atomic units]                         
}


help prefix -vartype character -helpfmt txt2html -helptext {
              prefix for file names                               
<p> ( default = 'pwscf' )
}

help outdir -vartype character -helpfmt txt2html -helptext {
              scratch directory                                   
<p> ( default = './' )
}
     
help fildyn -vartype character -helpfmt txt2html -helptext {
              the file where the derivative of the dynamical      
              matrix will be written
<p> ( default = 'd3dyn' )
}

help ethr_ph -vartype real -helpfmt txt2html -helptext {      
              threshold for iterative diagonalization             
              (accuracy in ryd of the calculated eigenvalues)
              For scf calculations:
<p> ( default = 1.0d-5 )
}
	      
help q0mode_todo -vartype integer -helpfmt txt2html -helptext { 
              This array contains the list of the q=0 mode that   
              will be computed. If q0mode_todo(1).eq.0 the
              program will compute every q=0 mode.
<p> ( default = 0 )
}

help wraux -vartype logical -helpfmt txt2html -helptext {
              if .true. the program will write different terms    
              of the matrix on different files.
<p> ( default = .false. )
}
	      
help recv -vartype logical -helpfmt txt2html -helptext {
              if .true. this will be a recover run.               
<p> ( default = .false. )
}
	      
help istop -vartype integer -helpfmt txt2html -helptext { 
              If this number is set different from zero the       
              program will stop after the specified routine
              and will write the partial result in the recover
              file.
<p> ( default = 0 )
}

help iverbosity -vartype integer -helpfmt txt2html -helptext {
    type of printing to stdout ( 0 few, 1 all )                   
<p> ( default = 0 )
}
	      

