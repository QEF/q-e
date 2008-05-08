module validate -title "Validation procs test" -script {

    # support for undefined variables in the input files
    set ::guib::settings(NAMELIST.variable_support_undefined) 1

    namelist validate -name "VALIDATION" {
	foreach type {
	    whatever          
	    nonnegint         
	    posint            
	    nonposint         
	    negint            
	    nonnegreal        
	    posreal	      
	    nonposreal	      
	    negreal	      
	    fortranreal	      
	    fortrannonnegreal 
	    fortranposreal    
	    fortrannonposreal 
	    fortrannegreal    
	    string
	} {
	    var $type -label [string totitle $type]: -validate $type
	}
    }
}

