module dimension -title "Test namelist's dimensions" -script {

    namelist dim -name "DIMENSION" {
	
	var howmany -label "Size of dimension:" -widget spinint -default 3
	
	dimension myDim1 -start 1 -end   3	
	
	dimension myDim2 {
	    -variable  textvalue-test
	    -widget    radiobox
	    -start     1
	    -end       3
	    -value     { first second third }	    
	    -textvalue { 1st 2nd 3rd }
	    -default   2nd
	}
    }

    tracevar howmany w {
	widgetconfigure myDim1 -end [varvalue howmany]
    }
}

    