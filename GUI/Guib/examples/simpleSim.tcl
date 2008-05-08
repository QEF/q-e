module simpleSim -title "Simulation Setup" -script {
    
    line title -name "Title" {
	var jobtitle -label "Job title:"  -fmt %S	
    }
    namelist pars -name Parameters {
	var type {
	    -label    "Job type:"
	    -value    {"'Single-point calculation'" "'Structural optimization'"}
	    -widget   optionmenu
	    -default  "'Single-point calculation'"
	}
	var functional {
	    -label    "DFT Functional:"
	    -value    {'PBE' 'RPBE' 'B3LYP' }
	    -widget   optionmenu
	    -default  'B3LYP'
	}
	var basis {
	    -label    "Gaussian basis set:"
	    -value    {'STO-3G' '3-21G' '6-21G' '6-311G' '6-311G*'}
	    -widget   optionmenu
	    -default  '6-311G'
	}
	var natoms -label "natoms"  -widget spinint -validate posint -default 1
    }
    
    keyword coord ATOMIC_COORDINATES\n
    table atoms {
	-caption  "Enter atomic coorditanes"
	-head     {"Atomic symbol" X-Coordinate Y-Coordinate Z-Coordinate}
	-cols     4
	-rows     1
	-outfmt   "%3s %15.10f %15.10f %15.10f"
    }
    keyword coord_end END\n    

    tracevar natoms w {
	widgetconfigure atoms -rows [varvalue natoms]
    }
}
