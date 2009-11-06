#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

proc ::main {argc argv} {

global reference_file
global output_file

	get_commandline $argc $argv
	puts -nonewline "Checking $output_file using $reference_file : "
	look_for_bads
	compare_initial_vectors
        compare_chains
	puts " \[OK\]"

}
##########################
#A procedure to read commandline
##########################
proc ::get_commandline {argc argv} {

global reference_file
global output_file


        if { $argc < 1} {
                puts stderr "Usage: check_tddfpt <pw output file> <reference file>"
		exit 1
        }

	if { $argc > 0 } {
		set state first
		foreach arg $argv {
   			switch -- $state {
      		first {
			set state second
			#check if the infile exists
			set output_file $arg
		        if { [file exists $output_file] } {
		        } else {
                		puts stderr "Input file $output_file does not exits"
		                exit 1
        		}

		}
      		second {
			set state flag
			#check if the infile exists
			set reference_file $arg
		        if { [file exists $reference_file] } {
		        } else {
                		puts stderr "Reference file $reference_file does not exits"
		                exit 1
        		}

		}

		flag {
         		switch -glob -- $arg {
            		default {
				puts stderr "unknown flag $arg" 
 				exit 2
				}
         		}
      		}
      		

                        }
                }
      }
}
#Compares the initial Lanczos vectors
proc ::compare_initial_vectors {} {

global reference_file
global output_file

global lc_pos
 
        set lc_pos 0 
       
	if [catch {open $reference_file r} fileId] {
                puts stderr "Cannot open $reference_file: $fileId"
                exit 1
        }
        set reffile [read $fileId]
        close $fileId
	if [catch {open $output_file r} fileId] {
                puts stderr "Cannot open $output_file: $fileId"
                exit 1
        }
        set outfile [read $fileId]
        close $fileId


         foreach {lc_pos} [regexp -all -line -inline -indices {Norm of initial Lanczos vectors=\s*-?[0-9]+\.?[0-9]*} $reffile] {
  		set lc_pos [lindex [split $lc_pos] 0]
                if {![regexp -nocase -start $lc_pos -- {\s*-?[0-9]+\.?[0-9]*\s*} $reffile norm_ref] } {
                	puts stderr "Error reading Initial Lanczos vector norm from reference file"
 	               exit 5
           	}
                if {![regexp -nocase -start $lc_pos -- {\s*-?[0-9]+\.?[0-9]*\s*} $outfile norm_out] } {
                	puts stderr "Error reading Initial Lanczos vector norm from output file"
 	               exit 5
           	}


		if { [ expr { abs( $norm_ref - $norm_out )} ] > 0.00001 } {
                        puts "\[NOT OK\]"
			puts stderr "Discrepancy in initial Lanczos vector norms"
                        exit 5
		}
         }


}

#Compares the lanczos chains
proc ::compare_chains {} {

global reference_file
global output_file

global loop_pos
global chain_pos
global steps
 
        set lc_pos 0 
       
	if [catch {open $reference_file r} fileId] {
                puts stderr "Cannot open $reference_file: $fileId"
                exit 1
        }
        set reffile [read $fileId]
        close $fileId
	if [catch {open $output_file r} fileId] {
                puts stderr "Cannot open $output_file: $fileId"
                exit 1
        }
        set outfile [read $fileId]
        close $fileId
        #Number of Lanczos steps performed
        if {![regexp -nocase -- {Number of Lanczos iterations\s*=\s*([0-9]+)} $reffile match steps_ref] } {
               	puts stderr "Error reading Number of Lanczos steps performed from reference file"
                exit 5
      	}
        if {![regexp -nocase -- {Number of Lanczos iterations\s*=\s*([0-9]+)} $outfile match steps_out] } {
               	puts stderr "Error reading Number of Lanczos steps performed from output file"
                exit 5
      	}
	if { [ expr { $steps_ref != $steps_out } ] } {
                puts "\[NOT OK\]"
		puts stderr "Discrepancy in Lanczos steps performed"
                exit 5
	}
	set steps $steps_ref
	
        #Loop on Lanczos loops
	foreach {loop loop_no} [regexp -all -line -inline -- {Starting Lanczos loop\s+([1-3])} $reffile] {
		set search "Starting Lanczos loop\\s+$loop_no"

 		 if {![regexp -indices -- $search $reffile loop_pos_ref] } {
                                puts stderr "Error reading loop $loop_no in reference"
                                exit 5
                 }
  		set loop_pos_ref [lindex [split $loop_pos_ref] 1]
		 if {![regexp -indices -- $search $outfile loop_pos_out] } {
                                puts stderr "Error reading loop $loop_no in output"
                                exit 5
                 }
  		set loop_pos_out [lindex [split $loop_pos_out] 1]

		for {set i 2} {$i<=$steps+1} {incr i} {
                        set search "beta\\s*\\\(0+$i\\\)=\(.*\)"
		        if {![regexp -start $loop_pos_ref -linestop -- $search $reffile match beta_ref] } {
		                puts stderr "Error reading beta $i in reference"
                		exit 5
		        }	
		        if {![regexp -start $loop_pos_out -linestop -- $search $outfile match beta_out] } {
		                puts stderr "Error reading beta $i in output"
                		exit 5
		        }
			if { [ expr { abs($beta_ref-$beta_out) } ] > 0.001} {
                		puts "\[NOT OK\]"
				puts stderr "Discrepancy in magnitude of p.q step:$i reference:$beta_ref read:$beta_out"
        	        	exit 5
			}
                        set search "gamma\\s*\\\(0+$i\\\)=\(.*\)"
		        if {![regexp -start $loop_pos_ref -linestop -- $search $reffile match gamma_ref] } {
		                puts stderr "Error reading gamma $i in reference"
                		exit 5
		        }	
		        if {![regexp -start $loop_pos_out -linestop -- $search $outfile match gamma_out] } {
		                puts stderr "Error reading gamma $i in output"
                		exit 5
		        }
			if { [ expr { $gamma_ref*$gamma_out } ] < 0.0} {
                		puts "\[NOT OK\]"
				puts stderr "Discrepancy in direction of p.q step:$i reference:$gamma_ref read:$gamma_out"
        	        	exit 5
			}



		}


	}
}
#Looks for bad signs in the execution
proc ::look_for_bads {} {

global reference_file
global output_file

global lc_pos
 
        set lc_pos 0 
       
	if [catch {open $reference_file r} fileId] {
                puts stderr "Cannot open $reference_file: $fileId"
                exit 1
        }
        set reffile [read $fileId]
        close $fileId
	if [catch {open $output_file r} fileId] {
                puts stderr "Cannot open $output_file: $fileId"
                exit 1
        }
        set outfile [read $fileId]
        close $fileId
	if {[regexp -nocase -- {linter: root not converged} $reffile] } {
                puts "\[NOT OK\]"
		puts stderr "DVPSI root not converged in reference file"
                exit 5
        }
	if {![regexp -nocase -- {response charge density does not sum to zero} $reffile] } {
	    if {[regexp -nocase -- {response charge density does not sum to zero} $outfile] } {
               	puts "\[NOT OK\]"
		puts stderr "Response charge leaks in output file whereas reference is normal"
       	        exit 5
            }
       	}

}


main $argc $argv
