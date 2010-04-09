#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

proc ::main {argc argv} {

global beta_gamma_z_file
global pw_output_file

	
	get_commandline $argc $argv
         
        extract_beta_gamma_z



}
##########################
#A procedure to read commandline
##########################
proc ::get_commandline {argc argv} {

global pw_output_file
global beta_gamma_z_file


        if { $argc < 1} {
                puts stderr "Usage: tddfpt_regen_bgz.tcl <input file> \[parameters\]"
                puts stderr "Parameters:"
		puts stderr "none at the moment "
		exit 1
        }

	if { $argc > 0 } {
		set state first
		foreach arg $argv {
   			switch -- $state {
      		first {
			set state flag
			#check if the infile exists
			set pw_output_file $arg
		        if { [file exists $pw_output_file] } {
                		puts "using output file $pw_output_file"
		        } else {
                		puts stderr "File $pw_output_file does not exits"
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
set beta_gamma_z_file "$pw_output_file.beta_gamma_z"

}
##########################
#Decides where to put the output
##########################

proc ::extract_beta_gamma_z {} {

global beta_gamma_z_file
global pw_output_file
global n_step
global n_loop
global lc_pos
global lanc_norm
global z1
#read the output file into memory
	if [catch {open $pw_output_file r} fileId] {
                puts stderr "Cannot open $pw_output_file: $fileId"
                exit 1
        }
        set tmpfile [read $fileId]
        close $fileId
#find the number of iterations
       if {![regexp -linestop -- {Number of Lanczos iterations =\s*([0-9]+)} $tmpfile match n_step]} {
		puts "Error reading Number of Lanczos steps performed, exiting"
		exit 2
 	}
        puts "The number of Lanczos steps for each loop : $n_step"  
#find number of loops
      if {![regexp -all -- {.*tarting Lanczos loop\s+([1-3])} $tmpfile match n_loop ]} {
		puts "Error reading Number of Lanczos loops performed, exiting"
		exit 2
 	}
        puts "Number of loops : $n_loop"
        incr n_loop 1
#loop over polarization directions
        set lc_pos 0
        for {set current_loop 1} {$current_loop<$n_loop} {incr current_loop} {
#find the initial norm of vectors
	        if {![regexp -start $lc_pos -- {Norm of initial Lanczos vectors=\s+([0-9]+\.[0-9]+)} $tmpfile match lanc_norm ]} {
        	        puts "Error reading Norm for loop $current_loop"
	                exit 2
	        }
                puts "Norm for loop $current_loop : $lanc_norm"
               
#open the file and put the preliminaries
       		 if [catch {open $beta_gamma_z_file.$current_loop w} fileId] {
        	         puts stderr "Cannot open $beta_gamma_z_file.$current_loop: $fileId"
                	 exit 1
	         }
                 puts $fileId [format "%12i" $n_step]
                 puts $fileId [format "%19.14f" $lanc_norm]
#locate the current loop
                 if {![regexp -indices -start $lc_pos -linestop -- {tarting Lanczos loop\s+[1-3]} $tmpfile match]} {
	        	puts "Error locating loop $current_loop"
		        exit 2
 		 }
                 set lc_pos [lindex [split $match] 1]
#start reading the indices
                 for {set current_step 1} {$current_step<=$n_step} {incr current_step} {
	                  if {![regexp -indices -start $lc_pos -linestop -- "Lanczos iteration:\\s*$current_step" $tmpfile match]} {
	        		puts "Error locating step $current_step in loop $current_loop"
			        exit 2
 			 }
	                 set lc_pos [lindex [split $match] 1]
                         #puts "Step no $current_step at $lc_pos"
#read beta 
		      	 if {![regexp -start $lc_pos -linestop -- {beta.*=(-?[0-9]*\.[0-9]+E[-+]?[0-9]+)} $tmpfile match beta]} {
                        	puts "Error locating beta in step $current_step, loop $current_loop"
	                        exit 2
        	          }
                          puts $fileId [format "%19.14f" $beta]

#read gamma	
		      	 if {![regexp -start $lc_pos -linestop -- {gamma.*=(-?[0-9]*\.[0-9]+E[-+]?[0-9]+)} $tmpfile match gamma]} {
                        	puts "Error locating beta in step $current_step, loop $current_loop"
	                        exit 2
        	          }
                          puts $fileId [format "%19.14f" $gamma]
#read zeta
                         for {set  i 1} {$i<$n_loop} {incr i} {
                         	 if {![regexp -start $lc_pos -linestop -- "z1=\\s+$i\\s+(-?\[0-9\]*\\.\[0-9\]+E\[-+\]?\[0-9\]+)\\s+(-?\[0-9\]*\\.\[0-9\]+E\[-+\]?\[0-9\]+)" $tmpfile match z1i z1r]} {
                                 	puts "Error locating z1_$i in step $current_step, loop $current_loop"
	                                exit 2
        	                 }
                                 set z1i [format "%21.15E" $z1i]
                                 regsub -- {(E[+-])} $z1i {\10} z1i
                                 set z1r [format "%21.15E" $z1r]
                                 regsub -- {(E[+-])} $z1r {\10} z1r
                                 puts $fileId " ($z1i,$z1r)"

                         } 

                
                 }
               

                 

                 close $fileId
        }
       
 
       
       
        

}
main $argc $argv
