#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

proc ::main {argc argv} {

global reference_file
global output_file

	get_commandline $argc $argv
	puts -nonewline "Checking $output_file using $reference_file : "
	compare_bands
        compare_total_energy
	puts " \[OK\]"

}
##########################
#A procedure to read commandline
##########################
proc ::get_commandline {argc argv} {

global reference_file
global output_file


        if { $argc < 1} {
                puts stderr "Usage: check_pw <pw output file> <reference file>"
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

#Compares the bands
proc ::compare_bands {} {

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

	if {![regexp -nocase -- {End of self-consistent calculation(.*)!    total energy} $outfile match bands_out_raw] } {
		puts stderr "Error reading bands from output file"
                exit 5
 	   }
        #puts $bands_out_raw
       	if {![regexp -nocase -- {End of self-consistent calculation(.*)!    total energy} $reffile match bands_ref_raw] } {
		puts stderr "Error reading bands from reference file"
                exit 5
  	  }
       #puts $bands_ref_raw

       if {![regexp -nocase -all  -line  -indices  {End of self-consistent calculation} $bands_ref_raw lc_pos] } {
	       set lc_pos 1
	   }
       set lc_pos [lindex $lc_pos end]
       set bands_ref_raw  [regexp -inline -start $lc_pos {.*} $bands_ref_raw]

       #puts $bands_ref_raw

       if {![regexp -nocase -all  -line  -indices  {End of self-consistent calculation} $bands_out_raw lc_pos] } {
	       set lc_pos 1
	   }
       set lc_pos [lindex $lc_pos end]
       set bands_out_raw  [regexp -inline -start $lc_pos {.*} $bands_out_raw]
       
       #puts $bands_out_raw

         foreach {lc_pos} [regexp -all -line -inline -indices {k =} $bands_out_raw] {
  		set lc_pos [lindex [split $lc_pos] 1]
                #puts $lc_pos
                #puts [regexp -inline -start $lc_pos {.*} $bands_out_raw]
                if {![regexp -nocase -start $lc_pos -- {\s*(-?[0-9]+\.?[0-9]*)\s*(-?[0-9]+\.?[0-9]*)\s*(-?[0-9]+\.?[0-9]*)\s*\([^k]*} $bands_out_raw band_kp_out k_out_x k_out_y k_out_z] } {
                	puts stderr "Error reading k-point from output file"
 	               exit 5
           	}
                #puts $band_kp_out 
                if {![regexp -nocase -start $lc_pos -- {\s*(-?[0-9]+\.?[0-9]*)\s*(-?[0-9]+\.?[0-9]*)\s*(-?[0-9]+\.?[0-9]*)\s*\([^k]*} $bands_ref_raw band_kp_ref k_ref_x k_ref_y k_ref_z] } {
                	puts stderr "Error reading k-point from reference file"
 	               exit 5
           	}

		if { [ expr {$k_out_x != $k_ref_x} ] && [ expr {$k_out_y != $k_ref_y} ] && [ expr {$k_out_z != $k_ref_z} ] } {
                        puts "\[NOT OK\]"
			puts stderr "Discrepancy in read K points"
                        exit 5
		}
                if {![regexp -nocase -line -indices -- {bands \(ev\):} $band_kp_out lc_pos] } {
	               puts stderr "Error reading bands for K point $k_out_x $_out_y $k_out_z"
 	               exit 5
           	}

		set lc_pos [lindex [split $lc_pos] 1]
		#puts $lc_pos
                 regexp -start $lc_pos -- {-?[0-9]+\.?[0-9]*} $band_kp_ref band_ref
                 regexp -start $lc_pos -- {-?[0-9]+\.?[0-9]*} $band_kp_out band_out
		 if {[ expr {abs($band_ref - $band_out)} ] > 0.001 } {
                        puts "\[NOT OK\]"
			puts stderr "Discrepancy in band values ref:$band_ref read:$band_out"
                        exit 5
		     }
         }


}
#Compare the total energy
proc ::compare_total_energy {} {

global reference_file
global output_file

 
       
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

	if {![regexp -nocase -all -- {total energy\s*=\s*(-?[0-9]+\.?[0-9]*)\s*Ry} $outfile match energy_out] } {
		puts stderr "Error reading total energy from output file"
                exit 5
 	   }

	if {![regexp -nocase -all -- {total energy\s*=\s*(-?[0-9]+\.?[0-9]*)\s*Ry} $reffile match energy_ref] } {
		puts stderr "Error reading total energy from reference file"
                exit 5
 	   }
	 if {[ expr {abs($energy_ref - $energy_out)} ] > 0.0001 } {
                puts "\[NOT OK\]"
		puts stderr "Discrepancy in total energy: ref:$energy_ref read:$energy_out"
                exit 5
	 }


}
main $argc $argv
