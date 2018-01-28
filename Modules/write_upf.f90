!
! Copyright (C) 2008-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!--------------------------------------------------------
MODULE  write_upf_module
   !-----------------------------------------------------
   !! this module collects the interfaces for the subroutines 
   !! writing pseudopotential usin XML. 
   !! write_upf_v2 (xf, upf, conf, u_input) writes using UPF v2.1 format 
   !! write_upf_schema(xf, upf, conf, u_input) writes according to qe_pseudo schema
   !! *xf* is a xmlf_f derived type from FoX library: is the xml file descriptor
   !! *upf* is  pseudo_upf derived type from qe Module pseudo_types: stores the whole pseudo data set
   !! *conf* (optional) is a pseudo_conf derived type from qe Module pseudo_types: stores the  info about 
   !!        the configuration used to generate the pseudo. 
   !! *u_input* is the Fortran unit number pointing to the input file of the pseudo creation program
   
   USE write_upf_v2_module, ONLY: write_upf_v2
   USE write_upf_schema_module, ONLY: write_upf_schema
   PUBLIC
END MODULE write_upf_module

