##NOTES##

* Routines or utilities existing in two copies, one for QE and one for upflib:
  - randy
    in module uspp, file uspp.f90
  - invmat
    simplified version, in module upf_invmat, file upf_invmat.f90
  - capital, lowercase, isnumeric, matches, version_compare
    in module upf_utils, file upf_utils.f90
  - errore and infomsg
    as upf_error, in file upf_error.f90

* Module variables that have been (partially) duplicated:
   - kinds      => upf_kinds  (only dp)
   - constants  => upf_const

* TO BE DONE: 
  - set the correct value of nsp in uspp_param when allocate_uspp is called,
    use it ONLY inside upflib, remove link of nsp in ions_base to uspp_param
  - nh(:) is allocated in init_uspp_dims, but maybe it should allocated
    together with upf(:), when upf is read? Or maybe nh should be part of upf?
    It is used in many many places, though!
  - Merge pseudopotential_indexes from CPV/src/pseudopot_sub.f90 with the
    uspp initialization in upflib (init_us_1 etc)
  - upf_ions now contains just a function n_atom_wfc: move somewhere else?
  - upf_spinorb contains just two variables: merge into uspp?

* upflib restructuring:
  - shall we keep just one src folder ? or structure it a bit more, such as

  upflib/baselib     all basic data structures and io  
  upflib/advlib      advanced initializations (init_us_0,1,2)
  upflib/tools       tools from upftools
  
