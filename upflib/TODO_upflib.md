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
  - nh should be dynamically allocated together with upf(:), or at the very
    beginning anyway. It seems to me that it belongs to the upf structure,
    but it is used in too many places. Then, npsx can be removed.
  - pre_init, in PW/src/init_run.f90, should be moved to upflib.
    Same for pseudopotential_indexes from CPV/src/pseudopot_sub.f90.
  - upf_ions now contains just a function n_atom_wfc: move somewhere else?
  - upf_spinorb contains just two variables: merge into uspp?

* upflib restructuring:
  - shall we keep just one src folder ? or structure it a bit more, such as

  upflib/baselib     all basic data structures and io  
  upflib/advlib      advanced initializations (init_us_0,1,2)
  upflib/tools       tools from upftools
  
