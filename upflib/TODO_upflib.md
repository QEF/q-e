
* Routines or utilities existing in two copies, one for QE and one for upflib:
  - randy
    in module uspp, file uspp.f90
  - invmat
    simplified version, in module upf_invmat, file upf_invmat.f90
  - capital, lowercase, isnumeric, matches, version_compare
    in module upf_utils, file upf_utils.f90
  - qe_erf
    as upf_erf, in file upf_erf.f90
  - errore and infomsg
    as upf_error, in file upf_error.f90
  The following modules have been (partially) duplicated:
    - kinds      => upf_kinds  (only dp)
    - parameters => upf_params
    - constants  => upf_const
    - ions_base  => upf_ions   (only ityp, nsp)
    - spinorb    => upf_spinorb
  Makefile simplified

* Note:
  - nsp is both in upf_ions and in uspp_param
  - lmaxx and lqmax are both in upf_params and in spinorb
  - upf_erf is likely obsolete: the erf function should be standard fortran
    and work for all compilers (wasn't true years ago)

* upflib restructuring:
  - shall we keep just one src folder ? or structure it a bit more, such as

  upflib/baselib     all basic data structures and io  
  upflib/advlib      advanced initializations (init_us_0,1,2)
  upflib/tools       tools from upftools
  
