
* State of the separation of upflib from the rest of QE:
  the following routines or utilities have been copied inside upflib
  - splinelib
    Given that its prevalent usage is for PPs, I have removed the one
    in Modules/
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
    - kinds      => upf_kinds
    - parameters => upf_params
    - constants  => upf_const
  Makefile simplified

* upflib restructuring:
  - shall we keep just one src folder ? or structure it a bit more, such as

  upflib/baselib     all basic data structures and io  
  upflib/advlib      advanced initializations (init_us_0,1,2)
  upflib/tools       tools from upftools
  
  AF: probably structuring is not the highest priority now, though it can
      be useful to organize the material.
  PG: agreed
