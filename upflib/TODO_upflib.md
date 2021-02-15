
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

* To Be Done:
  - The same code interpolating tab_at appears in atomic_wfc, gen_at_dy,
    plus_u_full. It should be extracted and put into u_fdata
  - The same code appears also in PP/src/atomic_wfc_nc_proj.f90, that should
    be merged with (or into) PW/src/atomic_wfc.f90
  - Very similar code in gen_at_dj should also be moved into uspp_data
  - nsp is both in upf_ions and in uspp_param
  - lmaxx and lqmax are both in upf_params and in spinorb

* upflib restructuring:
  - shall we keep just one src folder ? or structure it a bit more, such as

  upflib/baselib     all basic data structures and io  
  upflib/advlib      advanced initializations (init_us_0,1,2)
  upflib/tools       tools from upftools
  
  AF: probably structuring is not the highest priority now, though it can
      be useful to organize the material.
  PG: agreed
