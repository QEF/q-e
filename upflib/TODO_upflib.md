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
  - ions_base points to upf_ions for variables na(npsx), nax, nat, ityp, nsp.
    The latter points to uspp for nsp and to upf_params for npsx (this is
    similar to but not exactly the same as variables ntypx in ions_base)
    These variables are however used in a very few places inside upflib:
    init_at_1.f90:  USE upf_ions,     ONLY : ntyp => nsp
    init_us_1.f90:  USE upf_ions,     ONLY : ntyp => nsp, ityp, nat
    init_us_1.f90:  USE upf_ions,     ONLY : ntyp => nsp
    init_us_2_base.f90:  USE upf_ions,     ONLY : nat, ntyp => nsp, ityp
    init_us_2_base_gpu.f90:  USE upf_ions,     ONLY : nat, ntyp => nsp, ityp
    while na, nax are never used
  - first step: pass ntyp, nat, ityp as arguments, remove reference to upf_ions
    DONE
  - second step: nsp becomes an internal variable, not to be exported
  - finally, upf_ions can be deleted and its variables moved back to ions_base
  - upf_spinorb contains just two variables: merge into uspp?

* upflib restructuring:
  - shall we keep just one src folder ? or structure it a bit more, such as

  upflib/baselib     all basic data structures and io  
  upflib/advlib      advanced initializations (init_us_0,1,2)
  upflib/tools       tools from upftools
  
