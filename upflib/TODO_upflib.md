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
  - ions_base points to uspp_param for variable nsp: is this desirable?
    There are two variables, ntyp (number of types of atoms) and nsp
    (number of types of pseudopotentials), that are not exactly the same 
    but they are in practice. Same for ntypx/nspx. Should we always pass
    ntyp/nsp as argument, or always use the value in the module, or what?
    Now there is a mixture of the two cases.
  - upf_ions now contains just a function n_atom_wfc: move somewhere else?
  - upf_spinorb contains just two variables: merge into uspp?

* upflib restructuring:
  - shall we keep just one src folder ? or structure it a bit more, such as

  upflib/baselib     all basic data structures and io  
  upflib/advlib      advanced initializations (init_us_0,1,2)
  upflib/tools       tools from upftools
  
