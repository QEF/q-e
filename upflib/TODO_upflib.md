## NOTES ##

* Routines or utilities existing in two copies, one for QE and one for upflib:
  - randy
    in module uspp, file uspp.f90
  - invmat
    simplified version, in module upf_invmat, file upf_invmat.f90
  - capital, lowercase, isnumeric, matches, version_compare
    in module upf_utils, file upf_utils.f90
  - errore
    as upf_error, in file upf_error.f90

* Module variables that have been (partially) duplicated:
   - kinds      => upf_kinds  (only dp)
   - constants  => upf_const
   - nsp in ions_base points to nsp in uspp_param

* TO BE DONE: 

  - vkb (array holding beta functions in reciprocal space) is one of the
    variables of upflib, is computed and deallocated in upflib, but it is 
    allocated outside it. Not very logical. More in general:
    how to deal with objects like vkb that depend upon atomic positions? 
    What belongs to upflib and what doesn't? e.g. module upf_ions, 
    containing only function n_atom_wfc: does it belong here?

  - The usage of the full grid with "upf%mesh" points is useless and even
    dangerous in some cases (may enhance large-r numerical noise).
    The upf% structure contains shorter grids for nonlocal projectors
    (upf%kkbeta and upf%kbeta(:)). Another shorter grid is defined as
    "msh" with a rather arbitrary cutoff at Rmax=10 a.u. and used in
    several cases (e.g. integration of local part). Maybe we should
    move "msh" into the upf structure and allocate/read/write arrays
    using such a shorter grid

  - Simpson integration does not seem to be the best choice. 
    For testing purposes, it would be useful to select at compile time
    (with a preprocessing flag) the fully analytical or numerical form 
    for GTH pseudopotentials and experiments with those.

  - semilocal Vnl and human-readable sections are not read from PP files,
    but they should: converters may make a good usage of that information.

  - It would be VERY useful to read from a given record in fortran (without 
    direct I/O): how can one use fseek, ftell, pos= identifier, stream I/O?

  - set the correct value of nsp in uspp_param when allocate_uspp is called,
    use it ONLY inside upflib, remove link of nsp in ions_base to uspp_param

  - nh(:) is allocated in init_uspp_dims, but maybe it should be allocated
    together with upf(:), when upf is read. Or even better (but annoying
    to do): nh should be part of upf, since it is an atomic quantity?

  - lmaxq should be "the maximum value of L in Q functions", not "... + 1"
    and should be used to dimension arrays where l=0,...,L. 
    The dimension of spherical harmonics (lmaxq+1)^2 might be stored
    in another variable, something like ylmdim, or maxlm

  - Interpolation tables: rationalize names of variables and related routines
```
      CP     PW      better name     contains 	           computed in
    betagx   tab_beta              beta(G) functions	 compute_betagx,
    dbetagx             	   dbeta(G)/dG 		 compute_betagx
    dqradx              	   dQ(G)/dG  		 compute_qradx
             tab_at  tab_atwfc     atomic R_nl(G)	 init_tab_atwfc
    qradx    tab_qrad              Q(G) for  USPP/PAW	 qrad_mod
             tab_rho               atomic rho(G)	 rhoa_mod
             tab_rhc               pseudocore rho(G)	 rhoc_mod
	     tab_vloc              local pseudopotential vloc_mod
```
  - Merge pseudopotential_indexes from CPV/src/pseudopot_sub.f90 with the
    uspp initialization in upflib (init_us_1 etc); merge qvan2b and qvan2
    (requires merge of interpolation tables qrad and qradb)


* upflib restructuring:
  - shall we keep just one src folder ? or structure it a bit more, such as
```
    upflib/baselib     all basic data structures and io  
    upflib/advlib      advanced initializations (init_us_0,1,2)
    upflib/tools       tools from upftools
```

