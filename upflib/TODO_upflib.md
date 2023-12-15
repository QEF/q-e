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

  - The usage of the full grid with "upf%mesh" points is useless and even
    dangerous in some cases (may enhance large-r numerical noise).
    The upf% structure contains shorter grids for nonlocal projectors
    (upf%kkbeta and upf%kbeta(:)). Another shorter grid is defined as
    "msh" with a rather arbitrary cutoff at Rmax=10 a.u. and used in
    several cases (e.g. integration of local part). Maybe we should
    move "msh" into the upf structure and allocate/read/write arrays
    using such a shorter grid

  - Simpson integration does not seem to be the best choice. One might use 
    GTH analytical PPs converted to numerical form for testing purposes

  - semilocal Vnl and human-readable sections are not read from PP files,
    but they should: converters may make a good usage of that information.

  - It would be VERY useful to read from a given record in fortran (without 
    direct I/O): how can one use fseek, ftell, pos= identifier, stream I/O?

  - set the correct value of nsp in uspp_param when allocate_uspp is called,
    use it ONLY inside upflib, remove link of nsp in ions_base to uspp_param

  - nh(:) is allocated in init_uspp_dims, but maybe it should allocated
    together with upf(:), when upf is read. Or even better (but annoying
    to do): nh should be part of upf, since it is an atomic quantity?

  - Merge pseudopotential_indexes from CPV/src/pseudopot_sub.f90 with the
    uspp initialization in upflib (init_us_1 etc); merge qvan2b and qvan2
    (requires merge of interpolation tables qrad and qradb)

  - upf_ions now contains just a function n_atom_wfc: move somewhere else?
    same for upf_auxtools, that only contains upf_check_atwfc_norm
  - upf_spinorb contains just two variables: merge into uspp? add to it
    the two functions spinor and sph_ind used only for spin-orbit?
  - lmaxq should be "the maximum value of L in Q functions", not "... + 1"
    and should be used to dimension arrays where l=0,...,L. The dimension
    of spherical harmonics (2*lmaxkb+1)^2 is something different and should
    be stored in a different variable (something like ylmdim, or maxlm)

  - Interpolation tables: rationalize names of variables and related routines
```
      CP     PW      better name     contains 	           computed in
    betagx   tab     tab_beta      beta(G) functions	 compute_betagx,
    dbetagx             	   dbeta(G)/dG 		 compute_betagx
    dqradx              	   dQ(G)/dG  		 compute_qradx
             tab_at  tab_atwfc?    atomic R_nl(G)	 init_tab_atwfc
    qradx    tab_qrad              Q(G) for  USPP/PAW	 qrad_mod
             tab_rho               atomic rho(G)	 rhoa_mod
             tab_rhc               pseudocore rho(G)	 rhoc_mod
	     tab_vloc              local pseudopotential vloc_mod
```
  - Interpolation tables: rationalize the structure of the code.
    Move allocation  of interpolation tables into initialization routines,
    setting max |G| as input. Collect interpolation data and related routines
    into a module, one per variable. DONE: for Vloc, rhoc, rhoat
  - Interpolation tables: get rid of CUDA Fortran. Currently interpolation
    tables are computed on CPU, copied to GPU using OpenACC, used via OpenACC.
    Exception: init_us_2 still use CUDA Fortran, so tab_d has DEVICE
    attribute and ylmr2, dylmr2 have a CUDA Fortran version.

* upflib restructuring:
  - shall we keep just one src folder ? or structure it a bit more, such as
```
    upflib/baselib     all basic data structures and io  
    upflib/advlib      advanced initializations (init_us_0,1,2)
    upflib/tools       tools from upftools
```

