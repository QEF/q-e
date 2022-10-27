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
   - nsp in ions_base points to nsp in uspp_param

* TO BE DONE: 

  - PP files:
    semilocal Vnl and human-readable sections are not read, but they should:
    converters may make a good usage of that information
    How to read from a given record in fortran (without direct I/O) ?
    fseek, ftell, pos= identifier, stream I/O?

  - set the correct value of nsp in uspp_param when allocate_uspp is called,
    use it ONLY inside upflib, remove link of nsp in ions_base to uspp_param
  - nh(:) is allocated in init_uspp_dims, but maybe it should allocated
    together with upf(:), when upf is read? Or maybe nh should be part of upf?
    It is used in many many places, though!

  - Merge pseudopotential_indexes from CPV/src/pseudopot_sub.f90 with the
    uspp initialization in upflib (init_us_1 etc); merge qvan2b and qvan2
    (requires merge of interpolation tables qrad and qradb)

  - upf_ions now contains just a function n_atom_wfc: move somewhere else?
  - upf_spinorb contains just two variables: merge into uspp? add to it
    the two functions spinor and sph_ind used only for spin-orbit?
  - lmaxq should be "the maximum value of L in Q functions", not "... + 1"
    and should be used to dimension arrays where l=0,...,L. The dimension
    of spherical harmonics (2*lmaxkb+1)^2 is something different and should
    be stored in a different variable (something like ylmdim, or maxlm)

  - Names of interpolation tables and related routines are random:
      CP     PW      new name? 	     contains 	         computed in
    betagx   tab     tab_beta      beta(G) functions	compute_betagx,
    dbetagx             	   dbeta(G)/dG 		compute_betagx
    qradx    qrad    tab_q         Q(G) for  USPP/PAW	init_tab_qrad
    dqradx              	   dQ(G)/dG  		compute_qradx
             tab_at  tab_atwfc     atomic R_nl(G)	init_tab_atwfc
                    (tab_atrho     atomic rho(G)	to be done)
                    (tab_vloc      local potential	to be done)


* upflib restructuring:
  - shall we keep just one src folder ? or structure it a bit more, such as
    upflib/baselib     all basic data structures and io  
    upflib/advlib      advanced initializations (init_us_0,1,2)
    upflib/tools       tools from upftools

