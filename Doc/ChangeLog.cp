See file ChangeLog.old for changes after aug. 2004

24-jul-04  few changes in module usage, sort of workaround for ifc 7.1 (CC)

23-jul-04  inputs for string dynamics merged to CP input
           very preliminary sort of manual for FPMD/CP codes (CC)

19-jul-04 - further merging of low level subroutine between FPMD and CP 
          ( cell_move in Module/cell_base.f90 )
          - More input parameters check in Module/read_namelists
          - For CP, restart file is saved in working directory like in FPMD
            and not in output_dir where MD data are saved, this is because
            usually one keep MD trajectories in home dir.
          - added pseudopotential for wannier dynamics example
          - added Wannier postprocessing (from Manu Sharma )
          - fixed a small bug for FPMD and 'diis' electron dynamics (CC)

15-jul-04 added module cp_mass (for car-parrinello electronic mass)
          cpr.f90: lot of staff moved to subroutines (CC)

10-jul-04  Reference to nonexistent subroutines or variables removed
           from newly added code. Tabulators removed (PG).

07-jul-04  New kind of calculation cp-wf added - varius fix for CP with
           wannier functions, Now I'm able to run Sharma examples,
           but the code is still not fully tested.
           Fix in readpp for pseudo different from UPF (CC)

29-jun-04  Added to CPV the string dynamics as implemented by Yosuke
           (not fully tested yet) (CC)

12-Jun-04  deeq and dvan merged with those in uspp.f90
           ipp temporarily moved from module ions_base to cvan (PG)

02-Jun-04  oops, deeq and dvan are complex in uspp.f90, real in CP...
           To be fixed; for the time being, deeq and dvan not merged (PG)

01-Jun-04  deeq, betae merged with deeq, vkb in PW; order of indices
           in deeq and rhovan made compatible with order in PW (PG)

31-May-04  More USPP_related variables moved to Modules/uspp.f90
           Note that nhx => nhm for consistency with other names
           (those ending in x are static dimensioning)
           Parameter ipp no longer needed on input (still used internally):
           PP type assumed following the same logic as in PWscf (PG)

26-may-04  Most variables in module ncprm have been moved to a new module
           uspp_param, shared between PW and CP (in file Modules/uspp.f90)
           Remaining variables in ncprm moved to new module qrl_mod (PG)

28-apr-04  PP cleanup and merge: module "atom", common with PW,
           replaces "atomic_wfc" and part of "ncprm",
           ifpcor => nlcc, rscore => rho_atc as in PW

27-apr-04  PP cleanup and merge: vloc_at is v(r), not r*v(r)

26-apr-04  PP cleanup and merge: rucore => vloc_at

23-apr-04  PP cleanup and merge: mmaxx => ndmx

22-apr-04  Same logic (or lack of it) for DFT used as in PW

21-apr-04  Derivatives of ylm merged, variable cell works again
           (maybe). Indices of gx and gxb reversed, cleanup (PG)
           L=3 sort of implemented (untested). ng0 => gstart (PG)
	
19-apr-04  Next step in USPP harmonization: aainit, spherical
           harmonics merged - derivatives of ylm NOT YET,
           variable cell NOT WORKING (PG)

13-apr-04  First step in USPP harmonization: lx, lqx => lqmax, 
           lix => lmaxx+1, variables in module "uspp.f90", common 
           with PW, used (merge of aainit not yet done)
           invmat3 moved to flib/ and merged with invmat of PW
           Misc: dfloat => dble (PG)

29-mar-04  Various cleanup and code harmonization: 
           date_and_tim moved to flib and used by all code,
           tictac substituted by start_clock/stop_clock
           celldm/alat/at input parameters in FPMD/CP read
           and set as in PW . (CC)

15-mar-04  Almost all neb routines moved to Modules (CC)
           New module check_stop used by all codes
           to check for exit conditions ( maximum time
           or EXIT file ) (CC)

11-mar-04  NEB works for CP as well (CC)

07-mar-04  Cleanup in CPV: no more SSUM and CSUM
           Modules/smallbox.f90 should work again

26-feb-04  Martin Hilgeman, SGI:
   - support for the SGI Altix class of machines, with Intel Itanium2
     processors. These machines run Linux. Please find more information
     on http://www.sgi.com/servers/altix/. I have added an extra
     configure target named 'altix', as well as a '__ALTIX' pre-processor
     macro. The 'altix' target runs either serial or parallel with the
     SGI MPT MPI library, which is optimised for our low-latency,
     high-bandwidth NUMAflex interconnect which allows the use of shared
     memory.

   - modified Makeflags for the 'origin' target and added support for
     SCSL.

   - added support for 1-D, multiple 1-D and 3-D FFT routines from the
     SGI SCSL scientific library. SCSL is the successor of Complib (which
     is currently supported in CP). The two libraries have a different
     calling sequence.and the main advantage is that the same library is
     also supported (with the same calling sequence) on Altix systems. I
     have added a '__SCSL' macro for it and renamed the '__SGI' macro to
     '__COMPLIB' in 'Modules/fft_scalar.f90.

   - I also found a typo in 'CPV/cpr.f90', where all OPEN statements for
     external files had the same unit number. This bug was not in CP90
     v1.3.

   - I had to change the comment character in the scaLAPACK routines,
     because this was causing problems with the Intel Compilers. This
     isn't used anyway.

25-feb-04  merging FPMD/CP added common subroutines (wave_steepest
           wave_verlet ) to advance wave_functions . FPMD friction
           parameter for electrons "gdelt" substituted with "frice"

-------------------------------------------------------------------
Date: 24 Feb 2004    Version: 2.0
-------------------------------------------------------------------

18-feb-04  Initial support for NEB and meta dynamics.
           I do not include NEB dynamics modules in this version,
           because I want to wait for common neb modules, to be built
           as soon as this version has been released (CC)

17-feb-04  outdir added to the path of the output and restart files,
           pseudopotential reading moved out from cprmain subroutine (CC)

16-feb-04  CPV has been "subroutinized" and is ready for NEB like dynamics.
           Note that iosys has been split into two subroutines:
           read_input_file and iosys. 
           The first routine simply calls read_namelists and read_cards
           to read in the stdin, and does not perform any initialization.
           The second (iosys) does not read anythings but copies values
           from input_parameters to local variables. 
           read_input_file is called from the new main program.
           iosys is called from the cprmain subroutine (the old main program).
           This is the scheme used in FPMD.
i          Deallocation statements added to CPV for neb like dynamics. (CC)

09-nov-03  Unit 6 replaced by stdout (module io_global)
           Wavefunctions are in module wavefunction_module

31-jul-03  Major input restructuring, now common with all codes

01-jul-03  Variable-cell is working again (call to sph_bes fixed)

25-jun-03  More merging of common routines (CC)

19-may-03  some cleanup for occupancy and empty state calculation

14-may-03  Bug: namelist &ions must be read in all cases
           Write charge density (if required) only at last step
           Documentation updated

21-apr-03  fft restructuring (Carlo)
           Exch_corr: gradr not deallocated in some cases

12-apr-03  rsg in ortho => rs

27-feb-03  Misc. installation changes

21-feb-03  "error" renamed to "errore", "rnd" to "rndx"
           bug in io_base fixed

11-feb-03  pseudo_dir implemented

10-feb-03  Some cleanup (ibrav, tau written at the end)
           support for intel compiler and linux re-added

------------------------------------------------------------------
First release
------------------------------------------------------------------
 2-feb-03  Ultrasoft UPF bug fixed, more small changes related
           to cpv => cp

 1-feb-03  added check on dimension of pseudopotential arrays
           configure and example cpr.j fixed

10-jan-03  "make tar" or "make dist" produces a tar.gz file with
           a source distribution - Make.sample removed (PG)

05-jan-03  ggen: same ordering of PW and FPMD (using d(:) vector)
           interoperability with FPMD checked also in parallel

04-jan-03  file dimensions.f90 replaced by file parameters.f90
           changes to restart file (CC):
           - io_base.f90 mp.f90 mp_global.f90 mp_wave.f90 updated
           - directory "arch" replaced by "system", file Machine.*
             replaced by Make.*

20-dec-02  Spin-polarized calculation at fixed cell possible again
           Error in core corrections fixed 

16-dec-02  readpseudo.f90: yet another uninitialized variable fixed

11-dec-02  restart.f90: compilation warnings fixed
           readpseudo.f90: upf%tvanp always initialized

04-dec-02  __VARIABLECELL removed everywhere
           Small changes to UPF reading

01-dec-02  New writefile and readfile added
           same restart file layout as FPMD
           Program main alone in the file cpr.f90, all other subroutines
           moved to cprsub.f90 .  Subroutine matinv moved to cplib.f90
           para_mod.f90 compiled even if __PARA is not defined
           startup subroutine now appropriate also in the scalar code

30-nov-02  Module cell changed in cell_module 
           function and types added from FPMD
           mill_l, bi1, bi2, bi3 added
           erroneus usage of twmass corrected

22-nov-02  Minor glitches, documentation updated

21-nov-02  Input updated (final), cpv removed

15-nov-02  cpr.x as fast as cpv.x for fixed-cell calculation
           (useless calls to formf removed) - cpv.f90 is obsolete

14-nov-02  More input changes
           New installation procedure (like FPMD)
           Double underscore prepended to all the CPP macro
           Added modules from FPMD used in the new output format
           bug fix to mp_get and mp_put routines (module "mp")
           Added old "nbeg=-1" option ( suggested by Vittadini)
           Moved calculation of center of mass (suggested by Varadha)

06-nov-02  Compilation error on sp4

04-nov-02  Copyright corrected
           Added possibility to read UPF pseudopotentials

21 oct-02  Compilation problems for cpr on parallel machines,
           gnu license, Make.sample updated, misc. 

08-oct-02  More trouble from unitialized variables (variable-cell,
           intel compiler) fixed

11-sep-02  INPUT documentation updated

31-aug-02  New input layout with the namelists:
           CONTROL, SYSTEM, ELECTRONS, IONS, CELL .
           New ATOMIC_SPECIES card introduced, with the syntax:
             Label(is) pmass(is) psfile(is) ipp(is)
           with:
             character(len=2) label
             real(kind = 8) pmass
             character(len=*) psfile
             integer ipp
           New ATOMIC_POSITIONS card introduced, with the syntax
             label(ia) px(ia) py(ia) pz(ia) .....
           with:
             character(len=2) label 
               this label identify the atom and should match one of those
               present in ATOMIC_SPECIE, and could be optionally follewed
               by an index ( like Cu20 ), to be compliant with the 
               XYZ format.
             real( kind=8 ) px, py, pz  

16-aug-02  flag 'atomic_positions' properly (?) implemented
           fricp was incorrectly read
           more obvious format for units 77 and 78
           Units f77 and f78 are flushed (at least for some compilers)

12-aug-02  Misc. changes for compatibility with other codes:
           iforce for each component, may be specified on input as before
           in spin-polarized case, nbnd = number of spin up states = 
           number of spin down states, not their sum.
           Files are opened and closed during the run in order
           to preserve their content in case of crash;
           I/O-related useless crap removed

08-aug-02  New input - sort of working also in parallel
           PP files are now separated and called by name

06-aug-02  New input - sort of working (not in parallel)

17-jul-02  Start of the Grand Unification
-------------------------------------------------------------------------

24-apr-02  Readvan: check if nang=0 (Yudong)

-------------------------------------------------------------------------
tag:cpr11

 7-mar-02  Added check for consistency between US format and ipp
           (Seungwu)

28-feb-02  Format used in unit 78 increased (Andrea Trave)

27-feb-02  Initialization of Nose' variables not properly done in
           some cases (Xiaofei+Ralph)
           A few formats increased to avoid *** in the output

26-feb-02  More problems in variable-cell + Nose' in the parallel case:
           readpfile, writepfile modified (found by Andrea Trave)
           File format is once again not compatible with previous versions

25-feb-02  Serious (and stupid) bug in init1 if ibrav=0
           and first basis vector had a component along z
           Found by Balazs Hetenyi

22-feb-02  Nose' bug in cpr fixed also when using steepest descent on ions
           Box grid unit vectors are written on output
           (both suggested by Andrea Trave)

-------------------------------------------------------------------------
tag:cpr10

06-feb-02  fix problem with preprocessing on ibm introduced yesterday
           Remaining untyped variables explicitely typed

05-feb-02  added support for pgi compiler on a PC beowulf 
           (Andrea Vittadini): minor changes, documentation update.
           Intel compiler: cpu_time does not work, replaced by etime

01-feb-02  cplib: subroutine rhoset was using uninitialized variables
           in spin-polarized case (found by Yudong).

30-jan-02  cpv: in subroutine ggenb, gxnb(1,*) must be set to zero
           (found by Yudong)

23-jan-02  Default mmx changed to 5000 (500 was too small in most cases)
           (Ralph)

-------------------------------------------------------------------------
tag:cpr9

22-jan-02  More small changes for Compaq parallel machines (Yudong)
           Yet another serious Nose' bug in cpr (found by Ralph)

18-jan-02  Potential bug in Nose' dynamics fixed (some variables were
           not set to zero - the bug appeared with Intel compiler)
           More minor changes (timing routines, Make.sample)

17-jan-02  Added support for intel fortran compiler on linux PC
           (does not work for Nose') and for Compaq parallel machines
           (Thanks to Yudong Wu) (untested)
           Preprocessing simplified, documentation updated,
           minor changes here and there

15-jan-02  fixed bug in readpfile that caused serious trouble to 
           Nose' dynamics when restarting from file in the parallel
           case (xnhpm was not broadcast to all nodes in readpfile)
           Thanks to Xiaofei Wang for remarking the bug

09-nov-01  memory message for origin fixed

-------------------------------------------------------------------------
tag:cpr8

22-oct-01  serious bug in cpr when restarting from previous dynamics
           run fixed
18-oct-01  serious bug in drhov fixed (thanks to Ralph Gebauer): 
           stress was wrong if no ultrasoft atoms were present
27-aug-01  Added memory and file size estimator

-------------------------------------------------------------------------
tag:cpr7

25-aug-01  awful bug in newd (wrong forces in spin-polarized case)
14-aug-01  bug in new init for cpr fixed
           bug in parallel fft for boxes on ibm for n1rx=nr1+1
13-aug-01  more cleaning
           init1 for cpr heavily modified (calls other routines)
10-aug-01  cleaning of unused variables

-------------------------------------------------------------------------
tag:cpr6

09-aug-01  merged file format and related routines (readfile/writefile)
           between cpr and cpv. NOTA BENE:
           files produced by previous versions of the code
           cannot be read by this version.
           Scalar and parallel files still have different formats
           Documentation update
08-aug-01  cpr: major cleanup of nlinit and newnlinit
19-jul-01  First attempt of a parallelization for boxes
           (routines rhov, drhov, newd, set_cc, force_cc)

-------------------------------------------------------------------------
tag:cpr5

17-jul-01  Merge of vofrho in cpv and cpr
           More rhoofr and various other cleaning

-------------------------------------------------------------------------
tag:cpr4

16-jul-01  Variables rhovan, drhovan use compact indices like qgb
           cpr: rhoofr simplified and merged with cpv rhoofr

-------------------------------------------------------------------------
tag:cpr3

14-jul-01  Small box section heavily modified in order to make it
           parallel (parallelization to be finished): 
           - newd works now in real space instead of g-space: slower
             in scalar, in parallel reduces communications to minimum
           - newd, rhov, drhov, set_cc, force_cc:
             common code extracted and put into subroutines
             (box2grid, box2grid2, boxdotgrid)
           - two fft at a time implemented in force_cc
           Timing (hopefully) more readable
           Case ibrav=0 works (again)
           Documentation update

-------------------------------------------------------------------------
tag:cpr2

12-Jul-01  Yet another bug in force_cc for parallel execution
11-Jul-01  Rather serious bug in set_cc fixed
06-Jul-01  Added core corrections to cpv
           Documentation update
21-Jun-01  Documentation update
04-May-01  Out-of-bounds bug in atomic_wfc

-------------------------------------------------------------------------
tag:cpr1

27-Apr-01  First merge of variable-cell calculation, major changes
           There are two executable, "cpr.x" and "cpv.x"
           NOTA BENE: input data for cpv.x changed wrt preceding version

-------------------------------------------------------------------------
tag:cp90_16

19-Apr-01  Yet another bug in boxes (for nr odd) fixed
           printing of elapsed times on origin works (sort of)
           Bug in estimate of S(S+1) with Becke's formula
           in parallel case fixed
           dft is read from file in BHS pseudopotentials as well
           Minor changes to allow more than 64 processors
           Minor corrections here and there

07-Mar-01  Check on pseudopotential sanity added

21-feb-01  Added INPUT.HOWTO

-------------------------------------------------------------------------
tag:cp90_15

09-feb-01  bug in wavefunction write/read for the parallel case fixed
           Make.sample updated for NEC sx-5
           Estimate of S(S+1) added

27-jan-01  latgen modified (once again) so as to yield for ibrav=5
           right-handed axis triplets.

26-jan-01  pseudopotential format converter "pw2us.f90" updated

23-jan-01  latgen modified again to yield more accurate lattices for
           ibrav=5. Also: calculation of shells in ggen and ggenb 
           modified to be more numerically robust.

22-jan-01  latgen modified so as to yield for ibrav=7 and 10 right-handed 
           axis triplets. Boxes for US PPs do not seem to work with the
           original (left-handed) axis triplets. INPUT updated.
           TODO: find what is wrong with the logic of boxes.

18-jan-01  INPUT completed, Make.sample updated for t3e

16-jan-01  checks on nqlc and nang modified so that local PPs work

12-dec-00  nec bug in good_fft_dimension fixed
           added support for nec sx-5 and updated Make.sample
           redefinition of grid in BHS case removed
           added definition of variable f as array in all fft routines

21-nov-00  parallel case for nproc=1 and nr3x=nr3+1 fixed

-------------------------------------------------------------------------
tag:cp90_14

15-nov-00  added routine that reads PPs in Andrea Dal Corso's format

07-nov-00  deeq must be set to zero if non-us pp are to be used!
           Dynamical variables eigr, eigrb, ei1, ei2, ei3 are allocated 
           to the actual maximum number "nas" of atoms of the same kind
           and no longer to fixed parameter nax. 
           Static variables are still dimensioned as (nax,nsx)

06-nov-00  more energic stop in error for parallel case
           Removed hard-coded scratch directory for SP3 case: the scratch 
           directory is read from value of SCRDIR environment variable

25-oct-00  bug in PW91 spin-polarised (finally) found

-------------------------------------------------------------------------
tag:cp90_13

20-oct-00  added support for NEC SX-4

16-oct-00  fixed bug if number of atoms > numbers of states 
           (relevant only for two molecules of H2 or similar cases)

03-oct-00  Make.sample update
           naux increased to 15000 in ibmfft
           ndr=ndw is now allowed (had problem on origin)

26-sep-00  bug in initbox fixed: numerical rounding could lead to rather 
           large error for US pseudopotentials if an atom was very very 
           close to a grid point.
           Limitation on nr1b,nr2b,nr3b even removed.
           Latgen for ibrav=9,10,11,13, fixed
           Minor corrections.

-------------------------------------------------------------------------
tag:cp90_12

09-aug-00  slightly inconsistent calculation of box grid modified;
           exch-corr routines modified so as to be compatible
           with future introduction of cell dynamics. Note that
           the former version of PW91 is still present as "ggapwold".

28-jun-00  COPY is the real, not complex version: 
           needs factor 2 when COPYing complex wavefunctions
           mysterious line "emaec=73" removed

21-jun-00  reduce was missing in ggapw

19-jun-00  PW91 spin-polarised added.  NOTA BENE: since there are some
           differences wrt preceding (spin-unpolarised) results, the
           old routine "ggapwold" has been retained. Use "ggapw" instead
           (in exch_corr) for spin-polarized calculations.
           INPUT file updated

12-jun-00  ortho: test of floating-point error added

-------------------------------------------------------------------------
tag:cp90_11

10-jun-00  very serious bug in sigset for spin-polarized case

08-Jun-00  parallel I/O finally (?) correct (??)

01-Jun-00  parallel I/O better implemented
           some comments added or updated

-------------------------------------------------------------------------
tag:cp90_10

31-May-00  write wavefunctions on one file for parallel execution

29-May-00  write rho on one file for parallel execution

25-May-00  numerical problem in very special cases in LSDA fixed

22-May-00  lim2 in ggapbe was wrong

05-Apr-00  very stupid and serious bug with constraints fixed

--------------------------------------------------------------------------
tag:cp90_9

14-Mar-00  calculation of forces in vofrho is done in separate routines
           direct and reciprocal lattices moved into modules
           more logical names for rhet (=>rhog) and rhoe (=>rhor)
           obvious PBE bug fixed

05-Mar-00  added PBE (written by Michele Lazzeri)

--------------------------------------------------------------------------
tag:cp90_8

07-Feb-00  modules mass, pptype, rcmax_mod moved into ions
           module leng and spin moved into elct
           module control added
           many comments updated, added, displaced

--------------------------------------------------------------------------
tag: cp90_7

06-Feb-00  modules eigrb_mod, irb_mod, teigr removed

05-Feb-00  modules becdr_mod, betae_mod, wbeta_mod, forc removed
           tau0, sfac, deeq, rhovan removed from modules

--------------------------------------------------------------------------
tag: cp90_6

04-Feb-00  added support for absoft, Make.sample updated
           calphi, ortho cleaned

03-Feb-00  added support for origin
           prefor simplified

--------------------------------------------------------------------------
tag: cp90_5

03-Feb-00  added index ish for easier indexing of bec and becdr
           iterative orthonormalization: redundant variables removed

02-Feb-00  indices of becdr rearranged in the same way as for bec

--------------------------------------------------------------------------
tag: cp90_4

02-Feb-00  removed loop (no longer used) for constraints,
              gam, gamold => lambda, olambda
           eigs does no longer produce INF (produces 0.0 ...) on empty states
           major index rearrangements of bec and similar quantities:
              bec(nax,nx,nhx,nsp) => bec(nhsa,nx)

01-Feb-00  formf moved out of the main loop into initialization
           bec removed from modules and called explicitely
           some tictac's moved into subroutines

--------------------------------------------------------------------------
tag: cp90_3

01-Feb-00  Argh! serious bug in formf corrected

31-Jan-00  blypnum removed

29-Jan-00  reversed order of indexes in sfac, rhops, vps
           (should be faster and more logical)

--------------------------------------------------------------------------
tag: cp90_2

29-Jan-00  serious error fixed
           more extensive cleaning:
              phfac and nlpre merged
              strucf does no longer calculate eigr
              read, write, random initialization moved to separate routines
28-Jan-00  some minor cleaning

--------------------------------------------------------------------------
tag: cp90_1 

27-Jan-00  Initial release of f90 code. Main differences wrt f77 version:

- dynamic allocation of memory
- commons replaced by modules
- some general cleanup 
