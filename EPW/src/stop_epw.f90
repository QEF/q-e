  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !
  ! Copyright (C) 2001 PWSCF group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! Modified from stop_ph
  !
  !-----------------------------------------------------------------------
  SUBROUTINE stop_epw()
  !-----------------------------------------------------------------------
  !!
  !! Close all files and synchronize processes before stopping.
  !! Called at the end of the run
  !!
  USE mp,        ONLY : mp_end, mp_barrier
  USE mp_global, ONLY : inter_pool_comm, mp_global_end
  USE io_global, ONLY : stdout
  USE printing,  ONLY : print_clock_epw
  USE epwcom,    ONLY : eliashberg, plselfen, specfun_pl, scattering, iterative_bte, lpolar, lindabs
  USE elph2,     ONLY : adapt_smearing
  USE io_var,    ONLY : epwbib
  USE mp_world,  ONLY : mpime
  USE io_global, ONLY : ionode_id
  !
  IMPLICIT NONE
  !
  CALL print_clock_epw()
  !
  WRITE(stdout, '(a)') "     ==============================================================================="
  WRITE(stdout, '(a)') "     The functionality-dependent EPW.bib file was created with suggested citations. "
  WRITE(stdout, '(a)') "     Please consider citing the papers listed in EPW.bib.                           "
  WRITE(stdout, '(a)') "     ==============================================================================="
  WRITE(stdout, '(a)') "                                                                                    "
  ! 
  IF (mpime == ionode_id) THEN
    !
    OPEN(UNIT = epwbib, FILE = 'EPW.bib')
    ! 
    WRITE(epwbib, '(a)') " % Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino"
    WRITE(epwbib, '(a)') "                                                                                         "
    WRITE(epwbib, '(a)') " % Paper describing the method on which EPW relies                                       "
    WRITE(epwbib, '(a)') " @Article{Giustino2007                                                                   "
    WRITE(epwbib, '(a)') "   Title   = {Electron-phonon interaction using Wannier functions},                      "
    WRITE(epwbib, '(a)') "   Author  = {F. Giustino and M. L. Cohen and S. G. Louie},                              "
    WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                             "  
    WRITE(epwbib, '(a)') "   Year    = {2007},                                                                     "
    WRITE(epwbib, '(a)') "   Volume  = {76},                                                                       " 
    WRITE(epwbib, '(a)') "   Pages   = {165108},                                                                   "
    WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.76.165108}                                                " 
    WRITE(epwbib, '(a)') " }                                                                                       "
    WRITE(epwbib, '(a)') "                                                                                         "
    WRITE(epwbib, '(a)') " % Paper describing the EPW software                                                     "
    WRITE(epwbib, '(a)') " @Article{Ponce2016,                                                                     "
    WRITE(epwbib, '(a)') "   Title   = {EPW: Electron–phonon coupling, transport and superconducting properties &
                                        &using maximally localized Wannier functions},"
    WRITE(epwbib, '(a)') "   Author  = {S. Ponc\'e and E.R. Margine and C. Verdi and F. Giustino},                 "
    WRITE(epwbib, '(a)') "   Journal = {Computer Physics Communications},                                          "
    WRITE(epwbib, '(a)') "   Year    = {2016},                                                                     " 
    WRITE(epwbib, '(a)') "   Volume  = {209},                                                                      " 
    WRITE(epwbib, '(a)') "   Pages   = {116 - 133},                                                                "
    WRITE(epwbib, '(a)') "   Doi     = {https://doi.org/10.1016/j.cpc.2016.07.028}                                 "
    WRITE(epwbib, '(a)') " }                                                                                       "
    ! 
    ! Specific functionalities
    WRITE(epwbib, '(a)') "                                                                                         "
    ! 
    ! Eliashberg superconductivity
    IF (eliashberg) THEN
      WRITE(epwbib, '(a)') "                                                                                         "
      WRITE(epwbib, '(a)') " % Since you used the [eliashberg] input, please consider also citing                    "
      WRITE(epwbib, '(a)') " @Article{Margine2013,                                                                   " 
      WRITE(epwbib, '(a)') "   Title   = {Anisotropic Migdal-Eliashberg theory using Wannier functions},             " 
      WRITE(epwbib, '(a)') "   Author  = {E. R. Margine and F. Giustino},                                            "
      WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                             "
      WRITE(epwbib, '(a)') "   Year    = {2013},                                                                     "
      WRITE(epwbib, '(a)') "   Volume  = {87}                                                                        " 
      WRITE(epwbib, '(a)') "   Pages   = {024505},                                                                   "
      WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.87.024505}                                                "
      WRITE(epwbib, '(a)') " }                                                                                       "
    ENDIF         
    ! 
    ! Polar analytic interpolation
    IF (lpolar) THEN
      WRITE(epwbib, '(a)') "                                                                                         "
      WRITE(epwbib, '(a)') " % Since you used the [lpolar] input, please consider also citing                        "
      WRITE(epwbib, '(a)') " @Article{Verdi2015,                                                                     "
      WRITE(epwbib, '(a)') "   Title   = {Frohlich Electron-Phonon Vertex from First Principles},                    "
      WRITE(epwbib, '(a)') "   Author  = {C. Verdi and F. Giustino},                                                 "
      WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. Lett.},                                                         "
      WRITE(epwbib, '(a)') "   Year    = {2015},                                                                     "
      WRITE(epwbib, '(a)') "   Volume  = {115}                                                                       "
      WRITE(epwbib, '(a)') "   Pages   = {176401},                                                                   "
      WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevLett.115.176401}                                            "
      WRITE(epwbib, '(a)') " }                                                                                       "
    ENDIF
    ! 
    ! Plasmons
    IF (plselfen .OR. specfun_pl) THEN
      WRITE(epwbib, '(a)') "                                                                                         "
      WRITE(epwbib, '(a)') " % Since you used the [plselfen] or [specfun_pl] input, please consider also citing      "
      WRITE(epwbib, '(a)') " @Article{Caruso2018,                                                                    "
      WRITE(epwbib, '(a)') "   Title   = {Electron-plasmon and electron-phonon satellites in the angle-resolved &
                                         &photoelectron spectra of $n$-doped anatase ${\mathrm{TiO}}_{2}$},"
      WRITE(epwbib, '(a)') "   Author  = {F. Caruso and C. Verdi and S. Ponc\'e and F. Giustino},                    "
      WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                             "
      WRITE(epwbib, '(a)') "   Year    = {2018},                                                                     "
      WRITE(epwbib, '(a)') "   Volume  = {97}                                                                        "
      WRITE(epwbib, '(a)') "   Pages   = {165113},                                                                   "
      WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.97.165113}                                                "
      WRITE(epwbib, '(a)') " }                                                                                       "
    ENDIF
    ! 
    ! Transport module
    IF (scattering .AND. .NOT. iterative_bte) THEN
      WRITE(epwbib, '(a)') "                                                                                         "
      WRITE(epwbib, '(a)') " % Since you used the [scattering] input, please consider also citing                    "
      WRITE(epwbib, '(a)') " @Article{Ponce2018,                                                                     "
      WRITE(epwbib, '(a)') "   Title   = {Towards predictive many-body calculations of phonon-limited carrier &
                                         &mobilities in semiconductors},"
      WRITE(epwbib, '(a)') "   Author  = {S. Ponc\'e and E. R. Margine and F. Giustino},                             "
      WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                             "
      WRITE(epwbib, '(a)') "   Year    = {2018},                                                                     "
      WRITE(epwbib, '(a)') "   Volume  = {97}                                                                        "
      WRITE(epwbib, '(a)') "   Pages   = {121201},                                                                   "
      WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.97.121201}                                                "
      WRITE(epwbib, '(a)') " }                                                                                       "
    ENDIF    
    ! 
    IF (iterative_bte) THEN
      WRITE(epwbib, '(a)') "                                                                                         "
      WRITE(epwbib, '(a)') " % Since you used the [iterative_bte] input, please consider also citing                 "
      WRITE(epwbib, '(a)') " @Article{Ponce2018,                                                                     "
      WRITE(epwbib, '(a)') "   Title   = {Towards predictive many-body calculations of phonon-limited carrier &
                                         &mobilities in semiconductors},"
      WRITE(epwbib, '(a)') "   Author  = {S. Ponc\'e and E. R. Margine and F. Giustino},                             "
      WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                             "
      WRITE(epwbib, '(a)') "   Year    = {2018},                                                                     "
      WRITE(epwbib, '(a)') "   Volume  = {97}                                                                        "
      WRITE(epwbib, '(a)') "   Pages   = {121201},                                                                   "
      WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.97.121201}                                                "
      WRITE(epwbib, '(a)') " }                                                                                       "
      WRITE(epwbib, '(a)') "                                                                                         "
      WRITE(epwbib, '(a)') " @Article{Macheda2018,                                                                   "
      WRITE(epwbib, '(a)') "   Title   = {Magnetotransport phenomena in $p$-doped diamond from first principles},    "
      WRITE(epwbib, '(a)') "   Author  = {F. Macheda and N. Bonini},                                                 "
      WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                             "
      WRITE(epwbib, '(a)') "   Year    = {2018},                                                                     "
      WRITE(epwbib, '(a)') "   Volume  = {98}                                                                        "
      WRITE(epwbib, '(a)') "   Pages   = {201201},                                                                   "
      WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.98.201201}                                                "
      WRITE(epwbib, '(a)') " }                                                                                       "      
    ENDIF
    !  
    ! Improvements
    IF (adapt_smearing) THEN
      WRITE(epwbib, '(a)') "                                                                                         "
      WRITE(epwbib, '(a)') " % Since you used the [adapt_smearing] input, please consider also citing                "
      WRITE(epwbib, '(a)') " @Article{Macheda2018,                                                                   "
      WRITE(epwbib, '(a)') "   Title   = {Magnetotransport phenomena in $p$-doped diamond from first principles},    "
      WRITE(epwbib, '(a)') "   Author  = {F. Macheda and N. Bonini},                                                 "
      WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                             "
      WRITE(epwbib, '(a)') "   Year    = {2018},                                                                     "
      WRITE(epwbib, '(a)') "   Volume  = {98}                                                                        "
      WRITE(epwbib, '(a)') "   Pages   = {201201},                                                                   "
      WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.98.201201}                                                "
      WRITE(epwbib, '(a)') " }                                                                                       "
    ENDIF
    !
    ! Indirect optics
    IF (lindabs) THEN
      WRITE(epwbib, '(a)') "                                                                                         "
      WRITE(epwbib, '(a)') " % Since you used the [indabs] input, please consider also citing                        "
      WRITE(epwbib, '(a)') " @Article{Noffsinger2012,                                                                "
      WRITE(epwbib, '(a)') "   Title   = {Phonon-Assisted Optical Absorption in Silicon from First Principles},      "
      WRITE(epwbib, '(a)') "   Author  = {J. Noffsinger and E. Kioupakis and C. G. Van de Walle &
                                         &and S. G. Louie and M. L. Cohen}, "
      WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. Lett.},                                                         "
      WRITE(epwbib, '(a)') "   Year    = {2012},                                                                     "
      WRITE(epwbib, '(a)') "   Volume  = {108}                                                                       "
      WRITE(epwbib, '(a)') "   Pages   = {167402},                                                                   "
      WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevLett.108.167402}                                            "
      WRITE(epwbib, '(a)') " }                                                                                       "
      WRITE(epwbib, '(a)') "                                                                                         "
      WRITE(epwbib, '(a)') " @Article{zhang2022,                                                                     "
      WRITE(epwbib, '(a)') "   Title   = {Ab initio theory of free-carrier absorption in semiconductors},            "
      WRITE(epwbib, '(a)') "   Author  = {X. Zhang and G. Shi and J. A. Leveillee and F. Giustino and E. Kioupakis}, "
      WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                             "
      WRITE(epwbib, '(a)') "   Year    = {2022},                                                                     "
      WRITE(epwbib, '(a)') "   Volume  = {106}                                                                       "
      WRITE(epwbib, '(a)') "   Pages   = {205203},                                                                   "
      WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.106.205203}                                               "
      WRITE(epwbib, '(a)') " }                                                                                       "
    ENDIF
    ! 
    CLOSE(epwbib)
    ! 
  ENDIF
  !
  CALL mp_end(inter_pool_comm)
  CALL mp_global_end()
  !
  STOP
  !
  RETURN
  !-----------------------------------------------------------------------
  END SUBROUTINE stop_epw
  !-----------------------------------------------------------------------
