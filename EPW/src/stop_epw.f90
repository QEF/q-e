  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !
  ! Copyright (C) 2001 PWSCF group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE epw_stop
  !----------------------------------------------------------------------
  !!
  !!
  !! This module contains the routines related to k-point or q-point grid
  !! generation as well as selection of k/q points.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE stop_epw()
    !-----------------------------------------------------------------------
    !!
    !! Close all files and synchronize processes before stopping.
    !! Called at the end of the run
    !!
    USE mp,            ONLY : mp_end, mp_barrier
    USE mp_global,     ONLY : inter_pool_comm, mp_global_end
    USE io_global,     ONLY : stdout
    USE printing,      ONLY : print_clock_epw
    USE epwcom,        ONLY : eliashberg, plselfen, specfun_pl, scattering, iterative_bte, lpolar, lindabs, &
                              bfieldx, bfieldy, bfieldz, system_2d, ii_scattering, plrn, loptabs, wfpt
    USE elph2,         ONLY : adapt_smearing, qrpl
    USE io_var,        ONLY : epwbib
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE constants_epw, ONLY : eps40
    !
    IMPLICIT NONE
    !
    CALL print_clock_epw()
    !
    ! S.Tiwari , We write the contents of the bib file in .out file as well
    !
    !
    WRITE(stdout, '(a)') "     % Copyright (C) 2016-2023 EPW-Collaboration                                    "
    WRITE(stdout, '(a)') "                                                                                    "

    WRITE(stdout, '(a)') "     ==============================================================================="
    WRITE(stdout, '(a)') "     Please consider citing the following papers.                                   "
    !
    CALL write_citation()
    !
    WRITE(stdout, '(a)') "                                                                                    "
    WRITE(stdout, '(a)') "     For your convenience, this information is also reported in the                 "
    WRITE(stdout, '(a)') "     functionality-dependent EPW.bib file.                                          "
    WRITE(stdout, '(a)') "     ==============================================================================="
    WRITE(stdout, '(a)') "                                                                                    "
    !
    IF (mpime == ionode_id) THEN
      !
      OPEN(UNIT = epwbib, FILE = 'EPW.bib')
      !
      WRITE(epwbib, '(a)') " % Copyright (C) 2016-2023 EPW Collaboration                                             "
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
      WRITE(epwbib, '(a)') " % Papers describing the EPW software                                                    "
      WRITE(epwbib, '(a)') " @Article{Lee2023,                                                                       "
      WRITE(epwbib, '(a)') "   Title   = {Electron--phonon physics from first principles using the EPW code},        "
      WRITE(epwbib, '(a)') "   Author  = {H. Lee and S. Ponc\'e and K. Bushick and S. Hajinazar and &
                                         &J. Lafuente-Bartolome and J. Leveillee and C. Lian and J. Lihm and &
                                         &F. Macheda and H. Mori and H. Paudyal and W.H. Sio and S. Tiwari and &
                                         &M. Zacharias and X. Zhang and N. Bonini, Nicola and E. Kioupakis and &
                                         &E.R. Margine and F. Giustino},                                             "
      WRITE(epwbib, '(a)') "   Journal = {npj Computational Materials},                                              "
      WRITE(epwbib, '(a)') "   Year    = {2023},                                                                     "
      WRITE(epwbib, '(a)') "   Volume  = {9},                                                                        "
      WRITE(epwbib, '(a)') "   Pages   = {156},                                                                      "
      WRITE(epwbib, '(a)') "   Doi     = {10.1038/s41524-023-01107-3}                                                "
      WRITE(epwbib, '(a)') " }                                                                                       "
      WRITE(epwbib, '(a)') "                                                                                         "
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
        IF(qrpl) THEN
          WRITE(epwbib, '(a)') " % Since you used the [quadrupole] input, please consider also citing                  "
          WRITE(epwbib, '(a)') " @Article{Ponce2021,                                                                   "
          WRITE(epwbib, '(a)') "   Title   = {First-principles predictions of Hall and drift mobilities in &
                                       &semiconductors"
          WRITE(epwbib, '(a)') "   Author  = {S. Ponc\'e and F. Macheda and E.R. Margine and N. Marzari and N. Bonini &
                                       &and F. Giustino  },                                                      "
          WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. Res.},                                                        "
          WRITE(epwbib, '(a)') "   Year    = {2021},                                                                   "
          WRITE(epwbib, '(a)') "   Volume  = {4}                                                                       "
          WRITE(epwbib, '(a)') "   Pages   = {143022},                                                                 "
          WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevResearch.3.043022}                                        "
          WRITE(epwbib, '(a)') " }                                                                                     "
          WRITE(epwbib, '(a)') "                                                                                       "
          WRITE(epwbib, '(a)') " @Article{Ponce2023,                                                                   "
          WRITE(epwbib, '(a)') "   Title   = {Long-range electrostatic contribution to electron-phonon couplings and &
                                       &mobilities of two-dimensional and bulk materials},                       "
          WRITE(epwbib, '(a)') "   Author  = {S. Ponc\'e and M. Royo and M. Stengel and N. Marzari and M. Gibertini},  "
          WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                           "
          WRITE(epwbib, '(a)') "   Year    = {2023},                                                                   "
          WRITE(epwbib, '(a)') "   Volume  = {107}                                                                     "
          WRITE(epwbib, '(a)') "   Pages   = {155424},                                                                 "
          WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.107.155424}                                             "
          WRITE(epwbib, '(a)') " }                                                                                     "
        ENDIF
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
      IF (scattering .OR. iterative_bte) THEN
        WRITE(epwbib, '(a)') "                                                                                         "
        WRITE(epwbib, '(a)') " % Since you used the [scattering/iterative_bte] input, please consider also citing      "
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
      IF (ABS(bfieldx)+ABS(bfieldy)+ABS(bfieldz) > eps40) THEN
        WRITE(epwbib, '(a)') "                                                                                         "
        WRITE(epwbib, '(a)') " % Since you used the [bfield] input, please consider also citing                        "
        WRITE(epwbib, '(a)') " @Article{Macheda2018,                                                                   "
        WRITE(epwbib, '(a)') "   Title   = {Magnetotransport phenomena in $p$-doped diamond from first principles},    "
        WRITE(epwbib, '(a)') "   Author  = {F. Macheda and N. Bonini},                                                 "
        WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                             "
        WRITE(epwbib, '(a)') "   Year    = {2018},                                                                     "
        WRITE(epwbib, '(a)') "   Volume  = {98}                                                                        "
        WRITE(epwbib, '(a)') "   Pages   = {201201},                                                                   "
        WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.98.201201}                                                "
        WRITE(epwbib, '(a)') " }                                                                                       "
        WRITE(epwbib, '(a)') "                                                                                         "
        WRITE(epwbib, '(a)') " @Article{Ponce2021,                                                                     "
        WRITE(epwbib, '(a)') "   Title   = {First-principles predictions of Hall and drift mobilities in &
                                           &semiconductors},                                                           "
        WRITE(epwbib, '(a)') "   Author  = {S. Ponc\'e and F. Macheda and E. R. Margine and N. Marzari and &
                                           & N. Bonini and F. Giustino},                                               "
        WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. Res.},                                                          "
        WRITE(epwbib, '(a)') "   Year    = {2021},                                                                     "
        WRITE(epwbib, '(a)') "   Volume  = {3}                                                                         "
        WRITE(epwbib, '(a)') "   Pages   = {043022},                                                                   "
        WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevResearch.3.043022}                                          "
        WRITE(epwbib, '(a)') " }                                                                                       "
      ENDIF
      !
      ! 2D materials
      IF (system_2d /= 'no') THEN
        IF((system_2d .EQ. 'dipole_sp') .OR. (system_2d .EQ. 'quadrupole')) THEN
          WRITE(epwbib, '(a)') "                                                                                       "
          WRITE(epwbib, '(a)') " % Since you used the [system_2d==dipole_sp] input, please consider also citing        "
          WRITE(epwbib, '(a)') " @Article{Ponce2023,                                                                   "
          WRITE(epwbib, '(a)') "   Title   = {Long-range electrostatic contribution to electron-phonon couplings and &
                                       &mobilities of two-dimensional and bulk materials},                       "
          WRITE(epwbib, '(a)') "   Author  = {S. Ponc\'e and M. Royo and M. Stengel and N. Marzari and M. Gibertini},  "
          WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                           "
          WRITE(epwbib, '(a)') "   Year    = {2023},                                                                   "
          WRITE(epwbib, '(a)') "   Volume  = {107}                                                                     "
          WRITE(epwbib, '(a)') "   Pages   = {155424},                                                                 "
          WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.107.155424}                                             "
          WRITE(epwbib, '(a)') " }                                                                                     "
          WRITE(epwbib, '(a)') "                                                                                       "
          WRITE(epwbib, '(a)') "                                                                                       "
          WRITE(epwbib, '(a)') " @Article{Ponce2023,                                                                   "
          WRITE(epwbib, '(a)') "   Title   = {Accurate Prediction of Hall Mobilities in Two-Dimensional Materials &
                                             &through Gauge-Covariant Quadrupolar Contributions},                      "
          WRITE(epwbib, '(a)') "   Author  = {S. Ponc\'e and M. Royo and M. Gibertini and N. Marzari and M. Stengel},  "
          WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. Lett.},                                                       "
          WRITE(epwbib, '(a)') "   Year    = {2023},                                                                   "
          WRITE(epwbib, '(a)') "   Volume  = {130}                                                                     "
          WRITE(epwbib, '(a)') "   Pages   = {166301},                                                                 "
          WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevLett.130.166301}                                          "
          WRITE(epwbib, '(a)') " }                                                                                     "
          WRITE(epwbib, '(a)') "                                                                                       "
        ENDIF
        IF(system_2d .EQ. 'dipole_sh') THEN
          WRITE(epwbib, '(a)') "                                                                                       "
          WRITE(epwbib, '(a)') " % Since you used the [system_2d==dipole_sh] input, please consider also citing        "
          WRITE(epwbib, '(a)') " @Article{Sio2022,                                                                     "
          WRITE(epwbib, '(a)') "   Title   = {Unified ab initio description of Frohlich electron-phonon interactions&
                                               &in two-dimensional andthree-dimensional materials },"
          WRITE(epwbib, '(a)') "   Author  = {W.H. Sio and F. Giustino},                                               "
          WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                           "
          WRITE(epwbib, '(a)') "   Year    = {2022},                                                                   "
          WRITE(epwbib, '(a)') "   Volume  = {105}                                                                     "
          WRITE(epwbib, '(a)') "   Pages   = {115414},                                                                 "
          WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.105.115414}                                             "
          WRITE(epwbib, '(a)') " }                                                                                     "
          WRITE(epwbib, '(a)') "                                                                                       "
        ENDIF
        IF(system_2d .EQ. 'gaussian') THEN
          WRITE(epwbib, '(a)') "                                                                                       "
          WRITE(epwbib, '(a)') " % Since you used the [system_2d==gaussian] input, please consider also citing         "
          WRITE(epwbib, '(a)') " @Article{Sohier2016,                                                                  "
          WRITE(epwbib, '(a)') "   Title   = {Two-dimensional Fr\'ohlich interaction in transition-metal &
                                              &dichalcogenide monolayers: Theoretical modeling and &
                                              &first-principles calculations                                           "
          WRITE(epwbib, '(a)') "   Author  = {T. Sohier and M. Calandra and F. Mauri},                                 "
          WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                           "
          WRITE(epwbib, '(a)') "   Year    = {2016},                                                                   "
          WRITE(epwbib, '(a)') "   Volume  = {94}                                                                      "
          WRITE(epwbib, '(a)') "   Pages   = {085415},                                                                 "
          WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.94.085415}                                              "
          WRITE(epwbib, '(a)') " }                                                                                     "
          WRITE(epwbib, '(a)') "                                                                                       "
        ENDIF
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
        WRITE(epwbib, '(a)') " % Since you used the [lindabs] input, please consider also citing                       "
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
      ! Impurity scattering
      IF (ii_scattering) THEN
        WRITE(epwbib, '(a)') "                                                                                         "
        WRITE(epwbib, '(a)') " % Since you used the [ii_scattering] input, please consider also citing                 "
        WRITE(epwbib, '(a)') " @Article{Leveillee2023,                                                                 "
        WRITE(epwbib, '(a)') "   Title   = {Ab initio calculation of carrier mobility in semiconductors including &
                                      ionized-impurity scattering},    "
        WRITE(epwbib, '(a)') "   Author  = {J. Leveillee, and X. Zhang and E. Kioupakis and F. Giustino},              "
        WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                             "
        WRITE(epwbib, '(a)') "   Year    = {2023},                                                                     "
        WRITE(epwbib, '(a)') "   Volume  = {107}                                                                       "
        WRITE(epwbib, '(a)') "   Pages   = {125207},                                                                   "
        WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.107.125207}                                               "
        WRITE(epwbib, '(a)') " }                                                                                       "
      ENDIF
      !
      !Polaron
      IF (plrn) THEN
        WRITE(epwbib, '(a)') "                                                                                         "
        WRITE(epwbib, '(a)') " % Since you used the [plrn] input, please consider also citing                          "
        WRITE(epwbib, '(a)') " @Article{Sio2019,                                                                       "
        WRITE(epwbib, '(a)') "   Title   = {Ab initio theory of polarons: Formalism and applications},                 "
        WRITE(epwbib, '(a)') "   Author  = {W.H. Sio, C.Verdi, S. Ponc\'e, F. Giustino},                               "
        WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. B},                                                             "
        WRITE(epwbib, '(a)') "   Year    = {2019},                                                                     "
        WRITE(epwbib, '(a)') "   Volume  = {99}                                                                        "
        WRITE(epwbib, '(a)') "   Pages   = {235139},                                                                   "
        WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevB.99.235139}                                                "
        WRITE(epwbib, '(a)') " }                                                                                       "
        WRITE(epwbib, '(a)') " @Article{Sio2019,                                                                       "
        WRITE(epwbib, '(a)') "   Title   = {Polarons from First Principles, without Supercells},                       "
        WRITE(epwbib, '(a)') "   Author  = {W.H. Sio, C.Verdi, S. Ponc\'e, F. Giustino},                               "
        WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. Lett.},                                                         "
        WRITE(epwbib, '(a)') "   Year    = {2019},                                                                     "
        WRITE(epwbib, '(a)') "   Volume  = {122}                                                                       "
        WRITE(epwbib, '(a)') "   Pages   = {246403},                                                                   "
        WRITE(epwbib, '(a)') "   Doi     = {10.1103/PhysRevLett.122.246403}                                            "
        WRITE(epwbib, '(a)') " }                                                                                       "
      ENDIF
      !
      !Full optical absorption (not yet published)
      IF (loptabs) THEN
        WRITE(epwbib, '(a)') "                                                                                         "
        WRITE(epwbib, '(a)') " % Since you used the [loptabs] input, please consider also citing                       "
        WRITE(epwbib, '(a)') " @Article{Tiwari2023,                                                                    "
        WRITE(epwbib, '(a)') "   Title   = {Quasidegenerate manybody perturbation theory of direct and phonon-assisted&
                                            optical absorption and luminescence},                                      "
        WRITE(epwbib, '(a)') "   Author  = {S.Tiwari and E. Kioupakis and J. Menendez and F. Giustino},                "
        WRITE(epwbib, '(a)') "   Journal = {},                                                                         "
        WRITE(epwbib, '(a)') "   Year    = {2023},                                                                     "
        WRITE(epwbib, '(a)') "   Volume  = {},                                                                          "
        WRITE(epwbib, '(a)') "   Pages   = {},                                                                         "
        WRITE(epwbib, '(a)') "   Doi     = {}                                                                          "
        WRITE(epwbib, '(a)') " }                                                                                       "
        WRITE(epwbib, '(a)') " @Article{Tiwari2023,                                                                    "
        WRITE(epwbib, '(a)') "   Title   = {Ab-initio quasi-degenerate many-body perturbation theory of optical &
                                              absorption including both direct and phonon-assisted processes},         "
        WRITE(epwbib, '(a)') "   Author  = {S.Tiwari and F. Giustino},                                                 "
        WRITE(epwbib, '(a)') "   Journal = {},                                                                         "
        WRITE(epwbib, '(a)') "   Year    = {2023},                                                                     "
        WRITE(epwbib, '(a)') "   Volume  = {},                                                                         "
        WRITE(epwbib, '(a)') "   Pages   = {},                                                                         "
        WRITE(epwbib, '(a)') "   Doi     = {}                                                                          "
        WRITE(epwbib, '(a)') " }                                                                                       "
      ENDIF
      !
      ! WFPT
      IF (wfpt) THEN
        WRITE(epwbib, '(a)') "                                                                                         "
        WRITE(epwbib, '(a)') " % Since you used the [wfpt] input, please consider also citing                          "
        WRITE(epwbib, '(a)') " @Article{Lihm2021,                                                                      "
        WRITE(epwbib, '(a)') "   Title   = {Wannier Function Perturbation Theory: Localized Representation and &
                                            &Interpolation of Wave Function Perturbation},                 "
        WRITE(epwbib, '(a)') "   Author  = {Jae-Mo Lihm and Cheol-Hwan Park},                                          "
        WRITE(epwbib, '(a)') "   Journal = {Phys. Rev. X},                                                             "
        WRITE(epwbib, '(a)') "   Year    = {2021},                                                                     "
        WRITE(epwbib, '(a)') "   Volume  = {11}                                                                        "
        WRITE(epwbib, '(a)') "   Pages   = {041053},                                                                   "
        WRITE(epwbib, '(a)') "   Doi     = {10.1103/physrevx.11.041053}                                                "
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
    SUBROUTINE write_citation()
    !-----------------------------------------------------------------------
    !!
    !! Close all files and synchronize processes before stopping.
    !! Called at the end of the run.
    !!
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : eliashberg, plselfen, specfun_pl, scattering, iterative_bte, lpolar, lindabs, &
                              bfieldx, bfieldy, bfieldz, system_2d, loptabs, plrn, ii_scattering, wfpt
    USE elph2,         ONLY : adapt_smearing, qrpl
    USE io_global,     ONLY : ionode_id
    USE constants_epw, ONLY : eps40
    !
    IMPLICIT NONE
    !
    WRITE(stdout, '(a)') "                                                                                                 "
    WRITE(stdout, '(a)') "     % Paper describing the method on which EPW relies                                           "
    WRITE(stdout, '(a)') "       F. Giustino and M. L. Cohen and S. G. Louie, Phys. Rev. B 76, 165108 (2007)               "
    WRITE(stdout, '(a)') "                                                                                                 "
    WRITE(stdout, '(a)') "     % Papers describing the EPW software                                                        "
    WRITE(stdout, '(a)') "       H. Lee et al., npj Comput. Mater. 9, 156 (2023)                                           "
    WRITE(stdout, '(a)') "       S. Ponc\'e, E.R. Margine, C. Verdi and F. Giustino, Comput. Phys. Commun. 209, 116 (2016) "
    WRITE(stdout, '(a)') "       J. Noffsinger et al., Comput. Phys. Commun. 181, 2140 (2010)                              "
    !
    ! Specific functionalities
    WRITE(stdout, '(a)') "                                                                                                 "
    !
    ! Eliashberg superconductivity
    IF (eliashberg) THEN
      WRITE(stdout, '(a)') "                                                                                             "
      WRITE(stdout, '(a)') "     % Since you used the [eliashberg] input, please consider also citing                    "
      WRITE(stdout, '(a)') "       E. R. Margine and F. Giustino, Phys. Rev. B 87, 024505 (2013)                         "
    ENDIF
    !
    ! Polar analytic interpolation
    IF (lpolar) THEN
      WRITE(stdout, '(a)') "                                                                                             "
      WRITE(stdout, '(a)') "     % Since you used the [lpolar] input, please consider also citing                        "
      WRITE(stdout, '(a)') "       C. Verdi and F. Giustino, Phys. Rev. Lett. 115, 176401 (2015)                         "
      IF(qrpl) THEN
        WRITE(stdout, '(a)') "     % Since you used the quadrupole.fmt file, please consider also citing                 "
        WRITE(stdout, '(a)') "       S. Ponc\'e et al, Phys. Rev. Res. 4, 143022 (2021)                                  "
        WRITE(stdout, '(a)') "       S. Ponc\'e et al. Phys. Rev. B 107, 155424 (2023)                                   "
      ENDIF
    ENDIF
    !
    ! Plasmons
    IF (plselfen .OR. specfun_pl) THEN
      WRITE(stdout, '(a)') "                                                                                             "
      WRITE(stdout, '(a)') "     % Since you used the [plselfen] or [specfun_pl] input, please consider also citing      "
      WRITE(stdout, '(a)') "       F. Caruso et al, Phys. Rev. B 97, 165113 (2018)                                       "
    ENDIF
    !
    ! Transport module
    IF (scattering .OR. iterative_bte) THEN
      WRITE(stdout, '(a)') "                                                                                             "
      WRITE(stdout, '(a)') "     % Since you used the [scattering/iterative_bte] input, please consider also citing      "
      WRITE(stdout, '(a)') "       S. Ponc\'e, E. R. Margine and F. Giustino, Phys. Rev. B 97, 121201 (2018)             "
      WRITE(stdout, '(a)') "       F. Macheda and N. Bonini, Phys. Rev. B 98, 201201 (2018)                              "
    ENDIF
    !
    IF (ABS(bfieldx)+ABS(bfieldy)+ABS(bfieldz) > eps40) THEN
      WRITE(stdout, '(a)') "                                                                                             "
      WRITE(stdout, '(a)') "     % Since you used the [bfield] input, please consider also citing                        "
      WRITE(stdout, '(a)') "       F. Macheda and N. Bonini, Phys. Rev. B 98, 201201 (2018)                              "
      WRITE(stdout, '(a)') "       S. Ponc\'e et al, Phys. Rev. Res. 4, 143022 (2021)                                    "
    ENDIF
    !
    ! 2D materials
    IF (system_2d /= 'no') THEN
      IF((system_2d .EQ. 'dipole_sp') .OR. (system_2d .EQ. 'quadrupole'))THEN
        WRITE(stdout, '(a)') "                                                                                           "
        WRITE(stdout, '(a)') "     % Since you used the [system_2d=dipole_sp or quadrupole] input, please consider also &
                                   &citing                                                                               "
        WRITE(stdout, '(a)') "       S. Ponc\'e et al, Phys. Rev. B 107, 155424 (2023)                                   "
        WRITE(stdout, '(a)') "       S. Ponc\'e et al, Phys. Rev. Lett. 130, 166301 (2023)                               "
      ENDIF
      IF(system_2d .EQ. 'dipole_sh') THEN
        WRITE(stdout, '(a)') "                                                                                           "
        WRITE(stdout, '(a)') "     % Since you used the [system_2d==dipole_sh] input, please consider also citing        "
        WRITE(stdout, '(a)') "       W.H. Sio and F. Giustino, Phys. Rev. B 105, 115414 (2022)                           "
      ENDIF
      IF(system_2d .EQ. 'gaussian') THEN
        WRITE(stdout, '(a)') "                                                                                           "
        WRITE(stdout, '(a)') "     % Since you used the [system_2d==gaussian] input, please consider also citing         "
        WRITE(stdout, '(a)') "       T. Sohier and M. Calandra and F. Mauri, Phys. Rev. B 94, 085415 (2016)              "
      ENDIF
    ENDIF
    !
    ! Improvements
    IF (adapt_smearing) THEN
      WRITE(stdout, '(a)') "                                                                                             "
      WRITE(stdout, '(a)') "     % Since you used the [adapt_smearing] input, please consider also citing                "
      WRITE(stdout, '(a)') "       F. Macheda and N. Bonini, Phys. Rev. B 98, 201201 (2018)                              "
    ENDIF
    !
    ! Indirect optics
    IF (lindabs) THEN
      WRITE(stdout, '(a)') "                                                                                             "
      WRITE(stdout, '(a)') "     % Since you used the [lindabs] input, please consider also citing                       "
      WRITE(stdout, '(a)') "       J. Noffsinger et al, Phys. Rev. Lett. 108, 167402 (2012)                              "
      WRITE(stdout, '(a)') "       X. Zhang et al, Phys. Rev. B 106, 205203 (2022)                                       "
    ENDIF
    !
    ! Impurity scattering
    IF (ii_scattering) THEN
      WRITE(stdout, '(a)') "                                                                                             "
      WRITE(stdout, '(a)') "     % Since you used the [ii_scattering] input, please consider also citing                 "
      WRITE(stdout, '(a)') "       J. Leveillee et al, Phys. Rev. B 107, 125207 (2023)                                   "
    ENDIF
    !
    ! Polaron
    IF (plrn) THEN
      WRITE(stdout, '(a)') "                                                                                             "
      WRITE(stdout, '(a)') "     % Since you used the [plrn] input, please consider also citing                          "
      WRITE(stdout, '(a)') "       W. H. Sio et al, Phys. Rev. B 99, 235139 (2019)                                       "
      WRITE(stdout, '(a)') "       W. H. Sio et al, Phys. Rev. Lett. 122, 246403 (2019)                                  "
    ENDIF
    !
    ! Full optical absorption
    IF (loptabs) THEN
      WRITE(stdout, '(a)') "                                                                                             "
      WRITE(stdout, '(a)') "     % Since you used the [loptabs] input, please consider also citing                       "
      WRITE(stdout, '(a)') "       S. Tiwari et al, unpublished (2023)                                                   "
      WRITE(stdout, '(a)') "       S. Tiwari and F. Giustino, unpublished (2023)                                         "
    ENDIF
    !
    ! Wannier Function perturbation theory
    IF (wfpt) THEN
      WRITE(stdout, '(a)') "                                                                                             "
      WRITE(stdout, '(a)') "     % Since you used the [wfpt] input, please consider also citing                          "
      WRITE(stdout, '(a)') "       J.-M. Lihm and C.-H. Park, PRX 11, 041053 (2021)                                      "
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE  write_citation
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE ephf_deallocate(nrr_k, nrr_q, nrr_g, dims, dims2, epmatwef,      &
                               cufkk, cufkq, uf, w2, cfac, cfacq, rdotk, rdotk2,&
                               irvec_r, irvec_k, irvec_q, irvec_g, ndegen_k,    &
                               ndegen_q, ndegen_g, wslen_k, wslen_q, wslen_g,   &
                               etf_all, epmatf, eimpmatf)
    !-----------------------------------------------------------------------
    !!
    !! Deallocate variables at the end of fine grid interpolation
    !!
    USE kinds,        ONLY : DP
    USE modes,        ONLY : nmodes
    USE klist_epw,    ONLY : isk_dummy
    USE ions_base,    ONLY : nat
    USE epwcom,       ONLY : iterative_bte, ephwrite, mp_mesh_k, etf_mem, vme, &
                             epmatkqread, cumulant, eliashberg, assume_metal,  &
                             lindabs, carrier, ii_g, nbndsub, scattering
    USE elph2,        ONLY : map_rebal, map_rebal_inv, vmef, cvmew, cdmew,     &
                             dmef, epmatwp, chw, chw_ks, rdw, epsi, zstar,     &
                             wf, etf, etf_ks, nbndfst, nktotf, eps_rpa,        &
                             epstf_therm, bztoibz, s_bztoibz, gtemp, et_ks,    &
                             dos, Qmat, ef0_fca, partion, qtf2_therm,          &
                             epsilon2_abs, epsilon2_abs_lorenz,                &
                             epsilon2_abs_all, epsilon2_abs_lorenz_all
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_k
    !! Number of WS points for electrons
    INTEGER, INTENT(in) :: nrr_q
    !! Number of WS points for phonons
    INTEGER, INTENT(in) :: nrr_g
    !! Number of WS points for electron-phonons
    INTEGER, INTENT(in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(in) :: dims2
    !! Dims is either nat if use_ws or 1 if not
    INTEGER, ALLOCATABLE, INTENT(inout) :: irvec_k(:, :)
    !! Integer components of the ir-th Wigner-Seitz grid point in the basis
    !! of the lattice vectors for electrons
    INTEGER, ALLOCATABLE, INTENT(inout) :: irvec_q(:, :)
    !! INTEGER components of the ir-th Wigner-Seitz grid point for phonons
    INTEGER, ALLOCATABLE, INTENT(inout) :: irvec_g(:, :)
    !! INTEGER components of the ir-th Wigner-Seitz grid point for electron-phonon
    INTEGER, ALLOCATABLE, INTENT(inout) :: ndegen_k(:, :, :)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, ALLOCATABLE, INTENT(inout) :: ndegen_q(:, :, :)
    !! Wigner-Seitz weights for the phonon grid that depend on atomic positions $R + \tau(nb) - \tau(na)$
    INTEGER, ALLOCATABLE, INTENT(inout) :: ndegen_g(:, :, :)
    !! Wigner-Seitz weights for the electron-phonon grid that depend on
    !! atomic positions $R - \tau(na)$
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: w2(:)
    !! Interpolated phonon frequency
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: rdotk(:)
    !! $r\cdot k$
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: rdotk2(:)
    !! $r\cdot k$
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: irvec_r(:, :)
    !! Wigner-Size supercell vectors, store in real instead of integer
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: wslen_k(:)
    !! real-space length for electrons, in units of alat
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: wslen_q(:)
    !! real-space length for phonons, in units of alat
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: wslen_g(:)
    !! real-space length for electron-phonons, in units of alat
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: etf_all(:, :)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: epmatwef(:, :, :, :)
    !! e-p matrix  in el wannier - fine Bloch phonon grid
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: epmatf(:, :, :)
    !! e-p matrix  in smooth Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: eimpmatf(:, :)
    !! carrier-ionized impurity matrix in smooth Bloch basis
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: cufkk(:, :)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: cufkq(:, :)
    !! the same, for points k+q
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: uf(:, :)
    !! Rotation matrix for phonons
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: cfac(:, :, :)
    !! Used to store $e^{2\pi r \cdot k}$ exponential
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: cfacq(:, :, :)
    !! Used to store $e^{2\pi r \cdot k+q}$ exponential
    !
    ! Local variables
    INTEGER :: ierr
    !! Error status
    !
    IF ((iterative_bte .OR. ephwrite) .AND. mp_mesh_k .AND. etf_mem < 3) THEN
      DEALLOCATE(map_rebal, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating map_rebal', 1)
      DEALLOCATE(map_rebal_inv, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating map_rebal_inv', 1)
    ENDIF
    IF (vme == 'wannier') THEN
      DEALLOCATE(vmef, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating vmef', 1)
      DEALLOCATE(cvmew, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating cvmew', 1)
    ELSE
      DEALLOCATE(cdmew, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating cdmew', 1)
      DEALLOCATE(dmef, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating dmef', 1)
    ENDIF
    IF (etf_mem == 0) THEN
      DEALLOCATE(epmatwp, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating epmatwp', 1)
    ENDIF
    !
    DEALLOCATE(chw, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating chw', 1)
    DEALLOCATE(chw_ks, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating chw_ks', 1)
    DEALLOCATE(rdw, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating rdw', 1)
    DEALLOCATE(epsi, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating epsi', 1)
    DEALLOCATE(zstar, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating zstar', 1)
    !
    DEALLOCATE(epmatwef, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating epmatwef', 1)
    IF (.NOT. epmatkqread) THEN
      DEALLOCATE(wf, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating wf', 1)
    ENDIF
    DEALLOCATE(etf, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating etf', 1)
    DEALLOCATE(etf_ks, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating etf_ks', 1)
    DEALLOCATE(epmatf, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating epmatf', 1)
    DEALLOCATE(eimpmatf, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating eimpmatf', 1)
    DEALLOCATE(cufkk, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating cufkk', 1)
    DEALLOCATE(cufkq, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating cufkq', 1)
    DEALLOCATE(uf, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating uf', 1)
    DEALLOCATE(isk_dummy, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating isk_dummy', 1)
    DEALLOCATE(eps_rpa, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating eps_rpa', 1)
    DEALLOCATE(epstf_therm, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating epstf_therm', 1)
    DEALLOCATE(w2, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating w2', 1)
    DEALLOCATE(cfac, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating cfac', 1)
    DEALLOCATE(cfacq, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating cfacq', 1)
    DEALLOCATE(rdotk, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating rdotk', 1)
    DEALLOCATE(rdotk2, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating rdotk2', 1)
    DEALLOCATE(irvec_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating irvec_r', 1)
    DEALLOCATE(irvec_k, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating irvec_k', 1)
    DEALLOCATE(irvec_q, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating irvec_q', 1)
    DEALLOCATE(irvec_g, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating irvec_g', 1)
    DEALLOCATE(ndegen_k, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating ndegen_k', 1)
    DEALLOCATE(ndegen_q, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating ndegen_q', 1)
    DEALLOCATE(ndegen_g, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating ndegen_g', 1)
    DEALLOCATE(wslen_k, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating wslen_k', 1)
    DEALLOCATE(wslen_q, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating wslen_q', 1)
    DEALLOCATE(wslen_g, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating wslen_g', 1)
    DEALLOCATE(etf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating etf_all', 1)
    IF (mp_mesh_k .OR. (etf_mem == 3)) THEN
      DEALLOCATE(bztoibz, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating bztoibz', 1)
      DEALLOCATE(s_bztoibz, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating s_bztoibz', 1)
    ENDIF
    ! Deallocate temperature when no cumulant or supercond
    IF ((.NOT. cumulant) .AND. (.NOT. eliashberg)) THEN
      DEALLOCATE(gtemp, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating gtemp', 1)
    ENDIF
    DEALLOCATE(et_ks, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating et_ks', 1)
    IF (assume_metal) THEN
      DEALLOCATE(dos, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating dos', 1)
    ENDIF
    DEALLOCATE(Qmat, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating Qmat', 1)
    !
    IF (lindabs .AND. (.NOT. scattering)) THEN
      !
      IF (carrier) THEN
        DEALLOCATE(ef0_fca, STAT = ierr)
        IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating ef0_fca', 1)
        IF (ii_g) THEN
          DEALLOCATE(partion, STAT = ierr)
          IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating partion', 1)
          DEALLOCATE(qtf2_therm, STAT = ierr)
          IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating qtf2_therm', 1)
        ENDIF
      ENDIF
      DEALLOCATE(epsilon2_abs, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating epsilon2_abs', 1)
      DEALLOCATE(epsilon2_abs_lorenz, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating epsilon2_abs_lorenz', 1)
      DEALLOCATE(epsilon2_abs_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating epsilon2_abs_all', 1)
      DEALLOCATE(epsilon2_abs_lorenz_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephf_deallocate', 'Error deallocating epsilon2_abs_lorenz_all', 1)
    ENDIF

    !-----------------------------------------------------------------------
    END SUBROUTINE  ephf_deallocate
    !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  END MODULE epw_stop
  !-----------------------------------------------------------------------
