  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  !
  ! Copyright (C) 2002-2013 Quantum ESPRESSO group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! Code adapted from Modules/io_global.f90 - Quantum-ESPRESSO group
  !
  ! SP: I'm missing some that depend on QE like iuwfc, lrwfc .. Should ideally
  !     be included 
  !
  !----------------------------------------------------------------------------
  MODULE io_epw
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  PUBLIC :: lambda_phself, linewidth_phself, linewidth_elself, iospectral, &
            iua2ffil, iudosfil, iufillambda, iuqdos, iufe, iufilker, &
            iufilgap, iospectral_sup, iua2ftrfil, iufilgapFS, iufillambdaFS, &
            iuwanep, iuwane, iunukk, iudvscf
  PUBLIC :: epwdata, iundmedata, iunvmedata, iunksdata, iudyn, iukgmap, iuepb,&
            iurecover, iufilfreq, iufilegnv, iufileph, iufilkqmap, &
            iufilikmap, iueig, iunepmatwp, iunepmatwe, iunkf, iunqf, &
            iufileig, iukmap, crystal, iunifc
  PUBLIC :: iuwinfil, iun_plot, iuukk, iuprojfil !, iummn
  PUBLIC :: iufilsigma, iufilseebeck, iufilkappael, iufilkappa, iufilscatt_rate,&
            iufilFi_all, iufilsigma_all, iufiltau_all
  !
  ! Output of physically relevant quantities (60-100)
  !    
  INTEGER :: lambda_phself   = 60  ! Lambda factor of the phonon self-energy
                                   ! [lambda.phself] 
  INTEGER :: linewidth_phself= 61  ! Imaginary part of the phonon self-energy
                                   ! [linewidth.phself]
  INTEGER :: linewidth_elself= 62  ! Imaginary part of the electron self-energy
                                   ! [linewidth.elself]
  INTEGER :: iospectral      = 63  ! Electronic spectral function [specfun.elself]
  INTEGER :: iospectral_sup  = 64  ! Support data for the spectral function 
                                   ! [specfun_sup.elself]
  INTEGER :: iua2ffil        = 65  ! Eliashberg a2f function [.a2f]
  INTEGER :: iudosfil        = 66  ! Phonon density of state [.phdos]
  INTEGER :: iufillambda     = 67  ! Electron-phonon coupling strength lambda
                                   ! [.lambda_X]
  INTEGER :: iuqdos          = 68  ! Quasiparticle density of states in the
                                   ! superconducting state [.qdos]
  INTEGER :: iufe            = 69  ! Free energy [.fe]
  INTEGER :: iufilker        = 70  ! Eliashberg kernel [.ker]
  INTEGER :: iufilgap        = 71  ! Eliashberg superconducting gap [.gapr]
  INTEGER :: iua2ftrfil      = 72  ! Eliashberg transport a2f function [.a2f_tr]
  INTEGER :: iufilgapFS      = 73  ! Eliashberg superconducting gap on FS with k-points  
  INTEGER :: iufillambdaFS   = 74  ! Electron-phonon coupling strength on FS with k-points
!DBSP : iukgmap was 96. Should be the same as set_kplusq.f90. 
  INTEGER :: iunukk          = 77  ! Unit with rotation matrix U(k) from wannier code
  INTEGER :: iudvscf         = 80  ! Unit for the dvscf_q file
  INTEGER :: iudyn           = 81  ! Unit for the dynamical matrix file
  INTEGER :: iufilkqmap      = 82  ! Map of k+q
  INTEGER :: iukgmap         = 96  ! Map of folding G-vector indexes [.kgmap]
  INTEGER :: iuwanep         = 97  ! Spatial decay of e-p matrix elements in wannier basis 
                                   ! Electrons + phonons [epmat_wanep]
  INTEGER :: iuwane          = 98  ! Spatial decay of matrix elements in Wannier basis    
                                   ! [.epwane]  

  !
  ! Output of quantity for restarting purposes (101-200)
  ! 
  INTEGER :: epwdata         = 101  ! EPW data [epwdata.fmt] 
  INTEGER :: iundmedata      = 102  ! Dipole matrix in wannier basis [dmedata.fmt] 
  INTEGER :: iunvmedata      = 103  ! Velocity matrix in wannier basis [vmedata.fmt]
  INTEGER :: iunksdata       = 104  ! Hamiltonian in wannier basis
  INTEGER :: iuepb           = 105  ! Electron-phonon matrix in Bloch 
                                    ! representation [.epb]
  INTEGER :: iurecover       = 107  ! Dvanqq2 recovery file
  INTEGER :: iufilfreq       = 108  ! Phonon frequency from a previous epw run
                                    ! [.freq]
  INTEGER :: iufilegnv       = 109  ! Eigenvalues from a previous epw run [.egnv]
  INTEGER :: iufileph        = 110  ! Electron-phonon matrix element in the
                                    ! Bloch representation on the fine mesh
                                    ! [.ephmat]
  INTEGER :: iufilikmap      = 112  ! Index of k+(sign)q on the irreducible k-mesh
                                    ! [.ikmap]
!  INTEGER :: iuetf           = 113  ! Interpolated hamiltonian eigenvalues
  INTEGER :: iueig           = 114  ! Temporary eig for interpolation    

  INTEGER :: iunepmatwp      = 115  ! The unit with the e-ph matrix in Wannier-Wannier representation
  INTEGER :: iunepmatwe      = 116  ! The unit with the e-ph matrix in Wannier-Bloch representation
  INTEGER :: iunkf           = 117  ! The unit with the fine k-point mesh in crystal coord.
  INTEGER :: iunqf           = 118  ! The unit with the fine q-point mesh in crystal coord. 
  INTEGER :: iufileig        = 119  ! The unit with eigenenergies [band.eig]
  INTEGER :: iukmap          = 120  ! Unit for the k-point map generation
  INTEGER :: crystal         = 121  ! Unit for crystal data
  INTEGER :: iunifc          = 122  ! Unit for the IFC produced by q2r.x

  !
  ! Output quantites related to Wannier (201-250)
  !  
  INTEGER :: iuwinfil        = 201  ! Wannier projectors and other quantities
! SP : Not used for now but could be in the future. Would require the amn as well.
!  INTEGER :: iummn           = 202  ! Overlap of the cell periodic part of the Bloch 
                                    ! states <u_nmk|u_nk+b>
  INTEGER :: iun_plot        = 203  ! UNK file (needed by Wannier90 for plotting the 
                                    ! real space Wannier functions)
  INTEGER :: iuukk           = 204  ! Final ukk rotation matrix (the big U!)
  INTEGER :: iuprojfil       = 205  ! Unit for projector [.projw90]
  
  !
  ! Output quantites related to transport (251-300)
  
  INTEGER :: iufilsigma      = 251 ! Electrical conductivity
  INTEGER :: iufilseebeck    = 252 ! Seebeck coefficient
  INTEGER :: iufilkappael    = 253 ! Electronic contribution to thermal conductivity
  INTEGER :: iufilkappa      = 254 ! Electronic contribution to thermal conductivity
  INTEGER :: iufilscatt_rate = 255 ! scattering rate
  INTEGER :: iufilFi_all     = 256 ! Fi_all file to retart at X iteration
  INTEGER :: iufilsigma_all  = 257 ! Sigmar_all and Sigmai_all file to retart an interpolation
  INTEGER :: iufiltau_all    = 258 ! inv_tau_all file to retart an interpolation
  !
END MODULE io_epw
