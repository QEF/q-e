  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !
  ! Copyright (C) 2002-2013 Quantum ESPRESSO group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  !
  !----------------------------------------------------------------------------
  MODULE io_var
  !----------------------------------------------------------------------------
  !!
  !! EPW I/O variables
  !! Each integer range has a defined meaning (see below).
  !!
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  PUBLIC :: lambda_phself, linewidth_phself, linewidth_elself, iospectral,    &
            iua2ffil, iudosfil, iufillambda, iuqdos, iufe, iufilker, iuquad,  &
            iufilgap, iospectral_sup, iua2ftrfil, iufilgapFS, iufillambdaFS,  &
            iospectral_cum, iuwanep, iuwane, iunukk, iudvscf, iuqpeig, iures, &
            iuint3paw, iufildos, iufilmat
  PUBLIC :: epwdata, iundmedata, iunvmedata, iunksdata, iudyn, iukgmap, iuepb, &
            iufilfreq, iufilegnv, iufileph, iufilkqmap, iunpattern, iufilmu_q, &
            iufilikmap, iueig, iunepmatwp, iunepmatwe, iunkf, iunqf, iufilFS,  &
            iufileig, iukmap, crystal, iunifc, iunimem, iunepmatwp2, &
            iudwwe, iudgwe, iusthwe, iupmwe, iuxqc, iuqmap, iudmat
  PUBLIC :: iuwinfil, iun_plot, iuprojfil, iudecayH, iudecayP, &
            iudecaydyn, iudecayv, iunnkp, iuamn, iummn, iubvec
  PUBLIC :: iufilsigma, iufilseebeck, iufilkappael, iufilkappa, iufilscatt_rate,   &
            iufilFi_all, iufilsigma_all, iufiltau_all, iuindabs, iuntau, iuntaucb, &
            iufilesigma_all, epwbib, iuindabs_all, iudirabs
  PUBLIC :: iunsparseq, iunsparsek, iunsparsei, iunsparsej, iunsparset, iunselecq, &
            iunsparseqcb, iunsparsekcb, iunsparseicb, iunsparsejcb, iunsparsetcb,  &
            iunrestart, iufilibtev_sup, iunepmat, iunepmatcb, iufilF, iufilmu_nk,  &
            iunsparseq_merge, iunsparsek_merge, iunsparsei_merge,iunsparsej_merge, &
            iunsparset_merge, iunepmatcb_merge, iunsparseqcb_merge,                &
            iunsparseicb_merge, iunsparsejcb_merge, iunsparsetcb_merge,            &
            iunsparsekcb_merge, iunepmat_merge
  PUBLIC :: iunRpscell, iunkgridscell, iunpsirscell, iepfall, ihamil, iMmn, irho,  &
            iUmn, iekanu

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
  INTEGER :: iospectral_cum  = 75  ! Electronic spectral function with the cumulant method
                                   ! [specfun_cum##.elself]
  INTEGER :: iunukk          = 77  ! Unit with rotation matrix U(k) from wannier code
  INTEGER :: iuquad          = 78  ! Unit to read the quadrupole tensor from file
  INTEGER :: iudvscf         = 80  ! Unit for the dvscf_q file
  INTEGER :: iudyn           = 81  ! Unit for the dynamical matrix file
  INTEGER :: iufilkqmap      = 82  ! Map of k+q
  INTEGER :: iufilmat        = 87  ! Matsubara indices
  INTEGER :: iufildos        = 88  ! electronic DOS in Fermi windows [prefix.dos]
  INTEGER :: iukgmap         = 96  ! Map of folding G-vector indexes [.kgmap]
  INTEGER :: iuwanep         = 97  ! Spatial decay of e-p matrix elements in wannier basis
                                   ! Electrons + phonons [epmat_wanep]
  INTEGER :: iuwane          = 98  ! Spatial decay of matrix elements in Wannier basis [.epwane]
  INTEGER :: iuint3paw       = 99  ! Unit for the dvscf_paw_q file
  !
  ! Output of quantity for restarting purposes (101-200)
  ! Note that 100-102 are reserved Cray unit and cannot be used.
  !
  INTEGER :: iunvmedata      = 103  ! Velocity matrix in wannier basis [vmedata.fmt]
  INTEGER :: iunksdata       = 104  ! Hamiltonian in wannier basis
  INTEGER :: iuepb           = 105  ! Electron-phonon matrix in Bloch
                                    ! representation [.epb]
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
  INTEGER :: iunimem         = 123  ! Unit for reading memory information from the system status file
  INTEGER :: epwdata         = 124  ! EPW data [epwdata.fmt]
  INTEGER :: iundmedata      = 125  ! Dipole matrix in wannier basis [dmedata.fmt]
  INTEGER :: iunepmatwp2     = 126  ! Opening the epmatwp file
  INTEGER :: iufilibtev_sup  = 127  ! Files containing velocities for IBTE
  INTEGER :: iunsparseq      = 128  ! Q-mapping for IBTE
  INTEGER :: iunsparsek      = 129  ! K-mapping for IBTE
  INTEGER :: iunsparsei      = 130  ! i band mapping for IBTE
  INTEGER :: iunsparsej      = 131  ! j band mapping for IBTE
  INTEGER :: iunsparset      = 132  ! temperature mapping for IBTE
  INTEGER :: iunsparseqcb    = 133  ! Q-mapping for IBTE of conduction band
  INTEGER :: iunsparsekcb    = 134  ! K-mapping for IBTE for conduction band
  INTEGER :: iunsparseicb    = 135  ! i band mapping for IBTE for conduction band
  INTEGER :: iunsparsejcb    = 136  ! j band mapping for IBTE for conduction band
  INTEGER :: iunsparsetcb    = 137  ! temperature mapping for IBTE for conduction band
  INTEGER :: iunselecq       = 138  ! file containing q-point inside the fsthick windows
  INTEGER :: iunrestart      = 139  ! restart file during writing of IBTE scattering elements
  INTEGER :: iunepmat        = 140  ! Opening the epmatkq files
  INTEGER :: iunepmatcb      = 141  ! Opening the epmatkqcb file
  INTEGER :: iuntau          = 142  ! Opening the tau file
  INTEGER :: iuntaucb        = 143  ! Opening the taucb file
  INTEGER :: iuqpeig         = 144  ! Reading quasi-particle eigenenergies from file
  INTEGER :: iunpattern      = 145  ! Unit for reading the pattern files.
  INTEGER :: iufilFS         = 146  ! Unit for Fermi surface files
  INTEGER :: iudwwe          = 147  ! Debye-Waller matrix in Wannier representation (.dwmatwe)
  INTEGER :: iudgwe          = 148  ! delta g (hopping correction) matrix in Wannier representation (.dgmatwe)
  INTEGER :: iusthwe         = 149  ! Sternheimer matrix in Wannier representation (.sthmatwe)
  INTEGER :: iupmwe          = 150  ! Momentum matrix element (.cpmew)
  INTEGER :: iuxqc           = 151  ! coarse q points (.xqc)
  INTEGER :: iuqmap          = 152  ! File for the symmetry relation of q-points.
                                    ! Needed for unfolding the upper Fan term.
  INTEGER :: iudmat          = 153  ! Symmetry matrix elements
  !
  ! Output quantites related to Wannier (201-250)
  !
  INTEGER :: iuwinfil        = 201  ! Wannier projectors and other quantities
  INTEGER :: iunnkp          = 202  !
! SP : Not used for now but could be in the future. Would require the amn as well.
  INTEGER :: iuamn           = 203  !
  INTEGER :: iummn           = 204  ! Overlap of the cell periodic part of the Bloch
                                    ! states <u_nmk|u_nk+b>
  INTEGER :: iun_plot        = 205  ! UNK file (needed by Wannier90 for plotting the
                                    ! real space Wannier functions)
  INTEGER :: iuprojfil       = 206  ! Unit for projector [.projw90]
  INTEGER :: iudecayH        = 207  ! Hamiltonian decay in real space
  INTEGER :: iudecayP        = 208  ! Dipole decay in real space
  INTEGER :: iudecaydyn      = 209  ! Dynamical matrix decay in real space
  INTEGER :: iudecayv        = 210  ! Velocity matrix decay in real space
  INTEGER :: iubvec          = 211  ! b-vectors and their weight wb
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
  INTEGER :: iufilF          = 259 ! $\partial_E f_{nk}$ in .fmt mode
  INTEGER :: iufilesigma_all = 260 ! Sigmar_all and Sigmai_all file to retart spectral calculation
  INTEGER :: iufilmu_nk      = 261 ! $\mu_{nk}^{\alpha\beta}$ in mobility_nk.fmt file
  INTEGER :: iufilmu_q       = 262 ! $\mu_{\nu q}^{\alpha\beta}$ in mobility_nuq.fmt mode
  INTEGER :: iures           = 263 ! Resistivity in metals using Ziman formula [.res]
  !
  ! Output quantities related to Indirect absorption (301-325)
  INTEGER :: iuindabs        = 301 ! Indirect absorption data
  INTEGER :: iuindabs_all    = 302 ! Indirect absorption read/write
  INTEGER :: iudirabs        = 303 ! Direct absorption data
  !
  ! Miscellaneous (326-350)
  INTEGER :: epwbib          = 326 ! EPW bibliographic file.
  !
  ! Output quantities related to polaron (350-400)
  INTEGER :: iunRpscell      = 351 ! Rp unit cell list within polaron supercell
  INTEGER :: iunkgridscell   = 352 ! Gs k-grid used for transformed supercell
  INTEGER :: iunpsirscell    = 353 ! Polaron wf in real space for transformed supercell
  INTEGER :: iepfall         = 354 ! Polaron
  INTEGER :: ihamil          = 355 ! Polaron
  INTEGER :: iMmn            = 356 ! Polaron
  INTEGER :: irho            = 357 ! Polaron
  INTEGER :: iUmn            = 358 ! Polaron
  INTEGER :: iekanu          = 359 ! Polaron
  !
  ! Merging of files (400-450)
  INTEGER :: iunepmat_merge    = 400
  INTEGER :: iunsparseq_merge  = 401
  INTEGER :: iunsparsek_merge  = 402
  INTEGER :: iunsparsei_merge  = 403
  INTEGER :: iunsparsej_merge  = 404
  INTEGER :: iunsparset_merge  = 405
  !
  INTEGER :: iunepmatcb_merge    = 406
  INTEGER :: iunsparseqcb_merge  = 407
  INTEGER :: iunsparsekcb_merge  = 408
  INTEGER :: iunsparseicb_merge  = 409
  INTEGER :: iunsparsejcb_merge  = 410
  INTEGER :: iunsparsetcb_merge  = 411
  !
  !
  PUBLIC :: io_error
  !
  CONTAINS
    !----------------------------------------------------------------------------
    SUBROUTINE io_error(error_msg)
    !----------------------------------------------------------------------------
    !!
    !! Abort the code and gives an error message
    !!
    !! This routine is adapted from wannier90-3.0.0/src/io.F90
    !!
    !----------------------------------------------------------------------------
    !
    USE io_global, ONLY : stdout
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: error_msg
    !! Error message
    !
    ! Local variables
#ifdef MPI
    CHARACTER(LEN = 50) :: filename
    !! name of the file
    INTEGER :: stderr
    !! Standard error
    INTEGER :: ierr
    !! Error number
    INTEGER :: whoami
    !! Returns node number
    INTEGER :: num_nodes
    !! Number of nodes
    !
    CALL mpi_comm_rank(mpi_comm_world, whoami, ierr)
    CALL mpi_comm_size(mpi_comm_world, num_nodes, ierr)
    !
    IF (num_nodes > 1) THEN
      IF (whoami > 99999) THEN
        WRITE(filename, '(a,a,I0,a)') TRIM(seedname), '.node_', whoami, '.werr'
      ELSE
        WRITE(filename, '(a,a,I5.5,a)') TRIM(seedname), '.node_', whoami, '.werr'
      ENDIF
      stderr = io_file_unit()
      OPEN(UNIT = stderr, FILE = TRIM(filename), FORM = 'formatted', ERR = 105)
      WRITE(stderr, '(1x,a)') TRIM(error_msg)
      WRITE(stderr)
    ENDIF
    !
105 WRITE(*, '(1x,a)') TRIM(error_msg)
106 WRITE(*, '(1x,a,I0,a)') "Error on node ", whoami, ": examine the output/error files for details"
    !
    IF (whoami == 0) THEN
      WRITE(stdout, *) 'Exiting.......'
      WRITE(stdout, '(1x,a)') TRIM(error_msg)
      CLOSE(stdout)
    ENDIF
    !
    CALL MPI_abort(MPI_comm_world, 1, ierr)
    !
#else
    !
    WRITE(stdout, *) 'Exiting.......'
    WRITE(stdout, '(1x,a)') TRIM(error_msg)
    !
    CLOSE(stdout)
    !
    WRITE(*, '(1x,a)') TRIM(error_msg)
    WRITE(*, '(A)') "Error: examine the output/error file for details"
#endif
    !
#ifdef EXIT_FLAG
    CALL EXIT(1)
#else
    STOP
#endif
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE io_error
    !----------------------------------------------------------------------------
    !
  !----------------------------------------------------------------------------
  END MODULE io_var
  !----------------------------------------------------------------------------
