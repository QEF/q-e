! Copyright (C) 2016-2023 Marios Zacharias, Feliciano Giustino 
!                                                                            
! This file is distributed under the terms of the GNU General Public         
! License. See the file `LICENSE' in the root directory of the               
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
Module ifconstants
  ! This code generates ZG displacements
  !
  ! All variables read from file that need dynamical allocation
  !
  USE kinds, ONLY: DP
  REAL(DP), ALLOCATABLE :: frc(:, :, :, :, :, :, :), tau_blk(:, :),  zeu(:, :, :), &
               m_loc(:, :)
  ! frc : interatomic force constants in REAL space
  ! tau_blk : atomic positions for the original cell
  ! zeu : effective charges for the original cell
  ! m_loc: the magnetic moments of each atom
  INTEGER, ALLOCATABLE  :: ityp_blk(:)
  ! ityp_blk : atomic types for each atom of the original cell
  !
  CHARACTER(LEN=3), ALLOCATABLE :: atm(:)
end Module ifconstants
!
!-------------------------------------------------------------------------
PROGRAM ZG
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !! authors: Marios Zacharias, Feliciano Giustino 
  !! acknowledgement: Sabyasachi Tiwari for help packaging this release
  !! version: v1.1
  !! license: GNU
  !
  !  This program generates the ZG displacements for a list of
  !  q vectors that are comensurate to the supercell size to be employed for 
  !  the special displacement method. The first part of the code is obtained by 
  !  modifying the matdyn.f90 subroutine of QE. The program starts from the 
  !  interatomic force constants generated from the DFPT phonon code through
  !  the companion program q2r.
  !
  !  ZG_displacements generate a supercell of the original cell with the atoms 
  !  displaced via Eq. (2) of https://doi.org/10.1103/PhysRevResearch.2.013357. 
  !  Data required for the ZG displacements is read from the force constant 
  !  file "*.fc" and.
  !
  !  Input cards for ZG.in: namelist &input (first six as for matdyn.f90)
  !     flfrc     file produced by q2r containing force constants (needed)
  !               It is the same as in the input of q2r.x (+ the .xml extension
  !               IF the dynamical matrices produced by ph.x were in xml
  !               format). No default value: must be specified.
  !      asr      (character) indicates the type of Acoustic Sum Rule imposed
  !               - 'no': no Acoustic Sum Rules imposed (default)
  !               - 'simple':  previous implementation of the asr used
  !                  (3 translational asr imposed by correction of
  !                  the diagonal elements of the force constants matrix)
  !               - 'crystal': 3 translational asr imposed by optimized
  !                  correction of the force constants (projection).
  !               - 'one-dim': 3 translational asr + 1 rotational asr
  !                  imposed by optimized correction of the force constants
  !                  (the rotation axis is the direction of periodicity;
  !                   it will work only IF this axis considered is one of
  !                   the cartesian axis).
  !               - 'zero-dim': 3 translational asr + 3 rotational asr
  !                  imposed by optimized correction of the force constants
  !               Note that in certain cases, not all the rotational asr
  !               can be applied (e.g. IF there are only 2 atoms in a
  !               molecule or IF all the atoms are aligned, etc.).
  !               In these cases the supplementary asr are cancelled
  !               during the orthonormalization procedure (see below).
  !               using tetrahedra and a uniform q-point grid (see below)
  !               NB: may not work properly in noncubic materials
  !               IF .FALSE. calculate phonon bands from the list of q-points
  !               supplied in input (default)
  !     amass     masses of atoms in the supercell (a.m.u.), one per atom type
  !               (default: use masses read from file flfrc)
  !    "q_in_band_form" and "q_in_cryst_coord" meaningful if "q_external" 
  !    (see below) is set to .true.
  !     q_in_band_form IF .TRUE. the q points are given in band form:
  !               Only the first and last point of one or more lines 
  !               are given. See below. (default: .FALSE.).
  !     q_in_cryst_coord IF .TRUE. input q points are in crystalline 
  !              coordinates (default: .FALSE.)
  !     fd         (logical) if .t. the ifc come from the finite displacement calculation
  !     na_ifc     (logical) add non analitic contributions to the interatomic force 
  !                constants if finite displacement method is used (as in Wang et al.
  !                Phys. Rev. B 85, 224303 (2012)) 
  !     loto_2d  set to .true. to activate two-dimensional treatment of LO-TO
  !              siplitting for q= 0.  (default: .false.)
  !
  !  IF (q_in_band_form) THEN
  !     nq     ! number of q points
  !     (q(i, n), i = 1, 3), nptq   nptq is the number of points between this point
  !                            and the next. These points are automatically
  !                            generated. the q points are given in Cartesian
  !                            coordinates, 2pi/a units (a=lattice parameters)
  !  ELSE  :
  !     nq         number of q-points
  !     (q(i, n), i = 1, 3)    nq q-points in cartesian coordinates, 2pi/a units
  !  If q = 0, the direction qhat (q=>0) for the non-analytic part
  !  is extracted from the sequence of q-points as follows:
  !     qhat = q(n) - q(n- 1)  or   qhat = q(n) - q(n+ 1)
  !  depending on which one is available and nonzero.
  !  For low-symmetry crystals, specify twice q = 0 in the list
  !  IF you want to have q = 0 results for two different directions
  !  ----------------------------------------------------------------------------------
  ! 
  !  Input cards to control "ZG_configuration" subroutine:
  !
  !     "ZG_conf"            : Logical flag that enables the creation of the ZG-displacement. 
  !                            (default .true.) 
  !     "T"                  : Real number indicating the temperature at which the calculations will be performed. 
  !                            "T" essentially defines the amplitude of the normal coordinates. 
  !                            (default 0.00)
  !     "dim1","dim2","dim3" : Integers corresponding to the dimensionality of the supercell.
  !                            (default 0,0,0)
  !     "atm_zg(1), etc.."   : String describing the element of each atomic species
  !                            (default "Element")
  !     "flscf"              : String for the name of the scf input file used to calculate the phonons. The 
  !                            code will read information for preparing the input file for the supercell 
  !                            scf calculation. If left empty the code will not generate the input file for
  !                            the supercell scf calculation.
  !                            (default ' ')
  !     "qhat_in"            : Vector with three real entries for specifying the direction qhat 
  !                            for the non-analytic part when dim1=dim2=dim3=1.
  !                            Use for example "qhat_in(1) = 0.1, qhat_in(2) =0.0, qhat_in(3) = 0.0"
  !                            to account for LO-TO splitting from the direction [1 0 0]. 
  !                            (default 0.1,0.1,0.1)
  !     "synch"              : Logical flag that enables the synchronization of the modes. 
  !                            (default .false.)
  !     "niters"             : Integer for the number of iterations the algorithm needs to 
  !                            go through for finding the optimum configuration. The algorithm 
  !                            generates a set of "+,-,+,-" signs and its possible permutations, 
  !                            trying to minimize the error coming from the coupling of modes with 
  !                            the same q-wavevector but at different branch. For a finite supercell
  !                            size the order of using the "+,-,+,-" set and its permutations is  
  !                            important giving different results. Therefore the algorithm checks 
  !                            the combination that brings the error lower than a threshold.
  !                            (default 15000)
  !     "compute_error"      : Logical flag: if set to .true. allows the code to find the optimal ZG configuration 
  !                            by minimizing the error based on the "threshold" flag (see below). Set it
  !                            to .false. if speed up is required. Setting it to .false. is useful when 
  !                            (i) large supercell sizes are considered for which the error is minimized by 
  !                            the first set of signs, and (ii) only single phonon displacements are of interest (see below) 
  !                            (default .true.)
  !     "error_thresh"        : Real number indicating the error at which the algorithm stops while it's 
  !                            looking for possible combinations of signs. Once this limit is reached 
  !                            the ZG-displacement is constructed. The threshold is usually chosen 
  !                            to be less than 30% of the diagonal terms, i.e. those terms that contribute 
  !                            to the calculation of temperature-dependent properties. 
  !                            (default 0.05)
  !     "incl_qA"            : Logical flag, to decide whether to include phonon modes in set A or not. 
  !                            (default .true.)
  !     "ASDM"               : Logical flag, that enables the iterative procedure for evaluating anharmonic IFCs.
  !                            (default .false.)
  !     "apply_fd"           : Logical flag, that generates input files for the calculation of the IFCs. Finite displacements 
  !                            of amplitude "fd_displ" along each direction and atom of the ZG configuration are applied 
  !                            (default .true.)
  !     "fd_displ"           : Real number representing the displacement applied along each direction and atom of the 
  !                            ZG configuration.
  !                            (default 0.01058) 
  !                            Meaningful only if ASDM and ZG_conf are true.
  !     "read_fd_forces"     : Logical flag that enables reading of the forces from folder "fd_forces" after
  !                            performing finite differences on the polymorphous or the ZG-configuration.
  !                            Meaningful only if ASDM is true and IFCs have been calculated.
  !                            (default .false.)
  !     "iter_idx"           : Integer indicating the index of the current iteration. If "poly" flag is true then the code sets 
  !                            iter_idx = 0. If not specified the code stops. Based on this index, the code performs 
  !                            iterative mixing of IFCs between iter_idx0 and iter_idx.
  !                            (default 1) 
  !     "iter_idx0"          : Integer indicating the index of the initial IFCs matrix for the iterative mixing. Mixing is
  !                            performed between IFCs of iterations "iter_idx0" and  "iter_idx". If "poly" flag is true then 
  !                            the code sets iter_idx0 = 0. 
  !                            (default 0) 
  !     "mixing"             : Logical flag that enables mixing of IFCs between "iter_idx0" and "iter_idx". 
  !                            (default .false.) 
  !     "update_equil"       : Logical flag that enables the update of the atomic positions at each iteration. This flag 
  !                            serves for minimizing the free energy with respect to atomic positions based on the 
  !                            Newton-Raphson method. Meaningful only if atomic coordinates are specified by general 
  !                            Wyckoff positions or symetries are disabled.
  !                            (default .false.) 
  !     "incl_epsil"         : Logical flag that allows the code to print the information of the dielectric constant and 
  !                            Born effective charge tensors in the ".fc" file printed at each iteration. This flag enables 
  !                            essentially the inclusion of long-range dipole-dipole interactions in the computation of 
  !                            temperature-dependent anharmonic phonons. Meaningful if "read_fd_forces" = .true.
  !                            (default .false.)
  !     "poly"               : Logical flag that enables the starting point of ASDM using the polymorphous structure.
  !                            The code searches for ZG-relax.out whose final atomic coordinates represent the ground state 
  !                            polymorphous network. If "ZG-relax.out" exists the code will apply finite displacements and
  !                            generate "ZG-scf_poly_iter_00" input files. If "ZG-relax.out" does not exist, then the code 
  !                            will generate "ZG-relax.in" containing an initial ZG configuration for performing
  !                            optimization of the nuclei coordinates and thus obtain "ZG-relax.out". 
  !     "poly_fd_forces"     : Logical flag that enables reading of the forces from files "ZG-scf_poly_iter_00" in 
  !                            folder "fd_forces" after performing finite differences on the polymorphous configuration 
  !                            After completion of this starting step set "poly_fd_forces = .false." and "poly = .false.". 
  !                            (default .false.)
  !                            (starting point of ASDM). Meaningful only if ASDM is true and IFCs have been calculated.
  !     "single_ph_displ"    : Logical flag that allows to displace the nuclei along single phonon modes. 
  !                            Use output configurations to compute electron-phonon matrix elements with a direct 
  !                            supercell calculation. Set the displacement to the zero point by "T = 0". 
  !                            This finite displacement should carry precisely the effect of diagonal elements of [g(q)+g(-q)].
  !                            Output files: "single_phonon-displacements.dat" and 
  !                            "single_phonon-velocities.dat".
  !                            (default .false.)
  !     "q_external"         : Logical flag that allows the use of a q-point list specified by the user in the input file. 
  !                            If .false. the q-point list is specified by the supercell dimensions dim1, dim2, and dim3. 
  !                            If .false. any q-point list after the input flags is ignored.
  !                            If .true. the q-point list must be provided by the user (or see "qlist_AB.txt").
  !                            IF ph_unfold = .true. then q_external = .true. automatically and the q-path is provided as 
  !                            in a standard phonon dispersion calculation. 
  !                            (default .false.)
  !     "qlist_AB.txt"       : This file contains the external q-list in crystal coordinates as in the "ZG_444.in" example,
  !                            after the input flags. It corresponds to the q-points commensurate to the supercell size. 
  !                            Only one of the q-point time-reversal partners is kept for the construction of the 
  !                            ZG-displacement. The calculations, for the moment, assume systems with time-reversal symmetry. 
  !                            For the generation of the "qlist_AB.txt" set the q-gird in file 
  !                            "example/silicon/input/qlist.in" and run "../../../src/create_qlist.x < qlist.in > qlist.out".
  !                            One can modify the "create_qlist.f90" to generate a different path for consecutive q-points.
  !                            Paste the output of "qlist_AB.txt" to "ZG.in" after namelist &input. Set the flag 
  !                            q_external = .true. for the code to read the list.  
  ! 
  !     "ph_unfold"          : Logical flag to activate phonon unfolding procedure. (default: .false.). To perform phonon 
  !                            unfolding ZG_conf must be set to .false.. If ph_unfold = .true. then q_external = .true. 
  !
  !     "flfrq"              : Output file for frequencies to printed with unfolding weights (default: 'frequencies.dat')
  !
  !     "flweights"          : Output file for unfolding weights to printed with frequncies (default: 'unfold_weights.dat')
  !
  !    "ng1","ng2","ng3"     : Integers corresponding to the (h k l) indices of the reciprocal lattice vector g.  
  !                            Increase their values to check convergence. Default is a good starting point. 
  !                            (default 10,10,10)
  !
  USE kinds,            ONLY : DP
  USE mp,               ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_world,         ONLY : world_comm 
  USE mp_global,        ONLY : mp_startup, mp_global_end, inter_pool_comm
  USE environment,      ONLY : environment_start, environment_end
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE io_dyn_mat,       ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                               read_ifc_param, read_ifc
  USE cell_base,        ONLY : at, bg, celldm
  USE constants,        ONLY : RY_TO_THZ, RY_TO_CMM1, amu_ry
  USE symm_base,        ONLY : set_sym
  USE rap_point_group,  ONLY : code_group
  USE bz_form,          ONLY : transform_label_coord
  USE parser,           ONLY : read_line
  USE rigid,            ONLY : dyndiag, nonanal, nonanal_ifc

  USE ifconstants,      ONLY : frc, atm, zeu, tau_blk, ityp_blk, m_loc
  !
  IMPLICIT NONE
  !
  !
  ! variables *_blk refer to the original cell, other variables
  ! to the (super)cell (which may coincide with the original cell)
  !
  INTEGER, PARAMETER  :: ntypx= 10, nrwsx=200
  REAL(DP), PARAMETER :: eps = 1.0d-6
  INTEGER             :: nr1, nr2, nr3, nsc, ibrav
  CHARACTER(LEN=256)  :: flfrc, flscf, flfrq, flweights
  CHARACTER(LEN= 10)  :: asr
  LOGICAL             :: has_zstar, q_in_cryst_coord, loto_disable
  COMPLEX(DP), ALLOCATABLE :: dyn(:, :, :, :), dyn_blk(:, :, :, :), frc_ifc(:, :, :, :)
  COMPLEX(DP), ALLOCATABLE :: z(:, :) 
  REAL(DP), ALLOCATABLE    :: tau(:, :), q(:, :), w2(:, :), wq(:)
  INTEGER, ALLOCATABLE     :: ityp(:), itau_blk(:)
  REAL(DP) ::     omega, alat, &! cell parameters and volume
                  at_blk(3, 3), bg_blk(3, 3),  &! original cell
                  omega_blk,                 &! original cell volume
                  epsil(3, 3),                &! dielectric tensor
                  amass(ntypx),              &! atomic masses
                  amass_blk(ntypx),          &! original atomic masses
                  atws(3, 3),      &! lattice vector for WS initialization
                  rws(0 : 3, nrwsx)   ! nearest neighbor list, rws(0,*) = norm^2
  !
  INTEGER  :: nat, nat_blk, ntyp, ntyp_blk, &
              l1, l2, l3,                   &! supercell dimensions
              nrws,                         &! number of nearest neighbor
              code_group_old

  INTEGER  :: nspin_mag, nqs, ios
  !
  LOGICAL  :: xmlifc, lo_to_split, loto_2d, na_ifc, fd, nosym
  !
  REAL(DP) :: qhat(3), qh, E, qq 
  REAL(DP) :: delta
  REAL(DP), ALLOCATABLE :: xqaux(:, :)
  INTEGER, ALLOCATABLE :: nqb(:)
  INTEGER :: n, i, j, it, nq, nqx, na, nb, nqtot
  INTEGER :: lower_bnd, upper_bnd ! For parallelization 
  LOGICAL, EXTERNAL :: has_xml
  INTEGER, ALLOCATABLE :: num_rap_mode(:, :)
  LOGICAL, ALLOCATABLE :: high_sym(:)
  LOGICAL :: q_in_band_form
  ! .... variables for band plotting based on similarity of eigenvalues
  COMPLEX(DP), ALLOCATABLE :: f_of_q(:, :, :, :)
  INTEGER :: location(1), isig
  CHARACTER(LEN=6) :: int_to_char
  INTEGER            :: npk_label, nch
  CHARACTER(LEN=3), ALLOCATABLE :: letter(:)
  INTEGER, ALLOCATABLE :: label_list(:)
  LOGICAL :: tend, terr
  CHARACTER(LEN=256) :: input_line, buffer
  CHARACTER(LEN= 10) :: point_label_type
  CHARACTER(len=80) :: k_points = 'tpiba'
  ! 
  CHARACTER(LEN=3)         :: atm_zg(ntypx)
  LOGICAL                  :: ZG_conf, synch, incl_qA, q_external
  LOGICAL                  :: ZG_strf, compute_error, single_ph_displ
  INTEGER                  :: dim1, dim2, dim3, niters, qpts_strf
  REAL(DP)                 :: error_thresh, T, qhat_in(3)
  REAL(DP)                 :: atmsf_a(ntypx,5), atmsf_b(ntypx,5) 
  REAL(DP),    ALLOCATABLE :: q_nq(:, :) ! 3, nq
  COMPLEX(DP), ALLOCATABLE :: z_nq(:, :, :) ! nomdes, nmodes, nq
  ! 
  INTEGER                  :: nrots, kres1, kres2, col1, col2, Np
  REAL(DP)                 :: kmin, kmax
  !
  ! for A-SDM, A_ZG
  LOGICAL                  :: ASDM, apply_fd, read_fd_forces, incl_epsil
  LOGICAL                  :: poly, poly_fd_forces, mixing, update_equil
  INTEGER                  :: iter_idx, iter_idx0
  REAL(DP)                 :: fd_displ
  ! for phonon_unfolding
  LOGICAL                  :: ph_unfold
  INTEGER                  :: ng1, ng2, ng3
  !
  NAMELIST /input/ flfrc, amass, asr, at, ntyp, loto_2d, loto_disable, &
       &           q_in_band_form, q_in_cryst_coord, point_label_type,  &
       &           na_ifc, fd, &
! we add the inputs for generating the ZG-configuration
       &           ZG_conf, dim1, dim2, dim3, niters, error_thresh, q_external, & 
       &           compute_error, synch, atm_zg, T, incl_qA, single_ph_displ, & 
       &           ZG_strf, flscf, ph_unfold, qhat_in, ASDM
!
  NAMELIST / A_ZG / apply_fd, read_fd_forces, incl_epsil, iter_idx, iter_idx0, & 
                    fd_displ, poly, poly_fd_forces, mixing, update_equil
!
  NAMELIST /strf_ZG/ atmsf_a, atmsf_b, qpts_strf, &
                     nrots, kres1, kres2, kmin, kmax, col1, col2, Np
!
  NAMELIST / phonon_unfold / ng1, ng2, ng3, dim1, dim2, dim3, & 
                             flfrq, flweights
! Last line of inputs are for the ZG structure factor calculation
  !
  CALL mp_startup()
  CALL environment_start('ZG')
  !
  l1 = 1
  l2 = 1
  l3 = 1
  IF (ionode) CALL input_from_file ( )
     !
     ! ... all calculations are done by the first cpu
     !
     ! set namelist default
     !
     asr = 'no'
     flfrc = ' '
     amass(:) = 0.d0
     amass_blk(:) = 0.d0
     at(:, :) = 0.d0
     ntyp = 0
     q_in_band_form   = .FALSE.
     q_in_cryst_coord = .FALSE.
     point_label_type = 'SC'
     na_ifc           = .FALSE.
     fd               = .FALSE.
     loto_2d          = .FALSE.
     loto_disable     = .FALSE.
     ! 
     flscf         = ' '
     ZG_conf         = .TRUE.
     compute_error   = .TRUE.
     synch           = .FALSE.
     q_external      = .FALSE.
     incl_qA         = .TRUE.
     single_ph_displ = .FALSE.
     T               = 0
     error_thresh    = 3.0E-01
     dim1            = 0
     dim2            = 0
     dim3            = 0
     niters          = 15000 
     atm_zg          = "Element"
     ZG_strf         = .FALSE.
     ph_unfold       = .FALSE.
     ASDM            = .FALSE.
     read_fd_forces  = .FALSE.
     update_equil    = .FALSE.
     apply_fd        = .FALSE.
     poly            = .FALSE.
     poly_fd_forces  = .FALSE.
     mixing          = .FALSE.
     fd_displ        = 0.01058353 ! 0.01058353 in Ang = 0.02 Bohr 
     incl_epsil       = .FALSE. ! for read_fd_forces
     iter_idx        = 1
     iter_idx0       = 0
     qhat_in         = 0.1
     !
     nrots           = 1
     kres1           = 250
     kres2           = 250
     kmin            = -5
     kmax            = 10
     col1            = 1
     col2            = 2
     Np              = 100
     qpts_strf       = 0
     atmsf_a         = 0.d0
     atmsf_b         = 0.d0
     ! 
     ng1             = 10
     ng2             = 10
     ng3             = 10 
     flweights       = 'unfold_weights.dat'
     flfrq           = 'frequencies.dat'
     ! 
     !
     !
     IF (ionode) READ (5, input, IOSTAT = ios)
     CALL mp_bcast(ios, ionode_id, world_comm) 
     CALL errore('ZG', 'reading input namelist', ABS(ios))
     IF ((ionode) .AND. (ZG_strf)) READ (5, strf_ZG , IOSTAT = ios)
     CALL mp_bcast(ios, ionode_id, world_comm)
     CALL errore('strf_ZG', 'reading strf_ZG namelist', ABS(ios))
     IF ((ionode) .AND. (ph_unfold)) READ (5, phonon_unfold , IOSTAT = ios)
     CALL mp_bcast(ios, ionode_id, world_comm)
     CALL errore('ph_unfold', 'reading phonon_unfold namelist', ABS(ios))
     IF ((ionode) .AND. (ASDM)) READ (5, A_ZG , IOSTAT = ios)
     CALL mp_bcast(ios, ionode_id, world_comm)
     CALL errore('ASDM', 'reading A_ZG namelist', ABS(ios))
     !
     CALL mp_bcast(asr, ionode_id, world_comm)
     CALL mp_bcast(flfrc, ionode_id, world_comm)
     CALL mp_bcast(amass, ionode_id, world_comm)
     CALL mp_bcast(amass_blk, ionode_id, world_comm)
     CALL mp_bcast(at, ionode_id, world_comm)
     CALL mp_bcast(ntyp, ionode_id, world_comm)
     CALL mp_bcast(na_ifc,ionode_id, world_comm) 
     CALL mp_bcast(fd,ionode_id, world_comm)
     CALL mp_bcast(q_in_band_form, ionode_id, world_comm)
     CALL mp_bcast(q_in_cryst_coord, ionode_id, world_comm)
     CALL mp_bcast(point_label_type, ionode_id, world_comm)
     CALL mp_bcast(loto_2d, ionode_id, world_comm) 
     CALL mp_bcast(loto_disable,ionode_id, world_comm)
     ! 
     CALL mp_bcast(ZG_conf, ionode_id, world_comm)
     CALL mp_bcast(flscf, ionode_id, world_comm)
     CALL mp_bcast(compute_error, ionode_id, world_comm)
     CALL mp_bcast(synch, ionode_id, world_comm)
     CALL mp_bcast(q_external, ionode_id, world_comm)
     CALL mp_bcast(incl_qA, ionode_id, world_comm)
     CALL mp_bcast(single_ph_displ, ionode_id, world_comm)
     CALL mp_bcast(T, ionode_id, world_comm)
     CALL mp_bcast(error_thresh, ionode_id, world_comm)
     CALL mp_bcast(dim1, ionode_id, world_comm)
     CALL mp_bcast(dim2, ionode_id, world_comm)
     CALL mp_bcast(dim3, ionode_id, world_comm)
     CALL mp_bcast(niters, ionode_id, world_comm)
     CALL mp_bcast(atm_zg, ionode_id, world_comm)
     CALL mp_bcast(qhat_in, ionode_id, world_comm)
     CALL mp_bcast(ZG_strf, ionode_id, world_comm)
     CALL mp_bcast(ph_unfold, ionode_id, world_comm)
     CALL mp_bcast(ASDM, ionode_id, world_comm)
     CALL mp_bcast(apply_fd, ionode_id, world_comm)
     CALL mp_bcast(fd_displ, ionode_id, world_comm)
     CALL mp_bcast(read_fd_forces, ionode_id, world_comm)
     CALL mp_bcast(update_equil, ionode_id, world_comm)
     CALL mp_bcast(incl_epsil, ionode_id, world_comm)
     CALL mp_bcast(iter_idx, ionode_id, world_comm)
     CALL mp_bcast(iter_idx0, ionode_id, world_comm)
     CALL mp_bcast(poly, ionode_id, world_comm)
     CALL mp_bcast(poly_fd_forces, ionode_id, world_comm)
     CALL mp_bcast(mixing, ionode_id, world_comm)
     !
     CALL mp_bcast(qpts_strf, ionode_id, world_comm)
     CALL mp_bcast(atmsf_a, ionode_id, world_comm)
     CALL mp_bcast(atmsf_b, ionode_id, world_comm)
     CALL mp_bcast(nrots, ionode_id, world_comm)
     CALL mp_bcast(kres1, ionode_id, world_comm)
     CALL mp_bcast(kres2, ionode_id, world_comm)
     CALL mp_bcast(kmin, ionode_id, world_comm)
     CALL mp_bcast(kmax, ionode_id, world_comm)
     CALL mp_bcast(col1, ionode_id, world_comm)
     CALL mp_bcast(col2, ionode_id, world_comm)
     CALL mp_bcast(Np, ionode_id, world_comm)
     ! 
     CALL mp_bcast(ng1, ionode_id, world_comm)
     CALL mp_bcast(ng2, ionode_id, world_comm)
     CALL mp_bcast(ng3, ionode_id, world_comm)
     CALL mp_bcast(flfrq, ionode_id, world_comm)
     CALL mp_bcast(flweights, ionode_id, world_comm)
     !
     IF (loto_2d .AND. loto_disable) CALL errore('ZG', &
         'loto_2d and loto_disable cannot be both true', 1)
     ! 
     IF (ZG_conf .AND. ph_unfold) CALL errore('ZG', &
         'ZG_conf and ph_unfold cannot be both true', 1)
     ! 
     ! To check that user specifies supercell dimensions
     IF (ZG_conf) THEN 
        IF (dim1 < 1) CALL errore('ZG', 'reading supercell size, check dim1', 1)
        IF (dim2 < 1) CALL errore('ZG', 'reading supercell size, check dim2', 1)
        IF (dim3 < 1) CALL errore('ZG', 'reading supercell size, check dim3', 1)
        IF ( single_ph_displ .AND. compute_error  ) CALL errore('ZG', & 
             "for single phonon displacements set 'compute_error' to false", 1)
     ENDIF
     !
     IF (ASDM .AND. apply_fd .AND. read_fd_forces) CALL errore('ASDM', &
         '"apply_fd" and "read_fd_forces" cannot be both true', 1)
     IF (ASDM .AND. iter_idx < 0 ) CALL errore('ASDM', &
         'provide "iter_idx" >= 0', 1)
     IF (ASDM .AND. iter_idx0 < 0 ) CALL errore('ASDM', &
         'provide "iter_idx0" >= 0', 1)
     IF (ASDM .AND. (iter_idx .LE. iter_idx0)) CALL errore('ASDM', &
         'provide "iter_idx" > iter_idx0', 1)
     IF (ASDM .AND. poly) THEN 
       iter_idx = 0
       iter_idx0 = 0
       apply_fd = .FALSE.
       read_fd_forces = .FALSE.
     ENDIF
     IF (ASDM .AND. poly_fd_forces .AND. .NOT. poly) CALL errore('ASDM', &
         'for "poly_fd_forces", turn flag "poly" true', 1)
     !
     IF (ph_unfold) THEN 
       IF ((dim1 < 1)  .OR. (dim2 < 1) .OR. (dim3 < 1)) CALL errore('ph_unfold', &
            'reading supercell size', 1) 
       IF ((ng1 < 1)  .OR. (ng2 < 1) .OR. (ng3 < 1)) CALL errore('ph_unfold', & 
            'wrong g-vectors', 1)
       q_external = .true.
       !IF (.NOT. q_external) CALL errore('ph_unfold', 'set q_external = .true.', 1)
     ENDIF 
     !
     ! read force constants
     !
     ntyp_blk = ntypx ! avoids fake out-of-bound error
     xmlifc=has_xml(flfrc)
     IF (xmlifc) THEN
        CALL read_dyn_mat_param(flfrc, ntyp_blk, nat_blk)
        ALLOCATE (m_loc(3, nat_blk))
        ALLOCATE (tau_blk(3, nat_blk))
        ALLOCATE (ityp_blk(nat_blk))
        ALLOCATE (atm(ntyp_blk))
        ALLOCATE (zeu(3, 3, nat_blk))
        CALL read_dyn_mat_header(ntyp_blk, nat_blk, ibrav, nspin_mag, &
                 celldm, at_blk, bg_blk, omega_blk, atm, amass_blk, &
                 tau_blk, ityp_blk,  m_loc, nqs, has_zstar, epsil, zeu )
        alat= celldm(1)
        call volume(alat, at_blk(1, 1), at_blk(1, 2), at_blk(1, 3), omega_blk)
        CALL read_ifc_param(nr1, nr2, nr3)
        ALLOCATE(frc(nr1, nr2, nr3, 3, 3, nat_blk, nat_blk))
        CALL read_ifc(nr1, nr2, nr3, nat_blk,frc)
     ELSE
        CALL readfc ( flfrc, nr1, nr2, nr3, epsil, nat_blk, &
            ibrav, alat, at_blk, ntyp_blk, &
            amass_blk, omega_blk, has_zstar)
     ENDIF
     !
     CALL recips ( at_blk(1, 1), at_blk(1, 2), at_blk(1, 3),  &
          bg_blk(1, 1), bg_blk(1, 2), bg_blk(1, 3) )
     !
     ! set up (super)cell
     !
     IF (ntyp < 0) THEN
        call errore ('ZG','wrong ntyp ', ABS(ntyp))
     ELSE IF (ntyp == 0) THEN
        ntyp =ntyp_blk
     ENDIF
     !
     ! masses (for mass approximation)
     !
     DO it= 1, ntyp
        IF (amass(it) < 0.d0) THEN
           CALL errore ('ZG','wrong mass in the namelist', it)
        ELSE IF (amass(it) == 0.d0) THEN
           IF (it.LE.ntyp_blk) THEN
              WRITE (stdout,'(a, i3, a, a)') ' mass for atomic type ', it,      &
                   &                     ' not given; uses mass from file ',flfrc
              amass(it) = amass_blk(it)
           ELSE
              CALL errore ('ZG','missing mass in the namelist', it)
           ENDIF
        ENDIF
     ENDDO
     !
     ! lattice vectors
     !
     IF (SUM(ABS(at(:, :))) == 0.d0) THEN
        IF (l1.LE.0 .OR. l2.LE.0 .OR. l3.LE.0) CALL                    &
             &             errore ('ZG',' wrong l1,l2 or l3', 1)
        at(:, 1) = at_blk(:, 1) * DBLE(l1)
        at(:, 2) = at_blk(:, 2) * DBLE(l2)
        at(:, 3) = at_blk(:, 3) * DBLE(l3)
     ENDIF
     !
     CALL check_at(at, bg_blk, alat, omega)
     !
     ! the supercell contains "nsc" times the original unit cell
     !
     nsc = NINT(omega / omega_blk)
     IF (ABS(omega / omega_blk-nsc) > eps) &
          CALL errore ('ZG', 'volume ratio not integer', 1)
     !
     ! read/generate atomic positions of the (super)cell
     !
     nat = nat_blk * nsc
     !
     ALLOCATE ( tau (3, nat), ityp(nat), itau_blk(nat) )
     !
     ! read atomic positions from IFC file
     CALL set_tau  &
             (nat, nat_blk, at, at_blk, tau, tau_blk, ityp, ityp_blk, itau_blk)
     !
     !
     ! reciprocal lattice vectors
     !
     CALL recips (at(1, 1), at(1, 2), at(1, 3), bg(1, 1), bg(1, 2), bg(1, 3))
     !
     ! build the WS cell corresponding to the force constant grid
     !
     atws(:, 1) = at_blk(:, 1) * DBLE(nr1)
     atws(:, 2) = at_blk(:, 2) * DBLE(nr2)
     atws(:, 3) = at_blk(:, 3) * DBLE(nr3)
     ! initialize WS r-vectors
     CALL wsinit(rws, nrwsx, nrws, atws)
     !
     ! end of (super)cell setup
     !
     !
     ! read q-point list
     !
     ! 
     IF (.NOT. q_external) THEN 
       CALL qpoint_gen1(dim1, dim2, dim3, nq) 
       ! nq = ctrAB
       CALL mp_bcast(nq, ionode_id, world_comm)
       !
       ALLOCATE ( q(3, nq) )        
       CALL qpoint_gen2(dim1, dim2, dim3, nq, q) 
       !
       CALL mp_bcast(q, ionode_id, world_comm)
       !
       CALL cryst_to_cart(nq, q, bg, +1) ! convert them to Cartesian
     ELSE
     ! 
       IF (ionode) READ (5, *) nq
       CALL mp_bcast(nq, ionode_id, world_comm)
       ALLOCATE ( q(3, nq) )
       IF (.NOT.q_in_band_form) THEN
             DO n = 1, nq
                IF (ionode) READ (5, *) (q(i, n), i = 1, 3)
             ENDDO
             CALL mp_bcast(q, ionode_id, world_comm)
             !
             IF (q_in_cryst_coord)  CALL cryst_to_cart(nq, q, bg, +1)
       ELSE
             ALLOCATE( nqb(nq) )
             ALLOCATE( xqaux(3, nq) )
             ALLOCATE( letter(nq) )
             ALLOCATE( label_list(nq) )
             npk_label= 0
             DO n = 1, nq
                CALL read_line( input_line, end_of_file = tend, error = terr )
                IF (tend) CALL errore('ZG','Missing lines', 1)
                IF (terr) CALL errore('ZG','Error reading q points', 1)
                DO j = 1, 256   ! loop over all characters of input_line
                   IF ( (ICHAR(input_line(j:j)) < 58 .AND. &   ! a digit
                         ICHAR(input_line(j:j)) > 47)      &
                     .OR.ICHAR(input_line(j:j)) == 43 .OR. &   ! the + sign
                         ICHAR(input_line(j:j)) == 45 .OR. &   ! the - sign
                         ICHAR(input_line(j:j)) == 46 ) THEN   ! a dot .
  !
  !   This is a digit, therefore this line contains the coordinates of the
  !   k point. We read it and EXIT from the loop on characters
  !
                       READ(input_line,*) xqaux(1, n), xqaux(2, n), xqaux(3, n), &
                                                      nqb(n)
                       EXIT
                   ELSEIF ((ICHAR(input_line(j:j)) < 123 .AND. &
                            ICHAR(input_line(j:j)) > 64))  THEN
  !
  !   This is a letter, not a space character. We read the next three 
  !   characters and save them in the letter array, save also which k point
  !   it is
  !
                      npk_label=npk_label+ 1
                      READ(input_line(j:),'(a3)') letter(npk_label)
                      label_list(npk_label) =n
  !
  !  now we remove the letters from input_line and read the number of points
  !  of the line. The next two line should account for the case in which
  !  there is only one space between the letter and the number of points.
  !
                      nch=3
                      IF ( ICHAR(input_line(j+ 1:j+ 1)) ==32 .OR. &
                           ICHAR(input_line(j+2:j+2)) ==32 ) nch=2
                      buffer =input_line(j+nch:)
                      READ(buffer,*, err =20, iostat=ios) nqb(n)
  20                  IF (ios /= 0) CALL errore('ZG',&
                                        'problem reading number of points', 1)
                      EXIT
                   ENDIF
                ENDDO
             ENDDO
             IF (q_in_cryst_coord) k_points ='crystal'
             IF ( npk_label > 0 ) &
                CALL transform_label_coord(ibrav, celldm, xqaux, letter, &
                     label_list, npk_label, nq, k_points, point_label_type )
  
             DEALLOCATE(letter)
             DEALLOCATE(label_list)
  
             CALL mp_bcast(xqaux, ionode_id, world_comm)
             CALL mp_bcast(nqb, ionode_id, world_comm)
             IF (q_in_cryst_coord)  CALL cryst_to_cart(nq,xqaux, bg,+ 1)
             nqtot=SUM(nqb(1 : nq - 1)) + 1
             DO i = 1, nq - 1
                IF (nqb(i) == 0) nqtot=nqtot+ 1
             ENDDO
             DEALLOCATE(q)
             ALLOCATE(q(3, nqtot))
             ALLOCATE(wq(nqtot))
             CALL generate_k_along_lines(nq, xqaux, nqb, q, wq, nqtot)
             nq = nqtot
             DEALLOCATE(xqaux)
             DEALLOCATE(nqb)
          ENDIF
          ! 
     ENDIF ! q_external, q-list
     !
     IF (asr /= 'no') THEN
        CALL set_asr (asr, nr1, nr2, nr3, frc, zeu, &
             nat_blk, ibrav, tau_blk)
     ENDIF
     !
     ALLOCATE ( dyn(3, 3, nat, nat), dyn_blk(3, 3, nat_blk, nat_blk) )
     ALLOCATE ( z(3 * nat, 3 * nat), w2(3 * nat, nq), f_of_q(3, 3, nat, nat) )
     ! 
     ! Have to initialize w2
     w2 = 0.d0
     IF (ZG_conf .OR. ph_unfold) THEN
       ALLOCATE ( z_nq(3 * nat, 3 * nat, nq), q_nq(3, nq))
       z_nq(:, :, :) = (0.d0, 0.d0)
       q_nq(:, :) = 0.d0
     ENDIF
     ! 

     IF (xmlifc) CALL set_sym(nat, tau, ityp, nspin_mag, m_loc )

     ALLOCATE(num_rap_mode(3 * nat, nq))
     ALLOCATE(high_sym(nq))
     num_rap_mode=- 1
     high_sym=.TRUE.
     !
     CALL fkbounds( nq, lower_bnd, upper_bnd )
     !
     DO n = lower_bnd, upper_bnd ! 1, nq
     !DO n = 1, nq
        dyn(:, :, :, :) = (0.d0, 0.d0)

        lo_to_split = .FALSE.
        f_of_q(:, :, :, :) = (0.d0, 0.d0)

        IF(na_ifc) THEN

           qq=SQRT(q(1,n)**2+q(2,n)**2+q(3,n)**2)
           if(ABS(qq) < 1d-8) qq= 1.0
           qhat(1)=q(1,n)/qq
           qhat(2)=q(2,n)/qq
           qhat(3)=q(3,n)/qq

           CALL nonanal_ifc (nat,nat_blk,itau_blk,epsil,qhat,zeu,omega,dyn, &
                           nr1, nr2, nr3,f_of_q)
        ENDIF
        
        CALL setupmat (q(1, n), dyn, nat, at, bg, tau, itau_blk, nsc, alat, &
             dyn_blk, nat_blk, at_blk, bg_blk, tau_blk, omega_blk,  &
                   loto_2d, &
             epsil, zeu, frc, nr1,nr2,nr3, has_zstar, rws, nrws, na_ifc, f_of_q, fd)

        IF (.not.loto_2d) THEN
        qhat(1) = q(1, n) * at(1, 1) + q(2, n) * at(2, 1) + q(3, n) * at(3, 1)
        qhat(2) = q(1, n) * at(1, 2) + q(2, n) * at(2, 2) + q(3, n) * at(3, 2)
        qhat(3) = q(1, n) * at(1, 3) + q(2, n) * at(2, 3) + q(3, n) * at(3, 3)
        IF ( ABS( qhat(1) - NINT (qhat(1) ) ) <= eps .AND. &
             ABS( qhat(2) - NINT (qhat(2) ) ) <= eps .AND. &
             ABS( qhat(3) - NINT (qhat(3) ) ) <= eps ) THEN
           !
           ! q = 0 : we need the direction q => 0 for the non-analytic part
           !
           IF ( n == 1 ) THEN
              ! IF q is the first point in the list
              IF ( nq > 1 ) THEN
                 ! one more point
                 qhat(:) = q(:, n) - q(:, n + 1)
              ELSE
                 ! no more points
                 qhat(:) = 0.d0
              ENDIF
           ELSE IF ( n > 1 ) THEN
              ! IF q is not the first point in the list
              IF ( q(1, n - 1) == 0.d0 .AND. &
                   q(2, n - 1) == 0.d0 .AND. &
                   q(3, n - 1) == 0.d0 .AND. n < nq ) THEN
                 ! IF the preceding q is also 0 :
                 qhat(:) = q(:, n) - q(:, n + 1)
              ELSE
                 ! IF the preceding q is npt 0 :
                 qhat(:) = q(:, n) - q(:, n - 1)
              ENDIF
           ENDIF
           ! consider qhat from input or default when dim1=dim2=dim3=1 and
           ! LO-TO splitting needs to be included.
           IF (dim1 == 1 .AND. dim2 == 1 .AND. dim3 == 1) THEN
               qhat = qhat_in
               IF (ionode) THEN
                  WRITE(*,*) "=============================================="
                  WRITE(*,*)
                  WRITE(*,*) "Direction for q is specified from input:"
                  WRITE(*,'(3F8.4)') qhat
                  WRITE(*,*)
                  WRITE(*,*) "=============================================="
               ENDIF
           ENDIF
           qh = SQRT(qhat(1)**2 + qhat(2)**2 + qhat(3)**2)
           ! WRITE(*,*) ' qh,  has_zstar ',qh,  has_zstar
           IF (qh /= 0.d0) qhat(:) = qhat(:) / qh
           IF (qh /= 0.d0 .AND. .NOT. has_zstar) THEN
                IF (ionode) WRITE(*,*)
                CALL infomsg  &
                ('ZG','Z* not found in file '//TRIM(flfrc)// &
                          ', TO-LO splitting at q = 0 will be absent!')
           ELSEIF (loto_disable) THEN
                CALL infomsg('ZG', &
                    'loto_disable is true. Disable LO-TO splitting at q=0.')
           ELSE
              lo_to_split=.TRUE.
           ENDIF
           !
           IF (lo_to_split) CALL nonanal (nat, nat_blk, itau_blk, epsil, qhat, zeu, omega, dyn)
           !
        ENDIF
        !
        END IF
        !
        CALL dyndiag(nat, ntyp, amass, ityp, dyn, w2(1, n), z)
        !
        ! Fill a 3D matrix with all eigenvectors
        !
        IF (ZG_conf .OR. ph_unfold) THEN
           z_nq(:, :, n) = z(:, :)               
           q_nq(:, n) = q(:, n)
        ENDIF
        !
        ! Cannot use the small group of \Gamma to analize the symmetry
        ! of the mode IF there is an electric field.
        !
        IF (xmlifc.AND..NOT.lo_to_split) THEN
             WRITE(stdout,'(10x,"xq=", 3F8.4)') q(:, n)
             CALL find_representations_mode_q(nat, ntyp, q(:, n), &
                       w2(:, n), z,tau, ityp, amass, num_rap_mode(:, n), nspin_mag)
            IF (code_group == code_group_old.OR.high_sym(n- 1)) high_sym(n) =.FALSE.
            code_group_old= code_group
        ENDIF
        !
        !
        !
     ENDDO  !nq
     !
     !
     CALL mp_sum(z_nq, inter_pool_comm)
     CALL mp_sum(q_nq, inter_pool_comm)
     CALL mp_sum(w2, inter_pool_comm)
     CALL mp_barrier(inter_pool_comm)
     !
     !  Here is the main subroutine for generating ZG displacements.
     !
     IF ( ZG_conf .AND. .NOT. ph_unfold) & 
          CALL ZG_configuration(nq, nat, ntyp, amass, ityp, q_nq, w2, z_nq, & 
                                dim1, dim2, dim3, niters, error_thresh, &
                                synch, tau, alat, atm_zg, ntypx, at, &
                                q_in_cryst_coord, q_external, T, incl_qA, & 
                                compute_error, single_ph_displ, & 
                                ZG_strf, qpts_strf, atmsf_a, atmsf_b, &
                                nrots, kres1, kres2, kmin, kmax, col1, col2, Np, & 
                                flscf, ASDM, apply_fd, read_fd_forces, incl_epsil, & 
                                iter_idx, iter_idx0, fd_displ, ibrav, epsil, zeu, & 
                                poly, poly_fd_forces, mixing, update_equil)
     ! 
     !
     IF (ph_unfold .AND. .NOT. ZG_conf) & 
         CALL phonon_unfolding(nq, tau, nat, ntyp, amass, ityp, & 
                               dim1, dim2, dim3, ng1, ng2, ng3, q_nq, z_nq, &
                               flfrq, flweights, w2) 
     !
     DEALLOCATE (z, w2, dyn, dyn_blk)
     ! 
     IF (ZG_conf .OR. ph_unfold) DEALLOCATE (z_nq, q_nq) 
     ! 
     !
     !    for a2F
     !
     DEALLOCATE(num_rap_mode)
     DEALLOCATE(high_sym)
  !

  CALL environment_end('ZG')
  !
  CALL mp_global_end()
  !
  STOP
  !
END PROGRAM ZG
!
!-----------------------------------------------------------------------
SUBROUTINE readfc ( flfrc, nr1, nr2, nr3, epsil, nat,    &
                    ibrav, alat, at, ntyp, amass, omega, has_zstar )
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE ifconstants,ONLY : tau => tau_blk, ityp => ityp_blk, frc, zeu
  USE cell_base,  ONLY : celldm
  USE io_global,  ONLY : ionode, ionode_id, stdout
  USE mp,         ONLY : mp_bcast 
  USE mp_world,   ONLY : world_comm 
  USE constants,  ONLY : amu_ry
  !
  IMPLICIT NONE
  ! I/O variable
  CHARACTER(LEN=256) :: flfrc
  INTEGER :: ibrav, nr1, nr2, nr3, nat, ntyp
  REAL(DP) :: alat, at(3, 3), epsil(3, 3)
  LOGICAL :: has_zstar
  ! local variables
  INTEGER :: i, j, na, nb, m1,m2,m3
  INTEGER :: ibid, jbid, nabid, nbbid, m1bid,m2bid,m3bid
  REAL(DP) :: amass(ntyp), amass_from_file, omega
  INTEGER :: nt
  CHARACTER(LEN=3) :: atm
  !
  !
  IF (ionode) OPEN (unit= 1,file=flfrc, status ='old',form='formatted')
  !
  !  read cell data
  !
  IF (ionode)THEN
     READ(1,*) ntyp, nat, ibrav,(celldm(i), i = 1,6)
     IF (ibrav== 0) THEN
        read(1,*) ((at(i, j), i = 1, 3), j =1, 3)
     ENDIF
  ENDIF
  CALL mp_bcast(ntyp, ionode_id, world_comm)
  CALL mp_bcast(nat, ionode_id, world_comm)
  CALL mp_bcast(ibrav, ionode_id, world_comm)
  CALL mp_bcast(celldm, ionode_id, world_comm)
  IF (ibrav== 0) THEN
     CALL mp_bcast(at, ionode_id, world_comm)
  ENDIF
  !
  CALL latgen(ibrav, celldm, at(1, 1), at(1, 2), at(1, 3), omega)
  alat = celldm(1)
  at = at / alat !  bring at in units of alat
  CALL volume(alat, at(1, 1), at(1, 2), at(1, 3), omega)
  !
  !  read atomic types, positions and masses
  !
  DO nt = 1, ntyp
     IF (ionode) READ(1,*) i, atm, amass_from_file
     CALL mp_bcast(i, ionode_id, world_comm)
     CALL mp_bcast(atm, ionode_id, world_comm)
     CALL mp_bcast(amass_from_file, ionode_id, world_comm)
     IF (i.NE.nt) CALL errore ('readfc','wrong data read', nt)
     IF (amass(nt).EQ.0.d0) THEN
        amass(nt) = amass_from_file/amu_ry
     ELSE
        WRITE(stdout,*) 'for atomic type', nt,' mass from file not used'
     ENDIF
  ENDDO
  !
  ALLOCATE (tau(3, nat), ityp(nat), zeu(3, 3, nat))
  !
  DO na= 1, nat
     IF (ionode) READ(1,*) i, ityp(na),(tau(j, na), j = 1, 3)
     CALL mp_bcast(i, ionode_id, world_comm)
     IF (i.NE.na) CALL errore ('readfc','wrong data read', na)
  ENDDO
  CALL mp_bcast(ityp, ionode_id, world_comm)
  CALL mp_bcast(tau, ionode_id, world_comm)
  !
  !  read macroscopic variable
  !
  IF (ionode) READ (1,*) has_zstar
  CALL mp_bcast(has_zstar, ionode_id, world_comm)
  IF (has_zstar) THEN
     IF (ionode) READ(1,*) ((epsil(i, j), j = 1, 3), i =1, 3)
     CALL mp_bcast(epsil, ionode_id, world_comm)
     IF (ionode) THEN
        DO na= 1, nat
           READ(1,*)
           READ(1,*) ((zeu(i, j, na), j = 1, 3), i =1, 3)
        ENDDO
     ENDIF
     CALL mp_bcast(zeu, ionode_id, world_comm)
  ELSE
     zeu  (:, :, :) = 0.d0
     epsil(:, :) = 0.d0
  ENDIF
  !
  IF (ionode) READ (1,*) nr1, nr2, nr3
  CALL mp_bcast(nr1, ionode_id, world_comm)
  CALL mp_bcast(nr2, ionode_id, world_comm)
  CALL mp_bcast(nr3, ionode_id, world_comm)
  !
  !  read REAL-space interatomic force constants
  !
  ALLOCATE ( frc(nr1, nr2, nr3, 3, 3, nat, nat) )
  frc(:, :, :, :, :, :, :) = 0.d0
  DO i = 1, 3
     DO j = 1, 3
        DO na= 1, nat
           DO nb= 1, nat
              IF (ionode) READ (1,*) ibid, jbid, nabid, nbbid
              CALL mp_bcast(ibid, ionode_id, world_comm)
              CALL mp_bcast(jbid, ionode_id, world_comm)
              CALL mp_bcast(nabid, ionode_id, world_comm)
              CALL mp_bcast(nbbid, ionode_id, world_comm)
              IF(i .NE.ibid  .OR. j .NE.jbid .OR.                   &
                 na.NE.nabid .OR. nb.NE.nbbid)                      &
                 CALL errore  ('readfc','error in reading', 1)
              IF (ionode) READ (1,*) (((m1bid, m2bid, m3bid,        &
                          frc(m1,m2,m3, i, j, na, nb),                  &
                           m1= 1, nr1),m2 =1, nr2),m3=1, nr3)
               
              CALL mp_bcast(frc(:, :, :, i, j, na, nb), ionode_id, world_comm)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  IF (ionode) CLOSE(unit= 1)
  !
  RETURN
END SUBROUTINE readfc
!
!-----------------------------------------------------------------------
SUBROUTINE frc_blk(dyn,q,tau, nat, nr1, nr2, nr3, frc, & 
                    at, bg, rws, nrws, f_of_q, fd)
  !-----------------------------------------------------------------------
  ! calculates the dynamical matrix at q from the (short-range part of the)
  ! force constants
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi
  USE io_global,  ONLY : ionode, ionode_id, stdout ! mz adds ionode and ionode_id
  USE mp_world,   ONLY : world_comm
  USE mp,         ONLY : mp_bcast ! mz adds mp_bcast and world_comm
  !
  IMPLICIT NONE
  INTEGER nr1, nr2, nr3, nat, n1, n2, n3, nr1_, nr2_, nr3_, &
          ipol, jpol, na, nb, m1, m2, m3, NINT, i, j, nrws 
  COMPLEX(DP) dyn(3, 3, nat, nat), f_of_q(3, 3, nat, nat)
  REAL(DP) frc(nr1, nr2, nr3, 3, 3, nat, nat), tau(3, nat), q(3), arg, &
               at(3, 3), bg(3, 3), r(3), weight, r_ws(3),  &
               total_weight, rws(0:3, nrws), alat
  REAL(DP), EXTERNAL :: wsweight
  REAL(DP),SAVE,ALLOCATABLE :: wscache(:, :, :, :, :)
  REAL(DP), ALLOCATABLE :: ttt(:, :, :, :, :), tttx(:, :)
  LOGICAL,SAVE :: first=.TRUE.
  LOGICAL      :: fd
  !
  nr1_=2*nr1
  nr2_=2*nr2
  nr3_=2*nr3
  FIRST_TIME : IF (first) THEN
    first=.FALSE.
    ALLOCATE( wscache(-nr3_:nr3_, -nr2_:nr2_, -nr1_:nr1_, nat, nat) )
    DO na= 1, nat
       DO nb= 1, nat
          total_weight= 0.0d0
          !
          DO n1=-nr1_, nr1_
             DO n2 =-nr2_, nr2_
                DO n3=-nr3_, nr3_
                   DO i = 1, 3
                      r(i) = n1*at(i, 1) +n2*at(i, 2) +n3*at(i, 3)
                      r_ws(i) = r(i) + tau(i, na) -tau(i, nb)
                      if (fd) r_ws(i) = r(i) + tau(i,nb)-tau(i,na)
                   ENDDO
                   wscache(n3, n2, n1, nb, na) = wsweight(r_ws, rws, nrws)
                   CALL mp_bcast(wscache(n3,n2,n1,nb,na),ionode_id, world_comm)
                ENDDO
             ENDDO
          ENDDO
      ENDDO
    ENDDO
  ENDIF FIRST_TIME
  !
  ALLOCATE(ttt(3, nat, nr1, nr2, nr3))
  ALLOCATE(tttx(3, nat*nr1*nr2*nr3))
  ttt(:, :, :, :, :) = 0.d0

  DO na = 1, nat
     DO nb= 1, nat
        total_weight= 0.0d0
        DO n1=-nr1_, nr1_
           DO n2 =-nr2_, nr2_
              DO n3=-nr3_, nr3_
                 !
                 ! SUM OVER R VECTORS IN THE SUPERCELL - VERY VERY SAFE RANGE!
                 !
                 DO i = 1, 3
                    r(i) = n1*at(i, 1) +n2*at(i, 2) +n3*at(i, 3)
                 ENDDO

                 weight = wscache(n3, n2, n1, nb, na) 
                 IF (weight > 0.0d0) THEN
                    !
                    ! FIND THE VECTOR CORRESPONDING TO R IN THE ORIGINAL CELL
                    !
                    m1 = MOD(n1 + 1, nr1)
                    IF(m1.LE.0) m1=m1 +nr1
                    m2 = MOD(n2 + 1, nr2)
                    IF(m2.LE.0) m2 =m2 +nr2
                    m3 = MOD(n3 + 1, nr3)
                    IF(m3.LE.0) m3=m3 +nr3
                 !   WRITE(*,'(6i4)') n1, n2, n3,m1,m2,m3
                    !
                    ! FOURIER TRANSFORM
                    !
                    DO i = 1, 3
                       ttt(i, na,m1,m2,m3) =tau(i, na) +m1*at(i, 1) +m2*at(i, 2) +m3*at(i, 3)
                       ttt(i, nb,m1,m2,m3) =tau(i, nb) +m1*at(i, 1) +m2*at(i, 2) +m3*at(i, 3)
                    ENDDO

                    arg = tpi* (q(1) *r(1) + q(2) *r(2) + q(3) *r(3))
                    DO ipol= 1, 3
                       DO jpol= 1, 3
                          dyn(ipol, jpol, na, nb) =                 &
                               dyn(ipol, jpol, na, nb) +            &
                               (frc(m1,m2,m3, ipol, jpol, na, nb) +f_of_q(ipol, jpol, na, nb))     &
                               *CMPLX(COS(arg),-SIN(arg), kind= DP) *weight
                       ENDDO
                    ENDDO
                 ENDIF
                 total_weight=total_weight + weight
              ENDDO
           ENDDO
        ENDDO
        IF (ABS(total_weight-nr1*nr2*nr3).GT.1.0d-8) THEN
          ! WRITE(stdout,*) total_weight
          ! CALL errore ('frc_blk','wrong total_weight', 1)
        ENDIF
     ENDDO
  ENDDO
  !
  RETURN
END SUBROUTINE frc_blk
!
!-----------------------------------------------------------------------
SUBROUTINE setupmat (q,dyn,nat,at,bg,tau,itau_blk,nsc,alat, &
     &         dyn_blk,nat_blk,at_blk,bg_blk,tau_blk,omega_blk, &
     &         loto_2d, &
     &         epsil,zeu,frc,nr1,nr2,nr3,has_zstar,rws,nrws,na_ifc,f_of_q,fd)
  !-----------------------------------------------------------------------
  ! compute the dynamical matrix (the analytic part only)
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi
  USE cell_base,  ONLY : celldm
  USE rigid,      ONLY : rgd_blk
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER:: nr1, nr2, nr3, nat, nat_blk, nsc, nrws, itau_blk(nat)
  REAL(DP) :: q(3), tau(3, nat), at(3, 3), bg(3, 3), alat,      &
                  epsil(3, 3), zeu(3, 3, nat_blk), rws(0:3, nrws),   &
                  frc(nr1, nr2, nr3, 3, 3, nat_blk, nat_blk)
  REAL(DP) :: tau_blk(3, nat_blk), at_blk(3, 3), bg_blk(3, 3), omega_blk
  COMPLEX(DP) dyn_blk(3, 3, nat_blk, nat_blk), f_of_q(3, 3, nat, nat)
  COMPLEX(DP) ::  dyn(3, 3, nat, nat)
  LOGICAL :: has_zstar, na_ifc, fd
  !
  ! local variables
  !
  REAL(DP) :: arg
  COMPLEX(DP) :: cfac(nat)
  INTEGER :: i, j, k, na, nb, na_blk, nb_blk, iq
  REAL(DP) :: qp(3), qbid(3, nsc) ! automatic array
  LOGICAL ::  loto_2d
  !
  !
  CALL q_gen(nsc,qbid, at_blk, bg_blk, at, bg)
  !
  DO iq= 1, nsc
     !
     DO k = 1, 3
        qp(k) = q(k) + qbid(k, iq)
     ENDDO
     !
     dyn_blk(:, :, :, :) = (0.d0,0.d0)
     CALL frc_blk (dyn_blk, qp,tau_blk, nat_blk,              &
          &              nr1, nr2, nr3,frc, at_blk, bg_blk, rws, nrws, f_of_q, fd)
      IF (has_zstar .and. .not. na_ifc) &
           CALL rgd_blk(nr1, nr2, nr3, nat_blk, dyn_blk, qp,tau_blk,   &
                        epsil, zeu, bg_blk, omega_blk, celldm(1), loto_2d, +1.d0)
     !
     DO na= 1, nat
        na_blk = itau_blk(na)
        DO nb= 1, nat
           nb_blk = itau_blk(nb)
           !
           arg=tpi* ( qp(1) * ( (tau(1, na) -tau_blk(1, na_blk)) -   &
                                (tau(1, nb) -tau_blk(1, nb_blk)) ) + &
                      qp(2) * ( (tau(2, na) -tau_blk(2, na_blk)) -   &
                                (tau(2, nb) -tau_blk(2, nb_blk)) ) + &
                      qp(3) * ( (tau(3, na) -tau_blk(3, na_blk)) -   &
                                (tau(3, nb) -tau_blk(3, nb_blk)) ) )
           !
           cfac(nb) = CMPLX(COS(arg),SIN(arg), kind= DP)/nsc
           !
        ENDDO ! nb
        !
        DO i = 1, 3
           DO j = 1, 3
              !
              DO nb= 1, nat
                 nb_blk = itau_blk(nb)
                 dyn(i, j, na, nb) = dyn(i, j, na, nb) + cfac(nb) * &
                      dyn_blk(i, j, na_blk, nb_blk)
              ENDDO ! nb
              !
           ENDDO ! j
        ENDDO ! i
     ENDDO ! na
     !
  ENDDO ! iq
  !
  RETURN
END SUBROUTINE setupmat
!
!
!----------------------------------------------------------------------
SUBROUTINE set_asr (asr, nr1, nr2, nr3, frc, zeu, nat, ibrav, tau)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : ionode, stdout
  !
  IMPLICIT NONE
  CHARACTER (LEN= 10), INTENT(in) :: asr
  INTEGER, INTENT(in) :: nr1, nr2, nr3, nat, ibrav
  REAL(DP), INTENT(in) :: tau(3, nat)
  REAL(DP), INTENT(inout) :: frc(nr1, nr2, nr3, 3, 3, nat, nat), zeu(3, 3, nat)
  !
  INTEGER :: axis, n, i, j, na, nb, n1, n2, n3, m, p, k,l,q, r, i1, j1, na1
  REAL(DP) :: zeu_new(3, 3, nat)
  REAL(DP), ALLOCATABLE :: frc_new(:, :, :, :, :, :, :)
  type vector
     REAL(DP), pointer :: vec(:, :, :, :, :, :, :)
  end type vector
  !
  type (vector) u(6*3*nat)
  ! These are the "vectors" associated with the sum rules on force-constants
  !
  integer :: u_less(6*3*nat), n_less, i_less
  ! indices of the vectors u that are not independent to the preceding ones,
  ! n_less = number of such vectors, i_less = temporary parameter
  !
  integer, allocatable :: ind_v(:, :, :)
  REAL(DP), allocatable :: v(:, :)
  ! These are the "vectors" associated with symmetry conditions, coded by
  ! indicating the positions (i.e. the seven indices) of the non-zero elements (there
  ! should be only 2 of them) and the value of that element. We DO so in order
  ! to limit the amount of memory used.
  !
  REAL(DP), allocatable :: w(:, :, :, :, :, :, :), x(:, :, :, :, :, :, :)
  ! temporary vectors and parameters
  REAL(DP) :: scal, norm2, sum
  !
  REAL(DP) :: zeu_u(6 * 3, 3, 3, nat)
  ! These are the "vectors" associated with the sum rules on effective charges
  !
  integer :: zeu_less(6 * 3), nzeu_less, izeu_less
  ! indices of the vectors zeu_u that are not independent to the preceding ones,
  ! nzeu_less = number of such vectors, izeu_less = temporary parameter
  !
  REAL(DP) :: zeu_w(3, 3, nat), zeu_x(3, 3, nat)
  ! temporary vectors

  ! Initialization. n is the number of sum rules to be considered (IF asr.ne.'simple')
  ! and 'axis' is the rotation axis in the case of a 1D system
  ! (i.e. the rotation axis is (Ox) IF axis ='1', (Oy) if axis='2' and (Oz) if axis='3')
  !
  if((asr.ne.'simple').and.(asr.ne.'crystal').and.(asr.ne.'one-dim') &
                      .and.(asr.ne.'zero-dim')) THEN
     call errore('set_asr','invalid Acoustic Sum Rule:' // asr, 1)
  ENDIF
  !
  if(asr.eq.'simple') THEN
     !
     ! Simple Acoustic Sum Rule on effective charges
     !
     DO i = 1, 3
        DO j = 1, 3
           sum= 0.0d0
           DO na= 1, nat
              sum = sum + zeu(i, j, na)
           ENDDO
           DO na= 1, nat
              zeu(i, j, na) = zeu(i, j, na) - sum/nat
           ENDDO
        ENDDO
     ENDDO
     !
     ! Simple Acoustic Sum Rule on force constants in REAL space
     !
     DO i = 1, 3
        DO j = 1, 3
           DO na= 1, nat
              sum= 0.0d0
               DO nb= 1, nat
                  DO n1= 1, nr1
                     DO n2 = 1, nr2
                        DO n3= 1, nr3
                           sum=sum+frc(n1, n2, n3, i, j, na, nb)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
               frc(1, 1, 1, i, j, na, na) = frc(1, 1, 1, i, j, na, na) - sum
               !               WRITE(6,*) ' na, i, j, sum = ', na, i, j, sum
            ENDDO
         ENDDO
      ENDDO
      !
      return
      !
   ENDIF
  if(asr.eq.'crystal') n=3
  if(asr.eq.'one-dim') THEN
     ! the direction of periodicity is the rotation axis
     ! It will work only IF the crystal axis considered is one of
     ! the cartesian axis (typically, ibrav= 1, 6 or 8, or 4 along the
     ! z-direction)
     IF (nr1*nr2*nr3.eq.1) axis =3
     IF ((nr1.ne.1).and.(nr2*nr3.eq.1)) axis = 1
     IF ((nr2.ne.1).and.(nr1*nr3.eq.1)) axis =2
     IF ((nr3.ne.1).and.(nr1*nr2.eq.1)) axis =3
     IF (((nr1.ne.1).and.(nr2.ne.1)).or.((nr2.ne.1).and. &
          (nr3.ne.1)).or.((nr1.ne.1).and.(nr3.ne.1))) THEN
        call errore('set_asr','too many directions of &
             & periodicity in 1D system', axis)
     ENDIF
     IF ((ibrav.ne.1).and.(ibrav.ne.6).and.(ibrav.ne.8).and. &
          ((ibrav.ne.4).or.(axis.ne.3)) ) THEN
        WRITE(stdout,*) 'asr: rotational axis may be wrong'
     ENDIF
     WRITE(stdout,'("asr rotation axis in 1D system= ",I4)') axis
     n=4
  ENDIF
  if(asr.eq.'zero-dim') n=6
  !
  ! Acoustic Sum Rule on effective charges
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the effective charges matrix on
  !
  zeu_u(:, :, :, :) = 0.0d0
  DO i = 1, 3
     DO j = 1, 3
        DO na= 1, nat
           zeu_new(i, j, na) =zeu(i, j, na)
        ENDDO
     ENDDO
  ENDDO
  !
  p = 0
  DO i = 1, 3
     DO j = 1, 3
        ! These are the 3*3 vectors associated with the
        ! translational acoustic sum rules
        p = p + 1
        zeu_u(p, i, j,:) = 1.0d0
        !
     ENDDO
  ENDDO
  !
  IF (n.eq.4) THEN
     DO i = 1, 3
        ! These are the 3 vectors associated with the
        ! single rotational sum rule (1D system)
        p = p + 1
        DO na= 1, nat
           zeu_u(p, i, MOD(axis, 3) + 1, na) = -tau(MOD(axis + 1, 3) + 1, na)
           zeu_u(p, i, MOD(axis + 1, 3) + 1, na) = tau(MOD(axis, 3) + 1, na)
        ENDDO
        !
     ENDDO
  ENDIF
  !
  IF (n.eq.6) THEN
     DO i = 1, 3
        DO j = 1, 3
           ! These are the 3*3 vectors associated with the
           ! three rotational sum rules (0D system - typ. molecule)
           p = p + 1
           DO na= 1, nat
              zeu_u(p, i, MOD(j, 3) + 1, na) = -tau(MOD(j + 1, 3) + 1, na)
              zeu_u(p, i, MOD(j + 1, 3) + 1, na) = tau(MOD(j, 3) + 1, na)
           ENDDO
           !
        ENDDO
     ENDDO
  ENDIF
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  !
  nzeu_less = 0
  DO k = 1, p
     zeu_w(:, :, :) =zeu_u(k, :, :, :)
     zeu_x(:, :, :) =zeu_u(k, :, :, :)
     DO q= 1, k- 1
        r = 1
        DO izeu_less = 1, nzeu_less
           IF (zeu_less(izeu_less).eq.q) r = 0
        ENDDO
        IF (r.ne.0) THEN
           call sp_zeu(zeu_x, zeu_u(q, :, :, :), nat, scal)
           zeu_w(:, :, :) = zeu_w(:, :, :) - scal* zeu_u(q, :, :, :)
        ENDIF
     ENDDO
     call sp_zeu(zeu_w, zeu_w, nat, norm2)
     IF (norm2.gt.1.0d-16) THEN
        zeu_u(k, :, :, :) = zeu_w(:, :, :) / DSQRT(norm2)
     ELSE
        nzeu_less =nzeu_less+ 1
        zeu_less(nzeu_less) =k
     ENDIF
  ENDDO
  !
  ! Projection of the effective charge "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules
  !
  zeu_w(:, :, :) = 0.0d0
  DO k = 1, p
     r = 1
     DO izeu_less = 1, nzeu_less
        IF (zeu_less(izeu_less).eq.k) r = 0
     ENDDO
     IF (r.ne.0) THEN
        zeu_x(:, :, :) = zeu_u(k, :, :, :)
        call sp_zeu(zeu_x, zeu_new, nat, scal)
        zeu_w(:, :, :) = zeu_w(:, :, :) + scal * zeu_u(k, :, :, :)
     ENDIF
  ENDDO
  !
  ! Final substraction of the former projection to the initial zeu, to get
  ! the new "projected" zeu
  !
  zeu_new(:, :, :) = zeu_new(:, :, :) - zeu_w(:, :, :)
  call sp_zeu(zeu_w, zeu_w, nat, norm2)
  IF (ionode) WRITE(*,*)
  WRITE(stdout,'("Norm of the difference between old and new effective ", &
       & "charges: ", F25.20)') SQRT(norm2)
  !
  ! Check projection
  !
  !WRITE(6,'("Check projection of zeu")')
  !DO k = 1, p
  !  zeu_x(:, :, :) =zeu_u(k,:, :, :)
  !  call sp_zeu(zeu_x, zeu_new, nat, scal)
  !  IF (DABS(scal).gt.1d- 10) WRITE(6,'("k = ",I8," zeu_new|zeu_u(k) = ", F15.10)') k, scal
  !ENDDO
  !
  DO i = 1, 3
     DO j = 1, 3
        DO na= 1, nat
           zeu(i, j, na) = zeu_new(i, j, na)
        ENDDO
     ENDDO
  ENDDO
  !
  ! Acoustic Sum Rule on force constants
  !
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the force-constants matrix on
  !
  DO k = 1, 18*nat
     ALLOCATE(u(k) % vec(nr1, nr2, nr3, 3, 3, nat, nat))
     u(k) % vec (:, :, :, :, :, :, :) = 0.0d0
  ENDDO
  ALLOCATE (frc_new(nr1, nr2, nr3, 3, 3, nat, nat))
  DO i = 1, 3
     DO j = 1, 3
        DO na= 1, nat
           DO nb= 1, nat
              DO n1= 1, nr1
                 DO n2 = 1, nr2
                    DO n3= 1, nr3
                       frc_new(n1, n2, n3, i, j, na, nb) =frc(n1, n2, n3, i, j, na, nb)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  p = 0
  DO i = 1, 3
     DO j = 1, 3
        DO na= 1, nat
           ! These are the 3*3*nat vectors associated with the
           ! translational acoustic sum rules
           p = p + 1
           u(p) % vec (:, :, :, i, j, na,:) = 1.0d0
           !
        ENDDO
     ENDDO
  ENDDO
  !
  IF (n.eq.4) THEN
     DO i = 1, 3
        DO na= 1, nat
           ! These are the 3*nat vectors associated with the
           ! single rotational sum rule (1D system)
           p = p + 1
           DO nb= 1, nat
              u(p) % vec (:, :, :, i, MOD(axis, 3) + 1, na, nb) =-tau(MOD(axis+ 1, 3) + 1, nb)
              u(p) % vec (:, :, :, i, MOD(axis+ 1, 3) + 1, na, nb) =tau(MOD(axis, 3) + 1, nb)
           ENDDO
           !
        ENDDO
     ENDDO
  ENDIF
  !
  IF (n.eq.6) THEN
     DO i = 1, 3
        DO j = 1, 3
           DO na= 1, nat
              ! These are the 3*3*nat vectors associated with the
              ! three rotational sum rules (0D system - typ. molecule)
              p = p + 1
              DO nb= 1, nat
                 u(p) % vec (:, :, :, i, MOD(j, 3) + 1, na, nb) =-tau(MOD(j+ 1, 3) + 1, nb)
                 u(p) % vec (:, :, :, i, MOD(j+ 1, 3) + 1, na, nb) =tau(MOD(j, 3) + 1, nb)
              ENDDO
              !
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  !
  ALLOCATE (ind_v(9*nat*nat*nr1*nr2*nr3, 2,7), v(9*nat*nat*nr1*nr2*nr3, 2) )
  m= 0
  DO i = 1, 3
     DO j = 1, 3
        DO na= 1, nat
           DO nb= 1, nat
              DO n1= 1, nr1
                 DO n2 = 1, nr2
                    DO n3= 1, nr3
                       ! These are the vectors associated with the symmetry constraints
                       q= 1
                       l= 1
                       DO while((l.le.m).and.(q.ne.0))
                          IF ((ind_v(l, 1, 1).eq.n1).and.(ind_v(l, 1, 2).eq.n2).and. &
                               (ind_v(l, 1, 3).eq.n3).and.(ind_v(l, 1,4).eq.i).and. &
                               (ind_v(l, 1,5).eq.j).and.(ind_v(l, 1,6).eq.na).and. &
                               (ind_v(l, 1,7).eq.nb)) q= 0
                          IF ((ind_v(l, 2, 1).eq.n1).and.(ind_v(l, 2, 2).eq.n2).and. &
                               (ind_v(l, 2, 3).eq.n3).and.(ind_v(l, 2,4).eq.i).and. &
                               (ind_v(l, 2,5).eq.j).and.(ind_v(l, 2,6).eq.na).and. &
                               (ind_v(l, 2,7).eq.nb)) q= 0
                          l=l+ 1
                       ENDDO
                       IF ((n1.eq.MOD(nr1 + 1-n1, nr1) + 1).and.(n2.eq.MOD(nr2 + 1-n2, nr2) + 1) &
                            .and.(n3.eq.MOD(nr3 + 1-n3, nr3) + 1).and.(i.eq.j).and.(na.eq.nb)) q= 0
                       IF (q.ne.0) THEN
                          m=m+ 1
                          ind_v(m, 1, 1) =n1
                          ind_v(m, 1, 2) =n2
                          ind_v(m, 1, 3) =n3
                          ind_v(m, 1,4) =i
                          ind_v(m, 1,5) =j
                          ind_v(m, 1,6) =na
                          ind_v(m, 1,7) =nb
                          v(m, 1) = 1.0d0/DSQRT(2.0d0)
                          ind_v(m, 2, 1) =MOD(nr1 + 1-n1, nr1) + 1
                          ind_v(m, 2, 2) =MOD(nr2 + 1-n2, nr2) + 1
                          ind_v(m, 2, 3) =MOD(nr3 + 1-n3, nr3) + 1
                          ind_v(m, 2,4) =j
                          ind_v(m, 2,5) =i
                          ind_v(m, 2,6) =nb
                          ind_v(m, 2,7) =na
                          v(m, 2) =- 1.0d0/DSQRT(2.0d0)
                       ENDIF
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  ! Note that the vectors corresponding to symmetry constraints are already
  ! orthonormalized by construction.
  !
  n_less = 0
  ALLOCATE (w(nr1, nr2, nr3, 3, 3, nat, nat), x(nr1, nr2, nr3,3,3, nat, nat))
  DO k = 1, p
     w(:, :, :, :, :, :, :) =u(k) % vec (:, :, :, :, :, :, :)
     x(:, :, :, :, :, :, :) =u(k) % vec (:, :, :, :, :, :, :)
     DO l= 1,m
        !
        call sp2(x,v(l,:), ind_v(l,:, :), nr1, nr2, nr3, nat, scal)
        DO r = 1, 2
           n1=ind_v(l, r, 1)
           n2 =ind_v(l, r, 2)
           n3=ind_v(l, r, 3)
           i =ind_v(l, r,4)
           j =ind_v(l, r,5)
           na=ind_v(l, r,6)
           nb=ind_v(l, r,7)
           w(n1, n2, n3, i, j, na, nb) =w(n1, n2, n3, i, j, na, nb) -scal*v(l, r)
        ENDDO
     ENDDO
     IF (k.le.(9*nat)) THEN
        na1=MOD(k, nat)
        IF (na1.eq.0) na1=nat
        j1=MOD((k-na1)/nat, 3) + 1
        i1=MOD((((k-na1)/nat) -j1 + 1)/3, 3) + 1
     ELSE
        q=k-9*nat
        IF (n.eq.4) THEN
           na1=MOD(q, nat)
           IF (na1.eq.0) na1=nat
           i1=MOD((q-na1)/nat, 3) + 1
        ELSE
           na1=MOD(q, nat)
           IF (na1.eq.0) na1=nat
           j1=MOD((q-na1)/nat, 3) + 1
           i1=MOD((((q-na1)/nat) -j1 + 1)/3, 3) + 1
        ENDIF
     ENDIF
     DO q= 1, k- 1
        r = 1
        DO i_less = 1, n_less
           IF (u_less(i_less).eq.q) r = 0
        ENDDO
        IF (r.ne.0) THEN
           call sp3(x,u(q) % vec (:, :, :, :, :, :, :), i1, na1, nr1, nr2, nr3, nat, scal)
           w(:, :, :, :, :, :, :) = w(:, :, :, :, :, :, :) - scal* u(q) % vec (:, :, :, :, :, :, :)
        ENDIF
     ENDDO
     call sp1(w,w, nr1, nr2, nr3, nat, norm2)
     IF (norm2.gt.1.0d-16) THEN
        u(k) % vec (:, :, :, :, :, :, :) = w(:, :, :, :, :, :, :) / DSQRT(norm2)
     ELSE
        n_less =n_less+ 1
        u_less(n_less) =k
     ENDIF
  ENDDO
  !
  ! Projection of the force-constants "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules and symmetry contraints
  !
  w(:, :, :, :, :, :, :) = 0.0d0
  DO l= 1,m
     call sp2(frc_new,v(l,:), ind_v(l,:, :), nr1, nr2, nr3, nat, scal)
     DO r = 1, 2
        n1=ind_v(l, r, 1)
        n2 =ind_v(l, r, 2)
        n3=ind_v(l, r, 3)
        i =ind_v(l, r,4)
        j =ind_v(l, r,5)
        na=ind_v(l, r,6)
        nb=ind_v(l, r,7)
        w(n1, n2, n3, i, j, na, nb) =w(n1, n2, n3, i, j, na, nb) +scal*v(l, r)
     ENDDO
  ENDDO
  DO k = 1, p
     r = 1
     DO i_less = 1, n_less
        IF (u_less(i_less).eq.k) r = 0
     ENDDO
     IF (r.ne.0) THEN
        x(:, :, :, :, :, :, :) =u(k) % vec (:, :, :, :, :, :, :)
        call sp1(x,frc_new, nr1, nr2, nr3, nat, scal)
        w(:, :, :, :, :, :, :) = w(:, :, :, :, :, :, :) + scal*u(k)%vec(:, :, :, :, :, :, :)
     ENDIF
     DEALLOCATE(u(k) % vec)
  ENDDO
  !
  ! Final substraction of the former projection to the initial frc, to get
  ! the new "projected" frc
  !
  frc_new(:, :, :, :, :, :, :) =frc_new(:, :, :, :, :, :, :) - w(:, :, :, :, :, :, :)
  call sp1(w,w, nr1, nr2, nr3, nat, norm2)
  WRITE(stdout,'("Norm of the difference between old and new force-constants:",&
       &     F25.20)') SQRT(norm2)
  !
  ! Check projection
  !
  !WRITE(6,'("Check projection IFC")')
  !DO l= 1,m
  !  call sp2(frc_new,v(l,:), ind_v(l,:, :), nr1, nr2, nr3, nat, scal)
  !  IF (DABS(scal).gt.1d- 10) WRITE(6,'("l= ",I8," frc_new|v(l) = ", F15.10)') l, scal
  !ENDDO
  !DO k = 1, p
  !  x(:, :, :, :, :, :, :) =u(k) % vec (:, :, :, :, :, :, :)
  !  call sp1(x,frc_new, nr1, nr2, nr3, nat, scal)
  !  IF (DABS(scal).gt.1d- 10) WRITE(6,'("k = ",I8," frc_new|u(k) = ", F15.10)') k, scal
  !  DEALLOCATE(u(k) % vec)
  !ENDDO
  !
  DO i = 1, 3
     DO j = 1, 3
        DO na= 1, nat
           DO nb= 1, nat
              DO n1= 1, nr1
                 DO n2 = 1, nr2
                    DO n3= 1, nr3
                       frc(n1, n2, n3, i, j, na, nb) =frc_new(n1, n2, n3, i, j, na, nb)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  DEALLOCATE (x, w)
  DEALLOCATE (v, ind_v)
  DEALLOCATE (frc_new)
  !
  return
end subroutine set_asr
!
!----------------------------------------------------------------------
subroutine sp_zeu(zeu_u, zeu_v, nat, scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two effective charges matrices zeu_u and zeu_v
  ! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
  !
  USE kinds, ONLY: DP
  implicit none
  integer i, j, na, nat
  REAL(DP) zeu_u(3, 3, nat)
  REAL(DP) zeu_v(3, 3, nat)
  REAL(DP) scal
  !
  !
  scal= 0.0d0
  DO i = 1, 3
    DO j = 1, 3
      DO na= 1, nat
        scal=scal+zeu_u(i, j, na) *zeu_v(i, j, na)
      ENDDO
    ENDDO
  ENDDO
  !
  return
  !
end subroutine sp_zeu
!
!
!----------------------------------------------------------------------
subroutine sp1(u,v, nr1, nr2, nr3, nat, scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two force-constants matrices u and v (considered as
  ! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space, and coded in the usual way)
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1, nr2, nr3, i, j, na, nb, n1, n2, n3, nat
  REAL(DP) u(nr1, nr2, nr3, 3, 3, nat, nat)
  REAL(DP) v(nr1, nr2, nr3, 3, 3, nat, nat)
  REAL(DP) scal
  !
  !
  scal= 0.0d0
  DO i = 1, 3
    DO j = 1, 3
      DO na= 1, nat
        DO nb= 1, nat
          DO n1= 1, nr1
            DO n2 = 1, nr2
              DO n3= 1, nr3
                scal=scal+u(n1, n2, n3, i, j, na, nb) *v(n1, n2, n3, i, j, na, nb)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  return
  !
end subroutine sp1
!
!----------------------------------------------------------------------
subroutine sp2(u,v, ind_v, nr1, nr2, nr3, nat, scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two force-constants matrices u and v (considered as
  ! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space). u is coded in the usual way
  ! but v is coded as explained when defining the vectors corresponding to the
  ! symmetry constraints
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1, nr2, nr3, i, nat
  REAL(DP) u(nr1, nr2, nr3, 3, 3, nat, nat)
  integer ind_v(2,7)
  REAL(DP) v(2)
  REAL(DP) scal
  !
  !
  scal= 0.0d0
  DO i = 1, 2
    scal=scal+u(ind_v(i, 1), ind_v(i, 2), ind_v(i, 3), ind_v(i,4), ind_v(i,5), ind_v(i,6), &
         ind_v(i,7)) *v(i)
  ENDDO
  !
  return
  !
end subroutine sp2
!
!----------------------------------------------------------------------
subroutine sp3(u,v, i, na, nr1, nr2, nr3, nat, scal)
  !-----------------------------------------------------------------------
  !
  ! like sp1, but in the particular case when u is one of the u(k)%vec
  ! defined in set_asr (before orthonormalization). In this case most of the
  ! terms are zero (the ones that are not are characterized by i and na), so
  ! that a lot of computer time can be saved (during Gram-Schmidt).
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1, nr2, nr3, i, j, na, nb, n1, n2, n3, nat
  REAL(DP) u(nr1, nr2, nr3, 3, 3, nat, nat)
  REAL(DP) v(nr1, nr2, nr3, 3, 3, nat, nat)
  REAL(DP) scal
  !
  !
  scal= 0.0d0
  DO j = 1, 3
    DO nb= 1, nat
      DO n1= 1, nr1
        DO n2 = 1, nr2
          DO n3= 1, nr3
            scal=scal+u(n1, n2, n3, i, j, na, nb) *v(n1, n2, n3, i, j, na, nb)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  return
  !
end subroutine sp3
!
!-----------------------------------------------------------------------
SUBROUTINE q_gen(nsc,qbid, at_blk, bg_blk, at, bg)
  !-----------------------------------------------------------------------
  ! generate list of q (qbid) that are G-vectors of the supercell
  ! but not of the bulk
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  INTEGER :: nsc
  REAL(DP) qbid(3, nsc), at_blk(3, 3), bg_blk(3, 3), at(3,3), bg(3,3)
  !
  INTEGER, PARAMETER:: nr1=4, nr2 =4, nr3=4, &
                       nrm=(2*nr1 + 1) * (2*nr2 + 1) * (2*nr3 + 1)
  REAL(DP), PARAMETER:: eps = 1.0d-7
  INTEGER :: i, j, k, i1, i2, i3, idum(nrm), iq
  REAL(DP) :: qnorm(nrm), qbd(3, nrm) ,qwork(3), delta
  LOGICAL lbho
  !
  i = 0
  DO i1=-nr1, nr1
     DO i2 =-nr2, nr2
        DO i3=-nr3, nr3
           i = i + 1
           DO j = 1, 3
              qwork(j) = i1*bg(j, 1) + i2*bg(j, 2) + i3*bg(j, 3)
           ENDDO ! j
           !
           qnorm(i)  = qwork(1)**2 + qwork(2)**2 + qwork(3) **2
           !
           DO j = 1, 3
              !
              qbd(j, i) = at_blk(1, j) *qwork(1) + &
                         at_blk(2, j) *qwork(2) + &
                         at_blk(3, j) *qwork(3)
           ENDDO ! j
           !
           idum(i) = 1
           !
        ENDDO ! i3
     ENDDO ! i2
  ENDDO ! i1
  !
  DO i = 1, nrm- 1
     IF (idum(i).EQ.1) THEN
        DO j =i + 1, nrm
           IF (idum(j).EQ.1) THEN
              lbho=.TRUE.
              DO k = 1, 3
                 delta = qbd(k, i) -qbd(k, j)
                 lbho = lbho.AND. (ABS(NINT(delta) -delta)<eps)
              ENDDO ! k
              IF (lbho) THEN
                 IF(qnorm(i).GT.qnorm(j)) THEN
                    qbd(1, i) = qbd(1, j)
                    qbd(2, i) = qbd(2, j)
                    qbd(3, i) = qbd(3, j)
                    qnorm(i) = qnorm(j)
                 ENDIF
                 idum(j) = 0
              ENDIF
           ENDIF
        ENDDO ! j
     ENDIF
  ENDDO ! i
  !
  iq = 0
  DO i = 1, nrm
     IF (idum(i).EQ.1) THEN
        iq=iq+ 1
        qbid(1, iq) = bg_blk(1, 1) *qbd(1, i) +  &
                    bg_blk(1, 2) *qbd(2, i) +  &
                    bg_blk(1, 3) *qbd(3, i)
        qbid(2, iq) = bg_blk(2, 1) *qbd(1, i) +  &
                    bg_blk(2, 2) *qbd(2, i) +  &
                    bg_blk(2, 3) *qbd(3, i)
        qbid(3, iq) = bg_blk(3, 1) *qbd(1, i) +  &
                    bg_blk(3, 2) *qbd(2, i) +  &
                    bg_blk(3, 3) *qbd(3, i)
     ENDIF
  ENDDO ! i
  !
  IF (iq.NE.nsc) CALL errore('q_gen',' probably nr1, nr2, nr3 too small ', iq)
  RETURN
END SUBROUTINE q_gen
!
!-----------------------------------------------------------------------
SUBROUTINE check_at(at, bg_blk, alat, omega)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP) :: at(3, 3), bg_blk(3, 3), alat, omega
  REAL(DP) :: work(3, 3)
  INTEGER :: i, j
  REAL(DP), PARAMETER :: small= 1.d-6
  !
  work(:, :) = at(:, :)
  CALL cryst_to_cart(3,work, bg_blk,-1)
  !
  DO j = 1, 3
     DO i = 1, 3
        IF ( ABS(work(i, j) -NINT(work(i, j))) > small) THEN
           WRITE (stdout,'(3f9.4)') work(:, :)
           CALL errore ('check_at','at not multiple of at_blk', 1)
        ENDIF
     ENDDO
  ENDDO
  !
  omega =alat**3 * ABS(at(1, 1) * (at(2, 2) *at(3, 3) -at(3, 2) *at(2, 3)) - &
                       at(1, 2) * (at(2, 1) *at(3, 3) -at(2, 3) *at(3, 1)) + &
                       at(1, 3) * (at(2, 1) *at(3, 2) -at(2, 2) *at(3, 1)))
  !
  RETURN
END SUBROUTINE check_at
!
!-----------------------------------------------------------------------
SUBROUTINE set_tau (nat, nat_blk, at, at_blk, tau, tau_blk, &
     ityp, ityp_blk, itau_blk)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  INTEGER nat, nat_blk, ityp(nat), ityp_blk(nat_blk), itau_blk(nat)
  REAL(DP) at(3, 3), at_blk(3, 3),tau(3, nat),tau_blk(3, nat_blk)
  !
  REAL(DP) bg(3, 3), r(3) ! work vectors
  INTEGER i, i1, i2, i3, na, na_blk
  REAL(DP) small
  INTEGER NN1,NN2,NN3
  PARAMETER (NN1=8, NN2 =8, NN3=8, small= 1.d-8)
  !
  CALL recips (at(1, 1), at(1, 2), at(1, 3), bg(1, 1), bg(1, 2), bg(1, 3))
  !
  na = 0
  !
  DO i1 = -NN1,NN1
     DO i2 = -NN2,NN2
        DO i3 = -NN3,NN3
           r(1) = i1*at_blk(1, 1) + i2*at_blk(1, 2) + i3*at_blk(1, 3)
           r(2) = i1*at_blk(2, 1) + i2*at_blk(2, 2) + i3*at_blk(2, 3)
           r(3) = i1*at_blk(3, 1) + i2*at_blk(3, 2) + i3*at_blk(3, 3)
           CALL cryst_to_cart(1, r, bg,- 1)
           !
           IF ( r(1).GT.-small .AND. r(1)<1.d0-small .AND.          &
                r(2).GT.-small .AND. r(2)<1.d0-small .AND.          &
                r(3).GT.-small .AND. r(3)<1.d0-small ) THEN
              CALL cryst_to_cart(1, r, at,+ 1)
              !
              DO na_blk = 1, nat_blk
                 na = na + 1
                 IF (na.GT.nat) CALL errore('set_tau','too many atoms', na)
                 tau(1, na)    = tau_blk(1, na_blk) + r(1)
                 tau(2, na)    = tau_blk(2, na_blk) + r(2)
                 tau(3, na)    = tau_blk(3, na_blk) + r(3)
                 ityp(na)     = ityp_blk(na_blk)
                 itau_blk(na) = na_blk
              ENDDO
              !
           ENDIF
           !
        ENDDO
     ENDDO
  ENDDO
  !
  IF (na.NE.nat) CALL errore('set_tau','too few atoms: increase NNs', na)
  !
  RETURN
END SUBROUTINE set_tau
!
!-----------------------------------------------------------------------
SUBROUTINE read_tau &
     (nat, nat_blk, ntyp, bg_blk, tau, tau_blk, ityp, itau_blk)
  !---------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : ionode_id, ionode
  USE mp,         ONLY : mp_bcast
  USE mp_world,   ONLY : world_comm
  !
  IMPLICIT NONE
  !
  INTEGER nat, nat_blk, ntyp, ityp(nat), itau_blk(nat)
  REAL(DP) bg_blk(3, 3),tau(3, nat),tau_blk(3, nat_blk)
  !
  REAL(DP) r(3) ! work vectors
  INTEGER i, na, na_blk
  !
  REAL(DP) small
  PARAMETER ( small = 1.d-6 )
  !
  DO na= 1, nat
     IF (ionode) READ(5,*) (tau(i, na), i = 1, 3), ityp(na)
     CALL mp_bcast(tau(:, na), ionode_id, world_comm)
     CALL mp_bcast(ityp(na), ionode_id, world_comm)
     IF (ityp(na).LE.0 .OR. ityp(na) > ntyp) &
          CALL errore('read_tau',' wrong atomic type', na)
     DO na_blk = 1, nat_blk
        r(1) = tau(1, na) - tau_blk(1, na_blk)
        r(2) = tau(2, na) - tau_blk(2, na_blk)
        r(3) = tau(3, na) - tau_blk(3, na_blk)
        CALL cryst_to_cart(1, r, bg_blk,- 1)
        IF (ABS( r(1) -NINT(r(1)) ) < small .AND.                 &
            ABS( r(2) -NINT(r(2)) ) < small .AND.                 &
            ABS( r(3) -NINT(r(3)) ) < small ) THEN
           itau_blk(na) = na_blk
           go to 999
        ENDIF
     ENDDO
     CALL errore ('read_tau',' wrong atomic position ', na)
999  CONTINUE
  ENDDO
  !
  RETURN
END SUBROUTINE read_tau
!
!-----------------------------------------------------------------------
!
!---------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine setgam (q, gam, nat, at, bg,tau, itau_blk, nsc, alat, &
     &             gam_blk, nat_blk, at_blk, bg_blk,tau_blk, omega_blk, &
     &             frcg, nr1, nr2, nr3, rws, nrws, fd)
  !-----------------------------------------------------------------------
  ! compute the dynamical matrix (the analytic part only)
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : tpi
  implicit none
  !
  ! I/O variables
  !
  integer        :: nr1, nr2, nr3, nat, nat_blk,  &
                    nsc, nrws, itau_blk(nat)
  REAL(DP)       :: q(3), tau(3, nat), at(3, 3), bg(3, 3), alat, rws(0:3, nrws)
  REAL(DP)       :: tau_blk(3, nat_blk), at_blk(3, 3), bg_blk(3, 3), omega_blk, &
                    frcg(nr1, nr2, nr3, 3, 3, nat_blk, nat_blk)
  COMPLEX(DP)    :: gam_blk(3, 3, nat_blk, nat_blk),f_of_q(3, 3, nat, nat)
  COMPLEX(DP)    ::  gam(3, 3, nat, nat)
  LOGICAL        :: fd
  !
  ! local variables
  !
  REAL(DP)        :: arg
  complex(DP)     :: cfac(nat)
  integer         :: i, j, k, na, nb, na_blk, nb_blk, iq
  REAL(DP)        :: qp(3), qbid(3, nsc) ! automatic array
  !
  !
  call q_gen(nsc,qbid, at_blk, bg_blk, at, bg)
  !
  f_of_q=(0.0_DP,0.0_DP)
  DO iq= 1, nsc
     !
     DO k = 1, 3
        qp(k) = q(k) + qbid(k, iq)
     ENDDO
     !
     gam_blk(:, :, :, :) = (0.d0,0.d0)
     CALL frc_blk (gam_blk, qp,tau_blk, nat_blk,              &
                   nr1, nr2, nr3,frcg, at_blk, bg_blk, rws, nrws,f_of_q, fd)
     !
     DO na= 1, nat
        na_blk = itau_blk(na)
        DO nb= 1, nat
           nb_blk = itau_blk(nb)
           !
           arg = tpi * ( qp(1) * ( (tau(1, na) -tau_blk(1, na_blk)) -   &
                                (tau(1, nb) -tau_blk(1, nb_blk)) ) + &
                      qp(2) * ( (tau(2, na) -tau_blk(2, na_blk)) -   &
                                (tau(2, nb) -tau_blk(2, nb_blk)) ) + &
                      qp(3) * ( (tau(3, na) -tau_blk(3, na_blk)) -   &
                                (tau(3, nb) -tau_blk(3, nb_blk)) ) )
           !
           cfac(nb) = CMPLX(cos(arg), sin(arg), kind=dp)/nsc
           !
        ENDDO ! nb
        DO nb= 1, nat
           DO i = 1, 3
              DO j = 1, 3
                 nb_blk = itau_blk(nb)
                 gam(i, j, na, nb) = gam(i, j, na, nb) + cfac(nb) * &
                     gam_blk(i, j, na_blk, nb_blk)
              ENDDO ! j
              ENDDO ! i
        ENDDO ! nb
     ENDDO ! na
     !
  ENDDO ! iq
  !
  return
end subroutine setgam
!
!--------------------------------------------------------------------
function dos_gam (nbndx, nq, jbnd, gamma, et, ef)
  !--------------------------------------------------------------------
  ! calculates weights with the tetrahedron method (Bloechl version)
  ! this subroutine is based on tweights.f90 belonging to PW
  ! it calculates a2F on the surface of given frequency <=> histogram
  ! Band index means the frequency mode here
  ! and "et" means the frequency(mode,q-point)
  !
  USE kinds,       ONLY: DP
  USE parameters
  USE ktetra, ONLY : ntetra, tetra
  implicit none
  !
  integer :: nq, nbndx, jbnd
  REAL(DP) :: et(nbndx, nq), gamma(nbndx, nq), func

  REAL(DP) :: ef
  REAL(DP) :: e1, e2, e3, e4, c1, c2, c3, c4, etetra(4)
  integer      :: ik, ibnd, nt, nk, ns, i, ik1, ik2, ik3, ik4, itetra(4)

  REAL(DP) ::   f12,f13,f14,f23,f24,f34, f21,f31,f41,f42,f32,f43
  REAL(DP) ::   P1,P2,P3,P4, G, o13, Y1,Y2,Y3,Y4, eps,vol, Tint
  REAL(DP) :: dos_gam

  Tint = 0.0d0
  o13 = 1.0_dp/3.0_dp
  eps  = 1.0d-14
  vol  = 1.0d0/ntetra
  P1 = 0.0_dp
  P2 = 0.0_dp
  P3 = 0.0_dp
  P4 = 0.0_dp
  DO nt = 1, ntetra
     ibnd = jbnd
     !
     ! etetra are the energies at the vertexes of the nt-th tetrahedron
     !
     DO i = 1, 4
        etetra(i) = et(ibnd, tetra(i, nt))
     ENDDO
     itetra(1) = 0
     call hpsort (4, etetra, itetra)
     !
     ! ...sort in ascending order: e1 < e2 < e3 < e4
     !
     e1 = etetra (1)
     e2 = etetra (2)
     e3 = etetra (3)
     e4 = etetra (4)
     !
     ! kp1-kp4 are the irreducible k-points corresponding to e1- e4
     !
     ik1 = tetra(itetra(1), nt)
     ik2 = tetra(itetra(2), nt)
     ik3 = tetra(itetra(3), nt)
     ik4 = tetra(itetra(4), nt)
     Y1  = gamma(ibnd, ik1)/et(ibnd, ik1)
     Y2  = gamma(ibnd, ik2)/et(ibnd, ik2)
     Y3  = gamma(ibnd, ik3)/et(ibnd, ik3)
     Y4  = gamma(ibnd, ik4)/et(ibnd, ik4)

     IF ( e3 < ef .and. ef < e4) THEN

        f14 = (ef- e4)/(e1- e4)
        f24 = (ef- e4)/(e2- e4)
        f34 = (ef- e4)/(e3- e4)

        G  =  3.0_dp * f14 * f24 * f34 / (e4- ef)
        P1 =  f14 * o13
        P2 =  f24 * o13
        P3 =  f34 * o13
        P4 =  (3.0_dp - f14 - f24 - f34 ) * o13

     ELSE IF ( e2 < ef .and. ef < e3 ) THEN

        f13 = (ef- e3)/(e1- e3)
        f31 = 1.0_dp - f13
        f14 = (ef- e4)/(e1- e4)
        f41 = 1.0_dp-f14
        f23 = (ef- e3)/(e2- e3)
        f32 = 1.0_dp - f23
        f24 = (ef- e4)/(e2- e4)
        f42 = 1.0_dp - f24

        G   =  3.0_dp * (f23*f31 + f32*f24)
        P1  =  f14 * o13 + f13*f31*f23 / G
        P2  =  f23 * o13 + f24*f24*f32 / G
        P3  =  f32 * o13 + f31*f31*f23 / G
        P4  =  f41 * o13 + f42*f24*f32 / G
        G   =  G / (e4- e1)

     ELSE IF ( e1 < ef .and. ef < e2 ) THEN

        f12 = (ef- e2)/(e1- e2)
        f21 = 1.0_dp - f12
        f13 = (ef- e3)/(e1- e3)
        f31 = 1.0_dp - f13
        f14 = (ef- e4)/(e1- e4)
        f41 = 1.0_dp - f14

        G  =  3.0_dp * f21 * f31 * f41 / (ef- e1)
        P1 =  o13 * (f12 + f13 + f14)
        P2 =  o13 * f21
        P3 =  o13 * f31
        P4 =  o13 * f41

     ELSE

        G = 0.0_dp

     ENDIF

     Tint = Tint + G * (Y1*P1 + Y2*P2 + Y3*P3 + Y4*P4) * vol

  ENDDO   ! ntetra


  dos_gam = Tint  !2 because DOS_ee is per 1 spin

  return
end function dos_gam
!
!
!-----------------------------------------------------------------------
subroutine readfg ( ifn, nr1, nr2, nr3, nat, frcg )
  !-----------------------------------------------------------------------
  !
  USE kinds,       ONLY : DP
  USE io_global,   ONLY : ionode, ionode_id, stdout
  USE mp,          ONLY : mp_bcast
  USE mp_world,    ONLY : world_comm
  implicit none
  ! I/O variable
  integer, INTENT(in) ::  nr1, nr2, nr3, nat
  REAL(DP), INTENT(out) :: frcg(nr1, nr2, nr3, 3, 3, nat, nat)
  ! local variables
  integer i, j, na, nb, m1,m2,m3, ifn
  integer ibid, jbid, nabid, nbbid, m1bid,m2bid,m3bid
  !
  !
  IF (ionode) READ (ifn,*) m1, m2, m3
  CALL mp_bcast(m1, ionode_id, world_comm)
  CALL mp_bcast(m2, ionode_id, world_comm)
  CALL mp_bcast(m3, ionode_id, world_comm)
  IF ( m1 /= nr1 .or. m2 /= nr2 .or. m3 /= nr3) &
       call errore('readfg','inconsistent nr1, nr2, nr3 read', 1)
  DO i = 1, 3
     DO j = 1, 3
        DO na= 1, nat
           DO nb= 1, nat
              IF (ionode) read (ifn, *) ibid, jbid, nabid, nbbid
              CALL mp_bcast(ibid, ionode_id, world_comm)
              CALL mp_bcast(jbid, ionode_id, world_comm)
              CALL mp_bcast(nabid, ionode_id, world_comm)
              CALL mp_bcast(nbbid, ionode_id, world_comm)
              
              if(i.ne.ibid.or.j.ne.jbid.or.na.ne.nabid.or.nb.ne.nbbid)  THEN
                  WRITE(stdout,*) i, j, na, nb,'  <>  ', ibid, jbid, nabid, nbbid
                  call errore  ('readfG','error in reading', 1)
              ELSE
                  IF (ionode) read (ifn,*) (((m1bid, m2bid, m3bid,     &
                                 frcg(m1,m2,m3, i, j, na, nb), &
                                 m1= 1, nr1),m2 =1, nr2),m3=1, nr3)
              ENDIF
              CALL mp_bcast(frcg(:, :, :, i, j, na, nb), ionode_id, world_comm)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  IF (ionode) CLOSE(ifn)
  !
  return
end subroutine readfg
!
!
SUBROUTINE find_representations_mode_q ( nat, ntyp, xq, w2, u, tau, ityp, &
                  amass, num_rap_mode, nspin_mag )

  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE symm_base,  ONLY : s, sr, ft, irt, nsym, nrot, t_rev, time_reversal,&
                         sname, copy_sym, s_axis_to_cart

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat, ntyp, nspin_mag
  REAL(DP), INTENT(IN) :: xq(3), amass(ntyp), tau(3, nat)
  REAL(DP), INTENT(IN) :: w2(3*nat)
  INTEGER, INTENT(IN) :: ityp(nat)
  COMPLEX(DP), INTENT(IN) :: u(3*nat, 3*nat)
  INTEGER, INTENT(OUT) :: num_rap_mode(3*nat)
  REAL(DP) :: gi (3, 48), gimq (3), sr_is(3, 3,48), rtau(3,48, nat)
  INTEGER :: irotmq, nsymq, nsym_is, isym, i, ierr
  LOGICAL :: minus_q, search_sym, sym(48), magnetic_sym
!
!  find the small group of q
!
  time_reversal=.TRUE.
  IF (.NOT.time_reversal) minus_q=.FALSE.

  sym(1:nsym) =.TRUE.
  call smallg_q (xq, 0, at, bg, nsym, s, ft, sym, minus_q)
  nsymq= copy_sym(nsym, sym )
  call s_axis_to_cart ()
  CALL set_giq (xq, s, nsymq, nsym, irotmq,minus_q,gi,gimq)
!
!  IF the small group of q is non symmorphic,
!  search the symmetries only IF there are no G such that Sq -> q+G
!
  search_sym=.TRUE.
   IF ( ANY ( ABS(ft(:,1:nsymq)) > 1.0d-8 ) ) THEN
     DO isym= 1, nsymq
        search_sym=( search_sym.and.(ABS(gi(1, isym))<1.d-8).and.  &
                                    (ABS(gi(2, isym))<1.d-8).and.  &
                                    (ABS(gi(3, isym))<1.d-8) )
     ENDDO
  ENDIF
!
!  Set the representations tables of the small group of q and
!  find the mode symmetry
!
  IF (search_sym) THEN
     magnetic_sym=(nspin_mag==4)
     CALL prepare_sym_analysis(nsymq, sr,t_rev,magnetic_sym)
     sym (1:nsym) = .TRUE.
     CALL sgam_lr (at, bg, nsym, s, irt, tau, rtau, nat)
     CALL find_mode_sym_new (u, w2, tau, nat, nsymq, s, sr, irt, xq,    &
             rtau, amass, ntyp, ityp, 1, .FALSE., .FALSE., num_rap_mode, ierr)

  ENDIF
  RETURN
  END SUBROUTINE find_representations_mode_q
!
SUBROUTINE  qpoint_gen1(dim1, dim2, dim3, ctrAB) 
!
  use kinds, only: dp
  USE mp,   ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_global,  ONLY : inter_pool_comm
  
  IMPLICIT NONE
  ! input
  INTEGER, INTENT(in)             :: dim1, dim2, dim3
  INTEGER, INTENT(out)            :: ctrAB
!!  REAL(DP), INTENT(out)           :: q_AB(:,:)
  ! local
  INTEGER                         :: i, j, k, n, nqs
  INTEGER                         :: lower_bnd, upper_bnd
  REAL(DP), ALLOCATABLE           :: q_all(:, :)
  REAL(DP)                        :: q_B(3), q_A(3), eps
  !
  nqs = dim1 * dim2 * dim3  
  eps = 1.0E-06
  !
  ALLOCATE(q_all(3, nqs))
  q_all = 0.d0
  !
  CALL fkbounds( dim1, lower_bnd, upper_bnd )
!  DO i = 1, dim1
  DO i = lower_bnd, upper_bnd
      DO j = 1, dim2
         DO k = 1, dim3
            !  this is nothing but consecutive ordering
            n = (k - 1) + (j - 1) * dim3 + (i - 1) * dim2 * dim3 + 1
            !  q_all are the components of the complete grid in crystal axis
            q_all(1, n) = dble(i - 1) / dim1 ! + dble(k1)/2/dim1
            q_all(2, n) = dble(j - 1) / dim2 ! + dble(k2)/2/dim2
            q_all(3, n) = dble(k - 1) / dim3 ! + dble(k3)/2/dim3 ! k1 , k2 , k3 is for the shift
         ENDDO
      ENDDO
  ENDDO
  CALL mp_sum(q_all, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
  ctrAB = 0
  CALL fkbounds( nqs, lower_bnd, upper_bnd )
     !
  DO i = lower_bnd, upper_bnd !1, nqs
    q_A = q_all(:, i) + q_all(:, i) ! q_A to find if q belongs in A 
    IF (((ABS(q_A(1)) .LT. eps) .OR. (abs(abs(q_A(1)) - 1) .LT. eps)) .AND. &
        ((ABS(q_A(2)) .LT. eps) .OR. (abs(abs(q_A(2)) - 1) .LT. eps)) .AND. &
        ((ABS(q_A(3)) .LT. eps) .OR. (abs(abs(q_A(3)) - 1) .LT. eps))) THEN
        ctrAB = ctrAB + 1
    ELSE
     DO j = i + 1, nqs
        q_B = q_all(:, i) + q_all(:, j)       
       IF (((ABS(q_B(1)) .LT. eps) .OR. (abs(abs(q_B(1)) - 1) .LT. eps)) .AND. &
           ((ABS(q_B(2)) .LT. eps) .OR. (abs(abs(q_B(2)) - 1) .LT. eps)) .AND. &
           ((ABS(q_B(3)) .LT. eps) .OR. (abs(abs(q_B(3)) - 1) .LT. eps))) THEN
           ctrAB = ctrAB + 1
       END IF 
      END DO
    END IF
  END DO
  CALL mp_sum(ctrAB, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
  DEALLOCATE(q_all)
  ! 
  RETURN
  !
END SUBROUTINE qpoint_gen1


SUBROUTINE  qpoint_gen2(dim1, dim2, dim3, ctrAB, q_AB) 
!
  use kinds, only: dp
  USE mp,         ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_global,  ONLY : inter_pool_comm
  
  IMPLICIT NONE
  ! input
  INTEGER, INTENT(in)             :: dim1, dim2, dim3, ctrAB
  REAL(DP), INTENT(out)           :: q_AB(3, ctrAB)
  ! local
  INTEGER                         :: i, j, k, n, ctr, nqs
  INTEGER                         :: lower_bnd, upper_bnd
  REAL(DP), ALLOCATABLE           :: q_all(:, :), q_AB_TMP(:, :)
  REAL(DP)                        :: q_B(3), q_A(3), eps
  !
  nqs = dim1 * dim2 * dim3  
  eps = 1.0E-06
  !
  ALLOCATE(q_all(3, nqs), q_AB_TMP(3, nqs))
  ! 
  q_all = 0.d0
  !
  CALL fkbounds( dim1, lower_bnd, upper_bnd )
  DO i = lower_bnd, upper_bnd
  !DO i = 1, dim1
      DO j = 1, dim2
         DO k = 1, dim3
            !  this is nothing but consecutive ordering
            n = (k - 1) + (j - 1) * dim3 + (i - 1) * dim2 * dim3 + 1
            !  q_all are the components of the complete grid in crystal axis
            q_all(1, n) = dble(i - 1) / dim1 ! + dble(k1)/2/dim1
            q_all(2, n) = dble(j - 1) / dim2 ! + dble(k2)/2/dim2
            q_all(3, n) = dble(k - 1) / dim3 ! + dble(k3)/2/dim3 ! k1 , k2 , k3 is for the shift
         ENDDO
      ENDDO
  ENDDO
  CALL mp_sum(q_all, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
  CALL fkbounds( nqs, lower_bnd, upper_bnd )
     !
  q_AB_TMP = 0.d0
  DO i = lower_bnd, upper_bnd !1, nqs
    !DO i = 1, nqs
    q_A = q_all(:, i) + q_all(:, i) ! q_A to find if q belongs in A 
    IF (((ABS(q_A(1)) .LT. eps) .OR. (abs(abs(q_A(1)) - 1) .LT. eps)) .AND. &
        ((ABS(q_A(2)) .LT. eps) .OR. (abs(abs(q_A(2)) - 1) .LT. eps)) .AND. &
        ((ABS(q_A(3)) .LT. eps) .OR. (abs(abs(q_A(3)) - 1) .LT. eps))) THEN
        q_AB_TMP(:, i) = q_all(:, i)
  !      write(*,*) "A", q_AB(:, ctr)
    ELSE
     DO j = i + 1, nqs
        q_B = q_all(:, i) + q_all(:, j)       
       IF (((ABS(q_B(1)) .LT. eps) .OR. (abs(abs(q_B(1)) - 1) .LT. eps)) .AND. &
           ((ABS(q_B(2)) .LT. eps) .OR. (abs(abs(q_B(2)) - 1) .LT. eps)) .AND. &
           ((ABS(q_B(3)) .LT. eps) .OR. (abs(abs(q_B(3)) - 1) .LT. eps))) THEN
           q_AB_TMP(:, i) = q_all(:, i)
  !         write(*,*) q_AB(:, ctr)
       END IF 
      END DO
    END IF
  ! 
  END DO
  CALL mp_sum(q_AB_TMP, inter_pool_comm)
  CALL mp_barrier(inter_pool_comm)
  !
  ctr = 1 ! so that Gamma is the first entry
  q_AB = 0.d0
  DO i = 1, nqs
    IF ((SUM(ABS(q_AB_TMP(:, i)))) .GT. eps ) THEN
        ctr = ctr + 1
        q_AB(:, ctr) = q_AB_TMP(:, i)
    ENDIF
  ENDDO
  !
  DEALLOCATE(q_all, q_AB_TMP)
  !
  RETURN
  !
END SUBROUTINE qpoint_gen2
!!
SUBROUTINE phonon_unfolding(nq, tau, nat, ntyp, amass, ityp, & 
                            dim1, dim2, dim3, ng1, ng2, ng3, q_nq, z_nq, &
                            flfrq, flweights, w2)
! author: Marios Zacharias
! version: v0.1
! In this subroutine phonon unfolding technique is implemented as described in
! Refs. [Phys. Rev. B 87, 085322 (2013)] and [Comput. Mat. Sci. 125, 218-223, (2016)].
!
  USE cell_base,  ONLY : bg, at
  USE io_global,  ONLY : ionode, ionode_id, stdout
  USE kinds,      ONLY : dp
  USE mp_world,   ONLY : world_comm 
  USE mp_global,  ONLY : inter_pool_comm
  USE mp,         ONLY : mp_bcast, mp_barrier, mp_sum
  USE constants,  ONLY : amu_ry, RY_TO_CMM1, tpi
!
  !
  IMPLICIT NONE
  ! input
  CHARACTER(LEN=256), INTENT(in) :: flfrq, flweights
  INTEGER, INTENT(in)            :: nq, nat, dim1, dim2, dim3, ng1, ng2, ng3 
  INTEGER, INTENT(in)            :: ntyp, ityp(nat)
  REAL(DP), INTENT(in)           :: q_nq(3, nq), w2(3 * nat, nq), amass(ntyp), tau(3, nat)
  COMPLEX(DP), INTENT(in)        :: z_nq(3 * nat, 3 * nat, nq)
  ! local
  INTEGER                        :: lower_bnd, upper_bnd ! For parallelization 
  INTEGER                        :: d1, d2, d3, ctr, ipol, n, i, na
  REAL(DP)                       :: dotp, pathL, q(3, nq) 
  REAL(DP)                       :: freq(3 * nat, nq), abc(3, nat), spctrl_weight2(3 * nat,nq)
  COMPLEX(DP)                    :: spctrl_weight(3 * nat,nq,3)
  COMPLEX(DP)                    :: imagi  
  !
  imagi = (0.0d0, 1.0d0) !imaginary unit
  abc   = tau
  q     = q_nq
  ! to convert tau/abc in crystal coordinates
  !
  IF (ionode) WRITE(stdout, *) "=============================================="
  IF (ionode) WRITE(stdout, *) "Performing phonon unfolding ..."
  IF (ionode) WRITE(stdout, *) "=============================================="
  !
  CALL cryst_to_cart(nat, abc, bg, -1)
  ! Assuming that we have a one to one correspondence for Q to q
  CALL cryst_to_cart(nq, q, at, -1) ! convert to crystal
  !
  spctrl_weight2(:, :) = 0.d0
  CALL fkbounds( nq, lower_bnd, upper_bnd )
  !
  DO n = lower_bnd, upper_bnd ! 1, nq 
  ! DO n = 1, nq
    DO i = 1, 3 * nat ! modes nu 
     ctr = 0
     DO d1 = -ng1 * dim1, ng1 * dim1, dim1
       DO d2 = -ng2 * dim2, ng2 * dim2, dim2
         DO d3 = -ng3 * dim3, ng3 * dim3, dim3
           spctrl_weight(i, n, :) = (0.d0, 0.d0)
             DO ipol = 1, 3
              DO na = 1, nat
              !
                dotp = ( q(1, n) + d1 ) * abc(1, na) &     !
                     + ( q(2, n) + d2 ) * abc(2, na) &     !
                     + ( q(3, n) + d3 ) * abc(3, na)  !    !
                spctrl_weight(i, n, ipol) = spctrl_weight(i, n, ipol) + z_nq((na - 1) * 3 + ipol, i, n) &
                                      * SQRT(amu_ry * amass(ityp(na))) / SQRT(DBLE(dim1 * dim2 * dim3)) & ! & 
                                      * EXP( - imagi * tpi * dotp)
              ENDDO ! na
              spctrl_weight2(i, n) = spctrl_weight2(i, n) + &
                                     DBLE(spctrl_weight(i, n, ipol) * CONJG(spctrl_weight(i, n, ipol)))
            ENDDO ! ipol
          ctr = ctr + 1
          ENDDO ! d3
        ENDDO ! d2
      ENDDO ! d1
    ENDDO ! i    
  ENDDO ! nq
  CALL mp_sum(spctrl_weight2, inter_pool_comm)
  CALL cryst_to_cart(nq, q, bg, +1) ! convert back to cartessian
!
  DO n=1,nq
    ! freq(i,n) = frequencies in cm^(-1), with negative sign if omega^2 is
    ! negative
    DO i=1,3*nat
       freq(i,n)= SQRT(ABS(w2(i, n))) * RY_TO_CMM1
       IF (w2(i, n) < 0.0d0) freq(i, n) = -freq(i,n)
    END DO
  END DO
  !
  IF(flfrq.NE.' '.and.ionode) THEN
        OPEN (unit=2,file=flfrq ,status='unknown',form='formatted')
        WRITE(2, '(" &plot nbnd=",i4,", nks=",i4," /")') 3*nat, nq
        DO n=1, nq
           WRITE(2, '(10x,3f10.6)')  q(1, n) / DBLE(dim1), q(2, n) / DBLE(dim2),q(3, n) / DBLE(dim3)
           WRITE(2,'(6f10.4)') (freq(i,n), i=1,3*nat)
        END DO
        CLOSE(unit=2)

        OPEN (unit=2,file=trim(flfrq)//'.gp' ,status='unknown',form='formatted')
        pathL = 0._dp
        WRITE(2, '(f10.6,3x,999f10.4)')  pathL,  (freq(i, 1), i = 1, 3 * nat)
        DO n=2, nq
           pathL=pathL+(SQRT(SUM(  (q(:,n)-q(:,n-1))**2 )))
           WRITE(2, '(f10.6,3x,999f10.4)')  pathL,  (freq(i, n), i = 1, 3 * nat)
        END DO
        CLOSE(unit=2)

  END IF
     ! Print weights 
  IF(flweights.NE.' '.and.ionode) THEN
        OPEN (unit=2,file=flweights ,status='unknown',form='formatted')
        WRITE(2, '(" &plot nbnd=",i4,", nks=",i4," /")') 3*nat, nq
        DO n=1, nq
           WRITE(2, '(10x,3f10.6)')  q(1, n) / DBLE(dim1), q(2, n) / DBLE(dim2),q(3, n) / DBLE(dim3)
           WRITE(2,'(6f10.4)') (spctrl_weight2(i, n) / ctr, i=1,3 * nat)
        END DO
        CLOSE(unit=2)

        OPEN (unit=2,file=trim(flweights)//'.gp', status='unknown',form='formatted')
        pathL = 0._dp
        WRITE(2, '(f10.6,3x,999f10.4)')  pathL,  (spctrl_weight2(i, 1) / ctr, i = 1, 3 * nat)
        DO n=2, nq
           pathL=pathL+(SQRT(SUM(  (q(:, n) - q(:, n - 1))**2 )))
           WRITE(2, '(f10.6,3x,999f10.4)')  pathL,  (spctrl_weight2(i, n) / ctr, i = 1, 3 * nat)
        END DO
        CLOSE(unit=2)
     !
  END IF 
  !
  RETURN
  !
END SUBROUTINE phonon_unfolding
!!
!!
SUBROUTINE ZG_configuration(nq, nat, ntyp, amass, ityp, q, w2, z_nq, & 
                      dim1, dim2, dim3, niters, error_thresh, synch, tau, alat, atm, &
                      ntypx, at, q_in_cryst_coord, q_external, T, incl_qA, & 
                      compute_error, single_ph_displ, &
                      ZG_strf, qpts_strf, atmsf_a, atmsf_b, &
                      nrots, kres1, kres2, kmin, kmax, col1, col2, Np, &
                      flscf, ASDM, apply_fd, read_fd_forces, incl_epsil, & 
                      iter_idx, iter_idx0, fd_displ, ibrav, epsil, zeu, & 
                      poly, poly_fd_forces, mixing, update_equil)
  !        
  USE kinds,      ONLY : dp
  USE constants,  ONLY : amu_ry, ry_to_thz, ry_to_cmm1, H_PLANCK_SI, &  
                         K_BOLTZMANN_SI, AMU_SI, pi, tpi, BOHR_RADIUS_ANGS, &
                         RYDBERG_SI
  USE cell_base,  ONLY : bg, celldm
  USE io_global,  ONLY : ionode, ionode_id, stdout
  USE mp_world,   ONLY : world_comm 
  USE mp_global,  ONLY : inter_pool_comm
  USE mp,         ONLY : mp_bcast, mp_barrier, mp_sum
  USE read_input, ONLY : read_input_file
  USE input_parameters, ONLY : ecutwfc, calculation, prefix, restart_mode, &
                               pseudo_dir, outdir, diagonalization, mixing_mode, &
                               mixing_beta, conv_thr, atom_pfile, &
                               nk1, nk2, nk3, k1, k2, k3, &
                               occupations, degauss, smearing 
  !
  IMPLICIT NONE
  ! input
  CHARACTER(LEN=3),   INTENT(in) :: atm(ntypx)
  CHARACTER(LEN=256), INTENT(in) :: flscf
  LOGICAL, INTENT(in)          :: synch, q_in_cryst_coord, q_external, ZG_strf
  LOGICAL, INTENT(in)          :: incl_qA, compute_error, single_ph_displ, ASDM, mixing
  LOGICAL, INTENT(in)          :: apply_fd, read_fd_forces, incl_epsil, poly, poly_fd_forces, update_equil
  INTEGER, INTENT(in)          :: dim1, dim2, dim3, niters, qpts_strf
  INTEGER, INTENT(in)          :: nq, nat, ntyp, ntypx, iter_idx, iter_idx0
  ! nq is the number of qpoints in sets A and B
  INTEGER, INTENT(in)          :: ityp(nat)
  INTEGER, INTENT(in)          :: nrots, kres1, kres2, col1, col2, Np
  INTEGER, INTENT(in)          :: ibrav ! for read_fd_forces in ASDM
  REAL(DP), INTENT(in)         :: kmin, kmax, fd_displ ! for apply_fd
  REAL(DP), INTENT(in)         :: error_thresh, alat, T
  REAL(DP), INTENT(in)         :: at(3, 3), atmsf_a(ntypx,5), atmsf_b(ntypx,5), zeu(3, 3, nat)
  ! zeu for read_fd_forces so we can print it
  REAL(DP), INTENT(in)         :: q(3, nq), w2(3 * nat, nq), amass(ntyp), tau(3, nat)
  COMPLEX(DP), INTENT(in)      :: z_nq(3 * nat, 3 * nat, nq)
  ! 
  CHARACTER(len=256)           :: filename, pt_T, pt_1, pt_2, pt_3
  !
  INTEGER                  :: nat3, na, nta, ipol, i, j, k, qp, ii, p, kk
  INTEGER                  :: nq_tot, pn, combs, combs_all, sum_zg
  INTEGER                  :: lower_bnd, upper_bnd
  INTEGER                  :: ctr, ctr2, ctrA, ctrB, ctrAB
  ! nq_tot total number of q-points (including sets A, B, C)
  ! pn combinations 
  INTEGER, ALLOCATABLE     :: Mx_mat(:, :), Mx_mat_or(:, :), M_mat(:, :), V_mat(:)
  INTEGER, ALLOCATABLE     :: Rlist(:, :)
  ! M matrices : sign matrices 
  !
  REAL(DP)                 :: freq(3 * nat, nq), ph_w(3 * nat, nq), l_q(3 * nat, nq)
  REAL(DP)                 :: q_A(3), q_B(3), p_q(3 * nat, nq), FE_q(3 * nat, nq), epsil(3, 3)
  ! epsil for read_fd_forces
  REAL(DP)                 :: hbar, ang, u_rand, dotp, PE_nq, KE_nq, FE_nq
  !
  REAL(DP), PARAMETER      :: eps = 1.0d-6, eps2 = 1.0d-15
  ! l_q --> amplitude \sigma at temperature T
  ! PE_nq --> Potential enrgy: 1/2 Mp \omega_\nu^2 x_\nu^2
  ! KE_nq --> Kinetic enrgy: 1/2 Mp \omega_\nu^2 x_\nu^2
  ! p_q is the momentum on the nuclei \hbar\2\l_\bq\nu \SQRT(n_{q\nu, T}+ 1/2)
  !  
  ! ALLOCATE TABLES
  REAL(DP), ALLOCATABLE    :: equil_p(:, :, :), crystal_pos(:, :, :), qA(:, :), qB(:, :)
  REAL(DP), ALLOCATABLE    :: T_fact(:, :), DW_fact(:, :), DWp_fact(:, :), Tp_fact(:, :) 
  ! for displacements
  REAL(DP), ALLOCATABLE    :: Cx_matA(:, :), Cx_matB(:, :), Cx_matAB(:, :), Bx_vect(:)
  ! for momenta/velocities 
  REAL(DP), ALLOCATABLE    :: Cpx_matA(:, :), Cpx_matB(:, :), Cpx_matAB(:, :)
  ! matrices to account for the coupling terms between different phonon branches ! 
  REAL(DP), ALLOCATABLE    :: sum_er_D(:, :), sum_diag_D(:, :), sum_er_B(:) 
  REAL(DP), ALLOCATABLE    :: sum_diag_B(:), sum_er_B2(:), sum_diag_B2(:) 
  REAL(DP), ALLOCATABLE    :: D_tau(:, :, :), P_tau(:, :, :), ratio_zg(:)! displacements and velocities
  REAL(DP), ALLOCATABLE    :: R_mat(:, :), E_vect(:, :), D_vect(:, :), F_vect(:, :)
  ! D_tau  : atomic displacements
  ! z_nq_A : eigenvectors for q-points in set A 
  ! z_nq_B : eigenvectors for q-points in set B
  ! R_mat, E_vect, D_vect, F_vect : are used to compute the minimization of the 
  ! error coming from the off diagonal terms --> sum_er_D ! 
  ! sum_diag_D : the sum of diagonal terms contributing to the T-dependent properties  
  !
  COMPLEX(DP)              :: z_zg(3 * nat, 3 * nat, nq)
  COMPLEX(DP)              :: imagi 
  COMPLEX(DP), ALLOCATABLE :: z_nq_synch(:, :, :), z_nq_A(:, :, :), z_nq_B(:, :, :)
  ! singular value decomposition matrices U = R*conj(L)  
  !
  INTEGER                         :: INFO, N_dim, M_dim, K_dim, L_dim, LWORK
  REAL(DP),       ALLOCATABLE     :: RWORK(:), S_svd(:)
  COMPLEX(DP),    ALLOCATABLE     :: M_over(:, :, :), U_svd(:, :, :), U_svd_d(:, :), dotp_mat(:, :)
  COMPLEX(DP),    ALLOCATABLE     :: L_svd(:, :), R_svd(:, :), WORK(:), U_svd_d_new(:, :)
  COMPLEX*16 dum( 1 )    ! for the ZGEEV
  !
  !  
  ! constants to be used
  hbar    = H_PLANCK_SI / tpi ! reduce Plnack constant
  ang     = 1.0E-10            ! angstrom units
  imagi   = (0.0d0, 1.0d0) !imaginary unit
  ! Set intitial values
  nq_tot  = dim1 * dim2 * dim3 
  nat3    = 3 * nat
  pn      = 2**(nat3 - 1)
  ! pointless to allocate more signs for a large number of branches
  IF ( nat3 > 12) pn = 2**(12 - 1) 
  !
  !
  ! create equilibrium configuration
  ALLOCATE(equil_p(nq_tot, nat, 3), crystal_pos(nq_tot, nat, 3))
  !
  CALL create_supercell(at, tau, alat, dim1, dim2, dim3, nat, equil_p, crystal_pos)
  !
  CALL mp_bcast(equil_p, ionode_id, world_comm)
  !
  CALL mp_bcast(crystal_pos, ionode_id, world_comm)
  !
  IF ( TRIM(flscf) /= ' ' .AND. ionode) WRITE(*,*) "=============================================="
  IF ( TRIM(flscf) /= ' ') CALL read_input_file( 'PW', flscf )
  IF ( TRIM(flscf) /= ' ' .AND. ionode) WRITE(*,*) "=============================================="
  !
  ! 
  IF (ionode) THEN
  !
    IF (dim1 < 10) WRITE(pt_1,'(1i1)') dim1
    IF (dim2 < 10) WRITE(pt_2,'(1i1)') dim2
    IF (dim3 < 10) WRITE(pt_3,'(1i1)') dim3
    IF (dim1 >= 10 .AND. dim1 < 100 ) WRITE(pt_1,'(1i2)') dim1
    IF (dim2 >= 10 .AND. dim2 < 100 ) WRITE(pt_2,'(1i2)') dim2
    IF (dim3 >= 10 .AND. dim3 < 100 ) WRITE(pt_3,'(1i2)') dim3
    IF (dim1 >= 100) WRITE(pt_1,'(1i3)') dim1
    IF (dim2 >= 100) WRITE(pt_2,'(1i3)') dim2
    IF (dim3 >= 100) WRITE(pt_3,'(1i3)') dim3
    !
    IF ( TRIM(flscf) /= ' ') THEN
          OPEN (unit = 83, file = 'equil-scf_' //  TRIM( pt_1 ) // TRIM( pt_2 ) & 
                // TRIM( pt_3 ) // '.in', status = 'unknown', form = 'formatted')
          WRITE(83,*) "&control"
          WRITE(83,'(100A)') "  calculation = '", TRIM(calculation),"'"
          WRITE(83,'(100A)') "  restart_mode = '", TRIM(restart_mode),"'"
          WRITE(83,'(100A)') "  prefix = 'equil-", TRIM(prefix),"'"
          WRITE(83,'(100A)') "  pseudo_dir = '", TRIM(pseudo_dir),"'"
          WRITE(83,'(100A)') "  outdir = '", TRIM(outdir),"'"
          WRITE(83,'(100A)') "/"
          WRITE(83,'(100A)') "&system"
          WRITE(83,'(100A)') "  ibrav = 0"
          WRITE(83,'(A7,1i5)') "  nat =", nat * nq_tot
          WRITE(83,'(A8,1i5)') "  ntyp =", ntyp
          WRITE(83,'(A11,1F7.2)') "  ecutwfc =", ecutwfc
          WRITE(83,'(100A)') "  occupations = '", TRIM(occupations),"'"
          WRITE(83,'(100A)') "  smearing = '", TRIM(smearing),"'"
          WRITE(83,'(A11,1D10.1)') "  degauss =", degauss
          WRITE(83,'(100A)') "/"
          WRITE(83,'(100A)') "&electrons"
          WRITE(83,'(100A)') "  diagonalization = '",TRIM(diagonalization),"'"
          WRITE(83,'(100A)') "  mixing_mode= '",TRIM(mixing_mode),"'"
          WRITE(83,'(A16,1F4.2)') "  mixing_beta = ", mixing_beta
          WRITE(83,'(A12,1D10.1)') "  conv_thr = ", conv_thr
          WRITE(83,'(100A)') "/"
          !
          WRITE(83,'(100A)') "ATOMIC_SPECIES"
          DO k = 1, ntyp ! type of atom
            WRITE(83,'(A6, 1F8.3, A, A)') atm(k), amass(k), ' ', TRIM(atom_pfile(k))
          ENDDO
          !
          WRITE(83,'(100A)') "K_POINTS automatic"
          WRITE(83,'(6i4)') NINT(nk1/DBLE(dim1)), NINT(nk2/DBLE(dim2)), &
                            NINT(nk3/DBLE(dim3)), k1, k2, k3
          !
          WRITE(83,'(100A)') "CELL_PARAMETERS (angstrom)"
          WRITE(83,'(3F16.8)') DBLE(at(:, 1) * alat * dim1 * BOHR_RADIUS_ANGS )
          WRITE(83,'(3F16.8)') DBLE(at(:, 2) * alat * dim2 * BOHR_RADIUS_ANGS )
          WRITE(83,'(3F16.8)') DBLE(at(:, 3) * alat * dim3 * BOHR_RADIUS_ANGS )
          !
          WRITE(83,'(100A)') "ATOMIC_POSITIONS (angstrom)"
         !
    ENDIF ! flscf
    !
    filename = 'equil_pos.dat'
    OPEN (unit = 70, file = filename, status = 'unknown', form = 'formatted')
    WRITE(70,*) "Number of atoms", nat * dim1 * dim2 * dim3
    WRITE(70,*) 'equilibrium positions, (Ang):'
    DO k = 1, nat
      DO i = 1, nq_tot
        WRITE(70,'(A6, 3F13.8)') atm(ityp(k)), equil_p(i, k, :)
        IF ( TRIM(flscf) /= ' ') WRITE(83,'(A6, 3F13.8)') atm(ityp(k)), equil_p(i, k, :)
      ENDDO
    ENDDO
    CLOSE(70)
    IF ( TRIM(flscf) /= ' ') close(83)
  !
  ENDIF ! ionode
  !
  ! Inititialize eigenvectors matrix
  z_zg    = (0.d0, 0.d0)
  ! convert eigenvectors to mass-unscalled
  DO i = 1, nat3
    DO na = 1, nat
      nta = ityp(na)
      DO ipol = 1, 3
        DO qp = 1, nq
          z_zg((na - 1) * 3 + ipol, i, qp) = z_nq((na - 1) * 3 + ipol, i, qp) * SQRT(amu_ry*amass(nta))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
!
!
! Frequency check
  freq = 0.0d0
  dotp = 0.0d0
  IF (ionode) THEN
    WRITE(*,*) 
    DO qp = 1, nq
      DO i = 1, nat3 
     !   IF (w2(i, qp) .lt. 1.0E-8) THEN
        IF (w2(i, qp) .lt. 0.0d0) THEN
            dotp = -SQRT(ABS(w2(i, qp))) * ry_to_thz * tpi
            WRITE(*,'(A, F14.6, A, I2, A, I2)') "  WARNING: Negative freq. (Thz):", dotp, & 
                                                "  mode: ",i, "  q-point: ", qp 
            WRITE(*,*)
            freq(i, qp) = SQRT(ABS(w2(i, qp)))
        ELSE
            freq(i, qp) = SQRT(ABS(w2(i, qp)))
        ENDIF
      ENDDO
    ENDDO
    IF (dotp .lt. 0.0d0) WRITE(*,'(100A)') "  We freeze them, but dynamically stable phonons are recommended ..."
    WRITE(*,*)
  ENDIF
  CALL mp_bcast(freq, ionode_id, world_comm)
  !
  ph_w = freq * ry_to_thz * (1.0E12) * tpi ! correct frequency for phonons in SI
  !
  ! set amplitudes of displacements l_q = \sigma_\bq\nu and momenta 
  p_q = 0.0d0
  l_q = 0.0d0
  FE_q = 0.0d0
  DO qp = 1, nq
    DO i = 1, nat3 
       IF (w2(i, qp) .lt. 0.0d0 + eps2) THEN
        l_q(i, qp) = 0.d0
        p_q(i, qp) = 0.d0
      ELSE
      ! σν,T**2 = (2nν,T + 1) lν**2 ,
      ! V 0 + \sum \hbar ων /2 − kBT ln[1 + nB(ων ,T )]
        ! Free energy
        FE_q(i,qp) = hbar * ph_w(i, qp) / 2.0d0 - K_BOLTZMANN_SI * T * & 
                     LOG( 1.0d0 + 1.0d0 / (EXP(hbar * ph_w(i, qp) / (K_BOLTZMANN_SI * T)) - 1.0d0))
        !
        l_q(i, qp) = SQRT(hbar / ph_w(i, qp) / 2.0d0 / AMU_SI / ang**2.0d0) * & 
                     SQRT(DBLE(1.0d0 + 2.0d0 / (EXP(hbar * ph_w(i, qp) / (K_BOLTZMANN_SI * T)) - 1.0d0))) 
        p_q(i, qp) = hbar / SQRT(2.0d0) / (SQRT(hbar / ph_w(i, qp) / 2.0d0 /AMU_SI)) / ang * & !*1.0E-12& 
                     SQRT(DBLE(0.5d0 + 1.0d0 / (EXP(hbar * ph_w(i, qp) / (K_BOLTZMANN_SI * T)) - 1.0d0))) 
      ENDIF
                   ! we can multiply by 1.0E-12 to get 'picos'
    ENDDO
  ENDDO
  !     
  IF (q_external) THEN
      IF (q_in_cryst_coord .EQV. .FALSE.) THEN
      ! in both cases convert them to crystal 
          CALL cryst_to_cart(nq, q, at, -1)
      ELSE
          CALL cryst_to_cart(nq, q, at, -1)
      ENDIF
    ELSE 
    CALL cryst_to_cart(nq, q, at, -1)
  ENDIF
  !
  ! for accoustic modes put l_q\nu = 0 and p_q\nu = 0 so we freeze them
  !
  DO qp = 1, nq
    q_A = q(:, qp)  ! q_A to find IF q belongs in A
    IF (((ABS(q_A(1)) < eps)) .AND. & 
        ((ABS(q_A(2)) < eps)) .AND. &
        ((ABS(q_A(3)) < eps)))  THEN
        !
        l_q(1, qp) = 0.0d0
        l_q(2, qp) = 0.0d0
        l_q(3, qp) = 0.0d0
        p_q(1, qp) = 0.0d0
        p_q(2, qp) = 0.0d0
        p_q(3, qp) = 0.0d0
        !
    ENDIF
  ENDDO
  !
  ! To distinguish between different sets of qpoints, A, B, C
  ! to find how many points belong to set A and then allocate matrix accordingly
  ! NOTE that we want the qpoints always in crystal coordinates
  !
  ctrA  = 0
  ctrAB = 0
  PE_nq = 0.0d0
  KE_nq = 0.0d0
  FE_nq = 0.0d0
  !
  DO qp = 1, nq
    q_A = q(:, qp) + q(:, qp) ! q_A to find IF q belongs in A
    IF (((ABS(q_A(1)) < eps) .OR. (ABS(ABS(q_A(1)) - 1) < eps)) .AND. &
        ((ABS(q_A(2)) < eps) .OR. (ABS(ABS(q_A(2)) - 1) < eps)) .AND. &
        ((ABS(q_A(3)) < eps) .OR. (ABS(ABS(q_A(3)) - 1) < eps))) THEN
  !     WRITE(*,*) "set A", qp, q(:, qp)
       ctrA  = ctrA + 1
       ctrAB = ctrAB + 1
       DO i = 1, nat3 
         FE_nq = FE_nq + FE_q(i, qp) / RYDBERG_SI
         PE_nq = PE_nq + 0.5d0 * AMU_SI * ph_w(i, qp)**2.d0 * l_q(i, qp)**2.d0 * ang**2.0d0 / RYDBERG_SI
         KE_nq = KE_nq + 0.5d0 / AMU_SI * p_q(i, qp)**2.d0 * ang**2.0d0 / RYDBERG_SI
       ENDDO
    ELSE
       ctrAB = ctrAB + 1
       DO i = 1, nat3 
         FE_nq = FE_nq + 2.0d0 * FE_q(i, qp) / RYDBERG_SI
         PE_nq = PE_nq + 2.0d0 * 0.5d0 * AMU_SI * ph_w(i, qp)**2.d0 * l_q(i, qp)**2.d0 * ang**2.0d0 / RYDBERG_SI 
         KE_nq = KE_nq + 2.0d0 * 0.5d0 / AMU_SI * p_q(i, qp)**2.d0 * ang**2.0d0 / RYDBERG_SI
       ENDDO
    ENDIF
  ENDDO 
  !
  ctrB = ctrAB - ctrA
  IF (ionode) THEN
    WRITE(*,*) "=============================================="
    WRITE(*,*)
    WRITE(*,'(A30, 1F13.8, A4)') "Total vibrational energy: ", PE_nq + KE_nq, "Ry" 
    WRITE(*,*)
    WRITE(*,'(A30, 1F13.8, A4)') "Potential energy: ", PE_nq, "Ry" 
    WRITE(*,*)
    WRITE(*,'(A30, 1F13.8, A4)') "Kinetic energy: ", KE_nq, "Ry" 
    WRITE(*,*)
    WRITE(*,'(A30, 1F13.8, A4)') "Vibrational free energy: ", FE_nq, "Ry" 
    WRITE(*,*)
    WRITE(*,*) "Note that the total energy output from a DFT-ZG calculation" 
    WRITE(*,*) "calculation accounts for half the total vibrational energy"
    WRITE(*,*)
    WRITE(*,*) "=============================================="
    !
    WRITE(*,*) 
    WRITE(*,'(A26, 3I6)') "Points in sets AB, A, B :", ctrAB, ctrA, ctrB
    WRITE(*,*) 
    WRITE(*,*) "=============================================="
  ENDIF
  ! 
  ALLOCATE(qA(ctrA, 3), qB(ctrB, 3), z_nq_A(nat3, nat3, ctrA), z_nq_B(nat3, nat3, ctrB))  
  ALLOCATE(Cx_matAB(nat3, ctrAB), Cx_matA(nat3, ctrA), Cx_matB(nat3, ctrB))
  ALLOCATE(Cpx_matAB(nat3, ctrAB), Cpx_matA(nat3, ctrA), Cpx_matB(nat3, ctrB))
  ALLOCATE(D_tau(nq_tot, nat, 3), P_tau(nq_tot, nat, 3), Rlist(nq_tot, 3))
  ALLOCATE(T_fact(nat, 3), DW_fact(nat, 3), DWp_fact(nat, 3), Tp_fact(nat, 3))
  !
  Cx_matAB  = 0
  Cpx_matAB = 0
  !
  ! Generate lattice vectors in crystal coordinates   
  ctr2 = 1
  !
  DO i = 0, dim3 - 1
    DO j = 0, dim2 - 1
      DO  k = 0, dim1 - 1
        Rlist(ctr2, 1) = k
        Rlist(ctr2, 2) = j
        Rlist(ctr2, 3) = i !(/ k, j, i /)
     !WRITE(*,*) Rlist(ctr2,:)
        ctr2 = ctr2 + 1
      ENDDO
    ENDDO
  ENDDO
  !
  ! Now main synch proc as in paper. Do this procedure for every pn. 
  ! 
  IF (synch) THEN
    ! SVD parameters
    M_dim = nat3
    N_dim = nat3
    K_dim = MIN(M_dim, N_dim)
    L_dim = MAX(M_dim, N_dim) 
    !
    ALLOCATE(M_over(nat3, nat3, pn - 1), U_svd(nat3, nat3, pn - 1), z_nq_synch(nat3, nat3, ctrAB)) 
    ALLOCATE(U_svd_d(nat3, pn - 1), dotp_mat(nat3, nat3), U_svd_d_new(nat3, pn - 1))
    ALLOCATE(L_svd(M_dim, K_dim), R_svd(K_dim, N_dim),S_svd(K_dim))
    z_nq_synch = (0.0d0 , 0.0d0)
    ! query workspace
    !
    LWORK = 5 * nat3 !MAX(1, 2*K_dim+L_dim)
    !
    ALLOCATE( RWORK( LWORK ) )
    ALLOCATE( WORK( LWORK ) )
    LWORK = -1

    call ZGESVD('A','A', nat3, nat3, M_over(:, :, 1), nat3, S_svd, L_svd, &
                       nat3, R_svd, nat3, WORK, LWORK, RWORK, INFO)
     !
    LWORK = INT(WORK(1)) + 1
     !
    IF( LWORK > SIZE( WORK ) ) THEN
      DEALLOCATE( WORK )
      ALLOCATE( WORK( LWORK ) )
    ENDIF
     !
     !
    DO i = 0, ctrAB - pn, pn
      z_nq_synch(:, :, i + 1) = z_zg(:, :, i + 1)
      ! z_nq_synch(:, :, ctrAB-pn -i + 1) = z_zg(:, :, i + 1)
      DO ii = 1, pn - 1
        M_over = 0.d0
        ! Construct the overlap matrix M_{\nu,\nu'}
        S_svd = 0.d0
        DO p = 1, nat3
          DO j = 1, nat3
            DO k = 1, nat3 ! sum over \k,\a
               M_over(j, p, ii) = M_over(j, p, ii) + (z_zg(k, j, ii + i + 1) * CONJG(z_nq_synch(k, p, ii + i)))
            ENDDO ! k-loop
          ENDDO ! j-loop
        ENDDO ! p-loop 
        ! perform singular value decomposition
        call ZGESVD('A', 'A', nat3, nat3, M_over(:, :, ii), nat3, S_svd, L_svd, &
                    nat3, R_svd, nat3, WORK, LWORK, RWORK, INFO)
        U_svd(:, :, ii) = MATMUL(TRANSPOSE(CONJG(R_svd)),TRANSPOSE(CONJG(L_svd)))
        call ZGEEV('N', 'N', nat3, U_svd(:, :, ii), nat3, U_svd_d(:, ii), dum, 1, dum, 1, &
                    WORK, LWORK, RWORK, INFO)
        !
        M_over = 0.0d0
        DO p = 1, nat3
          DO j = 1, nat3
            DO k = 1, nat3 ! sum over \k,\a
               M_over(j, p, ii) = M_over(j, p, ii) + (z_zg(k, j, ii + i + 1) * CONJG(z_nq_synch(k, p, ii + i)))
            ENDDO ! k-loop
          ENDDO ! j-loop
        ENDDO ! p-loop 
        DO qp = 1, nat3
         DO k = 1, nat3
            dotp_mat(qp, k) =  CONJG(M_over(qp, qp, ii)) * CONJG(U_svd_d(k, ii))
         ENDDO
        ENDDO
        dotp_mat = ABS(REAL(dotp_mat))
        DO qp = 1, nat3
           p = MAXLOC(REAL(dotp_mat(qp,:)), 1)
           U_svd_d_new(qp, ii) = U_svd_d(p, ii)
        ENDDO   
        DO qp = 1, nat3
          DO k = 1, nat3
               z_nq_synch(k, qp, ii + i + 1) = U_svd_d_new(qp, ii) * z_zg(k, qp, ii + i + 1)
           ENDDO
        ENDDO
        z_zg(:, :, ii + i + 1) = z_nq_synch(:, :, ii + i + 1)
      ENDDO ! ii-loop
    ENDDO ! i-loop
     ! Here we synchronize the remaining eigenvectors IF ctrAB is not divided by pn
    IF (mod(ctrAB, pn) > 0) THEN
      ctr = ctrAB - mod(ctrAB, pn)
      z_nq_synch(:, :, ctr + 1) = z_zg(:, :, ctr + 1)
      DO ii = 1, mod(ctrAB, pn) - 1
      M_over = 0.0d0
      ! Construct the overlap matrix M_{\nu,\nu'}
      S_svd = 0.0d0
      DO p = 1, nat3
        DO j = 1, nat3
          DO k = 1, nat3
             M_over(j, p, ii) = M_over(j, p, ii) + (z_zg(k, j, ii + ctr + 1) * CONJG(z_nq_synch(k, p, ii + ctr)))
          ENDDO ! k-loop
        ENDDO ! j-loop
      ENDDO ! p-loop 
      ! perform singular value decomposition
      call ZGESVD('S','S', nat3, nat3, M_over(:, :, ii), nat3, S_svd, L_svd, &
             nat3, R_svd, nat3, WORK, LWORK, RWORK,INFO)
      ! ZGESVD returns R_svd**H (the hermitian TRANSPOSE of R_svd)
      U_svd(:, :, ii) = MATMUL(TRANSPOSE(CONJG(R_svd)),TRANSPOSE(CONJG(L_svd)))
      call ZGEEV( 'N', 'N', nat3, U_svd(:, :, ii), nat3, U_svd_d(:, ii), dum, 1, dum, 1, &
                     WORK, LWORK, RWORK, INFO )
      M_over = 0.0d0
      DO p = 1, nat3
        DO j = 1, nat3
          DO k = 1, nat3 ! sum over \k,\a
             M_over(j, p, ii) = M_over(j, p, ii) + (z_zg(k, j, ii + ctr + 1) * CONJG(z_nq_synch(k, p, ii + ctr)))
          ENDDO ! k-loop
        ENDDO ! j-loop
      ENDDO ! p-loop 
      DO qp = 1, nat3
        DO k = 1, nat3
          dotp_mat(qp, k) =  CONJG(M_over(qp, qp, ii)) * CONJG(U_svd_d(k, ii))
        ENDDO
      ENDDO
      dotp_mat = ABS(DBLE(dotp_mat))
      DO qp = 1, nat3
        p = MAXLOC(DBLE(dotp_mat(qp, :)), 1)
        U_svd_d_new(qp, ii) = U_svd_d(p, ii)
      ENDDO   
      !
      DO qp = 1, nat3
        DO k = 1, nat3
          z_nq_synch(k, qp, ii + ctr + 1) =  U_svd_d_new(qp, ii) * z_zg(k, qp, ii + ctr + 1)
        ENDDO
      ENDDO
      ! overwrite z_zg
      z_zg(:, :, ii + ctr + 1) = z_nq_synch(:, :, ii + ctr + 1)
      ENDDO ! ii-loop
    ENDIF ! mod(ctrAB, pn)
    DEALLOCATE(WORK, RWORK)
  !
  ! DEALLOCATE matrices for synch proc.
    DEALLOCATE(M_over, U_svd, z_nq_synch, U_svd_d, U_svd_d_new, dotp_mat)
    DEALLOCATE(L_svd, R_svd, S_svd)
  !
  ENDIF ! end IF synch is TRUE
  !
  !
  ! Initialize sign matrices: 
  ! sign matrices Mx and My . Total entries of Mx are 2^nmodes/2. 
  ! Divided by two so we get only independent entries:
  ! i.e. [1 1 1 1 1 1] gives same result to [- 1 -1 -1 -1 -1 -1]. 
  ! 
  ALLOCATE(M_mat(2 * pn, nat3), Mx_mat(pn, nat3), Mx_mat_or(pn, nat3), V_mat(2))
  M_mat = 1 ! initialize M_mat
  !
  ! Initialize sign matrix
  IF (nat3 .LE. 12) THEN
    M_mat = 1
    ! ELSE Initialize with random entries of +1 or -1 that are only used 
    ! for systems with large unit-cells.
  ELSE
    DO i = 1, nat3
      DO p = 1, 2*pn
         call random_number(u_rand)
         IF (u_rand .LE. 0.5) M_mat(p, i) = -1
         IF (u_rand .GT. 0.5) M_mat(p, i) =  1
      ENDDO
    ENDDO
  ENDIF
  !
  ! initialize V_mat whose entries will generate the deterministic sign matrices
  V_mat = (/ 1, -1/) 
  DO i = 1, nat3
    ctr = 1
    DO p = 1, 2**(i - 1)
      DO qp = 1, 2
        DO k = 1, 2**(nat3 - i)
          IF (ctr > 2 * pn) EXIT ! in case there many branches in the system and 
                                 ! we do not need to allocate more signs              
          M_mat(ctr, i) = V_mat(qp)
          ctr = ctr + 1
          IF (ctr > 2 * pn) EXIT              
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  ! NOTE: In M_mat the first half entries are the independent set of signs (2** (nat3- 1)) !
  !       The rest entries are the antithetics !!  
  ! checks
  !IF (ionode) WRITE(*,*) "Sign matrix"
  !IF (ionode) WRITE(*,*) "-----------"
  DO j = 1, 2 * pn !2** (nat3)
    IF (MOD(j, 2) == 0) M_mat(j, :) = -1 * M_mat(j, :)
  !  IF (ionode) WRITE(*,'(100i3)') M_mat(j, :)
  ENDDO
  ! checks_done     
  !
  combs = 0! how many unique pairs x1x2, x1x3,... we have. 
  DO i = 1, nat3 - 1
    combs = combs + i
  ENDDO
  combs_all = 2 * combs + nat3 !; % with x1^2, x2^2 ...
  !
  ! combs_all refere also to  all possible pais ({\k,\a}, {\k' \a'})
  ! 
  ALLOCATE(ratio_zg(combs_all))
  ratio_zg = 0.d0
  IF (compute_error) THEN
    ALLOCATE(sum_er_D(combs_all, INT(ctrAB / pn) + 1), sum_er_B(combs_all)) 
    ALLOCATE(sum_er_B2(combs_all * NINT(nq_tot / 2.0d0)), sum_diag_B2(combs_all * NINT(nq_tot / 2.0d0)))
    ! I add one because of the reminder when ctrAB is not divided by pn
    ALLOCATE(sum_diag_D(combs_all, INT(ctrAB / pn) + 1), sum_diag_B(combs_all))
    ALLOCATE(R_mat(combs_all, combs_all), D_vect(combs_all, ctrAB), F_vect(combs_all, ctrAB))
    ALLOCATE(Bx_vect(combs_all), E_vect(combs_all, ctrAB))
  ENDIF
  !
  
  !
  !
  IF (ionode) THEN
    IF (T < 10000.d0 - eps .AND. T > 1000.d0 - eps) WRITE(pt_T,'(1F7.2)') T 
    IF (T < 1000.d0 - eps .AND. T > 100.d0 - eps) WRITE(pt_T,'(1F6.2)') T 
    IF (T < 100.d0 - eps .AND. T > 10.d0 - eps) WRITE(pt_T,'(1F5.2)') T 
    IF (T < 10.d0 - eps .AND. T > 0.d0 - eps) WRITE(pt_T,'(1F4.2)') T 
    !
    !
    filename = 'ZG-configuration_' // TRIM( pt_T ) // 'K.dat' 
    OPEN (unit = 80, file = filename, status = 'unknown', form = 'formatted')
    filename = 'ZG-velocities_' // TRIM( pt_T ) // 'K.dat'
    OPEN (unit = 81, file = filename, status = 'unknown', form = 'formatted')
    !
    IF ( TRIM(flscf) /= ' ') THEN
        OPEN (unit = 82, file = 'ZG-scf_' //  TRIM( pt_1 ) // TRIM( pt_2 ) // TRIM( pt_3 ) // & 
                                '_' // TRIM( pt_T ) // 'K.in', status = 'unknown', form = 'formatted')
        WRITE(82,*) "&control"
        WRITE(82,'(100A)') "  calculation = '", TRIM(calculation),"'"
        WRITE(82,'(100A)') "  restart_mode = '", TRIM(restart_mode),"'"
        WRITE(82,'(100A)') "  prefix = 'ZG-", TRIM(prefix),"'"  
        WRITE(82,'(100A)') "  pseudo_dir = '", TRIM(pseudo_dir),"'"  
        WRITE(82,'(100A)') "  outdir = '", TRIM(outdir),"'"  
        WRITE(82,'(100A)') "/"
        WRITE(82,'(100A)') "&system"
        WRITE(82,'(100A)') "  ibrav = 0"  
        WRITE(82,'(A7,1i5)') "  nat =", nat * nq_tot
        WRITE(82,'(A8,1i5)') "  ntyp =", ntyp
        WRITE(82,'(A11,1F7.2)') "  ecutwfc =", ecutwfc 
        WRITE(82,'(100A)') "  occupations = '", TRIM(occupations),"'"
        WRITE(82,'(100A)') "  smearing = '", TRIM(smearing),"'"
        WRITE(82,'(A11,1D10.1)') "  degauss =", degauss
        WRITE(82,'(100A)') "/"
        WRITE(82,'(100A)') "&electrons"
        WRITE(82,'(100A)') "  diagonalization = '",TRIM(diagonalization),"'"
        WRITE(82,'(100A)') "  mixing_mode= '",TRIM(mixing_mode),"'"
        WRITE(82,'(A16,1F4.2)') "  mixing_beta = ", mixing_beta
        WRITE(82,'(A12,1D10.1)') "  conv_thr = ", conv_thr
        WRITE(82,'(100A)') "/"
        !
        WRITE(82,'(100A)') "ATOMIC_SPECIES"
        DO k = 1, ntyp ! type of atom
          WRITE(82,'(A6, 1F8.3, A, A)') atm(k), amass(k), ' ', TRIM(atom_pfile(k))
        ENDDO
        !
        WRITE(82,'(100A)') "K_POINTS automatic"
        WRITE(82,'(6i4)') NINT(nk1/DBLE(dim1)), NINT(nk2/DBLE(dim2)), & 
                          NINT(nk3/DBLE(dim3)), k1, k2, k3
        !
        WRITE(82,'(100A)') "CELL_PARAMETERS (angstrom)"
        WRITE(82,'(3F16.8)') DBLE(at(:, 1) * alat * dim1 * BOHR_RADIUS_ANGS )
        WRITE(82,'(3F16.8)') DBLE(at(:, 2) * alat * dim2 * BOHR_RADIUS_ANGS )
        WRITE(82,'(3F16.8)') DBLE(at(:, 3) * alat * dim3 * BOHR_RADIUS_ANGS )
        !
        WRITE(82,'(100A)') "ATOMIC_POSITIONS (angstrom)"
       !
    ENDIF ! flscf
  ENDIF ! ionode
  !
  IF (single_ph_displ) THEN
    pn = 2
    ! when single_ph_displ = .true. we set pn = 2
    ! pointless to allocate signs for signle_phonon_displacements
    WRITE(*, *) " WARNING: 'single_ph_displ' flag is on, so error is not minimized" 
  ENDIF
  !
  ! Instead of taking all possible permutations which are pn = 2** (nmodes- 1)! 
  ! we just select possible permutations until the error is lower than a
  ! threshold. The lower the threshold the longer the algorithm can take.
  DO kk = 1, niters
  ! Allocate original matrices ! half the entries of M_mat
  ! We also make the inherent choice that each column of Mx_mat_or 
  ! has the same number of positive and negative signs 
    Mx_mat_or = 1
    DO i = 1, 2 * pn / 4, 2
      Mx_mat_or(i, :) = M_mat(i, :)
    ENDDO
    ctr = 1
    DO i = 2, 2 * pn / 4, 2
      Mx_mat_or(i, :) = M_mat(2 * pn + 1 - i, :)
      ctr = ctr + 1
    ENDDO
    DO i = 2 * pn / 4 + 1, 2 * pn / 2, 2
      Mx_mat_or(i, :) = M_mat(i + 1, :)
      ctr = ctr + 1
    ENDDO
    DO i = 2 * pn / 4 + 2, 2 * pn / 2, 2
      Mx_mat_or(i, :) = M_mat(2 * pn + 2 - i, :)
      ctr = ctr + 1
    ENDDO
    !
    !
    !
    DO i = pn, 1, -1
    ! To generate integer numbers from 1 to pn
    ! and ALLOCATE the M_x,y matrices ! so we DO not have specific order !
      call random_number(u_rand)
      ii = 1 + FLOOR(i * u_rand) 
      Mx_mat(i, :) = Mx_mat_or(ii, :)
      Mx_mat_or(ii, :) = Mx_mat_or(i, :) ! so I DO not repeat this entry
    ENDDO  
    ! DO q-points in sets A n B
    ! based on the sets of signs we genereated from the above loop
    DO ii = 1, ctrAB ! loop over qpoints
    ! change the signs of Mx in every 2^(nmodes-2) entries 
    ! I use Mx_mat_or to apply concatenate matrices and obtain set E
      IF (mod(ii, pn) .EQ. 1) THEN
        Mx_mat_or(1 : pn - 1, :) = Mx_mat(2 : pn, :) ! second element goes to the top
        Mx_mat_or(pn, :) = Mx_mat(1, :) ! first element goes to the bottom
      ENDIF
      Mx_mat = Mx_mat_or
      ! to take antithetics every pn
      ! Remember always the error goes to very small as long we have equal number of + and - signs
      ! and displacements remain around equilibrium ! 
      IF (mod(ii, pn).EQ.1) THEN
          Mx_mat = -Mx_mat
      ENDIF
      !
      ! 
      ! Cx_matAB contains all the sigmas with the appropriates signs 
      IF (MOD(ii, pn) > 0) THEN
        DO k = 1, nat3
          Cx_matAB(k, ii)  = l_q(k, ii) * Mx_mat(mod(ii, pn), k) ! mod so values of every pn q-points are repeated
          Cpx_matAB(k, ii) = p_q(k, ii) * Mx_mat(mod(ii, pn), k) ! mod so values of every pn q-points are repeated
        ENDDO
      ELSE
        DO k = 1, nat3
          Cx_matAB(k, ii)  = l_q(k, ii) * Mx_mat(pn, k)
          Cpx_matAB(k, ii) = p_q(k, ii) * Mx_mat(pn, k)
        ENDDO
      ENDIF
      ! create R matrix that contains Re[e_ka^nu(q) * e_k'a'^nu'* (q)]
      IF (compute_error) THEN
        R_mat = 0.0d0
        D_vect = 0.0d0
        ctr = 1
        DO p = 1, nat3
          DO qp = 1, nat3 !those are for \k \a, \k' \a'
            ctr2 = 1
          !------------------------------
            DO i = 1, nat3
              DO j = 1, nat3
                D_vect(ctr, ii) = D_vect(ctr, ii) + DBLE(z_zg(p, i, ii) * CONJG(z_zg(qp, j, ii))) * & 
                                                         Cx_matAB(i, ii) * Cx_matAB(j, ii)
                R_mat(ctr, ctr2) = DBLE(z_zg(p, i, ii) * CONJG(z_zg(qp, j, ii)))
                ctr2 = ctr2 + 1
              ENDDO
            ENDDO
          !----------------------------
            ctr = ctr + 1
          ENDDO
        ENDDO ! end p loop ! R_mat is filled 
        !
        Bx_vect = 0.d0
        ! Bx_vect will contain all the cross terms v .neq. v' and diagonal for each q
        ctr = 1
        DO i = 1, nat3
          DO j = 1, nat3
            Bx_vect(ctr) = Cx_matAB(i, ii) * Cx_matAB(j, ii)
            ctr = ctr + 1
          ENDDO
        ENDDO
    !
        E_vect(:, ii) = 0.d0
    !   E_vect contains only the diagonal terms (i.e. v = v')
        ctr = 1
        DO p = 1, nat3 ! those are for \k \a, \k' \a'
          DO qp = 1, nat3
            DO i = 1, nat3
              DO j = 1, nat3
                IF (j == i) THEN
                  E_vect(ctr, ii) = E_vect(ctr, ii) + DBLE(z_zg(p, i, ii) * CONJG(z_zg(qp, j, ii))) * Cx_matAB(i, ii)**2
                ENDIF
              ENDDO
            ENDDO
            ctr = ctr + 1
          ENDDO
        ENDDO ! p loop
    !   D_vect(:, ii) = MATMUL(R_mat, Bx_vect)
    !   D_vect contains the diagonal and the non-diagonal terms (i.e. v = v' and v .neq. v')
    !   E_vect contains only the diagonal terms (i.e. v = v')
    !   F_vect contains the error (i.e. each entry is the contribution from v \neq v') at each q point to minimize
    !   checks
        F_vect(:, ii) = D_vect(:, ii) - E_vect(:, ii)
    !
    ENDIF ! compute_error
    ENDDO ! ii loop over qpoints
    !
    !
    !Compute error 
    !
    !
    IF (compute_error) THEN
      sum_er_D = 0.0d0
      sum_er_B = 0.0d0
      ! sum_er_D : contains the error from \nu and \nu' every pn
      !!!!!!!
      IF (ionode) THEN
        WRITE(*,*) 
        WRITE(*,*) "Searching for optimum configuration..."
        WRITE(*,'(A11, i8)') "  Attempt #", kk 
      ENDIF
      !
      DO p = 1, combs_all
        ctr = 1
        DO i = 0, INT(ctrAB / pn) - 1 
          sum_er_D(p, ctr) =SUM(F_vect(p, pn * i + 1 : pn * (i + 1))) ! pn) is the length of Mx
          ctr = ctr + 1
        ENDDO 
        ! Here we add the reminder IF ctrAB is not divided exactly by pn
        IF (mod(ctrAB, pn) > 0) THEN
          sum_er_D(p, ctr) = SUM(F_vect(p, ctrAB - mod(ctrAB, pn) + 1 : ctrAB)) ! add the remaining terms
        ENDIF
        ! evaluate also error from all q-points in B
        DO i = 1, ctrAB
          sum_er_B(p) = sum_er_B(p) + F_vect(p, i) ! 
        ENDDO 
        !
      ENDDO ! end p-loop
      ! 
      sum_er_B2 = 0.0d0
      ctr2 = 1
      DO j = 1, NINT(nq_tot / 2.0d0)
        DO p = 1, combs_all ! for ever \k,\a,\k',\a'
          DO i = 1, ctrAB ! over qpoints
            dotp = 0.0d0
            DO ii = 1, 3
              !dotp = dotp + q(i, ii) * Rlist(j, ii)!
              ! q(3, nq)
              dotp = dotp + q(ii, i) * Rlist(j, ii)
            ENDDO ! ii
            sum_er_B2(ctr2) = sum_er_B2(ctr2) + cos(tpi * dotp) * F_vect(p, i)
            !
          ENDDO ! i
          ctr2 = ctr2 + 1
        ENDDO ! p     
      ENDDO ! j
      !
      sum_diag_D = 0.0d0
      sum_diag_B = 0.0d0
      ! sum_diag_D : contains the diagonal terms 
      DO p = 1, combs_all
      ctr = 1
        DO i = 0, INT(ctrAB / pn) - 1 
          sum_diag_D(p, ctr) =SUM(E_vect(p, pn * i + 1 : pn * (i + 1))) ! pn) is the length of Mx
          ctr = ctr + 1
        ENDDO
        ! Here we add the reminder IF ctrAB is not divided exactly by pn
        IF (mod(ctrAB, pn) > 0) THEN
          sum_diag_D(p, ctr) = SUM(E_vect(p, ctrAB - mod(ctrAB, pn) + 1 : ctrAB)) ! add the remaining terms
        ENDIF
        DO i = 1, ctrAB
          sum_diag_B(p) = sum_diag_B(p) + E_vect(p, i) ! pn) is the length of Mx
        ENDDO
      ENDDO ! end p-loop
      sum_diag_B2 = 0.0d0
      ctr = 1
      sum_zg = 0
      ratio_zg = 0.0
      !
      DO p = 1, nat3 ! those are for \k \a, \k' \a'
        DO qp = 1, nat3
          IF (p == qp) THEN
            ratio_zg(ctr) = sum_er_B(ctr) / sum_diag_B(ctr)
            !!! WRITE(*,*) "Error from each branch", sum_diag_B(ctr), sum_er_B(ctr), p , qp, ratio_zg(ctr)
            IF (ABS(ratio_zg(ctr)) < error_thresh) THEN
              sum_zg = sum_zg + 1
            ENDIF   
          ENDIF
          ctr = ctr + 1
        ENDDO
      ENDDO
      !
      IF (ionode) WRITE(*,'(A, 1F12.6)') "      Total error:", SUM(ABS(ratio_zg)) / nat3
      IF (ionode) WRITE(*,*) 
      IF (ionode .AND. SUM(ABS(ratio_zg)) / nat3 < error_thresh ) & 
                  WRITE(*,*) "Optimum configuration found !"
    ENDIF ! compute_error
    !
    !IF (sum_zg == nat3) THEN
    IF (SUM(ABS(ratio_zg)) / nat3 < error_thresh ) THEN
      ctrA = 0
      ctrB = 0
      !
      DO qp = 1, ctrAB
        q_A = q(:, qp) +  q(:, qp) ! q_A to find IF q belongs in A
        IF (((ABS(q_A(1)) < eps) .OR. (ABS(ABS(q_A(1)) - 1) < eps)) .AND. &
            ((ABS(q_A(2)) < eps) .OR. (ABS(ABS(q_A(2)) - 1) < eps)) .AND. &
            ((ABS(q_A(3)) < eps) .OR. (ABS(ABS(q_A(3)) - 1) < eps))) THEN
              ctrA = ctrA + 1
              Cx_matA(:, ctrA) = Cx_matAB(:, qp)
              Cpx_matA(:, ctrA) = Cpx_matAB(:, qp)
              z_nq_A(:, :, ctrA) =  z_zg(:, :, qp)
              qA(ctrA, :) =  q(:, qp)
              IF (ABS(qA(ctrA, 1)) < eps) qA(ctrA, 1) = 0.0
              IF (ABS(qA(ctrA, 2)) < eps) qA(ctrA, 2) = 0.0
              IF (ABS(qA(ctrA, 3)) < eps) qA(ctrA, 3) = 0.0
        ELSE
              ctrB = ctrB + 1
              Cx_matB(:, ctrB) = Cx_matAB(:, qp)
              Cpx_matB(:, ctrB) = Cpx_matAB(:, qp)
              z_nq_B(:, :, ctrB) =  z_zg(:, :, qp)
              qB(ctrB,:) = q(:, qp)
              IF (ABS(qB(ctrB, 1)) < eps) qB(ctrB, 1) = 0.0
              IF (ABS(qB(ctrB, 2)) < eps) qB(ctrB, 2) = 0.0
              IF (ABS(qB(ctrB, 3)) < eps) qB(ctrB, 3) = 0.0
            !
        ENDIF
      ENDDO
      !
      IF (ionode) THEN
        IF (single_ph_displ) THEN
            WRITE(*,*) "Print single phonon displacements" 
            CALL single_phonon(nq_tot, nat, ctrB, ctrA, nat3, ityp, ntyp, & 
                               ntypx, qA, qB, amass, atm, equil_p, & 
                               Rlist, z_nq_B, z_nq_A, Cx_matB, & 
                               Cx_matA, Cpx_matB, Cpx_matA)
        ENDIF
        !     
        IF (compute_error) THEN
          WRITE(*,*)
          WRITE(*,'(A, 1F12.6)') " Sum of diagonal terms per q-point:", DBLE(SUM(sum_diag_B) / ctrAB)
          WRITE(*,'(A, 1F12.6,i8)') " Error and niter index:", SUM(ABS(ratio_zg)) / nat3, kk !
          WRITE(*,*)
        ENDIF
        WRITE(*,*) "=============================================="
        WRITE(*,*) 
        WRITE(*,*) "Print ZG configuration"
        WRITE(*,*) 
        !WRITE(80,*) "Sum of error per q-point and loop index:", SUM(sum_er_B)/ctrAB, kk !
        WRITE(80,'(A20, 1F8.2,A3)') 'Temperature is: ' , T ,' K'
        WRITE(80,*) "Atomic positions", nat * nq_tot
        WRITE(81,*) "ZG-Velocities (Ang/ps)"
      ENDIF
      ! Generate displacements and velocities.
      ! Remember nq_tot is also equal to the number of cells
      ! Here the displacements are generated according to 
      ! Np^(- 1/2)(Mo/Mk)^(1/2)[\sum_{q \in B} e^{1qR_p}e^v_{ka}(q)(x_{qv}+y_{q\nu})
      ! z_zg(nat3, nat3, nq)) 
      !
      D_tau = 0.0d0
      P_tau = 0.0d0
      ! 
      ! Main loop to construct ZG configuration
      !
      CALL fkbounds( nq_tot, lower_bnd, upper_bnd )
      !
      DO p = lower_bnd, upper_bnd !1, nq_tot
        ctr = 1
        DO k = 1, nat ! k represents the atom
          nta = ityp(k)
          DO i = 1, 3  ! i is for cart directions
          !
            DO qp = 1, ctrB
              dotp = 0.0d0
              DO ii = 1, 3
                dotp = dotp + qB(qp, ii) * Rlist(p, ii) * tpi ! dot product between q and R 
              ENDDO
              DO j = 1, nat3        
                D_tau(p, k, i) = D_tau(p, k, i) + SQRT(2.0d0 / nq_tot / amass(nta)) * DBLE(CMPLX(COS(dotp), SIN(dotp), kind=dp) &
                               * z_nq_B(ctr, j, qp) * (1.d0 + imagi) * Cx_matB(j, qp)) !) 
                P_tau(p, k, i) = P_tau(p, k, i) + SQRT(2.0d0 / nq_tot * amass(nta)) * DBLE(CMPLX(COS(dotp), SIN(dotp), kind=dp) &
                               * z_nq_B(ctr, j, qp) * (1.d0 + imagi) * Cpx_matB(j, qp)) / (amass(nta) * AMU_SI)
                ! Here we calculate the momenta of the nuclei and finally 
                !we divide by (amass(nta) *AMU_SI) to get the velocities.
              ENDDO
            ENDDO ! qp loop
            !
            IF (incl_qA) THEN ! If we want to include modes in set A. Those modes are known
                            ! to break the degeneracy for finite size systems
              DO qp = 1, ctrA
                dotp = 0.0d0
                DO ii = 1, 3
                  dotp = dotp + qA(qp, ii) * Rlist(p, ii)! 
                ENDDO
                DO j = 1, nat3
                  D_tau(p, k, i) = D_tau(p, k, i) + SQRT(1.0d0 / nq_tot / amass(nta)) * cos(tpi * dotp) &
                                           * DBLE(z_nq_A(ctr, j, qp) * Cx_matA(j, qp)) ! 
                  P_tau(p, k, i) = P_tau(p, k, i) + SQRT(1.0d0 / nq_tot * amass(nta)) * cos(tpi * dotp) &
                                           * DBLE(z_nq_A(ctr, j, qp) * Cpx_matA(j, qp)) / (amass(nta) * AMU_SI) !
                ENDDO
              ENDDO
            ENDIF ! IF incl_qA
            !
            ctr = ctr + 1 ! for k and i
            !IF (ABS(D_tau(p, k, i)) .GT. 5) CALL errore('ZG', 'Displacement very large', D_tau(p, k, i) )
            IF (ABS(D_tau(p, k, i)) .GT. 5) CALL errore('ZG', 'Displacement very large', 1)
            D_tau(p, k, i) = equil_p(p, k, i) + D_tau(p, k, i) ! add equil structure
          ENDDO ! end i for cart directions
        ENDDO ! end k loop over nat
      ENDDO ! end p loop over unit cells
      ! print displacements
      !
      CALL mp_sum(D_tau, inter_pool_comm)
      CALL mp_sum(P_tau, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
      IF (ionode) THEN
        DO k = 1, nat ! k represents the atom
          DO p = 1, nq_tot 
            WRITE(80,'(A6, 3F16.8)') atm(ityp(k)), D_tau(p, k, :) 
            IF (flscf /= ' ') WRITE(82,'(A6, 3F16.8)') atm(ityp(k)), D_tau(p, k, :) 
            WRITE(*,'(A6, 3F16.8)') atm(ityp(k)), D_tau(p, k, :) 
            WRITE(81,'(A6, 3F15.8)') atm(ityp(k)), P_tau(p, k, :) * 1.0E-12 ! multiply to obtain picoseconds 
          ENDDO
        ENDDO 
        !WRITE(80,*) "sign matrices" 
        !DO i = 1, pn
        !    WRITE(80,*) Mx_mat(i,:)
        !ENDDO
        WRITE(*,*) 
        WRITE(*,*) "=============================================="
        WRITE(*,*) 
        WRITE(*,'(100A)') ' Anisotropic mean-squared displacement tensor vs exact values (Ang^2):'
        WRITE(81,'(100A)') ' ZG-velocities vs exact velocities from momentum operator in second quantization:'
      ENDIF ! ionode
      !
      ! Exact anisotropic displacement parameter
      DW_fact  = 0.0d0
      DWp_fact = 0.0d0
      ctr = 1
      DO k = 1, nat ! k represents the atom
        nta = ityp(k)
        DO i = 1, 3  ! i is for cart directions
        !
          DO qp = 1, ctrB
            DO j = 1, nat3
              DW_fact(k, i) = DW_fact(k, i) + 2.0d0 / DBLE(nq_tot * amass(nta)) * z_nq_B(ctr, j, qp) & 
                                         * CONJG(z_nq_B(ctr, j, qp)) * Cx_matB(j, qp)**2
              DWp_fact(k, i) = DWp_fact(k, i) + 2.0d0 * amass(nta) / DBLE(nq_tot) * z_nq_B(ctr, j, qp) & 
                                        * CONJG(z_nq_B(ctr, j, qp)) * Cpx_matB(j, qp)**2 & 
                                        / ((amass(nta) * AMU_SI)**2) * 1.0E-24
            ENDDO
          ENDDO
          !
          IF (incl_qA) THEN ! 
            DO qp = 1, ctrA
              DO j = 1, nat3
                DW_fact(k, i) = DW_fact(k, i) + 1.0d0 / DBLE(nq_tot * amass(nta)) * z_nq_A(ctr, j, qp) & 
                                * CONJG(z_nq_A(ctr, j, qp)) * Cx_matA(j, qp)**2
                DWp_fact(k, i) = DWp_fact(k, i) + amass(nta) / DBLE(nq_tot) * z_nq_A(ctr, j, qp) * CONJG(z_nq_A(ctr, j, qp)) & 
                                * Cpx_matA(j, qp)**2 / ((amass(nta) * AMU_SI)**2) * 1.0E-24
              ENDDO
            ENDDO
          ENDIF
         !
          ctr = ctr + 1 ! for k and i
        ENDDO ! end i for cart directions
      ENDDO ! end k loop over nat
      !


      T_fact(:,:) = 0.d0
      Tp_fact(:,:) = 0.d0
      DO k = 1, nat
        DO p = 1, nq_tot
          T_fact(k, 1) =  T_fact(k, 1) + (D_tau(p, k, 1) - equil_p(p, k, 1))**2 / nq_tot
          T_fact(k, 2) =  T_fact(k, 2) + (D_tau(p, k, 2) - equil_p(p, k, 2))**2 / nq_tot
          T_fact(k, 3) =  T_fact(k, 3) + (D_tau(p, k, 3) - equil_p(p, k, 3))**2 / nq_tot
          Tp_fact(k, 1) =  Tp_fact(k, 1) + (P_tau(p, k, 1))**2 / nq_tot * 1.0E-24
          Tp_fact(k, 2) =  Tp_fact(k, 2) + (P_tau(p, k, 2))**2 / nq_tot * 1.0E-24
          Tp_fact(k, 3) =  Tp_fact(k, 3) + (P_tau(p, k, 3))**2 / nq_tot * 1.0E-24
        ENDDO
      ENDDO
      IF (ionode) THEN 
        WRITE(*,*)
        DO k = 1, nat
          nta = ityp(k)
          WRITE(*,'(A6, 2i2)') " Atom: ", k 
          WRITE(*,'(A8)') " --------"
          WRITE(*,'(A10, A6, 3F11.6)') " ZG_conf:", atm(nta), T_fact(k, 1), T_fact(k, 2), T_fact(k, 3)
      !    WRITE(*,'(A20, 3F11.6)') "  Exact values" 
          WRITE(*,'(A10, A6, 3F11.6)') " Exact:", atm(nta), DW_fact(k, 1), DW_fact(k, 2), DW_fact(k, 3)
          WRITE(*,*)
          !Here we print the mean velocities
          WRITE(81,'(A6, 3F12.8)') atm(nta), SQRT(Tp_fact(k, 1)), SQRT(Tp_fact(k, 2)), SQRT(Tp_fact(k, 3))
          WRITE(81,'(A6, 3F12.8)') atm(nta), SQRT(DWp_fact(k, 1)), SQRT(DWp_fact(k, 2)), SQRT(DWp_fact(k, 3))
          !
        ENDDO
      ENDIF
      !
      ! off-diagonal terms of tensor
      IF (ionode) WRITE(*,*) 
      IF (ionode) WRITE(*,*) "  off-diagonal terms"

      T_fact(:,:) = 0.d0
      Tp_fact(:,:) = 0.d0
      DO k = 1, nat
        DO p = 1, nq_tot
          T_fact(k, 1) =  T_fact(k, 1) + (D_tau(p, k, 1) - equil_p(p, k, 1)) * & 
                                         (D_tau(p, k, 2) - equil_p(p, k, 2)) / nq_tot
          T_fact(k, 2) =  T_fact(k, 2) + (D_tau(p, k, 1) - equil_p(p, k, 1)) * & 
                                         (D_tau(p, k, 3) - equil_p(p, k, 3)) / nq_tot
          T_fact(k, 3) =  T_fact(k, 3) + (D_tau(p, k, 2) - equil_p(p, k, 2)) * & 
                                         (D_tau(p, k, 3) - equil_p(p, k, 3)) / nq_tot
        ENDDO
      ENDDO
      IF (ionode) THEN 
        DO k = 1, nat
          nta = ityp(k)
          WRITE(*,'(A6, 3F11.6)') atm(nta), T_fact(k, 1), T_fact(k, 2), T_fact(k, 3) 
        ENDDO
      ENDIF ! (ionode)
      EXIT ! exit kk-loop if the error is less than a threshold
    ENDIF
  !
  ENDDO ! end kk for niters
  !
  IF (ionode) CLOSE(80) ! close ZG-configuration file
  IF (ionode) CLOSE(81) ! close ZG-velocities file
  IF (TRIM(flscf) /= ' ' .AND. ionode) CLOSE(82) ! close ZG-scf file
  !
  IF (ionode) WRITE(*,*) 
  IF ( ZG_strf .AND. ( SUM(ABS(ratio_zg)) / nat3 < error_thresh ) .AND. ionode ) & 
  WRITE(*,*) "Computing ZG structure factor ..."
  IF ( ZG_strf .AND. ( SUM(ABS(ratio_zg)) / nat3 < error_thresh) ) &
             call ZG_structure_factor(qpts_strf, D_tau, equil_p, nq_tot, &
                          nat, alat, ityp, ntypx, atmsf_a, atmsf_b, & 
                          nrots, kres1, kres2, kmin, kmax, col1, col2, Np)
  !
  IF (ionode) WRITE(*,*)
  IF (SUM(ABS(ratio_zg)) / nat3 > error_thresh .AND. ionode ) & 
      CALL errore('ZG', 'Error is not less than threshold', 5 )
  !
  ! Perform ASDM if ASDM is true
  IF (ionode .AND. ASDM) & 
    CALL ASDM_IFC(nat, nq_tot, D_tau, flscf, pt_T, ntyp, ntypx, atm, ityp, amass, & 
                  dim1, dim2, dim3, at, crystal_pos, alat, tau, fd_displ, & 
                  apply_fd, read_fd_forces, incl_epsil, iter_idx, iter_idx0, poly, & 
                  ibrav, epsil, zeu, poly_fd_forces, mixing, update_equil)
    
  ! 
  DEALLOCATE(T_fact, Tp_fact, DW_fact, DWp_fact)
  DEALLOCATE(equil_p, Rlist, D_tau, P_tau, qA, qB, z_nq_A, z_nq_B)
  DEALLOCATE(Cx_matA, Cx_matB, Cx_matAB)
  DEALLOCATE(Cpx_matA, Cpx_matB, Cpx_matAB)
  DEALLOCATE(Mx_mat_or, Mx_mat, M_mat, V_mat, ratio_zg)
  IF (compute_error) THEN
    DEALLOCATE(R_mat, D_vect, F_vect, Bx_vect, E_vect)
    DEALLOCATE(sum_er_D, sum_diag_D, sum_er_B, sum_diag_B, sum_er_B2, sum_diag_B2)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE ZG_configuration

SUBROUTINE ZG_structure_factor(qpts_strf, D_tau, equil_p, nq_tot, &
                               nat, alat, ityp, ntypx, atmsf_a, atmsf_b, &
                               nrots, kres1, kres2, kmin, kmax, col1, col2, Np)
!
 USE kinds,      ONLY : DP
 USE cell_base,  ONLY : bg
 USE constants,  ONLY : pi, tpi, BOHR_RADIUS_ANGS
 USE io_global,  ONLY : ionode, ionode_id, stdout
 USE mp_global,  ONLY : inter_pool_comm
 USE mp,         ONLY : mp_bcast, mp_barrier, mp_sum
 USE mp_world,   ONLY : world_comm
 !
 IMPLICIT NONE
 ! 
 INTEGER, INTENT(in)          :: nq_tot, nat, qpts_strf, ntypx
 INTEGER, INTENT(in)          :: nrots, kres1, kres2, col1, col2, Np
 REAL(DP), INTENT(in)         :: kmin, kmax
 REAL(DP), INTENT(in)         :: D_tau(nq_tot, nat, 3), equil_p(nq_tot, nat, 3), alat
 REAL(DP), INTENT(in)         :: atmsf_a(ntypx, 5), atmsf_b(ntypx, 5)
 INTEGER                      :: i, k, p, j, kk, pp, nta, ii
 INTEGER                      :: lower_bnd, upper_bnd, ctr
 INTEGER ityp(nat)
 REAL(DP)                     :: q_strf(3, qpts_strf), atomic_form_factor(nat, qpts_strf)
 REAL(DP)                     :: dotp, dotpp, eps
 REAL(DP), ALLOCATABLE        :: strf_rot(:, :)
 COMPLEX(DP)                  :: imagi, strf_map(qpts_strf)
 !
 eps   = 1d-5
 !
 q_strf = 0.d0
 IF (ionode) THEN
   OPEN (unit = 99, file = './qpts_strf.dat', status = 'unknown', form = 'formatted')
   !
   DO k = 1, qpts_strf
     READ(99,*) q_strf(:, k)
   !   WRITE(*,*) "aa", q_strf(k,:)
     CALL cryst_to_cart(1, q_strf(:, k), bg, +1)
     q_strf(:, k) = q_strf(:, k) * ( tpi / alat / BOHR_RADIUS_ANGS ) ! / 0.138933 * 0.073520
   ENDDO
   !
   CLOSE(99)
 ENDIF
 CALL mp_bcast(q_strf, ionode_id, world_comm)
 !
 !
 atomic_form_factor = 0.d0
 !
 DO k = 1, nat
   nta = ityp(k)
   DO i = 1, qpts_strf
      DO ii = 1, 5
        atomic_form_factor(k, i) = atomic_form_factor(k, i) + (atmsf_a(nta, ii) * &
                               EXP(-atmsf_b(nta, ii) * (NORM2(q_strf(:, i)) / 4.d0 / pi)**2))
      ENDDO
   ENDDO
 ENDDO
 !
 !
 CALL fkbounds( qpts_strf, lower_bnd, upper_bnd )
 strf_map = (0.d0, 0.d0)
 DO i = lower_bnd, upper_bnd !1, qpts_strf
   DO k = 1, nat
!     DO kk = 1, nat !! 
       DO p = 1, nq_tot
!         DO pp = 1, nq_tot  !!
           dotp = 0.0d0
           DO j = 1, 3
             dotp = dotp + q_strf(j, i) * D_tau(p, k, j)
!            dotp = dotp + q_strf(j, i) * (D_tau(p, k, j) - D_tau(pp, kk, j)) !! 
           ENDDO
           strf_map(i) = strf_map(i) + CMPLX(COS(dotp), SIN(dotp), kind=dp) * atomic_form_factor(k, i)
          !strf_map(i) = strf_map(i) + CMPLX(COS(dotp), SIN(dotp), kind=dp) * atomic_form_factor(k, i) * atomic_form_factor(kk, i)
!        ENDDO !!
      ENDDO ! p
!    ENDDO !!
  ENDDO ! k 
 ENDDO ! i
 CALL mp_sum(strf_map, inter_pool_comm)
 CALL mp_barrier(inter_pool_comm)
 !
 ! APPLY BROADENING and print outputs
 ctr = 0
 DO k = 1, qpts_strf
   IF ((q_strf(col1, k) .GT. 0.d0 - eps) .AND. & 
       (q_strf(col2, k) .GT. 0.d0 - eps)) THEN
      ctr = ctr + 1
   ENDIF
 ENDDO
 ALLOCATE(strf_rot(ctr * nrots, 4))
 !
 strf_rot = 0.d0
 IF (ionode) CALL rotate(DBLE(ABS(strf_map)**2), q_strf, qpts_strf, 0, nrots, &
                         ctr, strf_rot, col1, col2)
 !
 CALL mp_bcast(strf_rot, ionode_id, world_comm)
 CALL disca_broadening(strf_rot, ctr * nrots, kres1, kres2, alat, &
                       kmin, kmax, col1, col2, Np, 'strf_ZG_broad.dat')
 ! 
 DEALLOCATE(strf_rot)
 !
 IF (ionode) THEN
   OPEN (unit = 98, file = './structure_factor_ZG_raw.dat', status = 'unknown', form = 'formatted')
   !
   DO i = 1, qpts_strf
   !
     WRITE(98,'(4f26.6)') q_strf(:, i), ABS(strf_map(i))**2
   !abs(strf_map(i))**2  ! real(strf_map(i)) !
   !
   ENDDO
   !
   CLOSE(98)
 ENDIF
 !
 RETURN
 !
END SUBROUTINE

SUBROUTINE create_supercell(at, tau, alat, dim1, dim2, dim3, nat, equil_p, crystal_pos)
!
 USE kinds,      ONLY : DP
 USE constants,  ONLY : BOHR_RADIUS_ANGS
 USE cell_base,  ONLY : bg
 IMPLICIT NONE
!
!
 INTEGER,  INTENT(in)   :: dim1, dim2, dim3, nat
 REAL(DP), INTENT(in)   :: tau(3, nat), at(3, 3), alat
 REAL(DP), INTENT(out)  :: equil_p(dim1 * dim2 * dim3, nat, 3)
 INTEGER                :: i, j, k, ctr, p
 REAL(DP)               :: alat_ang, abc(3, nat)
 REAL(DP), INTENT(out)  :: crystal_pos(dim1 * dim2 * dim3, nat, 3)
 alat_ang = alat * BOHR_RADIUS_ANGS !bohr_to_angst ! to convert them in angstrom ! 
 abc = tau
 ! to convert tau/abc in crystal coordinates
 !
 call cryst_to_cart(nat, abc, bg, -1)
 !
 !
 crystal_pos = 0.d0
 ctr = 1
 DO i = 0, dim3 - 1
   DO j = 0, dim2 - 1
     DO k = 0, dim1 - 1
       DO p = 1, nat
         crystal_pos(ctr, p, 1) = (abc(1, p) + k) / FLOAT(dim1)
         crystal_pos(ctr, p, 2) = (abc(2, p) + j) / FLOAT(dim2)
         crystal_pos(ctr, p, 3) = (abc(3, p) + i) / FLOAT(dim3)
          ! WRITE(*, *) "cryspos:", crystal_pos(ctr, p, :)
       ENDDO
     ctr = ctr + 1
     ENDDO
   ENDDO
 ENDDO
 !
 !
 ctr = 1
 equil_p = 0.d0
 DO i = 0, dim3 - 1
   DO j = 0, dim2 - 1
     DO k = 0, dim1 - 1
       DO p = 1, nat
         equil_p(ctr, p, :) = (abc(1, p) + k) * at(:, 1) &
                            + (abc(2, p) + j) * at(:, 2) &
                            + (abc(3, p) + i) * at(:, 3)
       ENDDO
     ctr = ctr + 1
     ENDDO
   ENDDO
 ENDDO
 !
 equil_p = equil_p * alat_ang
 !
 RETURN
 !
 END SUBROUTINE create_supercell


 SUBROUTINE single_phonon(nq_tot, nat, ctrB, ctrA, nat3, ityp, ntyp, &
                          ntypx, qA, qB, amass, atm, equil_p, &
                          Rlist, z_nq_B, z_nq_A, Cx_matB, &
                          Cx_matA, Cpx_matB, Cpx_matA) 
 !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : AMU_SI, tpi
  IMPLICIT NONE
 !
 !
  CHARACTER(LEN=3), INTENT(in) :: atm(ntypx)
  INTEGER,  INTENT(in)         :: nq_tot, nat, ctrB, ctrA, nat3, ntyp, ntypx
  INTEGER,  INTENT(in)         :: Rlist(nq_tot, 3) 
  INTEGER,  INTENT(in)         :: ityp(nat)
  REAL(DP), INTENT(in)         :: qA(ctrA, 3), qB(ctrB, 3), amass(ntyp), equil_p(nq_tot, nat, 3)
  REAL(DP), INTENT(in)         :: Cx_matB(nat3, ctrB), Cx_matA(nat3, ctrA), Cpx_matA(nat3, ctrA), Cpx_matB(nat3, ctrB) 
  COMPLEX(DP), INTENT(in)      :: z_nq_A(nat3, nat3, ctrA), z_nq_B(nat3, nat3, ctrB)
  CHARACTER(len=256)           :: filename
 !
  INTEGER                      :: i, j, k, p, ii, ctr, nta, qp
  REAL(DP)                     :: dotp
  REAL(DP), ALLOCATABLE        :: D_tau(:, :, :), P_tau(:, :, :)
  COMPLEX(DP)                  :: imagi
  !
  !
  imagi = (0.0d0, 1.0d0) !imaginary unit
  !
  ALLOCATE(D_tau(nq_tot, nat, 3), P_tau(nq_tot, nat, 3))
  !
  ! 
  filename = 'single_phonon-displacements.dat' !'.fp'
  OPEN (unit = 85, file = filename, status = 'unknown', form = 'formatted')
  !
  filename = 'single_phonon-velocities.dat' !'.fp'
  OPEN (unit = 86, file = filename, status = 'unknown', form = 'formatted')
  ! Main loop to give single phonon displacements / velocities
  WRITE(85,'(A50)') "Displaced positions along phonon modes in set B"
   DO qp = 1, ctrB
     DO j = 1, nat3        
       WRITE(85,'(A25, 3F8.4, A15, I0)') "Phonon mode at q-point", qB(qp, :), " and branch:", j
       WRITE(86,'(A25, 3F8.4, A15, I0)') "Phonon mode at q-point", qB(qp, :), " and branch:", j
       D_tau = 0.0d0
       P_tau = 0.0d0
       DO p = 1, nq_tot
         dotp = 0.0d0
         DO ii = 1, 3
            dotp = dotp + qB(qp, ii) * Rlist(p, ii) * tpi! dot product between q and R 
         ENDDO
         ctr = 1
         DO k = 1, nat ! k represents the atom
           nta = ityp(k)
           DO i = 1, 3  ! i is for cart directions
            D_tau(p, k, i) = D_tau(p, k, i) + SQRT(2.0d0 / nq_tot / amass(nta)) * DBLE(CMPLX(COS(dotp), SIN(dotp), kind=dp) &
                                  * z_nq_B(ctr, j, qp) * (1.d0 + imagi) * ABS(Cx_matB(j, qp))) 
            P_tau(p, k, i) = P_tau(p, k, i) + SQRT(2.0d0 / nq_tot * amass(nta)) * DBLE(CMPLX(COS(dotp), SIN(dotp), kind=dp) &
                                  * z_nq_B(ctr, j, qp) * (1.d0 + imagi) * ABS(Cpx_matB(j, qp))) / (amass(nta) * AMU_SI)
            ! Here we calculate the momenta of the nuclei and finally 
            ! we divide by (amass(nta) *AMU_SI) to get the velocities.
            ctr = ctr + 1 ! for k and i
            !IF (ABS(D_tau(p, k, i)) .GT. 5) CALL errore('ZG', 'Displacement very large', D_tau(p, k, i) )
            IF (ABS(D_tau(p, k, i)) .GT. 5) CALL errore('ZG', 'Displacement very large', 1)
            D_tau(p, k, i) = equil_p(p, k, i) + D_tau(p, k, i) ! add equil structure
           ENDDO ! i loop
         ! write output data
         WRITE(85,'(A6, 3F13.8)') atm(ityp(k)), D_tau(p, k, :) 
         WRITE(86,'(A6, 3F15.8)') atm(ityp(k)), P_tau(p, k,:) * 1.0E-12 ! multiply to obtain picoseconds 
         !
         ENDDO ! k loop
       ENDDO ! p loop
     ENDDO ! j loop
   ENDDO ! qp loop
   !
   WRITE(85,'(A50)') "Displaced positions along phonon modes in set A"
   DO qp = 1, ctrA
     DO j = 1, nat3        
       WRITE(85,'(A25, 3F8.4, A15, I0)') "Phonon mode at q-point", qA(qp, :), " and branch:", j
       WRITE(86,'(A25, 3F8.4, A15, I0)') "Phonon mode at q-point", qA(qp, :), " and branch:", j
       D_tau = 0.0d0
       P_tau = 0.0d0
       DO p = 1, nq_tot
         dotp = 0.0d0
         DO ii = 1, 3
            dotp = dotp + qA(qp, ii) * Rlist(p, ii) * tpi ! dot product between q and R 
         ENDDO
         ctr = 1
         DO k = 1, nat ! k represents the atom
           nta = ityp(k)
           DO i = 1, 3  ! i is for cart directions
            D_tau(p, k, i) = D_tau(p, k, i) + SQRT(1.0d0 / nq_tot / amass(nta)) * DBLE(CMPLX(COS(dotp), SIN(dotp), kind=dp) &
                                  * z_nq_A(ctr, j, qp) * (1.d0 + imagi) * ABS(Cx_matA(j, qp))) 
            P_tau(p, k, i) = P_tau(p, k, i) + SQRT(1.0d0 / nq_tot * amass(nta)) * DBLE(CMPLX(COS(dotp), SIN(dotp), kind=dp) &
                                  * z_nq_A(ctr, j, qp) * (1.d0 + imagi) * ABS(Cpx_matA(j, qp))) / (amass(nta) * AMU_SI)
            ! Here we calculate the momenta of the nuclei and finally 
            !we divide by (amass(nta) *AMU_SI) to get the velocities.
            ctr = ctr + 1 ! for k and i
            !IF (ABS(D_tau(p, k, i)) .GT. 5) CALL errore('ZG', 'Displacement very large', D_tau(p, k, i) )
            IF (ABS(D_tau(p, k, i)) .GT. 5) CALL errore('ZG', 'Displacement very large', 1)
            D_tau(p, k, i) = equil_p(p, k, i) + D_tau(p, k, i) ! add equil structure
           ENDDO ! i loop
         ! write output data
         WRITE(85,'(A6, 3F13.8)') atm(ityp(k)), D_tau(p, k, :) 
         WRITE(86,'(A6, 3F15.8)') atm(ityp(k)), P_tau(p, k,:) * 1.0E-12 ! multiply to obtain picoseconds 
         !
         ENDDO ! k loop
       ENDDO ! p loop
     ENDDO ! j loop
   ENDDO ! qp loop
       !
 !      
 !
  DEALLOCATE(D_tau, P_tau)
  CLOSE(85)
  CLOSE(86)
 !
 !
 RETURN
 !
 END SUBROUTINE single_phonon
 !
SUBROUTINE fkbounds( nktot, lower_bnd, upper_bnd )
  !-----------------------------------------------------------------------
  !!
  !!   Subroutine from EPW finds the lower and upper bounds a k-grid in parallel
  !!
  !! @ Note: 
  !!    If you have 19 kpts and 2 pool, this routine will return
  !!    lower_bnd= 1 and upper_bnd=10 for the first pool
  !!    lower_bnd= 1 and upper_bnd=9 for the second pool
  !-----------------------------------------------------------------------
  !
  USE mp_global,    ONLY: my_pool_id, npool
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (in) :: nktot
  !! nktot k-points splited over pools
  INTEGER, INTENT (out) :: lower_bnd
  !! Lower kpt bounds for that image pool 
  INTEGER, INTENT (out) :: upper_bnd
  !! Upper kpt for that image pool
  !
#if defined(__MPI)
  !
  INTEGER :: nkl, nkr
  !
  ! find the bounds of k-dependent arrays in the parallel case
  ! number of kpoint blocks, kpoints per pool and reminder
  !
  nkl = nktot / npool
  nkr = nktot - nkl * npool
  !
  ! the reminder goes to the first nkr pools (0... nkr - 1)
  !
  IF (my_pool_id < nkr ) nkl = nkl + 1
  !
  ! the index of the first k point in this pool
  !
  lower_bnd = my_pool_id * nkl + 1
  IF ( my_pool_id >= nkr ) lower_bnd = my_pool_id * nkl + 1 + nkr
  !
  ! the index of the last k point in this pool
  !
  upper_bnd = lower_bnd + nkl - 1
  !
#else  
  !     
  ! In serial the definitions are much easier 
  !     
  lower_bnd = 1
  upper_bnd = nktot
  !     
#endif 
  !
  RETURN
  !
END SUBROUTINE fkbounds
!
SUBROUTINE rotate(strf, q, nq, nq_super, nrots, & 
                   ctr, strf_rot, col1, col2)
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in)    :: nq, nrots, col1, col2, ctr, nq_super
  REAL(DP), INTENT(in)   :: strf(nq), q(3, nq)
  REAL(DP), INTENT(out)  :: strf_rot(ctr * nrots, 4)
  INTEGER                :: i, p, ctr2
  REAL(DP)               :: str_f(ctr,4), str_fp(ctr, 4)
  !
  REAL(DP)               :: Rmat(2,2), theta, eps
  ! 
  eps        = 1d-5
  !
  !WRITE(*,*) "Number of rotations provided:", nrots
  theta      = tpi / FLOAT(nrots)
  !
  !
  str_f = 0.d0
  ctr2 = 0
  DO p = nq_super + 1, nq
    IF ((q(col1, p) .GT. 0.0 - eps) .AND. (q(col2, p) .GT. 0.0 - eps)) THEN
      ctr2 = ctr2 + 1
      str_f(ctr2, 1:3) = q(:, p)
      str_f(ctr2, 4) = strf(p)
    ENDIF
  ENDDO
  !
  ! To remove double contribution upon rotation
  !
  str_fp = 0.d0
  str_fp(1, :) = str_f(1, :)
  DO p = 2, ctr
    IF (ATAN(str_f(p, 1) / str_f(p, 2)) .LT. ( tpi / FLOAT(nrots) - eps) ) THEN
    str_fp(p, :) = str_f(p, :)
    ENDIF
  ENDDO
  !
  !
  strf_rot = 0.d0
  ctr2 = 1
  DO i = 0, nrots - 1
    Rmat(1, :) = (/ COS(i * theta), -SIN(i * theta) /)
    Rmat(2, :) = (/ SIN(i * theta),  COS(i * theta) /)
    DO p = 1, ctr
      strf_rot(ctr2, 1 : 2) = MATMUL(Rmat, str_fp(p, 1 : 2))
      strf_rot(ctr2, 3) = str_fp(p, 3)
      strf_rot(ctr2, 4) = str_fp(p, 4)
      ctr2 = ctr2 + 1
    END DO
  ENDDO
END SUBROUTINE
!
SUBROUTINE disca_broadening(strf_rot, steps, kres1, kres2, alat, &
                            kmin, kmax, col1, col2, Np, flstrfout)
!-------------------------------------------------------------------------
!! authors: Marios Zacharias, Feliciano Giustino 
!! acknowledgement: Hyungjun Lee for help packaging this release
!! version: v0.1
!! license: GNU
!
USE kinds,       ONLY : dp
USE mp_global,   ONLY : inter_pool_comm
USE mp_world,    ONLY : world_comm
USE mp,          ONLY : mp_bcast, mp_barrier, mp_sum
USE io_global,   ONLY : ionode, ionode_id, stdout
USE constants,   ONLY : tpi
!
IMPLICIT NONE
!
 CHARACTER(LEN=256), INTENT(IN)  :: flstrfout
 INTEGER, INTENT(IN)             :: steps, kres1, kres2, Np, col1, col2
 REAL(DP), INTENT(IN)            :: kmin, kmax, alat
 REAL(DP), INTENT(IN)            :: strf_rot(steps, 4)
 INTEGER                         :: ii, lower_bnd, upper_bnd, ik, iky
 REAL(DP)                        :: jump, sf_smearingx, sf_smearingy, maxv !, pi 
 REAL(DP), ALLOCATABLE           :: kgridx(:), kgridy(:), strf_out(:, :)
!
!
!
ALLOCATE( kgridx(kres1), kgridy(kres2))
ALLOCATE(strf_out(kres1,kres2))
!kmin = -10.0
!kmax = 10.0

jump = (kmax - kmin) / DBLE(kres1 - 1)
DO ik = 1, kres1
  kgridx(ik) = kmin + (ik - 1) * jump
ENDDO
sf_smearingx = (kmax - kmin) / DBLE(kres1)
!
jump = (kmax - kmin) / DBLE(kres2 - 1)
DO ik = 1, kres2
  kgridy(ik) = kmin + (ik- 1)*jump
ENDDO
sf_smearingy = (kmax - kmin) / DBLE(kres2)
!! 
!
strf_out = 0.d0
!
CALL fkbounds( steps, lower_bnd, upper_bnd )
!
DO ii = lower_bnd, upper_bnd
  DO ik = 1, kres1 !
    DO iky = 1, kres2
    !
    strf_out(ik, iky) =  strf_out(ik, iky) +  &
                          strf_rot(ii, 4) / sf_smearingx / SQRT(tpi) / sf_smearingy / SQRT(tpi) * &
                          (EXP(-(strf_rot(ii, col1) - kgridx(ik))**2.d0 / sf_smearingx**2.d0 / 2.d0))*&
                          (EXP(-(strf_rot(ii, col2) - kgridy(iky))**2.d0 / sf_smearingy**2.d0 / 2.d0))
    !
    ENDDO
  ENDDO
ENDDO
!
CALL mp_sum(strf_out, inter_pool_comm)
CALL mp_barrier(inter_pool_comm)
!
IF (ionode) THEN
  OPEN(46,FILE=flstrfout)
  maxv =  maxval(strf_out)
  WRITE(46,*) "#", maxv, maxval(strf_out)
  DO ik = 1, kres1
    DO iky = 1, kres2
      WRITE(46,'(3F28.12)') kgridx(ik), kgridy(iky), strf_out(ik,iky) * Np**(-2.0d0) ! 
                                         !Np**(-2.0d0) ! / maxv
    ENDDO
    WRITE(46,*)
  ENDDO
  CLOSE(46)
ENDIF
!
DEALLOCATE(strf_out, kgridx, kgridy)
!
RETURN
!
END SUBROUTINE
! 
!
  SUBROUTINE ASDM_IFC(nat, nq_tot, D_tau, flscf, pt_T, ntyp, ntypx, atm, ityp, amass, &
                      dim1, dim2, dim3, at, crystal_pos, alat, tau, fd_displ, & 
                      apply_fd, read_fd_forces, incl_epsil, iter_idx, iter_idx0, poly, &
                      ibrav, epsil, zeu, poly_fd_forces, mixing, update_equil)

  USE kinds,      ONLY : dp
  USE constants,  ONLY : amu_ry, BOHR_RADIUS_ANGS
  USE cell_base,  ONLY : celldm, bg
  USE read_input, ONLY : read_input_file
  USE input_parameters, ONLY : ecutwfc, calculation, prefix, restart_mode, &
                               pseudo_dir, outdir, diagonalization, mixing_mode, &
                               mixing_beta, conv_thr, atom_pfile, &
                               nk1, nk2, nk3, k1, k2, k3, & 
                               occupations, degauss, smearing 
  USE symm_base,  ONLY : sr, nsym, find_sym, find_sym_ifc, nrot, irt, & 
                         set_sym_bl ! for read_fd_forces
  USE noncollin_module, ONLY : noncolin, domag, m_loc ! for read_fd_forces
  
  IMPLICIT NONE
  ! 
  CHARACTER(LEN=3), INTENT(in) :: atm(ntypx)
  CHARACTER(LEN=256), INTENT(in) :: flscf, pt_T
  LOGICAL, INTENT(in)          :: read_fd_forces, incl_epsil, poly, poly_fd_forces
  LOGICAL, INTENT(in)          :: mixing, update_equil
  INTEGER, INTENT(in)          :: nq_tot, nat, ntyp, ntypx, iter_idx, iter_idx0
  INTEGER, INTENT(in)          :: dim1, dim2, dim3
  ! nq is the number of qpoints in sets A and B
  INTEGER, INTENT(in)          :: ityp(nat), ibrav ! for read_fd_forces
  REAL(DP), INTENT(in)         :: amass(ntyp), crystal_pos(nq_tot, nat, 3) 
  REAL(DP), INTENT(in)         :: epsil(3, 3), zeu(3, 3, nat), tau(3, nat), alat, fd_displ
  !
  CHARACTER (len=34)           :: A1 ! for read_fd_forces
  CHARACTER (len=*), PARAMETER :: search_str = "Forces" ! for read_fd_forces
  CHARACTER (len=*), PARAMETER :: search_str2 = "Begin" ! for read_fd_forces
  CHARACTER (len=50)           :: text ! for read_fd_forces
  CHARACTER (len=20)           :: word ! for read_fd_forces
  CHARACTER(len=256)           :: filename, pt_1, pt_2, pt_3
  CHARACTER (len=13)           :: filename2
  !
  LOGICAL                      :: apply_fd, file_exists, magnetic_sym  ! for read_fd_forces
  !
  ! local
  INTEGER, PARAMETER           :: iu = 84  ! for read_fd_forces
  INTEGER                      :: iosf, Rlist(nq_tot, 3), info ! for read_fd_forces
  INTEGER                      :: ctr, ctr2, i, j, p, k, ii, kk, qp
  INTEGER, ALLOCATABLE         :: ityp_super(:), irt_trans(:,:)  ! for read_fd_forces
  ! irt_trans --> atom index connected by translations, 
  ! irt_trans(dim1 * dim2 * dim3 * nat, dim1 * dim2 * dim3)
  REAL(DP)                     :: D_tau(nq_tot, nat, 3), P_tau(nq_tot, nat, 3)
  REAL(DP)                     :: at_or(3, 3), at(3, 3), atx(3, 3), inv_xkg(3, 3), alatx
  !atx --> supercell lattice param, for read_fd_forces
  REAL(DP), ALLOCATABLE        :: equil_p2(:, :), eqp_tmp(:, :) ! for read_fd_forces 
  REAL(DP), ALLOCATABLE        :: xkg(:, :, :, :), xkg_tmp(:, :), xkg_av(:, :, :, :) ! for read_fd_forces 
  REAL(DP), ALLOCATABLE        :: xkg_new(:, :, :, :, :), xkg_new_av(:, :, :, :) ! for read_fd_forces 
  REAL(DP), ALLOCATABLE        :: fd_forces(:, :, :), force_consts(:, :, :)! for read_fd_forces
  REAL(DP), ALLOCATABLE        :: force_consts_av(:, :, :), force_consts_iters(:, :, :, :) ! for read_fd_forces
  REAL(DP), ALLOCATABLE        :: eqf(:, :), eqf_tmp(:), eqf_av(:, :) ! for update_equil (rotations)
  REAL(DP), ALLOCATABLE        :: eqf_new(:, :, :), eqf_new_av(:, :) ! for update_equil (translations)
  REAL(DP), ALLOCATABLE        :: eqf_iters(:, :, :), eqf_av_iters(:, :), new_tau(:, :) ! for update_equil 
! 
  REAL(DP), PARAMETER      :: eps = 1.0d-6, eps2 = 1.0d-15
!
  WRITE(pt_2,'(1I2.2)') iter_idx 
  ! 
  ! Generate file for initial relaxation
  !
  
  filename2 = 'ZG-relax.out'
  INQUIRE(FILE=filename2, EXIST=file_exists)
  !
  IF (poly .AND. .NOT. file_exists) THEN 
    filename = 'ZG-relax.in'
    OPEN (unit = 82, file = filename, status = 'unknown', form = 'formatted')
    !IF ( TRIM(flscf) == ' ') CALL errore('ZG', 'No flscf specified', flscf )
    IF ( TRIM(flscf) == ' ') CALL errore('ZG', 'No flscf specified', 1)
    IF ( TRIM(flscf) /= ' ') THEN
      WRITE(82,*) "&control"
      WRITE(82,'(100A)') "  calculation = 'relax'"
      WRITE(82,'(100A)') "  restart_mode = '", TRIM(restart_mode),"'"
      WRITE(82,'(100A)') "  prefix = 'ZG-relax'"
      WRITE(82,'(100A)') "  pseudo_dir = '", TRIM(pseudo_dir),"'"
      WRITE(82,'(100A)') "  outdir = '", TRIM(outdir),"'"
      WRITE(82,'(100A)') "  tprnfor=.true."
      WRITE(82,'(100A)') "  disk_io = 'none'"
      WRITE(82,'(100A)') "  forc_conv_thr = 1.0D-5"
      WRITE(82,'(100A)') "  nstep = 300"
      WRITE(82,'(100A)') "/"
      !
      WRITE(82,'(100A)') "&system"
      WRITE(82,'(100A)') "  ibrav = 0"
      WRITE(82,'(A7,1i5)') "  nat =", nat * nq_tot
      WRITE(82,'(A8,1i5)') "  ntyp =", ntyp
      WRITE(82,'(A11,1F7.2)') "  ecutwfc =", ecutwfc
      WRITE(82,'(100A)') "  occupations = '", TRIM(occupations),"'"
      WRITE(82,'(100A)') "  smearing = '", TRIM(smearing),"'"
      WRITE(82,'(A11,1D10.1)') "  degauss =", degauss
      WRITE(82,'(100A)') "/"
      !
      WRITE(82,'(100A)') "&electrons"
      WRITE(82,'(100A)') "  diagonalization = '",TRIM(diagonalization),"'"
      WRITE(82,'(100A)') "  mixing_mode= '",TRIM(mixing_mode),"'"
      WRITE(82,'(A16,1F4.2)') "  mixing_beta = ", mixing_beta
      WRITE(82,'(A12,1D10.1)') "  conv_thr = ", conv_thr
      WRITE(82,'(100A)') "/"
      !
      WRITE(82,'(100A)') "&ions"
      WRITE(82,'(A23,1D10.1)') "  ion_dynamics = 'bfgs'"
      WRITE(82,'(100A)') "/"
      !
      WRITE(82,'(100A)') "ATOMIC_SPECIES"
      DO kk = 1, ntyp ! type of atom
        WRITE(82,'(A6, 1F8.3, A, A)') atm(kk), amass(kk), ' ', TRIM(atom_pfile(kk))
      ENDDO
      !
      WRITE(82,'(100A)') "K_POINTS automatic"
      WRITE(82,'(6i4)') NINT(nk1/DBLE(dim1)), NINT(nk2/DBLE(dim2)), &
                        NINT(nk3/DBLE(dim3)), k1, k2, k3
      !
      WRITE(82,'(100A)') "CELL_PARAMETERS (angstrom)"
      WRITE(82,'(3F16.8)') at(:, 1) * alat * dim1 * BOHR_RADIUS_ANGS
      WRITE(82,'(3F16.8)') at(:, 2) * alat * dim2 * BOHR_RADIUS_ANGS
      WRITE(82,'(3F16.8)') at(:, 3) * alat * dim3 * BOHR_RADIUS_ANGS
      WRITE(82,'(100A)') "ATOMIC_POSITIONS (angstrom)"
      !
      DO kk = 1, nat
        DO qp = 1, nq_tot
           WRITE(82,'(A6, 3F16.8)') atm(ityp(kk)), D_tau(qp, kk, :) 
        ENDDO
      ENDDO
      !
      CLOSE(82)
    ENDIF ! flscf
  ENDIF ! poly
  ! 
  ! apply fd for polymorphous configuration
  IF (poly .AND. file_exists) THEN
    D_tau = 0.d0
    WRITE(*, *) "  ======================================"
    WRITE(*, *) "  WARNING: ", filename2, "exists, code applies fd on polymorphous network ..."
    WRITE(*, *) "  Make sure that this is what you want otherwise remove file ZG-relax.out"
    WRITE(*, *) "  ======================================"
    WRITE(*, *) 
    OPEN (unit = iu, file = filename2, status = 'old', form = 'formatted')
    DO
      READ (iu, *, iostat = iosf) text ! iosf to check if we reach the
                                      ! end of file without finding forces
      READ (text, *) word 
      IF (word == search_str2) THEN
        WRITE(*, *) "  ", search_str2, " found in ", filename2
        WRITE(*, *) "  ======================================"
        READ (iu, *) ! skip empty lines
        READ (iu, *) ! skip empty lines
        DO kk = 1, nat
          DO qp = 1, nq_tot
            READ(iu,'(A3,3X,3F20.10)') A1, D_tau(qp, kk, :)
            ! Format should be the same with output_tau.f90 in PW/src
            !WRITE( stdout,'(A3,3X,3F20.10)') atm(ityp(na)), tau(:,na)*rescale
          ENDDO ! qp
        ENDDO ! kk
        EXIT ! exit do loop                
      ELSE
      IF (iosf /= 0) CALL errore('A-ZG',& 
                               'Final coords in file are missing ...', 1)
      ENDIF ! word
      ! 
    ENDDO ! do loop
    CLOSE(unit = iu)
  ! 
    apply_fd = .true.
    IF (poly_fd_forces) apply_fd = .false.
  !
  ENDIF ! poly .AND. file_exists
  !  
  ! Finite difference input files
  ! use P_tau to save D_tau
  P_tau = 0.d0
  IF (apply_fd) THEN
    ! Generate files with displaced coordinates
    ctr = 1
    DO k = 1, nat
      DO p = 1, nq_tot
        DO j = 1, 3
          DO i = 1, 2 ! for positive or negative displacements
            P_tau = D_tau
            IF (i == 1) P_tau(p, k, j) = P_tau(p, k, j) + fd_displ ! in Ang (default 0.02 Bohr)
            IF (i == 2) P_tau(p, k, j) = P_tau(p, k, j) - fd_displ ! is the finite displacement
            WRITE(pt_1,'(1I4.4)') ctr
            !IF ( TRIM(flscf) == ' ') CALL errore('ZG', 'No flscf specified', flscf )
            IF ( TRIM(flscf) == ' ') CALL errore('ZG', 'No flscf specified', 1)
            IF ( TRIM(flscf) /= ' ') THEN
              filename = 'ZG-scf_' // TRIM( pt_T ) // 'K_iter_' // TRIM(pt_2) // '_' // TRIM(pt_1) // '.in'
              IF (poly) filename = 'ZG-scf_poly_iter_' // TRIM(pt_2) // '_' // TRIM(pt_1) // '.in'
              OPEN (unit = 82, file = filename, status = 'unknown', form = 'formatted')
              WRITE(82,*) "&control"
              WRITE(82,'(100A)') "  calculation = '", TRIM(calculation),"'"
              WRITE(82,'(100A)') "  restart_mode = '", TRIM(restart_mode),"'"
              WRITE(82,'(100A)') "  prefix = 'ZG-", TRIM(prefix),"'"
              WRITE(82,'(100A)') "  pseudo_dir = '", TRIM(pseudo_dir),"'"
              WRITE(82,'(100A)') "  outdir = '", TRIM(outdir),"'"
              WRITE(82,'(100A)') "  tprnfor=.true."
              WRITE(82,'(100A)') "  disk_io = 'none'"
              WRITE(82,'(100A)') "/"
              WRITE(82,'(100A)') "&system"
              WRITE(82,'(100A)') "  ibrav = 0"
              WRITE(82,'(A7,1i5)') "  nat =", nat * nq_tot
              WRITE(82,'(A8,1i5)') "  ntyp =", ntyp
              WRITE(82,'(A11,1F7.2)') "  ecutwfc =", ecutwfc
              WRITE(82,'(100A)') "  occupations = '", TRIM(occupations),"'"
              WRITE(82,'(100A)') "  smearing = '", TRIM(smearing),"'"
              WRITE(82,'(A11,1D10.1)') "  degauss =", degauss
              WRITE(82,'(100A)') "/"
              WRITE(82,'(100A)') "&electrons"
              WRITE(82,'(100A)') "  diagonalization = '",TRIM(diagonalization),"'"
              WRITE(82,'(100A)') "  mixing_mode= '",TRIM(mixing_mode),"'"
              WRITE(82,'(A16,1F4.2)') "  mixing_beta = ", mixing_beta
              WRITE(82,'(A12,1D10.1)') "  conv_thr = ", conv_thr
              WRITE(82,'(100A)') "/"
              !
              WRITE(82,'(100A)') "ATOMIC_SPECIES"
              DO kk = 1, ntyp ! type of atom
                WRITE(82,'(A6, 1F8.3, A, A)') atm(kk), amass(kk), ' ', TRIM(atom_pfile(kk))
              ENDDO
              !
              WRITE(82,'(100A)') "K_POINTS automatic"
              WRITE(82,'(6i4)') NINT(nk1/DBLE(dim1)), NINT(nk2/DBLE(dim2)), &
                                NINT(nk3/DBLE(dim3)), k1, k2, k3
              !
              WRITE(82,'(100A)') "CELL_PARAMETERS (angstrom)"
              WRITE(82,'(3F16.8)') at(:, 1) * alat * dim1 * BOHR_RADIUS_ANGS
              WRITE(82,'(3F16.8)') at(:, 2) * alat * dim2 * BOHR_RADIUS_ANGS
              WRITE(82,'(3F16.8)') at(:, 3) * alat * dim3 * BOHR_RADIUS_ANGS
              WRITE(82,'(100A)') "ATOMIC_POSITIONS (angstrom)"
              !
              DO kk = 1, nat
                DO qp = 1, nq_tot
                   WRITE(82,'(A6, 3F16.8)') atm(ityp(kk)), P_tau(qp, kk, :) 
                ENDDO
              ENDDO
              !
              ctr = ctr + 1 
              !
              CLOSE(82)
            ENDIF ! flscf
          ENDDO ! i
        ENDDO ! j
      ENDDO ! k
    ENDDO ! p
    ! 
    ! Generate also the file without displacements (i.e. the ZG-conf)
    P_tau = D_tau
    !
    filename = 'ZG-scf_' // TRIM( pt_T ) // 'K_iter_' // TRIM(pt_2) // '.in'
    IF (poly) filename = 'ZG-scf_poly_iter_' // TRIM(pt_2) // '.in'
    OPEN (unit = 82, file = filename, status = 'unknown', form = 'formatted')
    WRITE(82,*) "&control"
    WRITE(82,'(100A)') "  calculation = '", TRIM(calculation),"'"
    WRITE(82,'(100A)') "  restart_mode = '", TRIM(restart_mode),"'"
    WRITE(82,'(100A)') "  prefix = 'ZG-", TRIM(prefix),"'"
    WRITE(82,'(100A)') "  pseudo_dir = '", TRIM(pseudo_dir),"'"
    WRITE(82,'(100A)') "  outdir = '", TRIM(outdir),"'"
    WRITE(82,'(100A)') "  tprnfor=.true."
    WRITE(82,'(100A)') "  disk_io = 'none'"
    WRITE(82,'(100A)') "/"
    WRITE(82,'(100A)') "&system"
    WRITE(82,'(100A)') "  ibrav = 0"
    WRITE(82,'(A7,1i5)') "  nat =", nat * nq_tot
    WRITE(82,'(A8,1i5)') "  ntyp =", ntyp
    WRITE(82,'(A11,1F7.2)') "  ecutwfc =", ecutwfc
    WRITE(82,'(100A)') "  occupations = '", TRIM(occupations),"'"
    WRITE(82,'(100A)') "  smearing = '", TRIM(smearing),"'"
    WRITE(82,'(A11,1D10.1)') "  degauss =", degauss
    WRITE(82,'(100A)') "/"
    WRITE(82,'(100A)') "&electrons"
    WRITE(82,'(100A)') "  diagonalization = '",TRIM(diagonalization),"'"
    WRITE(82,'(100A)') "  mixing_mode= '",TRIM(mixing_mode),"'"
    WRITE(82,'(A16,1F4.2)') "  mixing_beta = ", mixing_beta
    WRITE(82,'(A12,1D10.1)') "  conv_thr = ", conv_thr
    WRITE(82,'(100A)') "/"
    !
    WRITE(82,'(100A)') "ATOMIC_SPECIES"
    DO kk = 1, ntyp ! type of atom
      WRITE(82,'(A6, 1F8.3, A, A)') atm(kk), amass(kk), ' ', TRIM(atom_pfile(kk))
    ENDDO
    !
    WRITE(82,'(100A)') "K_POINTS automatic"
    WRITE(82,'(6i4)') NINT(nk1/DBLE(dim1)), NINT(nk2/DBLE(dim2)), &
                      NINT(nk3/DBLE(dim3)), k1, k2, k3
    !
    WRITE(82,'(100A)') "CELL_PARAMETERS (angstrom)"
    WRITE(82,'(3F16.8)') at(:, 1) * alat * dim1 * BOHR_RADIUS_ANGS
    WRITE(82,'(3F16.8)') at(:, 2) * alat * dim2 * BOHR_RADIUS_ANGS
    WRITE(82,'(3F16.8)') at(:, 3) * alat * dim3 * BOHR_RADIUS_ANGS
    WRITE(82,'(100A)') "ATOMIC_POSITIONS (angstrom)"
    !
    DO kk = 1, nat
      DO qp = 1, nq_tot
         WRITE(82,'(A6, 3F16.8)') atm(ityp(kk)), P_tau(qp, kk, :) 
      ENDDO
    ENDDO
    !
    !
    CLOSE(82)
    !
  ENDIF ! apply_fd
  !
  ! Read fd_forces in polymorphous structure
  !
  IF (poly_fd_forces) THEN
    ALLOCATE(fd_forces(nat * nq_tot, 3, nat * nq_tot * 2 * 3))
    ctr = 1
    DO k = 1, nat
      DO p = 1, nq_tot
        DO j = 1, 3
          DO i = 1, 2 ! for positive or negative displacements
            WRITE(pt_1,'(1I4.4)') ctr
           ! filename = 'fd_forces' // '/' // 'ZG-scf_' // TRIM( pt_T ) // 'K_' // TRIM(pt_1) // '.out'
            filename = 'fd_forces' // '/' // 'ZG-scf_poly_iter_' // & 
                        TRIM(pt_2) // '_' // TRIM(pt_1) // '.out'
            INQUIRE(FILE = filename, EXIST = file_exists)
            IF (.NOT. file_exists) CALL errore('ZG', 'Missing poly file in fd_forces/ ...', ctr )
            OPEN (unit = iu, file = filename, status = 'old', form = 'formatted')
            DO
              READ (iu, *, iostat = iosf) text ! iosf to check if we reach the
                                              ! end of file without finding forces
              READ (text, *) word 
              IF (word == search_str) THEN
                !WRITE(*,*) search_str, " found in ", filename
                READ (iu, *) ! skip empty line
                DO ii = 1, nat * nq_tot
                  READ(iu,'(A34, 3F14.8)') A1, fd_forces(ii, :, ctr)
                  ! Format should be the same with forces.f90 in PW/src
                  !WRITE(*, '(A34, 3F14.8)') A1, fd_forces(ii, :, ctr)
                ENDDO ! ii
                EXIT ! exit do loop                
              ELSE 
                IF (iosf /= 0) CALL errore('ZG',& 
                                      'forces in file are missing ...', ctr)
              ENDIF ! word
              ! 
            ENDDO ! do loop
            CLOSE(unit = iu)
            ctr = ctr + 1 ! to go to the next file
          ENDDO ! DO i 
        ENDDO ! DO j
      ENDDO ! DO p
    ENDDO ! DO k
    ! 
    ! Now perform finite differences to get the force constants
    ALLOCATE(force_consts(nat * nq_tot, 3, nat * nq_tot * 3))
    ctr = 1
    ctr2 = 1
    DO k = 1, nat
      DO p = 1, nq_tot
        DO j = 1, 3
          force_consts(:, :, ctr2) = (fd_forces(:, :, ctr + 1) - fd_forces(:, :, ctr)) & 
                                     / (2.0d0 * fd_displ / BOHR_RADIUS_ANGS)
          ! 13.6057039763 / BOHR_RADIUS_ANGS convert to eV / Ang
          ctr = ctr + 2
          ctr2 = ctr2 + 1
        ENDDO ! j
      ENDDO ! p
    ENDDO ! k
    !
    ! xkg: to convert forces into a format for performing mat multipl.
    ALLOCATE(xkg(3, 3,  nat * nq_tot,  nat * nq_tot))
    xkg = 0.d0
    !
    !WRITE(pt_1,'(1I2.2)') iter_idx 
    ! to print them in a format similar to phonopy
    filename = 'FORCE_CONSTANTS_poly_iter_' // TRIM(pt_2)
    OPEN (unit = 85, file = filename, status = 'unknown', form = 'formatted')
    j = 1
    ctr2 = 1
    WRITE(85,'(2I4)') nat * nq_tot, nat * nq_tot
    DO WHILE (j <= nat * nq_tot * 3) 
      DO p = 1, nat * nq_tot ! atom index
        !
        WRITE(85,'(2I4)') ctr2, p
        WRITE(85,'(3F14.8)')  force_consts(p, :, j)
        WRITE(85,'(3F14.8)')  force_consts(p, :, j + 1)
        WRITE(85,'(3F14.8)')  force_consts(p, :, j + 2)
        
        xkg(1, :, ctr2, p) =  force_consts(p, :, j)
        xkg(2, :, ctr2, p) =  force_consts(p, :, j + 1)
        xkg(3, :, ctr2, p) =  force_consts(p, :, j + 2)
        !
      ENDDO ! p
      j = j + 3
      ctr2 = ctr2 + 1
    ENDDO ! j
    CLOSE(85)
    !
    ! 
    ! Find symmetries using equilibirum structure
    ! convert first equil_pos in required format
    ALLOCATE(equil_p2(3, nat * nq_tot), ityp_super(nat * nq_tot))
    ! 
    at_or = at ! need to save original at to print the unit-cell fc info below
    !
    atx(:, 1) = dim1 * at(:, 1)
    atx(:, 2) = dim2 * at(:, 2)
    atx(:, 3) = dim3 * at(:, 3)
    ! alat = SQRT ( at(1,1)**2+at(2,1)**2+at(3,1)**2 ) as in cell_base.f90
    alatx = SQRT ( atx(1, 1)**2 + atx(2, 1)**2 + atx(3, 1)**2 )
    atx(:, :) = atx(:, :) / alatx
    ! 
    at = atx ! at accounts for supercell now
    !
    CALL recips( at(1, 1), at(1, 2), at(1, 3), bg(1, 1), bg(1, 2), bg(1, 3) )
    !
    ! generate Cart. coordinates in alat
    ctr = 1
    DO kk = 1, nat
      DO qp = 1, nq_tot
         !equil_p2(:, ctr) = crystal_pos(qp, kk, :) 
         equil_p2(:, ctr) = crystal_pos(qp, kk, 1) * at(:, 1)  & 
                          + crystal_pos(qp, kk, 2) * at(:, 2)  & 
                          + crystal_pos(qp, kk, 3) * at(:, 3) 
          !equil_p2 = nint(dble(equil_p2) * 1e6) / 1e6
          WRITE(*,'(A20, 3F14.8)' ) "Equil alat", equil_p2(:, ctr)
         ityp_super(ctr) = ityp(kk)
         ctr = ctr + 1
      ENDDO
    ENDDO
    !
    ! irt is the index mat that maps the atoms after symmetry oper. 
    IF ( ALLOCATED( irt ) ) DEALLOCATE( irt )
    ALLOCATE( irt(48, nat * nq_tot)) 
    !
    CALL set_sym_bl ( ) ! generate sym. environment
    !
    magnetic_sym = noncolin .AND. domag
    CALL find_sym(nat * nq_tot, equil_p2, ityp_super, magnetic_sym, m_loc)
    !
    !
    WRITE(*, *) "  ======================================"
    CALL print_symmetries ( 1, noncolin, domag )
    WRITE(*, *) "  ======================================"
    WRITE(*, *) 
    !
    ! Now apply symmetries on force constants.
    ! WRITE(*,*) "Num of sym and rots", nsym, nrot
    !
    ALLOCATE(xkg_av(3, 3, nat * nq_tot, nat * nq_tot), xkg_tmp(3, 3))
    xkg_av = 0.d0
    xkg_tmp = 0.d0
    DO k = 1, nsym
      DO i = 1, nat * nq_tot
        DO j = 1, nat * nq_tot
          xkg_tmp(:, :) =  MATMUL(TRANSPOSE(sr(:, :, k)), MATMUL(xkg(:, :, irt(k, i), irt(k, j)), sr(:, :, k)))
          xkg_av(:, :, i, j) = xkg_av(:, :, i, j) + xkg_tmp(:, :) / DBLE(nsym)
        ENDDO
      ENDDO
    ENDDO
    !
    ! now check translations
    !
    ! convert to crystal coordinates
    ctr = 1
    DO kk = 1, nat
      DO qp = 1, nq_tot
         equil_p2(:, ctr) = crystal_pos(qp, kk, :) 
         IF (equil_p2(1, ctr) < 0.d0) equil_p2(1, ctr) = equil_p2(1, ctr) + 1 ! only positive coord.
         IF (equil_p2(2, ctr) < 0.d0) equil_p2(2, ctr) = equil_p2(2, ctr) + 1 ! only positive coord.
         IF (equil_p2(3, ctr) < 0.d0) equil_p2(3, ctr) = equil_p2(3, ctr) + 1 ! only positive coord.
         ctr = ctr + 1
      ENDDO
    ENDDO
    !
    ! Apply transl and fill irt_trans
    ALLOCATE(eqp_tmp(3, nat * nq_tot), irt_trans(nat * nq_tot, nq_tot))
    !
    ctr = 1
    DO i = 0, dim1 - 1
      DO j = 0, dim2 - 1
        DO k = 0, dim3 - 1
          DO kk = 1, nat * nq_tot
              eqp_tmp(1, kk) = equil_p2(1, kk) - i / DBLE(dim1)
              IF (eqp_tmp(1, kk) < 0.d0) eqp_tmp(1, kk) = eqp_tmp(1, kk) + 1
              !
              eqp_tmp(2, kk) = equil_p2(2, kk) - j / DBLE(dim2)
              IF (eqp_tmp(2, kk) < 0.d0) eqp_tmp(2, kk) = eqp_tmp(2, kk) + 1
              !
              eqp_tmp(3, kk) = equil_p2(3, kk) - k / DBLE(dim3)
              IF (eqp_tmp(3, kk) < 0.d0) eqp_tmp(3, kk) = eqp_tmp(3, kk) + 1
          ENDDO ! kk
          ! 
          ! check distances if equal
          DO ii = 1, nat * nq_tot
            DO kk = 1, nat * nq_tot
              IF (NORM2(eqp_tmp(:, ii) - equil_p2(:, kk)) < eps ) irt_trans(ii, ctr) = kk
            ENDDO ! kk
          ENDDO ! ii
          ctr = ctr + 1
        ENDDO ! k
      ENDDO ! j
    ENDDO ! i
    ! 
    WRITE(*,*)  
    WRITE(*,*) "  Atom indices connected by translations"
    WRITE(*,*) "  ======================================"
    DO kk = 1, nat * nq_tot
      WRITE(*,'(1000I4)') irt_trans(kk, :)
    ENDDO
    !
    ! Now take average over force constants related by translation
    ALLOCATE(xkg_new(3, 3, nat * nq_tot, nat * nq_tot, nq_tot), xkg_new_av(3, 3, nat * nq_tot, nat * nq_tot))
    !
    xkg_new = 0.d0
    DO k = 1, nq_tot
      DO i = 1, nat * nq_tot
        DO j = 1, nat * nq_tot
          xkg_new(:, :, irt_trans(i, k), irt_trans(j, k), k) = xkg_av(:, :, i, j)
         ENDDO
      ENDDO
    ENDDO
    !
    xkg_new_av = 0.d0
    DO i = 1, nq_tot
      xkg_new_av(:, :, :, :) = xkg_new_av(:, :, :, :) + DBLE(xkg_new(:, :, :, :, i) / nq_tot)
    ENDDO
    !
    ! Impose hermiticity
    !!DO i = 1, nat * nq_tot
    !!  DO j = 1, nat * nq_tot
    !!    DO k = 1, 3
    !!      !DO kk = 1, 3
    !!      kk = 1
    !!      DO WHILE (kk < k) 
    !!        xkg_new_av(k, kk, i, j) = (xkg_new_av(k, kk, i, j) + xkg_new_av(kk, k, i, j))/2
    !!        xkg_new_av(kk, k, i, j) = xkg_new_av(k, kk, i, j) 
    !!        kk = kk + 1
    !!      ENDDO
    !!    ENDDO
    !!  ENDDO
    !!ENDDO
    ! 
    ! 
    ! Print them in phonopy format
    !WRITE(pt_1,'(1I2.2)') iter_idx 
    filename = 'FORCE_CONSTANTS_sym_poly_iter_' // TRIM(pt_2)
    OPEN (unit = 85, file = filename, status = 'unknown', form = 'formatted')
    WRITE(85, '(2i4)')  nat * nq_tot, nat * nq_tot 
    DO i = 1, nat * nq_tot
      DO j = 1, nat * nq_tot
      !
      WRITE(85, '(2i4)') i, j
      WRITE(85, '(3F14.8)') xkg_new_av(1, :, i, j)
      WRITE(85, '(3F14.8)') xkg_new_av(2, :, i, j)
      WRITE(85, '(3F14.8)') xkg_new_av(3, :, i, j)
      !
      ENDDO
    ENDDO
    CLOSE(85)
    ! To print in q2r.fc format
    !
    ctr2 = 1
    DO i = 1, dim3 
      DO j = 1, dim2 
        DO  k = 1, dim1
          Rlist(ctr2, 1) = k
          Rlist(ctr2, 2) = j
          Rlist(ctr2, 3) = i
          ctr2 = ctr2 + 1
        ENDDO
      ENDDO
    ENDDO
    !
    filename = 'poly_iter_' // TRIM(pt_2) // '.fc'
    OPEN (unit = 85, file = filename, status = 'unknown', form = 'formatted')
    !copy from do_q2r routine
    !Write initial part of fc file (info of the cell, dielectric const, BEC
    WRITE(85,'(i3,i5,i4,6f11.7)') ntyp, nat, ibrav, celldm
    IF (ibrav==0) then
      WRITE (85,'(2x,3f15.9)') ((at_or(i, j), i = 1, 3),j = 1, 3)
    ENDIF
    !
    DO kk = 1, ntyp
      WRITE(85, *) kk, " '", atm(kk), "' ", amass(kk) * amu_ry
    ENDDO
    !
    DO kk = 1, nat
      WRITE(85, '(2i5,3f18.10)') kk, ityp(kk), (tau(j, kk),j = 1, 3)
    ENDDO
    !
    WRITE (85, *) incl_epsil
    IF (incl_epsil) THEN
      WRITE(85, '(3f24.12)') ((epsil(i, j), j = 1, 3), i = 1, 3)
      DO kk = 1, nat
        WRITE(85, '(i5)') kk
        WRITE(85, '(3f15.7)') ((zeu(i, j, kk), j = 1, 3), i = 1, 3)
      ENDDO
    ENDIF ! incl_epsil
    WRITE(85, '(3I4)') dim1, dim2, dim3
    !
    ! WRITE FORCE CONSTANTS
    DO i = 1, 3
      DO j = 1, 3
        DO ii = 1, nat
          ctr = 1
          DO kk = 1, nat 
            WRITE(85, '(4I4)') i, j, ii, kk
            DO k = 1, nq_tot
              WRITE(85, '(3I3, 1F14.8)') Rlist(k, :), xkg_new_av(i, j, (ii - 1) * nq_tot + 1, ctr) 
              ctr = ctr + 1
            ENDDO ! k
          ENDDO ! kk
        ENDDO ! ii        
      ENDDO ! j
    ENDDO ! i
    CLOSE(85)
  !
    DEALLOCATE(fd_forces, force_consts, equil_p2, ityp_super)
    DEALLOCATE(xkg, xkg_av, xkg_tmp, eqp_tmp)
    DEALLOCATE(xkg_new, xkg_new_av)
  !
  ENDIF ! poly_fd_forces
  !
  ! Read fd_forces
  !
  IF (read_fd_forces) THEN
    ALLOCATE(fd_forces(nat * nq_tot, 3, nat * nq_tot * 2 * 3))
    ctr = 1
    DO k = 1, nat
      DO p = 1, nq_tot
        DO j = 1, 3
          DO i = 1, 2 ! for positive or negative displacements
            WRITE(pt_1,'(1I4.4)') ctr
            filename = 'fd_forces' // '/' // 'ZG-scf_' // TRIM( pt_T ) & 
                        // 'K_iter_' // TRIM(pt_2) // '_' // TRIM(pt_1) // '.out'
            !filename = 'ZG-scf_' // TRIM( pt_T ) // 'K_iter_' // TRIM(pt_2) // '_' // TRIM(pt_1) // '.in'
            INQUIRE(FILE = filename, EXIST = file_exists)
            IF (.NOT. file_exists) CALL errore('ZG', 'Missing file in fd_forces/ ...', ctr )
            OPEN (unit = iu, file = filename, status = 'old', form = 'formatted')
            DO
              READ (iu, *, iostat = iosf) text ! iosf to check if we reach the
                                              ! end of file without finding forces
              READ (text, *) word 
              IF (word == search_str) THEN
                !WRITE(*,*) search_str, " found in ", filename
                READ (iu, *) ! skip empty line
                DO ii = 1, nat * nq_tot
                  READ(iu,'(A34, 3F14.8)') A1, fd_forces(ii, :, ctr)
                  ! Format should be the same with forces.f90 in PW/src
                  !WRITE(*, '(A34, 3F14.8)') A1, fd_forces(ii, :, ctr)
                ENDDO ! ii
                EXIT ! exit do loop                
              ! 
              ELSE 
                IF (iosf /= 0) CALL errore('ZG',& 
                                      'forces in file are missing ...', ctr)
              ENDIF ! word
            ENDDO ! do loop
            CLOSE(unit = iu)
            ctr = ctr + 1 ! to go to the next file
          ENDDO ! DO i 
        ENDDO ! DO j
      ENDDO ! DO p
    ENDDO ! DO k
    ! 
    ! Now perform finite differences to get the force constants
    ALLOCATE(force_consts(nat * nq_tot, 3, nat * nq_tot * 3)) 
    ctr = 1
    ctr2 = 1
    DO k = 1, nat
      DO p = 1, nq_tot
        DO j = 1, 3
          force_consts(:, :, ctr2) = (fd_forces(:, :, ctr + 1) - fd_forces(:, :, ctr)) & 
                                     / (2.0d0 * fd_displ / BOHR_RADIUS_ANGS)
          ! I use BOHR_RADIUS_ANGS to convert them fd_displ to BOHR
          ! 13.6057039763 / BOHR_RADIUS_ANGS convert to eV / Ang
          ctr = ctr + 2
          ctr2 = ctr2 + 1
        ENDDO ! j
      ENDDO ! p
    ENDDO ! k
    !
    ALLOCATE(force_consts_av(nat * nq_tot, 3, nat * nq_tot * 3))
    force_consts_av = 0.d0
    ! perform mixing
    !
    IF (mixing) THEN
      ALLOCATE(force_consts_iters(nat * nq_tot, 3, nat * nq_tot * 3, iter_idx - iter_idx0 + 1)) 
      !
      force_consts_iters = 0.d0
      ctr = 1
      !
      DO i = iter_idx0, iter_idx - 1
        WRITE(pt_3,'(1I2.2)') i
        IF (i == 0) THEN 
          filename = 'FORCE_CONSTANTS_poly_iter_' // TRIM(pt_3)
          INQUIRE(FILE = filename, EXIST = file_exists)
          IF (.NOT. file_exists) CALL errore('ASDM', 'Missing iteration file with IFCs', i )
          OPEN (unit = 85, file = filename, status = 'old', form = 'formatted')
          j = 1
          READ(85,'(2I4)') 
          DO WHILE (j <= nat * nq_tot * 3) 
            DO p = 1, nat * nq_tot ! atom index
              !
              READ(85,'(2I4)') 
              READ(85,'(3F14.8)')  force_consts_iters(p, :, j, ctr)
              READ(85,'(3F14.8)')  force_consts_iters(p, :, j + 1, ctr)
              READ(85,'(3F14.8)')  force_consts_iters(p, :, j + 2, ctr)
              !
            ENDDO ! p
            j = j + 3
          ENDDO ! j
          CLOSE(85)
        ELSE 
          filename = 'FORCE_CONSTANTS_' // TRIM( pt_T ) // 'K_iter_' // TRIM(pt_3)
          INQUIRE(FILE = filename, EXIST = file_exists)
          IF (.NOT. file_exists) CALL errore('ASDM', 'Missing iteration file with IFCs', i )
          OPEN (unit = 85, file = filename, status = 'old', form = 'formatted')
          j = 1
          READ(85,'(2I4)') 
          DO WHILE (j <= nat * nq_tot * 3) 
            DO p = 1, nat * nq_tot ! atom index
              !
              READ(85,'(2I4)') 
              READ(85,'(3F14.8)')  force_consts_iters(p, :, j, ctr)
              READ(85,'(3F14.8)')  force_consts_iters(p, :, j + 1, ctr)
              READ(85,'(3F14.8)')  force_consts_iters(p, :, j + 2, ctr)
              !
            ENDDO ! p
            j = j + 3
          ENDDO ! j
          CLOSE(85)
        ENDIF ! i
        ! take average
        force_consts_av = force_consts_av + force_consts_iters(:, :, :, ctr) / (iter_idx - iter_idx0 + 1) 
        ctr = ctr + 1
      ENDDO ! i
      ! add also current IFCs
      force_consts_av = force_consts_av + force_consts / (iter_idx - iter_idx0 + 1) 
    ELSE 
    !
      force_consts_av = force_consts
    !
    ENDIF ! mixing
    ! 
    !
    ! xkg: to convert forces into a format for performing mat multipl.
    ALLOCATE(xkg(3, 3,  nat * nq_tot,  nat * nq_tot))
    xkg = 0.d0
    !
    !WRITE(pt_1,'(1I2.2)') iter_idx 
    ! to print them in a format similar to phonopy
    filename = 'FORCE_CONSTANTS_' // TRIM( pt_T ) // 'K_iter_' // TRIM(pt_2)
    OPEN (unit = 85, file = filename, status = 'unknown', form = 'formatted')
    j = 1
    ctr2 = 1
    WRITE(85,'(2I4)') nat * nq_tot, nat * nq_tot
    DO WHILE (j <= nat * nq_tot * 3) 
      DO p = 1, nat * nq_tot ! atom index
        !
        WRITE(85,'(2I4)') ctr2, p
        WRITE(85,'(3F14.8)')  force_consts(p, :, j)
        WRITE(85,'(3F14.8)')  force_consts(p, :, j + 1)
        WRITE(85,'(3F14.8)')  force_consts(p, :, j + 2)
        ! We pass to xkg the IFCs after mixing
        xkg(1, :, ctr2, p) =  force_consts_av(p, :, j)
        xkg(2, :, ctr2, p) =  force_consts_av(p, :, j + 1)
        xkg(3, :, ctr2, p) =  force_consts_av(p, :, j + 2)
        !
      ENDDO ! p
      j = j + 3
      ctr2 = ctr2 + 1
    ENDDO ! j
    CLOSE(85)
    ! 
    ! Find symmetries using equilibirum structure
    ! convert first equil_pos in required format
    ALLOCATE(equil_p2(3, nat * nq_tot), ityp_super(nat * nq_tot))
    !
    at_or = at ! need to save original at to print the unit-cell fc info below
    !
    atx(:, 1) = dim1 * at(:, 1)
    atx(:, 2) = dim2 * at(:, 2)
    atx(:, 3) = dim3 * at(:, 3)
    ! alat = SQRT ( at(1,1)**2+at(2,1)**2+at(3,1)**2 ) as in cell_base.f90
    alatx = SQRT ( atx(1, 1)**2 + atx(2, 1)**2 + atx(3, 1)**2 )
    atx(:, :) = atx(:, :) / alatx 
    !
    at = atx ! at accounts for supercell now
    !
    CALL recips( at(1, 1), at(1, 2), at(1, 3), bg(1, 1), bg(1, 2), bg(1, 3) )
    !
    ! generate Cart. coordinates in alat
    ctr = 1
    DO kk = 1, nat
      DO qp = 1, nq_tot
         !equil_p2(:, ctr) = crystal_pos(qp, kk, :) 
         equil_p2(:, ctr) = crystal_pos(qp, kk, 1) * at(:, 1)  & 
                          + crystal_pos(qp, kk, 2) * at(:, 2)  & 
                          + crystal_pos(qp, kk, 3) * at(:, 3) 
         !WRITE(*,'(A20, 3F14.8)' ) "Equil alat", equil_p2(:, ctr)
         ityp_super(ctr) = ityp(kk)
         ctr = ctr + 1
      ENDDO
    ENDDO
    !
    ! irt is the index mat that maps the atoms after symmetry oper. 
    IF ( ALLOCATED( irt ) ) DEALLOCATE( irt )
    ALLOCATE( irt(48, nat * nq_tot)) 
    !
    CALL set_sym_bl ( ) ! generate sym. environment
    !
    magnetic_sym = noncolin .AND. domag
    CALL find_sym(nat * nq_tot, equil_p2, ityp_super, magnetic_sym, m_loc)
    !
    WRITE(*, *) "  ======================================"
    CALL print_symmetries ( 1, noncolin, domag )
    WRITE(*, *) "  ======================================"
    WRITE(*, *) 
    !
    ! Now apply symmetries on force constants.
    ! WRITE(*,*) "Num of sym", nsym, nrot
    !
    ALLOCATE(xkg_av(3, 3, nat * nq_tot, nat * nq_tot), xkg_tmp(3, 3))
    xkg_av = 0.d0
    xkg_tmp = 0.d0
    DO k = 1, nsym
      DO i = 1, nat * nq_tot
        DO j = 1, nat * nq_tot
          xkg_tmp(:, :) =  MATMUL(TRANSPOSE(sr(:, :, k)), MATMUL(xkg(:, :, irt(k, i), irt(k, j)), sr(:, :, k)))
          xkg_av(:, :, i, j) = xkg_av(:, :, i, j) + xkg_tmp(:, :) / DBLE(nsym)
        ENDDO
      ENDDO
    ENDDO
    !
    ! now check translations
    !
    ! convert to crystal coordinates
    ctr = 1
    DO kk = 1, nat
      DO qp = 1, nq_tot
         equil_p2(:, ctr) = crystal_pos(qp, kk, :) 
         IF (equil_p2(1, ctr) < 0.d0) equil_p2(1, ctr) = equil_p2(1, ctr) + 1 ! only positive coord.
         IF (equil_p2(2, ctr) < 0.d0) equil_p2(2, ctr) = equil_p2(2, ctr) + 1 ! only positive coord.
         IF (equil_p2(3, ctr) < 0.d0) equil_p2(3, ctr) = equil_p2(3, ctr) + 1 ! only positive coord.
         ctr = ctr + 1
      ENDDO
    ENDDO
    !
    ! Apply transl and fill irt_trans
    ALLOCATE(eqp_tmp(3, nat * nq_tot), irt_trans(nat * nq_tot, nq_tot))
    !
    ctr = 1
    DO i = 0, dim1 - 1
      DO j = 0, dim2 - 1
        DO k = 0, dim3 - 1
          DO kk = 1, nat * nq_tot
              eqp_tmp(1, kk) = equil_p2(1, kk) - i / DBLE(dim1)
              IF (eqp_tmp(1, kk) < 0.d0) eqp_tmp(1, kk) = eqp_tmp(1, kk) + 1
              !
              eqp_tmp(2, kk) = equil_p2(2, kk) - j / DBLE(dim2)
              IF (eqp_tmp(2, kk) < 0.d0) eqp_tmp(2, kk) = eqp_tmp(2, kk) + 1
              !
              eqp_tmp(3, kk) = equil_p2(3, kk) - k / DBLE(dim3)
              IF (eqp_tmp(3, kk) < 0.d0) eqp_tmp(3, kk) = eqp_tmp(3, kk) + 1
          ENDDO ! kk
          ! 
          ! check distances if equal
          DO ii = 1, nat * nq_tot
            DO kk = 1, nat * nq_tot
              IF (NORM2(eqp_tmp(:, ii) - equil_p2(:, kk)) < eps ) irt_trans(ii, ctr) = kk
            ENDDO ! kk
          ENDDO ! ii
          ctr = ctr + 1
        ENDDO ! k
      ENDDO ! j
    ENDDO ! i
    ! 
    WRITE(*,*)  
    WRITE(*,*) "  Atom indices connected by translations"
    WRITE(*,*) "  ======================================"
    DO kk = 1, nat * nq_tot
      WRITE(*,'(1000I4)') irt_trans(kk, :)
    ENDDO
    !
    ! Now take average over force constants related by translation
    ALLOCATE(xkg_new(3, 3, nat * nq_tot, nat * nq_tot, nq_tot), xkg_new_av(3, 3, nat * nq_tot, nat * nq_tot))
    !
    xkg_new = 0.d0
    DO k = 1, nq_tot
      DO i = 1, nat * nq_tot
        DO j = 1, nat * nq_tot
          xkg_new(:, :, irt_trans(i, k), irt_trans(j, k), k) = xkg_av(:, :, i, j)
         ENDDO
      ENDDO
    ENDDO
    !
    xkg_new_av = 0.d0
    DO i = 1, nq_tot
      xkg_new_av(:, :, :, :) = xkg_new_av(:, :, :, :) + DBLE(xkg_new(:, :, :, :, i) / nq_tot)
    ENDDO
    ! 
    !!! Impose hermiticity
    !!DO i = 1, nat * nq_tot
    !!  DO j = 1, nat * nq_tot
    !!    DO k = 1, 3
    !!      !DO kk = 1, 3
    !!      kk = 1
    !!      DO WHILE (kk < k) 
    !!        xkg_new_av(k, kk, i, j) = (xkg_new_av(k, kk, i, j) + xkg_new_av(kk, k, i, j))/2
    !!        xkg_new_av(kk, k, i, j) = xkg_new_av(k, kk, i, j) 
    !!        kk = kk + 1
    !!      ENDDO
    !!    ENDDO
    !!  ENDDO
    !!ENDDO
    !
    ! Now determine new thermal equilibrium positions using Newton-Raphson
    ! minimization approach
    ALLOCATE(new_tau(3, nat))
    !
    IF ( update_equil ) THEN ! update_equil
    !
      ALLOCATE(eqf(nat * nq_tot, 3))
      filename = 'fd_forces' // '/' // 'ZG-scf_' // TRIM( pt_T ) // 'K_iter_' // TRIM(pt_2) // '.out'
      ! ZG-scf_222_1188.00K.out
      INQUIRE(FILE = filename, EXIST = file_exists)
      IF (.NOT. file_exists) CALL errore('ZG', 'Missing ZG-scf file in fd_forces/ ...', 1 )
      OPEN (unit = iu, file = filename, status = 'old', form = 'formatted')
      DO
        READ (iu, *, iostat = iosf) text ! iosf to check if we reach the
                                        ! end of file without finding forces
        READ (text, *) word 
        IF (word == search_str) THEN
          !WRITE(*,*) search_str, " found in ", filename
          READ (iu, *) ! skip empty line
          DO ii = 1, nat * nq_tot
            READ(iu,'(A34, 3F14.8)') A1, eqf(ii, :)
            ! Format should be the same with forces.f90 in PW/src
            !WRITE(*, '(A34, 3F14.8)') A1, eqf(ii, :)
          ENDDO ! ii
          EXIT ! exit do loop                
        ELSE 
          IF (iosf /= 0) CALL errore('ASDM',& 
                                'forces in ZG-scf file are missing ...', 1)
        ENDIF ! word
        ! 
      ENDDO ! do loop
      CLOSE(unit = iu)
      !
      ALLOCATE(eqf_av_iters(nat * nq_tot, 3))
      eqf_av_iters = 0.d0
      ! perform mixing over forces
      IF (mixing) THEN
        !
        ALLOCATE(eqf_iters(nat * nq_tot, 3, iter_idx - iter_idx0 + 1)) 
        !
        eqf_iters = 0.d0
        ctr = 1
        !
        DO i = iter_idx0, iter_idx - 1
          WRITE(pt_3,'(1I2.2)') i
          IF (i == 0) THEN 
            filename = 'fd_forces' // '/' // 'ZG-scf_poly_iter_' // TRIM(pt_3) // '.out'
            INQUIRE(FILE = filename, EXIST = file_exists)
            IF (.NOT. file_exists) CALL errore('ASDM', 'Missing iteration file with forces', i )
            OPEN (unit = iu, file = filename, status = 'old', form = 'formatted')
            DO
              READ (iu, *, iostat = iosf) text ! iosf to check if we reach the
                                              ! end of file without finding forces
              READ (text, *) word 
              IF (word == search_str) THEN
                !WRITE(*,*) search_str, " found in ", filename
                READ (iu, *) ! skip empty line
                DO ii = 1, nat * nq_tot
                  READ(iu,'(A34, 3F14.8)') A1, eqf_iters(ii, :, ctr)
                  ! Format should be the same with forces.f90 in PW/src
                  !WRITE(*, '(A34, 3F14.8)') A1, eqf(ii, :)
                ENDDO ! ii
                EXIT ! exit do loop                
              ELSE 
                IF (iosf /= 0) CALL errore('ASDM',& 
                                      'forces in ZG-scf file are missing ...', 1)
              ENDIF ! word
              ! 
            ENDDO ! do loop
            CLOSE(unit = iu)
            !
          ELSE 
            filename = 'fd_forces' // '/' // 'ZG-scf_' // TRIM( pt_T ) // 'K_iter_' // TRIM(pt_3) // '.out'
            INQUIRE(FILE = filename, EXIST = file_exists)
            IF (.NOT. file_exists) CALL errore('ASDM', 'Missing iteration file with forces', i )
            OPEN (unit = iu, file = filename, status = 'old', form = 'formatted')
            DO
              READ (iu, *, iostat = iosf) text ! iosf to check if we reach the
                                              ! end of file without finding forces
              READ (text, *) word 
              IF (word == search_str) THEN
                !WRITE(*,*) search_str, " found in ", filename
                READ (iu, *) ! skip empty line
                DO ii = 1, nat * nq_tot
                  READ(iu,'(A34, 3F14.8)') A1, eqf_iters(ii, :, ctr)
                  ! Format should be the same with forces.f90 in PW/src
                  !WRITE(*, '(A34, 3F14.8)') A1, eqf(ii, :)
                ENDDO ! ii
                EXIT ! exit do loop                
              ELSE 
                IF (iosf /= 0) CALL errore('ASDM',& 
                                      'forces in ZG-scf file are missing ...', 1)
              ENDIF ! word
              ! 
            ENDDO ! do loop
            CLOSE(unit = iu)
            !
          ENDIF ! i
          ! take average
          eqf_av_iters = eqf_av_iters + eqf_iters(:, :, ctr) / (iter_idx - iter_idx0 + 1) 
          ctr = ctr + 1
        ENDDO ! i
        ! add also current IFCs
        eqf_av_iters = eqf_av_iters + eqf / (iter_idx - iter_idx0 + 1) 
        !
        DEALLOCATE(eqf_iters)
        !
      ELSE 
      !
        eqf_av_iters = eqf
      !
      ENDIF ! mixing
      ! Impose symmetries on eqf
      ALLOCATE(eqf_av(nat * nq_tot, 3), eqf_tmp(3))
      eqf_av = 0.d0
      eqf_tmp = 0.d0
      DO k = 1, nsym
        DO i = 1, nat * nq_tot
            eqf_tmp(:) =  MATMUL(TRANSPOSE(sr(:, :, k)), eqf_av_iters(irt(k, i), :))
            eqf_av(i, :) = eqf_av(i, :) + eqf_tmp(:) / DBLE(nsym)
        ENDDO
      ENDDO
      !
      ! Now impose translational invariance on eqf_av
      ! 
      ALLOCATE(eqf_new( nat * nq_tot, 3, nq_tot), eqf_new_av(nat * nq_tot, 3))
      !
      eqf_new = 0.d0
      DO k = 1, nq_tot
        DO i = 1, nat * nq_tot
            eqf_new(irt_trans(i, k), :, k) = eqf_av(i, :)
        ENDDO
      ENDDO
      !
      eqf_new_av = 0.d0
      DO i = 1, nq_tot
        eqf_new_av(:, :) = eqf_new_av(:, :) + DBLE(eqf_new(:, :, i) / nq_tot)
      ENDDO
      !
      !!DO i = 1, nat * nq_tot
      !!  WRITE(*, *) "ooo", eqf_new_av(i, :), eqf(i, :)
      !!ENDDO
      !
      ! Now evaluate new equilibrium positions 
      ! using Newton-Raphson minimization
      !
      new_tau = 0.0d0
      ctr = 1
      !
      DO kk = 1, nat 
        inv_xkg = xkg_new_av(:, :, ctr, ctr)
        !
        CALL MatInv( 'G', 3, inv_xkg )
        ! MatInv from utils.f90 returns the inverse of a real matrix
        new_tau(1, kk) = tau(1, kk) + inv_xkg(1, 1) * eqf_new_av(ctr, 1) / alat / dim1 
        new_tau(2, kk) = tau(2, kk) + inv_xkg(2, 2) * eqf_new_av(ctr, 2) / alat / dim2
        new_tau(3, kk) = tau(3, kk) + inv_xkg(3, 3) * eqf_new_av(ctr, 3) / alat / dim3        
        ! We convert everything in unitcell alat
        ! Note that we have plus instead of minus since F = -dU/Dτ
        DO qp = 1, nq_tot
           ctr = ctr + 1 ! to skip entries in repeated unit cells
        ENDDO
      ENDDO
      ! Print new thermal equilibrium positions
      !DO kk = 1, nat
      !  WRITE(*,'(A20, 3F14.8)' ) "Equil alat", new_tau(1, kk) , & 
      !                                          new_tau(2, kk) , &
      !                                          new_tau(3, kk)  
     
      !ENDDO
      !
      DEALLOCATE(eqf, eqf_av, eqf_tmp, eqf_new, eqf_new_av, eqf_av_iters)
    ELSE 
    !
      new_tau = tau
    !
    ENDIF ! update_equil 
    !
    ! Print them in phonopy format
    !WRITE(pt_1,'(1I2.2)') iter_idx 
    filename = 'FORCE_CONSTANTS_sym_' // TRIM( pt_T ) // 'K_iter_' // TRIM(pt_2)
    OPEN (unit = 85, file = filename, status = 'unknown', form = 'formatted')
    WRITE(85, '(2i4)')  nat * nq_tot, nat * nq_tot 
    DO i = 1, nat * nq_tot
      DO j = 1, nat * nq_tot
      !
      WRITE(85, '(2i4)') i, j
      WRITE(85, '(3F14.8)') xkg_new_av(1, :, i, j)
      WRITE(85, '(3F14.8)') xkg_new_av(2, :, i, j)
      WRITE(85, '(3F14.8)') xkg_new_av(3, :, i, j)
      !
      ENDDO
    ENDDO
    CLOSE(85)
    ! To print in q2r.fc format
    !
    ctr2 = 1
    DO i = 1, dim3 
      DO j = 1, dim2 
        DO  k = 1, dim1
          Rlist(ctr2, 1) = k
          Rlist(ctr2, 2) = j
          Rlist(ctr2, 3) = i
          ctr2 = ctr2 + 1
        ENDDO
      ENDDO
    ENDDO
    !
    filename = TRIM( pt_T ) // 'K_iter_' // TRIM(pt_2) // '.fc'
    OPEN (unit = 85, file = filename, status = 'unknown', form = 'formatted')
    !copy from do_q2r routine
    !Write initial part of fc file (info of the cell, dielectric const, BEC
    WRITE(85,'(i3,i5,i4,6f11.7)') ntyp, nat, ibrav, celldm
    IF (ibrav==0) then
      WRITE (85,'(2x,3f15.9)') ((at_or(i, j), i = 1, 3),j = 1, 3)
    ENDIF
    !
    DO kk = 1, ntyp
      WRITE(85, *) kk, " '", atm(kk), "' ", amass(kk) * amu_ry
    ENDDO
    !
    DO kk = 1, nat
      WRITE(85, '(2i5,3f18.10)') kk, ityp(kk), (new_tau(j, kk),j = 1, 3)
    ENDDO
    !
    WRITE (85, *) incl_epsil
    IF (incl_epsil) THEN
      WRITE(85, '(3f24.12)') ((epsil(i, j), j = 1, 3), i = 1, 3)
      DO kk = 1, nat
        WRITE(85, '(i5)') kk
        WRITE(85, '(3f15.7)') ((zeu(i, j, kk), j = 1, 3), i = 1, 3)
      ENDDO
    ENDIF ! incl_epsil
    WRITE(85, '(3I4)') dim1, dim2, dim3
    !
    ! WRITE FORCE CONSTANTS
    DO i = 1, 3
      DO j = 1, 3
        DO ii = 1, nat
          ctr = 1
          DO kk = 1, nat 
            WRITE(85, '(4I4)') i, j, ii, kk
            DO k = 1, nq_tot
              WRITE(85, '(3I3, 1F14.8)') Rlist(k, :), xkg_new_av(i, j, (ii - 1) * nq_tot + 1, ctr) 
              ctr = ctr + 1
            ENDDO ! k
          ENDDO ! kk
        ENDDO ! ii        
      ENDDO ! j
    ENDDO ! i
    CLOSE(85)
  !
    DEALLOCATE(fd_forces, force_consts, equil_p2, ityp_super)
    DEALLOCATE(xkg, xkg_av, xkg_tmp, eqp_tmp, irt_trans)
    DEALLOCATE(xkg_new, xkg_new_av, force_consts_av, new_tau)    
    IF (mixing) DEALLOCATE (force_consts_iters)
  !
  ENDIF ! read_fd_forces
  ! 
  ! 
  RETURN
END SUBROUTINE 
!
