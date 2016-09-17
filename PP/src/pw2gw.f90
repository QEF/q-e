
! Copyright (C) 2005-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 13Aprile2005
! GENERATES INPUT for GW code
!tested on: Silicon bulk, Germanium Bulk, Na4, InP bulk
! Please note just symmorphic symm. op. have to be used
! Use input option of pw.x: force_symmorphic=.TRUE.

!-----------------------------------------------------------------------
PROGRAM pw2gw
  !-----------------------------------------------------------------------

  ! This subroutine writes files containing plane wave coefficients
  ! and other stuff needed by GW codes

  USE io_files,   ONLY : prefix, tmp_dir
  USE io_global,  ONLY : ionode, ionode_id
  USE mp,         ONLY : mp_bcast
  USE mp_world,   ONLY : world_comm, nproc
  USE mp_global,  ONLY : mp_startup
  USE mp_pools,   ONLY : kunit
  USE environment,ONLY : environment_start, environment_end
  USE us,         ONLY : spline_ps
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  CHARACTER(LEN=256) :: outdir
  !
  INTEGER :: ios
  INTEGER :: kunittmp
  LOGICAL :: use_gmaps
  CHARACTER(len=20) :: what
  CHARACTER(len=30) :: when

  NAMELIST / inputpp / prefix, outdir, what, use_gmaps
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'PW2GW' )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  what   = 'gw'
  use_gmaps = .false.

  ios = 0
  IF ( ionode )  THEN
     !
     READ (5, inputpp, iostat=ios)
     tmp_dir = trimcheck (outdir)
     !
  ENDIF
  !
  CALL mp_bcast( ios, ionode_id, world_comm )
  IF (ios /= 0)   CALL errore('pw2gw', 'reading inputpp namelist', abs(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast(tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( what, ionode_id, world_comm )
  CALL mp_bcast( use_gmaps, ionode_id, world_comm )
  !

  spline_ps = .false.

  CALL read_file
  CALL openfil_pp
  !
  CALL mp_bcast(spline_ps, ionode_id, world_comm)
#if defined __MPI
  kunittmp = kunit
#else
  kunittmp = 1
#endif
  !
  IF( trim( what ) == 'gw' ) THEN
    CALL compute_gw  ( use_gmaps )
  ELSE
    CALL write_gmaps ( kunittmp )
  ENDIF
  !
  CALL environment_end ( 'PW2GW' )
  !
  CALL stop_pp

END PROGRAM pw2gw


SUBROUTINE compute_gw( use_gmaps )

  ! This routine creates the QPLDA and the matrixelements
  ! tform = .false. UNFORMATTED QPLDA
  ! tform = .true.  FORMATTED QPLDA
  ! tsingle must be always true

  USE kinds,     ONLY : DP, sgl
  USE constants, ONLY : eps8, pi, AUTOEV, rytoev
  USE cell_base, ONLY : alat, tpiba2, at, bg, omega
  USE symm_base, ONLY : s, nsym
  USE wvfct,     ONLY : npwx, nbnd, wg, et
  USE gvecw,     ONLY : gcutw
  USE control_flags, ONLY : gamma_only
  USE gvect,         ONLY : ngm, g, gg, ig_l2g, nl
  USE fft_base,  ONLY: dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE klist ,        ONLY : nks, xk, wk, ngk, igk_k
  USE lsda_mod,      ONLY : nspin
  USE io_files,      ONLY : nwordwfc, iunwfc
  USE wavefunctions_module, ONLY : evc, psic
  USE mp_global, ONLY : intra_image_comm, npool
  USE io_global, ONLY : ionode, ionode_id
  USE mp,        ONLY : mp_sum , mp_max
  USE mp_world,  ONLY : world_comm, mpime, nproc
  USE mp_wave,   ONLY : mergewf
  USE parallel_include
  USE scf,       ONLY : rho, rho_core, rhog_core
  USE ener,      ONLY : etxc, vtxc

  USE uspp_param, ONLY : upf, nh
  USE uspp,       ONLY : nhtol
  USE us,         ONLY : tab, tab_d2y, spline_ps
  USE ions_base,  ONLY : ntyp => nsp
  USE klist,      ONLY : ngk

  IMPLICIT NONE

  LOGICAL, INTENT(in) :: use_gmaps

  INTEGER :: ii(16), ngw, nkpt, ig, ik, ir, n, i,j,k, io = 98, iband1, iband2
  INTEGER :: npw, omax, o, iproc
  INTEGER, ALLOCATABLE :: in1(:), in2(:), in3(:)
  INTEGER, ALLOCATABLE :: in1_tmp(:), in2_tmp(:), in3_tmp(:)
  INTEGER, ALLOCATABLE :: inx_rcv(:), ig_l2g_rcv(:)
  LOGICAL :: t_form = .false., t_single = .true.
  REAL(kind=sgl) :: a1_s(3), a2_s(3), a3_s(3)
  REAL(kind=sgl), ALLOCATABLE :: xk_s(:,:), eig_s(:,:), focc_s(:,:)
  REAL(kind=DP):: g2max, a1(3), a2(3), a3(3),norm, xkgk(3), rrhotwx(3), delta
  REAL(kind=DP):: alpha, egap, halfalpha, Df, const, dummy
  REAL(kind=DP), PARAMETER :: omegamax = 30.0
  REAL(kind=DP), ALLOCATABLE:: gsort(:), eig(:,:), focc(:,:), kpg(:,:), omegatt(:), omeg(:)
  REAL(kind=DP), ALLOCATABLE:: pp1(:,:), pp2(:,:), pp3(:,:)
  REAL(kind=DP), ALLOCATABLE:: epsx(:,:), epsy(:,:), epsz(:,:)
  REAL(kind=DP), ALLOCATABLE:: epstx(:), epsty(:), epstz(:)
  REAL(kind=DP) :: epsxx, epsyy, epszz
  REAL(kind=DP) :: vxcdiag
  REAL(kind=DP), ALLOCATABLE :: vxc(:,:)
  COMPLEX(kind=DP):: rhotwx(3), ctemp, dasomma(3)
  COMPLEX(kind=DP),  ALLOCATABLE:: c0(:), c0_m(:,:), c0_tmp_dp(:) !, c0_tmp(:) !, c0_gamma(:)
  COMPLEX(kind=sgl), ALLOCATABLE:: c0_s(:), c0_tmp(:) !, c0_gamma_s(:)
  CHARACTER(len=80) :: titleo(2)
  INTEGER :: igwx, igwxx, comm, ierr, ig_max, igwx_r
  INTEGER :: igwx_p(nproc)
  INTEGER, ALLOCATABLE :: igk_l2g(:)
  !
  REAL(kind=DP), ALLOCATABLE :: vkb0(:), djl(:), vec_tab(:), vec_tab_d2y(:)
  INTEGER :: nb, nt, size_tab, size_tab_d2y, ipw, l
  !
  ! REAL(kind=DP) :: norma ! Variable needed only for DEBUG
  !
#if defined __MPI
  INTEGER :: istatus( MPI_STATUS_SIZE )
#endif
  !
  IF( nspin > 1 ) CALL errore('pw2gw','Spin polarization not implemented',1)
  IF( npool > 1 ) CALL errore('pw2gw','parallel run with pools not allowed yet',1)
  !
  !
  IF( mpime == 0 ) THEN
     IF (t_form) THEN
        WRITE (6,'(//" writing LDA info on unit 98 FORMATTED")')
        OPEN (io, FILE='QPLDA',STATUS='unknown',FORM='FORMATTED')
     ELSE
        WRITE (6,'(//" writing LDA info on unit io UNFORMATTED")')
        OPEN (io, FILE='QPLDA',STATUS='unknown',FORM='UNFORMATTED')
     ENDIF
     WRITE (6,'(//" writing matrixelements on unit 98 FORMATTED")')
     OPEN (90, FILE='matrixelements',STATUS='unknown',FORM='FORMATTED')
  ENDIF
  !
  !  file's title [2 lines]
  !
  titleo(1)='pw2gw'
  titleo(2)='test version'
  IF( mpime == 0 ) THEN
     IF (t_form) THEN
        WRITE (io,'(A80/A80)') titleo(1), titleo(2)
     ELSE
        WRITE (io) titleo(1)
        WRITE (io) titleo(2)
     ENDIF
     !
     WRITE(6,*) 'qplda title'
     WRITE(6,*) titleo(1)
     WRITE(6,*) titleo(2)
  ENDIF
  !
  !  Read 16 integers (reserved for future flags)
  !  Flags used so far:
  !   I1 = 0 if QPLDA file is formatted, 1 if unformatted
  !   I2 = 0 if RWG format, 1 if BF format
  !   I3 = 1 if non-symmorphic operations (+vectors) included, otherwise 0
  !
  ii(:) = 0
  IF (t_form) THEN
     ii(1)=0
     IF( mpime == 0 ) WRITE (io,'(16I5)') ii
  ELSE
     ii(1)=1
     IF( mpime == 0 ) WRITE (io) ii
  ENDIF
  !
  WRITE(6,'(16I5)') ii
  !
  !  write real-space lattice vectors (Cartesian, in au) [3 lines]
  !
  a1(:)=at(:,1)*alat
  a2(:)=at(:,2)*alat
  a3(:)=at(:,3)*alat
  a1_s(:) = a1(:)
  a2_s(:) = a2(:)
  a3_s(:) = a3(:)

  IF( mpime == 0 ) THEN
     !
     IF (t_form) THEN
        WRITE (io,'(3E26.18)') a1, a2, a3
     ELSE
        IF (t_single) THEN
           WRITE (io) a1_s, a2_s, a3_s
        ELSE
           WRITE (io) a1, a2, a3
        ENDIF
     ENDIF
     !
     WRITE(6,*) 'Vettori di reticolo diretto'
     WRITE(6,'(a,3E26.18)') 'a1', a1_s
     WRITE(6,'(a,3E26.18)') 'a2', a2_s
     WRITE(6,'(a,3E26.18)') 'a3', a3_s
     !
  ENDIF
  !
  ! Write symmetry operations.
  ! The matrix s is the transpose of the symmetry matrix in direct space,
  ! in units of a_i. But the transpose of the symmetry matrix in real space
  ! is the symmetry matrix in reciprocal space so "s" is already the symmetry
  ! matrix in reciprocal space in units of b_i
  ! The gw code will read row by row a matrix and will treat it as symmetry
  ! matrix in reciprocal space in units of b_i
  ! In other words, the gw code expects as input the direct space symmetry
  ! matrix, in units of a_i, written columnwise
  !
  IF( mpime == 0 ) THEN
     WRITE(6,*)'nrot=',nsym
     WRITE(6,'(3E26.18)') (((float(s(i,j,k)),j=1,3),i=1,3),k=1,nsym)
     IF (t_form) THEN
        WRITE (io,'(I2)') nsym
        WRITE (io,'(3E26.18)') (((float(s(i,j,k)),j=1,3),i=1,3),k=1,nsym)
        IF (ii(3) == 1) THEN
           ! READ (10,1020) ((VOFFSET(I,J),I=1,3),J=1,NOP)
           ! WRITE (6,'(//" Run program CNVNSY to convert QPLDA file first.")')
           CALL errore('pw2gw','non-symmorphic translation vectors',ii(3))
        ENDIF
     ELSE
        WRITE (io) nsym
        IF (t_single) THEN
           WRITE (io) (((float(s(i,j,k)),j=1,3),i=1,3),k=1,nsym)
        ELSE
           WRITE (io) (((dble(s(i,j,k)),j=1,3),i=1,3),k=1,nsym)
        ENDIF
        IF (ii(3) == 1) THEN
           ! READ (10,1020) ((VOFFSET(I,J),I=1,3),J=1,NOP)
           CALL errore('pw2gw','non-symmorphic translation vectors',ii(3))
        ENDIF
     ENDIF
  ENDIF
  !
  !  write reciprocal lattice vectors (in reciprocal lattice units;
  !  ie in the basis of the reciprocal lattice basis vectors)
  !
  !  PWscf stores psi(k+G), using |k+G| to order the components;
  !  GW codes require on input psi_k(G), using the same set of G
  !
  g2max = 0.0d0
  !DEBUG
  IF (ionode) WRITE(6,*) ' nks ', nks
  IF (ionode) WRITE(6,*) ' k points in  cartesian coordinates'
  IF (ionode) WRITE(6,'(1x,3f10.6)') ( (xk(i,ik),i=1,3), ik=1,nks)
  !DEBUG
  igwx  = 0  !  maximum G vector index
  g2max = gcutw ! RAGGIO DELLA SFERA |G+k|<cut
  DO ik = 1, nks
     npw = ngk(ik)
     ! WRITE( 6, * ) 'DEBUG g2max ', g2max
     ! g2max, g2kin = RAGGIO DELLA SFERA |G+k|<cut, non MASSIMO |G| nella sfera
     ! g2max <= gcutw   PER COSTRUZIONE
     igwx = max( igwx, maxval( igk_k(1:npw,ik) ) )
  ENDDO
  !IF (ionode) write(*,*) "igwx = ", igwx
  !
  !  ngw =  number of G-vectors (complete shells) such that G2 <= G2max
  !  ovvero <= RAGGIO della SFERA, in pratica trova i G2 relativi a GAMMA
  !
  DO ngw = 1, ngm
     IF ( gg(ngw) > g2max + eps8) GOTO 100
  ENDDO
  CALL errore ( 'pw2gw','max G in PW not found?!?',ngw)
100 ngw = ngw - 1

  ! Pongo NGW pari al massimo indice tra i vettori G che fanno parte delle
  ! sfere |G+k|<cut per qualsiasi k
  !
  !IF (ionode) write( 6, * ) ' igwx= ', igwx

  ngw = igwx

  !
  !  PWscf stores G in order of increasing module
  !
  ALLOCATE (in1(ngw), in2(ngw), in3(ngw))
  DO ig=1,ngw
     in1(ig) = nint ( at(1,1)*g(1,ig) + at(2,1)*g(2,ig) + at(3,1)*g(3,ig) )
     in2(ig) = nint ( at(1,2)*g(1,ig) + at(2,2)*g(2,ig) + at(3,2)*g(3,ig) )
     in3(ig) = nint ( at(1,3)*g(1,ig) + at(2,3)*g(2,ig) + at(3,3)*g(3,ig) )
  ENDDO

  igwxx = maxval( ig_l2g( 1:ngw ) )
  CALL mp_max( igwxx, world_comm )
  IF (ionode) WRITE(*,*) "NDIMCP = ", igwxx

  igwx_p = 0
  igwx_p( mpime + 1 ) = igwx
  CALL mp_sum( igwx_p, world_comm )

  IF( mpime == 0 ) THEN
     !
     ! allocate arrays
     !
     ALLOCATE ( in1_tmp(igwxx), in2_tmp(igwxx), in3_tmp(igwxx) )
     ! copy local data of the root proc into global vector
     DO ig = 1, ngw
        in1_tmp( ig_l2g(ig) ) = in1(ig)
        in2_tmp( ig_l2g(ig) ) = in2(ig)
        in3_tmp( ig_l2g(ig) ) = in3(ig)
     ENDDO
     !
#if defined __MPI
     ALLOCATE( ig_l2g_rcv( igwxx ) )
     ALLOCATE( inx_rcv( igwxx ) )
     !
     DO iproc = 2, nproc
        CALL MPI_RECV( ig_l2g_rcv, igwx_p( iproc ), MPI_INTEGER, (iproc-1), iproc,       intra_image_comm, istatus, IERR )
        !
        CALL MPI_RECV( inx_rcv,    igwx_p( iproc ), MPI_INTEGER, (iproc-1), iproc+NPROC, intra_image_comm, istatus, IERR )
        DO ig = 1, igwx_p( iproc )
           in1_tmp( ig_l2g_rcv( ig ) ) = inx_rcv( ig )
        ENDDO
        CALL MPI_RECV( inx_rcv,    igwx_p( iproc ), MPI_INTEGER, (iproc-1), iproc+2*NPROC, intra_image_comm, istatus, IERR )
        DO ig = 1, igwx_p( iproc )
           in2_tmp( ig_l2g_rcv( ig ) ) = inx_rcv( ig )
        ENDDO
        CALL MPI_RECV( inx_rcv,    igwx_p( iproc ), MPI_INTEGER, (iproc-1), iproc+3*NPROC, intra_image_comm, istatus, IERR )
        DO ig = 1, igwx_p( iproc )
           in3_tmp( ig_l2g_rcv( ig ) ) = inx_rcv( ig )
        ENDDO
     ENDDO
     !
     DEALLOCATE( ig_l2g_rcv )
     DEALLOCATE( inx_rcv )
#endif
     !
  ELSE
     !
#if defined __MPI
     CALL MPI_SEND( ig_l2g, igwx, MPI_INTEGER, 0, mpime+1,         intra_image_comm, IERR )
     CALL MPI_SEND( in1(1), igwx, MPI_INTEGER, 0, mpime+1+NPROC,   intra_image_comm, IERR )
     CALL MPI_SEND( in2(1), igwx, MPI_INTEGER, 0, mpime+1+2*NPROC, intra_image_comm, IERR )
     CALL MPI_SEND( in3(1), igwx, MPI_INTEGER, 0, mpime+1+3*NPROC, intra_image_comm, IERR )
#endif
     !
  ENDIF

  IF (mpime == 0) WRITE(*,*) "fine debug sui punti g"



  IF (t_form) THEN
     IF( mpime == 0 ) THEN
        WRITE (io,'(I12)') igwxx
        WRITE (io,'(3I5)') (in1_tmp(ig),in2_tmp(ig),in3_tmp(ig),ig=1,igwxx)
     ENDIF
  ELSE
     IF( mpime == 0 ) THEN
        WRITE (io) igwxx
        WRITE (io) (in1_tmp(ig),in2_tmp(ig),in3_tmp(ig),ig=1,igwxx)
     ENDIF
  ENDIF
  !
  DEALLOCATE ( in1, in2, in3 )
  IF( mpime == 0 ) THEN
     DEALLOCATE ( in1_tmp, in2_tmp, in3_tmp )
  ENDIF
  !
  !  WRITE k-points (in RL units)
  !
  ! transformation in relative coordinates with respect to b1,b2,b3
  !
  CALL cryst_to_cart (nks, xk, at, -1)
  ! xk(3,nkpt) in input deve essere in coordinate cartesiane!
  nkpt = nks
  ALLOCATE (xk_s(3,nkpt))
  xk_s(:,:) = xk(:,1:nkpt)
  !
  IF( mpime == 0 ) THEN
     OPEN(65,file='k.dat')
     WRITE(65,'(1x,3f10.6,1x,f10.6)')  ( xk_s(:,ik), wk(ik)*0.5, ik=1,nks )
     CLOSE(unit=65)
  ENDIF

  IF( mpime == 0 ) THEN
     IF (t_form) THEN
        WRITE (io,'(I12)') nkpt
        WRITE (io,'(3E26.18)') ((xk_s(i,ik),i=1,3),ik=1,nkpt)
     ELSE
        WRITE (io) nkpt
        WRITE(6,*) 'nkpt',nkpt
        IF(t_single) THEN
           WRITE (io) ((xk_s(i,ik),i=1,3),ik=1,nkpt)
        ELSE
           WRITE (io) ((xk(i,ik),i=1,3),ik=1,nkpt)
        ENDIF
     ENDIF
     WRITE(6,'(1x,3f10.6)') ( (xk_s(i,ik),i=1,3), ik=1,nkpt)
  ENDIF
! --------------------------
! vkb0
! --------------------------
  DO ik=1,nkpt
    npw = ngk(ik)
    WRITE(15,*) "npw", npw
    ALLOCATE(vkb0(1:npw))

    size_tab=size(tab,1)
    size_tab_d2y=size(tab_d2y,1)

    ALLOCATE(vec_tab(1:size_tab))
    if(.not.allocated(vec_tab_d2y)) ALLOCATE(vec_tab_d2y(1:size_tab_d2y))

    DO nt = 1, ntyp
      DO nb = 1, upf(nt)%nbeta
        vkb0(:) = 0.0_dp
        vec_tab(:) = 0.0_dp
        vec_tab_d2y(:) = 0.0_dp
        vec_tab(:) = tab(:,nb,nt)
        IF(spline_ps) vec_tab_d2y(:) = tab_d2y(:,nb,nt)
        CALL gen_us_vkb0(ik,npw,vkb0,size_tab,vec_tab,spline_ps,vec_tab_d2y)
        WRITE(15,*) "---------------DEBUG-VKB0----------------------"
        WRITE(15,*) "ik= ", ik
        WRITE(15,*) "nt= ", nt
        WRITE(15,*) "nb= ", nb
        WRITE(15,*) "l= ", upf(nt)%lll(nb)
        WRITE (15,'(8f15.9)') vkb0
        WRITE(15,*) "--------------END-DEBUG------------------------"
!        WRITE(io) vkb0
      ENDDO
    ENDDO

   DEALLOCATE(vkb0)
   DEALLOCATE(vec_tab)
   IF(allocated(vec_tab_d2y))  DEALLOCATE(vec_tab_d2y)

  ENDDO
!---------------------------
! djl
!---------------------------
  DO ik=1,nkpt
    npw = ngk(ik)

    ALLOCATE(djl(1:npw))

    size_tab=size(tab,1)
    size_tab_d2y=size(tab_d2y,1)

    ALLOCATE(vec_tab(1:size_tab))
    IF(.not. allocated(vec_tab_d2y)) ALLOCATE(vec_tab_d2y(1:size_tab_d2y))
    DO nt = 1, ntyp
      DO nb = 1, upf(nt)%nbeta
        djl(:) = 0.0_dp
        vec_tab(:) = 0.0_dp
        vec_tab_d2y(:) = 0.0_dp
        vec_tab(:) = tab(:,nb,nt)
        IF(spline_ps) vec_tab_d2y(:) = tab_d2y(:,nb,nt)
        CALL gen_us_djl(ik,npw,djl,size_tab,vec_tab,spline_ps,vec_tab_d2y)
!        WRITE(0,*) "---------------DEBUG-----------------------"
!        WRITE(0,*) "spline: ", spline_ps
!        WRITE(0,*) "ik= ", ik
!        WRITE(0,*) "nt= ", nt
!        WRITE(0,*) "nb= ", nb
!        WRITE(0,*) "l= ", upf(nt)%lll(nb)
!        WRITE (0,'(8f15.9)') djl
!        WRITE(0,*) "--------------END-DEBUG------------------------"
!        WRITE(io) djl
      ENDDO
    ENDDO

   DEALLOCATE(djl)
   DEALLOCATE(vec_tab)
   IF(allocated(vec_tab_d2y))  DEALLOCATE(vec_tab_d2y)

  ENDDO
!-----------------------
!-----------------------

  !
  !  WRITE energies (Hartrees) (in ascending order, hopefully they are ordered)
  !
  n = nbnd
  ALLOCATE (eig(n,nkpt), eig_s(n,nkpt))
  eig(:,:)   = et(:,1:nkpt)*0.5d0

  IF (t_form) THEN
     IF( mpime == 0 ) WRITE (io,'(I12)') n
     IF( mpime == 0 ) WRITE (io,'(3E26.18)') ((eig(i,ik),ik=1,nkpt),i=1,n)
  ELSE
     IF( mpime == 0 ) WRITE (io) n
     IF(t_single) THEN
        DO ik=1,nkpt
           DO i=1,n
              eig_s(i,ik)=eig(i,ik)
           ENDDO
        ENDDO
        WRITE(6,*) 'nbndsi=',n
        IF( mpime == 0 ) WRITE (io) ((eig_s(i,ik),ik=1,nkpt),i=1,n)
     ELSE
        WRITE(6,*) 'nbndsi=',n
        IF( mpime == 0 ) WRITE (io) ((eig(i,ik),ik=1,nkpt),i=1,n)
     ENDIF
  ENDIF

!  write(6,*) 'autovalori energia per 10bande e tutti kpt'
!  WRITE(6,'(10F10.7)') ( ( eig(i,ik)*27.21 ,ik=1,nkpt), i=1,10 )
!  DEALLOCATE (eig_s, eig)
  !
  ! occupation numbers
  !
  ALLOCATE (focc(n,nkpt), focc_s(n,nkpt))
!  focc(:,:)   = wg(:,1:nkpt)
!  focc_s(:,:) = wg(:,1:nkpt)
  DO j=1,n
   DO ik=1,nkpt
    focc(j,ik)=wg(j,ik)*2.0d0/wk(ik)
    ENDDO
  ENDDO

  focc_s(:,:) = focc(:,:)

  IF( mpime == 0 ) THEN
     IF (t_form) THEN
        WRITE (io,'(3E26.18)') ((focc(i,ik), ik=1,nkpt), i=1,n)
     ELSE
        IF(t_single) THEN
           WRITE (io) ((focc_s(i,ik),ik=1,nkpt),i=1,n)
        ELSE
           WRITE (io) ((focc(i,ik),ik=1,nkpt),i=1,n)
        ENDIF
     ENDIF
  ENDIF

  WRITE (6,*) nkpt
  WRITE (6,*) 'weights:'
  WRITE (6,'(10f10.7)') (wk(ik), ik=1,nkpt)

  DO ik = 1, nkpt
     WRITE (6,*) 'ik=', ik
     WRITE (6,'(10f10.7)') (focc_s(j,ik), j=1,n)
  ENDDO

!  DEALLOCATE (focc_s, focc)
  !

  !  omax = nbnd*6
  omax = 600
  WRITE(6,*) 'io sono omax = ', omax
  alpha = omegamax/omax
  WRITE(6,*) 'alpha = ', alpha
  halfalpha= alpha*.5
  WRITE(6,*) 'halfalpha = ', halfalpha
  const = 4.d0*pi**2*AUTOEV**3/omega
  WRITE(6,*) 'const = ', const

  WRITE(*,*) "sono qui 6"
  ALLOCATE( omeg(omax+1))
  ALLOCATE( epsx(nkpt,omax+1), epsy(nkpt,omax+1), epsz(nkpt,omax+1) )
  ALLOCATE( epstx(omax+1), epsty(omax+1), epstz(omax+1) )
  ALLOCATE( pp1(nkpt,omax+1), pp2(nkpt,omax+1), pp3(nkpt,omax+1), omegatt(omax+1) )

  DO o = 1, omax + 1
     omeg(o) = (o-1)*alpha
  ENDDO

  epstx(:)=0.0
  epsty(:)=0.0
  epstz(:)=0.0
  epsx(:,:)=0.0
  epsy(:,:)=0.0
  epsz(:,:)=0.0
  pp1(:,:)=0.0
  pp2(:,:)=0.0
  pp3(:,:)=0.0

  WRITE(6,*) pp1(1,1), epstx(1)

  ALLOCATE ( c0(igwx), c0_s(igwx), kpg(3,igwx), c0_m(igwx,n), c0_tmp(igwxx) )
  !IF (gamma_only) ALLOCATE ( c0_gamma(2*igwx-1), c0_gamma_s(2*igwx-1) )
  CALL cryst_to_cart (nks, xk, bg, +1)
  IF (ionode) WRITE(6,*) 'Costruisco le psi ed il matrixelements'
  IF (ionode) WRITE(6,*) 'Controllo: I punti k ora devo essere in coordinate cartesiane!'
  IF (ionode) WRITE(6,'(1x,3f10.6)') ( (xk(i,ik),i=1,3), ik=1,nkpt)


  DO ik = 1, nkpt
    !
    npw = ngk(ik)
    ALLOCATE( igk_l2g( npw ) )
    !
    DO ig = 1, npw
       !
       igk_l2g(ig) = ig_l2g(igk_k(ig,ik))
       !
    ENDDO
    !
    IF( use_gmaps ) THEN
      !
      c0_m = 0.0d0
      !
      CALL read_and_collect( c0_m, size( c0_m, 1 ), n, ik )
      !
    ELSE
      !
      CALL davcio( evc, 2*nwordwfc, iunwfc, ik, -1 )
      !
      ! copy coefficient from array evc(:,n) (ordered as |k+G|)
      ! into array c0 with |G| ordering
      !
      DO ig = 1, npw
         IF( igk_k(ig,ik) < 1 .or. igk_k(ig,ik) > size( c0 ) ) &
            CALL errore(' pw2gw ', ' c0 too small ', 1 )
      ENDDO

      ! read wavefunctions and write the matrixelemnts

      DO i = 1, n

        ALLOCATE( c0_tmp_dp( igwxx ) )

        CALL mergewf( evc(:,i), c0_tmp_dp, npw, igk_l2g(:), mpime, nproc, 0, intra_image_comm )
        !
        ! important: missing components must be set to zero
        c0 (:) = 0.d0
        DO ig=1,npw
          c0(igk_k(ig,ik)) = evc(ig,i)
        ENDDO
        c0_m(:,i)=c0(:)

        c0_tmp = c0_tmp_dp
        IF( mpime == 0 ) WRITE(io) c0_tmp ! c0_s

        DEALLOCATE( c0_tmp_dp )

      ENDDO
    ENDIF

    DEALLOCATE( igk_l2g )

     ! k + g thet must be in 2piba units
     kpg(:,:) = 0.d0
     DO ig=1,npw
        kpg(:,igk_k(ig,ik))= xk_s(:,ik)+g(:,igk_k(ig,ik))
     ENDDO

     DO iband1 = 1,n
        IF ( focc(iband1,ik)>=1e-4) THEN
           DO iband2 = 1,n
              delta=2.0d0-focc(iband2,ik)
              IF (delta>1e-4) THEN

                 rhotwx = 0.0
                 DO ig=1,igwx
                    xkgk(1)= kpg(1,ig)
                    xkgk(2)= kpg(2,ig)
                    xkgk(3)= kpg(3,ig)
                    ctemp= conjg(c0_m(ig,iband1))*c0_m(ig,iband2)
                    rhotwx(1) = rhotwx(1) + xkgk(1) * ctemp
                    rhotwx(2) = rhotwx(2) + xkgk(2) * ctemp
                    rhotwx(3) = rhotwx(3) + xkgk(3) * ctemp
                 ENDDO

                 CALL mp_sum( rhotwx, world_comm )

                 IF (mpime == 0) THEN
                    rrhotwx(1)=tpiba2* real(rhotwx(1)*conjg(rhotwx(1)))
                    rrhotwx(2)=tpiba2* real(rhotwx(2)*conjg(rhotwx(2)))
                    rrhotwx(3)=tpiba2* real(rhotwx(3)*conjg(rhotwx(3)))
                    WRITE (90,'(1x,3i5,3e16.8,2f8.4)') ik,iband1,iband2,rrhotwx(1),rrhotwx(2), &
                    rrhotwx(3),(eig(iband2,ik)-eig(iband1,ik))*AUTOEV, (focc(iband1,ik)-focc(iband2,ik))
                    egap = (eig(iband2,ik)-eig(iband1,ik))*AUTOEV
                    Df = focc(iband1,ik)-focc(iband2,ik)
                    IF (egap>1e-3.and.Df>1e-4) THEN
                       DO o=1, omax+1
                          dummy = abs(egap - omeg(o))
                          IF (dummy<halfalpha) THEN
                             pp1(ik,o)=pp1(ik,o)+rrhotwx(1)*Df/egap**2
                             pp2(ik,o)=pp2(ik,o)+rrhotwx(2)*Df/egap**2
                             pp3(ik,o)=pp3(ik,o)+rrhotwx(3)*Df/egap**2
                          ENDIF
                       ENDDO
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
     ENDDO

  ENDDO

! CALCULATE POTENTIAL MATRIX ELEMNTS

   OPEN (UNIT=313,FILE="vxcdiag.dat",STATUS="UNKNOWN")
   WRITE(313,*) "#         K            BND          <Vxc>"
   ALLOCATE ( vxc(dfftp%nnr,nspin) )
   CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxc)
   DO ik=1,nkpt
      npw = ngk(ik)
      CALL davcio( evc, 2*nwordwfc, iunwfc, ik, -1 )
      DO iband1 = 1, nbnd
         psic(:) = (0.d0, 0.d0)
         DO ig = 1, npw
            psic(nl(igk_k(ig,ik)))  = evc(ig,iband1)
         ENDDO

         CALL invfft ('Dense', psic, dfftp)
         vxcdiag = 0.0d0
         !norma = 0.0d0
         DO ir = 1, dfftp%nnr
            vxcdiag = vxcdiag + vxc(ir,nspin) * &
                      ( dble(psic (ir) ) **2 + aimag(psic (ir) ) **2)
         !   norma = norma + ( DBLE(psic (ir) ) **2 + AIMAG(psic (ir) ) **2) / (nr1*nr2*nr3)
         ENDDO
 ! PG: this is the correct integral - 27/8/2010
         vxcdiag = vxcdiag * rytoev / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
         CALL mp_sum( vxcdiag, world_comm ) !, intra_pool_comm )
         ! ONLY FOR DEBUG!
         !IF (norma /= 1.0) THEN
         !   WRITE(*,*) "norma =", norma
         !   WRITE(*,*) "nrxx  =", nrxx
         !   STOP
         !ENDIF
         WRITE(313,"(i1,2x,i1,4x,f18.14)") ik, iband1, vxcdiag
      ENDDO
   ENDDO
   DEALLOCATE ( vxc )
   CLOSE (313)

!
  !
  IF ( mpime == 0 ) THEN

     WRITE(6, * )  ' out from k-points loop'
     WRITE(6, * )  ' Starting writing epsx,y,z,tot'

     WRITE(6,*) pp1(1,100)
     WRITE(6,*) pp1(1,350)

     OPEN (91, FILE='epsX.dat',STATUS='unknown',FORM='FORMATTED')
     OPEN (92, FILE='epsY.dat',STATUS='unknown',FORM='FORMATTED')
     OPEN (93, FILE='epsZ.dat',STATUS='unknown',FORM='FORMATTED')
     OPEN (94, FILE='epsTOT.dat',STATUS='unknown',FORM='FORMATTED')

     DO ik = 1, nkpt
        DO o =2, omax+1
            epsx(ik,o) = const * pp1(ik,o)*wk(ik)*0.5/ alpha
            epsy(ik,o) = const * pp2(ik,o)*wk(ik)*0.5/ alpha
            epsz(ik,o) = const * pp3(ik,o)*wk(ik)*0.5/ alpha
        ENDDO
     ENDDO

     WRITE(6, * ) epsx(1,150),epsx(1,300)

     DO o = 2, omax + 1
        omegatt(o) = (omeg(o-1)+omeg(o))*0.5
        DO ik = 1, nkpt
           epsxx= (epsx(ik,o-1)+epsx(ik,o))*0.5
           epsyy= (epsy(ik,o-1)+epsy(ik,o))*0.5
           epszz= (epsz(ik,o-1)+epsz(ik,o))*0.5
           epstx(o)=epstx(o)+epsxx
           epsty(o)=epsty(o)+epsyy
           epstz(o)=epstz(o)+epszz
        ENDDO
        WRITE(91,"(f15.6,1x,f15.6)") omegatt(o), epstx(o)
        WRITE(92,"(f15.6,1x,f15.6)") omegatt(o), epsty(o)
        WRITE(93,"(f15.6,1x,f15.6)") omegatt(o), epstz(o)
        WRITE(94,"(f15.6,1x,f15.6)") omegatt(o), (epstx(o)+ epsty(o)+ epstz(o))/3.0
     ENDDO

     WRITE(6, * )  ' Hey bello sto a fini'

     CLOSE(91)
     CLOSE(92)
     CLOSE(93)
     CLOSE(94)

  ENDIF

  DEALLOCATE (xk_s)
  DEALLOCATE (eig_s, eig)
  DEALLOCATE (focc_s, focc)
  DEALLOCATE (c0_s, c0, kpg, c0_m)
  DEALLOCATE (omeg, pp1,pp2, pp3, omegatt)
  DEALLOCATE ( epsx, epsy, epsz )
  DEALLOCATE ( epstx, epsty, epstz )
  !
  IF( mpime == 0 ) CLOSE(io)
  IF( mpime == 0 ) CLOSE(90)
  !
END SUBROUTINE compute_gw


!-----------------------------------------------------------------------
SUBROUTINE write_gmaps ( kunit)
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout
  USE cell_base, ONLY : at, bg, alat
  USE ions_base, ONLY : atm, nat
  USE gvect,     ONLY : ngm, ngm_g, ig_l2g, g
  USE lsda_mod,  ONLY : nspin, isk
  USE ions_base, ONLY : ntyp => nsp, tau, ityp
  USE wvfct,     ONLY : nbnd, npwx, et
  USE gvecw,     ONLY : gcutw
  USE klist,     ONLY : nkstot, ngk, nks, xk
  USE wavefunctions_module,  ONLY : evc
  USE io_files,  ONLY : nd_nmbr, tmp_dir, prefix, iunwfc, nwordwfc
  USE io_global, ONLY : ionode
  USE mp_images, ONLY : my_image_id
  USE mp_global, ONLY : nproc_pool, my_pool_id, my_image_id, intra_pool_comm
  USE mp,        ONLY : mp_sum, mp_max
  USE mp_world,  ONLY : world_comm, nproc, mpime


  IMPLICIT NONE
  INTEGER :: kunit

  INTEGER :: npw, i, j, k, ig, ik, ibnd, na, ngg, ikw
  INTEGER, ALLOCATABLE :: kisort(:)
  INTEGER :: ike, iks, npw_g, npwx_g, ispin
  INTEGER, EXTERNAL :: global_kpoint_index
  INTEGER, ALLOCATABLE :: ngk_g( : )
  INTEGER, ALLOCATABLE :: ngk_gw( : )
  INTEGER, ALLOCATABLE :: itmp( :, : )
  INTEGER, ALLOCATABLE :: igwk( : )
  INTEGER, ALLOCATABLE :: igk_l2g( :, : )
  REAL(kind=8), ALLOCATABLE :: gk(:)

  real(kind=8) :: wfc_scal
  LOGICAL :: twf0, twfm, twrite_wfc

  !
  !
  IF( ionode ) WRITE( stdout, fmt="(//,'WRITING G-MAPS for each processor' )" )

  IF( nkstot > 0 ) THEN

     IF( ( kunit < 1 ) .or. ( mod( nkstot, kunit ) /= 0 ) ) &
       CALL errore( ' write_wannier ',' wrong kunit ', 1 )

     IF( ( nproc_pool > nproc ) .or. ( mod( nproc, nproc_pool ) /= 0 ) ) &
       CALL errore( ' write_wannier ',' nproc_pool ', 1 )

     iks = global_kpoint_index (nkstot, 1)
     ike = iks + nks - 1

  ENDIF

  ! find out the global number of G vectors: ngm_g
  ngm_g = ngm
  CALL mp_sum( ngm_g, intra_pool_comm )


  ! build the G+k array indexes
  ALLOCATE ( igk_l2g( npwx, ik ) )
  ALLOCATE ( kisort( npwx ) )
  ALLOCATE ( gk( npwx ) )
  DO ik = 1, nks
     kisort = 0
     CALL gk_sort (xk (1, ik+iks-1), ngm, g, gcutw, npw, kisort(1), gk)
     DO ig = 1, npw
        igk_l2g(ig,ik) = ig_l2g(kisort(ig))
     ENDDO
     ngk (ik) = npw
  ENDDO
  DEALLOCATE (gk, kisort)

  ! compute the global number of G+k vectors for each k point
  ALLOCATE( ngk_g( nkstot ) )
  ALLOCATE( ngk_gw( nkstot/nspin ) )
  ngk_g = 0
  ngk_g( iks:ike ) = ngk( 1:nks )
  CALL mp_sum( ngk_g, world_comm )

  ! compute the Maximum G vector index among all G+k an processors
  npw_g = maxval( ig_l2g(:) ) ! ( igk_l2g(:,:) )
  CALL mp_max( npw_g, world_comm )

  ! compute the Maximum number of G vector among all k points
  npwx_g = maxval( ngk_g( 1:nkstot ) )


  ALLOCATE( igwk( npwx_g ) )

  DO ik = 1, nkstot
    igwk = 0
    ALLOCATE( itmp( npw_g, 1 ) )
    itmp = 0
    IF( ik >= iks .and. ik <= ike ) THEN
      DO  ig = 1, ngk( ik-iks+1 )
        itmp( ig_l2g( ig ), 1 ) = ig_l2g( ig )
      ENDDO
    ENDIF
    CALL mp_sum( itmp, world_comm )
    ngg = 0
    DO  ig = 1, npw_g
      IF( itmp( ig, 1 ) == ig ) THEN
        ngg = ngg + 1
        igwk( ngg ) = ig
      ENDIF
    ENDDO
    IF( ngg /= ngk_g( ik ) ) THEN
      WRITE( stdout,*) ' ik, ngg, ngk_g = ', ik, ngg, ngk_g( ik )
    ENDIF
    DEALLOCATE( itmp )
    IF( ionode ) THEN
        ! write (40)( igwk(ig), ig = 1, npwx_g )
    ENDIF
  ENDDO

  DEALLOCATE( igwk )

  DO ik = 1, nkstot
     IF( (ik >= iks) .and. (ik <= ike) ) THEN
        ispin = isk( ik )
        WRITE( 100 + mpime ) ik, iks, ike, nkstot, kunit, nproc, ispin, nspin, npw_g, &
                             nbnd, ngk(ik-iks+1), 2*nwordwfc, npwx, iunwfc, nd_nmbr
        WRITE( 100 + mpime ) ( igk_l2g( i, ik-iks+1 ), i = 1, ngk(ik-iks+1) )
     ENDIF
  ENDDO

  DEALLOCATE ( ngk_g )
  DEALLOCATE ( ngk_gw )
  DEALLOCATE (igk_l2g)

END SUBROUTINE write_gmaps


SUBROUTINE read_and_collect( c, ldc, n, ik )
  USE io_global,      ONLY : stdout
  USE io_files,       ONLY : prefix
  USE kinds,          ONLY : DP, sgl

  IMPLICIT NONE

  INTEGER     :: ldc, n, ik
  COMPLEX(DP) :: c( ldc, n )
  INTEGER     :: ik_ , iks, ike, nkstot, kunit, nproc_ , ispin, nspin, npw_g , nbnd , ngk
  INTEGER     :: nwordwfc, npwx, iunwfc
  INTEGER     :: nfile, ip, i, j
  COMPLEX(DP), ALLOCATABLE :: evc( :, : )
  INTEGER, ALLOCATABLE :: igk_l2g( : )
  LOGICAL     :: exst
  CHARACTER(len=3)  :: nd_nmbr

  READ( 100 ) ik_ , iks, ike, nkstot, kunit, nproc_ , ispin, nspin, npw_g , &
              nbnd , ngk, nwordwfc, npwx, iunwfc, nd_nmbr

  REWIND( 100 )

  nfile = nproc_

  CLOSE( iunwfc )

  DO ip = 0, nfile - 1
    READ( 100 + ip ) ik_ , iks, ike, nkstot, kunit, nproc_ , ispin, nspin, npw_g , &
                     nbnd , ngk, nwordwfc, npwx, iunwfc, nd_nmbr
    WRITE( stdout, * ) 'DEBUG nd_nmbr ', nd_nmbr
    IF( ( ik_ == ik ) .and. ( ik_ >= iks ) .and. ( ik_ <= ike ) ) THEN
      ALLOCATE( evc( npwx, nbnd ) )
      ALLOCATE( igk_l2g( ngk ) )
      READ( 100 + ip )  ( igk_l2g( i ), i = 1, ngk )
      CALL diropn_gw ( 99, trim( prefix )//'.wfc', 2*nwordwfc, exst, ip, nd_nmbr )
      CALL davcio ( evc, 2*nwordwfc, 99, (ik-iks+1), - 1 )
      CLOSE( 99 )
      DO j = 1, n
        DO i = 1, ngk
          c( igk_l2g( i ), j ) = evc( i, j )
        ENDDO
      ENDDO
      DEALLOCATE( evc )
      DEALLOCATE( igk_l2g )
    ENDIF
    REWIND( 100 + ip )
  ENDDO

  RETURN
END SUBROUTINE

!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE diropn_gw (unit, filename, recl, exst, mpime, nd_nmbr_ )
  !-----------------------------------------------------------------------
  !
  !     this routine opens a file in tmp_dir for direct I/O access
  !     If appropriate, the node number is added to the file name
  !
  USE kinds
  USE io_files
  IMPLICIT NONE

  !
  !    first the input variables
  !
  CHARACTER(len=*) :: filename
  ! input: name of the file to open
  INTEGER :: unit, recl
  ! input: unit of the file to open
  ! input: length of the records
  LOGICAL :: exst
  ! output: if true the file exists
  INTEGER :: mpime
  ! input: processor index
  CHARACTER(len=3) :: nd_nmbr_
  !
  !    local variables
  !
  CHARACTER(len=256) :: tempfile
  ! complete file name
  CHARACTER(len=80) :: assstr
  INTEGER :: ios, unf_recl, ierr
  ! used to check I/O operations
  ! length of the record
  ! error code
  LOGICAL :: opnd
  ! if true the file is already opened


  IF (unit < 0) CALL errore ('diropn', 'wrong unit', 1)
  !
  !    we first check that the file is not already openend
  !
  ios = 0
  INQUIRE (unit = unit, opened = opnd)
  IF (opnd) CALL errore ('diropn', "can't open a connected unit", abs(unit))
  !
  !      then we check the filename
  !

  IF (filename == ' ') CALL errore ('diropn', 'filename not given', 2)
  tempfile = trim(tmp_dir) // trim(filename) // trim( nd_nmbr_ )

  INQUIRE (file = tempfile, exist = exst)
  !
  !      the unit for record length is unfortunately machine-dependent
  !
#define DIRECT_IO_FACTOR 8
  unf_recl = DIRECT_IO_FACTOR * recl
  IF (unf_recl <= 0) CALL errore ('diropn', 'wrong record length', 3)
  !
  OPEN ( unit, file = trim(tempfile), iostat = ios, form = 'unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl )

  IF (ios /= 0) CALL errore ('diropn', 'error opening '//filename, unit)
  RETURN
END SUBROUTINE diropn_gw

!----------------------------------------------------------------------
subroutine gen_us_djl (ik,npw,djl,size_tab,vec_tab, spline_ps, vec_tab_d2y)
  !----------------------------------------------------------------------
  !
  !  Calculates the kleinman-bylander pseudopotentials with the
  !  derivative of the spherical harmonics projected on vector u
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE constants,  ONLY : tpi
  USE cell_base,  ONLY : tpiba
  USE klist,      ONLY : xk, igk_k
  USE gvect,      ONLY : g
  USE us,         ONLY : nqx, dq
  USE splinelib,  ONLY : splint_deriv
  USE uspp_param, ONLY : upf
  !
  implicit none
  !
  integer, intent(in) :: ik, npw
  real(DP), intent(inout) ::djl(1:npw)
  integer, intent(in) :: size_tab
  real(DP), intent(in) :: vec_tab(1:size_tab)
  real(DP), intent(in) :: vec_tab_d2y(1:size_tab)
  logical :: spline_ps
  !
  integer :: i0, i1, i2, &
       i3, ig
  real(DP), allocatable :: gk(:,:), q (:)
  real(DP) :: px, ux, vx, wx
  complex(DP), allocatable :: sk (:)

  integer :: iq
  real(DP), allocatable :: xdata(:)
  real(DP) :: qt


  allocate ( gk(3,npw) )
  allocate ( q(npw) )

  do ig = 1, npw
     gk (1, ig) = xk (1, ik) + g (1, igk_k (ig,ik) )
     gk (2, ig) = xk (2, ik) + g (2, igk_k (ig,ik) )
     gk (3, ig) = xk (3, ik) + g (3, igk_k (ig,ik) )
     q (ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo

  do ig = 1, npw
     q (ig) = sqrt ( q(ig) ) * tpiba
  end do

  if (spline_ps) then
    allocate(xdata(nqx))
    do iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    enddo
  endif

  ! calculate beta in G-space using an interpolation table
  do ig = 1, npw
    qt = sqrt(q(ig)) * tpiba
    if (spline_ps) then
        djl(ig) = splint_deriv(xdata, vec_tab(:), &
                                vec_tab_d2y(:), qt)
    else
        px = qt / dq - int (qt / dq)
        ux = 1.d0 - px
        vx = 2.d0 - px
        wx = 3.d0 - px
        i0 = qt / dq + 1
        i1 = i0 + 1
        i2 = i0 + 2
        i3 = i0 + 3
        djl (ig) = vec_tab (i0) * (-vx*wx-ux*wx-ux*vx) / 6.d0 + &
                          vec_tab (i1) * (+vx*wx-px*wx-px*vx) / 2.d0 - &
                          vec_tab (i2) * (+ux*wx-px*wx-px*ux) / 2.d0 + &
                          vec_tab (i3) * (+ux*vx-px*vx-px*ux) / 6.d0
    endif
  enddo

  deallocate (q)
  deallocate ( gk )
  if (spline_ps) deallocate(xdata)
  return
end subroutine gen_us_djl
!
!----------------------------------------------------------------------
subroutine gen_us_vkb0 (ik,npw,vkb0,size_tab,vec_tab, spline_ps, vec_tab_d2y)
  !----------------------------------------------------------------------
  !
  !  Calculates the kleinman-bylander pseudopotentials with the
  !  derivative of the spherical harmonics projected on vector u
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE constants,  ONLY : tpi
  USE cell_base,  ONLY : tpiba
  USE klist,      ONLY : xk, igk_k
  USE gvect,      ONLY : g
  USE us,         ONLY : nqx, dq
  USE splinelib,  ONLY : splint
  USE uspp_param, ONLY : upf
  !
  implicit none
  !
  integer, intent(in) :: ik, npw
  real(DP), intent(inout) ::vkb0(1:npw)
  integer, intent(in) :: size_tab
  real(DP), intent(in) :: vec_tab(1:size_tab)
  real(DP), intent(in) :: vec_tab_d2y(1:size_tab)
  logical :: spline_ps
  !
  integer :: na, nt, nb, ikb,i0, i1, i2, &
       i3, ig
  real(DP), allocatable :: gk(:,:), q (:)
  real(DP) :: px, ux, vx, wx
  complex(DP), allocatable :: sk (:)

  integer :: iq
  real(DP), allocatable :: xdata(:)

  allocate ( gk(3,npw) )
  allocate ( q(npw) )

  do ig = 1, npw
     gk (1, ig) = xk (1, ik) + g (1, igk_k (ig,ik) )
     gk (2, ig) = xk (2, ik) + g (2, igk_k (ig,ik) )
     gk (3, ig) = xk (3, ik) + g (3, igk_k (ig,ik) )
     q (ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo

  do ig = 1, npw
     q (ig) = sqrt ( q(ig) ) * tpiba
  end do

  if (spline_ps) then
    allocate(xdata(nqx))
    do iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    enddo
  endif

  ! calculate beta in G-space using an interpolation table
  do ig = 1, npw
    if (spline_ps) then
        vkb0(ig) = splint(xdata, vec_tab(:), &
                                vec_tab_d2y(:), q(ig))
    else
        px = q (ig) / dq - int (q (ig) / dq)
        ux = 1.d0 - px
        vx = 2.d0 - px
        wx = 3.d0 - px
        i0 = q (ig) / dq + 1
        i1 = i0 + 1
        i2 = i0 + 2
        i3 = i0 + 3
        vkb0 (ig) = vec_tab (i0) * ux * vx * wx / 6.d0 + &
                          vec_tab (i1) * px * vx * wx / 2.d0 - &
                          vec_tab (i2) * px * ux * wx / 2.d0 + &
                          vec_tab (i3) * px * ux * vx / 6.d0
    endif
  enddo

  deallocate (q)
  deallocate ( gk )
  if (spline_ps) deallocate(xdata)
  return
end subroutine gen_us_vkb0
