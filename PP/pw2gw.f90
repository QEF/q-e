
! Copyright (C) 2004-2007 Quantum-Espresso group 
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 
! 13Aprile2005 
! GENERATES INPUT for GW code
!tested on: Silicon bulk, Germanium Bulk, Na4, InP bulk
! Please note just symmorphic symm. op. have to be used
! to get rid of non symmorphic symm. op. set fractional_translations =.false.
! in sgama_at.f90



!----------------------------------------------------------------------- 
PROGRAM pw2gw
  !----------------------------------------------------------------------- 

  ! This subroutine writes files containing plane wave coefficients
  ! and other stuff needed by GW codes

  USE io_files,  ONLY : nd_nmbr, prefix, outdir, tmp_dir, trimcheck
  USE io_global, ONLY : ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  USE mp_global, ONLY : kunit, nproc
  !
  IMPLICIT NONE
  INTEGER :: ios
  INTEGER :: kunittmp
  LOGICAL :: use_gmaps
  CHARACTER(LEN=20) :: what
  CHARACTER(LEN=30) :: when

  NAMELIST / inputpp / prefix, outdir, what, use_gmaps

  CALL start_postproc(nd_nmbr)
  ! 
  !   set default values for variables in namelist 
  ! 
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  what   = 'gw'
  use_gmaps = .false.

  ios = 0
  IF ( ionode )  THEN 
     !
     READ (5, inputpp, iostat=ios)
     tmp_dir = trimcheck (outdir)
     !
  END IF
  ! 
  CALL mp_bcast( ios, ionode_id ) 
  IF (ios /= 0)   CALL errore('pw2gw', 'reading inputpp namelist', ABS(ios))
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( prefix, ionode_id ) 
  CALL mp_bcast(tmp_dir, ionode_id ) 
  CALL mp_bcast( what, ionode_id )
  CALL mp_bcast( use_gmaps, ionode_id )
  !
  CALL read_file
  CALL openfil_pp
  !
#if defined __PARA
  kunittmp = kunit
#else
  kunittmp = 1
#endif
  !
  IF( TRIM( what ) == 'gw' ) THEN
    CALL compute_gw  ( use_gmaps )
  ELSE
    CALL write_gmaps ( kunittmp )
  END IF
  !
  CALL stop_pp

END PROGRAM pw2gw


SUBROUTINE compute_gw( use_gmaps )

  ! This routine creates the QPLDA and the matrixelements
  ! tform = .false. UNFORMATTED QPLDA
  ! tform = .true.  FORMATTED QPLDA
  ! tsingle must be always true 

  USE kinds,     ONLY : DP, sgl
  USE constants, ONLY : eps8, pi, AUTOEV
  USE cell_base, ONLY : alat, tpiba2, at, bg, omega
  USE char,      ONLY : title
  USE symme,     ONLY : s, nsym
  USE wvfct,     ONLY : npw, npwx, nbnd, igk, g2kin, wg, et
  USE control_flags, ONLY : gamma_only
  USE gvect,         ONLY : ngm, g, gg, ig_l2g, ecutwfc
  USE klist ,        ONLY : nks, xk, wk
  USE lsda_mod,      ONLY : nspin
  USE io_files,      ONLY : nwordwfc, iunwfc
  USE wavefunctions_module, ONLY : evc
  use mp_global, ONLY : mpime, kunit, nproc, intra_image_comm, npool
  USE io_global, ONLY : ionode, ionode_id
  USE mp,        ONLY : mp_sum , mp_max
  USE mp_wave,   ONLY : mergewf
  USE parallel_include

  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: use_gmaps

  INTEGER :: ii(16), ngw, nkpt, ig, ik, n, i,j,k, io = 98, iband1, iband2
  INTEGER :: omax, o, iproc
  INTEGER, ALLOCATABLE :: in1(:), in2(:), in3(:)
  INTEGER, ALLOCATABLE :: in1_tmp(:), in2_tmp(:), in3_tmp(:)
  INTEGER, ALLOCATABLE :: inx_rcv(:), ig_l2g_rcv(:)
  LOGICAL :: t_form = .false., t_single = .true.
  REAL(kind=sgl) :: a1_s(3), a2_s(3), a3_s(3)
  REAL(kind=sgl), ALLOCATABLE :: xk_s(:,:), eig_s(:,:), focc_s(:,:)
  REAL(kind=DP):: g2max, a1(3), a2(3), a3(3),norm, xkgk(3), rrhotwx(3), delta
  REAL(kind=DP):: alpha, egap, halfalpha, Df, const, dummy
  REAL(kind=DP), parameter :: omegamax = 30.0
  REAL(kind=DP), ALLOCATABLE:: gsort(:), eig(:,:), focc(:,:), kpg(:,:), omegatt(:), omeg(:)
  REAL(kind=DP), ALLOCATABLE:: pp1(:,:), pp2(:,:), pp3(:,:)
  REAL(kind=DP), ALLOCATABLE:: epsx(:,:), epsy(:,:), epsz(:,:)
  REAL(kind=DP), ALLOCATABLE:: epstx(:), epsty(:), epstz(:)
  REAL(kind=DP) :: epsxx, epsyy, epszz
  COMPLEX(kind=DP):: rhotwx(3), ctemp, dasomma(3)
  COMPLEX(kind=DP),  ALLOCATABLE:: c0(:), c0_m(:,:), c0_tmp_dp(:) !, c0_tmp(:) !, c0_gamma(:)
  COMPLEX(kind=sgl), ALLOCATABLE:: c0_s(:), c0_tmp(:) !, c0_gamma_s(:)
  CHARACTER(LEN=80) :: titleo(2)
  INTEGER :: igwx, igwxx, comm, ierr, ig_max, igwx_r
  INTEGER :: igwx_p(nproc)
  INTEGER, ALLOCATABLE :: igk_l2g(:)
  !
#if defined __PARA
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
  END IF
  !
  !  file's title [2 lines]
  !
  titleo(1)='pw2gw'
  titleo(2)='test version'
  IF( mpime == 0 ) THEN
     IF (t_form) THEN
        write (io,'(A80/A80)') titleo(1), titleo(2)
     ELSE
        write (io) titleo(1)
        write (io) titleo(2)
     ENDIF
     !
     write(6,*) 'qplda title'
     write(6,*) titleo(1)
     write(6,*) titleo(2)
  END IF
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
     IF( mpime == 0 ) write (io,'(16I5)') ii
  ELSE
     ii(1)=1
     IF( mpime == 0 ) write (io) ii
  END IF
  !
  write(6,'(16I5)') ii
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
     write(6,'(a,3E26.18)') 'a1', a1_s
     write(6,'(a,3E26.18)') 'a2', a2_s
     write(6,'(a,3E26.18)') 'a3', a3_s
     !
  END IF
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
     write(6,*)'nrot=',nsym
     write(6,'(3E26.18)') (((float(s(i,j,k)),j=1,3),i=1,3),k=1,nsym)
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
           WRITE (io) (((dfloat(s(i,j,k)),j=1,3),i=1,3),k=1,nsym)
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
  g2kin(:) = 0.0d0
  !DEBUG
  IF (ionode) WRITE(6,*) ' nks ', nks
  IF (ionode) WRITE(6,*) ' k points in  cartesian coordinates'
  IF (ionode) WRITE(6,'(1x,3f10.6)') ( (xk(i,ik),i=1,3), ik=1,nks)
  !DEBUG
  igwx  = 0  !  maximum G vector index
  DO ik = 1, nks
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     g2max = MAX ( g2max, MAXVAL (g2kin(1:npw)) )  
     ! WRITE( 6, * ) 'DEBUG g2max ', g2max 
     ! g2max, g2kin = RAGGIO DELLA SFERA |G+k|<cut, non MASSIMO |G| nella sfera
     ! g2max <= ecutwfc / tpiba2   PER COSTRUZIONE
     igwx = MAX( igwx, MAXVAL( igk(1:npw) ) ) 
  END DO
  IF (ionode) write(*,*) "igwx = ", igwx
  !
  !  ngw =  number of G-vectors (complete shells) such that G2 <= G2max
  !  ovvero <= RAGGIO della SFERA, in pratica trova i G2 relativi a GAMMA
  !
  do ngw = 1, ngm
     if ( gg(ngw) > g2max + eps8) go to 100
  end do
  call errore ( 'pw2gw','max G in PW not found?!?',ngw)
100 ngw = ngw - 1

  ! Pongo NGW pari al massimo indice tra i vettori G che fanno parte delle 
  ! sfere |G+k|<cut per qualsiasi k
  !
  IF (ionode) write( 6, * ) ' igwx= ', igwx

  ngw = igwx

  !
  !  PWscf stores G in order of increasing module
  !
  ALLOCATE (in1(ngw), in2(ngw), in3(ngw))
  DO ig=1,ngw
     in1(ig) = NINT ( at(1,1)*g(1,ig) + at(2,1)*g(2,ig) + at(3,1)*g(3,ig) )
     in2(ig) = NINT ( at(1,2)*g(1,ig) + at(2,2)*g(2,ig) + at(3,2)*g(3,ig) )
     in3(ig) = NINT ( at(1,3)*g(1,ig) + at(2,3)*g(2,ig) + at(3,3)*g(3,ig) )
  END DO

  igwxx = MAXVAL( ig_l2g( 1:ngw ) )
  CALL mp_max( igwxx )

  igwx_p = 0
  igwx_p( mpime + 1 ) = igwx
  CALL mp_sum( igwx_p )

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
     END DO
     !
#if defined __PARA
     ALLOCATE( ig_l2g_rcv( igwxx ) )
     ALLOCATE( inx_rcv( igwxx ) )
     !
     DO iproc = 2, nproc
        CALL MPI_RECV( ig_l2g_rcv, igwx_p( iproc ), MPI_INTEGER, (iproc-1), iproc,       intra_image_comm, istatus, IERR )
        !
        CALL MPI_RECV( inx_rcv,    igwx_p( iproc ), MPI_INTEGER, (iproc-1), iproc+NPROC, intra_image_comm, istatus, IERR )
        DO ig = 1, igwx_p( iproc )
           in1_tmp( ig_l2g_rcv( ig ) ) = inx_rcv( ig )
        END DO
        CALL MPI_RECV( inx_rcv,    igwx_p( iproc ), MPI_INTEGER, (iproc-1), iproc+2*NPROC, intra_image_comm, istatus, IERR )
        DO ig = 1, igwx_p( iproc )
           in2_tmp( ig_l2g_rcv( ig ) ) = inx_rcv( ig )
        END DO
        CALL MPI_RECV( inx_rcv,    igwx_p( iproc ), MPI_INTEGER, (iproc-1), iproc+3*NPROC, intra_image_comm, istatus, IERR )
        DO ig = 1, igwx_p( iproc )
           in3_tmp( ig_l2g_rcv( ig ) ) = inx_rcv( ig )
        END DO
     END DO
     !
     DEALLOCATE( ig_l2g_rcv )
     DEALLOCATE( inx_rcv )
#endif
     !
  ELSE
     !
#if defined __PARA
     CALL MPI_SEND( ig_l2g, igwx, MPI_INTEGER, 0, mpime+1,         intra_image_comm, IERR )
     CALL MPI_SEND( in1(1), igwx, MPI_INTEGER, 0, mpime+1+NPROC,   intra_image_comm, IERR )
     CALL MPI_SEND( in2(1), igwx, MPI_INTEGER, 0, mpime+1+2*NPROC, intra_image_comm, IERR )
     CALL MPI_SEND( in3(1), igwx, MPI_INTEGER, 0, mpime+1+3*NPROC, intra_image_comm, IERR )
#endif
     !
  END IF

  IF (mpime == 0) write(*,*) "fine debug sui punti g"



  IF (t_form) THEN
     IF( mpime == 0 ) THEN
        WRITE (io,'(I12)') igwxx
        WRITE (io,'(3I5)') (in1_tmp(ig),in2_tmp(ig),in3_tmp(ig),ig=1,igwxx)
     END IF
  ELSE
     IF( mpime == 0 ) THEN
        WRITE (io) igwxx 
        WRITE (io) (in1_tmp(ig),in2_tmp(ig),in3_tmp(ig),ig=1,igwxx)
     END IF
  ENDIF
  !
  DEALLOCATE ( in1, in2, in3 )
  IF( mpime == 0 ) THEN
     DEALLOCATE ( in1_tmp, in2_tmp, in3_tmp )
  END IF
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
     WRITE(65,'(1x,3f10.6,x,f10.6)')  ( xk_s(:,ik), wk(ik)*0.5, ik=1,nks )
     CLOSE(unit=65)
  END IF

  IF( mpime == 0 ) THEN
     IF (t_form) THEN
        WRITE (io,'(I12)') nkpt
        WRITE (io,'(3E26.18)') ((xk_s(i,ik),i=1,3),ik=1,nkpt)
     ELSE
        WRITE (io) nkpt
        WRITE(6,*) 'nkpt',nkpt
        if(t_single) then
           WRITE (io) ((xk_s(i,ik),i=1,3),ik=1,nkpt)
        else
           WRITE (io) ((xk(i,ik),i=1,3),ik=1,nkpt)
        endif
     ENDIF
     WRITE(6,'(1x,3f10.6)') ( (xk_s(i,ik),i=1,3), ik=1,nkpt)
  ENDIF

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
     if(t_single) then
        do ik=1,nkpt
           do i=1,n
              eig_s(i,ik)=eig(i,ik)
           enddo
        enddo
        WRITE(6,*) 'nbndsi=',n
        IF( mpime == 0 ) WRITE (io) ((eig_s(i,ik),ik=1,nkpt),i=1,n)
     else
        WRITE(6,*) 'nbndsi=',n
        IF( mpime == 0 ) WRITE (io) ((eig(i,ik),ik=1,nkpt),i=1,n)
     endif
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
  do j=1,n
   do ik=1,nkpt
    focc(j,ik)=wg(j,ik)*2.0d0/wk(ik)
    enddo
  enddo
  
  focc_s(:,:) = focc(:,:)

  IF( mpime == 0 ) THEN
     IF (t_form) THEN
        WRITE (io,'(3E26.18)') ((focc(i,ik), ik=1,nkpt), i=1,n)
     ELSE
        if(t_single) then
           WRITE (io) ((focc_s(i,ik),ik=1,nkpt),i=1,n)
        else
           WRITE (io) ((focc(i,ik),ik=1,nkpt),i=1,n)
        endif
     ENDIF
  END IF

  WRITE (6,*) nkpt
  WRITE (6,*) 'weights:'
  WRITE (6,'(10f10.7)') (wk(ik), ik=1,nkpt)

  do ik = 1, nkpt
     WRITE (6,*) 'ik=', ik
     WRITE (6,'(10f10.7)') (focc_s(j,ik), j=1,n)
  enddo

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

  write(*,*) "sono qui 6"
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
    CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
    !
    ALLOCATE( igk_l2g( npw ) )
    !
    DO ig = 1, npw
       !
       igk_l2g(ig) = ig_l2g(igk(ig))
       !
    END DO
    !
    IF( use_gmaps ) THEN
      !
      c0_m = 0.0d0
      !
      CALL read_and_collect( c0_m, SIZE( c0_m, 1 ), n, ik )
      !
    ELSE
      !
      CALL davcio( evc, nwordwfc, iunwfc, ik, -1 )
      !
      ! copy coefficient from array evc(:,n) (ordered as |k+G|)
      ! into array c0 with |G| ordering
      !
      DO ig = 1, npw
         IF( igk(ig) < 1 .OR. igk(ig) > SIZE( c0 ) ) &
            CALL errore(' pw2gw ', ' c0 too small ', 1 )
      END DO

      ! read wavefunctions and write the matrixelemnts

      DO i = 1, n

        allocate( c0_tmp_dp( igwxx ) )

        CALL mergewf( evc(:,i), c0_tmp_dp, npw, igk_l2g(:), mpime, nproc, 0, intra_image_comm )
        !
        ! important: missing components must be set to zero
        c0 (:) = 0.d0
        DO ig=1,npw
          c0(igk(ig)) = evc(ig,i)
        END DO
        c0_m(:,i)=c0(:)

        c0_tmp = c0_tmp_dp
        IF( mpime == 0 ) WRITE(io) c0_tmp ! c0_s

        deallocate( c0_tmp_dp )

      ENDDO
    ENDIF

    DEALLOCATE( igk_l2g )
    
  ENDDO


  DO ik = 1, nkpt

     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)

     ! k + g thet must be in 2piba units
     kpg(:,:) = 0.d0
     DO ig=1,npw
        kpg(:,igk(ig))= xk_s(:,ik)+g(:,igk(ig))
     END DO

     DO iband1 = 1,n
        IF ( focc(iband1,ik).ge.1e-4) THEN
           DO iband2 = 1,n
              delta=2.0d0-focc(iband2,ik)
              IF (delta.gt.1e-4) THEN
   
                 rhotwx = 0.0
                 DO ig=1,igwx
                    xkgk(1)= kpg(1,ig)
                    xkgk(2)= kpg(2,ig)
                    xkgk(3)= kpg(3,ig)
                    ctemp= CONJG(c0_m(ig,iband1))*c0_m(ig,iband2)
                    rhotwx(1) = rhotwx(1) + xkgk(1) * ctemp
                    rhotwx(2) = rhotwx(2) + xkgk(2) * ctemp
                    rhotwx(3) = rhotwx(3) + xkgk(3) * ctemp
                 ENDDO
   
                 CALL mp_sum( rhotwx )
   
                 IF (mpime == 0) THEN
                    rrhotwx(1)=tpiba2* real(rhotwx(1)*conjg(rhotwx(1)))
                    rrhotwx(2)=tpiba2* real(rhotwx(2)*conjg(rhotwx(2)))
                    rrhotwx(3)=tpiba2* real(rhotwx(3)*conjg(rhotwx(3)))
                    WRITE (90,'(1x,3i5,3e16.8,2f8.4)') ik,iband1,iband2,rrhotwx(1),rrhotwx(2), &
                    rrhotwx(3),(eig(iband2,ik)-eig(iband1,ik))*AUTOEV, (focc(iband1,ik)-focc(iband2,ik))
                    egap = (eig(iband2,ik)-eig(iband1,ik))*AUTOEV
                    Df = focc(iband1,ik)-focc(iband2,ik)
                    IF (egap.gt.1e-3.and.Df.gt.1e-4) THEN
                       DO o=1, omax+1
                          dummy = abs(egap - omeg(o))
                          IF (dummy.lt.halfalpha) THEN
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
        write(91,"(f15.6,x,f15.6)") omegatt(o), epstx(o)
        write(92,"(f15.6,x,f15.6)") omegatt(o), epsty(o)
        write(93,"(f15.6,x,f15.6)") omegatt(o), epstz(o)
        write(94,"(f15.6,x,f15.6)") omegatt(o), (epstx(o)+ epsty(o)+ epstz(o))/3.0
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
subroutine write_gmaps ( kunit) 
  !----------------------------------------------------------------------- 
  ! 
  USE io_global, ONLY : stdout 
  USE cell_base, ONLY : at, bg, tpiba2, alat 
  USE ions_base, ONLY : atm, nat 
  USE gvect,     ONLY : ngm, ngm_g, ig_l2g, ig1, ig2, ig3, ecutwfc, & 
       nr1, nr2, nr3, g 
  USE lsda_mod,  ONLY : nspin, isk 
  USE ions_base, ONLY : ntyp => nsp, tau, ityp 
  USE wvfct,     ONLY : nbnd, npw, npwx, et, g2kin 
  USE gvect,     ONLY : ig_l2g
  USE klist,     ONLY : nkstot, ngk, nks, xk 
  USE wavefunctions_module,  ONLY : evc 
  use io_files,  only : nd_nmbr, tmp_dir, prefix, iunwfc, nwordwfc 
  use io_global, only : ionode 
  use mp_global, only : nproc, nproc_pool, mpime 
  use mp_global, only : my_pool_id, my_image_id, intra_pool_comm 
  use mp,        only : mp_sum, mp_max 
 
 
  implicit none 
  integer :: kunit 
 
  integer :: i, j, k, ig, ik, ibnd, na, ngg, ikw 
  integer, allocatable :: kisort(:) 
  integer :: npool, nkbl, nkl, nkr, npwx_g 
  integer :: ike, iks, npw_g, ispin 
  integer, allocatable :: ngk_g( : ) 
  integer, allocatable :: ngk_gw( : ) 
  integer, allocatable :: itmp( :, : ) 
  integer, allocatable :: igwk( : ) 
  integer, allocatable :: igk_l2g( :, : ) 
 
  real(kind=8) :: wfc_scal  
  logical :: twf0, twfm, twrite_wfc 
 
  ! 
  ! 
  IF( ionode ) WRITE( stdout, fmt="(//,'WRITING G-MAPS for each processor' )" ) 
 
  IF( nkstot > 0 ) THEN 
 
     IF( ( kunit < 1 ) .OR. ( MOD( nkstot, kunit ) /= 0 ) ) & 
       CALL errore( ' write_wannier ',' wrong kunit ', 1 ) 
 
     IF( ( nproc_pool > nproc ) .OR. ( MOD( nproc, nproc_pool ) /= 0 ) ) & 
       CALL errore( ' write_wannier ',' nproc_pool ', 1 ) 
 
     !  find out the number of pools 
     npool = nproc / nproc_pool 
 
     !  find out number of k points blocks 
     nkbl = nkstot / kunit 
 
     !  k points per pool 
     nkl = kunit * ( nkbl / npool ) 
 
     !  find out the reminder 
     nkr = ( nkstot - nkl * npool ) / kunit 
 
     !  Assign the reminder to the first nkr pools 
     IF( my_pool_id < nkr ) nkl = nkl + kunit 
 
     !  find out the index of the first k point in this pool 
     iks = nkl * my_pool_id + 1 
     IF( my_pool_id >= nkr ) iks = iks + nkr * kunit 
 
     !  find out the index of the last k point in this pool 
     ike = iks + nkl - 1 
 
  END IF 
 
  ! find out the global number of G vectors: ngm_g   
  ngm_g = ngm 
  call mp_sum( ngm_g, intra_pool_comm ) 
 
 
  ! build the G+k array indexes 
  allocate ( kisort( npwx ) ) 
  allocate ( igk_l2g( npwx, ik ) ) 
  do ik = 1, nks 
     kisort = 0 
     call gk_sort (xk (1, ik+iks-1), ngm, g, ecutwfc / tpiba2, npw, kisort(1), g2kin) 
     DO ig = 1, npw
        igk_l2g(ig,ik) = ig_l2g(kisort(ig))
     END DO
     ngk (ik) = npw 
  end do 
  deallocate (kisort) 
 
  ! compute the global number of G+k vectors for each k point 
  allocate( ngk_g( nkstot ) ) 
  allocate( ngk_gw( nkstot/nspin ) ) 
  ngk_g = 0 
  ngk_g( iks:ike ) = ngk( 1:nks ) 
  CALL mp_sum( ngk_g ) 
 
  ! compute the Maximum G vector index among all G+k an processors 
  npw_g = MAXVAL( ig_l2g(:) ) ! ( igk_l2g(:,:) ) 
  CALL mp_max( npw_g ) 
 
  ! compute the Maximum number of G vector among all k points 
  npwx_g = MAXVAL( ngk_g( 1:nkstot ) ) 
 
 
  allocate( igwk( npwx_g ) ) 
 
  do ik = 1, nkstot 
    igwk = 0 
    allocate( itmp( npw_g, 1 ) ) 
    itmp = 0 
    if( ik >= iks .AND. ik <= ike ) then  
      do  ig = 1, ngk( ik-iks+1 ) 
        itmp( ig_l2g( ig ), 1 ) = ig_l2g( ig )
      end do 
    end if 
    call mp_sum( itmp ) 
    ngg = 0 
    do  ig = 1, npw_g 
      if( itmp( ig, 1 ) == ig ) then 
        ngg = ngg + 1 
        igwk( ngg ) = ig 
      end if 
    end do 
    if( ngg /= ngk_g( ik ) ) then 
      WRITE( stdout,*) ' ik, ngg, ngk_g = ', ik, ngg, ngk_g( ik ) 
    end if 
    deallocate( itmp ) 
    if( ionode ) then 
        ! write (40)( igwk(ig), ig = 1, npwx_g ) 
    end if 
  end do 
 
  deallocate( igwk ) 
 
  do ik = 1, nkstot 
     IF( (ik >= iks) .AND. (ik <= ike) ) THEN 
        ispin = isk( ik ) 
        WRITE( 100 + mpime ) ik, iks, ike, nkstot, kunit, nproc, ispin, nspin, npw_g, & 
                             nbnd, ngk(ik-iks+1), nwordwfc, npwx, iunwfc, nd_nmbr 
        WRITE( 100 + mpime ) ( igk_l2g( i, ik-iks+1 ), i = 1, ngk(ik-iks+1) ) 
     END IF 
  end do 
 
  deallocate ( ngk_g ) 
  deallocate ( ngk_gw ) 
  deallocate (igk_l2g) 
 
end subroutine write_gmaps 
 
 
subroutine read_and_collect( c, ldc, n, ik ) 
  USE io_global,      ONLY : stdout 
  USE io_files,       ONLY : prefix 
  USE kinds,          ONLY : DP, sgl 
 
  implicit none 
 
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
    IF( ( ik_ == ik ) .AND. ( ik_ >= iks ) .AND. ( ik_ <= ike ) ) THEN 
      ALLOCATE( evc( npwx, nbnd ) ) 
      ALLOCATE( igk_l2g( ngk ) ) 
      READ( 100 + ip )  ( igk_l2g( i ), i = 1, ngk ) 
      CALL diropn_gw ( 99, TRIM( prefix )//'.wfc', nwordwfc, exst, ip, nd_nmbr ) 
      call davcio ( evc, nwordwfc, 99, (ik-iks+1), - 1 ) 
      CLOSE( 99 ) 
      DO j = 1, n 
        DO i = 1, ngk 
          c( igk_l2g( i ), j ) = evc( i, j ) 
        END DO 
      END DO 
      DEALLOCATE( evc ) 
      DEALLOCATE( igk_l2g ) 
    END IF 
    REWIND( 100 + ip ) 
  END DO 
 
  return 
end subroutine 
 
! 
! Copyright (C) 2001-2003 PWSCF group 
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 
! 
! 
!----------------------------------------------------------------------- 
subroutine diropn_gw (unit, filename, recl, exst, mpime, nd_nmbr_ ) 
  !----------------------------------------------------------------------- 
  ! 
  !     this routine opens a file in tmp_dir for direct I/O access 
  !     If appropriate, the node number is added to the file name 
  ! 
  USE kinds 
  use io_files 
  implicit none 
 
#if defined(__ALPHA)
#  define DIRECT_IO_FACTOR 2
#else
#  define DIRECT_IO_FACTOR 8
#endif

  ! 
  !    first the input variables 
  ! 
  character(len=*) :: filename 
  ! input: name of the file to open 
  integer :: unit, recl 
  ! input: unit of the file to open 
  ! input: length of the records 
  logical :: exst 
  ! output: if true the file exists 
  integer :: mpime 
  ! input: processor index 
  CHARACTER(LEN=3) :: nd_nmbr_ 
  ! 
  !    local variables 
  ! 
  character(len=256) :: tempfile 
  ! complete file name 
  character(len=80) :: assstr 
  integer :: ios, unf_recl, ierr 
  ! used to check I/O operations 
  ! length of the record 
  ! error code 
  logical :: opnd 
  ! if true the file is already opened 
 
 
  if (unit < 0) call errore ('diropn', 'wrong unit', 1) 
  ! 
  !    we first check that the file is not already openend 
  ! 
  ios = 0 
  inquire (unit = unit, opened = opnd) 
  if (opnd) call errore ('diropn', "can't open a connected unit", abs(unit)) 
  ! 
  !      then we check the filename 
  ! 
 
  if (filename == ' ') call errore ('diropn', 'filename not given', 2) 
  tempfile = trim(tmp_dir) // trim(filename) // trim( nd_nmbr_ ) 
 
  inquire (file = tempfile, exist = exst) 
  ! 
  !      the unit for record length is unfortunately machine-dependent 
  ! 
  unf_recl = DIRECT_IO_FACTOR * recl 
  if (unf_recl <= 0) call errore ('diropn', 'wrong record length', 3) 
  ! 
  !     on T3E reduce the size of the buffer if it is too large 
  ! 
#ifdef __T3E 
  if (unf_recl.gt.5000000) then 
     if (unit < 10) then 
        write (assstr, '("assign -b 1 u:",i1)') unit 
     else if(unit < 100) then 
        write (assstr, '("assign -b 1 u:",i2)') unit 
     else 
        call errore ('diropn', 'unit too large', 1) 
     endif 
     call assign (assstr, ierr) 
  endif 
#endif 
 
  open ( unit, file = TRIM(tempfile), iostat = ios, form = 'unformatted', & 
       status = 'unknown', access = 'direct', recl = unf_recl ) 
 
  if (ios /= 0) call errore ('diropn', 'error opening '//filename, unit) 
  return 
end subroutine diropn_gw 
