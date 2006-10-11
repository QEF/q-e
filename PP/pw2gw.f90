
! Copyright (C) 2004 PWSCF group 
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
  !
  IMPLICIT NONE
  INTEGER :: ios

  NAMELIST / inputpp / prefix, outdir

  CALL start_postproc(nd_nmbr)
  ! 
  !   set default values for variables in namelist 
  ! 
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'

  IF ( ionode )  THEN 
     !
     READ (5, inputpp, err=200, iostat=ios)
200  CALL errore('pw2gw', 'reading inputpp namelist', ABS(ios))
  tmp_dir = trimcheck (outdir)
     !
  END IF
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( prefix, ionode_id ) 
  CALL mp_bcast(tmp_dir, ionode_id ) 
  !
  CALL read_file
  CALL openfil_pp
  !
  CALL compute_gw
  !
  CALL stop_pp
  STOP

END PROGRAM pw2gw


SUBROUTINE compute_gw

  ! This routine creates the QPLDA and the matrixelements
  ! tform = .false. UNFORMATTED QPLDA
  ! tform = .true.  FORMATTED QPLDA
  ! tsingle must be always true 

  USE kinds, ONLY: DP, sgl
  USE constants, ONLY : eps8
  USE cell_base, ONLY: alat, tpiba2, at, bg
  USE char, ONLY: title
  USE symme, ONLY: s, nsym
  USE wvfct, ONLY: npw, npwx, nbnd, gamma_only, igk, g2kin, wg, et
  USE gvect, ONLY: ngm, g, gg, ecutwfc
  USE klist , ONLY: nks, xk, wk
  USE lsda_mod, ONLY: nspin
  USE io_files, ONLY: nwordwfc, iunwfc
  USE wavefunctions_module, ONLY : evc

  IMPLICIT NONE
  INTEGER :: ii(16), ngw, nkpt, ig, ik, n, i,j,k, io = 98, iband1, iband2
  INTEGER, ALLOCATABLE :: in1(:), in2(:), in3(:)
  LOGICAL :: t_form = .false., t_single = .true.
  REAL(kind=sgl) :: a1_s(3), a2_s(3), a3_s(3)
  REAL(kind=sgl), ALLOCATABLE :: xk_s(:,:), eig_s(:,:), focc_s(:,:)
  REAL(kind=DP):: g2max, a1(3), a2(3), a3(3),norm, xkgk(3), rrhotwx(3), delta
  REAL(kind=DP), ALLOCATABLE:: gsort(:), eig(:,:), focc(:,:), kpg(:,:)
  COMPLEX(kind=DP):: rhotwx(3), ctemp
  COMPLEX(kind=DP),  ALLOCATABLE:: c0(:), c0_m(:,:)
  COMPLEX(kind=sgl), ALLOCATABLE:: c0_s(:)
  CHARACTER(LEN=80) :: titleo(2) 
  INTEGER :: igwx
  !
  !
  IF( nspin > 1 ) CALL errore('pw2gw','Spin polarization not implemented',1)
  !
  !!! CALL init_us_1
  !!! CALL newd
  !
  IF (t_form) THEN
     WRITE (6,'(//" writing LDA info on unit 98 FORMATTED")')
     OPEN (io, FILE='QPLDA',STATUS='new',FORM='FORMATTED')
  ELSE
     WRITE (6,'(//" writing LDA info on unit io UNFORMATTED")')
     OPEN (io, FILE='QPLDA',STATUS='new',FORM='UNFORMATTED')
  ENDIF
     WRITE (6,'(//" writing matrixelements on unit 98 FORMATTED")')
     OPEN (90, FILE='matrixelements',STATUS='new',FORM='FORMATTED')
  !
  !  file's title [2 lines]
  !
  titleo(1)='pw2gw'
  titleo(2)='test version'
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
     write (io,'(16I5)') ii
  ELSE
     ii(1)=1
     write (io) ii
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
  IF (t_form) THEN
     WRITE (io,'(3E26.18)') a1, a2, a3
  ELSE
  IF (t_single) THEN
     WRITE (io) a1_s, a2_s, a3_s
  ELSE
     WRITE (io) a1, a2, a3
  ENDIF
  ENDIF

  WRITE(6,*) 'Vettori di reticolo diretto'
  write(6,'(a,3E26.18)') 'a1', a1_s
  write(6,'(a,3E26.18)') 'a2', a2_s
  write(6,'(a,3E26.18)') 'a3', a3_s
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
  write(6,*)'nrot=',nsym
  write(6,'(3E26.18)') (((REAL(s(i,j,k), kind=sgl),j=1,3),i=1,3),k=1,nsym) 
  IF (t_form) THEN
     WRITE (io,'(I2)') nsym
     WRITE (io,'(3E26.18)') (((REAL(s(i,j,k), KIND=sgl),j=1,3),i=1,3),k=1,nsym)
   IF (ii(3) == 1) THEN 
        ! READ (10,1020) ((VOFFSET(I,J),I=1,3),J=1,NOP)
        ! WRITE (6,'(//" Run program CNVNSY to convert QPLDA file first.")')
        CALL errore('pw2gw','non-symmorphic translation vectors',ii(3))
   ENDIF
  ELSE
     WRITE (io) nsym
   IF (t_single) THEN
      WRITE (io) (((REAL(s(i,j,k), KIND=sgl),j=1,3),i=1,3),k=1,nsym)
   ELSE
      WRITE (io) (((REAL(s(i,j,k), KIND=DP ), j=1,3),i=1,3),k=1,nsym) 
   ENDIF
   IF (ii(3) == 1) THEN
       ! READ (10,1020) ((VOFFSET(I,J),I=1,3),J=1,NOP)
       CALL errore('pw2gw','non-symmorphic translation vectors',ii(3))
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
  WRITE(6,*) ' nks ', nks
  WRITE(6,*) ' k points in  cartesian coordinates'
  WRITE(6,'(1x,3f10.6)') ( (xk(i,ik),i=1,3), ik=1,nks)
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
  !
  !  ngw =  number of G-vectors (complete shells) such that G2 <= G2max
  !  ovvero <= RAGGIO della SFERA, in pratica trova i G2 relativi a GAMMA
  !
  do ngw = 1, ngm
     if ( gg(ngw) > g2max + eps8) go to 100
  end do
  call errore ( 'pw2gw','max G in PW not found?!?',ngw)
100 ngw = ngw - 1

!  write( 6, * ) ' ngw = ', ngw
!  write( 6, * ) '  modules   = ', gg(ngw), gg(ngw+1), g2max

  ! Pongo NGW pari al massimo indice tra i vettori G che fanno parte delle 
  ! sfere |G+k|<cut per qualsiasi k
  !

  write( 6, * ) ' igwx= ', igwx

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
  IF (t_form) THEN
     WRITE (io,'(I12)') ngw
     WRITE (io,'(3I5)') (in1(ig),in2(ig),in3(ig),ig=1,ngw)
  ELSE
     WRITE (io) ngw
     WRITE (io) (in1(ig),in2(ig),in3(ig),ig=1,ngw)
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
  OPEN(65,file='k.dat')
  WRITE(65,'(1x,3f10.6,1x,f10.6)')  ( xk_s(:,ik), wk(ik)*0.5, ik=1,nks )
  CLOSE(unit=65)
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
!  DEALLOCATE (xk_s)
  !
  !  WRITE energies (Hartrees) (in ascending order, hopefully they are ordered)
  !
  n = nbnd
  ALLOCATE (eig(n,nkpt), eig_s(n,nkpt))
  eig(:,:)   = et(:,1:nkpt)*0.5
!  eig_s(:,:) = et(:,1:nkpt)*0.5
  IF (t_form) THEN
     WRITE (io,'(I12)') n
     WRITE (io,'(3E26.18)') ((eig(i,ik),ik=1,nkpt),i=1,n)
  ELSE
     WRITE (io) n
     if(t_single) then
        do ik=1,nkpt
           do i=1,n
              eig_s(i,ik)=eig(i,ik)
           enddo
        enddo
        WRITE(6,*) 'nbndsi=',n
        WRITE (io) ((eig_s(i,ik),ik=1,nkpt),i=1,n)
     else
        WRITE(6,*) 'nbndsi=',n
        WRITE (io) ((eig(i,ik),ik=1,nkpt),i=1,n)
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
    focc(j,ik)=wg(j,ik)*2.0/wk(ik)
    enddo
  enddo
  
  focc_s(:,:) = focc(:,:)

  IF (t_form) THEN
     WRITE (io,'(3E26.18)') ((focc(i,ik), ik=1,nkpt), i=1,n)
  ELSE
     if(t_single) then
        WRITE (io) ((focc_s(i,ik),ik=1,nkpt),i=1,n)
     else
        WRITE (io) ((focc(i,ik),ik=1,nkpt),i=1,n)
     endif
  ENDIF
  WRITE (6,*) nkpt
  WRITE (6,*) 'weights:'
  WRITE (6,'(10f10.7)') (wk(ik), ik=1,nkpt)

  do ik=1,nkpt
   WRITE (6,*) 'ik=', ik
   WRITE (6,'(10f10.7)') (focc_s(j,ik), j=1,n)
  enddo

!  DEALLOCATE (focc_s, focc)
  !
  ! wavefunctions
  !  Read wavefunctions (a(G), normalised to 1)
  !
  ALLOCATE ( c0(igwx), c0_s(igwx), kpg(3,igwx),c0_m(igwx,n) )
  CALL cryst_to_cart (nks, xk, bg, +1)
  WRITE(6,*) 'Costruisco le psi ed il matrixelements'
  WRITE(6,*) 'Controllo: I punti k ora devo essere in coordinate cartesiane!'
  WRITE(6,'(1x,3f10.6)') ( (xk(i,ik),i=1,3), ik=1,nkpt)

  DO ik = 1, nkpt     
    CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
    CALL davcio( evc, nwordwfc, iunwfc, ik, -1 )
    !
    ! copy coefficient from array evc(:,n) (ordered as |k+G|) 
    ! into array c0 with |G| ordering
    !
!DEBUG     WRITE(6, * )  'npw, n    = ', npw, n
!DEBUG     WRITE(6, * )  'SIZE(evc) = ', SIZE(evc,1), SIZE(evc,2)
     
     DO ig=1,npw
       IF( igk(ig) < 1 .OR. igk(ig) > SIZE( c0 ) ) &
         CALL errore(' pw2gw ', ' c0 too small ', 1 )
     END DO

    ! read wavefunctions and write the matrixelemnts

     DO i = 1, n
       ! important: missing components must be set to zero
       c0 (:) = 0.d0
       DO ig=1,npw
         c0(igk(ig)) = evc(ig,i)
       END DO
       c0_m(:,i)=c0(:)  
 
       WRITE(6,*) c0_m(1,1)
       WRITE(6,*) c0_m(3,1)
      
       IF (t_form) THEN
         WRITE (io,'(3E26.18)') c0
        ELSE
        IF (t_single) THEN
          c0_s(:) = c0(:)
          WRITE(io) c0_s
        ELSE
          WRITE(io) c0
        ENDIF
      ENDIF
     END DO
     WRITE(6, * )  ' out from bands loop of wfc'
    
    ! k + g thet must be in 2piba units

    kpg(:,:) = 0.d0
    DO ig=1,npw
      kpg(:,igk(ig))= xk_s(:,ik)+g(:,igk(ig))
    END DO

! START MATRIXELEMENTS

  DO iband1 = 1,n
  IF ( focc(iband1,ik).ge.1e-4) THEN
  DO iband2 = 1,n
  delta=2.0d0-focc(iband2,ik)
  IF (delta.gt.1e-4) THEN
                rhotwx = 0.0
                DO  ig=1,igwx
                    xkgk(1)= kpg(1,ig)
                    xkgk(2)= kpg(2,ig)
                    xkgk(3)= kpg(3,ig)
                    ctemp= CONJG(c0_m(ig,iband1))*c0_m(ig,iband2) 
                    rhotwx(1) = rhotwx(1) + xkgk(1) * ctemp
                    rhotwx(2) = rhotwx(2) + xkgk(2) * ctemp
                    rhotwx(3) = rhotwx(3) + xkgk(3) * ctemp
                ENDDO ! STOP END igwx
                rrhotwx(1)=tpiba2* real(rhotwx(1)*conjg(rhotwx(1)))
                rrhotwx(2)=tpiba2* real(rhotwx(2)*conjg(rhotwx(2)))
                rrhotwx(3)=tpiba2* real(rhotwx(3)*conjg(rhotwx(3)))
                WRITE (90,'(1x,3i5,3e16.8,2f8.4)') ik,iband1,iband2,rrhotwx(1),rrhotwx(2), &
                rrhotwx(3),(eig(iband2,ik)-eig(iband1,ik))*27.21, &
                       (focc(iband1,ik)-focc(iband2,ik))
           ENDIF 
       ENDDO
   ENDIF
  ENDDO
!end matrixelements
!---------------------------------------------------------------------

 ENDDO
     WRITE(6, * )  ' out from k-points loop'

  DEALLOCATE (xk_s)
  DEALLOCATE (eig_s, eig)
  DEALLOCATE (focc_s, focc)
  DEALLOCATE (c0_s, c0, kpg, c0_m)
  !
  CLOSE(90)
  CLOSE(io)
  !
END SUBROUTINE compute_gw


