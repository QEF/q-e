!-----------------------------------------------------------------------
program fd_raman
!-----------------------------------------------------------------------
  use constants
  use io_files,   ONLY : prefix, tmp_dir, psfile, pseudo_dir
  use io_global,  ONLY : stdout, ionode, ionode_id
  USE mp_global,  ONLY : mp_startup
  USE environment,ONLY : environment_start
  USE mp,         ONLY : mp_bcast
  USE cell_base,  ONLY : tpiba2, alat,omega, at, bg, ibrav, celldm
  USE ions_base,  ONLY : amass, nat, atm, zv, tau, ntyp => nsp, ityp
  USE kinds,      ONLY : dp 
  USE gvecw,      ONLY : ecutwfc
  USE symm_base,       ONLY : nsym, nsym_ns, nsym_na, invsym, s, sr, &
                              t_rev, ftau, sname
  USE symme
  USE rap_point_group, ONLY : code_group, nclass, nelem, elem, &
       which_irr, char_mat, name_rap, name_class, gname, ir_ram
  USE rap_point_group_so, ONLY : nrap, nelem_so, elem_so, has_e, &
       which_irr_so, char_mat_so, name_rap_so, name_class_so, d_spin, &
       name_class_so1
  USE rap_point_group_is, ONLY : nsym_is, sr_is, ftau_is, d_spin_is, &
       gname_is, sname_is, code_group_is
  USE fft_base, ONLY : dfftp

  USE parser,    ONLY : field_count, read_line, get_field, parse_unit

  implicit none
  character(len=9) :: code = 'FD_RAMAN'
  integer :: ios
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  character(len=200) :: pp_file
  logical :: uspp_spsi, ascii, single_file, raw, disp_only

  INTEGER :: apol, na, nt
  integer           :: nrx1,nrx2,nrx3,nr1,nr2,nr3,nb,nax,natx,inn
  real(kind=dp)     :: r1(3),r2(3),r3(3),rr(3,3)
  INTEGER :: nclass_ref   ! The number of classes of the point group
  INTEGER :: isym, ipol
  REAL (dp) :: ft1, ft2, ft3

  integer :: npol, npol_rm, npol1_rm, npol_eps, npol_zeu
  integer :: ndiag,noffd,nmod,npol1
  integer :: i,j,p,k,ii,jj,n
  real*8, allocatable :: F0(:,:),dechi(:,:,:,:),Fd(:,:,:,:),dechi_u(:,:,:,:)
  real*8, allocatable :: Fij(:,:,:,:),alpha(:,:,:),u(:,:,:),ui(:,:,:)
  real*8 :: de, de_raman, de_zeu, de_eps,conv,ry,a0
  CHARACTER(50)    :: filemodes
  logical :: lalpha,lpuma

  real*8, allocatable :: pol0(:),pol(:,:,:),eps(:,:)
  real*8, allocatable :: zeta(:,:,:)
  real*8 :: sum

  CHARACTER(len=2)           :: prog   ! calling program ( PW, CP, WA )
  CHARACTER(len=256)         :: input_line
  CHARACTER(len=80)          :: card
  CHARACTER(len=1), EXTERNAL :: capital
  LOGICAL                    :: tend, verbose

  NAMELIST /inputfd/ prefix,npol_rm, npol_eps, npol_zeu, ndiag,noffd,nmod,npol1,lpuma, &
                     de_raman, de_eps, de_zeu, filemodes, verbose

  lalpha = .true.
  lpuma = .false.
  npol_rm=4
  npol1 = 2
  npol_eps=2
  npol_zeu=2
  ndiag=3
  noffd=3
  de_raman=1.0
  de_eps=1.0
  de_zeu=1.0
  filemodes=' '
  verbose=.false.

  ! define conversion constant
  conv=BOHR_RADIUS_ANGS**2

  CALL mp_startup ( )
  CALL environment_start ( code )

  IF ( ionode ) THEN
    CALL input_from_file ( )
    READ(5,inputfd,IOSTAT=ios)
  endif
  if (filemodes .eq. ' ') lalpha=.false.

  !reading the xml file
  call read_file

  if (ionode) then
    write(6,*) '**************************************************'
    write(6,*) '* Info from the preceding pw.x run:              *'
    write(6,*) '**************************************************'
    write(6,*) ''
    write(6,*) '    prefix=  ',trim(prefix)
    write(6,*) '    outdir=  ',trim(tmp_dir)
    write(6,*) '    ecutwfc= ',ecutwfc, 'Ry'

    WRITE( stdout, 199) ibrav, alat, omega, nat, ntyp
    199 FORMAT(5X, &
      &     'bravais-lattice index     = ',I12,/,5X, &
      &     'lattice parameter (alat)  = ',F12.4,'  a.u.',/,5X, &
      &     'unit-cell volume          = ',F12.4,' (a.u.)^3',/,5X, &
      &     'number of atoms/cell      = ',I12,/,5X, &
      &     'number of atomic types    = ',I12)
    !
    WRITE( stdout, '(/2(3X,3(2X,"celldm(",I1,")=",F11.6),/))' ) &
      ( i, celldm(i), i = 1, 6 )
    !

    !lattice vectors in Angs
    rr = at*alat*bohr_radius_angs
    r1(:) = rr(:,1)
    r2(:) = rr(:,2)
    r3(:) = rr(:,3)
    WRITE( stdout, '(5X, &
      &     "Lattice vectors: (cart. coord. in Angs)",/, &
      &       3(15x,"a(",i1,") = (",3f11.6," )  ",/ ) )')  (apol,  &
      (rr (ipol, apol) , ipol = 1, 3) , apol = 1, 3)

    !atoms
    WRITE( stdout, '(/,5x,"Atomic coordiantes")')                                                  
    WRITE( stdout, '(5x,"site n.     atom                  positions (Angs)")')
    WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")') &
      (na, atm(ityp(na)), na, (tau(ipol,na)*alat*0.5291772, ipol=1,3), na=1,nat)

    !symmetries
    IF (verbose) THEN
    write(6,*) 
     
     WRITE( stdout, '(36x,"s",24x,"frac. trans.")')
     nsym_is=0
     DO isym = 1, nsym
        WRITE( stdout, '(/6x,"isym = ",i2,5x,a45/)') isym, sname(isym)
        IF ( ftau(1,isym).NE.0 .OR. ftau(2,isym).NE.0 .OR. &
             ftau(3,isym).NE.0) THEN
           ft1 = at(1,1)*ftau(1,isym)/dfftp%nr1 + at(1,2)*ftau(2,isym)/dfftp%nr2 + &
                at(1,3)*ftau(3,isym)/dfftp%nr3
           ft2 = at(2,1)*ftau(1,isym)/dfftp%nr1 + at(2,2)*ftau(2,isym)/dfftp%nr2 + &
                at(2,3)*ftau(3,isym)/dfftp%nr3
           ft3 = at(3,1)*ftau(1,isym)/dfftp%nr1 + at(3,2)*ftau(2,isym)/dfftp%nr2 + &
                at(3,3)*ftau(3,isym)/dfftp%nr3
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                &        " )    f =( ",f10.7," )")') &
                isym, (s(1,ipol,isym),ipol=1,3), DBLE(ftau(1,isym))/DBLE(dfftp%nr1)
           WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )")') &
                (s(2,ipol,isym),ipol=1,3), DBLE(ftau(2,isym))/DBLE(dfftp%nr2)
           WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )"/)') &
                (s(3,ipol,isym),ipol=1,3), DBLE(ftau(3,isym))/DBLE(dfftp%nr3)
           WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7, &
                &        " )    f =( ",f10.7," )")') &
                isym, (sr(1,ipol,isym),ipol=1,3), ft1
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )")') &
                (sr(2,ipol,isym),ipol=1,3), ft2
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )"/)') &
                (sr(3,ipol,isym),ipol=1,3), ft3
        ELSE
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), " )")') &
                isym,  (s (1, ipol, isym) , ipol = 1,3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )")')  (s(2,ipol,isym), ipol=1,3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )"/)') (s(3,ipol,isym), ipol=1,3)
           WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7," )")') &
                isym,  (sr (1, ipol,isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )")')  (sr (2, ipol,isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )"/)') (sr (3, ipol,isym) , ipol = 1, 3)
        END IF
     END DO
    END IF
  end if

100   CALL read_line( input_line, end_of_file=tend )
      !
      IF( tend ) GOTO 120
      IF( input_line == ' ' .OR. input_line(1:1) == '#' .OR. &
                                 input_line(1:1) == '!' ) GOTO 100
      !
      READ (input_line, *) card
      !
      DO i = 1, len_trim( input_line )
         input_line( i : i ) = capital( input_line( i : i ) )
      ENDDO

  IF ( trim(card) == 'RAMAN_TENSOR') THEN

! read forces from input card "RAMAN_TENSOR"
! Arrigo Calzolari's convention (to be automated)
!
! npol (lpuma) = 1,2,3,4, --> -2h, -h, +h, +2h
! npol (else)  = 1,2,3,4, --> -1-1, +1+1, +1,-1, -1+1
! npol1 (else) = 1,2 --> -h, h
! ndiag = 1,2,3 --> Ex, Ey, Ez
! noffd = 1,2,3 --> Exy, Exz, Eyz

  npol=npol_rm
  allocate(F0(3,nat))
  allocate(dechi(3,3,3,nat))
  if (lpuma) then
     allocate(Fd(npol,ndiag,3,nat))
  else
     allocate(Fd(npol1,ndiag,3,nat))
  end if
  allocate(Fij(npol,noffd,3,nat))

  F0(:,:)=0.0d0
  dechi(:,:,:,:)=0.0d0
  Fd(:,:,:,:)=0.0d0
  Fij(:,:,:,:)=0.0d0

  ! read data from input

  do i=1,nat
   read(5,*) (F0(k,i), k=1,3)
  end do

  do ii=1,ndiag
   do p=1,npol1
    do i=1,nat
     read(5,*) (Fd(p,ii,k,i), k=1,3)
    end do
   enddo
  end do
  
  do ii=1,noffd
   do p=1,npol
    do i=1,nat
     read(5,*) (Fij(p,ii,k,i), k=1,3)
    end do
   enddo
  end do
  
  dechi(:,:,:,:)=0.0d0
  de=de_raman
  do i=1,nat
   do ii=1,3
     do k=1,3
     if (lpuma)  then
       dechi(ii,ii,k,i)=(-1.0*Fd(1,ii,k,i)+16.0*Fd(2,ii,k,i)-30.0*F0(i,k)+16.0*Fd(3,ii,k,i)   &
                         -1.0*Fd(4,ii,k,i))/(12.*de**2)
     else
       dechi(ii,ii,k,i)=(Fd(1,ii,k,i)-2*F0(i,k)+Fd(2,ii,k,i))/(de**2)
     end if
     end do
   end do
  end do

! construct d chi/dE1dE2
  IF (verbose) THEN
  write(6,*) '**************************************************'
  WRITE( stdout, '("unsymmetrized Raman tensor")')
  write(6,*) '**************************************************'
  END IF

  do i=1,nat
    do k=1,3
     if (lpuma)  then
      dechi(1,2,k,i) = (-1.0*Fij(1,1,k,i)+16.0*Fij(2,1,k,i)-30.0*F0(i,k)+16.0*Fij(3,1,k,i)    &
                        -1.0*Fij(4,1,k,i))/(12.0*de**2)
      dechi(1,2,k,i) = 0.5*dechi(1,2,k,i)-0.5*dechi(1,1,k,i)-0.5*dechi(2,2,k,i)
      dechi(2,1,k,i) = dechi(1,2,k,i)
  
      dechi(1,3,k,i) = (-1.0*Fij(1,2,k,i)+16.0*Fij(2,2,k,i)-30.0*F0(i,k)+16.0*Fij(3,2,k,i)    &
                        -1.0*Fij(4,2,k,i))/(12.0*de**2)
      dechi(1,3,k,i) = 0.5*dechi(1,3,k,i)-0.5*dechi(1,1,k,i)-0.5*dechi(3,3,k,i)
      dechi(3,1,k,i) = dechi(1,3,k,i)
  
      dechi(2,3,k,i) = (-1.0*Fij(1,3,k,i)+16.0*Fij(2,3,k,i)-30.0*F0(i,k)+ 16.0*Fij(3,3,k,i)   &
                        -1.0*Fij(4,3,k,i))/(12.0*de**2)
      dechi(2,3,k,i) = 0.5*dechi(2,3,k,i)-0.5*dechi(2,2,k,i)-0.5*dechi(3,3,k,i)
      dechi(3,2,k,i) = dechi(2,3,k,i)
     else
      dechi(1,2,k,i) = (Fij(1,1,k,i)+Fij(2,1,k,i)-Fij(3,1,k,i)-Fij(4,1,k,i))/(4*de**2)
      dechi(2,1,k,i) = dechi(1,2,k,i)
  
      dechi(1,3,k,i) = (Fij(1,2,k,i)+Fij(2,2,k,i)-Fij(3,2,k,i)-Fij(4,2,k,i))/(4*de**2)
      dechi(3,1,k,i) = dechi(1,3,k,i)
  
      dechi(2,3,k,i) = (Fij(1,3,k,i)+Fij(2,3,k,i)-Fij(3,3,k,i)-Fij(4,3,k,i))/(4*de**2)
      dechi(3,2,k,i) = dechi(2,3,k,i)
     end if
    end do
  end do

  do i=1,nat
    do ii=1,3
      do jj=1,3
        do k=1,3
          dechi(ii,jj,k,i)=dechi(ii,jj,k,i)/omega ! *(-1.0d0)
        end do
      end do
    end do
  end do
  
    IF (verbose) THEN
  do i=1,nat
    do k = 1, 3
      write (6,'(5x,"atom # ",i4,"    pol.",i3)') i,k
      do ii =1,3
        write(6,43) (dechi(ii,jj,k,i)*omega, jj=1,3)
      end do
    end do
    write(6,*)'  '
  end do
    END IF
  write(6,*) '**************************************************'
  WRITE( stdout, '("Raman tensor (A^2)")')
  write(6,*) '**************************************************'

  ! convert in crystal coordinates

  do na = 1,nat
     call cart_to_crys_mat3 ( dechi(1,1,1,na) )
  end do

  call symtensor3( nat, dechi)

  do i=1,nat
    do ii=1,3
      do jj=1,3
        do k=1,3
          dechi(ii,jj,k,i)=conv*omega*dechi(ii,jj,k,i)
        end do
      end do
    end do
  end do
  
  do i=1,nat
    do k = 1, 3
      write (6,'(5x,"atom # ",i4,"    pol.",i3)') i,k
      do ii =1,3
        write(6,34) (dechi(ii,jj,k,i), jj=1,3)
      end do
    end do
    write(6,*)'  '
  end do

if (lalpha)  then
  write(6,*) '**************************************************'
  WRITE( stdout, '("Raman alpha tensor")')
  write(6,*) '**************************************************'

  nmod=3*nat
  allocate (u(nmod,3,nat))
  allocate (ui(nmod,3,nat))
  allocate(alpha(3,3,nmod))
  allocate(dechi_u(3,3,3,nat))

  alpha(:,:,:)=0.0d0
  u(:,:,:)=0.0d0

! read normalized eigenmodes from matdyn.modes

  open (2,file=TRIM(filemodes),form='formatted')
  
  read(2,*)
  read(2,*)
  read(2,*)
  read(2,*)

  do n=1,nmod
    read (2,*)
    do i=1,nat
      read(2,'(1x,1x,3 (f10.6,1x,f10.6,3x),1x)') (u(n,k,i),ui(n,k,i),k=1,3)
      do k=1,3
        u(n,i,k)=u(n,k,i)/Sqrt(amass(ityp(i)))
      end do
    end do
  end do
  close(2)

  do n=1,nmod
   do ii=1,3
     do jj=1,3
       do i=1,nat
         do k=1,3
           dechi_u(ii,jj,k,i)=dechi(ii,jj,k,i)*u(n,k,i)
           alpha(ii,jj,n)=alpha(ii,jj,n)+dechi_u(ii,jj,k,i)
         end do
       end do
       alpha(ii,jj,n)=Sqrt(omega)*alpha(ii,jj,n)
     end do
   end do
  end do

  write(6,*)''
  do n=1,nmod
  write(6,*) n
   do ii=1,3
     write(6,43) (alpha(ii,jj,n), jj=1,3)
   end do
  end do
  deallocate(u,alpha,dechi_u)
 end if


32 format(a,i5,a,i5)
43 format(3f12.6)
34 format(3e24.12)
deallocate(F0,dechi,Fd,Fij)

    ELSEIF ( trim(card) == 'DIELECTRIC_TENSOR') THEN

  write(6,*) '**************************************************'
  WRITE( stdout, '("Dielectric tensor")')
  write(6,*) '**************************************************'

   npol=npol_eps

   allocate (pol0(3))
   allocate (pol(npol,3,3))
   allocate (eps(3,3))
!
   pol0=0.0d0
   de=de_eps*omega
   read(5,*) (pol0(i),i=1,3)
   do i=1,3
     if (Abs (pol0(i)) .lt. conv) pol0(i)=0.0d0
   end do
!
   pol=0.0d0
   do p=1,npol
     do i=1,3
       read(5,*) (pol(p,i,j), j=1,3)
     end do
   end do
!
   do i=1,3
     do j=1,3
       if (npol==2) then
         eps(i,j)=(pol(1,i,j)-pol0(j))-(pol(2,i,j)-pol0(j))
         eps(i,j)=2*pi*eps(i,j)/de
       else
         eps(i,j)=(pol(1,i,j)-pol0(j))
         eps(i,j)=4*pi*eps(i,j)/de
       end if
       if (i==j) eps(i,j)=eps(i,j)+1.0
     end do
   end do
!
   call symmatrix ( eps )

   do i=1,3
     write(6,33) (eps(i,j),j=1,3)
   end do
   33 format(3f10.4)
   deallocate(pol0,pol,eps)

    ELSEIF ( trim(card) == 'BORN_CHARGES') THEN
  write(6,*) '**************************************************'
  WRITE( stdout, '("Born effective charges")')
  write(6,*) '**************************************************'

! npol=1 --> code reads forces due to Efield along
!            x,y,z directions. In this exact order.
! npol=2 --> code reads forces due to Efield along
!            x,y,z,-x,-y,-z (in this order)

  npol=npol_zeu
  de=de_zeu*sqrt(2.0d0)

  allocate (F0(nat,3))
  allocate (Fij(npol,nat,3,3))
  allocate (zeta(3,3,nat))

  F0(:,:)=0.0d0
  Fij(:,:,:,:)=0.0d0
  zeta=0.0d0

  do k=1,nat
    read(5,*) (F0(k,j),j=1,3)
  end do
  
  do p=1,npol
   do i=1,3   ! 3 --> x,y,z
     do k=1,nat
       read(5,*) (Fij(p,k,i,j), j=1,3)
     end do
    end do
  end do
  do k=1,nat
    do i=1,3
      do j=1,3
        if (npol==2) then
          zeta(i,j,k)=(Fij(1,k,i,j)-F0(k,j))-(Fij(2,k,i,j)-F0(k,j))
          zeta(i,j,k)=0.5*zeta(i,j,k)/de
        else
          zeta(i,j,k)=(Fij(1,k,i,j)-F0(k,j))
          zeta(i,j,k)=zeta(i,j,k)/de
        end if
      end do
    end do
  end do
  ! impose ASR
  do i=1,3
    do j=1,3
      sum=0.0d0
        do k=1,nat
          sum = sum + zeta(i,j,k) 
        end do
        do k=1,nat
          zeta(i,j,k) = zeta(i,j,k) - sum/nat
        end do
    end do
  end do
  ! symmetrize
  call symtensor (nat, zeta)

  do k=1,nat
  write(6,54) 'atom ', k
    do i=1,3
      write(6,53) (zeta(i,j,k),j=1,3)
    end do
  end do
  53 format(3f14.8)
  54 format(a,2x,i5)
  deallocate(F0,Fij,zeta)

    ELSE
         !
         IF ( ionode ) &
            WRITE( stdout,'(A)') 'Warning: card '//trim(input_line)//' ignored'
         !

    END IF
    GO TO 100

120 CONTINUE

end program fd_raman


   SUBROUTINE cart_to_crys_mat3 ( mat3 )
     !-----------------------------------------------------------------------
     !
     USE kinds,      ONLY : dp
     USE cell_base,  ONLY : at

     IMPLICIT NONE
     !
     REAL(DP), intent(INOUT) :: mat3(3,3,3)
     !
     REAL(DP) :: work(3,3,3)
     INTEGER :: i,j,k,l,m,n
     !
     work(:,:,:) = 0.0d0
     DO i = 1, 3
        DO j = 1, 3
           DO k = 1, 3
              DO l = 1, 3
                 DO m = 1, 3
                    DO n = 1, 3
                       work (i, j, k) = work (i, j, k) +  &
                          mat3 (l, m, n) * at (l, i) * at (m, j) * at (n, k)
                    END DO
                 END DO
              END DO
           END DO
        END DO
     END DO
     mat3(:,:,:) = work (:,:,:)
     !
   END SUBROUTINE cart_to_crys_mat3
