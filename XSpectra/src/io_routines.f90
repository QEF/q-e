

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE read_core_abs(filename,core_wfn,nl_init)
  !--------------------------------------------------------------------------
  USE kinds,   ONLY: DP
  USE atom,    ONLY: rgrid
  !USE atom,        ONLY : mesh     !mesh(ntypx) number of mesh points
  USE xspectra,ONLY: xiabs
  USE io_global,       ONLY : ionode, stdout

  IMPLICIT NONE

  INTEGER :: i, ierr, nbp, iblind
  INTEGER, dimension(2) :: nl_init
  CHARACTER (LEN=80) :: filename
  REAL(KIND=dp):: x
  REAL(KIND=dp):: core_wfn(*)

  !WRITE(6,*) 'xmesh=',rgrid(xiabs)%mesh
  open(unit=33,file=filename,form='formatted',iostat=ierr,status='old',err=123)

123   if( ierr == 29 )  then
          if( ionode ) write(stdout, &
          '("ERROR: core wavefunction file ",A,&
         & " does not exist or is not located in the right folder !")') &
           trim(filename)
          call stop_xspectra
        else if( ierr /= 0 ) then
          if( ionode ) write(stdout, '("ERROR reading the core wavefunction file")')
          call stop_xspectra
        else
          continue
        end if

  rewind(33)

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !<OB>           Brute force identification of the core state:
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  nbp = rgrid(xiabs)%mesh
! Determine how many lines to skip before reading the good core function
  select case( nl_init(1) )
    case(1)
      iblind = 1                       !1s
    case(2)
      if( nl_init(2) == 0 ) then       !2s
        iblind = nbp + 2
      else                             !2p
        iblind = 2*nbp + 3
      end if
    case(3)
      if( nl_init(2) == 0 ) then       !3s
        iblind = 3*nbp + 4
      else if( nl_init(2) == 1 ) then  !3p
        iblind = 4*nbp + 5
      else
        iblind = 5*nbp + 6             !3d
      end if
  end select

  do i = 1, iblind
    READ(33,*) 
  end do

  DO i=1,nbp
   READ(33 ,*) x,core_wfn(i)
  ENDDO

  close(33)

END SUBROUTINE read_core_abs


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE read_save_file(a,b,xnorm,ncalcv,x_save_file,core_energy)
  !----------------------------------------------------------------------------
  ! This routine reads the x_save_file (default name: xanes.sav)
  ! (automatically named prefix_K.save for a K edge, in the next version). 
  ! This routine is used only if xonly_plot=.true.
  !----------------------------------------------------------------------------
  USE kinds,       ONLY: DP
  USE constants,   ONLY: rytoev
  USE klist,       ONLY: nks, nkstot
  USE xspectra,    ONLY: xnitermax, xang_mom, xiabs,   &
                         n_lanczos, save_file_version, &
                         save_file_kind, calculated,   &
                         xe0, xe0_default, xe0_ry
  USE ener,        ONLY: ef
  USE io_global,   ONLY: stdout, ionode
  USE lsda_mod,    ONLY: nspin, lsda

  IMPLICIT NONE
  ! Arguments
  REAL(dp), INTENT (INOUT) ::  a(xnitermax,n_lanczos,nks)
  REAL(dp), INTENT (INOUT) ::  b(xnitermax,n_lanczos,nks)     
  REAL(dp), INTENT (INOUT) ::  xnorm(n_lanczos,nks) 
  INTEGER, INTENT (INOUT)  ::  ncalcv(n_lanczos,nks)
  CHARACTER(LEN=256), INTENT (IN) :: x_save_file
  REAL(dp), INTENT (OUT)   ::  core_energy
  ! Local variables
  ! NB: The '_r' suffix means 'read in x_save_file'
  INTEGER  :: ierr, nkstot_r
  INTEGER  :: xm_r, nc_r, ncomp_max
  INTEGER  :: i, j, k, ncalcv_max
  INTEGER  :: calculated_all(n_lanczos,nkstot)
  INTEGER, ALLOCATABLE :: ncalcv_all(:,:)
  REAL(dp) :: xepsilon_r(3), xkvec_r(3)
  REAL(dp), ALLOCATABLE :: a_all(:,:), b_all(:,:), xnorm_all(:,:), aux(:,:)

  ALLOCATE(a_all(xnitermax,nkstot))
  ALLOCATE(b_all(xnitermax,nkstot))
  ALLOCATE(xnorm_all(n_lanczos,nkstot))
  ALLOCATE(ncalcv_all(n_lanczos,nkstot))
  ALLOCATE(aux(xnitermax,nks))

  a(:,:,:) = 0.d0
  b(:,:,:) = 0.d0
  xnorm(:,:) = 0.d0
  ncalcv(:,:) = 0
  calculated_all(:,:) = 0


  OPEN ( UNIT = 10, FILE = x_save_file, FORM = 'FORMATTED', &
       STATUS = 'UNKNOWN', IOSTAT = ierr )
  CALL errore( 'iosys', 'x_save_file ' // TRIM( x_save_file ) // &
       & ' not found' , ierr )
  WRITE(stdout,'(5x,"x_save_file name: ",a)') TRIM( x_save_file )
  REWIND(10)

  IF (save_file_version == 0) then
     WRITE(stdout,'(5x,a)') 'x_save_file version: old'
  !ELSEIF(save_file_version == 1) then
  !   WRITE(stdout,'(5x,a)') 'x_save_file version: 1'
  ELSE
     WRITE(stdout,'(5x,a,i3)') 'x_save_file version: ', save_file_version
     DO i = 1, 6 
        READ(10,*)      
     ENDDO
  ENDIF
  WRITE(stdout,*)

  READ(10,*) lsda, nspin
  WRITE(stdout,'(5x,a,i2)') 'nspin:',nspin
  READ(10,*) xm_r, nkstot_r, xnitermax
  WRITE(stdout,'(5x,a,i4)') 'number of k-points:',nkstot
  WRITE(stdout,*)
  WRITE(stdout,'(5x,a,i4)') 'final-state angular momentum (xm_r): ',xm_r
  IF (xm_r==1) THEN
     WRITE(stdout,'(5x,a)') ' => electric-dipole approximation'
  ELSEIF (xm_r==2) THEN
     WRITE(stdout,'(5x,a)') ' => electric-quadrupole approximation'
  ELSE
     WRITE(stdout,'(5x,a)') 'Wrong value of xm_r: STOP'
     CALL stop_xspectra
  ENDIF
  IF(xm_r.NE.xang_mom) & 
     CALL errore('read_save_file','xm_r is different from xang_mom=',xang_mom)

  READ(10,*) ncalcv_max
  IF(ncalcv_max.GT.xnitermax) THEN
     WRITE(stdout,'(5x,a,i5)') 'ncalcv_max=',ncalcv_max
     CALL errore('read_save_file','ncalcv_max is grater than xnitermax=', &
                 xnitermax)
  ENDIF

  WRITE(stdout,*)
  IF (save_file_version < 2) THEN
     READ(10,*) core_energy
  ELSE
     READ(10,*) core_energy, ef
     WRITE(stdout,'(5x,a,f9.4)') 'Fermi level [eV]:', ef 
  ENDIF
  WRITE(stdout,'(5x,a,f10.3,/)') 'core energy [eV]:', core_energy

  READ(10,*) (xkvec_r(i),i=1,3)
  !WRITE(stdout,*) '---------------------------------------------------------'
  !WRITE(stdout,*) 'xkvec read from savefile'
  !WRITE(stdout,*) (xkvec_r(i),i=1,3)
  READ(10,*) (xepsilon_r(i),i=1,3)
  !WRITE(stdout,*) 'xepsilon read from file'
  !WRITE(stdout,*) (xepsilon_r(i),i=1,3)
  !WRITE(stdout,*) '---------------------------------------------------------'
  !write(stdout,*) 'n_lanczos=', n_lanczos
  WRITE(stdout,'(5x,a,1x,3(f10.6,1x))') 'xepsilon [Cartesian frame]:',&
       (xepsilon_r(i),i=1,3)
  IF (xm_r==2) WRITE(stdout,'(5x,a,1x,3(f10.6,1x))') &
       'xkvec [Cartesian frame]:', (xkvec_r(i),i=1,3)
  WRITE(stdout,*)

  DO i = 1, n_lanczos
     IF (TRIM(save_file_kind).eq.'unfinished') THEN
       READ(10,*) (calculated_all(i,j),j=1,nkstot)
     ENDIF
     READ(10,*) (xnorm_all(i,j),j=1,nkstot)
     READ(10,*) (ncalcv_all(i,j),j=1,nkstot)
     READ(10,*) ((a_all(j,k),j=1,ncalcv_max),k=1,nkstot)
     READ(10,*) ((b_all(j,k),j=1,ncalcv_max),k=1,nkstot)
#if defined(__MPI)
     CALL poolscatter(xnitermax,nkstot,a_all,nks,aux)
     a(1:xnitermax,i,1:nks)=aux(1:xnitermax,1:nks)
     CALL poolscatter(xnitermax,nkstot,b_all,nks,aux)
     b(1:xnitermax,i,1:nks)=aux(1:xnitermax,1:nks)
#else
     a(1:xnitermax,i,1:nkstot)=a_all(1:ncalcv_max,1:nkstot)
     b(1:xnitermax,i,1:nkstot)=b_all(1:ncalcv_max,1:nkstot)
#endif

  ENDDO
  CLOSE(10)

#if defined(__MPI) 
  CALL poolscatter(n_lanczos,nkstot,xnorm_all,nks,xnorm)
  CALL ipoolscatter(n_lanczos,nkstot,ncalcv_all,nks,ncalcv)
  CALL ipoolscatter(n_lanczos,nkstot,calculated_all,nks,calculated)
#else
  IF(nks.NE.nkstot) THEN
     CALL errore('read_save_file','nks\=nkstot',1)
  ENDIF
  xnorm(1:n_lanczos,1:nkstot) = xnorm_all(1:n_lanczos,1:nkstot)
  ncalcv(1:n_lanczos,1:nkstot) = ncalcv_all(1:n_lanczos,1:nkstot)
  calculated(1:n_lanczos,1:nkstot) = calculated_all(1:n_lanczos,1:nkstot)
#endif
  DEALLOCATE(a_all)
  DEALLOCATE(b_all)
  DEALLOCATE(xnorm_all)
  DEALLOCATE(ncalcv_all)
  DEALLOCATE(aux)

END SUBROUTINE read_save_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE read_header_save_file(x_save_file)
  !----------------------------------------------------------------------------
  USE kinds,    ONLY: DP
  USE klist,    ONLY: nkstot
  USE lsda_mod, ONLY: nspin,lsda
  USE xspectra, ONLY: save_file_version, save_file_kind, n_lanczos
  USE io_global,ONLY: stdout

  IMPLICIT NONE
  CHARACTER(LEN=256), INTENT (IN) :: x_save_file
  INTEGER    :: ierr, nkstot_r
  INTEGER    :: xm_r
  CHARACTER  :: c

  OPEN ( UNIT = 10, FILE = x_save_file, FORM = 'FORMATTED', &
       STATUS = 'UNKNOWN', IOSTAT = ierr )
  CALL errore( 'iosys', 'x_save_file ' // TRIM( x_save_file ) // &
       & ' not found' , ierr )
  REWIND(10)

  READ(10, '(a1)') c
  REWIND(10)
  IF (c == '#') then
     READ(10, '(20x,i8)') save_file_version
     READ(10, '(20x,a32)') save_file_kind
     READ(10,*) 
     READ(10,'(27x,i4)') n_lanczos
     READ(10,*) 
     READ(10,*) 
  ELSE
     save_file_version=0
     save_file_kind='xanes_old'
     n_lanczos=1
  ENDIF

  READ(10,*) lsda,nspin
  READ(10,*) xm_r,nkstot_r
  nkstot=nkstot_r
  CLOSE(10)

  RETURN 
END SUBROUTINE read_header_save_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE write_save_file(a,b,xnorm,ncalcv,x_save_file)
  !----------------------------------------------------------------------------
  USE kinds,      ONLY: DP
  USE io_global,       ONLY : stdout
  USE constants,  ONLY: rytoev
  USE klist,      ONLY: nks, nkstot
  USE xspectra,   ONLY: xnitermax, xang_mom, xkvec, xepsilon, xiabs, &
                        save_file_version, save_file_kind,           &
                        n_lanczos, calculated, edge
  USE ener,       ONLY: ef
  USE io_global,  ONLY: ionode
  !*apsi  USE uspp_param, ONLY : psd
  USE lsda_mod,   ONLY: nspin,lsda
  USE uspp_param, ONLY: upf
  USE edge_energy, ONLY: getE
  !
  IMPLICIT NONE
  ! Arguments
  REAL(dp), INTENT (IN) :: a(xnitermax,n_lanczos,nks)
  REAL(dp), INTENT (IN) :: b(xnitermax,n_lanczos,nks)     
  REAL(dp), INTENT (IN) :: xnorm(n_lanczos,nks)
  INTEGER,  INTENT (IN) :: ncalcv(n_lanczos,nks)
  CHARACTER(LEN=256), INTENT (IN) :: x_save_file ! could be removed, 
                                                 ! exists in module
  ! Local variables
  INTEGER :: ierr
  INTEGER :: i, j, k 
  INTEGER :: ncalcv_max
  INTEGER :: calculated_all(n_lanczos,nkstot)
  INTEGER, ALLOCATABLE :: ncalcv_all(:,:)
  REAL(dp), ALLOCATABLE :: a_all(:,:), b_all(:,:), xnorm_all(:,:)
  CHARACTER(LEN=8) :: dte

  IF (ionode)  CALL DATE_AND_TIME(date=dte)

  ALLOCATE(a_all(xnitermax,nkstot))
  ALLOCATE(b_all(xnitermax,nkstot))
  ALLOCATE(xnorm_all(n_lanczos,nkstot))
  ALLOCATE(ncalcv_all(n_lanczos,nkstot))

  ncalcv_all(:,:) = 0
  xnorm_all(:,:)  = 0.d0

  ncalcv_all(1:n_lanczos,1:nks) = ncalcv(1:n_lanczos,1:nks)
  xnorm_all(1:n_lanczos,1:nks)  = xnorm(1:n_lanczos,1:nks)
  calculated_all(1:n_lanczos,1:nks) = calculated(1:n_lanczos,1:nks)

#if defined(__MPI)
  CALL poolrecover(xnorm_all,n_lanczos,nkstot,nks)
  CALL ipoolrecover(ncalcv_all,n_lanczos,nkstot,nks)
  CALL ipoolrecover(calculated_all,n_lanczos,nkstot,nks)
#endif

  ncalcv_max = 0
  DO j = 1, n_lanczos
     DO i = 1, nkstot       
        IF (ncalcv_all(j,i).GT.ncalcv_max) ncalcv_max = ncalcv_all(j,i)
     ENDDO
  ENDDO

  IF ( ionode ) THEN
     OPEN ( UNIT = 10, FILE = x_save_file, FORM = 'FORMATTED', &
            STATUS = 'UNKNOWN', IOSTAT = ierr )
     REWIND(10)
     WRITE(10, '(a20,i8)') '# save_file_version=',save_file_version
     WRITE(10, '(a20,a32)') '# save_file_kind   =',save_file_kind
     WRITE(10,'(a7,a8)') '# date=', dte
     WRITE(10,'(a27,i4)') '# number of lanczos stored=',n_lanczos
     WRITE(10,'(a1)') '#'
     WRITE(10,'(a1)') '#'
     WRITE(10,*) lsda,nspin
     WRITE(10,*) xang_mom,nkstot,xnitermax
     WRITE(10,*) ncalcv_max
!     WRITE(10,*) mygetK(upf(xiabs)%psd), ef ! ef here is in eV
     WRITE(10,*) getE(upf(xiabs)%psd,edge), ef ! ef here is in eV
     WRITE(10,*) (xkvec(i),i=1,3)
     WRITE(10,*) (xepsilon(i),i=1,3)
  ENDIF
  DO i=1, n_lanczos
     a_all(:,:) = 0.d0
     b_all(:,:) = 0.d0
     a_all(1:xnitermax,1:nks) =  a(1:xnitermax,i,1:nks)
     b_all(1:xnitermax,1:nks) =  b(1:xnitermax,i,1:nks)
#if defined(__MPI)
     CALL poolrecover(a_all,xnitermax,nkstot,nks)
     CALL poolrecover(b_all,xnitermax,nkstot,nks)
#endif
     IF ( ionode) THEN
        IF (TRIM(save_file_kind).EQ.'unfinished') THEN
          WRITE(10,*) (calculated_all(i,j),j=1,nkstot)
        ENDIF
        WRITE(10,*) (xnorm_all(i,j),j=1,nkstot)
        WRITE(10,*) (ncalcv_all(i,j),j=1,nkstot)
        WRITE(10,*) ((a_all(j,k),j=1,ncalcv_max),k=1,nkstot)
        WRITE(10,*) ((b_all(j,k),j=1,ncalcv_max),k=1,nkstot)
     ENDIF
  ENDDO
  CLOSE(10)

  WRITE(stdout,'(/,5x,a)') &
       'Results of STEP 1 successfully written in x_save_file'
  WRITE(stdout,'(5x,a18,/,5x,a2,2x,a65)') 'x_save_file name: ',&
       '->', x_save_file 
  WRITE(stdout,'(5x,a21,i2)') 'x_save_file version: ', save_file_version
  
  WRITE(stdout,'(/,5x,"... End STEP 1 ...",/)')


  DEALLOCATE(a_all)
  DEALLOCATE(b_all)
  DEALLOCATE(xnorm_all)
  DEALLOCATE(ncalcv_all)

END SUBROUTINE write_save_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE read_gamma_file
  !----------------------------------------------------------------------------
  USE gamma_variable_mod

  IMPLICIT NONE
  INTEGER :: nl, ierr, i

  OPEN ( UNIT = 21, FILE = gamma_file, FORM = 'formatted',&
         STATUS = 'unknown', IOSTAT = ierr)
  CALL errore('io ', 'gamma file '//TRIM(gamma_file)//' not found', abs (ierr))
  REWIND(21)

  nl = 0

  DO
     READ (21,'(a1)',iostat=ierr)
     IF (ierr.NE.0) EXIT
     nl = nl + 1
  ENDDO
  CLOSE(21)

  gamma_lines = nl
  ALLOCATE(gamma_points(nl,2))

  OPEN ( UNIT = 21, FILE = gamma_file, FORM = 'formatted',&
         STATUS = 'unknown', IOSTAT = ierr)
  REWIND(21)

  DO i=1,nl
     READ(21,*) gamma_points(i,1), gamma_points(i,2)
  ENDDO

  CLOSE(21)

END SUBROUTINE read_gamma_file

