! Copyright (C) 2006-2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.d0,0.d0)
#define ONE (1.d0,0.d0)

!-----------------------------------------------------------------------
PROGRAM wannier_plot
!-----------------------------------------------------------------------
!
! This program plots charge density of selected wannier function in
! IBM Data Explorer format

  USE io_global, ONLY: stdout, ionode, ionode_id
  USE kinds,         ONLY : DP
  USE io_files,      ONLY : prefix, tmp_dir
  USE wannier_new,   ONLY : nwan, plot_wan_num, plot_wan_spin
  USE mp,            ONLY : mp_bcast
  USE mp_world,      ONLY : world_comm
  USE io_global,     ONLY : ionode, stdout
  USE mp_global,     ONLY : mp_startup
  USE environment,   ONLY : environment_start, environment_end

  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER(len=256) :: outdir
  INTEGER :: ios,nc(3),n0(3)
  NAMELIST /inputpp/ outdir, prefix, nwan, plot_wan_num, plot_wan_spin, nc, n0
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'WANNIER_PLOT' )

  ios = 0
  !
  IF ( ionode ) THEN
     !
     !   set default values for variables in namelist
     !
     CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
     IF ( trim( outdir ) == ' ' ) outdir = './'
     prefix ='pwscf'
     nwan = 0
     plot_wan_spin=1

     nc(1) = 3
     nc(2) = 3
     nc(3) = 3
     n0(1) = -1
     n0(2) = -1
     n0(3) = -1
     !
     CALL input_from_file ( )
     !
     READ (5, inputpp, iostat=ios )
     !
     tmp_dir = trimcheck (outdir)
  ENDIF
  !
  CALL mp_bcast( ios, ionode_id, world_comm )
  IF ( ios /= 0 ) CALL errore('wannier_ham','reading inputpp namelist',abs(ios))
  CALL read_file
  CALL openfil_pp

  CALL wannier_init(.true.)

  !debug
  WRITE(stdout,'(5x,"Calling plot_wannier for wannier",i3)') plot_wan_num
  !end of debug
  CALL plot_wannier(nc,n0)
  !debug
  WRITE(stdout,'(5x,"Calling plot_atoms")')
  !end of debug
  CALL plot_atoms()

  CALL stop_pp

  CALL environment_end ( 'WANNIER_PLOT' )

  CALL wannier_clean()

END PROGRAM wannier_plot

SUBROUTINE plot_wannier(nc,n0)

  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE io_files
  USE kinds,         ONLY : DP
  USE wannier_new,   ONLY : nwan,plot_wan_num,plot_wan_spin
  USE klist,         ONLY : nks, xk, wk, ngk, igk_k
  USE lsda_mod,      ONLY : isk, lsda, nspin
  USE wvfct,         ONLY : nbnd, npwx, g2kin
  USE constants,     ONLY : rytoev , tpi
  USE buffers
  USE symm_base,     ONLY : nsym
  USE basis,         ONLY : swfcatom
  USE fft_base,      ONLY : dffts, dfftp
  USE fft_interfaces,ONLY : invfft
  USE gvect
  USE gvecs
  USE cell_base
  USE ions_base, ONLY : nat, ntyp=>nsp, ityp, tau, atm, zv
  USE vlocal,    ONLY : strf


  IMPLICIT NONE
  INTEGER, INTENT(in) :: nc(3), n0(3)
  INTEGER :: npw, is, i,j, k, ik, n, ir, ios, n1, n2, n3,i1,j1,k1
  COMPLEX(DP) :: phase
  COMPLEX(DP), ALLOCATABLE :: wan_func(:,:), pp_ort(:,:), psic(:), psic3(:,:,:), psic3_0(:,:,:), psic_sum(:,:,:,:), paux(:,:)
  real(DP), ALLOCATABLE :: rho(:,:,:,:), raux(:)
  real(DP) :: r(3)

  IF (nsym>1) THEN
     CALL errore('wannier_cmptn','k-points set is in the irreducible brillouin zone - not implemented',1)
  ENDIF

  ALLOCATE(wan_func(npwx,nwan))
  ALLOCATE(psic(dffts%nnr))
  ALLOCATE(psic3(dffts%nr1x,dffts%nr2x,dffts%nr3x))
  ALLOCATE(psic3_0(dffts%nr1x,dffts%nr2x,dffts%nr3x))
  ALLOCATE(psic_sum(nc(1)*dffts%nr1x,nc(2)*dffts%nr2x,nc(3)*dffts%nr3x,nspin))
  ALLOCATE(rho(nc(1)*dffts%nr1x,nc(2)*dffts%nr2x,nc(3)*dffts%nr3x,nspin))

  CALL init_us_1
  CALL init_at_1

  CALL struc_fact (nat, tau, ntyp, ityp, ngm, g, bg, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
       strf, eigts1, eigts2, eigts3)

  wan_func = ZERO
  psic3 = ZERO
  psic3_0 = ZERO
  psic_sum = ZERO

  DO ik = 1, nks
     npw = ngk(ik)
     is  = isk(ik)

     wan_func = ZERO
     CALL get_buffer( wan_func, nwordwf, iunwf, ik)

     psic(1:dffts%nnr) = ZERO
     rho = ZERO
     DO j = 1, npw
        psic (nls (igk_k(j,ik) ) ) = wan_func (j, plot_wan_num)
     ENDDO

     CALL invfft ('Wave', psic, dffts)

     DO k=1, dffts%nr3x
        DO j=1,dffts%nr2x
           DO i=1,dffts%nr1x
              n = i + (j-1)*dffts%nr1x + (k-1)*dffts%nr2x*dffts%nr1x
              psic3_0(i,j,k) = psic(n)
           ENDDO
        ENDDO
     ENDDO

     DO k=1, (dffts%nr3x-1)*nc(3)
        DO j=1, (dffts%nr2x-1)*nc(2)
           DO i=1, (dffts%nr1x-1)*nc(1)
              r = n0(1)*at(:,1)+n0(2)*at(:,2)+n0(3)*at(:,3)
              r = r + dble(i-1)*at(:,1)/dble(dffts%nr1x-1) + &
                      dble(j-1)*at(:,2)/dble(dffts%nr2x-1) + &
                      dble(k-1)*at(:,3)/dble(dffts%nr3x-1)
              phase = cos(tpi*(xk(1,ik)*r(1)+xk(2,ik)*r(2)+xk(3,ik)*r(3))) + &
          (0.d0,1.d0)*sin(tpi*(xk(1,ik)*r(1)+xk(2,ik)*r(2)+xk(3,ik)*r(3)))

              i1 = i - floor(dble(i-0.01)/dble(dffts%nr1x-1))*(dffts%nr1x-1)
              j1 = j - floor(dble(j-0.01)/dble(dffts%nr2x-1))*(dffts%nr2x-1)
              k1 = k - floor(dble(k-0.01)/dble(dffts%nr3x-1))*(dffts%nr3x-1)
              psic_sum(i,j,k,is) = psic_sum(i,j,k,is)+ &
                   cmplx(wk(ik),0.d0,kind=DP)*psic3_0(i1,j1,k1)*phase
           ENDDO
        ENDDO
     ENDDO

  ENDDO !ik

  rho = 0.d0

  DO is=1, nspin
     DO i=1, dffts%nr1x*nc(1)
        DO j=1, dffts%nr2x*nc(2)
           DO k=1,dffts%nr3x*nc(3)
              rho(i,j,k,is) = dreal(psic_sum(i,j,k,is))**2 +  &
                              aimag(psic_sum(i,j,k,is))**2
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  OPEN (10, file='wannier.plot.dx', err = 100, iostat = ios)
100 CALL errore ('plot_wannier', 'Opening out file', abs (ios) )

  ! I want to write .dx file for dataexplorer
  WRITE(10,'(a36,3i6)') 'object 1 class gridpositions counts ', &
      dffts%nr3x*nc(3), dffts%nr2x*nc(2), dffts%nr1x*nc(1)
  WRITE(10,*) 'origin', n0(1)*at(:,1)+n0(2)*at(:,2)+n0(3)*at(:,3)
  WRITE(10,'(a5, 3f9.5)') 'delta', (at(i,1)/(1.d0*(dffts%nr3x-1)),i=1,3)
  WRITE(10,'(a5, 3f9.5)') 'delta', (at(i,2)/(1.d0*(dffts%nr2x-1)),i=1,3)
  WRITE(10,'(a5, 3f9.5)') 'delta', (at(i,3)/(1.d0*(dffts%nr1x-1)),i=1,3)
  WRITE(10,'(a38,3i6)') 'object 2 class gridconnections counts ', &
      dffts%nr3x*nc(3), dffts%nr2x*nc(2), dffts%nr1x*nc(1)
  WRITE(10,*) 'attribute "element type" string "cubes"'
  WRITE(10,*) 'attribute "ref" string "positions"'
  WRITE(10,'(a44,i10,a13)') 'object 3 class array type float rank 0 items', &
       dffts%nr3x*nc(3)*dffts%nr2x*nc(2)*dffts%nr1x*nc(1), 'data follows'

  DO i=1, dffts%nr3x*nc(3)
     DO j=1,dffts%nr2x*nc(2)
        DO k=1,dffts%nr1x*nc(1)
           WRITE(10,'(f13.7)') rho(k,j,i,plot_wan_spin)
           ! write(10,'(f13.7)') aimag(psic_sum(k,j,i,plot_wan_spin))
        ENDDO
     ENDDO
  ENDDO

  WRITE(10,'(a34)') 'attribute "dep" string "positions"'
  WRITE(10,*) 'object "regular positions regular connections" class field'
  WRITE(10,*) 'component "positions" value 1'
  WRITE(10,*) 'component "connections" value 2'
  WRITE(10,*) 'component "data" value 3'
  WRITE(10,*) 'end'

  CLOSE(10)

  DEALLOCATE(wan_func)
  DEALLOCATE(psic)
  DEALLOCATE(psic3)
  DEALLOCATE(psic3_0)
  DEALLOCATE(psic_sum)
  DEALLOCATE(rho)


END SUBROUTINE plot_wannier


SUBROUTINE plot_atoms
  USE io_global, ONLY: stdout
  USE kinds, ONLY: DP
  USE ions_base, ONLY: tau, nat, ityp, zv
  IMPLICIT NONE
  INTEGER :: i,na, ios

  OPEN (20, file='atoms.plot.dx', err = 200, iostat = ios)
200 CALL errore ('plot_wannier', 'Opening out atoms file', abs (ios) )

  WRITE(20,*) 'object 1 class array type float rank 1 shape 3 items', nat,' data follows'
  DO na = 1, nat
     WRITE(20,'(3f9.5)') (tau(i,na),i=1,3)
  ENDDO
  WRITE(20,*) 'object 2 class array type float rank 0 items', nat,' data follows'
  DO na = 1, nat
     WRITE(20,*) zv(ityp(na))
  ENDDO
  WRITE(20,*) 'attribute "dep" string "positions"'
  WRITE(20,*) 'object "irregular positions" class field'
  WRITE(20,*) 'component "positions" value 1'
  WRITE(20,*) 'component "data" value 2'
  WRITE(20,*) 'end'
  CLOSE(20)
END SUBROUTINE plot_atoms
