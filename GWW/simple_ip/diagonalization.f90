
SUBROUTINE diagonalization(q,sh,input,eig,ik,kptns)
  !------------------------------------------------------------------
  !
  ! This routine calculates H(q)_ij in the optimal basis for a given q and then diagonalizes it.
  ! It reconstructs V_nonloc_ij from the interpolated projectors
  ! h_level = 0 only kinetic part
  ! h_level = 1 kinetic part + local
  ! h_level > 1 kinetic part + local + non-local
  !
  USE kinds, ONLY : DP 
  USE constants, ONLY : pi , rytoev
  USE simple_ip_objects
  USE input_simple_ip
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(kind=DP), INTENT(in) :: q(3)  ! k-point (units 2*pi/alat)
  INTEGER, INTENT(in)  :: ik         ! k-point in the local dense grid
  TYPE(shirley) :: sh
  TYPE(eigen) :: eig
  TYPE(kpoints) :: kptns
  TYPE(input_options_simple_ip) :: input
  COMPLEX(kind=DP), DIMENSION(:,:), ALLOCATABLE  :: h_mat 
  COMPLEX(kind=DP), DIMENSION(:), ALLOCATABLE :: work
  REAL(kind=DP), DIMENSION(:), ALLOCATABLE :: rwork , energies_tmp
  INTEGER, DIMENSION(:), ALLOCATABLE :: iwork , ifail
  COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: b_nc , csca_nc 
  COMPLEX(kind=DP), DIMENSION(:,:), ALLOCATABLE :: b_c  ,  csca_mat , csca_mat2 , deeaux , csca_nc2 
  REAL(kind=DP) :: k2  , abstol , vl , vu 
  INTEGER :: i, ndim, nt, na, j, ih, jh, ikb, jkb, lwork, info, low_eig, high_eig ,ntot_e , num_bands
  !
  !call start_clock('diagonalization')
  !
  if(associated(eig%energy)) deallocate(eig%energy)
  allocate(eig%energy(sh%num_bands))
  if(associated(eig%wave_func)) deallocate(eig%wave_func)
  allocate(eig%wave_func(sh%ntot_e,sh%num_bands))
  !
  allocate(h_mat(sh%ntot_e,sh%ntot_e),energies_tmp(sh%ntot_e))
  h_mat(1:sh%ntot_e,1:sh%ntot_e) = 0.d0
  energies_tmp(1:sh%ntot_e) = 0.d0
  !
  ! Construction of the matrix H(q) = (1/2)*[q**2 + q*K1 + K0] + Vloc + Vnonloc(q) in the Shirley basis
  ! The factor 1/2 in front of the kinetic term is not needed in the QE units
  !
  ! Diagonal part: q^2  
  k2 = q(1)*q(1) + q(2)*q(2) + q(3)*q(3)
  k2 = k2 * (2.0*pi/sh%alat)**2
  do i=1,sh%ntot_e
   h_mat(i,i) = k2  
  enddo
  !
  ! Kinetic term: q*K1
  do i=1,3
   h_mat(1:sh%ntot_e,1:sh%ntot_e) = h_mat(1:sh%ntot_e,1:sh%ntot_e) + &
         & q(i)*sh%h1(1:sh%ntot_e,1:sh%ntot_e,i) * (2.0*pi/sh%alat)  !*0.5
  enddo  
  !
  ! Kinetic term: K0
  h_mat(1:sh%ntot_e,1:sh%ntot_e) = h_mat(1:sh%ntot_e,1:sh%ntot_e) + sh%h0(1:sh%ntot_e,1:sh%ntot_e)
  !
  ! Vloc
  if (input%h_level>0) then
     h_mat(1:sh%ntot_e,1:sh%ntot_e) = h_mat(1:sh%ntot_e,1:sh%ntot_e) + sh%Vloc(1:sh%ntot_e,1:sh%ntot_e)
  endif
  !
  ! Vnonloc
  !call start_clock('diago_vnloc')
  if (input%h_level>1) then
   !
   if (sh%noncolin) then
     ndim = sh%nkb * sh%npol * sh%ntot_e
     allocate(b_nc(sh%nkb,sh%npol,sh%ntot_e))
     if(input%nonlocal_interpolation) then
       call trilinear_parallel_nc(kptns,ik,sh,b_nc)
     else
       b_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb, 1:sh%npol, 1:sh%ntot_e,ik)
     endif
   else
     ndim = sh%nkb * sh%ntot_e
     allocate(b_c(sh%nkb,sh%ntot_e))
     if(input%nonlocal_interpolation) then
       call trilinear_parallelc(kptns,ik,sh,b_c)
     else
       b_c(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb, 1:sh%ntot_e,ik)
     endif
   endif
   !
   ! Build Vnonloc
   if (sh%noncolin) then
     !
     allocate(csca_nc(sh%nkb,sh%npol,sh%ntot_e),csca_nc2(sh%ntot_e,sh%ntot_e))
     csca_nc(:,:,:) = 0.d0
     do nt=1,sh%ntyp
      if (sh%nh(nt)==0) CYCLE
      do na=1,sh%nat
       if ( sh%ityp(na) == nt ) then
           do j=1,sh%ntot_e
            do jh=1,sh%nh(nt)
             jkb = sh%indv_ijkb0(na)+jh
             do ih=1,sh%nh(nt)
               ikb = sh%indv_ijkb0(na)+ih  
               ! 
               csca_nc(ikb,1,j) = csca_nc(ikb,1,j) +    & 
                              sh%deeq_nc(ih,jh,na,1)*b_nc(jkb,1,j)+ & 
                              sh%deeq_nc(ih,jh,na,2)*b_nc(jkb,2,j) 
               csca_nc(ikb,2,j) = csca_nc(ikb,2,j)  +   & 
                              sh%deeq_nc(ih,jh,na,3)*b_nc(jkb,1,j)+&
                              sh%deeq_nc(ih,jh,na,4)*b_nc(jkb,2,j) 
               !
               enddo
            enddo
           enddo
       endif    
      enddo
     enddo
     !
     call ZGEMM ('C', 'N', sh%ntot_e, sh%ntot_e, sh%npol*sh%nkb, ( 1.d0, 0.d0 ) , b_nc, &
                   sh%npol*sh%nkb, csca_nc, sh%npol*sh%nkb, ( 0.d0, 0.d0 ) , csca_nc2, sh%ntot_e)
     h_mat(1:sh%ntot_e,1:sh%ntot_e) = h_mat(1:sh%ntot_e,1:sh%ntot_e) + csca_nc2(1:sh%ntot_e,1:sh%ntot_e)
     !
     deallocate(csca_nc,csca_nc2)
     !
   else
     !
     allocate(csca_mat(sh%nkb,sh%ntot_e), csca_mat2(sh%ntot_e,sh%ntot_e))
     !
     do nt=1,sh%ntyp
      if (sh%nh(nt)==0) cycle
      allocate ( deeaux(sh%nh(nt),sh%nh(nt)) )
      do na=1,sh%nat
       if ( sh%ityp(na) == nt ) then
          !
          deeaux(1:sh%nh(nt),1:sh%nh(nt)) = sh%deeqc(1:sh%nh(nt),1:sh%nh(nt),na) 
          call ZGEMM('N','N', sh%nh(nt), sh%ntot_e, sh%nh(nt), (1.d0,0.d0), &
                     deeaux, sh%nh(nt), b_c(sh%indv_ijkb0(na)+1,1), sh%nkb, &
                     (0.d0, 0.d0), csca_mat(sh%indv_ijkb0(na)+1,1), sh%nkb )
          !
       endif    
      enddo
      deallocate(deeaux)
     enddo
     !
     CALL ZGEMM( 'C', 'N', sh%ntot_e, sh%ntot_e, sh%nkb, ( 1.d0, 0.d0 ) , b_c, &
                   sh%nkb, csca_mat, sh%nkb, ( 0.d0, 0.d0 ) , csca_mat2, sh%ntot_e )
     h_mat(1:sh%ntot_e,1:sh%ntot_e) = h_mat(1:sh%ntot_e,1:sh%ntot_e) + csca_mat2(1:sh%ntot_e,1:sh%ntot_e)
     !
     deallocate(csca_mat, csca_mat2)
     !
   endif
   !
   if (sh%noncolin) then
     deallocate(b_nc)
   else
     deallocate(b_c)
   endif
   !
  endif
  !call stop_clock('diago_vnloc')
  !
  ! Hamiltonian diagonalization 
  allocate(rwork(7*sh%ntot_e),iwork(5*sh%ntot_e),ifail(sh%ntot_e))
  allocate(work(1))
  lwork = -1
  abstol = -1.d0
  vl = 0.d0
  vu = 0.d0
  low_eig = sh%num_nbndv(1)-sh%num_val+1
  high_eig = sh%num_nbndv(1)+sh%num_cond
  num_bands = sh%num_bands
  !
  CALL ZHEEVX( 'V', 'I', 'U', sh%ntot_e, h_mat, sh%ntot_e, vl, vu, low_eig, high_eig, abstol, &
             & num_bands, energies_tmp, eig%wave_func, sh%ntot_e, work, lwork, rwork, iwork, ifail, info)
  !
  lwork = int(work(1))
  deallocate(work)
  allocate(work(1:lwork))
  !
  !Diagonalization of H(k) using ZHEEVX (only the first sh%num_bands energies/wavefunctions are explicitly calculated)
  !call start_clock('diago_zheevx')
  call ZHEEVX( 'V', 'I', 'U', sh%ntot_e, h_mat, sh%ntot_e, vl, vu, low_eig, high_eig, abstol, &
             & num_bands, energies_tmp, eig%wave_func, sh%ntot_e, work, lwork, rwork, iwork, ifail, info)
  !call stop_clock('diago_zheevx')
  !
  if (info /= 0) write(stdout,*) 'ERROR in the diagonalization of the k-dependent Hamiltonian', info
  eig%energy(1:sh%num_bands) = energies_tmp(1:sh%num_bands)
  !
  deallocate(h_mat, rwork, work, iwork, ifail, energies_tmp)
  !
  !call stop_clock('diagonalization')
  !
END SUBROUTINE diagonalization

