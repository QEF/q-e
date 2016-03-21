subroutine check_basis(numwprod,npw)
! checking if the polarizability basis is orthonormal

USE fft_custom_gwl
!USE io_files,             ONLY : find_free_unit, prefix, diropn
USE io_files,             ONLY : prefix, diropn
USE wavefunctions_module, ONLY :  psic
USE mp,          ONLY :mp_barrier
use io_global, ONLY : stdout, ionode 
USE kinds, ONLY : DP
USE mp,             ONLY : mp_sum
use mp_world, ONLY : mpime
USE mp_world,             ONLY : world_comm
USE gvect,          ONLY : gstart,ngm_g


implicit none
INTEGER, EXTERNAL :: find_free_unit
REAL(kind=DP), EXTERNAL :: ddot

integer numwprod
integer npw

COMPLEX(kind=DP), ALLOCATABLE :: p_basis(:,:)
real(kind=DP) :: prod
INTEGER ::iungprod
INTEGER ::iunnorm
LOGICAL       :: exst
integer :: ii,jj



iungprod = find_free_unit()
allocate(p_basis(npw,numwprod))
CALL diropn( iungprod, 'wiwjwfc_red', npw*2, exst )

do ii=1,numwprod
   call davcio(p_basis(:,ii),npw*2,iungprod,ii,-1)
enddo

call mp_barrier(world_comm)
close(iungprod)


! check normalization
if(ionode) then
   iunnorm = find_free_unit()
   open(iunnorm, file='pol_basis_norm.dat',status='unknown',form='formatted')
   write(iunnorm,*) '# Pol_vector_i, Norm'
endif

do ii=1,numwprod   
   prod=2.d0*ddot(2*npw,p_basis(:,ii),1,p_basis(:,ii),1)
   if (gstart==2) prod=prod-p_basis(1,ii)*p_basis(1,ii)
   call mp_sum(prod,world_comm)
!   prod=prod/ngm_g
   if(ionode) write(iunnorm,*) ii,prod
enddo   

if(ionode) close(iunnorm)

!check orthogonality

if(ionode) then
   iunnorm = find_free_unit()
   open(iunnorm, file='pol_basis_ortho.dat',status='unknown',form='formatted')
   write(iunnorm,*) '# Pol_vector_i, #Polarization vector j, Product'
endif

do ii=1,numwprod
   do jj=ii+1,numwprod   
      prod=2.d0*ddot(2*npw,p_basis(:,ii),1,p_basis(:,jj),1)
      if (gstart==2) prod=prod-p_basis(1,ii)*p_basis(1,jj)
      call mp_sum(prod,world_comm)
!      prod=prod/ngm_g
      if(ionode) write(iunnorm,*) ii,jj,prod
   enddo
enddo   

if(ionode) close(iunnorm)
return
end subroutine
