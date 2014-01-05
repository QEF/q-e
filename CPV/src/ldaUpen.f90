!
! Copyright (C) 2011-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
MODULE step_penalty
!-------------------------------------------------------------------------
  !
  ! LDA+U with occupation constraint 
  !
  USE kinds
  implicit none
  integer :: natx
  real(DP) :: E_pen = 0.d0
  real(DP), allocatable :: A_pen(:,:), sigma_pen(:), alpha_pen(:)
  logical :: step_pen
  PRIVATE
  PUBLIC :: ldaUpen_init, deallocate_step_pen, write_pen, penalty_e, penalty_f

CONTAINS
  !
  subroutine ldaUpen_init ( natx_, step_pen_, sigma_pen_, alpha_pen_, A_pen_ )
  !-----------------------------------------------------------------------
  !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: natx_
      LOGICAL, INTENT(IN) :: step_pen_
      REAL(DP),INTENT(IN) :: sigma_pen_(natx_), alpha_pen_(natx_), A_pen_(natx_,2)
      
      step_pen=step_pen_
      natx = natx_
      IF ( step_pen ) THEN
         allocate (A_pen(natx,2), sigma_pen(natx), alpha_pen(natx) )
         sigma_pen=sigma_pen_
         alpha_pen=alpha_pen_
         A_pen=A_pen_
      END IF
  END SUBROUTINE ldaUpen_init
  !
  subroutine deallocate_step_pen()
  !-----------------------------------------------------------------------
  !
     IF( ALLOCATED( alpha_pen ) ) DEALLOCATE( alpha_pen )
     IF( ALLOCATED( sigma_pen ) ) DEALLOCATE( sigma_pen )
     IF( ALLOCATED( A_pen ) ) DEALLOCATE( A_pen )
     !
  end subroutine deallocate_step_pen
  !-----------------------------------------------------------------------
  subroutine write_pen (nsp, nspin)
  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  INTEGER, intent(in) :: nsp, nspin
  INTEGER :: is, isp
  !
  if (step_pen) then
     do isp=1,nspin
        write (6,'(6(a,i2,a,i2,a,f8.4,6x))') &
        ('A_pen(',is,',',isp,') =', A_pen(is,isp),is=1,nsp)
     enddo
     write (6,'(6(a,i2,a,f8.4,6x))') &
           ('sigma_pen(',is,') =', sigma_pen(is), is=1,nsp)
     write (6,'(6(a,i2,a,f8.4,6x))') &
        ('alpha_pen(',is,') =', alpha_pen(is), is=1,nsp)
  endif
  END subroutine write_pen
!
!-----------------------------------------------------------------------
  SUBROUTINE penalty_e ( offset, swfc, proj, e_hubbard, hpsi )
!-----------------------------------------------------------------------
!
!      Calculate the energy (added to e_hubbard) and the potential (added
!      to hpsi) due to constraint
!
      USE kinds,              ONLY: dp        
      USE ions_base,          ONLY: na, nat, nsp
      USE gvecw,              ONLY: ngw
      USE electrons_base,     ONLY: nspin, n => nbsp, nx => nbspx, ispin, f
      USE ldaU_cp,            ONLY: Hubbard_U, Hubbard_l, ldmx, nwfcU, ns
      USE dspev_module,       ONLY: dspev_drv
!
      IMPLICIT NONE
      INTEGER,     intent(in) :: offset(nsp,nat)
      REAL(dp),    intent(in) :: proj(nwfcU,n)
      COMPLEX(dp), intent(in) :: swfc(ngw,nwfcU)
      REAL(dp),    intent(inout) :: e_hubbard
      COMPLEX(dp), intent(inout) :: hpsi(ngw,nx)
!
      REAL(dp), allocatable   :: lambda(:), f1(:), vet(:,:)
      REAL(dp) :: x_value, g_value, step_value
      COMPLEX(dp) :: tempsi
      INTEGER :: is, ia, iat, isp, m1, m2, k, i
!
      E_pen=0
      IF ( .NOT. step_pen ) RETURN
      allocate(f1(ldmx*ldmx), vet(ldmx,ldmx), lambda(ldmx) )
      iat=0
      do is = 1,nsp
         do ia = 1, na(is)
            iat = iat + 1
            if (Hubbard_U(is).ne.0.0_dp) then
               do isp = 1, nspin
                  if (A_pen(iat,isp).ne.0.0_dp) then
                     k = 0
                     f1=0.0
                     do m1 = 1, 2 * Hubbard_l(is) + 1
                        do m2 = m1, 2 * Hubbard_l(is) + 1
                           k = k + 1
                           f1 (k) = ns (m2,m1,iat,isp)
                        enddo
                     enddo
                     CALL dspev_drv( 'V', 'L', 2*Hubbard_l(is)+1, f1, &
                                     lambda, vet, ldmx  )
                     x_value=alpha_pen(iat)-lambda(2*Hubbard_l(is)+1)
                     call stepfn(A_pen(iat,isp),sigma_pen(iat),x_value, &
     &                           g_value,step_value)
                     do i=1, n
                        do m1 = 1, 2 * Hubbard_l(is) + 1
                           do m2 = 1, 2 * Hubbard_l(is) + 1
                              tempsi=-1.d0*f(i)*proj (offset(is,ia)+m1,i) * &
                                      vet(m1,2*Hubbard_l(is)+1) * &
                                      vet(m2,2*Hubbard_l(is)+1) * g_value
                              ! add to hpsi
                              call ZAXPY (ngw,tempsi,swfc(1,offset(is,ia)+m2),&
                                          1,hpsi(1,i),1)
                           enddo
                        enddo
                     end do
                     E_pen=E_pen+step_value
                  end if
               enddo
            endif
         enddo
      enddo
      e_hubbard = e_hubbard + E_pen
      deallocate(f1, vet, lambda)
      !
   end subroutine penalty_e
!
!-----------------------------------------------------------------------
   SUBROUTINE penalty_f ( is, iat, dns, forceh )
!-----------------------------------------------------------------------
!
!      Calculate forces due to constraint (added to forceh)
!
      USE kinds,              ONLY: dp        
      USE ions_base,          ONLY: na, nat, nsp
      USE gvecw,              ONLY: ngw
      USE electrons_base,     ONLY: nspin, n => nbsp, nx => nbspx, ispin, f
      USE ldaU_cp,            ONLY: Hubbard_U, Hubbard_l, ldmx, nwfcU, ns
      USE dspev_module,       ONLY: dspev_drv
!
      IMPLICIT NONE
      INTEGER, intent(in) :: is, iat
      REAL(dp), intent(in) :: dns(ldmx,ldmx,nat,nspin)
      REAL(dp), intent(inout) :: forceh
!
      REAL(dp), allocatable   :: lambda(:), f1(:), vet(:,:)
      REAL(dp) :: x_value, g_value, step_value
      COMPLEX(dp) :: tempsi
      INTEGER :: isp, m1, m2, k
!
      IF ( .NOT. step_pen ) RETURN
      allocate(f1(ldmx*ldmx), vet(ldmx,ldmx), lambda(ldmx) )
      do isp = 1, nspin
         if ( (A_pen(iat,isp).ne.0.0) .and. (Hubbard_U(is).ne.0.d0)) then
            k = 0
            f1=0.0
            do m1 = 1, 2 * Hubbard_l(is) + 1
               do m2 = m1, 2 * Hubbard_l(is) + 1
                  k = k + 1
                  f1 (k) = ns (m2,m1,iat,isp)
               enddo
            enddo
            CALL dspev_drv( 'V', 'L', 2 * Hubbard_l(is) + 1,&
                            f1, lambda, vet, ldmx  )
            x_value=alpha_pen(iat)-lambda(2*Hubbard_l(is)+1)
            call stepfn(A_pen(iat,isp),sigma_pen(iat),x_value,g_value,&
                                   step_value)
            do m1 = 1,2*Hubbard_l(is) + 1
               do m2 = 1,2*Hubbard_l(is) + 1
                  forceh = forceh + g_value * dns(m1,m2,iat,isp)           &
                                     * vet(m1,2*Hubbard_l(is)+1)           &
                                     * vet(m2,2*Hubbard_l(is)+1)
               end do
            end do
         endif
      end do
      deallocate ( f1, vet, lambda )
      !
   end subroutine penalty_f
!-----------------------------------------------------------------------
      subroutine stepfn(A,sigma,x_value,g_value,step_value)
!-----------------------------------------------------------------------
!     This subroutine calculates the value of the gaussian and step
!     functions with a given x_value. A and sigma are given in the
!     input file. ... to be used in occupation_constraint...
!
      USE constants, ONLY : pi
      implicit none
      real(kind=8) A, sigma, x_value, g_value, step_value
      real(kind=8) x
      integer i
      step_value=0.0d0
      g_value=0.0d0
!
      do i=1,100000
         x=x_value + (i-100000)/100000.0d0*(x_value + 5.d0*sigma)
!
! Integrate from 5 sigma before the x_value
!
         g_value=A*dexp(-x*x/(2*sigma*sigma))/(sigma*dsqrt(2*pi))
!         write(6,*) 'step', step_value,'g',g_value
!         if (g_value.le.0.0) g_value=0.0
         if ((x_value+5*sigma).ge.0.0d0) then
         step_value=step_value+g_value/100000.0d0*(x_value+5.d0*sigma)
         end if
      end do
      return
      end subroutine stepfn
end module step_penalty

