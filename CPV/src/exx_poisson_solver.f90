SUBROUTINE CGMIC_STDCG(iter, n, eps, fbsscale, coemicf, coeke, rho, pot)
    IMPLICIT NONE
    !------------------------------------------------------------------------
    ! --- pass in variables ---
    !------------------------------------------------------------------------
    integer, intent(out) :: iter
    integer, intent(in)  :: n(3)
    real(8), intent(in)  :: eps
    real(8), intent(in)  :: fbsscale
    real(8), intent(in)  :: coemicf(-3:3,3,3)
    real(8), intent(in)  :: coeke(-3:3,3,3)
    real(8), intent(in)  :: rho(n(1),n(2),n(3))
    real(8), intent(inout) :: pot(n(1),n(2),n(3))
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! --- external function ---
    !------------------------------------------------------------------------
    REAL(8)                   :: DDOT
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! --- local variables ---
    !------------------------------------------------------------------------
    integer              :: itr
    integer              :: i, j, k
    integer              :: nord2
    integer              :: npt
    integer              :: nd(6)
    integer              :: nb(6)
    real(8)              :: nro
    real(8)              :: nr
    real(8)              :: mnr
    real(8)              :: alfa
    real(8), allocatable :: x (:,:,:)
    real(8), allocatable :: r (:,:,:)
    real(8), allocatable :: d0(:,:,:)
    real(8), allocatable :: d1(:,:,:)
    !------------------------------------------------------------------------


    !------------------------------------------------------------------------
    ! --- allocate arrays ---
    !------------------------------------------------------------------------
    nord2=3
    !
    nd(1)=1                                      
    nd(2)=1         
    nd(3)=1         
    nd(4)=n(1)      
    nd(5)=n(2)      
    nd(6)=n(3)      
    !
    nb(1)=1-nord2
    nb(2)=1-nord2
    nb(3)=1-nord2
    nb(4)=n(1)+nord2
    nb(5)=n(2)+nord2
    nb(6)=n(3)+nord2
    !
    npt=(n(1)+2*nord2)*(n(2)+2*nord2)*(n(3)+2*nord2)
    !
    ALLOCATE( x(1-nord2:n(1)+nord2, 1-nord2:n(2)+nord2, 1-nord2:n(3)+nord2))
    ALLOCATE( r(1-nord2:n(1)+nord2, 1-nord2:n(2)+nord2, 1-nord2:n(3)+nord2))
    ALLOCATE(d0(1-nord2:n(1)+nord2, 1-nord2:n(2)+nord2, 1-nord2:n(3)+nord2))
    ALLOCATE(d1(1-nord2:n(1)+nord2, 1-nord2:n(2)+nord2, 1-nord2:n(3)+nord2))
    !
    r = 0.d0; x = 0.d0

    do k = 1, n(3)
        do j = 1, n(2)
            do i = 1, n(1)
                r(i,j,k) = rho(i,j,k)
                x(i,j,k) = pot(i,j,k)
            end do
        end do
    end do
    !
    CALL PADX(nd,nb,coeke,x,d1)
    !
    nro = 0.d0
    do k = nb(3), nb(6)
      do j = nb(2), nb(5)
        do i = nb(1), nb(4)
          r(i,j,k) = r(i,j,k) - 1.0d0 * d1(i,j,k)
          d0(i,j,k) = r(i,j,k) !MCA
          nro = nro + r(i,j,k)*r(i,j,k) 
        end do
      end do
    end do   
    !
    DO itr=0,1000 ! loop over the dimension of the problem                         ! std CG_4 : for i = 1:length(b)
        !
        CALL PADX(nd,nb,coeke,d0,d1)                                               ! std CG_5 : Ap = A * p; % p = d0_d; Ap = d1_d
        !
        alfa = 0.d0
        do k = nb(3), nb(6)
          do j = nb(2), nb(5)
            do i = nb(1), nb(4)
              alfa = alfa + d0(i,j,k) * d1(i,j,k)
            end do
          end do
        end do   
        !
        nr = 0.d0
        !alfa=nro/Ddot(npt,d0_d,1,d1_d,1)                                          ! std CG_6 : alfa = rsold / (p' * Ap);
        do k = nb(3), nb(6)
          do j = nb(2), nb(5)
            do i = nb(1), nb(4)
              x(i,j,k) = x(i,j,k) + nro / alfa * d0(i,j,k) 
              r(i,j,k) = r(i,j,k) - nro / alfa * d1(i,j,k)
              nr = nr + r(i,j,k) * r(i,j,k)
            end do
          end do
        end do   
        !                                 ! std CG_7 : x = x + alfa * p;
        !                                 ! std CG_8 : r = r - alfa * Ap;
        !                                                                                          ! std CG_9 : rsnew = r' * r;
        IF (nr < eps*eps*(1.d0+nro)) EXIT                                        ! std CG_10-12 : if (converge) : break
        do k = nb(3), nb(6)
          do j = nb(2), nb(5)
            do i = nb(1), nb(4)
              d0(i,j,k) = r(i,j,k) + d0(i,j,k) * nr/nro
            end do
          end do
        end do   
        !                                  ! std CG_13 : p = r + (rsnew / rsold) * p;
        !                                    ! std CG_13 : p = r + (rsnew / rsold) * p;
        nro = nr
    END DO
    !------------------------------------------------------------------------

    !!------------------------------------------------------------------------
    !! print error
    !!------------------------------------------------------------------------
    !r_d(nd_d(1):nd_d(4),nd_d(2):nd_d(5),nd_d(3):nd_d(6))=rho(:,:,:)
    !CALL PADX_CUDA(nd_d,nb_d,coeke_d,x_d,d0_d) ! TODO
    !call cublasDaxpy(npt,-1.d0,d0_d,1, r_d, 1)
    !!------------------------------------------------------------------------
    !! WRITE(*,"(A, E15.7, A, I4, A)") "error: ", DSQRT(PDDOT(nd,nb,r,r)), "  in", iter, "  steps"
    !!------------------------------------------------------------------------
    do k = 1, n(3)
        do j = 1, n(2)
            do i = 1, n(1)
                pot(i,j,k) = x(i,j,k)
            end do
        end do
    end do
    !
    iter = itr

    DEALLOCATE(x)
    DEALLOCATE(r)
    !DEALLOCATE(z_d)
    DEALLOCATE(d0)
    DEALLOCATE(d1)
    !
END SUBROUTINE CGMIC_STDCG
!============================================================================
!       OP     PCG    FLOPS      MEM         DEP
! ------------------------------------------------
! [01] DOT              2N       2N
!
! [02] MV              37N       3N
! [03] DOT      *       2N       2N
! [04] DOT              2N       2N        02,03
! [05] AXPY             2N       3N        04
! [06] AXPY             2N       3N        04
! [07] CP               --       2N        06
! [08] SCAL     *        N       2N        07
! [09] FW       *      18N       2N        08 ! TODO
! [10] BW       *      18N       2N        09 ! TODO
! [11] CP       *       --       2N        10
! [12] DOT              2N       2N        09
! [13] AXPY             2N       3N        12

SUBROUTINE PADX(nd,nb,coeke,d,Ad)
    IMPLICIT NONE
    !------------------------------------------------------------------------
    ! --- pass in variables ---
    !------------------------------------------------------------------------
    INTEGER     :: nd(6)
    INTEGER     :: nb(6)
    INTEGER     :: itr, jtr, ktr
    REAL(8)    :: coeke(-3:3,3,3)
    REAL(8)    :: d(nb(1):nb(4), nb(2):nb(5), nb(3):nb(6))
    REAL(8)    :: Ad(nb(1):nb(4), nb(2):nb(5), nb(3):nb(6))
    !------------------------------------------------------------------------

    DO ktr=nd(3),nd(6)
        DO jtr=nd(2),nd(5)
            DO itr=nd(1),nd(4)
                Ad(itr,jtr,ktr)=(coeke(0,1,1) & 
                    +coeke(0,2,2)+coeke(0,3,3))*d(itr,jtr,ktr) & 
                    +coeke(1,1,1)*(d(itr-1,jtr,ktr)+d(itr+1,jtr,ktr)) &
                    +coeke(2,1,1)*(d(itr-2,jtr,ktr)+d(itr+2,jtr,ktr)) &
                    +coeke(3,1,1)*(d(itr-3,jtr,ktr)+d(itr+3,jtr,ktr)) &
                    +coeke(1,2,2)*(d(itr,jtr-1,ktr)+d(itr,jtr+1,ktr)) &
                    +coeke(2,2,2)*(d(itr,jtr-2,ktr)+d(itr,jtr+2,ktr)) &
                    +coeke(3,2,2)*(d(itr,jtr-3,ktr)+d(itr,jtr+3,ktr)) &
                    +coeke(1,3,3)*(d(itr,jtr,ktr-1)+d(itr,jtr,ktr+1)) &
                    +coeke(2,3,3)*(d(itr,jtr,ktr-2)+d(itr,jtr,ktr+2)) &
                    +coeke(3,3,3)*(d(itr,jtr,ktr-3)+d(itr,jtr,ktr+3))
                IF (ABS(coeke(1,1,2)) .GT. 1.0e-6) THEN
                 Ad(itr,jtr,ktr)=Ad(itr,jtr,ktr)+ &
                     coeke(1,1,2)*(d(itr+1,jtr+1,ktr)-d(itr+1,jtr-1,ktr) & 
                     -d(itr-1,jtr+1,ktr)+d(itr-1,jtr-1,ktr)) &
                     +coeke(2,1,2)*(d(itr+2,jtr+2,ktr) &
                     -d(itr+2,jtr-2,ktr)-d(itr-2,jtr+2,ktr) &
                     +d(itr-2,jtr-2,ktr)) &
                     +coeke(3,1,2)*(d(itr+3,jtr+3,ktr) &
                     -d(itr+3,jtr-3,ktr)-d(itr-3,jtr+3,ktr)+d(itr-3,jtr-3,ktr))

                END IF

                IF (ABS(coeke(1,1,3)) .GT. 1.0e-6) THEN
                 Ad(itr,jtr,ktr)=Ad(itr,jtr,ktr)+ &
                     coeke(1,1,3)*(d(itr+1,jtr,ktr+1)- &
                     d(itr+1,jtr,ktr-1)-d(itr-1,jtr,ktr+1) &
                     +d(itr-1,jtr,ktr-1)) &
                     +coeke(2,1,3)*(d(itr+2,jtr,ktr+2) &
                     -d(itr+2,jtr,ktr-2)-d(itr-2,jtr,ktr+2) &
                     +d(itr-2,jtr,ktr-2)) &
                     +coeke(3,1,3)*(d(itr+3,jtr,ktr+3) &
                     -d(itr+3,jtr,ktr-3)-d(itr-3,jtr,ktr+3) &
                     +d(itr-3,jtr,ktr-3))
                END IF

                IF (ABS(coeke(1,2,3)) .GT. 1.0e-6) THEN
                 Ad(itr,jtr,ktr)=Ad(itr,jtr,ktr) &
                     +coeke(1,2,3)*(d(itr,jtr+1,ktr+1) &
                     -d(itr,jtr+1,ktr-1)-d(itr,jtr-1,ktr+1) &
                     +d(itr,jtr-1,ktr-1)) &
                     +coeke(2,2,3)*(d(itr,jtr+2,ktr+2) &
                     -d(itr,jtr+2,ktr-2)-d(itr,jtr-2,ktr+2)+d(itr,jtr-2,ktr-2)) &
                     +coeke(3,2,3)*(d(itr,jtr+3,ktr+3) &
                     -d(itr,jtr+3,ktr-3)-d(itr,jtr-3,ktr+3)+d(itr,jtr-3,ktr-3))
                END IF
            END DO
        END DO
    END DO
END SUBROUTINE PADX
