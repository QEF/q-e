!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

   MODULE orthogonalize_base

      USE kinds
      USE parallel_toolkit, ONLY: matmulp, cmatmulp, pdspev_drv, dspev_drv, &
                                  pzhpev_drv, zhpev_drv

      IMPLICIT NONE

      SAVE

      PRIVATE

      REAL(DP) :: one, zero, two, minus_one, minus_two
      PARAMETER ( one = 1.0d0, zero = 0.0d0, two = 2.0d0, minus_one = -1.0d0 )
      PARAMETER ( minus_two = -2.0d0 )
      COMPLEX(DP) :: cone, czero, mcone
      PARAMETER ( cone = (1.0d0, 0.0d0), czero = (0.0d0, 0.0d0) )
      PARAMETER ( mcone = (-1.0d0, 0.0d0) )
      REAL(DP) :: small = 1.0d-14

#if defined __AIX
       INTEGER, PARAMETER :: nrlx_tune = 128
#else
       INTEGER, PARAMETER :: nrlx_tune = 4
#endif

       INTERFACE sqr_matmul
         MODULE PROCEDURE sqr_dmatmul, sqr_cmatmul
       END INTERFACE

       INTERFACE sigset
         MODULE PROCEDURE rsigset, csigset
       END INTERFACE
       INTERFACE rhoset
         MODULE PROCEDURE rrhoset, crhoset
       END INTERFACE

       INTERFACE diagonalize_rho
         MODULE PROCEDURE diagonalize_rrho, diagonalize_crho
       END INTERFACE

       PUBLIC :: sqr_matmul, sigset, rhoset, diagonalize_rho
       PUBLIC :: backrhoset2, sigrhoset2, backrhoset, sigrhoset
       PUBLIC :: ortho_iterate

     CONTAINS


       SUBROUTINE sqr_dmatmul(transa,transb,a,b,c)
! ...    Multiply square matrices A, B and return the result in C
         USE mp_global, ONLY: nproc
         REAL(DP) :: c(:,:), a(:,:), b(:,:)
         CHARACTER*1 :: transa, transb
         INTEGER :: n
         n = SIZE(c,1)
         IF ( ( nproc > 1 ) .AND. ( n >= nproc ) ) THEN
           CALL matmulp( transa, transb, a, b, c, n )
         ELSE
           CALL DGEMM( transa, transb, n, n, n, one, a(1,1), n, b(1,1), n, zero, c(1,1), n)
         END IF
         RETURN
       END SUBROUTINE sqr_dmatmul

       SUBROUTINE sqr_cmatmul(transa,transb,a,b,c)
! ...    Multiply square matrices A, B and return the result in C
         USE mp_global, ONLY: nproc
         COMPLEX(DP) :: c(:,:), a(:,:), b(:,:)
         CHARACTER*1 transa, transb
         INTEGER :: n
         n = SIZE(c,1)
         IF ((nproc > 1 ).AND. (n >= nproc)) THEN
           CALL cmatmulp(transa,transb,A,B,C,n)
         ELSE
           CALL ZGEMM(transa,transb,n,n,n,cone,a(1,1),n,b(1,1),n,czero,c(1,1),n)
         END IF
         RETURN
       END SUBROUTINE sqr_cmatmul


!.....DIAGONALIZATION OF RHOS

       SUBROUTINE diagonalize_rrho( temp, rhod, s, pwrk)

#if defined __SHMEM
         USE shmem_include
#endif

         USE mp_global, ONLY: nproc, mpime
         USE mp, ONLY: mp_sum
         REAL(DP) :: rhod(:)
         REAL(DP) :: temp(:,:), s(:,:), pwrk(:)

         REAL(DP),   ALLOCATABLE :: aux(:)
         REAL(DP),   ALLOCATABLE :: diag(:,:)
         REAL(DP),   ALLOCATABLE :: vv(:,:)
         INTEGER :: n, nrl

         n = SIZE(temp,1)

         IF ( ( nproc > 1 ) .AND. ( n / nproc ) >= nrlx_tune ) THEN

           nrl = n/nproc
           IF(mpime < MOD(n,nproc)) THEN
             nrl = nrl + 1
           end if
           ALLOCATE( diag(nrl,n), vv(nrl,n) )

           CALL prpack(diag, temp)
           CALL pdspev_drv( 'V', diag, nrl, rhod, vv, nrl, nrl, n, nproc, mpime)
           CALL prunpack(temp, vv)

           DEALLOCATE( diag, vv )

#if defined __SHMEM
           call shmem_barrier_all
           CALL SHMEM_REAL8_SUM_TO_ALL(S,TEMP,N*N,0,0,nproc, pWrk, pSync_sta)
           call shmem_barrier_all
#else
           CALL mp_sum(temp, s)
#endif

         ELSE

           ALLOCATE( aux( n * ( n + 1 ) / 2 ) )
           CALL rpack( aux, temp )
           CALL dspev_drv('V', 'L', n, aux, rhod, s, n)
           DEALLOCATE( aux )

         END IF

         RETURN
       END SUBROUTINE diagonalize_rrho

!  BEGIN manual

      SUBROUTINE diagonalize_crho(a,d,ev)

!  this routine calls the appropriate Lapack routine for diagonalizing a
!  complex Hermitian matrix
!  ----------------------------------------------
!  END manual

         USE mp_global, ONLY: nproc, mpime
         USE mp, ONLY: mp_sum
         IMPLICIT NONE

         REAL(DP)    :: d(:)
         COMPLEX(DP) :: a(:,:), ev(:,:)

         INTEGER :: n, nrl

         COMPLEX(DP), ALLOCATABLE :: aloc(:)
         COMPLEX(DP), ALLOCATABLE :: ap(:,:)
         COMPLEX(DP), ALLOCATABLE :: vp(:,:)

! ...   end of declarations
!  ----------------------------------------------

         n = SIZE(a, 1)

         IF((nproc.EQ.2) .OR. (n.LT.nproc) .OR. (n.LT.256)) THEN

           ALLOCATE(aloc(n*(n+1)/2))

! ...      copy the lower-diagonal part of the matrix according to the
! ...      Lapack packed storage scheme for Hermitian matrices
           CALL zpack(aloc, a)
! ...      call the Lapack routine
           CALL zhpev_drv('V', 'L', n, aloc, d, ev, n)

           DEALLOCATE(aloc)

         ELSE

           nrl = n/nproc
           IF(mpime.LT.MOD(n,nproc)) THEN
             nrl = nrl + 1
           END IF

           ALLOCATE(ap(nrl,n), vp(nrl,n))

           CALL pzpack(ap, a)
           CALL pzhpev_drv( 'V', ap, nrl, d, vp, nrl, nrl, n, nproc, mpime)
           CALL pzunpack(a, vp)
           CALL mp_sum(a, ev)

           DEALLOCATE(ap, vp)

         END IF

         RETURN
       END SUBROUTINE diagonalize_crho


!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE rsigset ( ngw, nb, cp, sig, pwrk)

!     SIG = REAL PART OF ONE-2.0*ADJ(CP)*CP+CP(*,1)*ADJ(CP(*,1))
!     WHERE CP(*,1) IS REAL, AND THEREFORE TRANS() IS USED IN PLACE OF ADJ()
!  ----------------------------------------------
!  END manual

      USE mp_global, ONLY: nproc, group
      USE mp, ONLY: mp_sum
#if defined __SHMEM
      USE shmem_include
#endif

      IMPLICIT NONE

      COMPLEX(DP)        :: CP(:,:)
      REAL(DP)           :: SIG(:,:), PWRK(1)
      INTEGER, INTENT(IN) :: nb, ngw
      INTEGER :: i, ldc, twongw, j, k, lds, n

      ldc = 2 * SIZE( cp, 1 )
      lds =     SIZE( sig, 1 )
      twongw = 2*ngw
      n      = nb

!      WRITE( stdout,*) ' SIGSET 1 ', SUM(sig), SUM(cp)  ! DEBUG
!      DO i = 1, nb
!        DO j = 1, nb
!          DO k = 1, ngw
!            sig(i,j) = - 2.0d0 * ( DBLE(cp(k,i))*DBLE(cp(k,j))+AIMAG(cp(k,i))*AIMAG(cp(k,j)) )
!          END DO
!        END DO
!      END DO

      CALL DGEMM('T','N', n, n, twongw, -2.0d0, cp(1,1), ldc, cp(1,1), ldc, zero, sig(1,1), lds)
      DO i = 1, n
        sig(i,i) = sig(i,i) + one / DBLE(nproc)
      END DO

!      WRITE( stdout,*) ' SIGSET 2 ', SUM(sig)  ! DEBUG


#if defined __SHMEM
      call shmem_barrier_all
      CALL SHMEM_REAL8_SUM_TO_ALL(SIG, SIG, SIZE(sig),0,0, nproc, pWrk, pSync_sta)
      call shmem_barrier_all
#else
      CALL mp_sum( sig, group )
#endif


      RETURN
      END SUBROUTINE rsigset

!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE csigset( ngw, nx, cp, sig, pwrk )

! SIG = REAL PART OF ONE-2.0*ADJ(CP)*CP+CP(*,1)*ADJ(CP(*,1))
!     (WHERE CP(*,1) IS REAL, AND THEREFORE TRANS() IS USED IN
!      PLACE OF ADJ()
!  ----------------------------------------------
!  END manual

      USE mp_global, ONLY: nproc, mpime
      USE mp, ONLY: mp_sum
#if defined __SHMEM
      USE shmem_include
#endif

      IMPLICIT NONE

      INTEGER, INTENT( IN ) :: nx, ngw
      COMPLEX(DP) :: cp(:,:), sig(:,:), pwrk(1)
      INTEGER :: i, j, ldc, lds

      ldc = SIZE( cp, 1 )
      lds = SIZE( sig, 1 )
      
      CALL ZGEMM('C','N',NX,NX,NGW,mcone,CP(1,1),ldc,CP(1,1),ldc,czero,sig(1,1),lds)
      DO i = 1, nx
        sig(i,i) = sig(i,i) + cone / DBLE(nproc)
      END DO
      
#if defined __SHMEM
      call shmem_barrier_all
      CALL SHMEM_REAL8_SUM_TO_ALL(SIG,SIG,2*NX*NX,0,0,nproc,pWrk,pSync_sta)
      call shmem_barrier_all
#else
      CALL mp_sum( sig )
#endif

      DO i=1,nx
        DO j=i,nx
          sig(j,i) = 0.5d0 * ( sig(j,i) + CONJG(sig(i,j)) )
        ENDDO
      ENDDO
      DO i=1,nx
        DO j=1,i-1
          sig(j,i) = CONJG(sig(i,j))
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE csigset

!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE rrhoset ( ngw, nb, c0, cp, rho, tmass, pwrk)

!     RHO   = REAL PART OF 2*ADJ(C0/PMSS)*CP + 
!             C0(*,1)/PMSS*TRANS(CP(*,1)) 
!             (CP(*,1) AND C0(*,1) REAL!)
!
!     TMASS = REAL PART OF 2*ADJ(C0/PMSS)*C0/PMSS + ...
!
!     RHO AND TMASS ARE PLACED IN COMMON /HOPE/
!  ----------------------------------------------
!  END manual
 
      USE mp_global, ONLY: nproc
      USE mp, ONLY: mp_sum
#if defined __SHMEM
      USE shmem_include
#endif

      IMPLICIT  NONE

      COMPLEX(DP) :: CP(:,:), C0(:,:)
      REAL(DP)    :: RHO(:,:), TMASS(:,:)
      REAL(DP)    :: pWrk(1)
      INTEGER, INTENT(IN) :: ngw, nb

      INTEGER ::  i, j, ldc, ldr, tngw, n

      ldc = 2*SIZE( cp, 1 )
      ldr = SIZE( rho, 1 )
      tngw = 2*ngw
      n    = nb

      CALL DGEMM('T','N',n,n,tngw,two,c0(1,1),ldc,cp(1,1),ldc,zero,rho(1,1),ldr)
      CALL DGEMM('T','N',n,n,tngw,two,c0(1,1),ldc,c0(1,1),ldc,zero,tmass(1,1),ldr)

#if defined __SHMEM
      call shmem_barrier_all
      CALL SHMEM_REAL8_SUM_TO_ALL(RHO, RHO, size(rho),0,0, nproc, &
        pWrk, pSync_sta)
      CALL SHMEM_REAL8_SUM_TO_ALL(TMASS,TMASS,size(tmass),0,0,nproc, &
        pWrk,pSync_sta)
      call shmem_barrier_all
#else
      CALL mp_sum( rho )
      CALL mp_sum( tmass )
#endif

      RETURN
      END SUBROUTINE rrhoset

!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE crhoset( ngw, nx, c0, cp, rho, tmass, pwrk )
!
!     RHO   = REAL PART OF 2*ADJ(C0/PMSS)*CP + 
!             C0(*,1)/PMSS*TRANS(CP(*,1)) 
!             (CP(*,1) AND C0(*,1) REAL!)
!
!     TMASS = REAL PART OF 2*ADJ(C0/PMSS)*C0/PMSS + ...
!
!  ----------------------------------------------
!  END manual

      USE mp_global, ONLY: nproc
      USE mp, ONLY: mp_sum
#if defined __SHMEM
      USE shmem_include
#endif

      IMPLICIT  NONE

      INTEGER :: nx, ngw
      COMPLEX(DP) :: cp(:,:), c0(:,:)
      COMPLEX(DP) :: rho(:,:), tmass(:,:)
      COMPLEX(DP) :: pwrk(1)
      INTEGER  :: i, j, ldc, ldr

      ldc = SIZE( c0, 1 )
      ldr = SIZE( rho, 1 )
      CALL ZGEMM('C','N',nx,nx,ngw,cone,c0(1,1),ldc,cp(1,1),ldc,czero,rho(1,1),ldr)
      CALL ZGEMM('C','N',nx,nx,ngw,cone,c0(1,1),ldc,c0(1,1),ldc,czero,tmass(1,1),ldr)

#if defined __SHMEM
      call shmem_barrier_all
      CALL SHMEM_REAL8_SUM_TO_ALL(TMASS,TMASS,2*NX*NX,0,0,nproc,pWrk,pSync_sta)
      CALL SHMEM_REAL8_SUM_TO_ALL(RHO,RHO,2*NX*NX,0,0,nproc,pWrk,pSync_sta)
      call shmem_barrier_all
#else
      CALL mp_sum( RHO )
      CALL mp_sum( TMASS )
#endif

      DO I=1,NX
        DO J=I,NX
          TMASS(J,I) = 0.5d0 * ( TMASS(J,I) + CONJG(TMASS(I,J)) )
        ENDDO
      ENDDO
      DO I=1,NX
        DO J=1,I-1
          TMASS(J,I) = CONJG(TMASS(I,J))
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE crhoset

!----------------------------------------------------------------------!


      SUBROUTINE sigrhoset( ngw, nx, cp, c0, sig, rho, tmass, pmss, emass, gzero )

!***************************************************************
! SIG = REAL PART OF ONE-2.0*ADJ(CP)*CP+CP(*,1)*ADJ(CP(*,1))
!     (WHERE CP(*,1) IS REAL, AND THEREFORE TRANS() IS USED IN
!      PLACE OF ADJ()
!***************************************************************
      USE parallel_types
      USE descriptors_module
      USE mp, ONLY: mp_sum
      USE processors_grid_module, ONLY: get_grid_coor, get_grid_dims

      IMPLICIT NONE

      REAL (DP) :: PMSS(:), EMASS
      TYPE (real_parallel_matrix) ::  sig
      TYPE (real_parallel_matrix) ::  rho
      TYPE (real_parallel_matrix) ::  tmass
      COMPLEX (DP) :: CP(:,:)
      COMPLEX (DP) :: C0(:,:)
      LOGICAL, INTENT(IN) :: gzero
      INTEGER, INTENT(IN) :: ngw, nx

      INTEGER I,J,NPROW, NPCOL, MYROW, MYCOL, CURROW, CURCOL
      INTEGER ME, NRL, NCL, N, NB, CURNR, CURNC, II, JJ,IB,JB
      INTEGER IP,JP,NIB, NJB, JP1, JP2, IOFF1, IOFF2, IOFF3
      INTEGER RSRC, CSRC, npz, mez
      INTEGER INDXL2G, INDXG2L, INDXG2P 
      REAL(DP) DDOT
      REAL(DP), allocatable :: SIGTMP(:)
      real(DP), allocatable :: CTMP(:,:)
      real(DP), allocatable :: ebpmss(:)
      real(DP) sqrtfact
      COMPLEX (DP) :: C0ji
      INTEGER :: ldc

!
!     SUBROUTINE BODY
!
      CALL get_grid_dims(sig%desc%grid, nprow, npcol, npz)
      CALL get_grid_coor(sig%desc%grid, myrow, mycol, mez)
      ldc   = 2 * SIZE(C0,1)
      N     = sig%desc%nx
      NB    = sig%desc%nxblk
      RSRC  = sig%desc%ipexs
      CSRC  = sig%desc%ipeys

      ALLOCATE( ctmp( nx, 2*ngw ) )
      ALLOCATE( ebpmss( ngw ) )
      CALL avrec( ngw, emass, pmss, ebpmss )
      DO I = 1, N
        DO J = 1, NGW
          jp1         = j+j-1
          jp2         = j+j
          c0ji        = c0(j,i)
          c0ji        = c0ji * ebpmss(j)
          c0(j,i)     = c0ji
          ctmp(i,jp1) =  DBLE(c0ji)
          ctmp(i,jp2) = AIMAG(c0ji)
        ENDDO
      ENDDO
      DEALLOCATE( ebpmss )


      sqrtfact = sqrt(0.5_DP)
      IF(gzero) THEN
        DO I = 1, N
          CTMP(I,1) = sqrtfact * CTMP(I,1)
          C0(1,I)  = sqrtfact * C0(1,I)
          CP(1,I)  = sqrtfact * CP(1,I)
        ENDDO
      ENDIF


!     LOOP OVER THE 2D PROCESSORS GRID
!
      DO CURCOL = 0, NPCOL-1
        DO CURROW = 0, NPROW-1

          CURNR = NUMROC(N, NB, CURROW, RSRC, NPROW)
          CURNC = NUMROC(N, NB, CURCOL, CSRC, NPCOL)

          allocate(sigtmp(CURNR*CURNC*3))

!         LOOP OVER THE BLOCKS OWNED BY PE (CURROW,CURCOL)
!
          do jb = 1,N,NB
            jp = INDXG2P( jb, NB, MYCOL, CSRC, NPCOL )
            njb = min(nb,n-jb+1)
            if(jp.eq.curcol) then
              do ib = 1,N,NB
                nib = min(nb,n-ib+1)
                ip = INDXG2P( IB, NB, MYROW, RSRC, NPROW )
                if(ip.eq.currow) then
                  I = INDXG2L( IB, NB, CURROW, RSRC, NPROW )
                  J = INDXG2L( JB, NB, CURCOL, CSRC, NPCOL )
                  IOFF1 = i+(j-1)*CURNR
                  IOFF2 = i+(j-1)*CURNR +   CURNR*CURNC
                  IOFF3 = i+(j-1)*CURNR + 2*CURNR*CURNC
                  ! SIG
                  call DGEMM('T','N',nib,njb,2*ngw,-2.0d0, &
     &                 CP(1,ib),ldc,CP(1,jb),ldc,0.0d0, &
     &                 sigtmp(IOFF1),CURNR)
                  ! RHO
                  call DGEMM('N','N',nib,njb,2*ngw,2.0d0, &
     &                 CTMP(ib,1),n,CP(1,jb),ldc,0.0d0, &
     &                 sigtmp(IOFF2),CURNR)
                  ! TMASS
                  call DGEMM('N','N',nib,njb,2*ngw,2.0d0, &
     &                 CTMP(ib,1),n,C0(1,jb),ldc,0.0d0, &
     &                 sigtmp(IOFF3),CURNR)

                endif
              enddo
            endif
          enddo

          call mp_sum( sigtmp )

          IF( (CURROW.eq.MYROW) .and. (CURCOL.eq.MYCOL) ) THEN 
            DO J=1,CURNC
              DO I=1,CURNR
                SIG%m(I,J)   = sigtmp(i+(j-1)*CURNR)
                RHO%m(I,J)   = sigtmp(i+(j-1)*CURNR+   CURNR*CURNC)
                TMASS%m(I,J) = sigtmp(i+(j-1)*CURNR+ 2*CURNR*CURNC)
              ENDDO
            ENDDO
          ENDIF

          DEALLOCATE(SIGTMP)

        ENDDO
      ENDDO

      deallocate(ctmp)

      DO I=1,N
        II =     INDXG2L( I, NB, MYROW, RSRC, NPROW )
        JJ =     INDXG2L( I, NB, MYCOL, CSRC, NPCOL )
        CURROW = INDXG2P( I, NB, MYROW, RSRC, NPROW )
        CURCOL = INDXG2P( I, NB, MYCOL, CSRC, NPCOL )
        IF( ( CURROW .eq. MYROW ) .and. ( CURCOL .eq. MYCOL ) ) THEN
            SIG%m(II,JJ) =  SIG%m(II,JJ) + 1.0
        END IF
      ENDDO

      IF(gzero) THEN
        DO I=1,N
          C0(1,I) = C0(1,I)/sqrtfact
          CP(1,I) = CP(1,I)/sqrtfact
        ENDDO
      ENDIF


      RETURN
      END SUBROUTINE sigrhoset


      SUBROUTINE sigrhoset2( ngw, nx, CP,C0,SIG,RHO,TMASS,PMSS,EMASS,gzero)

!***************************************************************
! SIG = REAL PART OF ONE-2.0*ADJ(CP)*CP+CP(*,1)*ADJ(CP(*,1))
!     (WHERE CP(*,1) IS REAL, AND THEREFORE TRANS() IS USED IN
!      PLACE OF ADJ()
!***************************************************************
      USE parallel_types
      USE descriptors_module
      USE mp, ONLY: mp_sum
      USE processors_grid_module, ONLY: get_grid_coor, get_grid_dims

      IMPLICIT NONE

      REAL (DP) :: PMSS(:), EMASS
      TYPE (real_parallel_matrix) ::  sig
      TYPE (real_parallel_matrix) ::  rho
      TYPE (real_parallel_matrix) ::  tmass
      COMPLEX (DP) :: CP(:,:)
      COMPLEX (DP) :: C0(:,:)
      LOGICAL, INTENT(IN) :: gzero
      INTEGER, INTENT(IN) :: ngw, nx

      INTEGER :: I,J, NPROW, NPCOL, MYROW, MYCOL
      INTEGER :: NRL, N, II, JJ
      INTEGER :: ip, ldc
      INTEGER :: nngw, npz, mez, nproc, mpime
      INTEGER :: nrl_ip, nrlx
      REAL(DP), ALLOCATABLE :: RTMP(:,:,:)
      REAL(DP), ALLOCATABLE :: ebpmss(:)
      REAL(DP) :: sqrtfact

!
!     SUBROUTINE BODY
!
      CALL get_grid_dims( sig%desc%grid, nprow, npcol, npz )
      CALL get_grid_coor( sig%desc%grid, myrow, mycol, mez )

      ldc   = 2 * SIZE(C0, 1)
      nngw  = 2 * ngw
      n     = nx
      nproc = nprow
      mpime = myrow
      nrl  = n/nproc
      nrlx = n/nproc + 1
      IF( mpime < MOD(n,nproc) ) THEN
        nrl = nrl + 1
      END IF

      IF( npcol > 1 .OR. npz > 1 ) THEN
        CALL errore(' sigrhoset2 ',' wrong grid dimension ', npcol)
      END IF
      IF( mycol > 0 .OR. mez > 0 ) THEN
        CALL errore(' sigrhoset2 ',' wrong grid coordinates ', mycol)
      END IF
      IF( nrl /= SIZE(sig%m,1) .OR. n /= SIZE(sig%m,2)) THEN
        CALL errore(' sigrhoset2 ',' wrong sizes for matrix SIG ', nrl)
      END IF
      IF( nrl /= SIZE(rho%m,1) .OR. n /= SIZE(rho%m,2)) THEN
        CALL errore(' sigrhoset2 ',' wrong sizes for matrix RHO ', nrl)
      END IF
      IF( nrl /= SIZE(tmass%m,1) .OR. n /= SIZE(tmass%m,2)) THEN
        CALL errore(' sigrhoset2 ',' wrong sizes for matrix TMASS ', nrl)
      END IF
      IF( n > SIZE(C0,2) ) THEN
        CALL errore(' sigrhoset2 ',' wrong number of states ', n)
      END IF

      ALLOCATE( ebpmss( ngw ) )
      CALL avrec( ngw, emass, pmss, ebpmss )

      sqrtfact = sqrt(0.5d0)
      IF( gzero ) THEN
        DO i = 1, n
          C0(1,i)  = sqrtfact * C0(1,i) 
          CP(1,i)  = sqrtfact * CP(1,i)
        ENDDO
      ENDIF
      DO j = 1, n
        DO i = 1, ngw
          C0(i,j) = C0(i,j) * ebpmss(i)
        END DO
      END DO

      DEALLOCATE(ebpmss)

!
!     LOOP OVER THE 1D PROCESSORS GRID
!

      IF( nrlx > nrlx_tune ) THEN

        DO ip = 0, nproc-1

          nrl_ip = n/nproc
          IF( ip < MOD(n,nproc) ) THEN
            nrl_ip = nrl_ip + 1
          end if

          allocate( rtmp( nrl_ip, n, 1 ) )

          ! SIG
          CALL DGEMM('T', 'N', nrl_ip, n, nngw, -2.0d0, &
     &      cp(1,ip+1), ldc * nproc, CP(1,1), ldc, 0.0d0, rtmp(1,1,1), nrl_ip)

          call mp_sum( rtmp(:,:,1), sig%m(:,:), root=ip  )

          ! RHO
          CALL DGEMM('T', 'N', nrl_ip, n, nngw, 2.0d0, &
     &      c0(1,ip+1), ldc * nproc, CP(1,1), ldc, 0.0d0, rtmp(1,1,1), nrl_ip)

          call mp_sum( rtmp(:,:,1), rho%m(:,:), root=ip  )

          ! TMASS
          CALL DGEMM('T', 'N', nrl_ip, n, nngw, 2.0d0, &
     &      c0(1,ip+1), ldc * nproc, C0(1,1), ldc, 0.0d0, rtmp(1,1,1), nrl_ip)

          call mp_sum( rtmp(:,:,1), tmass%m(:,:), root=ip  )

          DEALLOCATE(rtmp)

        ENDDO

      ELSE

        ALLOCATE( rtmp( n, n, 1 ) )

        ! SIG
        CALL DGEMM('T', 'N', n, n, nngw, -2.0d0, &
     &      cp(1,1), ldc, CP(1,1), ldc, 0.0d0, rtmp(1,1,1), n)
        call mp_sum( rtmp(:,:,1) )
        DO j = 1, n
          ii = mpime + 1
          DO i = 1, nrl
            sig%m(i,j) = rtmp(ii,j,1)
            ii = ii + nproc
          END DO
        END DO

        ! RHO
        CALL DGEMM('T', 'N', n, n, nngw, 2.0d0, &
     &      c0(1,1), ldc, CP(1,1), ldc, 0.0d0, rtmp(1,1,1), n)
        call mp_sum( rtmp(:,:,1)  )
        DO j = 1, n
          ii = mpime + 1
          DO i = 1, nrl
            rho%m(i,j) = rtmp(ii,j,1)
            ii = ii + nproc
          END DO
        END DO

        ! TMASS
        CALL DGEMM('T', 'N', n, n, nngw, 2.0d0, &
     &      c0(1,1), ldc, C0(1,1), ldc, 0.0d0, rtmp(1,1,1), n)
        call mp_sum( rtmp(:,:,1) )
        DO j = 1, n
          ii = mpime + 1
          DO i = 1, nrl
            tmass%m(i,j) = rtmp(ii,j,1)
            ii = ii + nproc
          END DO
        END DO

        DEALLOCATE(rtmp)

      END IF

! ... Add "1" to the diagonal elements
      ii = mpime + 1
      DO i = 1, nrl
        SIG%m(i,ii) =  SIG%m(i,ii) + 1.0d0
        ii = ii + nproc
      ENDDO

! ... Restore the coefficients of the G=0 reciprocal vector
      IF(gzero) THEN
        DO I=1,N
          C0(1,I) = C0(1,I)/sqrtfact
          CP(1,I) = CP(1,I)/sqrtfact
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE sigrhoset2

!
!----------------------------------------------------------------------
!

      SUBROUTINE backrhoset( ngw, nx, cp, c0, rho, pmss, emass )

      USE parallel_types
      USE descriptors_module
      USE mp, ONLY: mp_sum
      USE processors_grid_module, ONLY: get_grid_coor, get_grid_dims

      IMPLICIT NONE
      TYPE (real_parallel_matrix) :: RHO
      COMPLEX (DP) :: CP(:,:)
      COMPLEX (DP) :: C0(:,:)
      REAL (DP)  :: PMSS(:), EMASS
      INTEGER, INTENT(IN) :: ngw, nx

      INTEGER I,J,NPROW, NPCOL, MYROW, MYCOL, CURROW, CURCOL
      INTEGER ME, NRL, NCL, N, NB, CURNR, CURNC, II, JJ,IB,JB
      INTEGER IP,JP,NIB, NJB, JP1, JP2, IOFF1, IOFF2, IOFF3
      INTEGER RSRC, CSRC, npz, mez, ldc
      INTEGER INDXL2G, INDXG2L, INDXG2P 
      REAL (DP)  :: DDOT
      REAL (DP)  :: FACT,ONE_BY_EMASS
      REAL (DP), allocatable :: SIGTMP(:)

!
!     SUBROUTINE BODY
!

      CALL get_grid_dims(rho%desc%grid, nprow, npcol, npz)
      CALL get_grid_coor(rho%desc%grid, myrow, mycol, mez)
      ldc   = 2 * SIZE( C0, 1 )
      N     = nx
      NB    = rho%desc%nxblk
      RSRC  = rho%desc%ipexs
      CSRC  = rho%desc%ipeys


!     LOOP OVER THE 2D PROCESSORS GRID
!
      DO CURCOL = 0, NPCOL-1
        DO CURROW = 0, NPROW-1

          CURNR = NUMROC(N, NB, CURROW, RSRC, NPROW)
          CURNC = NUMROC(N, NB, CURCOL, CSRC, NPCOL)

          allocate(sigtmp(CURNR*CURNC))
          sigtmp = 0.0d0
          IF( (CURROW.eq.MYROW) .and. (CURCOL.eq.MYCOL) ) THEN 
            DO J=1,CURNC
              DO I=1,CURNR
                sigtmp(i+(j-1)*CURNR) = RHO%m(I,J)
              ENDDO
            ENDDO
          ENDIF
          call mp_sum( sigtmp )

!         LOOP OVER THE BLOCKS OWNED BY PE (CURROW,CURCOL)
!
          do jb = 1,N,NB
            jp = INDXG2P( jb, NB, MYCOL, CSRC, NPCOL )
            njb = min(nb,n-jb+1)
            if(jp.eq.curcol) then
              do ib = 1,N,NB
                nib = min(nb,n-ib+1)
                ip = INDXG2P( IB, NB, MYROW, RSRC, NPROW )
                if(ip.eq.currow) then
                  I = INDXG2L( IB, NB, CURROW, RSRC, NPROW )
                  J = INDXG2L( JB, NB, CURCOL, CSRC, NPCOL )
                  IOFF1 = i+(j-1)*CURNR
                  CALL DGEMM('N','N',2*NGW,nib,njb,1.0D0, &
     &                 C0(1,ib),ldc,sigtmp(IOFF1),CURNR,1.0D0, &
     &                 CP(1,jb),ldc )

                endif
              enddo
            endif
          enddo

          DEALLOCATE(SIGTMP)

        ENDDO
      ENDDO

!
!.....RESTORE C0
!
      ALLOCATE(SIGTMP(NGW))
      ONE_BY_EMASS = 1.0d0 / EMASS
      DO J=1,NGW
        SIGTMP(J) = PMSS(J) * ONE_BY_EMASS
      END DO
      DO I=1,N
        DO J=1,NGW
          C0(J,I) = C0(J,I) * SIGTMP(J)
        END DO
      END DO
      DEALLOCATE(SIGTMP)


      RETURN
      END SUBROUTINE backrhoset


      SUBROUTINE backrhoset2( ngw, nx, cp, c0, rho, pmss, emass )

      USE parallel_types
      USE descriptors_module
      USE mp, ONLY: mp_sum, mp_bcast
      USE processors_grid_module, ONLY: get_grid_coor, get_grid_dims

      IMPLICIT NONE
      TYPE (real_parallel_matrix) :: RHO
      COMPLEX (DP) :: CP(:,:)
      COMPLEX (DP) :: C0(:,:)
      REAL (DP)  :: PMSS(:), EMASS
      INTEGER, INTENT(IN) :: ngw, nx

      INTEGER I,J,NPROW, NPCOL, MYROW, MYCOL
      INTEGER NRL, nrl_ip, n, ii, jj 
      INTEGER ip, nngw, nrlx
      INTEGER npz, mez, mpime, nproc, ldc

      REAL (DP)  :: DDOT
      REAL (DP)  :: FACT,ONE_BY_EMASS
      REAL (DP), allocatable :: SIGTMP(:)
      REAL (DP), allocatable :: rtmp(:,:)
      COMPLEX (DP), allocatable :: cc(:,:)

!
!     SUBROUTINE BODY
!

      CALL get_grid_dims(rho%desc%grid, nprow, npcol, npz)
      CALL get_grid_coor(rho%desc%grid, myrow, mycol, mez)

      ldc   = 2 * SIZE( C0, 1 )
      nngw  = 2 * ngw
      n     = nx
      nproc = nprow
      mpime = myrow
      nrlx  = n/nproc+1
      nrl   = n/nproc
      IF( mpime < MOD(n,nproc) ) THEN
        nrl = nrl + 1
      END IF

      IF( npcol > 1 .OR. npz > 1 ) THEN
        CALL errore(' backrhoset2 ',' wrong grid dimension ', npcol)
      END IF
      IF( mycol > 0 .OR. mez > 0 ) THEN
        CALL errore(' backrhoset2 ',' wrong grid coordinates ', mycol)
      END IF
      IF( nrl /= SIZE(rho%m,1) .OR. n /= SIZE(rho%m,2)) THEN
        CALL errore(' backrhoset2 ',' wrong sizes for matrix RHO ', nrl)
      END IF
      IF( n > SIZE(C0,2) ) THEN
        CALL errore(' backrhoset2 ',' wrong number of states ', n)
      END IF

      IF( nrlx > nrlx_tune ) THEN

        DO ip = 0, nproc-1

          nrl_ip = n/nproc
          IF( ip < MOD(n,nproc) ) THEN
            nrl_ip = nrl_ip + 1
          end if

          allocate( rtmp( nrl_ip, n ) )
  
          IF ( mpime == ip ) THEN
            rtmp = rho%m( 1:nrl_ip, 1:n )
          END IF
          call mp_bcast( rtmp, ip )

          CALL DGEMM('N', 'N', nngw, n, nrl_ip, 1.0D0, &
     &      c0(1,ip+1), ldc * nproc, rtmp(1,1), nrl_ip, 1.0D0, cp(1,1), ldc )

          DEALLOCATE(rtmp)

        ENDDO

      ELSE

        allocate( rtmp( n, n ) )
        rtmp = 0.0d0
        DO j = 1, n
          ii = mpime + 1
          DO i = 1, nrl
            rtmp(ii,j) = rho%m(i,j)
            ii = ii + nproc
          END DO
        END DO
        call mp_sum( rtmp(:,:) )

        CALL DGEMM('N', 'N', nngw, n, n, 1.0D0, &
     &    c0(1,1), ldc, rtmp(1,1), n, 1.0D0, cp(1,1), ldc )

        deallocate( rtmp )

      END IF

!
!.....RESTORE C0
!
      ALLOCATE(SIGTMP(NGW))
      ONE_BY_EMASS = 1.0d0 / EMASS
      DO J=1,NGW
        SIGTMP(J) = PMSS(J) * ONE_BY_EMASS
      END DO
      DO I=1,N
        DO J=1,NGW
          C0(J,I) = C0(J,I) * SIGTMP(J)
        END DO
      END DO
      DEALLOCATE(SIGTMP)

      RETURN
      END SUBROUTINE backrhoset2


!=----------------------------------------------------------------------------=!


      SUBROUTINE prpack( ap, a)
        USE mp_global, ONLY: mpime, nproc
        REAL(DP), INTENT(IN) :: a(:,:)
        REAL(DP), INTENT(OUT) :: ap(:,:)
        INTEGER :: i, j, jl
        DO i = 1, SIZE( ap, 2)
           j = mpime + 1
           DO jl = 1, SIZE( ap, 1)
             ap(jl,i) = a(j,i)
             j = j + nproc
           END DO
        END DO
        RETURN
      END SUBROUTINE prpack

      SUBROUTINE pzpack( ap, a)
        USE mp_global, ONLY: mpime, nproc
        COMPLEX(DP), INTENT(IN) :: a(:,:)
        COMPLEX(DP), INTENT(OUT) :: ap(:,:)
        INTEGER :: i, j, jl
        DO i = 1, SIZE( ap, 2)
           j = mpime + 1
           DO jl = 1, SIZE( ap, 1)
             ap(jl,i) = a(j,i)
             j = j + nproc
           END DO
        END DO
        RETURN
      END SUBROUTINE pzpack

      SUBROUTINE prunpack( a, ap)
        USE mp_global, ONLY: mpime, nproc
        REAL(DP), INTENT(IN) :: ap(:,:)
        REAL(DP), INTENT(OUT) :: a(:,:)
        INTEGER :: i, j, jl
        DO i = 1, SIZE(a, 2)
          DO j = 1, SIZE(a, 1)
            a(j,i) = zero
          END DO
          j = mpime + 1
          DO jl = 1, SIZE(ap, 1)
            a(j,i) = ap(jl,i)
            j = j + nproc
          END DO
        END DO
        RETURN
      END SUBROUTINE prunpack

      SUBROUTINE pzunpack( a, ap)
        USE mp_global, ONLY: mpime, nproc
        COMPLEX(DP), INTENT(IN) :: ap(:,:)
        COMPLEX(DP), INTENT(OUT) :: a(:,:)
        INTEGER :: i, j, jl
        DO i = 1, SIZE(a, 2)
          DO j = 1, SIZE(a, 1)
            a(j,i) = zero
          END DO
          j = mpime + 1
          DO jl = 1, SIZE(ap, 1)
            a(j,i) = ap(jl,i)
            j = j + nproc
          END DO
        END DO
        RETURN
      END SUBROUTINE pzunpack

      SUBROUTINE rpack( ap, a)
        REAL(DP), INTENT(IN) :: a(:,:)
        REAL(DP), INTENT(OUT) :: ap(:)
        INTEGER :: i, j, k
        K = 0
        DO J = 1, SIZE(a, 2)
          DO I = J, SIZE(a, 1)
            K = K + 1
            ap( k ) = a( i, j )
          END DO
        END DO
        RETURN
      END SUBROUTINE rpack

      SUBROUTINE zpack( ap, a)
        COMPLEX(DP), INTENT(IN) :: a(:,:)
        COMPLEX(DP), INTENT(OUT) :: ap(:)
        INTEGER :: i, j, k
        K=0
        DO J = 1, SIZE(a, 2)
          DO I = J, SIZE(a, 1)
            K = K + 1
            ap(k) = a(i,j)
          END DO
        END DO
        RETURN
      END SUBROUTINE zpack

!=----------------------------------------------------------------------------=!


   SUBROUTINE ortho_iterate( u, diag, xloc, sig, rhor, rhos, tau, nx, nss, max, eps )

      USE kinds,     ONLY: DP
      USE io_global, ONLY: stdout

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: max
      INTEGER, INTENT(IN) :: nx, nss
      REAL(DP), INTENT(IN) :: eps
      REAL(DP) :: u( nx, nx )
      REAL(DP) :: diag( nx )
      REAL(DP) :: xloc( nx, nx )
      REAL(DP) :: rhor( nx, nx )
      REAL(DP) :: rhos( nx, nx )
      REAL(DP) :: tau( nx, nx )
      REAL(DP) :: sig( nx, nx )

      INTEGER :: iter, i, j
      REAL(DP), ALLOCATABLE :: tmp1(:,:), tmp2(:,:), dd(:,:)
      REAL(DP), ALLOCATABLE :: con(:,:), x1(:,:)
      REAL(DP) :: diff

      ALLOCATE( tmp1(nx,nx), tmp2(nx,nx), dd(nx,nx), x1(nx,nx), con(nx,nx) )


      DO iter = 1, max
         !
         !       the following 4 MXMA-calls do the following matrix
         !       multiplications:
         !                       tmp1 = x0*rhor    (1st call)
         !                       dd   = x0*tau*x0  (2nd and 3rd call)
         !                       tmp2 = x0*rhos    (4th call)
         !
         CALL MXMA( xloc,1,nx,rhor,1,nx,tmp1,1,nx,nss,nss,nss)
         CALL MXMA( tau ,1,nx,xloc,1,nx,tmp2,1,nx,nss,nss,nss)
         CALL MXMA( xloc,1,nx,tmp2,1,nx,  dd,1,nx,nss,nss,nss)
         CALL MXMA( xloc,1,nx,rhos,1,nx,tmp2,1,nx,nss,nss,nss)
         !
         DO i=1,nss
            DO j=1,nss
               x1(i,j) = sig(i,j)-tmp1(i,j)-tmp1(j,i)-dd(i,j)
               con(i,j)= x1(i,j)-tmp2(i,j)-tmp2(j,i)
            END DO
         END DO
         !
         !         x1      = sig      -x0*rho    -x0*rho^t  -x0*tau*x0
         !
         diff = 0.d0
         DO i=1,nss
            DO j=1,nss
               IF(ABS(con(i,j)).GT.diff) diff=ABS(con(i,j))
            END DO
         END DO

         IF( diff <= eps ) go to 20

         !
         !     the following two MXMA-calls do:
         !                       tmp1 = x1*u
         !                       tmp2 = ut*x1*u
         !
         CALL MXMA(x1,1,nx,   u,1,nx,tmp1,1,nx,nss,nss,nss)
         CALL MXMA(u ,nx,1,tmp1,1,nx,tmp2,1,nx,nss,nss,nss)
         !
         !       g=ut*x1*u/d  (g is stored in tmp1)
         !
         DO i=1,nss
            DO j=1,nss
               tmp1(i,j)=tmp2(i,j)/(diag(i)+diag(j))
            END DO
         END DO
         !
         !       the following two MXMA-calls do:
         !                       tmp2 = g*ut
         !                       x0 = u*g*ut
         !
         CALL MXMA(tmp1,1,nx,  u,nx,1,tmp2,1,nx,nss,nss,nss)
         CALL MXMA(   u,1,nx,tmp2,1,nx,xloc,1,nx,nss,nss,nss)
         !
      END DO

      WRITE( stdout,*) ' diff= ',diff,' iter= ',iter
      CALL errore('ortho','max number of iterations exceeded',iter)

20    CONTINUE

      DEALLOCATE( tmp1, tmp2, dd, x1, con )

      RETURN
   END SUBROUTINE ortho_iterate


   END MODULE orthogonalize_base
