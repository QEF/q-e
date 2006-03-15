!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
#define CI   ( 0.D0, 1.D0 )
!
! ... written by Manu Sharma ( 2001-2005 )
!
!----------------------------------------------------------------------------
SUBROUTINE wf( clwf, c, bec, eigr, eigrb, taub, irb, &
               b1, b2, b3, Uall, becdr, what1, wfc, jw, ibrav )
  !----------------------------------------------------------------------------
  !
  ! ... this routine calculates overlap matrices
  !
  ! ... routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
  !
  USE kinds,                    ONLY : DP
  USE constants,                ONLY : pi, tpi
  USE ions_base,                ONLY : nsp, na, nax, nat
  USE cvan,                     ONLY : nvb, ish
  USE cell_base,                ONLY : omega, a1, a2, a3, alat, h, ainv
  USE electrons_base,           ONLY : nbspx, nbsp, nupdwn, iupdwn, nspin
  USE gvecb,                    ONLY : npb, nmb, ngb
  USE gvecw,                    ONLY : ngw
  USE reciprocal_vectors,       ONLY : gstart
  USE smooth_grid_dimensions,   ONLY : nnrsx
  USE control_flags,            ONLY : iprsta
  USE qgb_mod,                  ONLY : qgb
  USE wannier_base,             ONLY : wfg, nw, weight, indexplus, indexplusz, &
                                       indexminus, indexminusz, tag, tagp,     &
                                       expo, wfsd
  USE grid_dimensions,          ONLY : nr1, nr2, nr3
  USE smallbox_grid_dimensions, ONLY : nnrbx, nr1b, nr2b, nr3b, &
                                       nr1bx, nr2bx, nr3bx
  USE uspp_param,               ONLY : nh, nhm
  USE uspp,                     ONLY : nkb
  USE io_global,                ONLY : ionode, stdout
  USE mp,                       ONLY : mp_barrier
  USE mp_global,                ONLY : nproc, mpime
  USE fft_module,               ONLY : invfft
  USE fft_base,                 ONLY : dfftp
  USE parallel_include
  !
  IMPLICIT NONE
  !
  INTEGER,           INTENT(IN)    :: irb(3,nat), jw, ibrav
  REAL(DP),    INTENT(INOUT) :: bec(nkb,nbsp), becdr(nkb,nbsp,3)
  REAL(DP),    INTENT(IN)    :: b1(3), b2(3), b3(3), taub(3,nax)
  COMPLEX(DP), INTENT(INOUT) :: c(ngw,nbspx)
  COMPLEX(DP), INTENT(IN)    :: eigr(ngw,nat), eigrb(ngb,nat)
  REAL(DP),    INTENT(INOUT) :: Uall(nbsp,nbsp)
  !
  REAL(DP)    :: becwf(nkb,nbsp), temp3(nkb,nbsp)
  COMPLEX(DP) :: cwf(ngw,nbspx), bec2(nbsp), bec3(nbsp), bec2up(nupdwn(1))
  COMPLEX(DP) :: bec2dw(nupdwn(2)), bec3up(nupdwn(1)), bec3dw(nupdwn(2))
  !
  COMPLEX(DP), ALLOCATABLE :: c_m(:,:), c_p(:,:), c_psp(:,:)
  COMPLEX(DP), ALLOCATABLE :: c_msp(:,:)
  INTEGER,           ALLOCATABLE :: tagz(:)
  REAL(DP),    ALLOCATABLE :: Uspin(:,:)
  COMPLEX(DP), ALLOCATABLE :: X(:,:), X1(:,:), Xsp(:,:), X2(:,:), X3(:,:)
  COMPLEX(DP), ALLOCATABLE :: O(:,:,:), Ospin(:,:,:), Oa(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: qv(:)
  REAL(DP),    ALLOCATABLE :: gr(:,:), mt(:), W(:,:), wr(:)
  INTEGER,           ALLOCATABLE :: f3(:)
  !
  LOGICAL           :: what1
  INTEGER           :: inl, jnl, iss, isa, is, ia, ijv, i, j, k, l, ig, &
                       ierr, ti, tj, tk, iv, jv, inw, iqv, ibig1, ibig2, &
                       ibig3, ir1, ir2, ir3, ir, clwf, m,  &
                       ib, jb, total, nstat, jj, ngpww, irb3
  REAL(DP)    :: t1, t2, t3, taup(3)
  COMPLEX(DP) :: qvt
  REAL (DP)   :: temp_vec(3)
  INTEGER           :: adjust,ini, ierr1,nnn, me
  COMPLEX(DP) :: U2(nbsp,nbsp)
  INTEGER           :: igx, igy, igz
  REAL(DP)    :: wfcx, wfcy, wfcz
  REAL(DP)    :: wfc(3,nbsp)
  REAL(DP)    :: te(6)
  
  COMPLEX(DP), EXTERNAL :: boxdotgridcplx
  !
#if defined (__PARA)
  !
  INTEGER :: proc, ntot, ncol, mc, ngpwpp(nproc)
  INTEGER :: ncol1,nz1, nz_1 
  INTEGER :: nmin(3), nmax(3), n1,n2,nzx,nz,nz_
  INTEGER :: nmin1(3), nmax1(3)
  INTEGER :: root
  INTEGER :: rdispls(nproc),  recvcount(nproc)
  INTEGER :: rdispls1(nproc), recvcount1(nproc)
  INTEGER :: rdispls2(nproc), recvcount2(nproc)
  INTEGER :: sendcount(nproc),  sdispls(nproc)
  INTEGER :: sendcount1(nproc), sdispls1(nproc)
  INTEGER :: sendcount2(nproc), sdispls2(nproc)
  !
  COMPLEX(DP), ALLOCATABLE :: psitot(:,:), psitot_pl(:,:)
  COMPLEX(DP), ALLOCATABLE :: psitot1(:), psitot_p(:)
  COMPLEX(DP), ALLOCATABLE :: psitot_mi(:,:)
  COMPLEX(DP), ALLOCATABLE :: psitot_m(:)
  INTEGER,           ALLOCATABLE :: ns(:)
  !
#endif
  !
  ! ... set up the weights and the G vectors for wannie-function calculation
  !
  me = mpime + 1
  !
  te = 0.D0
  !
  SELECT CASE( ibrav )
  CASE( 0 )
     !
     ! ... free cell used for cpr
     !
     nw = 6
     !
     ALLOCATE( wfg( nw, 3 ), weight( nw ) )
     ALLOCATE( tagz( nw ) )
     !
     tagz(:) = 1
     tagz(3) = 0
     !
     ! ... assigns weights for the triclinic/cpr case
     !
     CALL tric_wts( b1, b2, b3, alat, weight )
     !
     wfg(:,:) = 0
     wfg(1,1) = 1
     wfg(2,2) = 1
     wfg(3,3) = 1
     wfg(4,1) = 1
     wfg(4,2) = 1
     wfg(5,2) = 1
     wfg(5,3) = 1
     wfg(6,1) = 1
     wfg(6,3) = 1
     !
  CASE( 1 )
     !
     ! ... cubic P [sc]
     !
     nw = 3
     !
     ALLOCATE( wfg( nw, 3 ), weight( nw ) )
     ALLOCATE( tagz( nw ) )
     !
     weight = 1.D0
     !
     wfg(:,:) = 0
     wfg(1,1) = 1
     wfg(2,2) = 1
     wfg(3,3) = 1
     !
     tagz(:) = 1
     tagz(3) = 0
     !
  CASE( 2 )
     !
     ! ... cubic F [fcc]
     !
     nw = 4
     !
     ALLOCATE( wfg( nw, 3 ), weight( nw ) )
     ALLOCATE( tagz( nw ) )
     !
     weight = 0.25D0
     !
     wfg(:,:) =  0     
     wfg(1,1) =  1
     wfg(2,2) =  1
     wfg(3,3) =  1
     wfg(4,:) = -1
     !
     tagz(:) = 1
     tagz(3) = 0
     !
  CASE( 3 )
     !
     ! ... cubic I [bcc]
     !
     nw=6
     !
     ALLOCATE( wfg( nw, 3 ), weight( nw ) )
     ALLOCATE( tagz( nw ) )
     !
     weight = 0.25D0
     !
     tagz(:) = 1
     tagz(3) = 0
     !
     wfg(:,:) =  0
     wfg(1,1) =  1
     wfg(2,2) =  1
     wfg(3,3) =  1
     wfg(4,1) =  1
     wfg(4,2) =  1
     wfg(5,2) =  1
     wfg(5,3) =  1
     wfg(6,1) = -1
     wfg(6,3) =  1
     !
  CASE( 4 )
     !
     ! ... hexagonal and trigonal P
     !
     nw = 4
     !
     ALLOCATE( wfg( nw, 3 ), weight( nw ) )
     ALLOCATE( tagz( nw ) )
     !
     tagz(:) = 1
     tagz(3) = 0
     !
     weight(1) = 0.5D0
     weight(2) = 0.5D0
     weight(3) = 1.D0 / b3(3)**2
     weight(4) = 0.5D0
     !
     wfg(:,:) =  0
     wfg(1,1) =  1
     wfg(2,2) =  1
     wfg(3,3) =  1
     wfg(4,1) =  1
     wfg(4,2) = -1
     !
  CASE( 5 )
     !
     ! ... trigonal R
     !
     nw = 6
     !
     ALLOCATE( wfg( nw, 3 ), weight( nw ) )
     ALLOCATE( tagz( nw ) )
     !
     tagz(:) = 1
     tagz(3) = 0
     !
     t1 = 1.D0 - 1.D0 / ( b2(2)**2 * 1.5D0 )
     !
     weight(1) =  1.D0
     weight(2) =  1.D0 + 2.D0 * t1
     weight(3) =  1.D0
     weight(4) =  t1
     weight(5) = -t1
     weight(6) = -t1
     !
     wfg(:,:) =  0
     wfg(1,1) =  1
     wfg(2,2) =  1
     wfg(3,3) =  1
     wfg(4,1) =  1
     wfg(4,3) =  1
     wfg(5,2) = -1
     wfg(5,3) =  1
     wfg(6,1) =  1
     wfg(6,2) = -1
     !
  CASE( 6 )
     !
     ! ... tetragonal P[st]
     !
     nw=3
     !
     ALLOCATE( wfg( nw , 3 ), weight( nw ) )
     ALLOCATE( tagz(nw) )
     !
     tagz(:) = 1
     tagz(3) = 0
     !
     weight(1) = 1.D0
     weight(2) = 1.D0
     weight(3) = 1.D0 / b3(3)**2
     !
     wfg(:,:) = 0
     wfg(1,1) = 1
     wfg(2,2) = 1
     wfg(3,3) = 1
     !
  CASE( 7 )
     !
     ! ... tetragonal I [bct]
     !
     nw=6
     !
     ALLOCATE( wfg( nw, 3 ), weight( nw ) )
     ALLOCATE( tagz( nw ) )
     !
     tagz(:) = 1
     tagz(3) = 0  
     !
     t1 = 0.25D0 / b2(3)**2
     !
     weight(1) = 0.5D0 - t1
     weight(2) = t1
     weight(3) = t1
     weight(4) = t1
     weight(5) = weight(1)
     weight(6) = t1
     !
     wfg(:,:) =  0
     wfg(1,1) =  1
     wfg(2,2) =  1
     wfg(3,3) =  1
     wfg(4,1) =  1
     wfg(4,3) =  1
     wfg(5,2) =  1
     wfg(5,3) = -1
     wfg(6,1) =  1
     wfg(6,2) =  1
     !
  CASE( 8 )
     !
     ! ... orthorhombic P
     !
     nw = 3
     !
     ALLOCATE( wfg( nw, 3 ), weight( nw ) )
     ALLOCATE( tagz( nw ) )
     !
     tagz(:) = 1
     tagz(3) = 0
     !
     weight(1) = 1.D0
     weight(2) = 1.D0 / b2(2)**2
     weight(3) = 1.D0 / b3(3)**2
     !
     wfg(:,:) = 0
     wfg(1,1) = 1
     wfg(2,2) = 1
     wfg(3,3) = 1
     !
  CASE( 9 )
     !
     ! ... one face centered orthorhombic C
     !
     IF ( b1(2) == 1 ) THEN
        !
        nw = 3
        !
        ALLOCATE( wfg( nw, 3 ), weight( nw ) )
        ALLOCATE( tagz( nw ) )
        !
        wfg(:,:) = 0
        !
        weight(1) = 0.5D0
        weight(2) = 0.5D0
        !
     ELSE
        !
        IF ( b1(2) > 1 ) THEN
           !
           CALL errore( 'cp-wf', 'Please make celldm(2) not less than 1', 1 )
           !
        ELSE
           !
           nw = 4
           !
           ALLOCATE( wfg( nw, 3 ), weight( nw ) )
           ALLOCATE( tagz( nw ) )
           !
           weight(1) = 0.5D0
           weight(2) = 0.5D0
           weight(4) = 0.25D0 * ( 1.D0 / b1(2)**2 - 1.D0 )
           !
           wfg(:,:) = 0
           wfg(4,1) = 1
           wfg(4,2) = 1
        END IF
        !
     END IF
     !
     weight(3) = 1.D0 / b3(3)**2
     !
     tagz(:) = 1
     tagz(3) = 0
     !
     wfg(1,1) = 1
     wfg(2,2) = 1
     wfg(3,3) = 1
     !
  CASE( 10 )
     !
     ! ... all face centered orthorhombic F
     !
     IF ( b1(2) == -1 .AND. b1(3) == 1 ) &
        CALL errore( 'cp-wf', 'Please change ibrav to 2', 1 )
     !
     IF ( b1(1) > b1(3) .OR. b1(1) + b1(2) < 0 ) &
        CALL errore( 'cp-wf', 'Please make celldm(2) >= 1 >= celldm(3)', 1 )
     !
     IF ( b1(3) == 1 ) THEN
        !
        nw = 5
        !
     ELSE
        !
        nw = 6
        !
     END IF
     !
     ALLOCATE( wfg( nw, 3 ), weight( nw ) )
     ALLOCATE( tagz( nw ) )
     !
     weight(:) = 0.25D0 / b1(3)**2
     !
     tagz(:) = 1
     tagz(3) = 0
     !
     wfg(:,:) = 0
     wfg(1,1) = 1
     wfg(2,2) = 1
     wfg(3,3) = 1
     wfg(4,:) = 1
     wfg(5,2) = 1
     wfg(5,3) = 1
     !
     weight(5) = 0.25D0 / b1(2)**2 - weight(1)
     !
     IF ( b1(3) /= 1 ) THEN
        !
        weight(6) = 0.25D0 - weight(1)
        !
        wfg(6,1) = 1
        wfg(6,2) = 1
        !
     END IF
     !
  CASE( 11 )
     !
     ! ... body centered orthohombic I
     !
     IF ( b1(3) == 1 .AND. b2(2) == 1 ) &
        CALL errore( 'cp-wf', 'Please change ibrav to 3', 1 )
     !
     IF ( b1(3) == 1 .OR. b2(2) == 1 ) &
        CALL errore( 'cp-wf', 'Please change ibrav to 7', 1 )
     !
     nw = 6
     !
     ALLOCATE( wfg( nw, 3 ), weight( nw ) )
     ALLOCATE( tagz( nw ) )
     !
     tagz(:) = 1
     tagz(3) = 0
     !
     t1 = 0.25D0
     t2 = t1 / b2(2)**2
     t3 = t1 / b1(3)**2
     !
     weight(1) =  t1 - t2 + t3
     weight(2) =  t1 + t2 - t3
     weight(3) = -t1 + t2 + t3
     weight(4) = weight(3)
     weight(5) = weight(1)
     weight(6) = weight(2)
     !
     wfg(:,:) =  0
     wfg(1,1) =  1
     wfg(2,2) =  1
     wfg(3,3) =  1
     wfg(4,1) =  1
     wfg(4,2) =  1
     wfg(5,2) =  1
     wfg(5,3) =  1
     wfg(6,1) = -1
     wfg(6,3) =  1
     !
  CASE( 12 )
     !
     ! ... monoclinic P
     !
     nw=4
     !
     ALLOCATE( wfg( nw, 3 ), weight( nw ) )
     ALLOCATE( tagz( nw ) )
     !
     tagz(:) = 1
     tagz(3) = 0
     !
     t1 = - b1(2) / b2(2)
     !
     k = NINT( t1 )
     !
     IF ( k == 0 .AND. t1 >= 0 ) k =  1
     IF ( k == 0 .AND. t1 <= 0 ) k = -1
     !
     t2 = t1 / SQRT( 1.D0 - 1.D0 / ( 1.D0 + b1(2)**2 ) )
     !
     weight(4) = t1 / DBLE( k )
     weight(1) = 1.D0 - weight(4)
     weight(2) = t2**2 - t1 * k
     weight(3) = 1.D0 / b3(3)**2
     !
     wfg(:,:) = 0
     wfg(1,1) = 1
     wfg(2,2) = 1
     wfg(3,3) = 1
     wfg(4,1) = 1
     wfg(4,2) = k
     !
  CASE( 13 )
     !
     ! ... one face centered monoclinic C
     !
     nw = 6
     !
     ALLOCATE( wfg( nw, 3 ), weight( nw ) )
     ALLOCATE( tagz( nw ) )
     !
     tagz(:) = 1
     tagz(3) = 0
     !
     t1 = - b1(2) / b3(2)
     !
     k = NINT( t1 )
     !
     IF (k == 0 .AND. t1 >= 0 ) k =  1
     IF (k == 0 .AND. t1 <= 0 ) k = -1
     !
     t2 = 1.D0 / ( 4.D0 * b1(1)**2 )
     t3 = 1.D0 / b3(2)**2
     !
     weight(1) = 2.D0 * t2 * ( 1.D0 - t1 )
     weight(2) = weight(1)
     weight(3) = t3 + 4.D0 * t2 * t1 * ( t1 - k )
     weight(4) = 1.D0 - 2.D0 * t2
     weight(5) = 2.D0 * t1 * t2 / DBLE( k )
     weight(6) = weight(5)
     !
     wfg(:,:) =  0
     wfg(1,1) =  1
     wfg(2,2) =  1
     wfg(3,3) =  1
     wfg(4,1) =  1
     wfg(4,2) = -1
     wfg(5,1) =  1
     wfg(5,3) =  k
     wfg(6,2) =  1
     wfg(6,3) =  k
     !
  CASE default
     !
     ! ... free cell used for cpr
     !
     nw = 6
     !
     ALLOCATE( wfg( nw, 3 ), weight( nw ) )
     ALLOCATE( tagz( nw ) )
     !
     tagz(:) = 1
     tagz(3) = 0
     !
     ! ... assigns weights for the triclinic/cpr case
     !
     CALL tric_wts( b1, b2, b3, alat, weight )
     !
     wfg(:,:) = 0
     wfg(1,1) = 1
     wfg(2,2) = 1
     wfg(3,3) = 1
     wfg(4,1) = 1
     wfg(4,2) = 1
     wfg(5,2) = 1
     wfg(5,3) = 1
     wfg(6,1) = 1
     wfg(6,3) = 1
     !
  END SELECT
  !
  ! ... set up matrix O
  !
  ALLOCATE( O( nw, nbsp, nbsp ), X( nbsp, nbsp ), Oa( nw, nbsp, nbsp ) )
  !
  IF ( nspin == 2 .AND. nvb > 0 ) THEN
     !
     ALLOCATE( X2( nupdwn(1), nupdwn(1) ) )
     ALLOCATE( X3( nupdwn(2), nupdwn(2) ) )
     !
  END IF
  !
#if defined (__PARA)
  !
  ALLOCATE( ns( nproc ) )
  !
  ns(1:nproc) = 0
  !
  IF ( nbsp == nproc ) THEN
     !
     ns(1:nbsp)=1
     !
  ELSE
     i=0
1    DO j=1,nproc   
        ns(j)=ns(j)+1
        i=i+1
        IF(i.GE.nbsp) GO TO 2
     END DO
     IF(i.LT.nbsp) GO TO 1 
  END IF
2 IF(iprsta.GT.4) THEN
     DO j=1,nproc
        WRITE( stdout, * ) ns(j)
     END DO
  END IF

  total = 0   
  DO proc=1,nproc
     ngpwpp(proc)=(dfftp%nwl(proc)+1)/2
     total=total+ngpwpp(proc)
     !           nstat=ns(proc)
     IF(iprsta.GT.4) THEN
        WRITE( stdout, * ) "I am proceessor", proc, "and i have ",ns(me)," states."
     END IF
  END DO
  nstat=ns(me)

  ALLOCATE(psitot(total,nstat))
  ALLOCATE(psitot1(total*nstat))
  ALLOCATE(psitot_pl(total,nstat))
  ALLOCATE(psitot_p(total*nstat))
  ALLOCATE(psitot_mi(total,nstat))
  ALLOCATE(psitot_m(total*nstat))

  ALLOCATE(c_p(ngw,nbspx))
  ALLOCATE(c_m(ngw,nbspx))
  IF(iprsta.GT.4) THEN
     WRITE( stdout, * ) "All allocations done"
  END IF


  DO proc=1,nproc
     sendcount(proc)=ngpwpp(me)*ns(proc)
     recvcount(proc)=ngpwpp(proc)*ns(me)
  END DO
  sdispls(1)=0
  rdispls(1)=0

  DO proc=2,nproc
     sdispls(proc)=sdispls(proc-1)+sendcount(proc-1)
     rdispls(proc)=rdispls(proc-1)+recvcount(proc-1)
  END DO
  !
  ! ... Step 1. Communicate to all Procs so that each proc has all
  ! ... G-vectors and some states instead of all states and some
  ! ... G-vectors. This information is stored in the 1-d array 
  ! ... psitot1.
  !
  CALL mp_barrier()
  !
  CALL MPI_ALLTOALLV(c, sendcount, sdispls, MPI_DOUBLE_COMPLEX,             &
       &             psitot1, recvcount, rdispls, MPI_DOUBLE_COMPLEX,       &
       &            MPI_COMM_WORLD,ierr)
  IF (ierr.NE.0) CALL errore('WF','alltoallv 1',ierr)
  !
#endif   
  IF(clwf.EQ.5) THEN
     !
#if defined (__PARA)
     !
     CALL write_psi(c,jw)
     CALL MPI_finalize(ierr)
     WRITE( stdout, * ) "State written", jw
     STOP
  END IF
  !
#else
    !
    DO i=1,ngw
       WRITE(22,*) c(i,jw)
    END DO
    WRITE( stdout, * ) "State written", jw
    STOP
  END IF
  !
#endif
#if defined (__PARA)
  !
  !   Step 2. Convert the 1-d array psitot1 into a 2-d array consistent with the
  !   original notation c(ngw,nbsp). Psitot contains ntot = SUM_Procs(ngw) G-vecs
  !   and nstat states instead of all nbsp states
  !
  ngpww=0
  DO proc=1,nproc
     DO i=1,ns(me)
        DO j=1,ngpwpp(proc)
           psitot(j+ngpww,i)=psitot1(rdispls(proc)+j+(i-1)*ngpwpp(proc))
        END DO
     END DO
     ngpww=ngpww+ngpwpp(proc)
  END DO
  IF(iprsta.GT.4) THEN
     WRITE( stdout, * ) "Step 2. Convert the 1-d array psitot1 into a 2-d array... Done, wf"
  END IF
  !
  !   Step 3. do the translation of the 2-d array to get the transtalted
  !   arrays psitot_pl and psittot_mi, corresponding to G+G' and -G+G'
  !   
  DO inw=1,nw   
   !
   !   Intermediate Check. If the translation is only along the z-direction
   !   no interprocessor communication and data rearrangement is required 
   !   because each processor contains all the G- components in the z-dir.
   !
   IF(tagz(inw).EQ.0) THEN
     DO i=1,nbsp
        DO ig=1,ngw
           IF(indexplusz(ig).EQ.-1) THEN
              c_p(ig,i)=(0.D0,0.D0)
           ELSE
              c_p(ig,i)=c(indexplusz(ig),i)
           END IF
           IF(indexminusz(ig).EQ.-1) THEN
              c_m(ig,i)=(0.D0,0.D0)
           ELSE
              c_m(ig,i)=CONJG(c(indexminusz(ig),i))
           END IF
        END DO
     END DO
   ELSE
     DO i=1,ns(me)
        DO ig=1,total
           IF(indexplus(ig,inw).EQ.-1) THEN
              psitot_pl(ig,i)=(0.D0,0.D0)
           ELSE   
              IF(tagp(ig,inw).EQ.1) THEN
                 psitot_pl(ig,i)=CONJG(psitot(indexplus(ig,inw),i))
              ELSE
                 psitot_pl(ig,i)=psitot(indexplus(ig,inw),i)
              END IF
           END IF
           IF(indexminus(ig,inw).EQ.-1) THEN
              psitot_mi(ig,i)=(0.D0,0.D0)
           ELSE
              IF(tag(ig,inw).EQ.1) THEN
                 psitot_mi(ig,i)=CONJG(psitot(indexminus(ig,inw),i))
              ELSE
                 psitot_mi(ig,i)=psitot(indexminus(ig,inw),i)
              END IF
           END IF
        END DO
     END DO
     IF(iprsta.GT.4) THEN
        WRITE( stdout, * ) "Step 3. do the translation of the 2-d array...Done, wf"
     END IF
     !
     !   Step 4. Convert the 2-d arrays psitot_p and psitot_m into 1-d
     !   arrays
     !
     ngpww=0
     DO proc=1,nproc
        DO i=1,ns(me)
           DO j=1,ngpwpp(proc)
              psitot_p(rdispls(proc)+j+(i-1)*ngpwpp(proc))=psitot_pl(j+ngpww,i)
              psitot_m(rdispls(proc)+j+(i-1)*ngpwpp(proc))=psitot_mi(j+ngpww,i)
           END DO
        END DO
        ngpww=ngpww+ngpwpp(proc)
     END DO
     IF(iprsta.GT.4) THEN
        WRITE( stdout, * ) "Convert the 2-d arrays psitot_p and psitot_m into 1-d arrays...Done, wf"
     END IF
     !
     !   Step 5. Redistribute among processors. The result is stored in 2-d
     !   arrays c_p and c_m consistent with the notation c(ngw,nbsp), such that
     !   c_p(j,i) contains the coeffiCIent for c(j,i) corresponding to G+G'
     !       and c_m(j,i) contains the coeffiCIent for c(j,i) corresponding to -G+G'
     !
     c_p = 0.D0
     !
     CALL mp_barrier()
     !
     CALL MPI_alltoallv(psitot_p, recvcount, rdispls, MPI_DOUBLE_COMPLEX,          &
          &                       c_p, sendcount , sdispls, MPI_DOUBLE_COMPLEX,              &
          &                       MPI_COMM_WORLD,ierr)
     IF (ierr.NE.0) CALL errore('WF','alltoallv 2',ierr)

     c_m = 0.D0
     CALL mp_barrier()
     !
     CALL MPI_alltoallv(psitot_m, recvcount, rdispls, MPI_DOUBLE_COMPLEX,          &
          &                       c_m, sendcount, sdispls, MPI_DOUBLE_COMPLEX,               &
          &                       MPI_COMM_WORLD,ierr)
     IF (ierr.NE.0) CALL errore('WF','alltoallv 3',ierr)
     IF(iprsta.GT.4) THEN
        WRITE( stdout, * ) "Step 5. Redistribute among processors...Done, wf"
     END IF
  END IF
    !
#else
    !
  ALLOCATE(c_p(ngw,nbspx))
  ALLOCATE(c_m(ngw,nbspx))
  DO inw=1,nw
     IF(tagz(inw).EQ.0) THEN
        DO i=1,nbsp
           DO ig=1,ngw
              IF(indexplusz(ig).EQ.-1) THEN
                 c_p(ig,i)=(0.D0,0.D0)
              ELSE
                 c_p(ig,i)=c(indexplusz(ig),i)
              END IF
              IF(indexminusz(ig).EQ.-1) THEN
                 c_m(ig,i)=(0.D0,0.D0)
              ELSE
                 c_m(ig,i)=CONJG(c(indexminusz(ig),i))
              END IF
           END DO
        END DO
     ELSE
        DO i=1,nbsp
           DO ig=1,ngw
              IF(indexplus(ig,inw).EQ.-1) THEN
                 c_p(ig,i)=(0.D0,0.D0)
              ELSE
                 IF(tagp(ig,inw).EQ.1) THEN
                    c_p(ig,i)=CONJG(c(indexplus(ig,inw),i))
                 ELSE
                    c_p(ig,i)=c(indexplus(ig,inw),i)
                 END IF
              END IF
              IF(indexminus(ig,inw).EQ.-1) THEN
                 c_m(ig,i)=(0.D0,0.D0)
              ELSE
                 IF(tag(ig,inw).EQ.1) THEN
                    c_m(ig,i)=CONJG(c(indexminus(ig,inw),i))
                 ELSE
                    c_m(ig,i)=c(indexminus(ig,inw),i)
                 END IF
              END IF
           END DO
        END DO
     END IF
    !
#endif
     !
     ! ... Step 6. Calculate Overlaps
     !
     ! ... Augmentation Part first
     !
     ALLOCATE( qv( nnrbx ) )
     !
     X = ZERO
     !
     isa = 1
     DO is = 1, nvb
        DO ia =1, na(is)
           DO iv = 1, nh(is)
              inl = ish(is) + (iv-1)*na(is) + ia
              jv = iv 
              ijv=(jv-1)*jv/2 + iv
              qv( 1 : nnrbx ) = 0.D0 
              DO ig=1,ngb
                 qv(npb(ig))=eigrb(ig,isa)*qgb(ig,ijv,is)
                 qv(nmb(ig))=CONJG(eigrb(ig,isa)*qgb(ig,ijv,is))
              END DO
#ifdef __PARA
              irb3=irb(3,isa)
#endif
              CALL invfft('Box',qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,isa)
              iqv=1
              qvt=(0.D0,0.D0)
              qvt=boxdotgridcplx(irb(1,isa),qv,expo(1,inw))

#ifdef __PARA
              CALL reduce(2, qvt)
#endif
              !
              IF (nspin.EQ.1) THEN
                 bec2(1:nbsp)=(0.D0,0.D0)
                 bec2(1:nbsp)=bec(inl,1:nbsp)*ONE
                 CALL ZSYRK('U','T',nbsp,1,qvt,bec2,1,ONE,X,nbsp)
              ELSE
                 X2=(0.D0,0.D0)
                 X3=(0.D0,0.D0)
                 bec2up(1:nupdwn(1))=(0.D0,0.D0)
                 bec2up(1:nupdwn(1))=bec(inl,1:nupdwn(1))
                 CALL ZSYRK('U','T',nupdwn(1),1,qvt,bec2up,1,ONE,X2,nupdwn(1))
                 bec2dw(1:nupdwn(2))=(0.D0,0.D0)
                 bec2dw(1:nupdwn(2))=bec(inl,iupdwn(2):nbsp)
                 CALL ZSYRK('U','T',nupdwn(2),1,qvt,bec2dw,1,ONE,X3,nupdwn(2))
                 DO i = 1, nupdwn(1)
                    DO j=i, nupdwn(1)
                       X(i,j)=X(i,j)+X2(i,j)
                    END DO
                 END DO
                 DO i = 1,nupdwn(2)
                    DO j=i,nupdwn(2)
                       X(i+nupdwn(1),j+nupdwn(1)) =X(i+nupdwn(1),j+nupdwn(1)) + X3(i,j)
                    END DO
                 END DO
              END IF
              DO jv = iv+1, nh(is)
                 jnl = ish(is) + (jv-1)*na(is) + ia
                 ijv = (jv-1)*jv/2 + iv
                 qv( 1:nnrbx ) = 0.D0
                 DO ig=1,ngb
                    qv(npb(ig))=eigrb(ig,isa)*qgb(ig,ijv,is)
                    qv(nmb(ig))=CONJG(eigrb(ig,isa)*qgb(ig,ijv,is))
                 END DO
                 CALL invfft('Box',qv,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,isa)
                 iqv=1
                 qvt=0.D0
                 qvt=boxdotgridcplx(irb(1,isa),qv,expo(1,inw))
#ifdef __PARA
                 CALL reduce(2, qvt)
#endif
                 !
                 IF (nspin.EQ.1) THEN
                    bec2(1:nbsp)=(0.D0,0.D0)
                    bec3(1:nbsp)=(0.D0,0.D0)
                    bec2(1:nbsp)=bec(inl,1:nbsp)*ONE
                    bec3(1:nbsp)=bec(jnl,1:nbsp)*ONE
                    CALL ZSYR2K('U','T',nbsp,1,qvt,bec2,1,bec3,1,ONE,X,nbsp)
                 ELSE
                    X2=(0.D0,0.D0)
                    X3=(0.D0,0.D0)
                    bec2up(1:nupdwn(1))=(0.D0,0.D0)
                    bec3up(1:nupdwn(1))=(0.D0,0.D0)
                    bec2up(1:nupdwn(1))=bec(inl,1:nupdwn(1))*ONE
                    bec3up(1:nupdwn(1))=bec(jnl,1:nupdwn(1))*ONE
                    CALL ZSYR2K('U','T',nupdwn(1),1,qvt,bec2up,1,bec3up,1,ONE,X2,nupdwn(1))
                    bec2dw(1:nupdwn(2))=(0.D0,0.D0)
                    bec3dw(1:nupdwn(2))=(0.D0,0.D0)
                    bec2dw(1:nupdwn(2))=bec(inl,iupdwn(2):nbsp)*ONE
                    bec3dw(1:nupdwn(2))=bec(jnl,iupdwn(2):nbsp)*ONE
                    CALL ZSYR2K('U','T',nupdwn(2),1,qvt,bec2dw,1,bec3dw,1,ONE,X3,nupdwn(2))
                    DO i = 1, nupdwn(1)
                       DO j=i, nupdwn(1)
                          X(i,j)=X(i,j)+X2(i,j)
                       END DO
                    END DO
                    DO i = 1,nupdwn(2)
                       DO j=i,nupdwn(2)
                          X(i+nupdwn(1),j+nupdwn(1)) =X(i+nupdwn(1),j+nupdwn(1)) + X3(i,j)
                       END DO
                    END DO
                 END IF
              END DO
           END DO
           isa = isa + 1
        END DO
     END DO
     t1=omega/DBLE(nr1*nr2*nr3)
     X=X*t1
     DO i=1, nbsp
        DO j=i+1, nbsp
           X(j, i)=X(i, j)
        END DO
     END DO
     Oa(inw, :, :)=X(:, :)
     IF(iprsta.GT.4) THEN
        WRITE( stdout, * ) "Augmentation Part Done"
     END IF

     DEALLOCATE( qv )

     !   Then Soft Part
     IF(nspin.EQ.1) THEN
        !   Spin Unpolarized calculation
        X=0.D0   
        IF(gstart.EQ.2) THEN
           c_m(1,:)=0.D0
        END IF
        !           cwf(:,:)=ZERO
        !           cwf(:,:)=c(:,:)
        CALL ZGEMM('c','nbsp',nbsp,nbsp,ngw,ONE,c,ngw,c_p,ngw,ONE,X,nbsp)
        CALL ZGEMM('T','nbsp',nbsp,nbsp,ngw,ONE,c,ngw,c_m,ngw,ONE,X,nbsp)
#ifdef __PARA
        CALL reduce (2*nbsp*nbsp,X)
#endif
        O(inw,:,:)=Oa(inw,:,:)+X(:,:)
        IF(iprsta.GT.4) THEN
           WRITE( stdout, * ) "Soft Part Done"
        END IF

     ELSE
        !   Spin Polarized case
        !   Up Spin First
        ALLOCATE(Xsp(nbsp,nupdwn(1)))
        ALLOCATE(c_psp(ngw,nupdwn(1)))
        ALLOCATE(c_msp(ngw,nupdwn(1)))
        Xsp=0.D0
        c_psp=0.D0 
        c_msp=0.D0
        DO i=1,nupdwn(1)
           c_psp(:,i)=c_p(:,i)
           c_msp(:,i)=c_m(:,i)
        END DO
        IF(gstart.EQ.2) THEN
           c_msp(1,:)=0.D0
        END IF
        !           cwf(:,:)=ZERO
        !           cwf(:,:)=c(:,:,1,1)
        CALL ZGEMM('c','nbsp',nbsp,nupdwn(1),ngw,ONE,c,ngw,c_psp,ngw,ONE,Xsp,nbsp)
        CALL ZGEMM('T','nbsp',nbsp,nupdwn(1),ngw,ONE,c,ngw,c_msp,ngw,ONE,Xsp,nbsp)
#ifdef __PARA
        CALL reduce (2*nbsp*nupdwn(1),Xsp)
#endif
        DO i=1,nupdwn(1)
           DO j=1,nbsp
              X(j,i)=Xsp(j,i)
           END DO
        END DO
        DEALLOCATE(Xsp,c_psp,c_msp)
        !    Then Down Spin
        ALLOCATE(Xsp(nbsp,iupdwn(2):nbsp))
        ALLOCATE(c_psp(ngw,iupdwn(2):nbsp))
        ALLOCATE(c_msp(ngw,iupdwn(2):nbsp))
        Xsp=0.D0
        c_psp=0.D0
        c_msp=0.D0
        DO i=iupdwn(2),nbsp
           c_psp(:,i)=c_p(:,i)
           c_msp(:,i)=c_m(:,i)
        END DO
        IF(gstart.EQ.2) THEN
           c_msp(1,:)=0.D0
        END IF
        !           cwf(:,:)=ZERO
        !           cwf(:,:)=c(:,:,1,1)
        CALL ZGEMM('c','nbsp',nbsp,nupdwn(2),ngw,ONE,c,ngw,c_psp,ngw,ONE,Xsp,nbsp)
        CALL ZGEMM('T','nbsp',nbsp,nupdwn(2),ngw,ONE,c,ngw,c_msp,ngw,ONE,Xsp,nbsp)
#ifdef __PARA
        CALL reduce (2*nbsp*nupdwn(2),Xsp)
#endif
        DO i=iupdwn(2),nbsp
           DO j=1,nbsp
              X(j,i)=Xsp(j,i)
           END DO
        END DO
        DEALLOCATE(Xsp,c_psp,c_msp)
        O(inw,:,:)=Oa(inw,:,:)+X(:,:)
     END IF
  END DO
#ifdef __PARA
  DEALLOCATE(ns)
#endif

  IF(clwf.EQ.2) THEN
     !    output the overlap matrix to fort.38
#ifdef __PARA
     IF(me.EQ.1) THEN
#endif
        REWIND 38
        WRITE(38, '(i5, 2i2, i3, f9.5)') nbsp, nw, nspin, ibrav, alat
        IF (nspin.EQ.2) THEN
           WRITE(38, '(i5)') nupdwn(1)
        END IF
        WRITE(38, *) a1
        WRITE(38, *) a2
        WRITE(38, *) a3
        WRITE(38, *) b1
        WRITE(38, *) b2
        WRITE(38, *) b3
        DO inw=1, nw
           WRITE(38, *) wfg(inw, :), weight(inw)
        END DO
        DO inw=1, nw
           DO i=1, nbsp
              DO j=1, nbsp
                 WRITE(38, *) O(inw, i, j)
              END DO
           END DO
        END DO
        DO i=1, nbsp
           DO j=1, nbsp
              WRITE(38, *) Uall(i, j)
           END DO
        END DO
        CLOSE(38)
#ifdef __PARA
     END IF
#endif
#ifdef __PARA
     CALL Mpi_finalize(ierr)
#endif
     STOP
  END IF

  IF(clwf.EQ.3.OR.clwf.EQ.4) THEN
     IF(nspin.EQ.1) THEN
        IF(.NOT.what1) THEN
           IF(wfsd) THEN
              CALL wfsteep(nbsp,O,Uall,b1,b2,b3)
           ELSE
              CALL ddyn(nbsp,O,Uall,b1,b2,b3)
           END IF
        END IF
        IF(iprsta.GT.4) THEN
           WRITE( stdout, * ) "Out from DDYN"
        END IF
     ELSE
        ALLOCATE(Uspin(nupdwn(1), nupdwn(1)), Ospin(nw, nupdwn(1), nupdwn(1)))
        DO i=1, nupdwn(1)
           DO j=1, nupdwn(1)
              Uspin(i, j)=Uall(i, j)
              Ospin(:, i, j)=O(:, i, j)
           END DO
        END DO
        IF(.NOT.what1) THEN
           IF(wfsd) THEN
              CALL wfsteep(nupdwn(1), Ospin, Uspin,b1,b2,b3)
           ELSE 
              CALL ddyn(nupdwn(1), Ospin, Uspin,b1,b2,b3)
           END IF
        END IF
        DO i=1, nupdwn(1)
           DO j=1, nupdwn(1)
              Uall(i, j)=Uspin(i, j)
              O(:,i,j)  =Ospin(:,i,j)
           END DO
        END DO
        DEALLOCATE(Uspin, Ospin)
        ALLOCATE(Uspin(nupdwn(2), nupdwn(2)), Ospin(nw, nupdwn(2), nupdwn(2)))
        DO i=1, nupdwn(2)
           DO j=1, nupdwn(2)
              Uspin(i, j)=Uall(i+nupdwn(1), j+nupdwn(1))
              Ospin(:, i, j)=O(:, i+nupdwn(1), j+nupdwn(1))
           END DO
        END DO
        IF(.NOT.what1) THEN
           IF(wfsd) THEN
              CALL wfsteep(nupdwn(2), Ospin, Uspin,b1,b2,b3)
           ELSE
              CALL ddyn(nupdwn(2), Ospin, Uspin,b1,b2,b3)
           END IF
        END IF
        DO i=1, nupdwn(2)
           DO j=1, nupdwn(2)
              Uall(i+nupdwn(1), j+nupdwn(1))=Uspin(i, j)
              O(:,i+nupdwn(1),j+nupdwn(1))=Ospin(:,i,j)
           END DO
        END DO
        DEALLOCATE(Uspin, Ospin)
     END IF
  END IF

  !       Update C and bec
  cwf=ZERO
  !        cwf(:,:)=c(:,:,1,1)
  becwf=0.0
  U2=Uall*ONE
  CALL ZGEMM('nbsp','nbsp',ngw,nbsp,nbsp,ONE,c,ngw,U2,nbsp,ZERO,cwf,ngw)
  !           call ZGEMM('nbsp','nbsp',ngw,nbsp,nbsp,ONE,cwf,ngw,U2,nbsp,ZERO,cwf,ngw)
  CALL DGEMM('nbsp','nbsp',nkb,nbsp,nbsp,ONE,bec,nkb,Uall,nbsp,ZERO,becwf,nkb)
  U2=ZERO
  IF(iprsta.GT.4) THEN
     WRITE( stdout, * ) "Updating Wafefunctions and Bec"
  END IF

  c(:,:)=cwf(:,:)
  bec(:,:)=becwf(:,:)

  IF(iprsta.GT.4) THEN
     WRITE( stdout, * ) "Wafefunctions and Bec Updated"
  END IF

  !
  ! calculate wannier-function centers
  !
  ALLOCATE( wr( nw ), W( nw, nw ), gr( nw, 3 ), f3( nw ), mt( nw ) )
  DO inw=1, nw
     gr(inw, :)=wfg(inw,1)*b1(:)+wfg(inw,2)*b2(:)+wfg(inw,3)*b3(:)
  END DO
  !
  ! set up a matrix with the element (i,j) is G_i·G_j·weight(j)
  ! to check the correctness of choices on G vectors
  !
  DO i=1, nw
     DO j=1, nw
        W(i,j)=SUM(gr(i,:)*gr(j,:))*weight(j)
        !        write(6, *) i,j,W(i,j)
     END DO
  END DO
  !  write(24, *) "wannier function centers: (unit:\AA)"
  DO i=1, nbsp
     mt=-AIMAG(LOG(O(:,i,i)))/tpi
     wfc(1, i)=SUM(mt*weight*gr(:,1))
     wfc(2, i)=SUM(mt*weight*gr(:,2))
     wfc(3, i)=SUM(mt*weight*gr(:,3))
     DO inw=1, nw
        wr(inw)=SUM(wfc(:,i)*gr(inw,:))-mt(inw)
     END DO
     mt=wr
     f3=0
     adjust=0
     !
     ! ... balance the phase factor if necessary
     !
     wfc(1,i) = ( wfc(1,i) + SUM( mt * weight * gr(:,1) ) ) * alat
     wfc(2,i) = ( wfc(2,i) + SUM( mt * weight * gr(:,2) ) ) * alat
     wfc(3,i) = ( wfc(3,i) + SUM( mt * weight * gr(:,3) ) ) * alat
     !
  END DO
  !
  IF ( ionode ) THEN
     !
     IF ( .NOT. what1 ) THEN
        !
        ! ... pbc are imposed here in the range [0,1]
        !
        DO i = 1, nbsp
           !
           temp_vec(:) = MATMUL( ainv(:,:), wfc(:,i) )
           !
           WHERE( temp_vec(:) >= 1.D0 ) temp_vec(:) = temp_vec(:) - 1.D0
           WHERE( temp_vec(:) < 0.D0 )  temp_vec(:) = temp_vec(:) + 1.D0
           !
           temp_vec(:) = MATMUL( temp_vec(:), h(:,:) )
           !
           WRITE( 26, '(3f11.6)' ) temp_vec(:)
           !
        END DO
        !
     END IF
     !
  END IF
  !
  IF ( nspin == 2 .AND. nvb > 0 ) DEALLOCATE( X2, X3 )
  !
  DEALLOCATE( wr, W, mt, f3, gr )
  !
#if defined (__PARA)
  !
  DEALLOCATE( psitot )
  DEALLOCATE( psitot1 )
  DEALLOCATE( psitot_pl )
  DEALLOCATE( psitot_p )
  DEALLOCATE( psitot_mi )
  DEALLOCATE( psitot_m )
  DEALLOCATE( c_p )
  DEALLOCATE( c_m )
  !
#else
  !
  DEALLOCATE( c_p, c_m )
  !
#endif
  !
  DEALLOCATE( X, O, Oa, wfg, weight, tagz )
  !
  RETURN
  !
END SUBROUTINE wf
!
!----------------------------------------------------------------------------
SUBROUTINE ddyn( m, Omat, Umat, b1, b2, b3 )
  !----------------------------------------------------------------------------
  ! ... This part of the subroutine wf has been added by Manu. It performes
  ! ... Damped Dynamics on the A matrix to get the Unitary transformation to
  ! ... obtain the wannier function at time(t+delta). It also updates the
  ! ... quantities bec and becdr
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE wannier_base,     ONLY : wf_friction, nsteps, tolw, adapt, wf_q, &
                               weight, nw, wfdt
  USE cell_base,        ONLY : alat
  USE constants,        ONLY : tpi
  USE electrons_base,   ONLY : nbsp
  USE control_flags,    ONLY : iprsta
  USE mp_global,        ONLY : mpime
  USE parallel_include
  !
  IMPLICIT NONE
  !
  INTEGER :: f3(nw), f4(nw), i,j,inw
  INTEGER ,INTENT(in) :: m
  REAL(DP), INTENT(in) :: b1(3),b2(3),b3(3)
  REAL(DP), INTENT(inout) :: Umat(m,m)
  COMPLEX(DP), INTENT(inout) :: Omat(nw,m,m)
  COMPLEX(DP) :: U2(m,m),U3(m,m)
  INTEGER :: adjust,ini, ierr1,nnn
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: wr
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: W
  REAL(DP) :: t0, fric,U(m,m), t2
  REAL(DP) :: A(m,m),oldt0,Wm(m,m),U1(m,m)
  REAL(DP) :: Aminus(m,m), Aplus(m,m),f2(4*m)
  REAL(DP) :: temp(m,m)
  COMPLEX(DP) :: d(m,m)
  COMPLEX(DP) :: f1(2*m-1), wp(m*(m+1)/2),z(m,m)
  COMPLEX(DP), ALLOCATABLE, DIMENSION(:, :) :: X1
  COMPLEX(DP), ALLOCATABLE, DIMENSION(:, :, :) :: Oc
  REAL(DP) , ALLOCATABLE , DIMENSION(:) :: mt
  REAL(DP), PARAMETER :: autoaf=0.529177d0
  REAL(DP) :: spread, sp
  REAL(DP) :: wfc(3,nbsp), gr(nw,3)
  INTEGER  :: me
  !
  me = mpime + 1
  !
  ALLOCATE(mt(nw))
  ALLOCATE(X1(m,m))
  ALLOCATE(Oc(nw,m,m))

  fric=wf_friction
  ALLOCATE (W(m,m),wr(m))

  Umat=0.D0
  DO i=1,m
     Umat(i,i)=1.D0
  END DO

  U2=Umat*ONE

  !
  ! update Oc using the initial guess of Uspin
  !
  DO inw=1, nw
     X1(:, :)=Omat(inw, :, :)
     U3=ZERO
     !    call ZGEMUL(U2, m, 'T', X1, m, 'nbsp', U3, m, m,m,m) 
     CALL ZGEMM ('T', 'nbsp', m,m,m,ONE,U2,m,X1,m,ZERO,U3,m)
     X1=ZERO
     !    call ZGEMUL(U3, m, 'nbsp', U2, m, 'nbsp', X1, m, m,m,m) 
     CALL ZGEMM ('nbsp','nbsp', m,m,m,ONE,U3,m,U2,m,ZERO,X1,m)
     Oc(inw, :, :)=X1(:, :)
  END DO

  U2=ZERO
  U3=ZERO

  oldt0=0.D0
  A=0.D0
  Aminus=A
  temp=Aminus


  !   START ITERATIONS HERE

  DO ini=1, nsteps

     t0=0.D0     !use t0 to store the value of omega
     DO inw=1, nw
        DO i=1, m
           t0=t0+DBLE(CONJG(Oc(inw, i, i))*Oc(inw, i, i))
        END DO
     END DO

     IF(ABS(t0-oldt0).LT.tolw) THEN
#ifdef __PARA
        IF(me.EQ.1) THEN
#endif
           WRITE(27,*) "MLWF Generated at Step",ini
#ifdef __PARA
        END IF
#endif
        IF(iprsta.GT.4) THEN
           WRITE( stdout, * ) "MLWF Generated at Step",ini
        END IF
        GO TO 241
     END IF

     IF(adapt) THEN
        IF(oldt0.LT.t0) THEN
           fric=fric/2.
           A=Aminus
           Aminus=temp
        END IF
     END IF

     !   calculate d(omega)/dA and store result in W
     !   this is the force for the damped dynamics
     !

     W=0.D0
     DO inw=1, nw
        t2=weight(inw)
        DO i=1,m
           DO j=1,m
              W(i,j)=W(i,j)+t2*DBLE(Oc(inw,i,j)*CONJG(Oc(inw,i,i)        &
                   -Oc(inw,j,j))+CONJG(Oc(inw,j,i))*(Oc(inw,i,i)-Oc(inw,j,j)))
           END DO
        END DO
     END DO


     !   the verlet scheme to calculate A(t+wfdt)

     Aplus=0.D0

     DO i=1,m
        DO j=i+1,m
           Aplus(i,j)=Aplus(i,j)+(2*wfdt/(2*wfdt+fric))*(2*A(i,j)               &
                -Aminus(i,j)+(wfdt*wfdt/wf_q)*W(i,j)) + (fric/(2*wfdt+fric))*Aminus(i,j)
        ENDDO
     ENDDO

     Aplus=Aplus-TRANSPOSE(Aplus)
     Aplus=(Aplus-A)

     DO i=1, m
        DO j=i,m 
           wp(i + (j-1)*j/2) = CMPLX(0.d0, Aplus(i,j))
        END DO
     END DO

#if ! defined __AIX
     CALL zhpev('V','U',m,wp,wr,z,m,f1,f2,ierr1)
#else
     CALL zhpev(21, wp, wr, z, m, m, f2, 4*m)
     ierr1 = 0
#endif

     IF (ierr1.NE.0) THEN 
        WRITE( stdout, * ) "failed to diagonalize W!"
        STOP
     END IF

     d=0.D0
     DO i=1, m
        d(i, i)=EXP(CI*wr(i)*wfdt)
     END DO      !d=exp(d)

     !   U=z*exp(d)*z+
     !   
     U3=ZERO
     !    call ZGEMUL(z, m, 'nbsp', d, m, 'nbsp', U3, m, m,m,m)
     CALL ZGEMM ('nbsp', 'nbsp', m,m,m,ONE,z,m,d,m,ZERO,U3,m)  
     U2=ZERO
     !    call ZGEMUL(U3, m, 'nbsp', z, m, 'c', U2, m, m,m,m)
     CALL ZGEMM ('nbsp','c', m,m,m,ONE,U3,m,z,m,ZERO,U2,m)
     U=DBLE(U2)
     U2=ZERO
     U3=ZERO

     temp=Aminus
     Aminus=A
     A=Aplus


     !   update Umat
     !
     !    call DGEMUL(Umat, m, 'nbsp', U, m, 'nbsp', U1, m, m,m,m) 
     U1=ZERO
     CALL DGEMM ('nbsp', 'nbsp', m,m,m,ONE,Umat,m,U,m,ZERO,U1,m)

     Umat=U1 

     !   update Oc
     !
     U2=Umat*ONE
     U3=ZERO
     DO inw=1, nw
        X1(:, :)=Omat(inw, :, :)
        !    call ZGEMUL(U2, m, 'T', X1, m, 'nbsp', U3, m, m,m,m)
        CALL ZGEMM ('T', 'nbsp', m,m,m,ONE,U2,m,X1,m,ZERO,U3,m)
        X1=ZERO
        !    call ZGEMUL(U3, m, 'nbsp', U2, m, 'nbsp', X1, m, m,m,m)
        CALL ZGEMM ('nbsp','nbsp',m,m,m,ONE,U3,m,U2,m,ZERO,X1,m)
        Oc(inw, :, :)=X1(:, :)
     END DO
     U2=ZERO
     U3=ZERO

     IF(ABS(t0-oldt0).GE.tolw.AND.ini.GE.nsteps) THEN
#ifdef __PARA
        IF(me.EQ.1) THEN
#endif 
           WRITE(27,*) "MLWF Not generated after",ini,"Steps." 
#ifdef __PARA
        END IF
#endif
        IF(iprsta.GT.4) THEN
           WRITE( stdout, * ) "MLWF Not generated after",ini,"Steps." 
        END IF
        GO TO 241
     END IF

     oldt0=t0

  END DO

241 DEALLOCATE(wr, W)
  spread=0.0
  DO i=1, m
     mt=1.D0-DBLE(Oc(:,i,i)*CONJG(Oc(:,i,i)))
     sp= (alat*autoaf/tpi)**2*SUM(mt*weight)
#ifdef __PARA
     IF(me.EQ.1) THEN
#endif
        WRITE(25, '(f10.7)') sp
#ifdef __PARA
     END IF
#endif
     IF ( sp < 0.D0 ) &
        CALL errore( 'cp-wf', 'Something wrong WF Spread negative', 1 )
     !
     spread=spread+sp
     !
  END DO
  spread=spread/m

#ifdef __PARA
  IF(me.EQ.1) THEN
#endif
     WRITE(24, '(f10.7)') spread
     WRITE(27,*) "Average spread = ", spread
#ifdef __PARA
  END IF
#endif
  Omat=Oc
  IF(iprsta.GT.4) THEN
     WRITE( stdout, * ) "Average spread = ", spread
  END IF
  !
  DEALLOCATE (mt,X1,Oc)
  !
  IF(iprsta.GT.4) THEN
     WRITE( stdout, * ) "Leaving DDYN"
  END IF
  RETURN
END SUBROUTINE ddyn
!
!----------------------------------------------------------------------------
SUBROUTINE wfunc_init( clwf, b1, b2, b3, ibrav )
  !----------------------------------------------------------------------------
  !
  USE io_global,          ONLY : stdout
  USE kinds,              ONLY : DP
  USE reciprocal_vectors, ONLY : gx, mill_l, gstart
  USE gvecw,              ONLY : ngw
  USE electrons_base,     ONLY : nbsp
  USE wannier_base,       ONLY : gnx, gnn, indexplus, indexminus, &
                                 indexplusz, indexminusz, tag, tagp
  USE cvan,               ONLY : nvb
  USE mp,                 ONLY : mp_barrier, mp_bcast
  USE mp_global,          ONLY : nproc, mpime
  USE fft_base,           ONLY : dfftp
  USE parallel_include     
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(in) :: b1(3),b2(3),b3(3)
#ifdef __PARA
  INTEGER :: ntot, proc, ierr, root, i,j,inw,ngppp(nproc)
  INTEGER :: ii,ig,recvcount(nproc), sendcount(nproc),displs(nproc)
#else
  INTEGER :: ierr, i,j,inw, ntot
  INTEGER :: ii,ig
#endif
  REAL (DP), ALLOCATABLE:: bigg(:,:)
  INTEGER, ALLOCATABLE :: bign(:,:)
  INTEGER :: igcount,nw1,jj,nw2, in, kk, ibrav
  INTEGER, ALLOCATABLE :: i_1(:), j_1(:), k_1(:)
  INTEGER :: ti, tj, tk
  REAL(DP) ::t1, vt, err1, err2, err3
  INTEGER :: ti1,tj1,tk1, clwf
  INTEGER :: me
  !
  me = mpime + 1
  !
  IF ( nbsp < nproc ) &
     CALL errore( 'cp-wf', &
                & 'Number of Processors is greater than the number of states', 1 )
  !
  ALLOCATE(gnx(3,ngw))
  ALLOCATE(gnn(3,ngw))
#ifdef __PARA
  root=0
#endif
  vt=1.0d-4
  j=0
  DO i=1,ngw
     gnx(1,i)=gx(1,i)
     gnx(2,i)=gx(2,i)
     gnx(3,i)=gx(3,i)
     gnn(1,i)=mill_l(1,i)
     gnn(2,i)=mill_l(2,i)
     gnn(3,i)=mill_l(3,i)
  END DO
#ifdef __PARA
  ntot=0
  DO i=1,nproc
     ngppp(i)=(dfftp%nwl(i)+1)/2
  END DO

  DO proc=1,nproc
     recvcount(proc)=ngppp(proc)*3
     IF(proc.EQ.1) THEN
        displs(proc)=0
     ELSE
        displs(proc)=displs(proc-1)+recvcount(proc-1)
     END IF
     ntot=ntot+recvcount(proc)/3
  END DO

  IF(me.EQ.1) THEN
     ALLOCATE(bigg(3,ntot))
     ALLOCATE(bign(3,ntot))
  END IF
#else
  ntot=ngw
  ALLOCATE(bigg(3,ntot))
  ALLOCATE(bign(3,ntot))
  bigg(1:3,1:ntot)=gnx(1:3,1:ntot)
  bign(1:3,1:ntot)=gnn(1:3,1:ntot)
#endif

  SELECT CASE(ibrav)
  CASE(0)
     !       free cell for cpr

     nw1=6
     WRITE( stdout, * ) "Translations to be done", nw1
     ALLOCATE(indexplus(ntot,nw1))
     ALLOCATE(indexplusz(ngw))
     ALLOCATE(indexminus(ntot,nw1))
     ALLOCATE(indexminusz(ngw))
     ALLOCATE(tag(ntot,nw1))
     ALLOCATE(tagp(ntot,nw1))
     ALLOCATE(i_1(nw1))
     ALLOCATE(j_1(nw1))
     ALLOCATE(k_1(nw1))

     i_1(1)=1     ! 1
     j_1(1)=0     ! 0
     k_1(1)=0     ! 0

     i_1(2)=0     ! 0
     j_1(2)=1     ! 1
     k_1(2)=0     ! 0

     i_1(3)=0     ! 0
     j_1(3)=0     ! 0
     k_1(3)=1     ! 1

     i_1(4)=1     ! 1
     j_1(4)=1     ! 1
     k_1(4)=0     ! 0

     i_1(5)=0     ! 0
     j_1(5)=1     ! 1
     k_1(5)=1     ! 1

     i_1(6)=1     ! 1
     j_1(6)=0     ! 0
     k_1(6)=1     ! 1

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  CASE(1)
     !   cubic P [sc]

     nw1=3
     WRITE( stdout, * ) "Translations to be done", nw1
     ALLOCATE(indexplus(ntot,nw1))
     ALLOCATE(indexminus(ntot,nw1))
     ALLOCATE(tag(ntot,nw1))
     ALLOCATE(tagp(ntot,nw1))
     ALLOCATE(indexplusz(ngw))
     ALLOCATE(indexminusz(ngw))
     ALLOCATE(i_1(nw1))
     ALLOCATE(j_1(nw1))
     ALLOCATE(k_1(nw1))

     i_1(1)=1   ! 1
     j_1(1)=0   ! 0 
     k_1(1)=0   ! 0

     i_1(2)=0   ! 0 
     j_1(2)=1   ! 1
     k_1(2)=0   ! 0

     i_1(3)=0   ! 0
     j_1(3)=0   ! 0
     k_1(3)=1   ! 1

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  CASE(2)
     !   cubic F [FCC]
     nw1=4
     WRITE( stdout, * ) "Translations to be done", nw1
     ALLOCATE(i_1(nw1))
     ALLOCATE(j_1(nw1))
     ALLOCATE(k_1(nw1))
     ALLOCATE(indexplus(ntot,nw1))
     ALLOCATE(indexminus(ntot,nw1))
     ALLOCATE(tag(ntot,nw1))
     ALLOCATE(tagp(ntot,nw1))
     ALLOCATE(indexplusz(ngw))
     ALLOCATE(indexminusz(ngw))

     i_1(1)=1   ! 1
     j_1(1)=0   ! 0
     k_1(1)=0   ! 0

     i_1(2)=0   ! 0
     j_1(2)=1   ! 1
     k_1(2)=0   ! 0

     i_1(3)=0   ! 0
     j_1(3)=0   ! 0
     k_1(3)=1   ! 1 

     i_1(4)=-1    ! -1
     j_1(4)=-1    ! -1
     k_1(4)=-1    ! -1

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  CASE(3)
     !       cubic I [bcc]
     nw1=6
     WRITE( stdout, * ) "Translations to be done", nw1
     ALLOCATE(i_1(nw1))
     ALLOCATE(j_1(nw1))
     ALLOCATE(k_1(nw1))
     ALLOCATE(indexplus(ntot,nw1))
     ALLOCATE(indexminus(ntot,nw1))
     ALLOCATE(tag(ntot,nw1))
     ALLOCATE(tagp(ntot,nw1))
     ALLOCATE(indexplusz(ngw))
     ALLOCATE(indexminusz(ngw))

     i_1(1)=1     ! 1
     j_1(1)=0     ! 0
     k_1(1)=0     ! 0

     i_1(2)=0     ! 0
     j_1(2)=1     ! 1
     k_1(2)=0     ! 0

     i_1(3)=0   ! 0
     j_1(3)=0   ! 0
     k_1(3)=1   ! 1

     i_1(4)=1   ! 1
     j_1(4)=1   ! 1
     k_1(4)=0   ! 0   

     i_1(5)=1   ! 1
     j_1(5)=1   ! 1
     k_1(5)=1   ! 1

     i_1(6)=-1   ! -1
     j_1(6)=0   ! 0
     k_1(6)=1   ! 1

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  CASE(4)
     !   hexagonal and trigonal P

     nw1=4
     WRITE( stdout, * ) "Translations to be done", nw1
     ALLOCATE(i_1(nw1))
     ALLOCATE(j_1(nw1))
     ALLOCATE(k_1(nw1))
     ALLOCATE(indexplus(ntot,nw1))
     ALLOCATE(indexminus(ntot,nw1))
     ALLOCATE(tag(ntot,nw1))
     ALLOCATE(tagp(ntot,nw1))
     ALLOCATE(indexplusz(ngw))
     ALLOCATE(indexminusz(ngw))

     i_1(1)=1     ! 1
     j_1(1)=0     ! 0
     k_1(1)=0     ! 0

     i_1(2)=0     ! 0
     j_1(2)=1     ! 1
     k_1(2)=0     ! 0

     i_1(3)=0     ! 0
     j_1(3)=0     ! 0
     k_1(3)=1     ! 1

     i_1(4)=1     ! 1
     j_1(4)=-1    ! -1
     k_1(4)=0     ! 0

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  CASE(5)
     !   trigonal R

     nw1=6   
     WRITE( stdout, * ) "Translations to be done", nw1
     ALLOCATE(i_1(nw1))
     ALLOCATE(j_1(nw1))
     ALLOCATE(k_1(nw1))
     ALLOCATE(indexplus(ntot,nw1))
     ALLOCATE(indexminus(ntot,nw1))
     ALLOCATE(tag(ntot,nw1))
     ALLOCATE(tagp(ntot,nw1))
     ALLOCATE(indexplusz(ngw))
     ALLOCATE(indexminusz(ngw))

     i_1(1)=1     ! 1
     j_1(1)=0     ! 0
     k_1(1)=0     ! 0

     i_1(2)=0     ! 0
     j_1(2)=1     ! 1
     k_1(2)=0     ! 0

     i_1(3)=0   ! 0
     j_1(3)=0   ! 0
     k_1(3)=1   ! 1

     i_1(4)=1   ! 1
     j_1(4)=0   ! 0
     k_1(4)=1   ! 1

     i_1(5)=0   ! 0
     j_1(5)=-1   ! -1
     k_1(5)=1   ! 1

     i_1(6)=1   ! 1
     j_1(6)=-1   ! -1
     k_1(6)=0   ! 0

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  CASE(6)
     ! tetragonal P [st]

     nw1=3
     WRITE( stdout, * ) "Translations to be done", nw1
     ALLOCATE(indexplus(ntot,nw1))
     ALLOCATE(indexminus(ntot,nw1))
     ALLOCATE(tag(ntot,nw1))
     ALLOCATE(tagp(ntot,nw1))
     ALLOCATE(indexplusz(ngw))
     ALLOCATE(indexminusz(ngw))

     ALLOCATE(i_1(nw1))
     ALLOCATE(j_1(nw1))
     ALLOCATE(k_1(nw1))

     i_1(1)=1     ! 1
     j_1(1)=0     ! 0
     k_1(1)=0     ! 0

     i_1(2)=0     ! 0
     j_1(2)=1     ! 1
     k_1(2)=0     ! 0

     i_1(3)=0     ! 0
     j_1(3)=0     ! 0
     k_1(3)=1   ! 1

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  CASE(7)
     ! tetragonal I [bct]

     nw1=6
     WRITE( stdout, * ) "Translations to be done", nw1
     ALLOCATE(indexplus(ntot,nw1))
     ALLOCATE(indexminus(ntot,nw1))
     ALLOCATE(tag(ntot,nw1))
     ALLOCATE(tagp(ntot,nw1))
     ALLOCATE(indexplusz(ngw))
     ALLOCATE(indexminusz(ngw))
     ALLOCATE(i_1(nw1))
     ALLOCATE(j_1(nw1))
     ALLOCATE(k_1(nw1))

     i_1(1)=1     ! 1
     j_1(1)=0     ! 0
     k_1(1)=0     ! 0

     i_1(2)=0     ! 0
     j_1(2)=1     ! 1
     k_1(2)=0     ! 0

     i_1(3)=0     ! 0
     j_1(3)=0     ! 0
     k_1(3)=1     ! 1

     i_1(4)=1   ! 1
     j_1(4)=0   ! 0
     k_1(4)=1   ! 1

     i_1(5)=0   ! 0
     j_1(5)=1   ! 1
     k_1(5)=-1   ! -1

     i_1(6)=1   ! 1
     j_1(6)=1   ! 1
     k_1(6)=0   ! 0

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  CASE(8)
     ! Orthorhombic P
     nw1=3
     WRITE( stdout, * ) "Translations to be done", nw1
     ALLOCATE(indexplus(ntot,nw1))
     ALLOCATE(indexminus(ntot,nw1))
     ALLOCATE(tag(ntot,nw1))
     ALLOCATE(tagp(ntot,nw1))
     ALLOCATE(indexplusz(ngw))
     ALLOCATE(indexminusz(ngw))
     ALLOCATE(i_1(nw1))
     ALLOCATE(j_1(nw1))
     ALLOCATE(k_1(nw1))
     indexplus=0
     indexminus=0
     tag=0
     tagp=0
     indexplusz=0
     indexminusz=0

     i_1(1)=1     ! 1
     j_1(1)=0     ! 0
     k_1(1)=0     ! 0

     i_1(2)=0     ! 0
     j_1(2)=1     ! 1
     k_1(2)=0     ! 0

     i_1(3)=0     ! 0
     j_1(3)=0     ! 0
     k_1(3)=1     ! 1

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  CASE(9)
     ! One face centered Orthorhombic C
     IF(b1(2).EQ.1) THEN 

        nw1=3
        WRITE( stdout, * ) "Translations to be done", nw1
        ALLOCATE(indexplus(ntot,nw1))
        ALLOCATE(indexminus(ntot,nw1))
        ALLOCATE(tag(ntot,nw1))
        ALLOCATE(tagp(ntot,nw1))
        ALLOCATE(indexplusz(ngw))
        ALLOCATE(indexminusz(ngw))
        ALLOCATE(i_1(nw1))
        ALLOCATE(j_1(nw1))
        ALLOCATE(k_1(nw1))

     ELSE
        IF ( b1(2) > 1 ) THEN
           !
           CALL errore( 'cp-wf', 'Please make celldm(2) not less than 1', 1 )
           !
        ELSE
           nw1=4
           WRITE( stdout, * ) "Translations to be done", nw1
           ALLOCATE(indexplus(ntot,nw1))
           ALLOCATE(indexminus(ntot,nw1))
           ALLOCATE(tag(ntot,nw1))
           ALLOCATE(tagp(ntot,nw1))
           ALLOCATE(indexplusz(ngw))
           ALLOCATE(indexminusz(ngw))

           ALLOCATE(i_1(nw1))
           ALLOCATE(j_1(nw1))
           ALLOCATE(k_1(nw1))

           i_1(4)=1
           j_1(4)=1
           k_1(4)=0
        END IF
     END IF

     i_1(1)=1     ! 1
     j_1(1)=0     ! 0
     k_1(1)=0     ! 0

     i_1(2)=0     ! 0
     j_1(2)=1     ! 1
     k_1(2)=0     ! 0

     i_1(3)=0     ! 0
     j_1(3)=0     ! 0
     k_1(3)=1     ! 1

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  CASE(10)
     ! all face centered orthorhombic F

     IF ( b1(2) == -1 .AND. b1(3) == 1 ) &
        CALL errore( 'cp-wf', 'Please change ibrav to 2', 1 )
     !
     IF ( ( b1(1) > b1(3) ) .OR. ( b1(1) + b1(2) < 0 ) ) &
        CALL errore( 'cp-wf', 'Please make celldm(2)>=1>=celldm(3)', 1 )
     !
     IF(b1(3).EQ.1) THEN
        nw1=5
        WRITE( stdout, * ) "Translations to be done", nw1
        ALLOCATE(indexplus(ntot,nw1))
        ALLOCATE(indexminus(ntot,nw1))
        ALLOCATE(tag(ntot,nw1))
        ALLOCATE(tagp(ntot,nw1))
        ALLOCATE(indexplusz(ngw))
        ALLOCATE(indexminusz(ngw))

        ALLOCATE(i_1(nw1))
        ALLOCATE(j_1(nw1))
        ALLOCATE(k_1(nw1))
     ELSE
        nw1=6
        WRITE( stdout, * ) "Translations to be done", nw1
        ALLOCATE(indexplus(ntot,nw1))
        ALLOCATE(indexminus(ntot,nw1))
        ALLOCATE(tag(ntot,nw1))
        ALLOCATE(tagp(ntot,nw1))
        ALLOCATE(indexplusz(ngw))
        ALLOCATE(indexminusz(ngw))

        ALLOCATE(i_1(nw1))
        ALLOCATE(j_1(nw1))
        ALLOCATE(k_1(nw1))

        i_1(6)=1   ! 1
        j_1(6)=1   ! 1
        k_1(6)=0   ! 0

     END IF

     i_1(1)=1     ! 1
     j_1(1)=0     ! 0
     k_1(1)=0     ! 0

     i_1(2)=0     ! 0
     j_1(2)=1     ! 1
     k_1(2)=0     ! 0

     i_1(3)=0     ! 0
     j_1(3)=0     ! 0
     k_1(3)=1     ! 1

     i_1(4)=1   ! 1
     j_1(4)=1    ! 1
     k_1(4)=1   ! 1

     i_1(5)=0   ! 0
     j_1(5)=1   ! 1
     k_1(5)=0   ! 1

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  CASE(11) 
     ! Body centered orthorhombic I

     IF ( ( b1(3) == 1 ) .AND. ( b2(2) == 1 ) ) &
        CALL errore( 'cp-wf', 'Please change ibrav to 3', 1 )
     !
     IF ( ( b1(3) == 1 ) .OR. ( b2(2) == 1 ) ) &
        CALL errore( 'cp-wf', 'Please change ibrav to 7', 1 )
     !
     nw1=6
     WRITE( stdout, * ) "Translations to be done", nw1
     ALLOCATE(indexplus(ntot,nw1))
     ALLOCATE(indexminus(ntot,nw1))
     ALLOCATE(tag(ntot,nw1))
     ALLOCATE(tagp(ntot,nw1))
     ALLOCATE(indexplusz(ngw))
     ALLOCATE(indexminusz(ngw))

     ALLOCATE(i_1(nw1))
     ALLOCATE(j_1(nw1))
     ALLOCATE(k_1(nw1))

     i_1(1)=1     ! 1
     j_1(1)=0     ! 0
     k_1(1)=0     ! 0

     i_1(2)=0     ! 0
     j_1(2)=1     ! 1
     k_1(2)=0     ! 0

     i_1(3)=0     ! 0
     j_1(3)=0     ! 0
     k_1(3)=1     ! 1

     i_1(4)=1     ! 1
     j_1(4)=1     ! 1
     k_1(4)=0     ! 0

     i_1(5)=1     ! 1
     j_1(5)=0     ! 0
     k_1(5)=1     ! 1

     i_1(6)=-1    ! -1
     j_1(6)=0     ! 0
     k_1(6)=1     ! 1

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  CASE(12)
     ! monoclinic P

     nw1=4
     WRITE( stdout, * ) "Translations to be done", nw1
     ALLOCATE(indexplus(ntot,nw1))
     ALLOCATE(indexminus(ntot,nw1))
     ALLOCATE(tag(ntot,nw1))
     ALLOCATE(tagp(ntot,nw1))
     ALLOCATE(indexplusz(ngw))
     ALLOCATE(indexminusz(ngw))
     ALLOCATE(i_1(nw1))
     ALLOCATE(j_1(nw1))
     ALLOCATE(k_1(nw1))

     t1=-b1(2)/b2(2)
     kk=NINT(t1)
     IF((kk.EQ.0).AND.(t1.GE.0)) kk=1
     IF((kk.EQ.0).AND.(t1.LE.0)) kk=-1

     i_1(1)=1     ! 1
     j_1(1)=0     ! 0
     k_1(1)=0     ! 0

     i_1(2)=0     ! 0
     j_1(2)=1     ! 1
     k_1(2)=0     ! 0

     i_1(3)=0     ! 0
     j_1(3)=0     ! 0
     k_1(3)=1     ! 1

     i_1(4)=1     ! 1
     j_1(4)=kk    ! kk
     k_1(4)=0     ! 0

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  CASE(13)
     ! one face centered monoclinic C

     nw1=6
     WRITE( stdout, * ) "Translations to be done", nw1
     ALLOCATE(indexplus(ntot,nw1))
     ALLOCATE(indexminus(ntot,nw1))
     ALLOCATE(tag(ntot,nw1))
     ALLOCATE(tagp(ntot,nw1))
     ALLOCATE(indexplusz(ngw))
     ALLOCATE(indexminusz(ngw))
     ALLOCATE(i_1(nw1))
     ALLOCATE(j_1(nw1))
     ALLOCATE(k_1(nw1))

     t1=-b1(2)/b3(2)
     kk=NINT(t1)
     IF((kk.EQ.0).AND.(t1.GE.0)) kk=1
     IF((kk.EQ.0).AND.(t1.LE.0)) kk=-1

     i_1(1)=1     ! 1
     j_1(1)=0     ! 0
     k_1(1)=0     ! 0

     i_1(2)=0     ! 0
     j_1(2)=1     ! 1
     k_1(2)=0     ! 0

     i_1(3)=0     ! 0
     j_1(3)=0     ! 0
     k_1(3)=1     ! 1

     i_1(4)=1     ! 1
     j_1(4)=-1    ! -1
     k_1(4)=0     ! 0

     i_1(5)=1   ! 1
     j_1(5)=0   ! 0
     k_1(5)=kk   ! kk

     i_1(6)=0   ! 0
     j_1(6)=1   ! 1
     k_1(6)=kk   ! kk

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  CASE default
     !       free cell for cpr

     nw1=6
     WRITE( stdout, * ) "Translations to be done", nw1
     ALLOCATE(indexplus(ntot,nw1))
     ALLOCATE(indexplusz(ngw))
     ALLOCATE(indexminus(ntot,nw1))
     ALLOCATE(indexminusz(ngw))
     ALLOCATE(tag(ntot,nw1))
     ALLOCATE(tagp(ntot,nw1))
     ALLOCATE(i_1(nw1))
     ALLOCATE(j_1(nw1))
     ALLOCATE(k_1(nw1))

     i_1(1)=1     ! 1
     j_1(1)=0     ! 0
     k_1(1)=0     ! 0

     i_1(2)=0     ! 0
     j_1(2)=1     ! 1
     k_1(2)=0     ! 0

     i_1(3)=0     ! 0
     j_1(3)=0     ! 0
     k_1(3)=1     ! 1

     i_1(4)=1     ! 1
     j_1(4)=1     ! 1
     k_1(4)=0     ! 0

     i_1(5)=0     ! 0
     j_1(5)=1     ! 1
     k_1(5)=1     ! 1

     i_1(6)=1     ! 1
     j_1(6)=0     ! 0
     k_1(6)=1     ! 1

     indexplus(:,3)=0
     indexminus(:,3)=0
     tag(:,3)=0
     tagp(:,3)=0

  END SELECT

  WRITE( stdout, * ) "ibrav selected:", ibrav
  !
  IF(nvb.GT.0) CALL small_box_wf(i_1, j_1, k_1, nw1)
#ifdef __PARA
  CALL mp_barrier()
  !
  CALL MPI_GATHERV(gnx, recvcount(me), MPI_REAL8,              &
       bigg, recvcount, displs, MPI_REAL8,         &
       root, MPI_COMM_WORLD, ierr)
  IF (ierr.NE.0) CALL errore('wfunc_init','MPI_GATHERV' , ierr)
  !
  CALL mp_barrier()
  !
  CALL MPI_GATHERV(gnn, recvcount(me), MPI_INTEGER,              &
       bign, recvcount, displs, MPI_INTEGER,         &
       root, MPI_COMM_WORLD, ierr)
  IF (ierr.NE.0) CALL errore('wfunc_init','MPI_GATHERV' , ierr)
#endif
#ifdef __PARA
  IF(me.EQ.1) THEN
#endif
     IF(clwf.EQ.5) THEN
#ifdef __PARA
        DO ii=1,ntot
           WRITE(21,*) bigg(:,ii)
        END DO
#else
        DO ii=1,ngw
           WRITE(21,*) gx(1,ii), gx(2,ii), gx(3,ii)
        END DO
#endif
        CLOSE(21)
     END IF
#ifdef __PARA 
  END IF
#endif

  DO inw=1,nw1
     IF(i_1(inw).EQ.0.AND.j_1(inw).EQ.0) THEN
        DO ig=1,ngw
           IF(gstart.EQ.2) THEN
              indexminusz(1)=-1
           END IF
           !           ti=(gnn(1,ig)+i_1(inw))*b1(1)+(gnn(2,ig)+j_1(inw))*b2(1)+(gnn(3,ig)+k_1(inw))*b3(1)
           !           tj=(gnn(1,ig)+i_1(inw))*b1(2)+(gnn(2,ig)+j_1(inw))*b2(2)+(gnn(3,ig)+k_1(inw))*b3(2)
           !           tk=(gnn(1,ig)+i_1(inw))*b1(3)+(gnn(2,ig)+j_1(inw))*b2(3)+(gnn(3,ig)+k_1(inw))*b3(3)
           ti=(gnn(1,ig)+i_1(inw))
           tj=(gnn(2,ig)+j_1(inw))
           tk=(gnn(3,ig)+k_1(inw))
           DO ii=1,ngw
              err1=ABS(gnx(1,ii)-ti)
              err2=ABS(gnx(2,ii)-tj)
              err3=ABS(gnx(3,ii)-tk)
              IF(gnn(1,ii).EQ.ti.AND.gnn(2,ii).EQ.tj.AND.gnn(3,ii).EQ.tk) THEN
                 !             if(err1.lt.vt.and.err2.lt.vt.and.err3.lt.vt) then
                 indexplusz(ig)=ii
                 !               write (6,*) "Found +", ig,ii,inw, ti,tj,tk
                 !               write (6,*) "looking for", ti,tj,tk
                 GO TO 224
              ELSE
              END IF
           END DO
           indexplusz(ig)=-1
           !                write (6,*) "Not Found +", ig,-1,inw
           !               write (6,*) "looking for", ti,tj,tk
           !224        ti=(-gnn(1,ig)+i_1(inw))*b1(1)+(-gnn(2,ig)+j_1(inw))*b2(1)+(-gnn(3,ig)+k_1(inw))*b3(1)
           !           tj=(-gnn(1,ig)+i_1(inw))*b1(2)+(-gnn(2,ig)+j_1(inw))*b2(2)+(-gnn(3,ig)+k_1(inw))*b3(2)
           !           tk=(-gnn(1,ig)+i_1(inw))*b1(3)+(-gnn(2,ig)+j_1(inw))*b2(3)+(-gnn(3,ig)+k_1(inw))*b3(3)
224        ti=(-gnn(1,ig)+i_1(inw))
           tj=(-gnn(2,ig)+j_1(inw))
           tk=(-gnn(3,ig)+k_1(inw))
           ti1=-gnn(1,ig)+i_1(inw)
           tj1=-gnn(2,ig)+j_1(inw)
           tk1=-gnn(3,ig)+k_1(inw)
           IF(ti1.LT.0.OR.(ti1.EQ.0.AND.(tj1.LT.0.OR.(tj1.EQ.0.AND.tk1.LT.0)))) THEN
              DO ii=1,ngw
                 err1=ABS(gnx(1,ii)+ti)
                 err2=ABS(gnx(2,ii)+tj)
                 err3=ABS(gnx(3,ii)+tk)
                 IF(gnn(1,ii).EQ.-ti.AND.gnn(2,ii).EQ.-tj.AND.gnn(3,ii).EQ.-tk) THEN
                    !                    if(err1.lt.vt.and.err2.lt.vt.and.err3.lt.vt) then
                    indexminusz(ig)=ii
                    !                     tag(ig,inw)=1
                    !                     write (6,*) "Found -", ig,ii,inw
                    !               write (6,*) "looking for", -ti,-tj,-tk
                    GO TO 223
                 ELSE
                 END IF
              END DO
              indexminusz(ig)=-1
              !                tag(ig,inw)=1
              !                write (6,*) "Not Found -", ig,-1,inw
              !               write (6,*) "looking for", -ti,-tj,-tk
           ELSE
              DO ii=1,ngw
                 err1=ABS(gnx(1,ii)-ti)
                 err2=ABS(gnx(2,ii)-tj)
                 err3=ABS(gnx(3,ii)-tk)
                 IF(gnn(1,ii).EQ.ti.AND.gnn(2,ii).EQ.tj.AND.gnn(3,ii).EQ.tk) THEN
                    !                  if(err1.lt.vt.and.err2.lt.vt.and.err3.lt.vt) then
                    indexminusz(ig)=ii
                    !                   tag(ig,inw)=-1
                    !                   write (6,*) "Found -", ig,ii,inw
                    !               write (6,*) "looking for", ti,tj,tk
                    GO TO 223
                 ELSE
                 END IF
              END DO
              indexminusz(ig)=-1
              !              tag(ig,inw)=-1
              !              write (6,*) "Not Found -", ig,-1,inw
              !               write (6,*) "looking for", ti,tj,tk
           END IF
223        CONTINUE
        END DO
        WRITE( stdout, * ) "Translation", inw, "for", ngw, "G vectors"
     ELSE
#ifdef __PARA
        IF(me.EQ.1) THEN   
#endif
           DO ig=1,ntot
              IF(gstart.EQ.2) THEN
                 indexminus(1,inw)=-1
              END IF
              !           ti=(bign(1,ig)+i_1(inw))*b1(1)+(bign(2,ig)+j_1(inw))*b2(1)+(bign(3,ig)+k_1(inw))*b3(1)
              !           tj=(bign(1,ig)+i_1(inw))*b1(2)+(bign(2,ig)+j_1(inw))*b2(2)+(bign(3,ig)+k_1(inw))*b3(2)
              !           tk=(bign(1,ig)+i_1(inw))*b1(3)+(bign(2,ig)+j_1(inw))*b2(3)+(bign(3,ig)+k_1(inw))*b3(3)
              ti=(bign(1,ig)+i_1(inw))
              tj=(bign(2,ig)+j_1(inw))
              tk=(bign(3,ig)+k_1(inw))
              ti1=bign(1,ig)+i_1(inw)
              tj1=bign(2,ig)+j_1(inw)
              tk1=bign(3,ig)+k_1(inw)
              IF(ti1.LT.0.OR.(ti1.EQ.0.AND.(tj1.LT.0.OR.(tj1.EQ.0.AND.tk1.LT.0)))) THEN
                 DO ii=1,ntot
                    err1=ABS(bigg(1,ii)+ti)
                    err2=ABS(bigg(2,ii)+tj)
                    err3=ABS(bigg(3,ii)+tk)
                    !              if(err1.lt.vt.and.err2.lt.vt.and.err3.lt.vt) then
                    IF(bign(1,ii).EQ.-ti.AND.bign(2,ii).EQ.-tj.AND.bign(3,ii).EQ.-tk) THEN
                       indexplus(ig,inw)=ii
                       tagp(ig,inw)=1
                       !                write (6,*) "Found +", ig,ii,inw 
                       !               write (6,*) "looking for", -ti,-tj,-tk
                       GO TO 214
                    ELSE
                    END IF
                 END DO
                 indexplus(ig,inw)=-1
                 tagp(ig,inw)=1
                 !          write (6,*) "Not Found +", ig,-1,inw 
                 !               write (6,*) "looking for", -ti,-tj,-tk
              ELSE
                 DO ii=1,ntot
                    err1=ABS(bigg(1,ii)-ti)
                    err2=ABS(bigg(2,ii)-tj)
                    err3=ABS(bigg(3,ii)-tk)
                    !              if(err1.lt.vt.and.err2.lt.vt.and.err3.lt.vt) then
                    IF(bign(1,ii).EQ.ti.AND.bign(2,ii).EQ.tj.AND.bign(3,ii).EQ.tk) THEN
                       indexplus(ig,inw)=ii
                       tagp(ig,inw)=-1
                       !                write (6,*) "Found +", ig,ii,inw
                       !               write (6,*) "looking for", ti,tj,tk
                       GO TO 214
                    ELSE
                    END IF
                 END DO
                 indexplus(ig,inw)=-1
                 tagp(ig,inw)=-1
                 !          write (6,*) "Not Found +", ig,-1,inw
                 !               write (6,*) "looking for", ti,tj,tk
              END IF
              !214        ti=(-bign(1,ig)+i_1(inw))*b1(1)+(-bign(2,ig)+j_1(inw))*b2(1)+(-bign(3,ig)+k_1(inw))*b3(1)
              !           tj=(-bign(1,ig)+i_1(inw))*b1(2)+(-bign(2,ig)+j_1(inw))*b2(2)+(-bign(3,ig)+k_1(inw))*b3(2)
              !           tk=(-bign(1,ig)+i_1(inw))*b1(3)+(-bign(2,ig)+j_1(inw))*b2(3)+(-bign(3,ig)+k_1(inw))*b3(3)
214           ti=(-bign(1,ig)+i_1(inw))
              tj=(-bign(2,ig)+j_1(inw))
              tk=(-bign(3,ig)+k_1(inw))
              ti1=-bign(1,ig)+i_1(inw)
              tj1=-bign(2,ig)+j_1(inw)
              tk1=-bign(3,ig)+k_1(inw)
              IF(ti1.LT.0.OR.(ti1.EQ.0.AND.(tj1.LT.0.OR.(tj1.EQ.0.AND.tk1.LT.0)))) THEN
                 DO ii=1,ntot
                    err1=ABS(bigg(1,ii)+ti)
                    err2=ABS(bigg(2,ii)+tj)
                    err3=ABS(bigg(3,ii)+tk)
                    !                    if(err1.lt.vt.and.err2.lt.vt.and.err3.lt.vt) then
                    IF(bign(1,ii).EQ.-ti.AND.bign(2,ii).EQ.-tj.AND.bign(3,ii).EQ.-tk) THEN
                       indexminus(ig,inw)=ii
                       tag(ig,inw)=1
                       !                     write (6,*) "Found -", ig,ii,inw 
                       !               write (6,*) "looking for", -ti,-tj,-tk
                       GO TO 213
                    ELSE
                    END IF
                 END DO
                 indexminus(ig,inw)=-1
                 tag(ig,inw)=1
                 !                write (6,*) "Not Found -", ig,-1,inw 
                 !               write (6,*) "looking for", -ti,-tj,-tk
              ELSE 
                 DO ii=1,ntot
                    err1=ABS(bigg(1,ii)-ti)
                    err2=ABS(bigg(2,ii)-tj)
                    err3=ABS(bigg(3,ii)-tk)
                    !                 if(err1.lt.vt.and.err2.lt.vt.and.err3.lt.vt) then
                    IF(bign(1,ii).EQ.ti.AND.bign(2,ii).EQ.tj.AND.bign(3,ii).EQ.tk) THEN
                       indexminus(ig,inw)=ii
                       tag(ig,inw)=-1
                       !                   write (6,*) "Found -", ig,ii,inw 
                       !               write (6,*) "looking for", ti,tj,tk
                       GO TO 213
                    ELSE
                    END IF
                 END DO
                 indexminus(ig,inw)=-1
                 tag(ig,inw)=-1
                 !              write (6,*) "Not Found -", ig,-1,inw 
                 !               write (6,*) "looking for", ti,tj,tk
              END IF
213           CONTINUE
           END DO
           WRITE( stdout, * ) "Translation", inw, "for", ntot, "G vectors"
#ifdef __PARA
        END IF
#endif
     END IF
  END DO

#ifdef __PARA

  CALL mp_barrier()
  !
  CALL mp_bcast( indexplus,  root )
  CALL mp_bcast( indexminus, root )
  CALL mp_bcast( tag,        root )
  CALL mp_bcast( tagp,       root )

  IF (me.EQ.1) THEN
#endif
     DEALLOCATE(bigg)
     DEALLOCATE(bign)
#ifdef __PARA
  END IF
#endif
  DEALLOCATE(i_1,j_1,k_1)

  RETURN
END SUBROUTINE wfunc_init
!
!----------------------------------------------------------------------------
SUBROUTINE grid_map()
  !----------------------------------------------------------------------------
  !
  USE kinds,                  ONLY : DP
  USE efcalc,                 ONLY : xdist, ydist, zdist
  USE smooth_grid_dimensions, ONLY : nnrsx, nr1s, nr2s, nr3s, &
                                     nr1sx, nr2sx, nr3sx
  USE fft_base,               ONLY : dffts
  USE mp_global,              ONLY : mpime
  USE parallel_include
  !
  IMPLICIT NONE
  !
  INTEGER :: ir1, ir2, ir3, ibig3, me
  !
  me = mpime + 1
  !
  ALLOCATE(xdist(nnrsx))
  ALLOCATE(ydist(nnrsx))
  ALLOCATE(zdist(nnrsx))
  !
  DO ir3=1,nr3s
#ifdef __PARA
     ibig3 = ir3 - dffts%ipp( me )
     IF(ibig3.GT.0.AND.ibig3.LE.dffts%npp(me)) THEN
#else
        ibig3=ir3
#endif
        DO ir2=1,nr2s
           DO ir1=1,nr1s
              xdist(ir1+(ir2-1)*nr1sx+(ibig3-1)*nr1sx*nr2sx) =                     &
                   &                  ((ir1-1)/DBLE(nr1sx))
              ydist(ir1+(ir2-1)*nr1sx+(ibig3-1)*nr1sx*nr2sx) =                   &
                   &                  ((ir2-1)/DBLE(nr2sx))
              zdist(ir1+(ir2-1)*nr1sx+(ibig3-1)*nr1sx*nr2sx) =                     &
                   &                  ((ir3-1)/DBLE(nr3sx))
              !         
           END DO
        END DO
#ifdef __PARA
     END IF
#endif
  END DO
  RETURN
END SUBROUTINE grid_map
!
!----------------------------------------------------------------------------
SUBROUTINE tric_wts( rp1, rp2, rp3, alat, wts )
  !----------------------------------------------------------------------------
  !
  ! ... This subroutine computes the weights to be used for
  ! ... R.P. translations in the WF calculation in the case
  ! ... of ibrav=0 or ibrav=14
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : pi
  USE cell_base, ONLY : tpiba, tpiba2
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: rp1(3), rp2(3), rp3(3)
  REAL(DP), INTENT(IN)  :: alat
  REAL(DP), INTENT(OUT) :: wts(6)
  ! 
  REAL(DP) :: b1x, b2x, b3x, b1y, b2y, b3y, b1z, b2z, b3z
  INTEGER        :: i
  !
  !
  b1x = rp1(1)*tpiba
  b2x = rp2(1)*tpiba
  b3x = rp3(1)*tpiba
  b1y = rp1(2)*tpiba
  b2y = rp2(2)*tpiba
  b3y = rp3(2)*tpiba
  b1z = rp1(3)*tpiba
  b2z = rp2(3)*tpiba
  b3z = rp3(3)*tpiba
  !        WRITE( stdout, * ) 'COMPUTING WEIGHTS NOW ...'

  wts(1) = tpiba2*(-b1z*b2x*b2z*b3x + b2y**2*b3x**2 + b1z*b2z*b3x**2 + & 
       b2z**2*b3x**2 - b1z*b2y*b2z*b3y - 2.D0*b2x*b2y*b3x*b3y + & 
       b2x**2*b3y**2 + b1z*b2z*b3y**2 + b2z**2*b3y**2 + & 
       b1z*b2x**2*b3z + b1z*b2y**2*b3z - b1z*b2x*b3x*b3z - & 
       2.D0*b2x*b2z*b3x*b3z - b1z*b2y*b3y*b3z - &
       2.D0*b2y*b2z*b3y*b3z + b2x**2*b3z**2 + b2y**2*b3z**2 + &
       b1x*(b2y**2*b3x + b2z**2*b3x - b2y*(b2x + b3x)*b3y - &
       b2z*(b2x + b3x)*b3z +  b2x*(b3y**2 + b3z**2)) + &
       b1y*(b2x**2*b3y - b2x*b3x*(b2y + b3y) + &
       b2z*b3y*(b2z - b3z) + b2y*(b3x**2 - b2z*b3z +  b3z**2)))/ &
       ((b1z*b2y*b3x - b1y*b2z*b3x - b1z*b2x*b3y + b1x*b2z*b3y &
       + b1y*b2x*b3z - b1x*b2y*b3z)**2)

  wts(2) = tpiba2*(b1z**2*(b2x*b3x + b3x**2 + b3y*(b2y + b3y)) + &
       b1y**2*(b2x*b3x + b3x**2 + b3z*(b2z + b3z)) - &
       b1z*(-b2z*(b3x**2 + b3y**2) + (b2x*b3x + b2y*b3y)*b3z + &
       b1x*(b2z*b3x + (b2x + 2.D0* b3x)*b3z)) - &
       b1y*(b1x*(b2y*b3x + (b2x + 2.D0*b3x)*b3y) + &
       b3y*(b1z*b2z + b2x*b3x + 2.D0*b1z*b3z + b2z*b3z) - &
       b2y*(b3x**2 - b1z* b3z + b3z**2)) + &
       b1x*(-b2y*b3x*b3y + b2x*b3y**2 - b2z*b3x*b3z + b2x*b3z**2 + &
       b1x*(b2y*b3y + b3y**2 + b3z*(b2z + b3z))))/ &                        
       ((b1z*b2y*b3x - b1y*b2z*b3x - b1z*b2x*b3y + b1x*b2z*b3y + &
       b1y*b2x*b3z - b1x*b2y*b3z)**2)

  wts(3) = tpiba2*(b1z**2*(b2x**2 + b2x*b3x + b2y*(b2y + b3y)) - &
       b1y*(2.D0*b1z*b2y*b2z + b2x*b2y*b3x - b2x**2*b3y + &
       b1z*b2z*b3y - b2z**2*b3y + b1x*(2.D0*b2x*b2y + b2y*b3x + b2x*b3y) + &
       b1z*b2y*b3z + b2y*b2z*b3z) + b1y**2*(b2x**2 + b2x*b3x + b2z*(b2z + b3z)) - &
       b1z*(b2x*b2z*b3x + b2y*b2z*b3y - b2x**2*b3z - b2y**2*b3z + &
       b1x*(2.D0*b2x*b2z + b2z*b3x + b2x*b3z)) + &
       b1x*(b2y**2*b3x + b2z**2*b3x - b2x*b2y*b3y - b2x*b2z*b3z + &
       b1x*(b2y**2 + b2y*b3y + b2z*(b2z +   b3z))))/ &
       ((b1z*b2y*b3x - b1y*b2z*b3x - b1z*b2x*b3y + b1x*b2z*b3y + b1y*b2x*b3z - &
       b1x*b2y*b3z)**2) 

  wts(4) = tpiba2*(b1z*(-b2z*(b3x**2 + b3y**2) + (b2x*b3x + b2y*b3y)*b3z) + & 
       b1y*(b3y*(b2x*b3x + b2z*b3z) - b2y*(b3x**2 + b3z**2)) + &
       b1x*(b2y*b3x*b3y + b2z*b3x*b3z - b2x*(b3y**2 +  b3z**2)))/ &
       ((b1z*b2y*b3x - b1y*b2z*b3x - b1z*b2x*b3y + b1x*b2z*b3y + &
       b1y*b2x*b3z - b1x*b2y*b3z)**2)

  wts(5) =  -tpiba2*(b1z**2*(b2x*b3x + b2y*b3y) - b1x*b1z*(b2z*b3x + b2x*b3z) - &
       b1y*(b1x*b2y*b3x + b1x*b2x*b3y + b1z*b2z*b3y + b1z*b2y*b3z) + &
       b1y**2*(b2x*b3x + b2z*b3z) + b1x**2*(b2y*b3y + b2z*b3z))/ &
       ((b1z*b2y*b3x - b1y*b2z*b3x - b1z*b2x*b3y + b1x*b2z*b3y + &
       b1y*b2x*b3z - b1x*b2y*b3z)**2)

  wts(6) = -tpiba2*(b1z*(-b2x*b2z*b3x - b2y*b2z*b3y + b2x**2*b3z + b2y**2*b3z) + &
       b1x*(b2y**2*b3x + b2z**2*b3x - b2x*b2y*b3y -  b2x*b2z*b3z) + &
       b1y*(-b2x*b2y*b3x + b2x**2*b3y + b2z*(b2z*b3y -  b2y*b3z)))/ &
       ((b1z*b2y*b3x - b1y*b2z*b3x - b1z*b2x*b3y + b1x*b2z*b3y + &
       b1y*b2x*b3z - b1x*b2y*b3z)**2)
  !
END SUBROUTINE tric_wts
!
!----------------------------------------------------------------------------
SUBROUTINE small_box_wf( i_1, j_1, k_1, nw1 )
  !----------------------------------------------------------------------------
  !
  USE kinds,                  ONLY : DP
  USE io_global,              ONLY : stdout
  USE constants,              ONLY : fpi
  USE wannier_base,           ONLY : expo
  USE grid_dimensions,        ONLY : nr1, nr2, nr3, nr1x, nr2x, nr3x, nnrx
  USE fft_base,               ONLY : dfftp
  USE mp_global,              ONLY : mpime
  USE parallel_include
  !
  IMPLICIT NONE

  INTEGER ir1, ir2, ir3, ibig3 , inw
  REAL(DP) x
  INTEGER , INTENT(in) :: nw1, i_1(nw1), j_1(nw1), k_1(nw1)
  INTEGER :: me

  me = mpime + 1

  ALLOCATE(expo(nnrx,nw1))

  DO inw=1,nw1

     WRITE( stdout, * ) inw ,":", i_1(inw), j_1(inw), k_1(inw)

     DO ir3=1,nr3
#ifdef __PARA
        ibig3 = ir3 - dfftp%ipp( me )
        IF(ibig3.GT.0.AND.ibig3.LE.dfftp%npp(me)) THEN
#else
           ibig3=ir3
#endif
           DO ir2=1,nr2
              DO ir1=1,nr1
                 x =  (((ir1-1)/DBLE(nr1x))*i_1(inw) +                          &
                      &                  ((ir2-1)/DBLE(nr2x))*j_1(inw) +             &
                      &                  ((ir3-1)/DBLE(nr3x))*k_1(inw))*0.5d0*fpi
                 expo(ir1+(ir2-1)*nr1x+(ibig3-1)*nr1x*nr2x,inw) =  CMPLX(COS(x), -SIN(x))
              END DO
           END DO
#ifdef __PARA
        END IF
#endif
     END DO
  END DO
  RETURN
END SUBROUTINE small_box_wf
!
!-----------------------------------------------------------------------
FUNCTION boxdotgridcplx(irb,qv,vr)
  !-----------------------------------------------------------------------
  !
  ! Calculate \sum_i qv(r_i)*vr(r_i)  with r_i on box grid
  ! array qv(r) is defined on box grid, array vr(r)on dense grid
  ! irb   : position of the box in the dense grid
  ! Parallel execution: remember to sum the contributions from other nodes
  !
  !      use ion_parameters
  !
  USE kinds,                    ONLY : DP
  USE grid_dimensions,          ONLY : nnrx, nr1, nr2, nr3, nr1x, nr2x, nr3x
  USE smallbox_grid_dimensions, ONLY : nnrbx, nr1b, nr2b, nr3b, &
                                       nr1bx, nr2bx, nr3bx
  USE fft_base,                 ONLY : dfftp
  USE mp_global,                ONLY : mpime
  !
  IMPLICIT NONE
  !
  INTEGER,           INTENT(IN):: irb(3)
  COMPLEX(DP), INTENT(IN):: qv(nnrbx), vr(nnrx)
  COMPLEX(DP)            :: boxdotgridcplx
  !
  INTEGER :: ir1, ir2, ir3, ir, ibig1, ibig2, ibig3, ibig, me
  !
  me = mpime + 1
  !
  boxdotgridcplx = ZERO

  DO ir3=1,nr3b
     ibig3=irb(3)+ir3-1
     ibig3=1+MOD(ibig3-1,nr3)
#ifdef __PARA
     ibig3 = ibig3 - dfftp%ipp( me )
     IF (ibig3.GT.0.AND.ibig3.LE.dfftp%npp(me)) THEN
#endif
        DO ir2=1,nr2b
           ibig2=irb(2)+ir2-1
           ibig2=1+MOD(ibig2-1,nr2)
           DO ir1=1,nr1b
              ibig1=irb(1)+ir1-1
              ibig1=1+MOD(ibig1-1,nr1)
              ibig=ibig1 + (ibig2-1)*nr1x + (ibig3-1)*nr1x*nr2x
              ir  =ir1 + (ir2-1)*nr1bx + (ir3-1)*nr1bx*nr2bx
              boxdotgridcplx = boxdotgridcplx + qv(ir)*vr(ibig)
           END DO
        END DO
#ifdef __PARA
     ENDIF
#endif
  END DO
  !
  RETURN
  !
END FUNCTION boxdotgridcplx
!
!----------------------------------------------------------------------------
SUBROUTINE write_rho_g( rhog )
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE io_global,          ONLY : stdout
  USE gvecp,              ONLY : ngm
  USE reciprocal_vectors, ONLY : gx, mill_l
  USE electrons_base,     ONLY : nspin
  USE fft_base,           ONLY : dfftp
  USE mp_global,          ONLY : nproc, mpime
  USE parallel_include
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) ,INTENT(IN) :: rhog(ngm,nspin) 
  REAL(DP),   ALLOCATABLE:: gnx(:,:), bigg(:,:)
  COMPLEX(DP),ALLOCATABLE :: bigrho(:)
  COMPLEX(DP) :: rhotmp_g(ngm)
  INTEGER           :: ntot, i, j, me
#ifdef __PARA
  INTEGER proc, ierr, root, ngdens(nproc),recvcount(nproc), displs(nproc)
  INTEGER recvcount2(nproc), displs2(nproc)
#endif
  CHARACTER (LEN=6)  :: name
  CHARACTER (LEN=15) :: name2

  me = mpime + 1

  ALLOCATE(gnx(3,ngm))


#ifdef __PARA
  root=0
#endif

  DO i=1,ngm
     gnx(1,i)=gx(1,i)
     gnx(2,i)=gx(2,i)
     gnx(3,i)=gx(3,i)
  END DO

#ifdef __PARA
  ntot=0
  DO i=1,nproc
     ngdens(i)=(dfftp%ngl(i)+1)/2
  END DO

  DO proc=1,nproc
     recvcount(proc)=ngdens(proc)*3
     recvcount2(proc)=ngdens(proc)
     IF(proc.EQ.1) THEN
        displs(proc)=0
        displs2(proc)=0
     ELSE
        displs(proc)=displs(proc-1)+recvcount(proc-1)
        displs2(proc)=displs2(proc-1)+recvcount2(proc-1)
     END IF
     ntot=ntot+recvcount(proc)/3
  END DO

  IF(me.EQ.1) THEN
     ALLOCATE(bigg(3,ntot))
  END IF

  CALL mpi_barrier(MPI_COMM_WORLD,ierr)
  IF(ierr.NE.0) CALL errore('write_rho_g0','mpi_barrier', ierr)
  CALL MPI_GATHERV(gnx,recvcount(me),MPI_REAL8,                        &
       bigg,recvcount,displs,MPI_REAL8,                    &
       root,MPI_COMM_WORLD, ierr)
  IF(ierr.NE.0) CALL errore('write_rho_g0','MPI_GATHERV', ierr)
  DO i=1,nspin

     rhotmp_g(1:ngm)=rhog(1:ngm,i)

     IF(me.EQ.1) THEN 
        ALLOCATE (bigrho(ntot))
     END IF

     CALL mpi_barrier(MPI_COMM_WORLD,ierr)
     IF(ierr.NE.0) CALL errore('write_rho_g1','mpi_barrier', ierr)
     CALL MPI_GATHERV(rhotmp_g,recvcount2(me),MPI_DOUBLE_COMPLEX,                &
          bigrho,recvcount2,displs2,MPI_DOUBLE_COMPLEX,              &
          root,MPI_COMM_WORLD, ierr)
     IF(ierr.NE.0) CALL errore('write_rho_g1','MPI_GATHERV', ierr)


     IF(me.EQ.1) THEN
        IF(i.EQ.1) name2="CH_DEN_G_PARA.1"
        IF(i.EQ.2) name2="CH_DEN_G_PARA.2"
        OPEN(unit=57, file=name2) 
        DO j=1,ntot
           WRITE(57,*) bigrho(j)
        END DO
        CLOSE(57)
        DEALLOCATE(bigrho)
     END IF

     WRITE( stdout, * ) "Charge density written to ", name2

  END DO

  IF(me.EQ.1) THEN
     name="G_PARA"
     OPEN(unit=56, file=name) 
     DO i=1,ntot
        WRITE(56,*) bigg(:,i)
     END DO
     CLOSE(56)
     DEALLOCATE(bigg)
  END IF
  WRITE( stdout, * ) "G-vectors written to G_PARA"
#else
  ntot=ngm
  ALLOCATE(bigg(3,ntot))
  bigg(1:3,1:ntot)=gnx(1:3,1:ngm)
  DO i=1,nspin
     ALLOCATE(bigrho(ntot))
     bigrho(1:ngm)=rhog(1:ngm,i)

     IF(i.EQ.1) name2="CH_DEN_G_SERL.1"
     IF(i.EQ.2) name2="CH_DEN_G_SERL.2"

     OPEN(unit=57, file=name2) 
     DO j=1,ntot
        WRITE(57,*) bigrho(j)
     END DO
     CLOSE(57)
     DEALLOCATE(bigrho)

     WRITE( stdout, * ) "Charge density written to", name2

  END DO

  name="G_SERL"
  OPEN(unit=56, file=name) 
  DO i=1,ntot
     WRITE(56,*) bigg(:,i)
  END DO
  CLOSE(56)
  DEALLOCATE(bigg)
  WRITE( stdout, * ) "G-vectors written to G_SERL"
#endif
  !
  DEALLOCATE(gnx)
  !
  RETURN
  !
END SUBROUTINE write_rho_g
!
!----------------------------------------------------------------------------
SUBROUTINE macroscopic_average( rhog, tau0, e_tuned )
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE reciprocal_vectors, ONLY : gx
  USE gvecp,              ONLY : ngm
  USE electrons_base,     ONLY : nspin
  USE tune,               ONLY : npts, xdir, ydir, zdir, B, &
                                 shift, start, av0, av1
  USE cell_base,          ONLY : a1, a2, a3, tpiba, omega
  USE ions_base,          ONLY : nsp, na, zv, nsx, nax
  USE constants,          ONLY : pi, tpi
  USE mp,                 ONLY : mp_barrier, mp_bcast
  USE fft_base,           ONLY : dfftp
  USE mp_global,          ONLY : nproc, mpime
  USE parallel_include
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE:: gnx(:,:), bigg(:,:)
  COMPLEX(DP) ,INTENT(in) :: rhog(ngm,nspin)
  COMPLEX(DP),ALLOCATABLE :: bigrho(:)
  COMPLEX(DP) :: rhotmp_g(ngm)
  INTEGER ntot, i, j, ngz, l, isa
  INTEGER ,ALLOCATABLE :: g_red(:,:)
#ifdef __PARA
  INTEGER proc, ierr, root, ngdens(nproc),recvcount(nproc), displs(nproc)
  INTEGER recvcount2(nproc), displs2(nproc)
#endif
  REAL(DP) zlen,vtot, pos(3,nax,nsx), a_direct(3,3),a_trans(3,3), e_slp, e_int
  REAL(DP), INTENT(out) :: e_tuned(3)
  REAL(DP), INTENT(in) :: tau0(3,nax)
  REAL(DP),ALLOCATABLE :: v_mr(:), dz(:), gz(:), g_1(:,:), vbar(:), cd(:), v_final(:)
  REAL(DP), ALLOCATABLE:: cdion(:), cdel(:), v_line(:), dist(:)
  COMPLEX(DP),ALLOCATABLE :: rho_ion(:),v_1(:),vmac(:),rho_tot(:),rhogz(:), bigrhog(:)
  INTEGER :: me

  me = mpime + 1

  ALLOCATE(gnx(3,ngm))


#ifdef __PARA
  root=0
#endif

  DO i=1,ngm
     gnx(1,i)=gx(1,i)
     gnx(2,i)=gx(2,i)
     gnx(3,i)=gx(3,i)
  END DO

#ifdef __PARA
  ntot=0
  DO i=1,nproc
     ngdens(i)=(dfftp%ngl(i)+1)/2
  END DO

  DO proc=1,nproc
     recvcount(proc)=ngdens(proc)*3
     recvcount2(proc)=ngdens(proc)
     IF(proc.EQ.1) THEN
        displs(proc)=0
        displs2(proc)=0
     ELSE
        displs(proc)=displs(proc-1)+recvcount(proc-1)
        displs2(proc)=displs2(proc-1)+recvcount2(proc-1)
     END IF
     ntot=ntot+recvcount(proc)/3
  END DO

  ALLOCATE(bigg(3,ntot))
  ALLOCATE(g_1(3,2*ntot-1))
  ALLOCATE(g_red(3,2*ntot-1))

  CALL mp_barrier()
  !
  CALL MPI_GATHERV(gnx,recvcount(me),MPI_REAL8,                        &
       bigg,recvcount,displs,MPI_REAL8,                    &
       root,MPI_COMM_WORLD, ierr)
  IF(ierr.NE.0) CALL errore('macroscopic_avergae','MPI_GATHERV', ierr)
  !
  CALL mpi_bcast( bigg )
  !
  rhotmp_g(1:ngm)=rhog(1:ngm,1)

  ALLOCATE (bigrho(ntot))
  ALLOCATE (bigrhog(2*ntot-1))

  CALL mp_barrier()
  !
  CALL MPI_GATHERV(rhotmp_g,recvcount2(me),MPI_DOUBLE_COMPLEX,       &
       bigrho,recvcount2,displs2,MPI_DOUBLE_COMPLEX,       &
       root,MPI_COMM_WORLD, ierr)
  IF(ierr.NE.0) CALL errore('macroscopic_avergae','MPI_GATHERV', ierr)
  !
  CALL mp_bcast( bigrho, root )
  !
#else
  ntot=ngm
  ALLOCATE(bigg(3,ntot))
  ALLOCATE(g_1(3,2*ntot-1))
  ALLOCATE(g_red(3,2*ntot-1))
  bigg(1:3,1:ntot)=gnx(1:3,1:ngm)
  ALLOCATE(bigrho(ntot))
  ALLOCATE(bigrhog(2*ntot-1))
  bigrho(1:ngm)=rhog(1:ngm,1)
#endif

  ALLOCATE(v_mr(npts))
  ALLOCATE(v_final(npts))
  ALLOCATE(dz(npts))
  ALLOCATE(vbar(npts))
  ALLOCATE(cd(npts))
  ALLOCATE(cdel(npts))
  ALLOCATE(cdion(npts))

  !-- needed for non-orthogonal cells

  a_direct(1,1:3)=a1(1:3)
  a_direct(2,1:3)=a2(1:3)
  a_direct(3,1:3)=a3(1:3)

  a_trans=TRANSPOSE(a_direct)

  !--- Construct rho(-g) from rho(g). rgo(-g)=rho*(g)

  bigrhog(1:ntot)=bigrho(1:ntot)
  g_1(:,1:ntot)=bigg(:,1:ntot)
  DO i=2,ntot
     bigrhog(ntot+i-1)=CONJG(bigrho(i))
     g_1(:,ntot+i-1)=-bigg(:,i)
  END DO

  !--- needed fot non-orthogonal cells

  DO i=1,2*ntot-1
     g_red(:,i)=NINT(MATMUL(a_trans(:,:),g_1(:,i))*tpiba/tpi)
  END DO

  !--- define the direction of the line

  xdir=1
  ydir=2

  IF ((zdir).EQ.1) xdir=3
  IF ((zdir).EQ.2) ydir=3

  IF(zdir.EQ.1) zlen=DSQRT(a1(1)**2+a1(2)**2+a1(3)**2)
  IF(zdir.EQ.2) zlen=DSQRT(a2(1)**2+a2(2)**2+a2(3)**2)
  IF(zdir.EQ.3) zlen=DSQRT(a3(1)**2+a3(2)**2+a3(3)**2)


  !--- We need the potentiail only along zdir, so pick the appropriate G-vectors with Gxdir=Gydir=0

  ngz=0
  DO i=1,2*ntot-1
     IF((g_red(xdir,i).EQ.0).AND.(g_red(ydir,i).EQ.0)) ngz=ngz+1
  END DO

  ALLOCATE(gz(ngz))
  ALLOCATE(rhogz(ngz))
  ALLOCATE(rho_ion(ngz))
  ALLOCATE(rho_tot(ngz))
  ALLOCATE(vmac(ngz))
  ALLOCATE(v_1(ngz))

  !--- The G-vectors are output in units of 2*pi/a, so convert them to the correct values

  j=0
  DO i=1,2*ntot-1
     IF((g_red(xdir,i).EQ.0).AND.(g_red(ydir,i).EQ.0)) THEN
        j=j+1
        gz(j)=g_1(zdir,i)*tpiba
        rhogz(j)=bigrhog(i)
     END IF
  END DO

  isa = 0
  DO i=1,nsp
     DO j=1,na(i)
        isa = isa + 1
        pos(:,j,i)=tau0(:,isa)
     END DO
  END DO

  !--- Construct the ionic Charge density in G-space

  rho_ion = ZERO
  !
  DO j=1,ngz
     DO i=1,nsp
        DO l=1,na(i)
           rho_ion(j)=rho_ion(j)+zv(i)*EXP(-CI*gz(j)*pos(zdir,l,i))*EXP(-gz(j)**2/(4.D0*ONE))
        END DO
     END DO
  END DO

  rho_ion=rho_ion/omega

  !--- Construct the total Charge density in G-space

  rho_tot=rho_ion-rhogz

  !--- Construct the electrostatic potential and macroscopic average in G-space

  v_1(1)=ZERO
  vmac(1)=ZERO
  v_1(2:ngz)=4*pi*rho_tot(2:ngz)/gz(2:ngz)**2
  vmac(2:)=v_1(2:)*SIN(gz(2:)*b)/(gz(2:)*b)


  !--- Calculate planewise average in R-space and FFT V(Gz) ---> V(z) ... well not really FFT but FT

  vbar=0.D0
  v_mr=0.D0
  cdel=0.D0
  cdion=0.D0
  cd=0.D0
  DO j=1,npts
     dz(j)=(j-1)*zlen/(npts*1.D0)
     DO i=1,ngz
        vbar(j)=vbar(j)-DBLE(EXP(CI*gz(i)*dz(j))*v_1(i))
        v_mr(j)=v_mr(j)-DBLE(EXP(CI*gz(i)*dz(j))*vmac(i))
        cdel(j)=cdel(j)-DBLE(EXP(CI*gz(i)*dz(j))*rhogz(i))
        cdion(j)=cdion(j)+DBLE(EXP(CI*gz(i)*dz(j))*rho_ion(i))
        cd(j)=cd(j)+DBLE(EXP(CI*gz(i)*dz(j))*rho_tot(i))
     END DO
     !           WRITE( stdout, * ) vbar(j), v_mr(j), cdel(j), cdion(j)
  END DO
  IF (shift) THEN
     vtot=(v_mr(start)+v_mr(start-1))/2.D0
     v_final(1:npts-start+1)=v_mr(start:npts)-vtot
     v_final(npts-start+2:npts)=v_mr(1:start-1)-vtot
  ELSE
     vtot=(v_mr(1)+v_mr(npts))/2.D0
     v_final(1:npts)=v_mr(1:npts)-vtot
  END IF

  e_tuned=0.D0

  ALLOCATE(v_line(1:av1-av0+1))
  ALLOCATE(dist(1:av1-av0+1))


  v_line(1:av1-av0+1)=v_final(av0:av1)
  dist(1:av1-av0+1) =dz(av0:av1)

  e_tuned(zdir)=-(v_final(av1)-v_final(av0))/((av1-av0)*zlen/(npts*1.D0))

#ifdef __PARA

  DEALLOCATE(bigg,g_1,bigrho,bigrhog,g_red)     

#else
  DEALLOCATE(bigg,g_1,bigrho,bigrhog,g_red)     
#endif

  DEALLOCATE(gnx,v_mr,v_final,dz,vbar,cd,cdel,cdion)
  DEALLOCATE(v_line, dist)
  DEALLOCATE(gz,rhogz,rho_ion,rho_tot,vmac,v_1)

  RETURN
END SUBROUTINE macroscopic_average
!
!----------------------------------------------------------------------------
SUBROUTINE least_square( npts, x, y, slope, intercept )
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,        INTENT(IN) :: npts
  REAL(DP), INTENT(IN) :: x(npts), y(npts)
  REAL(DP), INTENT(OUT):: slope, intercept
  !
  INTEGER        :: i
  REAL(DP) :: sumx,sumy,sumx2,sumxy,sumsqx
  REAL(DP) :: xav,yav

  sumxy=0.D0
  sumx =0.D0
  sumy =0.D0
  sumx2=0.D0
  DO i=1,npts
     sumxy=sumxy+x(i)*y(i)
     sumx =sumx +x(i)
     sumy =sumy +y(i)
     sumx2=sumx2+x(i)*x(i)
  END DO
  sumsqx=sumx**2
  xav=sumx/DBLE(npts)
  yav=sumy/DBLE(npts)

  slope=(npts*sumxy - sumx*sumy)/(npts*sumx2 - sumsqx)

  intercept=yav-slope*xav

  RETURN

END SUBROUTINE least_square
!
!----------------------------------------------------------------------------
SUBROUTINE wfsteep( m, Omat, Umat, b1, b2, b3 )
  !----------------------------------------------------------------------------
  !
  USE kinds,                  ONLY : DP
  USE io_global,              ONLY : stdout
  USE wannier_base,           ONLY : nw, weight, nit, tolw, wfdt, maxwfdt, nsd
  USE control_flags,          ONLY : iprsta
  USE cell_base,              ONLY : alat
  USE constants,              ONLY : tpi
  USE smooth_grid_dimensions, ONLY : nr1s, nr2s, nr3s
  USE mp_global,              ONLY : mpime
  USE parallel_include
  !
  IMPLICIT NONE

  !    (m,m) is the size of the matrix Ospin.
  !    Ospin is input overlap matrix.
  !    Uspin is the output unitary transformation.
  !             Rough guess for Uspin can be carried in.
  !
  !     conjugated gradient to search maximization
  !
  REAL(DP), PARAMETER :: autoaf=0.529177d0
  INTEGER, INTENT(in) :: m
  REAL(DP), INTENT(in) :: b1(3),b2(3),b3(3)
  COMPLEX(DP), INTENT(inout) :: Omat(nw, m, m)
  REAL(DP), INTENT(inout) :: Umat(m,m)
  !
  INTEGER :: i, j, k, l, ig, ierr, ti, tj, tk, inw, ir, adjust
  INTEGER :: f3(nw), f4(nw), ierr1
  REAL(DP) :: slope, slope2, t1, t2, t3, mt(nw),t21,temp1,maxdt
  REAL(DP) :: U(m,m), wfc(3, m), Wm(m,m), schd(m,m), f2(4*m), gr(nw, 3)
  REAL(DP) :: Uspin2(m,m),temp2,wfdtold,oldt1,t01, d3(m,m), d4(m,m), U1(m,m)
  REAL(DP) :: spread, sp
  REAL(DP), ALLOCATABLE  :: wr(:)
  REAL(DP), ALLOCATABLE  :: W(:,:)
  COMPLEX(DP) :: ct1, ct2, ct3, z(m, m), X(m, m), d(m,m), d2(m,m)
  COMPLEX(DP) :: f1(2*m-1), wp(m*(m+1)/2), Oc(nw, m, m)
  COMPLEX(DP) ::  Oc2(nw, m, m),wp1(m*(m+1)/2), X1(m,m), U2(m,m), U3(m,m)
  INTEGER :: me
  !
  me = mpime + 1
  !
  ALLOCATE(W(m,m), wr(m))
  !
  Umat=0.D0
  DO i=1,m
     Umat(i,i)=1.D0
  END DO
  Oc=ZERO
  Oc2=ZERO
  X1=ZERO
  U2=Umat*ONE
  !
  ! update Oc using the initial guess of Uspin
  !
  DO inw=1, nw
     X1(:, :)=Omat(inw, :, :)
     U3=ZERO
     CALL ZGEMM ('T', 'nbsp', m,m,m,ONE,U2,m,X1,m,ZERO,U3,m)
     X1=ZERO
     CALL ZGEMM ('nbsp','nbsp', m,m,m,ONE,U3,m,U2,m,ZERO,X1,m)
     Oc(inw, :, :)=X1(:, :)
  END DO

  U2=ZERO
  U3=ZERO

  W=0.D0
  schd=0.D0
  oldt1=0.D0
  wfdtold=0.D0

  DO k=1, nit
     t01=0.D0     !use t1 to store the value of omiga
     DO inw=1, nw
        DO i=1, m
           t01=t01+DBLE(CONJG(Oc(inw, i, i))*Oc(inw, i, i))
        END DO
     END DO

     !    WRITE( stdout, * ) t01

     IF(ABS(oldt1-t01).LT.tolw) THEN 
#ifdef __PARA
        IF(me.EQ.1) THEN
#endif
           WRITE(27,*) "MLWF Generated at Step",k
#ifdef __PARA
        END IF
#endif
        IF(iprsta.GT.4) THEN
           WRITE( stdout, * ) "MLWF Generated at Step",k
        END IF
        GO TO 40
     END IF

     !    oldt1=t01

     !   calculate d(omiga)/dW and store result in W
     !   W should be a real symmetric matrix for gamma-point calculation
     !
     Wm=W
     W=0.D0
     DO inw=1, nw
        t2=weight(inw)
        DO i=1,m
           DO j=i+1,m
              W(i,j)=W(i,j)+t2*DBLE(Oc(inw,i,j)*CONJG(Oc(inw,i,i)        &
                   -Oc(inw,j,j))+CONJG(Oc(inw,j,i))*(Oc(inw,i,i)-Oc(inw,j,j)))
           END DO
        END DO
     END DO
     W=W-TRANSPOSE(W)

     !   calculate slope=d(omiga)/d(lamda)
     slope=SUM(W**2)

     !   calculate slope2=d2(omiga)/d(lamda)2
     slope2=0.D0
     DO ti=1, m
        DO tj=1, m
           DO tk=1, m
              t2=0.D0
              DO inw=1, nw
                 t2=t2+DBLE(Oc(inw,tj,tk)*CONJG(Oc(inw,tj,tj)+Oc(inw,tk,tk) &
                      -2.D0*Oc(inw,ti,ti))-4.D0*Oc(inw,ti,tk)          &
                      *CONJG(Oc(inw,ti,tj)))*weight(inw)
              END DO
              slope2=slope2+W(tk,ti)*W(ti,tj)*2.D0*t2
           END DO
        END DO
     END DO
     slope2=2.D0*slope2

     !   use parabola approximation. Defined by 1 point and 2 slopes
     IF (slope2.LT.0) wfdt=-slope/2.D0/slope2
     IF (maxwfdt.GT.0.AND.wfdt.GT.maxwfdt) wfdt=maxwfdt

     IF (k.LT.nsd) THEN
        schd=W    !use steepest-descent technique

        !   calculate slope=d(omiga)/d(lamda)
        slope=SUM(schd**2)

        !       schd=schd*maxwfdt
        DO i=1, m
           DO j=i, m
              wp1(i + (j-1)*j/2) = CMPLX(0.d0, schd(i,j))
           END DO
        END DO

#if defined (__AIX)
        !
        CALL zhpev(21, wp1, wr, z, m, m, f2, 4*m)
        !
        ierr1 = 0
        !
#else   
        !    
        CALL zhpev('V','U',m,wp1,wr,z,m,f1,f2,ierr)
        !  
#endif

        IF (ierr.NE.0) STOP 'failed to diagonalize W!'

     ELSE
        !
        CALL DGEMM ('T','nbsp', m,m,m,ONE,Wm,m,Wm,m,ZERO,d3,m)

        t1=0.D0
        DO i=1, m
           t1=t1+d3(i, i)
        END DO
        IF (t1.NE.0) THEN
           d4=(W-Wm)
           CALL DGEMM ('T','nbsp', m,m,m,ONE,W,m,d4,m,ZERO,d3,m)
           t2=0.D0
           DO i=1, m
              t2=t2+d3(i, i)
           END DO
           t3=t2/t1
           schd=W+schd*t3
        ELSE
           schd=W
        END IF
        !
        !   calculate the new d(Lambda) for the new Search Direction
        !   added by Manu. September 19, 2001
        !
        !   calculate slope=d(omiga)/d(lamda)
        slope=SUM(schd**2)
        !------------------------------------------------------------------------
        !   schd=schd*maxwfdt
        DO i=1, m
           DO j=i, m
              wp1(i + (j-1)*j/2) = CMPLX(0.d0, schd(i,j))
           END DO
        END DO

#if defined __AIX
        CALL zhpev(21, wp1, wr, z, m, m, f2, 4*m)
        ierr1 = 0
#else
        CALL zhpev('V','U',m,wp1,wr,z,m,f1,f2,ierr)
#endif
        IF (ierr.NE.0) STOP 'failed to diagonalize W!'

        maxdt=maxwfdt

11      d=0.D0
        DO i=1, m
           d(i, i)=EXP(CI*(maxwfdt)*wr(i))
        END DO

        U3=ZERO
        CALL ZGEMM ('nbsp', 'nbsp', m,m,m,ONE,z,m,d,m,ZERO,U3,m)
        U2=ZERO
        CALL ZGEMM ('nbsp','c', m,m,m,ONE,U3,m,z,m,ZERO,U2,m)
        U=DBLE(U2)
        U2=ZERO
        U3=ZERO
        !
        !   update Uspin
        U1=ZERO
        CALL DGEMM ('nbsp', 'nbsp', m,m,m,ONE,Umat,m,U,m,ZERO,U1,m)
        Umat=U1

        !
        !   update Oc
        !
        U2=Umat*ONE
        U3=ZERO
        DO inw=1, nw
           X1(:,:)=Omat(inw,:,:)
           CALL ZGEMM ('T', 'nbsp', m,m,m,ONE,U2,m,X1,m,ZERO,U3,m)
           X1=ZERO
           CALL ZGEMM ('nbsp','nbsp',m,m,m,ONE,U3,m,U2,m,ZERO,X1,m)
           Oc2(inw, :,:)=X(:,:)
        END DO
        U2=ZERO 
        U3=ZERO
        !
        t21=0.D0     !use t21 to store the value of omiga
        DO inw=1, nw
           DO i=1, m
              t21=t21+DBLE(CONJG(Oc2(inw, i, i))*Oc2(inw, i, i))
           END DO
        END DO

        temp1=-((t01-t21)+slope*maxwfdt)/(maxwfdt**2)
        temp2=slope
        wfdt=-temp2/(2*temp1)

        IF (wfdt.GT.maxwfdt.OR.wfdt.LT.0.D0) THEN
           maxwfdt=2*maxwfdt
           GO TO 11
        END IF

        maxwfdt=maxdt
        !
        !
        !   use parabola approximation. Defined by 2 point and 1 slopes
        !    if (slope2.lt.0) wfdt=-slope/2.D0/slope2
        !    if (maxwfdt.gt.0.and.wfdt.gt.maxwfdt) wfdt=maxwfdt
        !
        !    write(6, '(e12.5E2,1x,e11.5E2,1x,f6.2)') slope2, slope, wfdt
        !-------------------------------------------------------------------------
        !
        !      schd is the new searching direction
        !
     END IF

     d=0.D0
     DO i=1, m
        d(i, i)=EXP(CI*wfdt*wr(i))
     END DO          !d=exp(d)


     !   U=z*exp(d)*z+
     !
     U3=ZERO
     CALL ZGEMM ('nbsp', 'nbsp', m,m,m,ONE,z,m,d,m,ZERO,U3,m)
     U2=ZERO
     CALL ZGEMM ('nbsp','c', m,m,m,ONE,U3,m,z,m,ZERO,U2,m)
     U=DBLE(U2)
     U2=ZERO
     U3=ZERO

     !   update Uspin
     !
     U1=ZERO
     CALL DGEMM ('nbsp', 'nbsp', m,m,m,ONE,Umat,m,U,m,ZERO,U1,m)
     Umat=U1

     !   update Oc
     !
     U2=Umat*ONE
     U3=ZERO
     DO inw=1, nw
        X1(:, :)=Omat(inw, :, :)
        CALL ZGEMM ('T', 'nbsp', m,m,m,ONE,U2,m,X1,m,ZERO,U3,m)
        X1=ZERO
        CALL ZGEMM ('nbsp','nbsp',m,m,m,ONE,U3,m,U2,m,ZERO,X1,m)
        Oc(inw, :, :)=X1(:, :)
     END DO
     U2=ZERO
     U3=ZERO
     IF(ABS(t01-oldt1).GE.tolw.AND.k.GE.nit) THEN
#ifdef __PARA
        IF(me.EQ.1) THEN
#endif
           WRITE(27,*) "MLWF Not generated after",k,"Steps."
#ifdef __PARA
        END IF
#endif
        IF(iprsta.GT.4) THEN
           WRITE( stdout, * ) "MLWF Not generated after",k,"Steps."
        END IF
        GO TO 40
     END IF
     oldt1=t01
  END DO

40 DEALLOCATE(W, wr)

  !
  ! calculate the spread
  !
  !  write(24, *) "spread: (unit \AA^2)"
  DO i=1, m
     mt=1.D0-DBLE(Oc(:,i,i)*CONJG(Oc(:,i,i)))
     sp = (alat*autoaf/tpi)**2*SUM(mt*weight)
#ifdef __PARA
     IF(me.EQ.1) THEN
#endif
        WRITE(25, '(f10.7)') sp
#ifdef __PARA
     END IF
#endif
     IF( sp < 0.D0 ) &
        CALL errore( 'cp-wf', 'Something wrong WF Spread negative', 1 )
     !
     spread=spread+sp
     !
  END DO
  spread=spread/DBLE(m)

#ifdef __PARA
  IF(me.EQ.1) THEN
#endif
     WRITE(24, '(f10.7)') spread
     WRITE(27,*) "Average spread = ", spread
#ifdef __PARA
  END IF
#endif
  !
  RETURN
END SUBROUTINE wfsteep
!
!----------------------------------------------------------------------------
SUBROUTINE dforce_field( bec, deeq, betae, i, c, ca, df, da, v, v1 )
  !----------------------------------------------------------------------------
  !
  ! ... computes: the generalized force df=CMPLX(dfr,dfi) acting on the i-th
  ! ...           electron state at the gamma point of the brillouin zone
  ! ...           represented by the vector c=CMPLX(cr,CI)
  !
  ! ...    d_n(g) = f_n { 0.5 g^2 c_n(g) + [vc_n](g) +
  ! ...                   sum_i,ij d^q_i,ij (-i)**l beta_i,i(g) 
  ! ...                                e^-ig.r_i < beta_i,j | c_n > }
  !
  USE kinds,                  ONLY : DP
  USE control_flags,          ONLY : iprint, tbuff
  USE gvecs,                  ONLY : nms, nps
  USE gvecw,                  ONLY : ngw
  USE cvan,                   ONLY : ish
  USE smooth_grid_dimensions, ONLY : nr1s, nr2s, nr3s, &
                                     nr1sx, nr2sx, nr3sx, nnrsx
  USE electrons_base,         ONLY : nbspx, nbsp, nspin, f, ispin
  USE constants,              ONLY : pi, fpi
  USE ions_base,              ONLY : nsp, na, nat
  USE gvecw,                  ONLY : ggp
  USE uspp_param,             ONLY : nh, nhm
  USE uspp,                   ONLY : nkb, dvan
  USE cell_base,              ONLY : tpiba2
  USE fft_module,             ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: betae(ngw,nkb), c(ngw), ca(ngw), df(ngw), da(ngw)
  REAL(DP)    :: bec(nkb,nbsp), deeq(nhm,nhm,nat,nspin), v(nnrsx,nspin), v1(nnrsx,nspin)
  INTEGER           :: i
  ! local variables
  INTEGER             :: iv, jv, ia, is, isa, ism, ios, iss1, iss2, ir, ig, inl, jnl
  REAL(DP)      :: fi, fip, dd
  COMPLEX(DP)   :: fp,fm
  REAL(DP)      :: af(nkb), aa(nkb)
  COMPLEX(DP)   :: dtemp(ngw)    !
  COMPLEX(DP), ALLOCATABLE :: psi(:)
  !
  !
  CALL start_clock( 'dforce_field' )
  ALLOCATE( psi( nnrsx ) )
  !
  !     important: if nbsp is odd => c(*,nbsp+1)=0.
  ! 
  IF (MOD(nbsp,2).NE.0.AND.i.EQ.nbsp) THEN
     DO ig=1,ngw
        ca(ig)=(0.,0.)
     END DO
  ENDIF
  !
  IF (.NOT.tbuff) THEN
     !
     psi( 1:nnrsx ) = 0.D0
     DO ig=1,ngw
        psi(nms(ig))=CONJG(c(ig)-CI*ca(ig))
        psi(nps(ig))=c(ig)+CI*ca(ig)
     END DO
     !
     CALL invfft('Wave',psi,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
     !     
  ELSE
     !
     !     read psi from buffer 21
     !
     READ(21,iostat=ios) psi
     IF(ios.NE.0) CALL errore( 'dforce', 'error in reading unit 21', ios )
     !
  ENDIF
  ! 
  iss1=ispin(i)
  !
  ! the following avoids a potential out-of-bounds error
  !
  IF (i.NE.nbsp) THEN
     iss2=ispin(i+1)
  ELSE
     iss2=iss1
  END IF
  !
  DO ir=1,nnrsx
     psi(ir)=CMPLX(v(ir,iss1)* DBLE(psi(ir)), v1(ir,iss2)*AIMAG(psi(ir)) )
  END DO
  !
  CALL fwfft('Wave',psi,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
  !
  !     note : the factor 0.5 appears 
  !       in the kinetic energy because it is defined as 0.5*g**2
  !       in the potential part because of the logics
  !
  fi =-  f(i)*0.5
  fip=-f(i+1)*0.5
  DO ig=1,ngw
     fp= psi(nps(ig)) + psi(nms(ig))
     fm= psi(nps(ig)) - psi(nms(ig))
     df(ig)= fi*(tpiba2*ggp(ig)* c(ig)+CMPLX(DBLE(fp), AIMAG(fm)))
     da(ig)=fip*(tpiba2*ggp(ig)*ca(ig)+CMPLX(AIMAG(fp),-DBLE(fm)))
  END DO
  !
  !     aa_i,i,nbsp = sum_j d_i,ij <beta_i,j|c_n>
  ! 
  IF(nkb.GT.0)THEN
     DO inl=1,nkb
        af(inl)=0.
        aa(inl)=0.
     END DO
     !
     DO is=1,nsp
        DO iv=1,nh(is)
           DO jv=1,nh(is)
              isa=0
              DO ism=1,is-1
                 isa=isa+na(ism)
              END DO
              DO ia=1,na(is)
                 inl=ish(is)+(iv-1)*na(is)+ia
                 jnl=ish(is)+(jv-1)*na(is)+ia
                 isa=isa+1
                 dd = deeq(iv,jv,isa,iss1)+dvan(iv,jv,is)
                 af(inl)=af(inl)-  f(i)*dd*bec(jnl,  i)
                 dd = deeq(iv,jv,isa,iss2)+dvan(iv,jv,is)
                 IF (i.NE.nbsp) aa(inl)=aa(inl)-f(i+1)*dd*bec(jnl,i+1)
              END DO
           END DO
        END DO
     END DO
     !
     DO ig=1,ngw
        dtemp(ig)=(0.,0.)
     END DO
     CALL MXMA(betae,1,2*ngw,af,1,nkb,dtemp,1,2*ngw,2*ngw,nkb,1)
     DO ig=1,ngw
        df(ig)=df(ig)+dtemp(ig)
     END DO
     !
     DO ig=1,ngw
        dtemp(ig)=(0.,0.)
     END DO
     CALL MXMA(betae,1,2*ngw,aa,1,nkb,dtemp,1,2*ngw,2*ngw,nkb,1)
     DO ig=1,ngw
        da(ig)=da(ig)+dtemp(ig)
     END DO
  ENDIF

  DEALLOCATE( psi )
  !
  CALL stop_clock( 'dforce_field' )
  !
  RETURN
END SUBROUTINE dforce_field
!
#if defined (__PARA)
!
!----------------------------------------------------------------------------
SUBROUTINE write_psi( c, jw )
  !----------------------------------------------------------------------------
  ! ... for calwf 5             - M.S
  ! ... collect wavefunctions on first node and write to file
  !
  USE kinds,                  ONLY : DP
  USE io_global,              ONLY : stdout
  USE gvecs,                  ONLY : nps
  USE electrons_base,         ONLY : nbspx
  USE smooth_grid_dimensions, ONLY : nnrsx, nr3sx, nr1s, nr2s, nr3s
  USE gvecw ,                 ONLY : ngw
  USE reciprocal_vectors,     ONLY : mill_l
  USE mp,                     ONLY : mp_barrier
  USE fft_base,               ONLY : dfftp
  USE mp_global,              ONLY : nproc, mpime
  USE parallel_include
  !
  IMPLICIT NONE
  !
  INTEGER :: unit, jw
  COMPLEX(DP) :: c(ngw,nbspx)
  COMPLEX(DP), ALLOCATABLE :: psis(:)
  !
  INTEGER ::i, ii, ig, proc, ierr, ntot, ncol, mc,ngpwpp(nproc)
  INTEGER ::nmin(3), nmax(3), n1,n2,nzx,nz,nz_
  INTEGER ::root, displs(nproc), recvcount(nproc)
  COMPLEX(DP), ALLOCATABLE:: psitot(:), psiwr(:,:,:)
  INTEGER :: me

  me = mpime + 1
  !
  ! nmin, nmax are the bounds on (i,j,k) indexes of wavefunction G-vectors
  !
  CALL nrbounds( ngw, nr1s, nr2s, nr3s, mill_l, nmin, nmax )
  !
  ! nzx is the maximum length of a column along z
  !
  nzx=nmax(3)-nmin(3)+1
  !
  root = 0
  ! root is the first node
  ntot = 0
  !
  DO proc=1,nproc
     ngpwpp(proc)=(dfftp%nwl(proc)+1)/2
     ntot=ntot+ngpwpp(proc)
  END DO

  DO proc=1,nproc

     recvcount(proc) = ngpwpp(proc)
     !
     !!
     ! recvcount(proc) = size of data received from processor proc
     !                   (number of columns times length of each column)
     !
     IF (proc.EQ.1) THEN
        displs(proc)=0
     ELSE
        displs(proc)=displs(proc-1) + recvcount(proc-1)
     END IF
     !
     ! displs(proc) is the position of data received from processor proc
     !
     !         ntot = ntot + recvcount(proc)
     !
     ! ntot = total size of gathered data
     !
  END DO
  !
  ! allocate the needed work spaces
  !
  ALLOCATE( psis( nnrsx ) ) 
  psis( 1:nnrsx ) = 0.D0
  IF (me.EQ.1) THEN
     ALLOCATE(psitot(ntot))
     ALLOCATE(psiwr(nmin(3):nmax(3),nmin(1):nmax(1),nmin(2):nmax(2)))
     !         write(unit) nbsp, nmin, nmax
  END IF
  !
  ! ... fill array psis with c_i(G) (as packed columns along z)
  !
  DO ig=1,ngw
     !
     ! ... ncol+1 is the index of the column
     !
     ncol=(nps(ig)-1)/nr3sx
     !
     ! ... nz_ is the z component in FFT style (refolded between 1 and nr3s)
     !
     nz_ =nps(ig)-ncol*nr3sx
     !
     ! ... nz is the z component in "natural" style (between nmin(3) and nmax(3))
     !
     nz = nz_-1
     IF (nz.GE.nr3s/2) nz=nz-nr3s

     ! ... ncpw(me) columns along z are stored in contiguous order on each node
     !
     psis(ig)=c(ig,jw)
     !
  END DO
  !
  ! ... gather all psis arrays on the first node, in psitot
  !
  CALL mp_barrier()
  !
  CALL MPI_GATHERV (psis, recvcount(me),     MPI_DOUBLE_COMPLEX, &
       &                      psitot,recvcount, displs,MPI_DOUBLE_COMPLEX, &
       &                 root, MPI_COMM_WORLD, ierr)
  IF (ierr.NE.0) CALL errore('write_wfc','MPI_GATHERV',ierr)
  !
  ! write the node-number-independent array
  !
  IF(me.EQ.1) THEN
     DO i=1,ntot
        WRITE(22,*) psitot(i)
     END DO
     WRITE( stdout, * ) "State Written", jw
  END IF
  !
  IF (me.EQ.1) THEN
     DEALLOCATE(psiwr)
     DEALLOCATE(psitot)
  END IF

  DEALLOCATE( psis ) 
  !
  RETURN
  !
END SUBROUTINE write_psi
!
#endif
!
!----------------------------------------------------------------------------
SUBROUTINE rhoiofr( nfi, c, irb, eigrb, bec, &
                    rhovan, rhor, rhog, rhos, enl, ekin, ndwwf )
  !----------------------------------------------------------------------------
  !
  ! ... the normalized electron density rhor in real space
  ! ... the kinetic energy ekin
  ! ... subroutine uses complex fft so it computes two ft's
  ! ... simultaneously
  !
  ! ... rho_i,ij = sum_n < beta_i,i | psi_n >< psi_n | beta_i,j >
  !     
  ! ... < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
  !                            2 sum_g re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
  !
  !     e_v = sum_i,ij rho_i,ij d^ion_is,ji
  !
  USE constants,              ONLY : bohr_radius_angs
  USE kinds,                  ONLY : DP
  USE control_flags,          ONLY : iprint, tbuff, iprsta, thdyn, tpre, trhor
  USE ions_base,              ONLY : nax, nat, nsp, na
  USE cell_base,              ONLY : a1, a2, a3
  USE recvecs_indexes,        ONLY : np, nM
  USE gvecs,                  ONLY : nms, nps, ngs
  USE gvecp,                  ONLY : ngm
  USE gvecb,                  ONLY : ngb
  USE gvecw,                  ONLY : ngw
  USE reciprocal_vectors,     ONLY : gstart
  USE cvan,                   ONLY : ish
  USE grid_dimensions,        ONLY : nr1, nr2, nr3, nr1x, nr2x, nr3x, nnrx
  USE cell_base,              ONLY : omega
  USE smooth_grid_dimensions, ONLY : nr1s, nr2s, nr3s, &
                                     nr1sx, nr2sx, nr3sx, nnrsx
  USE electrons_base,         ONLY : nbspx, nbsp, nspin, f, ispin
  USE constants,              ONLY : pi, fpi
  USE wannier_base,           ONLY : iwf
  USE dener,                  ONLY : dekin, denl
  USE io_global,              ONLY : stdout, ionode
  USE uspp_param,             ONLY : nh, nhm
  USE uspp,                   ONLY : nkb
  USE fft_module,             ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  REAL(DP) bec(nkb,nbsp), rhovan(nhm*(nhm+1)/2,nat,nspin)
  REAL(DP) rhovanaux(nhm,nhm,nat,nspin)
  REAL(DP) rhor(nnrx,nspin), rhos(nnrsx,nspin)
  REAL(DP) enl, ekin
  COMPLEX(DP) eigrb(ngb,nat), c(ngw,nbspx), rhog(ngm,nspin)
  INTEGER irb(3,nat), nfi, ndwwf
  ! local variables
  INTEGER iss, isup, isdw, iss1, iss2, ios, i, ir, ig
  INTEGER is,iv,jv,isa,isn, jnl, j, k, inl, ism, ia
  REAL(DP) rsumr(2), rsumg(2), sa1, sa2, sums(2)
  REAL(DP) rnegsum, rmin, rmax, rsum
  REAL(DP) enkin, ennl
  COMPLEX(DP) fp,fm
  COMPLEX(DP), ALLOCATABLE :: psi(:), psis(:)
  EXTERNAL ennl, enkin
  !
  !
  CALL start_clock(' rhoiofr ')
  ALLOCATE( psi( nnrx ) )
  ALLOCATE( psis( nnrsx ) ) 
  !
  DO iss=1,nspin
     rhor( 1:nnrx, iss ) = 0.D0
     rhos( 1:nnrsx, iss ) = 0.D0
     rhog( 1:ngm, iss ) = 0.D0
  END DO
  !
  !     ==================================================================
  !     calculation of kinetic energy ekin
  !     ==================================================================
  ekin=enkin(c,ngw,f,nbsp)
  IF(tpre) CALL denkin(c,dekin)
  !
  !     ==================================================================
  !     calculation of non-local energy
  !     ==================================================================
  !      enl=ennl(rhovan,bec)
  DO is=1,nsp
     DO iv=1, nh(is)
        DO jv=iv,nh(is)
           isa=0
           DO ism=1,is-1
              isa=isa+na(ism)
           END DO
           DO ia=1,na(is)
              inl=ish(is)+(iv-1)*na(is)+ia
              jnl=ish(is)+(jv-1)*na(is)+ia
              isa=isa+1
              sums(1)=f(iwf)*bec(inl,iwf)*bec(jnl,iwf)
              rhovanaux(iv,jv,isa,1) = sums(1)
              rhovanaux(jv,iv,isa,1) = sums(1)
           END DO
        END DO
     END DO
  END DO

  k=1
  DO i=1,nhm
     DO j=i,nhm
        rhovan(k,:,:)=rhovanaux(j,i,:,:)
        k=k+1
     END DO
  END DO

  IF(tpre) CALL dennl(bec,denl)
  !    
  !    warning! trhor and thdyn are not compatible yet!   
  !
  IF(trhor.AND.(.NOT.thdyn))THEN
     ! 
     CALL read_rho( nspin, rhor )
     !
     IF(nspin.EQ.1)THEN
        iss=1
        DO ir=1,nnrx
           psi(ir)=CMPLX(rhor(ir,iss),0.d0)
        END DO
        CALL fwfft('Dense',psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
        DO ig=1,ngm
           rhog(ig,iss)=psi(np(ig))
        END DO
     ELSE
        isup=1
        isdw=2
        DO ir=1,nnrx
           psi(ir)=CMPLX(rhor(ir,isup),rhor(ir,isdw))
        END DO
        CALL fwfft('Dense',psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
        DO ig=1,ngm
           fp=psi(np(ig))+psi(nm(ig))
           fm=psi(np(ig))-psi(nm(ig))
           rhog(ig,isup)=0.5*CMPLX( DBLE(fp),AIMAG(fm))
           rhog(ig,isdw)=0.5*CMPLX(AIMAG(fp),-DBLE(fm))
        END DO
     ENDIF
     !
  ELSE
     !     ==================================================================
     !     self-consistent charge
     !     ==================================================================
     !
     !     important: if nbsp is odd then nbspx must be .ge.nbsp+1 and c(*,nbsp+1)=0.
     ! 
     !         if (mod(nbsp,2).ne.0) then
     !            do ig=1,ngw
     !               c(ig,nbsp+1)=(0.,0.)
     !            end do
     !         endif
     !
     !         do i=1,nbsp,2
     i=iwf
     psis( 1:nnrsx ) = 0.D0
     DO ig=1,ngw
        !               c(ig,i+1)=(0.,0.)
        psis(nms(ig))=CONJG(c(ig,i))
        psis(nps(ig))=c(ig,i)
     END DO
     !
     CALL invfft('Wave',psis,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
     !
     !     wavefunctions in unit 21
     !
     IF(tbuff) WRITE(21,iostat=ios) psis
     !            iss1=ispin(i)
     iss1=1
     sa1=f(i)/omega
     !            if (i.ne.nbsp) then
     !              iss2=ispin(i+1)
     !              sa2=f(i+1)/omega
     !            else
     iss2=iss1  ! carlo
     sa2=0.0
     !            end if
     DO ir=1,nnrsx
        rhos(ir,iss1)=rhos(ir,iss1) + sa1*( DBLE(psis(ir)))**2
        rhos(ir,iss2)=rhos(ir,iss2) + sa2*(AIMAG(psis(ir)))**2
     END DO
     !
     !       buffer 21
     !     
     IF(tbuff) THEN
        IF(ios.NE.0) CALL errore                                  &
             &              (' rhoofr',' error in writing unit 21',ios)
     ENDIF
     !
     !         end do
     !
     IF(tbuff) REWIND 21
     !
     !     smooth charge in g-space is put into rhog(ig)
     !
     IF(nspin.EQ.1)THEN
        iss=1
        DO ir=1,nnrsx
           psis(ir)=CMPLX(rhos(ir,iss),0.d0)
        END DO
        CALL fwfft('Smooth',psis,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
        DO ig=1,ngs
           rhog(ig,iss)=psis(nps(ig))
        END DO
     ELSE
        isup=1
        isdw=2
        DO ir=1,nnrsx
           psis(ir)=CMPLX(rhos(ir,isup),rhos(ir,isdw))
        END DO
        CALL fwfft('Smooth',psis,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
        DO ig=1,ngs
           fp= psis(nps(ig)) + psis(nms(ig))
           fm= psis(nps(ig)) - psis(nms(ig))
           rhog(ig,isup)=0.5*CMPLX( DBLE(fp),AIMAG(fm))
           rhog(ig,isdw)=0.5*CMPLX(AIMAG(fp),-DBLE(fm))
        END DO
     ENDIF
     !
     IF(nspin.EQ.1) THEN
        !     ==================================================================
        !     case nspin=1
        !     ------------------------------------------------------------------
        iss=1
        psi( 1:nnrx ) = 0.D0
        DO ig=1,ngs
           psi(nm(ig))=CONJG(rhog(ig,iss))
           psi(np(ig))=      rhog(ig,iss)
        END DO
        CALL invfft('Dense',psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
        DO ir=1,nnrx
           rhor(ir,iss)=DBLE(psi(ir))
        END DO
     ELSE 
        !     ==================================================================
        !     case nspin=2
        !     ------------------------------------------------------------------
        isup=1
        isdw=2
        psi( 1:nnrx ) = 0.D0
        DO ig=1,ngs
           psi(nm(ig))=CONJG(rhog(ig,isup))+CI*CONJG(rhog(ig,isdw))
           psi(np(ig))=rhog(ig,isup)+CI*rhog(ig,isdw)
        END DO
        CALL invfft('Dense',psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
        DO ir=1,nnrx
           rhor(ir,isup)= DBLE(psi(ir))
           rhor(ir,isdw)=AIMAG(psi(ir))
        END DO
     ENDIF
     !
     !         if(iprsta.ge.3)then
     WRITE( stdout,*) 'Smooth part of charge density :'
     DO iss=1,nspin
        rsumg(iss)=omega*DBLE(rhog(1,iss))
        rsumr(iss)=SUM(rhor(1:nnrx,iss))*omega/DBLE(nr1*nr2*nr3)
     END DO
#ifdef __PARA
     IF (gstart.NE.2) THEN
        ! in the parallel case, only one processor has G=0 ! 
        DO iss=1,nspin
           rsumg(iss)=0.0
        END DO
     END IF
     CALL reduce(nspin,rsumg)
     CALL reduce(nspin,rsumr)
#endif
     IF (nspin.EQ.1) THEN
        WRITE( stdout,1) rsumg(1),rsumr(1)
     ELSE
        WRITE( stdout,2) (rsumg(iss),iss=1,nspin),(rsumr(iss),iss=1,nspin)
     ENDIF
     !         endif
     !     ==================================================================
     !
     !     add vanderbilt contribution to the charge density
     !
     !     drhov called before rhov because input rho must be the smooth part
     !
     IF (tpre) CALL drhov(irb,eigrb,rhovan,rhog,rhor)
     !
     CALL rhov(irb,eigrb,rhovan,rhog,rhor)
  ENDIF
  !     ======================================endif for trhor=============
  REWIND ndwwf
  !
  IF ( ionode ) THEN
     !
     WRITE( ndwwf, '("3  2")' )
     !
     WRITE( ndwwf, '(3(2X,I3))' ) nr1x, nr2x, nr3x
     !
     WRITE( ndwwf, '(3(2X,"0",2X,F16.10))' ) &
         ( DBLE( nr1x - 1 ) / DBLE( nr1x ) ) * a1(1) * bohr_radius_angs, &
         ( DBLE( nr2x - 1 ) / DBLE( nr2x ) ) * a2(2) * bohr_radius_angs, &
         ( DBLE( nr3x - 1 ) / DBLE( nr3x ) ) * a3(3) * bohr_radius_angs
     !
  END IF
  !
  !     here to check the integral of the charge density
  !
  !
  IF(iprsta.GE.2) THEN
     CALL checkrho(nnrx,nspin,rhor,rmin,rmax,rsum,rnegsum)
     rnegsum=rnegsum*omega/DBLE(nr1*nr2*nr3)
     rsum=rsum*omega/DBLE(nr1*nr2*nr3)
     WRITE( stdout,'(a,4(1x,f12.6))')                                     &
          &     ' rhoofr: rmin rmax rnegsum rsum  ',rmin,rmax,rnegsum,rsum
  END IF
  !
  IF(nfi.EQ.0.OR.MOD(nfi-1,iprint).EQ.0) THEN

     WRITE( stdout, * ) 
     WRITE( stdout, * ) 'Smooth part + Augmentatio Part: '
     DO iss=1,nspin
        rsumg(iss)=omega*DBLE(rhog(1,iss))
        rsumr(iss)=SUM(rhor(1:nnrx,iss))*omega/DBLE(nr1*nr2*nr3)
     END DO
#ifdef __PARA
     IF (gstart.NE.2) THEN
        ! in the parallel case, only one processor has G=0 ! 
        DO iss=1,nspin
           rsumg(iss)=0.0
        END DO
     END IF
     CALL reduce(nspin,rsumg)
     CALL reduce(nspin,rsumr)
#endif
     IF (nspin.EQ.1) THEN
        WRITE( stdout,1) rsumg(1),rsumr(1)
     ELSE
        IF(iprsta.GE.3)                                             &
             &          WRITE( stdout,2)  rsumg(1),rsumg(2),rsumr(1),rsumr(2)
        WRITE( stdout,1) rsumg(1)+rsumg(2),rsumr(1)+rsumr(2)
     ENDIF
  ENDIF
  !
2 FORMAT(//' subroutine rhoofr: total integrated electronic',       &
       &     ' density'/' in g-space =',f10.6,2x,f10.6,4x,                &
       &     ' in r-space =',f10.6,2x,f10.6)
1 FORMAT(//' subroutine rhoofr: total integrated electronic',       &
       &     ' density'/' in g-space =',f10.6,4x,                         &
       &     ' in r-space =',f10.6)
  !

  WRITE( stdout, * ) 
  WRITE( stdout , * ) 'State Written : ' ,iwf, 'to unit',ndwwf
  WRITE( stdout, * ) 

  DEALLOCATE( psi  )
  DEALLOCATE( psis ) 

  CALL stop_clock(' rhoiofr ')
  !
  RETURN
END SUBROUTINE rhoiofr
