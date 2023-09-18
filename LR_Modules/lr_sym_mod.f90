! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module contains routines that rotate the charge density and the
! magnetization induced by a phonon or an electric field perturbation.
!
MODULE lr_sym_mod

IMPLICIT NONE
SAVE
PRIVATE

PUBLIC rotate_mesh, rotate_mesh_1s, find_mesh_ijk, compute_phase, &
       psymeq

CONTAINS
!
!---------------------------------------------------------------------
SUBROUTINE rotate_mesh(my_nrxx, nsym, rir)
!---------------------------------------------------------------------

USE fft_base,  ONLY : dfftp
USE symm_base, ONLY : s, ft

IMPLICIT NONE
INTEGER, INTENT(IN) :: my_nrxx, nsym
INTEGER, INTENT(INOUT) :: rir(my_nrxx, nsym)
INTEGER :: isym, ftau(3,48)

ftau(1,1:nsym) = NINT ( ft(1,1:nsym)*dfftp%nr1 ) 
ftau(2,1:nsym) = NINT ( ft(2,1:nsym)*dfftp%nr2 ) 
ftau(3,1:nsym) = NINT ( ft(3,1:nsym)*dfftp%nr3 )

DO isym=1, nsym
   CALL rotate_mesh_1s(my_nrxx, s(1,1,isym), ftau(1,isym), rir(1,isym))
ENDDO

RETURN
END SUBROUTINE rotate_mesh
!
!---------------------------------------------------------------------
SUBROUTINE rotate_mesh_1s(my_nrxx, s, ftau, rir)
!---------------------------------------------------------------------

USE fft_base,  ONLY : dfftp

IMPLICIT NONE
INTEGER, INTENT(IN) :: my_nrxx, s(3,3), ftau(3)
INTEGER, INTENT(INOUT) :: rir(my_nrxx)

INTEGER :: j0, k0, ir, idx, i, j, k, ri, rj, rk, nr1x, nr2x, nr12x, ss(3,3)
INTEGER :: nr1, nr2, nr3

nr1 = dfftp%nr1
nr2 = dfftp%nr2
nr3 = dfftp%nr3
nr1x = dfftp%nr1x
nr2x = dfftp%nr2x
nr12x = nr1x * nr2x

rir=0
ss (1, 1) = s (1, 1)
ss (2, 1) = s (2, 1) * nr1 / nr2
ss (3, 1) = s (3, 1) * nr1 / nr3
ss (1, 2) = s (1, 2) * nr2 / nr1
ss (2, 2) = s (2, 2)
ss (3, 2) = s (3, 2) * nr2 / nr3
ss (1, 3) = s (1, 3) * nr3 / nr1
ss (2, 3) = s (2, 3) * nr3 / nr2
ss (3, 3) = s (3, 3)
   
j0 = dfftp%my_i0r2p
k0 = dfftp%my_i0r3p
DO ir = 1, my_nrxx
   idx = ir - 1
   k   = idx / (nr1x*dfftp%my_nr2p)
   idx = idx - (nr1x*dfftp%my_nr2p)*k
   k   = k + k0
   j   = idx / nr1x
   idx = idx - nr1x * j
   j   = j + j0
   i   = idx
   i=i + 1
   j=j + 1
   k=k + 1
   IF (i > nr1 .OR. j > nr2 .OR. k > nr3) CYCLE
   CALL ruotaijk (ss, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk ) 
   rir(ir)=ri+(rj-1)*nr1x+(rk-1)*nr12x
ENDDO

RETURN
END SUBROUTINE rotate_mesh_1s

SUBROUTINE find_mesh_ijk(my_nrxx, iir, jir, kir)

USE fft_base, ONLY : dfftp

IMPLICIT NONE
INTEGER, INTENT(IN) :: my_nrxx
INTEGER, INTENT(INOUT) :: iir(my_nrxx), jir(my_nrxx), kir(my_nrxx)

INTEGER :: nr1x, i, j, k, idx, j0, k0, ir

nr1x = dfftp%nr1x

j0 = dfftp%my_i0r2p
k0 = dfftp%my_i0r3p

DO ir = 1, my_nrxx
   idx = ir - 1
   k   = idx / (nr1x*dfftp%my_nr2p)
   idx = idx - (nr1x*dfftp%my_nr2p)*k
   k   = k + k0
   j   = idx / nr1x
   idx = idx - nr1x * j
   j   = j + j0
   i   = idx
   i   = i + 1
   j   = j + 1
   k   = k + 1
   iir(ir) = i
   jir(ir) = j
   kir(ir) = k
ENDDO

RETURN
END SUBROUTINE find_mesh_ijk
!
!---------------------------------------------------------------------
SUBROUTINE compute_phase(phase1, phase2, phase3, nr1, nr2, nr3, nsym, gi, zb)
!---------------------------------------------------------------------
!
!  This routine computes the phases e^{i G r} for the vector gi, which
!  is supposed to be in cartesian coordinates. 
!
USE kinds,     ONLY : DP
USE constants, ONLY : tpi
USE cell_base, ONLY : at

IMPLICIT NONE
INTEGER, INTENT(IN) :: nr1, nr2, nr3, nsym
COMPLEX(DP), INTENT(INOUT) :: phase1(nr1, nsym), phase2(nr2, nsym), &
                                                 phase3(nr3, nsym)
REAL(DP), INTENT(IN) :: gi(3, nsym)
LOGICAL, INTENT(OUT) :: zb(nsym)

INTEGER  :: isym, i, j, k
INTEGER  :: gk(3, nsym)
REAL(DP) :: arg, wrk(3,nsym)
COMPLEX(DP) :: term
!
!  Bring gi to crystal coordinates (of the bg). 
!
wrk(:,:)=gi(:,:)
CALL cryst_to_cart(nsym, wrk, at, -1)
gk(:,:)=NINT(wrk(:,:))

DO isym=1, nsym
   zb(isym) = (gk(1,isym) /= 0 .OR. gk(2,isym) /= 0 .OR. gk(3,isym) /= 0 )

   IF (gk(1, isym) /= 0) THEN
      phase1(1,isym)=(1.0_DP,0.0_DP)
      arg = tpi*gk(1, isym)/DBLE(nr1)
      term= CMPLX(COS(arg), SIN(arg), KIND=DP)
      DO i = 2, nr1
         phase1(i,isym) = phase1(i-1,isym) * term
      ENDDO
   ELSE
      phase1(:,isym)=(1.0_DP,0.0_DP)
   ENDIF

   IF (gk(2, isym) /= 0) THEN
      phase2(1, isym)=(1.0_DP,0.0_DP)
      arg = tpi*gk(2, isym)/DBLE(nr2)
      term= CMPLX(COS(arg), SIN(arg), KIND=DP)
      DO j = 2, nr2
           phase2(j, isym) = phase2(j-1, isym) * term
      ENDDO
   ELSE
      phase2(:, isym)=(1.0_DP,0.0_DP)
   ENDIF

   IF (gk(3,isym) /= 0) THEN
      phase3(1,isym)=(1.0_DP,0.0_DP)
      arg = tpi*gk(3,isym)/DBLE(nr3)
      term= CMPLX(COS(arg), SIN(arg), KIND=DP)
      DO k = 2, nr3
         phase3(k, isym) = phase3(k-1, isym) * term
      ENDDO
   ELSE
      phase3(:, isym)=(1.0_DP,0.0_DP)
   ENDIF
ENDDO

RETURN
END SUBROUTINE compute_phase
!
!---------------------------------------------------------------------
SUBROUTINE psymeq (dvsym)
!---------------------------------------------------------------------
!
!     This routine symmetrizes the change of the charge/magnetization
!     or potential/magnetic field due to an electric field perturbation 
!     written as E e^iqr. The three components of the magnetization or
!     magnetic field are in cartesian coodinates.
!
USE kinds,     ONLY : DP
USE cell_base, ONLY : at, bg
USE fft_base,  ONLY : dfftp
USE symm_base, ONLY : s, sname, invs, t_rev
USE noncollin_module, ONLY : noncolin, domag, nspin_lsda, nspin_mag
USE lr_symm_base, ONLY : nsymq, gi
USE scatter_mod, ONLY : cgather_sym
IMPLICIT NONE

COMPLEX(DP) :: dvsym (dfftp%nnr, nspin_mag)
  ! the potential to symmetrize
COMPLEX(DP), ALLOCATABLE ::  aux(:), aux_nc(:,:)
  ! auxiliary quantity
COMPLEX(DP), ALLOCATABLE ::  phase1(:,:), phase2(:,:), phase3(:,:)
  ! additional phases
LOGICAL, ALLOCATABLE :: zb(:)
  ! if .true. the symmetry requires a G vector for giving Sq=q+G 
INTEGER :: is, isym, ipol, jpol, kpol, my_nrxx, nr1x, ir, i, j, k, &
           nr1, nr2, nr3
  ! counter on spin polarization
  ! counter on rotations
  ! counter on symmetries
  ! counter on polarizations
  ! number of points on this processor
  ! nrx1 size of the FFT mesh array along direction 1
  ! counters on mesh points
  ! sizes of the FFT mesh
INTEGER, ALLOCATABLE :: iir(:), jir(:), kir(:), rir(:,:)

COMPLEX(DP) :: dmags(3), mag(3), phase

INTEGER :: stilde(3,3,nsymq)

IF (nsymq == 1) RETURN

nr1=dfftp%nr1
nr2=dfftp%nr2
nr3=dfftp%nr3
nr1x=dfftp%nr1x
my_nrxx=nr1x*dfftp%my_nr2p*dfftp%my_nr3p

ALLOCATE(rir(my_nrxx, nsymq))
ALLOCATE(iir(my_nrxx))
ALLOCATE(jir(my_nrxx))
ALLOCATE(kir(my_nrxx))
ALLOCATE(phase1(nr1, nsymq))
ALLOCATE(phase2(nr2, nsymq))
ALLOCATE(phase3(nr3, nsymq))
ALLOCATE(zb(nsymq))
ALLOCATE (aux(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))

CALL rotate_mesh(my_nrxx, nsymq, rir)
CALL find_mesh_ijk(my_nrxx, iir, jir, kir)
CALL compute_phase(phase1, phase2, phase3, nr1, nr2, nr3, nsymq, gi, zb)

DO is = 1, nspin_lsda
   !
   !  collect all the quantity to symmetrize in all processors
   !
   CALL cgather_sym(dfftp, dvsym(:, is), aux(:))
   dvsym(:,is)=(0.0_DP, 0.0_DP)
   !
   !  symmmetrize. Each processor symmetrizes only the points that it has
   !
   DO isym = 1, nsymq
      IF (zb(isym)) THEN
         DO ir = 1, my_nrxx
            IF (rir(ir,isym)==0) CYCLE
            i=iir(ir)
            j=jir(ir)
            k=kir(ir)
            phase=phase1(i,isym)*phase2(j,isym)*phase3(k,isym)
            dvsym(ir,is) = dvsym(ir,is) + aux(rir(ir,isym))*phase
         ENDDO
      ELSE
         DO ir = 1, my_nrxx
            IF (rir(ir,isym)==0) CYCLE
            dvsym(ir,is) = dvsym(ir,is) + aux(rir(ir,isym))
         ENDDO
      ENDIF
   ENDDO
ENDDO

DEALLOCATE (aux)
!
!  Rotate the magnetization in the noncollinear magnetic case
!
IF (noncolin.AND.domag) THEN
   ALLOCATE (aux_nc(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, 3))
!
!  set the symmetry matrices that rotate the magnetization
!
   CALL set_stilde(s, stilde, sname, t_rev, invs, nsymq)  
!
!  bring the magnetization in crystal coordinates
!  collect all the quantity to symmetrize in all processors
!
   CALL ccryst_to_cart_t (dfftp%nnr, dvsym(:,2:4), bg, -1)
   DO is = 2, nspin_mag
      CALL cgather_sym( dfftp, dvsym(:, is), aux_nc(:, is-1))
      dvsym(:,is) = (0.0_DP, 0.0_DP)
   ENDDO
   !
   DO isym = 1, nsymq
      DO ir = 1, my_nrxx
         IF (rir(ir,isym)==0) CYCLE
         IF (zb(isym)) THEN
            i=iir(ir)
            j=jir(ir)
            k=kir(ir)
            phase=phase1(i,isym)*phase2(j,isym)*phase3(k,isym)
            dmags(:)=aux_nc(rir(ir,isym),:) * phase
         ELSE
            dmags(:)=aux_nc(rir(ir,isym),:) * phase
         ENDIF      

!
! rotate the magnetic moment
!
         DO kpol = 1, 3
            mag(kpol) = stilde(1,kpol,isym)*dmags(1) + &
                        stilde(2,kpol,isym)*dmags(2) + &
                        stilde(3,kpol,isym)*dmags(3)
         ENDDO
!
! and add to the complete field
!
         IF (t_rev(isym)==1) THEN
            DO kpol = 1, 3
               dvsym(ir,kpol+1) = dvsym(ir,kpol+1) + CONJG(mag(kpol))
            ENDDO
         ELSE
            DO kpol = 1, 3
               dvsym(ir,kpol+1) = dvsym(ir,kpol+1) + mag(kpol)
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   DEALLOCATE (aux_nc)
!
! go back to cartesian coordinates
!
   CALL ccryst_to_cart_t (dfftp%nnr, dvsym(:,2:4), at, 1)
ENDIF

dvsym = dvsym / DBLE(nsymq)

DEALLOCATE(iir)
DEALLOCATE(jir)
dEALLOCATE(kir)
DEALLOCATE(rir)
DEALLOCATE(phase1)
DEALLOCATE(phase2)
DEALLOCATE(phase3)
DEALLOCATE(zb)

RETURN
END SUBROUTINE psymeq

SUBROUTINE set_stilde(s, stilde, sname, t_rev, invs, nsym)
!
!  This routine sets the matrices needed to rotate the magnetization.
!  These matrices contain the proper part of the inverse of s 
!  (when t_rev==0) or minus the proper part of the inverse of s 
!  (when t_rev==1).
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nsym
INTEGER, INTENT(IN) :: s(3,3,nsym), t_rev(nsym), invs(nsym)
INTEGER, INTENT(INOUT) :: stilde(3,3,nsym)
CHARACTER(LEN=45), INTENT(IN) :: sname(nsym)

INTEGER :: isym

DO isym = 1, nsym
   stilde(:,:,isym)=s(:,:,invs(isym))
   IF (sname(isym)(1:3)=='inv') stilde(:,:,isym)=-stilde(:,:,isym)
   IF (t_rev(isym)==1) stilde(:,:,isym)=-stilde(:,:,isym)
ENDDO

RETURN
END SUBROUTINE set_stilde

!-----------------------------------------------------------------------
subroutine ccryst_to_cart_t (nvec, vec, trmat, iflag)
  !-----------------------------------------------------------------------
  !
  !     This routine transforms the atomic positions or the k-point
  !     components from crystallographic to cartesian coordinates
  !     ( iflag=1 ) and viceversa ( iflag=-1 ).
  !     Output cartesian coordinates are stored in the input ('vec') array
  !
  !
  USE kinds, ONLY : DP
  implicit none
  !
  integer, intent(in) :: nvec, iflag
  ! nvec:  number of vectors (atomic positions or k-points)
  !        to be transformed from crystal to cartesian and vice versa
  ! iflag: gives the direction of the transformation
  real(DP), intent(in) :: trmat (3, 3)
  ! trmat: transformation matrix
  ! if iflag=1:
  !    trmat = at ,  basis of the real-space lattice,       for atoms   or
  !          = bg ,  basis of the reciprocal-space lattice, for k-points
  ! if iflag=-1: the opposite
  COMPLEX(DP), intent(inout) :: vec (nvec,3)
  ! coordinates of the vector (atomic positions or k-points) to be
  ! transformed - overwritten on output
  !
  !    local variables
  !
  integer :: nv, kpol
  ! counter on vectors
  ! counter on polarizations
  !
  COMPLEX(DP) :: veu(nvec,3), trmatc(3,3)

  trmatc=CMPLX(trmat,0.0_DP)
  !
  veu=vec
  IF (iflag.EQ.1) THEN
     CALL zgemm('N','T',nvec,3,3,(1.0_DP,0.0_DP),veu,nvec,trmatc,&
                                           3,(0.0_DP,0.0_DP),vec,nvec)
  ELSE
     CALL zgemm('N','N',nvec,3,3,(1.0_DP,0.0_DP),veu,nvec,trmatc,&
                                           3,(0.0_DP,0.0_DP),vec,nvec)
  ENDIF

  return
end subroutine ccryst_to_cart_t

!----------------------------------------------------------------------
subroutine ruotaijk (ss, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk)
  !----------------------------------------------------------------------
  !
  !    This routine computes the rotated of the point i,j,k throught
  !    the symmetry (s,f). Then it computes the equivalent point
  !    on the original mesh
  !
  !
  USE kinds
  implicit none
  !
  !    first the dummy variables
  !
  integer :: ss (3, 3), ftau (3), i, j, k, nr1, nr2, nr3, ri, rj, rk
  ! input: the rotation matrix
  ! input: the fractionary translation
  !   !   input: the point to rotate

  ! /
  !   !   input: the dimension of the mesh

  ! /
  !  !  output: the rotated point

  !/
  !
  ri = ss (1, 1) * (i - 1) + ss (2, 1) * (j - 1) + ss (3, 1) &
       * (k - 1) - ftau (1)
  ri = mod (ri, nr1) + 1
  if (ri.lt.1) ri = ri + nr1
  rj = ss (1, 2) * (i - 1) + ss (2, 2) * (j - 1) + ss (3, 2) &
       * (k - 1) - ftau (2)
  rj = mod (rj, nr2) + 1
  if (rj.lt.1) rj = rj + nr2
  rk = ss (1, 3) * (i - 1) + ss (2, 3) * (j - 1) + ss (3, 3) &
       * (k - 1) - ftau (3)
  rk = mod (rk, nr3) + 1
  if (rk.lt.1) rk = rk + nr3

  return
end subroutine ruotaijk

END MODULE lr_sym_mod
