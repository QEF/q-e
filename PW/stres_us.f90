!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------

subroutine stres_us (ik, gk, sigmanlc)
  !----------------------------------------------------------------------
  ! nonlocal (separable pseudopotential) contribution to the stress
  !
#include "machine.h"
  use pwcom
  USE wavefunctions,    ONLY : evc
  use becmod
#ifdef __PARA
  use para
#endif
  implicit none
  integer :: ik

  real(kind=DP) :: sigmanlc (3, 3), gk (3, npw)
  real(kind=DP) :: eps

  parameter (eps = 1.0e-8)
  integer :: na, np, ibnd, ipol, jpol, l, i, ikb, jkb, ih, jh, &
       ijkb0
  real(kind=DP) :: fac, xyz (3, 3), q, evps, DDOT
  real(kind=DP) , allocatable :: qm1(:)
  complex(kind=DP) , allocatable :: work1 (:), work2 (:), dvkb (:,:)
  ! dvkb contains the derivatives of the kb potential
  complex(kind=DP) :: ps
  ! xyz are the three unit vectors in the x,y,z directions
  data xyz / 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /
  !
  if (nkb.eq.0) return
  !
  if (lsda) current_spin = isk (ik)
  if (nks.gt.1) call init_us_2 (npw, igk, xk (1, ik), vkb)
  call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
  allocate (work1(npwx), work2(npwx), qm1(npwx))

  do i = 1, npw
     q = sqrt (gk (1, i) **2 + gk (2, i) **2 + gk (3, i) **2)
     if (q.gt.eps) then
        qm1 (i) = 1.d0 / q
     else
        qm1 (i) = 0.d0
     endif
  enddo
#ifdef __PARA
  if (me.ne.1) goto 100
#endif
  !
  !     diagonal contribution
  !
  evps = 0.d0
  do ibnd = 1, nbnd
     fac = wg (ibnd, ik)
     ijkb0 = 0
     do np = 1, ntyp
        do na = 1, nat
           if (ityp (na) .eq.np) then
              do ih = 1, nh (np)
                 ikb = ijkb0 + ih
                 ps = deeq (ih, ih, na, current_spin) - &
                      et (ibnd, ik) * qq (ih, ih, np)
                 evps = evps + fac * ps * abs (becp (ikb, ibnd) ) **2

                 if (tvanp (np) .or.newpseudo (np) ) then
                    !    only in the US case there is a contribution for jh<>ih. We have use
                    !    here the symmetry in the interchange of ih and jh
                    !
                    do jh = ih + 1, nh (np)
                       jkb = ijkb0 + jh
                       ps = deeq (ih, jh, na, current_spin) - et (ibnd, ik) &
                            * qq (ih, jh, np)
                       evps = evps + ps * fac * 2.0d0 * &
                            DREAL(conjg (becp (ikb,ibnd) ) * becp (jkb, ibnd) )
                    enddo
                 endif
              enddo
              ijkb0 = ijkb0 + nh (np)
           endif
        enddo
     enddo

  enddo
  do l = 1, 3
     sigmanlc (l, l) = sigmanlc (l, l) - evps

  enddo
  !      write (6,*) ' non local energy ', evps, evps*uakbar/omega
100 continue
  !
  !     non diagonal contribution - derivative of the bessel function
  !

  allocate ( dvkb(npwx,nkb) )

  call gen_us_dj (ik, dvkb)
  do ibnd = 1, nbnd
     work2 (:) = (0.d0,0.d0)
     ijkb0 = 0
     do np = 1, ntyp
        do na = 1, nat
           if (ityp (na) .eq.np) then
              do ih = 1, nh (np)
                 ikb = ijkb0 + ih
                 if (.not. (tvanp (np) .or.newpseudo (np) ) ) then
                    ps = becp (ikb, ibnd) * (deeq (ih, ih, na, current_spin) &
                         - et (ibnd, ik) * qq (ih, ih, np) )
                 else
                    !
                    !    in the US case there is a contribution also for jh<>ih
                    !
                    ps = (0.d0, 0.d0)
                    do jh = 1, nh (np)
                       jkb = ijkb0 + jh
                       ps = ps + becp (jkb, ibnd) * (deeq (ih, jh, na, &
                            current_spin) - et (ibnd, ik) * qq (ih, jh, np) )
                    enddo
                 endif
                 call ZAXPY (npw, ps, dvkb (1, ikb), 1, work2, 1)
              enddo
              ijkb0 = ijkb0 + nh (np)
           endif
        enddo

     enddo
     do ipol = 1, 3
        do jpol = 1, ipol
           do i = 1, npw
              work1 (i) = evc(i, ibnd) * gk(ipol, i) * gk(jpol, i) * qm1(i)
           enddo
           sigmanlc (ipol, jpol) = sigmanlc (ipol, jpol) - &
                2.d0 * wg (ibnd, ik) * DDOT (2 * npw, work1, 1, work2, 1)
        enddo
     enddo
  enddo
  !
  !     non diagonal contribution - derivative of the spherical harmonics
  !     (no contribution from l=0)
  !

  if (lmaxkb.eq.0) goto 10

  do ipol = 1, 3
     call gen_us_dy (ik, xyz (1, ipol), dvkb)
     do ibnd = 1, nbnd
        work2(:) = (0.d0,0.d0)
        ijkb0 = 0
        do np = 1, ntyp
           do na = 1, nat
              if (ityp (na) .eq.np) then
                 do ih = 1, nh (np)
                    ikb = ijkb0 + ih
                    if (.not. (tvanp (np) .or.newpseudo (np) ) ) then
                       ps = becp(ikb, ibnd) * (deeq(ih,ih,na,current_spin) &
                            - et (ibnd, ik) * qq (ih, ih, np) )
                    else
                       !
                       ! in the US case there is a contribution also for jh<>ih
                       !
                       ps = (0.d0, 0.d0)
                       do jh = 1, nh (np)
                          jkb = ijkb0 + jh
                          ps = ps + becp (jkb, ibnd) * (deeq (ih, jh, na, &
                               current_spin) - et(ibnd, ik) * qq(ih, jh, np) )
                       enddo
                    endif
                    call ZAXPY (npw, ps, dvkb (1, ikb), 1, work2, 1)
                 enddo
                 ijkb0 = ijkb0 + nh (np)
              endif
           enddo
        enddo

        do jpol = 1, ipol
           do i = 1, npw
              work1 (i) = evc (i, ibnd) * gk (jpol, i)
           enddo
           sigmanlc (ipol, jpol) = sigmanlc (ipol, jpol) - &
                2.d0 * wg (ibnd,ik) * DDOT(2 * npw, work1, 1, work2, 1)
        enddo
     enddo
  enddo

10 continue
  deallocate (dvkb)
  deallocate (work1, work2, qm1)

  return
end subroutine stres_us

