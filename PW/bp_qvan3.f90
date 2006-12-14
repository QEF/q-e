!
! Copyright (C) 2004 Vanderbilt's group at Rutgers University, NJ
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
      subroutine qvan3(iv,jv,is,qg,ylm_k,qr)
!--------------------------------------------------------------------------
!
!     calculate qg = SUM_LM (-I)^L AP(LM,iv,jv) YR_LM QRAD(iv,jv,L,is)
      USE kinds, ONLY: DP
      USE ions_base,  ONLY : ntyp => nsp
      USE us, ONLY: dq, qrad
      USE uspp_param, ONLY: lmaxq, nbetam
      USE uspp, ONLY: nlx, lpl, lpx, ap, indv, nhtol, nhtolm

      implicit none
      integer :: iv,jv,is
      complex(DP) :: qg,sig
      real(DP) :: ylm_k(lmaxq*lmaxq)
      real(DP) :: qr(nbetam,nbetam,lmaxq,ntyp)
      
      integer ivs,jvs,ivl,jvl,lp,l,i
!       IV  = 1..8    ! s_1 p_x1 p_y1 p_z1 s_2 p_x2 p_z2 p_y2
!       IVS = 1..4    ! s_1 s_2 p_1 p_2 d_1 d_2
!       IVL = 1..4    ! s p_x p_y p_z
!
!  NOTE :   IV  = 1..8 (sppp sppp)   IVS = 1..4 (sspp) OR 1..2 (sp)
!           IVL = 1..4 (sppp)
!
      ivs = indv(iv,is)
      jvs = indv(jv,is)
      ivl = nhtolm(iv,is)
      jvl = nhtolm(jv,is)

      if (ivs > nbetam .OR. jvs > nbetam) &
           call errore (' qvan ', ' wrong dimensions (1)', MAX(ivs,jvs))
      if (ivl > nlx .OR. jvl > nlx) &
       call errore (' qvan ', ' wrong dimensions (2)', MAX(ivl,jvl))
 
      qg = (0.0d0,0.0d0)

!odl                  Write(*,*) 'QVAN3  --  ivs jvs = ',ivs,jvs
!odl                  Write(*,*) 'QVAN3  --  ivl jvl = ',ivl,jvl
      do i=1,lpx(ivl,jvl)
!odl                  Write(*,*) 'QVAN3  --  i = ',i
        lp = lpl(ivl,jvl,i)
!odl                  Write(*,*) 'QVAN3  --  lp = ',lp

!     EXTRACTION OF ANGULAR MOMENT L FROM LP:

        if (lp.eq.1) then
          l = 1
        else if ((lp.ge.2) .and. (lp.le.4)) then
          l = 2
        else if ((lp.ge.5) .and. (lp.le.9)) then
          l = 3
        else if ((lp.ge.10).and.(lp.le.16)) then
          l = 4
        else if ((lp.ge.17).and.(lp.le.25)) then
          l = 5
        else if (lp.ge.26) then
          call errore(' qvan3 ',' lp.ge.26 ',lp)
        end if

        sig = (0.d0,-1.d0)**(l-1)
        sig = sig * ap(lp,ivl,jvl)

!odl                  Write(*,*) 'QVAN3  --  sig = ',sig

!        WRITE( stdout,*) 'qvan3',ng1,LP,L,ivs,jvs

          qg = qg + sig * ylm_k(lp) * qr(ivs,jvs,l,is)

      end do
      return
      end subroutine qvan3
