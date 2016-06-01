!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
subroutine dvkb3(ik,dvkb)
!----------=========-------------------------------------------------------
!
!
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, tpiba
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE gvect,     ONLY : g
  USE lsda_mod,  ONLY : lsda, current_spin, isk
  USE klist,     ONLY : xk, ngk, igk_k
  USE wvfct,     ONLY : npwx
  USE wavefunctions_module,    ONLY : evc
  USE uspp,      ONLY: nkb
  USE uspp_param,ONLY: nh
  USE units_ph,  ONLY: this_dvkb3_is_on_file, lrdvkb3, iudvkb3
  implicit none

  integer, intent(in) :: ik
  complex(DP), intent(out) :: dvkb (npwx,nkb,3)

  integer :: npw, jpol,  nt, na, ikb, jkb, ig

  real(DP), allocatable  :: gk (:,:), g2kin(:)

  complex(DP), allocatable :: work (:,:)

  if (this_dvkb3_is_on_file(ik)) then
     call davcio (dvkb, lrdvkb3, iudvkb3, ik, -1)
  else
     npw = ngk(ik)
     allocate (work(npwx,nkb))
     allocate (gk(3, npw))
     allocate (g2kin(npw))
     !
     do ig = 1, npw
        gk (1, ig) = (xk (1, ik) + g (1, igk_k (ig,ik) ) ) * tpiba
        gk (2, ig) = (xk (2, ik) + g (2, igk_k (ig,ik) ) ) * tpiba
        gk (3, ig) = (xk (3, ik) + g (3, igk_k (ig,ik) ) ) * tpiba
        g2kin (ig) = gk (1, ig) **2 + gk (2, ig) **2 + gk (3, ig) **2
        if (g2kin (ig) .lt.1.0d-10) then
           gk (1, ig) = 0.d0
           gk (2, ig) = 0.d0
           gk (3, ig) = 0.d0
        else
           gk (1, ig) = gk (1, ig) / sqrt (g2kin (ig) )
           gk (2, ig) = gk (2, ig) / sqrt (g2kin (ig) )
           gk (3, ig) = gk (3, ig) / sqrt (g2kin (ig) )
        endif
     enddo

     if (lsda) current_spin = isk (ik)
     do jpol=1,3
        call gen_us_dy (ik, at (1, jpol), dvkb(1,1,jpol))
     end do
     call gen_us_dj (ik, work)

     jkb = 0
     do nt = 1, ntyp
        do na = 1, nat
           if (nt.eq.ityp (na)) then
              do ikb = 1, nh (nt)
                 jkb = jkb + 1
                 do jpol=1,3
                    do ig = 1, npw
                       dvkb(ig,jkb,jpol) = dvkb(ig,jkb,jpol) + work(ig, jkb) * &
                                                  (at (1, jpol) * gk (1, ig) + &
                                                   at (2, jpol) * gk (2, ig) + &
                                                   at (3, jpol) * gk (3, ig) )
                    enddo
                 enddo
              enddo
           endif
        enddo
     enddo

     deallocate(g2kin)
     deallocate(gk)
     deallocate(work)
     call davcio (dvkb, lrdvkb3, iudvkb3, ik, 1)
     this_dvkb3_is_on_file(ik) = .true.
  end if

  return

end subroutine dvkb3
