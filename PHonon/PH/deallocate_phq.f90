!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------
subroutine deallocate_phq
!----------========-----------------------
!
!  deallocates the variables allocated by allocate_phq
!
  USE noncollin_module, ONLY : m_loc
  USE becmod, ONLY: bec_type, becp, deallocate_bec_type
  USE wavefunctions_module,  ONLY: evc

  USE ramanm, ONLY: ramtns
  USE modes, ONLY : tmq, t, npert, u, name_rap_mode, num_rap_mode
  USE efield_mod, ONLY : zstareu, zstarue, zstarue0, zstareu0, &
                         zstarue0_rec
  USE phus, ONLY : int1, int1_nc, int2, int2_so, &
                   int4, int4_nc, int5, int5_so, becsum_nc, &
                   becsumort, alphasum, alphasum_nc, &
                   alphap
  USE gamma_gamma, ONLY : with_symmetry, has_equivalent, equiv_atoms, &
                   n_equiv_atoms
  USE nlcc_ph, ONLY : drc
  USE units_ph, ONLY : this_dvkb3_is_on_file, this_pcxpsi_is_on_file
  USE dynmat, ONLY : dyn00, dyn_rec, dyn, w2
  USE el_phon, ONLY : el_ph_mat
  USE freq_ph, ONLY : polar

  USE lrus,         ONLY : int3, int3_nc, int3_paw, becp1, dpqq, dpqq_so
  USE lr_symm_base, ONLY : rtau
  USE gc_lr,        ONLY : grho, gmag, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s, &
                           vsgga, segni
  USE qpoint,       ONLY : eigqts, ikks, ikqs, nksq, xk_col
  USE eqv,          ONLY : dmuxc, vlocq, dpsi, dvpsi, evq
  USE control_lr,   ONLY : lgamma, nbnd_occ

  IMPLICIT NONE
  INTEGER :: ik, ipol

  if(allocated(ramtns)) deallocate (ramtns)
  if (lgamma) then
     if(associated(evq)) nullify(evq)
  else
     if(associated(evq)) deallocate(evq)
  end if

  if(allocated(dvpsi)) deallocate (dvpsi)
  if(allocated(dpsi)) deallocate ( dpsi)
  !
  if(allocated(vlocq)) deallocate (vlocq)
  if(allocated(dmuxc)) deallocate (dmuxc)
  !
  if(allocated(ikks)) deallocate (ikks)
  if(allocated(ikqs)) deallocate (ikqs)
  if(allocated(eigqts)) deallocate (eigqts)
  if(allocated(rtau)) deallocate (rtau)
  if(associated(u)) deallocate (u)
  if(allocated(name_rap_mode)) deallocate (name_rap_mode)
  if(allocated(num_rap_mode)) deallocate (num_rap_mode)
  if(allocated(dyn)) deallocate (dyn)
  if(allocated(dyn_rec)) deallocate (dyn_rec)
  if(allocated(dyn00)) deallocate (dyn00)
  if(allocated(w2)) deallocate (w2)
  if(allocated(xk_col)) deallocate (xk_col)
  if(allocated(polar)) deallocate (polar)

  CALL deallocate_pert()

  if(allocated(npert)) deallocate (npert)
  if(allocated(zstareu)) deallocate (zstareu)
  if(allocated(zstareu0)) deallocate (zstareu0)
  if(allocated(zstarue)) deallocate (zstarue)
  if(allocated(zstarue0)) deallocate (zstarue0)
  if(allocated(zstarue0_rec)) deallocate (zstarue0_rec)

  if(allocated(int1)) deallocate (int1)
  if(allocated(int2)) deallocate (int2)
  if(allocated(int3)) deallocate (int3)
  if(allocated(int3_paw)) deallocate (int3_paw)
  if(allocated(int4)) deallocate (int4)
  if(allocated(int5)) deallocate (int5)
  if(allocated(dpqq)) deallocate (dpqq)
  if(allocated(int1_nc)) deallocate(int1_nc)
  if(allocated(int3_nc)) deallocate(int3_nc)
  if(allocated(int4_nc)) deallocate(int4_nc)
  if(allocated(becsum_nc)) deallocate(becsum_nc)
  if(allocated(becsumort)) deallocate(becsumort)
  if(allocated(alphasum_nc)) deallocate(alphasum_nc)
  if(allocated(int2_so)) deallocate(int2_so)
  if(allocated(int5_so)) deallocate(int5_so)
  if(allocated(dpqq_so)) deallocate(dpqq_so)

  if(allocated(alphasum)) deallocate (alphasum)
  if(allocated(this_dvkb3_is_on_file)) deallocate (this_dvkb3_is_on_file)

  if(allocated(this_pcxpsi_is_on_file)) deallocate (this_pcxpsi_is_on_file)
  if(allocated(alphap)) then
     do ik=1,nksq
        do ipol=1,3
           call deallocate_bec_type ( alphap(ipol,ik) )
        enddo
     end do
     deallocate (alphap)
  endif
  if(allocated(becp1))  then
     do ik=1,size(becp1)
        call deallocate_bec_type ( becp1(ik) )
     end do
     deallocate(becp1)
  end if
  call deallocate_bec_type ( becp )

  if(allocated(el_ph_mat)) deallocate (el_ph_mat)
  if(allocated(m_loc))     deallocate(m_loc)

  if(allocated(drc)) deallocate(drc)

  if(allocated(dvxc_rr)) deallocate (dvxc_rr)
  if(allocated(dvxc_sr)) deallocate (dvxc_sr)
  if(allocated(dvxc_ss)) deallocate (dvxc_ss)
  if(allocated(dvxc_s)) deallocate (dvxc_s)
  if(allocated(grho)) deallocate (grho)
  if(allocated(segni)) deallocate (segni)
  if(allocated(vsgga)) deallocate (vsgga)
  if(allocated(gmag))  deallocate (gmag)

  IF (allocated(has_equivalent))   DEALLOCATE(has_equivalent)
  IF (allocated(with_symmetry))    DEALLOCATE(with_symmetry)
  IF (allocated(n_equiv_atoms))    DEALLOCATE(n_equiv_atoms)
  IF (allocated(equiv_atoms))      DEALLOCATE(equiv_atoms)

  IF (allocated(nbnd_occ))         DEALLOCATE(nbnd_occ)

  return
end subroutine deallocate_phq
