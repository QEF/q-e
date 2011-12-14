!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine init_ns
   !-----------------------------------------------------------------------
   !
   ! This routine computes the starting ns (for lda+U calculation) filling
   ! up the d states (the only interested by the on-site potential for the
   ! moment) according to the Hund's rule (valid for the isolated atoms on
   ! which starting potential is built), and to the starting_magnetization:
   ! majority spin levels are populated first, then the remaining electrons
   ! are equally distributed among the minority spin states
   !
   USE kinds,     ONLY : DP
   USE ions_base, ONLY : nat, ityp
   USE lsda_mod,  ONLY : nspin, starting_magnetization
   USE ldaU,      ONLY : hubbard_u, hubbard_alpha, hubbard_l
   USE scf,       ONLY : rho
   USE uspp_param,ONLY : upf
   !
   implicit none

   real(DP) :: totoc
   real(DP), external :: hubbard_occ

   integer :: na, nt, is, m1, majs, mins
   logical :: nm        ! true if the atom is non magnetic

   rho%ns(:,:,:,:) = 0.d0

   do na = 1, nat
      nt = ityp (na)
      if (Hubbard_U(nt).ne.0.d0 .or. Hubbard_alpha(nt).ne.0.d0) then
         totoc = hubbard_occ ( upf(nt)%psd )
         nm=.true.
         if (nspin.eq.2) then
            if (starting_magnetization (nt) .gt.0.d0) then  
               nm=.false.
               majs = 1  
               mins = 2  
            elseif (starting_magnetization (nt) .lt.0.d0) then  
               nm=.false.
               majs = 2  
               mins = 1  
            endif  
         endif
         if (.not.nm) then  
            if (totoc.gt.2*Hubbard_l(nt)+1) then  
               do m1 = 1, 2*Hubbard_l(nt)+1  
                  rho%ns (m1, m1, majs, na) = 1.d0  
                  rho%ns (m1, m1, mins, na) = (totoc -(2*Hubbard_l(nt)+1) ) / &
                                                      (2*Hubbard_l(nt)+1)  
               enddo  
            else  
               do m1 = 1, 2*Hubbard_l(nt)+1  
                  rho%ns (m1, m1, majs, na) = totoc / (2*Hubbard_l(nt)+1)
               enddo  
            endif  
         else  
            do is = 1,nspin
               do m1 = 1, 2*Hubbard_l(nt)+1  
                  rho%ns (m1, m1, is, na) = totoc /  2.d0 / (2*Hubbard_l(nt)+1)
               enddo  
            enddo  
         endif  
      endif  
   enddo  
   return  
end subroutine init_ns
