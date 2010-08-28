!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine addusdynmat (dynwrk)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the additional terms which are contained in
  !     <psi|V"|psi> part of the dynamical matrix and which are due
  !     to the change of the self consistent D term in the pseudopotential
  !     There are four additional terms which we compute here.
  !

  USE kinds, only : DP
  USE ions_base, ONLY : nat, ityp
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp, ONLY: okvan, becsum
  USE uspp_param, only: upf, nh
  USE lsda_mod, ONLY : nspin
  USE spin_orb, ONLY : lspinorb
  USE noncollin_module, ONLY : nspin_lsda

  USE phus,    ONLY : int1, int1_nc, int2, int2_so, int4, int4_nc, &
                      int5, int5_so, alphasum, alphasum_nc, becsum_nc
  USE modes,   ONLY : nmodes

  implicit none

  complex(DP) :: dynwrk (3 * nat, 3 * nat)
  ! inp/out: the dynamical matrix

  integer :: ipol, jpol, np, na, nb, nu_i, nu_j, ih, jh, ijh, dim, &
       is, is1, is2, ijs
  ! counter on polarizations
  ! counter on pseudopotentials
  ! counter on atoms
  ! counter on modes
  ! counter on solid beta functions
  ! composed dimension of the beta
  ! counter on spin

  complex(DP) :: term (3, 3), dyn1 (3 * nat, 3 * nat)
  ! auxiliary space
  ! auxiliary dynamical matrix


  if (.not.okvan) return
  call start_clock ('addusdynmat')

  IF (noncolin) CALL set_int12_nc(1)

  dyn1 (:,:) = (0.d0, 0.d0)
  !
  !  We compute the four terms required
  !
  do na = 1, nat
     np = ityp (na)
     if (upf(np)%tvanp  ) then
        dim = (nh (np) * (nh (np) + 1) ) / 2
        do ipol = 1, 3
           nu_i = 3 * (na - 1) + ipol
           do jpol = 1, 3
              nu_j = 3 * (na - 1) + jpol
              IF (noncolin) THEN
                 ijh=1
                 DO ih=1,nh(np)
                    DO jh=ih,nh(np)
                       ijs=0
                       DO is1=1,npol
                          DO is2=1,npol
                             ijs=ijs+1
                             dynwrk(nu_i, nu_j)=dynwrk(nu_i, nu_j)  + &
                                  int4_nc(ih,jh,ipol,jpol,na,ijs)   * &
                                  becsum_nc(ijh,na,is1,is2)
                             IF (ih.NE.jh) THEN
                                dynwrk(nu_i, nu_j)=dynwrk(nu_i, nu_j) + &
                                     int4_nc(jh,ih,ipol,jpol,na,ijs)  * &
                                     CONJG(becsum_nc(ijh,na,is2,is1))
                             END IF
                          END DO
                       END DO
                       ijh=ijh+1
                    END DO
                 END DO
              ELSE
                 do is = 1, nspin
                    do ijh = 1, dim
                       dynwrk(nu_i, nu_j)=dynwrk(nu_i, nu_j)+ &
                           int4(ijh,ipol,jpol,na,is) * becsum(ijh,na,is)
                    enddo
                 enddo
              END IF
           enddo
        enddo
        !
        !   The second term requires an exchange of the components.
        !
        term (:,:) = (0.d0, 0.d0)
        do ipol = 1, 3
           do jpol = 1, 3
              ijh = 0
              do ih = 1, nh (np)
                 do jh = ih, nh (np)
                    ijh = ijh + 1
                    IF (noncolin) THEN
                       ijs=0
                       do is1 = 1, npol
                          do is2 = 1, npol
                             ijs=ijs+1
                             term(ipol,jpol) = term(ipol,jpol) + &
                             int1_nc(ih,jh,ipol,na,ijs)*     &
                             alphasum_nc(ijh,jpol,na,is1,is2)
                             IF (ih.ne.jh) THEN
                                term(ipol,jpol) = term(ipol,jpol) + &
                                    int1_nc(jh,ih,ipol,na,ijs)* &
                                    CONJG(alphasum_nc(ijh,jpol,na,is2,is1))
                             ENDIF
                          enddo
                       enddo
                    ELSE
                       do is = 1, nspin
                          term(ipol,jpol) = term(ipol,jpol) + &
                          CONJG(int1(ih,jh,ipol,na,is))*alphasum(ijh,jpol,na,is)
                       enddo
                    END IF
                 enddo
              enddo
           enddo
        enddo
        !
        !  And then we add the appropriate terms to the dynamical matrix
        !
        do ipol = 1, 3
           nu_i = 3 * (na - 1) + ipol
           do jpol = 1, 3
              nu_j = 3 * (na - 1) + jpol
              dynwrk (nu_i, nu_j) = dynwrk (nu_i, nu_j) + &
                                    term (ipol, jpol) + term (jpol, ipol)
           enddo
        enddo
        !
        !   the other two terms do not contain a delta ss'
        !
        do nb = 1, nat
           do ipol = 1, 3
              nu_i = 3 * (nb - 1) + ipol
              do jpol = 1, 3
                 nu_j = 3 * (na - 1) + jpol
                 ijh = 0
                 do ih = 1, nh (np)
                    do jh = ih, nh (np)
                       ijh = ijh + 1
                       IF (lspinorb) THEN
                          ijs=0
                          do is1 = 1, npol
                             do is2 = 1, npol
                                ijs=ijs+1
                                dyn1(nu_i,nu_j)=dyn1(nu_i,nu_j) + &
                                     int2_so(ih,jh,ipol,nb,na,ijs) * &
                                     alphasum_nc(ijh,jpol,na,is1,is2)   + &
                                     int5_so(ih,jh,ipol,jpol,nb,na,ijs) * &
                                     becsum_nc(ijh,na,is1,is2)
                                IF (ih.ne.jh) THEN
                                   dyn1(nu_i,nu_j)=dyn1(nu_i,nu_j) + &
                                    int2_so(jh,ih,ipol,nb,na,ijs) * &
                                    CONJG(alphasum_nc(ijh,jpol,na,is2,is1))+&
                                    int5_so(jh,ih,ipol,jpol,nb,na,ijs) * &
                                    CONJG(becsum_nc(ijh,na,is2,is1))
                                END IF
                             enddo
                          enddo
                       ELSE
                          do is = 1, nspin_lsda
                             dyn1(nu_i,nu_j)=dyn1(nu_i,nu_j) + &
                                          CONJG(int2(ih,jh,ipol,nb,na)) * &
                                               alphasum(ijh,jpol,na,is) + &
                                          int5(ijh,ipol,jpol,nb,na) *     &
                                                  becsum(ijh,na,is)
                          enddo
                       END IF
                    enddo
                 enddo
              enddo
           enddo
        enddo
     endif

  enddo
  do nu_i = 1, nmodes
     do nu_j = 1, nmodes
        dynwrk (nu_i, nu_j) = dynwrk (nu_i, nu_j) + &
                              dyn1 (nu_i, nu_j) + CONJG(dyn1 (nu_j, nu_i) )
    enddo
  enddo
  deallocate (int4)
  deallocate (int5)

  IF (noncolin) THEN
     call set_int12_nc(0)
     deallocate(int4_nc)
     if (lspinorb) deallocate(int5_so)
  END IF


  call stop_clock ('addusdynmat')
  return
end subroutine addusdynmat
