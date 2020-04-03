!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE add_nl_evc(c_zero)
  !----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE cell_base,            ONLY : at, bg, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE klist,                ONLY : nks, xk, ngk
  USE gvect,                ONLY : g
  USE uspp,                 ONLY : nkb, vkb, qq_at, deeq, qq_so, deeq_nc
  USE uspp_param,           ONLY : upf, nh, newpseudo, nhm
  USE wvfct,                ONLY : nbnd, npw, npwx, wg, et
  USE klist,                ONLY : igk_k
  USE lsda_mod,             ONLY : lsda, current_spin, isk, nspin
  USE symme,                ONLY : symvector
  USE wavefunctions_module, ONLY : evc
  USE noncollin_module,     ONLY : npol, noncolin
  USE spin_orb,             ONLY : lspinorb
  USE io_global,            ONLY : ionode
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE becmod,               ONLY : bec_type, becp, allocate_bec_type, deallocate_bec_type
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
  USE mp,                   ONLY : mp_sum, mp_get_comm_null
  USE becmod,  ONLY : calbec
  USE zero_mod, ONLY : ion_vel!,becpd
  !
  IMPLICIT NONE
  !
  ! ... the dummy variable
  !
  logical ::l_test
  complex(DP), allocatable  ::evc_mod(:,:,:)
  real(DP),intent(inout) :: c_zero(3)
  REAL(DP) :: forcenl(3,nat)
  ! output: the nonlocal contribution
  TYPE(bec_type) :: rdbecp (3),rdbecp_mod (nat,3), becp_mod(nat)
  ! auxiliary variable, contains <dbeta|psi>
  COMPLEX(DP), ALLOCATABLE :: vkb1(:,:)
  ! auxiliary variable contains g*|beta>
  REAL(DP) :: ps
  INTEGER       :: ik, ipol, ibnd, ibnd_loc, ig, ih, jh, na, nt, ikb, jkb, ijkb0,a
  integer ::axes,iata,b,iat
  ! counters
  !
!ciclo sulle componenti della corrente
  open(unit=16,file='rdbecp_mod(nl)->becprd(nc)',status='unknown')
  open(unit=17,file='becp_mod(nl)->becpr(nc)',status='unknown')
  open(unit=18,file='rdbecp(nl)->becpd(nc)',status='unknown')
  open(unit=19,file='becp(nl)->becp(nc)',status='unknown')
  open(unit=20,file='vkb1(nl)',status='unknown')
  open(unit=21,file='vkb(nl)',status='unknown')
  open(unit=22,file='evc(nl)',status='unknown')

  axes_cycle: do axes=1,3
!
  CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm )
  do iat=1,nat
     CALL allocate_bec_type ( nkb, nbnd, becp_mod(iat), intra_bgrp_comm )
  end do
  !
  forcenl(:,:) = 0.D0
       !
  DO ipol = 1, 3
!     CALL allocate_bec_type ( nkb, nbnd, becpd(ipol))
     CALL allocate_bec_type ( nkb, nbnd, rdbecp(ipol), intra_bgrp_comm )   
     do iat=1,nat
        CALL allocate_bec_type ( nkb, nbnd, rdbecp_mod(iat,ipol), intra_bgrp_comm )   
     end do
  END DO
  ALLOCATE( vkb1(  npwx, nkb ) ) 
!costoso per ora in termini di memory management
  ALLOCATE(evc_mod(npwx,nbnd,nat))

  CALL modify_evc(evc_mod,axes)
!
!da qui in poi evc_mod si suppone inizializzato.
       !
       ! ... the forces are a sum over the K points and the bands
       !
       DO ik = 1, nks
          IF ( lsda ) current_spin = isk(ik)
          !
          npw = ngk (ik)
!commentato io
!          IF ( nks > 1 ) THEN
!             READ( iunigk ) igk
!             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
!             IF ( nkb > 0 ) &
!                CALL init_us_2( npw, igk, xk(1,ik), vkb )
!                print*,'IGK',igk
!                print*,'XK',xk(1,ik)
!          END IF
          !
          CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
          !
!          print*,'BEFORE MODIFICATION',vkb(1:npw,1)
          
          !
          CALL calbec ( npw, vkb, evc, becp )
!
          do iat=1,nat
             CALL calbec ( npw, vkb, evc_mod(1:npwx,1:nbnd,iat), becp_mod(iat) )
          end do
!
!fin qui becp e becp_mod(iat) sono inizializzati.
!          do ig=1,npw
!             print*,ig,igk(ig)
!          end do  
!          print*,'GAMMA?',gamma_only
          DO ipol = 1, 3
             DO jkb = 1, nkb
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
                DO ig = 1, npw
                   vkb1(ig,jkb) = vkb(ig,jkb) * (0.D0,-1.D0) * g(ipol,igk_k(ig,ik))
                END DO
!$OMP END PARALLEL DO
             END DO
             !
!----------------------------------------------------------
             l_test=.true.
             if (l_test) then
                 if ((ionode).and.(axes==1)) then
                    write(20,*)'IPOL: ',ipol
                    do jkb=1,nkb
                       write(20,*) 'JKB: ',jkb
                       do ig=1,npw
                          write(20,*) ig,vkb1(ig,jkb)
                       end do 
                    end do
                 end if
                 if ((ionode).and.(axes==1).and.(ipol==1)) then
                    do jkb=1,nkb
                       write(21,*) 'JKB: ',jkb
                       do ig=1,npw
                         write(21,*) ig,vkb(ig,jkb)
                       end do
                    end do
                 end if
                  if ((ionode).and.(axes==1).and.(ipol==1)) then
                    do ibnd=1,nbnd
                       write(22,*) 'IBND: ',ibnd
                       do ig=1,npw
                         write(22,*) ig,evc(ig,ibnd)
                       end do
                    end do
                 end if

             end if
!----------------------------------------------------------
             !
             CALL calbec ( npw, vkb1, evc, rdbecp(ipol) )
!             CALL calbec ( npw, vkb1, evc, becpd(ipol) )

             do iat=1,nat
                CALL calbec ( npw, vkb1, evc_mod(1:npwx,1:nbnd,iat) , rdbecp_mod(iat,ipol) )
             end do
             !
          END DO
!!fin qui rdbecp e rdbecp_mod(iat,ipol) sono inizializzati.
!--------------------------------------------------------------------
     l_test=.true.
     if (l_test) then 
          if (ionode) then
              do ibnd=1,nbnd
                 do ikb=1,nkb
                    do b=1,3
                       write(16, "('BECPRD: A-B-IKB-IBND..value',4I4,F12.7)")&
axes,b,ikb,ibnd,rdbecp_mod(1,b)%r(ikb,ibnd)
                    end do      
                 end do
              end do
          end if
!confronto con becpr

           if (ionode) then
              do ibnd=1,nbnd
                 do ikb=1,nkb
                    write(17,"('BECPR: A-IKB-IBND..value',3I4,F12.7)")&
axes,ikb,ibnd,becp_mod(1)%r(ikb,ibnd)
                 end do
              end do
          end if
!confrondo per becpd
           if ((ionode).and.(axes==1)) then
              do ibnd=1,nbnd
                 do ikb=1,nkb
                    do ipol=1,3
                        write(18, "('BECPR: A-IKB-IBND..value',3I4,2F12.7)")&
ipol,ikb,ibnd,rdbecp(ipol)%r(ikb,ibnd)!,becpd(ipol)%r(ikb,ibnd)
                    end do  
                 end do
              end do
          end if
!confronto per becp 
          if ((ionode).and.(axes==1)) then
              do ibnd=1,nbnd
                 do ikb=1,nkb
                     write(19,"('BECPR: IKB-IBND..value',2I4,F12.7)")&
ikb,ibnd,becp%r(ikb,ibnd) 
                 end do
              end do
          end if
     end if
!--------------------------------------------------------------
          !
          ijkb0 = 0
          DO nt = 1, ntyp
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ih = 1, nh(nt)
                      ikb = ijkb0 + ih
                      DO ibnd_loc = 1, becp%nbnd_loc
                         ibnd = ibnd_loc + becp%ibnd_begin - 1
                         ps = deeq(ih,ih,na,current_spin) - &
                              et(ibnd,ik) * qq_at(ih,ih,na)
                         DO ipol = 1, 3
!                              forcenl(ipol,na) = forcenl(ipol,na) - &
!                                        ps * wg(ibnd,ik) * 2.D0 * tpiba * &
!                                        rdbecp(ipol)%r(ikb,ibnd_loc) *becp%r(ikb,ibnd_loc)
                                       
!cambiato così
                             forcenl(ipol,na) = forcenl(ipol,na) + &
                                       ps * wg(ibnd,ik) * tpiba * &
                                       rdbecp_mod(na,ipol)%r(ikb,ibnd_loc) *becp%r(ikb,ibnd_loc) + &                                        
                                       ps * wg(ibnd,ik) * tpiba * &
                                       rdbecp(ipol)%r(ikb,ibnd_loc) *becp_mod(na)%r(ikb,ibnd_loc)
!
                         END DO

                      END DO
                      !
                      IF ( upf(nt)%tvanp .OR. newpseudo(nt) ) THEN
                         !
                         ! ... in US case there is a contribution for jh<>ih. 
                         ! ... We use here the symmetry in the interchange 
                         ! ... of ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(nt)
                            jkb = ijkb0 + jh
                            DO ibnd_loc = 1, becp%nbnd_loc
                               ibnd = ibnd_loc + becp%ibnd_begin - 1
                               ps = deeq(ih,jh,na,current_spin) - &
                                    et(ibnd,ik) * qq_at(ih,jh,na)
                               DO ipol = 1, 3
                                  forcenl(ipol,na) = forcenl(ipol,na) - &
                                     ps * wg(ibnd,ik) * 2.d0 * tpiba * &
                                     (rdbecp(ipol)%r(ikb,ibnd_loc) *becp%r(jkb,ibnd_loc) + &
                                      rdbecp(ipol)%r(jkb,ibnd_loc) *becp%r(ikb,ibnd_loc) )
                               END DO
                            END DO
                         END DO
                      END IF
                   END DO
                   ijkb0 = ijkb0 + nh(nt)
                END IF
             END DO
          END DO
       END DO
       !
       IF( becp%comm /= mp_get_comm_null() ) CALL mp_sum( forcenl, becp%comm )
       !
       ! ... The total D matrix depends on the ionic position via the
       ! ... augmentation part \int V_eff Q dr, the term deriving from the 
       ! ... derivative of Q is added in the routine addusforce
       !
!         if (ionode) print*,'CZERO_P',c_zero(axes),axes,ion_vel,forcenl
!         do na=1,nat
!            c_zero(axes)=c_zero(axes)+DOT_PRODUCT(ion_vel(1:3,na),forcenl(1:3,na))
!         end do
!         if (ionode) print*, 'CZERO_D',c_zero(axes),axes


!cosa è? dovrebbe essere il termine prettamente ultrasoffice.
!        CALL addusforce( forcenl )
!tolto
         do na=1,nat
            c_zero(axes)=c_zero(axes)+DOT_PRODUCT(ion_vel(1:3,na),forcenl(1:3,na))
         end do
       !
       !
       ! ... collect contributions across pools
       !
       CALL mp_sum( forcenl, inter_pool_comm )
       !
       ! ... Since our summation over k points was only on the irreducible 
       ! ... BZ we have to symmetrize the forces
       !
       CALL symvector ( nat, forcenl )
       !
       DEALLOCATE( vkb1 )
       DEALLOCATE( evc_mod)
       DO ipol = 1, 3
          CALL deallocate_bec_type ( rdbecp(ipol) )
!          CALL deallocate_bec_type (becpd(ipol) )   
          do iat=1,nat
             CALL deallocate_bec_type ( rdbecp_mod(iat,ipol) )   
          end do
       END DO
       !
  !
  do iat=1,nat
     CALL deallocate_bec_type ( becp_mod(iat) )
  end do
  CALL deallocate_bec_type ( becp )

!calcolo ora la corrente di componente "axes" con il termine forcenl qui calcolato
!anche qui per ora funziona solo per cella cubica credo, perchè ion_vel viene dato in coordinate
!(x,y,z)
!  print*,'CZERO_P',c_zero(axes),axes,ion_vel,forcenl
!  do na=1,nat   
!     c_zero(axes)=c_zero(axes)+DOT_PRODUCT(ion_vel(1:3,na),forcenl(1:3,na))
!  end do
!  print*, 'CZERO_D',c_zero(axes),axes
!
  end do axes_cycle
  close(16)
  close(17)
  close(18)
  close(19)
  close(20)
  close(21)
  close(22)
!
END subroutine add_nl_evc
