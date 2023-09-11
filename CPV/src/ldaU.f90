!
! Copyright (C) 2011-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
      SUBROUTINE s_wfc(nwfc,becwfc,betae,wfc,swfc)
      !-----------------------------------------------------------------------
      !! Take as input the wfc, \(\text{becwfc}=\langle\text{wfc}|\text{beta}\rangle\),
      !! \(\text{betae}=|\text{beta}\rangle\).
      !! Gives as output: \(\text{swfc}=S|\text{wfc}\rangle\)
      !
      USE kinds, ONLY: DP
      USE ions_base, ONLY: nat, ityp
      USE uspp, ONLY: nkb, nkbus, qq_nt, ofsbeta
      USE uspp_param, ONLY: nh, upf
      USE gvecw, ONLY: ngw
      IMPLICIT NONE
! input
      INTEGER, INTENT(in)     :: nwfc
      COMPLEX(DP), INTENT(in) :: betae(ngw,nkb),                   &
     &                             wfc(ngw,nwfc)
      REAL(DP), INTENT(in)    :: becwfc(nkb,nwfc)
! output
      COMPLEX(DP), INTENT(out):: swfc(ngw,nwfc)
! local
      INTEGER :: is, iv, jv, ia, inl, jnl, i
      REAL(DP) :: qtemp(nkb,nwfc)
!
      swfc = wfc
!
      IF (nkbus > 0) THEN
         qtemp=0.d0
         DO ia=1,nat
            is=ityp(ia)
            IF( upf(is)%tvanp ) THEN
               DO iv=1,nh(is)
                  DO jv=1,nh(is)
                     IF(ABS(qq_nt(iv,jv,is)).GT.1.e-5) THEN
                        inl = ofsbeta(ia) + iv
                        jnl = ofsbeta(ia) + jv
                        DO i=1,nwfc
                           qtemp(inl,i) = qtemp(inl,i) + qq_nt(iv,jv,is)*becwfc(jnl,i)
                        END DO
                     ENDIF
                  END DO
               END DO
            END IF
         END DO
!
         CALL dgemm &
              ('N','N',2*ngw,nwfc,nkb,1.0d0,betae,2*ngw,&
               qtemp,nkb,1.0d0,swfc,2*ngw)
!
      END IF
!
      RETURN
      END SUBROUTINE s_wfc
!-----------------------------------------------------------------------
      subroutine ldaU_init
!-----------------------------------------------------------------------
!
      use ldaU_cp,          ONLY: nwfcU, lda_plus_u, Hubbard_U
      use ldaU_cp,          ONLY: Hubbard_lmax, Hubbard_l, ldmx, ns, vupsi
      use ions_base,        only: nsp, atm, nat
      use gvecw,            only: ngw
      use electrons_base,   only: nspin, nx => nbspx
      USE uspp_param,       ONLY: upf
      !
      implicit none
      integer is, nb, l

      IF ( .NOT.lda_plus_u ) RETURN
      ! FIXME: wasteful allocation, should be removed
      allocate(vupsi(ngw,nx))
      vupsi=(0.0d0,0.0d0)

      Hubbard_lmax = -1
      do is=1,nsp
         if (Hubbard_U(is).ne.0.d0) then 
            ! Hubbard_l is read from the HUBBARD card in the input file
            Hubbard_lmax = max(Hubbard_lmax,Hubbard_l(is))
            write (6,*) ' HUBBARD L FOR TYPE ',atm(is),' IS ', Hubbard_l(is)
         else
            Hubbard_l(is) = -1
         end if
      end do
      write (6,*) ' MAXIMUM HUBBARD L IS ', Hubbard_lmax
      if (Hubbard_lmax.eq.-1) call errore                            &
     &        ('setup','lda_plus_u calculation but Hubbard_l not set',1)
      !
      ldmx = 2 * Hubbard_lmax + 1
      allocate(ns(ldmx,ldmx,nspin,nat))
      !
      return
      end subroutine ldaU_init
!
!-----------------------------------------------------------------------
      subroutine new_ns( c, eigr, betae, hpsi, forceh )
!-----------------------------------------------------------------------
      !! This routine computes the on site occupation numbers of the Hubbard ions.
      !! It also calculates the contribution of the Hubbard Hamiltonian to the
      !! electronic potential and to the forces acting on ions.

      use kinds,              ONLY: DP        
      use control_flags,      ONLY: tfor, tprnfor
      use ions_base,          only: nat, nsp, ityp
      use gvecw,              only: ngw
      use gvect,              only: gstart
      USE uspp,               ONLY: nkb
      USE uspp_param,         ONLY: upf
      use electrons_base,     only: nspin, n => nbsp, nx => nbspx, ispin, f
      USE ldaU_cp,            ONLY: Hubbard_U, Hubbard_l, ldmx
      USE ldaU_cp,            ONLY: nwfcU, ns, e_hubbard
      USE step_penalty,       ONLY: penalty_e, penalty_f
      USE mp_pools,           ONLY: intra_pool_comm, me_pool, nproc_pool
      USE mp_bands,           only: nbgrp
      USE cp_interfaces,      only: calbec, nlsm2_bgrp
!
      implicit none
      complex(DP), intent(in) :: c(ngw,nx), eigr(ngw,nat)
      complex(DP), intent(inout) :: betae(ngw,nkb)
      complex(DP), intent(out) :: hpsi(ngw,nx)
      real(DP), INTENT(OUT) :: forceh(3,nat)
!
      complex(DP), allocatable:: wfcU(:,:), swfc(:,:), spsi(:,:)
      real(DP), allocatable   :: becwfc(:,:), bp(:,:), dbp(:,:,:), wdb(:,:,:)
      real(DP), allocatable   :: dns(:,:,:,:), proj(:,:), tempsi(:,:)
      integer is, ia, nb, isp, l, m, m1, m2, k, i, ldim, ig
      integer iv, jv, inl, jnl,alpha_a,alpha_s,ipol
      integer, allocatable ::  offset (:)
      INTEGER :: nb_s, nb_e, mykey
!
      if( nbgrp > 1 ) call errore(' new_ns ', &
         ' parallelization over bands not yet implemented ', 1 )
      call start_clock('new_ns')
!
      allocate(offset(nat))
      offset(:) = -1
      ! offset = -1 means "not a Hubbard wfc"
      nwfcU = 0
      do ia = 1, nat
         is = ityp(ia)
         do i = 1, upf(is)%nwfc
            l = upf(is)%lchi(i)
            if (l == Hubbard_l(is)) then 
               offset (ia) = nwfcU
               nwfcU = nwfcU + 2 * l + 1
            end if
         end do   
      end do
      !
      allocate(wfcU(ngw,nwfcU))
      allocate(becwfc(nkb,nwfcU))
      allocate(swfc(ngw,nwfcU))
      allocate(proj(nwfcU,n))
      !
      ! calculate proj = <wfcU|S|c>
      !
      CALL projwfc_hub( c, nx, eigr, betae, n, nwfcU, &
     &                  offset, Hubbard_l, wfcU, becwfc, swfc, proj )
      !
      ns(:,:,:,:) = 0.d0
      do ia = 1,nat
         is = ityp(ia)
         if (Hubbard_U(is).ne.0.d0) then 
            k = offset(ia)
            do m1 = 1, 2*Hubbard_l(is) + 1
               do m2 = m1, 2*Hubbard_l(is) + 1
                  do i = 1,n
                   ns(m1,m2,ispin(i),ia) = ns(m1,m2,ispin(i),ia) + &
     &                            f(i) * proj(k+m2,i) * proj(k+m1,i)
                  end do
               end do
               do m2 = m1+1, 2*Hubbard_l(is) + 1
                  ns(m2,m1,:,ia) = ns(m1,m2,:,ia)
               end do
            end do
         end if
      end do
      if (nspin.eq.1) ns = 0.5d0 * ns
! Contributions to total energy
      e_hubbard = 0.d0
      do ia = 1,nat
         is = ityp(ia)
         if (Hubbard_U(is).ne.0.d0) then
             do isp = 1,nspin
                do m1 = 1, 2*Hubbard_l(is) + 1
                  e_hubbard = e_hubbard + 0.5d0 * Hubbard_U(is) *    &
     &                        ns(m1,m1,isp,ia)
                  do m2 = 1, 2*Hubbard_l(is) + 1
                     e_hubbard = e_hubbard - 0.5d0 * Hubbard_U(is) * &
     &                           ns(m1,m2,isp,ia) * ns(m2,m1,isp,ia)
                  end do
               end do
            end do
         end if
      end do
      if (nspin.eq.1) e_hubbard = 2.d0*e_hubbard
!
!      Calculate the potential and forces on wavefunctions due to U
!
      hpsi(:,:)=(0.d0,0.d0)
      ALLOCATE ( tempsi(ldmx,n) )
      tempsi(:,:)=(0.d0,0.d0)
      do ia = 1,nat
         is = ityp(ia)
         if (Hubbard_U(is).ne.0.d0) then
            ldim = 2*Hubbard_l(is) + 1
            do i=1, n
               do m1 = 1, ldim
                  tempsi(m1,i) = proj (offset(ia)+m1,i)
                  do m2 = 1, ldim
                     tempsi(m1,i) = tempsi(m1,i) - &
                                    2.0_dp*ns(m1,m2,ispin(i),ia) * &
                                            proj (offset(ia)+m2,i)
                  enddo
                  tempsi(m1,i) = tempsi(m1,i) * Hubbard_U(is)/2.d0*f(i)
               enddo
            enddo
            !
            CALL dgemm ( 'N','N', 2*ngw, n, ldim, 1.0_dp, &
                          swfc(1,offset(ia)+1), 2*ngw, tempsi, &
                          ldmx, 1.0_dp, hpsi, 2*ngw )
         endif
      enddo
      DEALLOCATE ( tempsi )
!
!      Calculate the potential and energy due to constraint
!
      CALL  penalty_e ( offset, swfc, proj, e_hubbard, hpsi )
!
! Calculate the contribution to forces on ions due to U and constraint
!
      forceh=0.d0

      if ( tfor .or. tprnfor ) then
        call start_clock('new_ns:forc')
        allocate (bp(nkb,n), dbp(nkb,nx,3), wdb(nkb,nwfcU,3))
        allocate(dns(ldmx,ldmx,nspin,nat))
        allocate (spsi(ngw,n))
!
        call calbec ( n, betae, c, bp )
        call s_wfc ( n, bp, betae, c, spsi )
        call nlsm2_bgrp( ngw, nkb, betae, c, dbp, nx, n )
        call nlsm2_bgrp( ngw, nkb, betae, wfcU, wdb, nwfcU, nwfcU )
        !
        ! poor-man parallelization over bands
        ! - if nproc_pool=1   : nb_s=1, nb_e=n, mykey=0
        ! - if nproc_pool<=nbnd:each processor calculates band nb_s to nb_e; mykey=0
        ! - if nproc_pool>nbnd :each processor takes care of band nb_s=nb_e;
        !   mykey labels how many times each band appears (mykey=0 first time etc.)
        !
        CALL block_distribute( n, me_pool, nproc_pool, nb_s, nb_e, mykey )
        ! 
        do alpha_a = 1, nat
            alpha_s = ityp(alpha_a)
            do ipol = 1,3
               call dndtau(alpha_a,alpha_s,becwfc,spsi,bp,dbp,wdb,          &
                           offset,wfcU,eigr,proj,ipol,nb_s,nb_e,mykey,&
                           dns)
               do ia=1, nat
                  is = ityp(ia)
                  if (Hubbard_U(is).ne.0.d0) then
                     do isp = 1,nspin
                        do m2 = 1,2*Hubbard_l(is) + 1
                           forceh(ipol,alpha_a) = forceh(ipol,alpha_a) -   &
     &                     Hubbard_U(is) * 0.5d0 * dns(m2,m2,isp,ia)
                           do m1 = 1,2*Hubbard_l(is) + 1
                              forceh(ipol,alpha_a) = forceh(ipol,alpha_a) + &
     &                        Hubbard_U(is)*ns(m2,m1,isp,ia)*       &
     &                        dns(m1,m2,isp,ia)
                           end do
                        end do
                     end do
                  end if
! Occupation constraint added here to forceh(ipol,alpha)
                  CALL penalty_f ( is, ia, dns, forceh(ipol,alpha_a) )
               end do
            end do
        end do
        !
        ! I am not sure why the following instruction (present in PW)
        ! seems to yield a wrong factor here ... PG
        !if (nspin.eq.1) then
        !   forceh = 2.d0 * forceh
        !end if
        !
        deallocate ( spsi, dns, bp, dbp, wdb)
        call stop_clock('new_ns:forc')
      end if
      !
      deallocate ( wfcU, becwfc, proj, offset, swfc)
      !
      call stop_clock('new_ns')
      !
      return
      end subroutine new_ns
!
!-----------------------------------------------------------------------
      subroutine write_ns
!-----------------------------------------------------------------------
!
! This routine computes the occupation numbers on atomic orbitals.
! It also write the occupation number in the output file.
!
      USE kinds,            only: DP
      USE constants,        ONLY: autoev
      use electrons_base,   only: nspin
      use electrons_base,   only: n => nbsp 
      use ions_base,        only: nat, nsp, ityp
      use gvecw,            only: ngw
      USE ldaU_cp,          ONLY: Hubbard_U, Hubbard_l, ldmx
      USE ldaU_cp,          ONLY: nwfcU, ns, e_hubbard
      USE step_penalty,     ONLY: write_pen

      implicit none

      include 'laxlib.fh'

  integer :: is, isp, ia, m1, m2, err, k
  real(DP), allocatable   :: ftemp1(:), ftemp2(:), f1 (:), vet (:,:)

  real(DP) :: lambda (ldmx), nsum, nsuma

  CALL write_pen (nsp, nspin)

  write (6,'(6(a,i2,a,f8.4,6x))') &
        ('U(',is,') =', Hubbard_U(is) * autoev, is=1,nsp)
      nsum = 0.d0
      allocate( ftemp1(ldmx), ftemp2(ldmx), f1(ldmx*ldmx), vet(ldmx,ldmx) )
      write(6,*) 'nsp',nsp
      do ia = 1, nat
         is=ityp(ia)
         nsuma = 0.d0
         if (Hubbard_U(is).ne.0.d0) then
            do isp = 1, nspin
                do m1 = 1, 2 * Hubbard_l(is) + 1
                   nsuma = nsuma + ns (m1,m1,isp,ia)
                end do
            end do
            if (nspin.eq.1) nsuma = 2.d0 * nsuma
            write(6,'(a,1x,i2,2x,a,f11.7)') 'atom', ia,              &
     &                                      ' Tr[ns(na)]= ',nsuma
            nsum = nsum + nsuma
!
            do isp = 1, nspin

               k = 0
               do m1 = 1, 2 * Hubbard_l(is) + 1
                  do m2 = m1, 2 * Hubbard_l(is) + 1
                     k = k + 1
                     f1 ( k ) = ns (m2,m1,isp,ia)
                  enddo
               enddo

               CALL dspev_drv( 'V', 'L', 2 * Hubbard_l(is) + 1, &
                               f1, lambda, vet, ldmx  )

               write(6,'(a,1x,i2,2x,a,1x,i2)') 'atom', ia, 'spin', isp
               write(6,'(a,7f10.7)') 'eigenvalues: ',(lambda(m1),m1=1,&
     &                                2 * Hubbard_l(is) + 1)
               write(6,*) 'eigenvectors'
               do m2 = 1, 2*Hubbard_l(is)+1
                  write(6,'(i2,2x,7(f10.7,1x))') m2,(real(vet(m1,m2)),&
     &                            m1=1,2 * Hubbard_l(is) + 1)
               end do
               write(6,*) 'occupations'
               do m1 = 1, 2*Hubbard_l(is)+1
                  write (6,'(7(f6.3,1x))') (ns(m1,m2,isp,ia),m2=1,    &
     &                     2*Hubbard_l(is)+1)
               end do
            end do
         end if
      end do
      deallocate ( ftemp1, ftemp2,f1, vet )
      return
      end subroutine write_ns
!
!-------------------------------------------------------------------------
      subroutine dndtau(alpha_a,alpha_s,becwfc,spsi,bp,dbp,wdb,         &
                        offset,wfcU,eigr,proj,ipol,nb_s,nb_e,mykey,dns)
!-----------------------------------------------------------------------
      !! This routine computes the derivative of the ns with respect to the ionic
      !! displacement \(\text{tau}(\text{alpha},\text{ipol})\) used to obtain
      !! the Hubbard contribution to the atomic forces.
      !
      use ions_base, only: nat, nsp, ityp
      use gvecw, only: ngw
      use electrons_base, only: nspin, n => nbsp, nx => nbspx, ispin, f
      USE uspp,           ONLY: nkb
      USE ldaU_cp,        ONLY: Hubbard_U, Hubbard_l, ldmx
      USE ldaU_cp,        ONLY: nwfcU, ns
      USE kinds,          ONLY: DP
      USE mp,             ONLY: mp_sum
      USE mp_pools,       ONLY: intra_pool_comm
!
      implicit none
! input
      integer,      intent(in) :: offset(nat)
      integer,      intent(in) :: alpha_a,alpha_s,ipol
      INTEGER,      INTENT(in) :: nb_s, nb_e, mykey
      COMPLEX(dp),  INTENT(in) :: wfcU(ngw,nwfcU), eigr(ngw,nat)
      REAL(dp),     INTENT(IN) :: becwfc(nkb,nwfcU), bp(nkb,n), &
                                  dbp(nkb,nx,3), wdb(nkb,nwfcU,3)
      real(DP),     intent(in) :: proj(nwfcU,n)
      complex (DP), intent(in) :: spsi(ngw,n)
! output
      real (DP),   intent(out) :: dns(ldmx,ldmx,nspin,nat)
!     dns   derivative of ns(:,:,:,:) w.r.t. tau
!
      integer ibnd,is,i,ia, m1,m2, l, ldim
      real (DP),   allocatable :: dproj(:,:)
!     dproj(nwfcU,n)   derivative of proj(:,:) w.r.t. tau 
!
      CALL start_clock('dndtau')
      !
      allocate (dproj(nwfcU,nb_s:nb_e) )
      call dprojdtau(wfcU,becwfc,spsi,bp,dbp,wdb,eigr,alpha_a,     &
                     alpha_s,ipol,offset(alpha_a),nb_s,nb_e,mykey, &
                     dproj)
      !
      ! compute the derivative of occupation numbers (the quantities dn(m1,m2))
      ! of the atomic orbitals. They are real quantities as well as n(m1,m2)
      !
      dns(:,:,:,:) = 0.d0
      !
      ! band parallelization. If each band appears more than once
      ! compute its contribution only once (i.e. when mykey=0)
      !
      IF ( mykey /= 0 ) GO TO 10
      do ia = 1,nat
         is = ityp(ia) 
         if (Hubbard_U(is).ne.0.d0) then
            ldim = 2*Hubbard_l(is) + 1
            do m1 = 1, ldim
               do m2 = m1, ldim
                  do ibnd = nb_s,nb_e
                     dns(m1,m2,ispin(ibnd),ia) =                    &
     &               dns(m1,m2,ispin(ibnd),ia) +                    &
     &                f(ibnd)*REAL(  proj(offset(ia)+m1,ibnd) *   &
     &                             (dproj(offset(ia)+m2,ibnd))+   &
     &                              dproj(offset(ia)+m1,ibnd) *   &
     &                              (proj(offset(ia)+m2,ibnd)) )
                  end do
                  dns(m2,m1,:,ia) = dns(m1,m2,:,ia)
               end do
            end do
         end if
      end do
!
 10   deallocate (dproj)
      CALL mp_sum(dns, intra_pool_comm)
      CALL stop_clock('dndtau')
      return
      end subroutine dndtau
!
!-----------------------------------------------------------------------
      subroutine dprojdtau(wfcU,becwfc,spsi,bp,dbp,wdb,eigr,alpha_a,    &
                           alpha_s,ipol,offset,nb_s,nb_e,mykey,dproj)
!-----------------------------------------------------------------------
      !! This routine computes the first derivative of the projection
      !! \(\langle\phi^{at}_{I,m1}|S|\psi_{k,v,s}\rangle\) with respect to the
      !! atomic displacement \(u(\alpha,\text{ipol})\) (we remember that 
      !! \(ns_{m1,m2,s,I} = \sum_{k,v} f_{kv} \langle\phi^{at}_{I,m1}|S|\psi_{k,v,s}
      !! \rangle\langle \psi_{k,v,s}|S|\phi^{at}_{I,m2}\rangle\) ).
      !
      use ions_base, only: nat
      use gvecw, only: ngw
      use gvect, only: g, gstart
      use electrons_base, only: n => nbsp, nx => nbspx
      USE uspp,           ONLY: nkb, qq_nt, ofsbeta
      USE ldaU_cp,        ONLY: Hubbard_U, Hubbard_l
      USE ldaU_cp,        ONLY: nwfcU
      use cell_base,      ONLY: tpiba
      USE uspp_param,     only: nh
      use mp_global,      only: intra_bgrp_comm
      use mp,             only: mp_sum
      USE kinds,          ONLY: DP
!
       implicit none
       integer, INTENT(in) :: alpha_a, alpha_s,ipol, offset
! input: the displaced atom
! input: the component of displacement
! input: the offset of the wfcs of the atom "alpha_a,alpha_s"
       INTEGER, INTENT(in) :: nb_s, nb_e, mykey
       complex (DP), intent(in) :: spsi(ngw,n),                     &
     &                   eigr(ngw,nat)
! input: S|evc>, structure factors
       real(DP), intent(in) ::becwfc(nkb,nwfcU), &
     &            bp(nkb,n), dbp(nkb,nx,3), wdb(nkb,nwfcU,3)
       COMPLEX(dp), INTENT(IN) :: wfcU(ngw,nwfcU)
       real(DP), intent(out) :: dproj(nwfcU,nb_s:nb_e)
! output: the derivative of the projection
!
      integer i,ig,m1,ibnd,iwf,ia,is,iv,jv,ldim,alpha,l,m,k,inl
!
      real(dp), allocatable :: dproj0(:,:)
      real(dp) :: gvec
      complex (DP), allocatable :: dwfc(:,:)
      real (DP), allocatable :: betapsi(:,:), dbetapsi(:,:), &
     &                          wfcbeta(:,:),wfcdbeta(:,:), auxwfc(:,:)
!      dwfc(ngw,ldmx),      ! the derivative of the atomic Hubbard wfc
!      betapsi(nh,n),       ! <beta|evc>
!      dbetapsi(nh,n),      ! <dbeta|evc>
!      wfcbeta(nwfcU,nh),   ! <wfc|beta>
!      wfcdbeta(nwfcU,nh),  ! <wfc|dbeta>
 
      ldim = 2 * Hubbard_l(alpha_s) + 1
      dproj(:,:)=0.d0
!
! At first the derivative of the atomic wfc is computed
!
      if (Hubbard_U(alpha_s).ne.0.d0) then
         !
         allocate ( dwfc(ngw,ldim), dproj0(ldim,n) )
         !
         do ig=1,ngw
            gvec = g(ipol,ig)*tpiba 
            do m1=1,ldim
               dwfc(ig,m1) = CMPLX (gvec*AIMAG(wfcU(ig,offset+m1)),      &
     &                             -gvec* DBLE(wfcU(ig,offset+m1)), kind=dp )
            end do
         end do
         !
         ! no need to calculate the G=0 term: it is zero
         !
         CALL dgemm( 'C', 'N', ldim, n, 2*ngw, 2.0_DP, dwfc, 2*ngw, spsi, &
                    2*ngw, 0.0_DP, dproj0, ldim )
         call mp_sum( dproj0, intra_bgrp_comm )
         !
         ! copy to dproj results for the bands treated by this processor
         !
         dproj(offset+1:offset+ldim,:) = dproj0(:,nb_s:nb_e)
         deallocate (dproj0, dwfc)
         !
      end if
      !
      IF( nh(alpha_s) > 0 ) THEN
         !
         allocate (  wfcbeta(nwfcU,nh(alpha_s)) )
         allocate ( wfcdbeta(nwfcU,nh(alpha_s)) )
         allocate (   auxwfc(nwfcU,nh(alpha_s)) )
         !
         do iv=1,nh(alpha_s)
            inl=ofsbeta(alpha_a) + iv
            do m=1,nwfcU
               auxwfc(m,iv) = becwfc(inl,m)
            end do
         end do
         ! following dgemm performs (note that qq is symmetric)
         ! wfcbeta(m,iv) = sum_jv qq(iv,jv,alpha_s)*auxwfc(m,jv)
         CALL dgemm( 'N', 'N', nwfcU, nh(alpha_s), nh(alpha_s), 1.0_DP, &
                  auxwfc, nwfcU, qq_nt(1,1,alpha_s), nh(alpha_s), &
                  0.0_DP, wfcbeta, nwfcU )
         do iv=1,nh(alpha_s)
            inl=ofsbeta(alpha_a) + iv
            do m=1,nwfcU
               auxwfc(m,iv) = wdb(inl,m,ipol)
            end do
         end do
         ! as above with wfcbeta(m,iv) => wfcdbeta
         CALL dgemm( 'N', 'N', nwfcU, nh(alpha_s), nh(alpha_s), 1.0_DP, &
                  auxwfc, nwfcU, qq_nt(1,1,alpha_s), nh(alpha_s), &
                  0.0_DP, wfcdbeta, nwfcU )
         deallocate(auxwfc)
         !
         IF ( mykey == 0 ) THEN
            allocate (  betapsi(nh(alpha_s),nb_s:nb_e) )
            allocate ( dbetapsi(nh(alpha_s),nb_s:nb_e) )
            do iv=1,nh(alpha_s)
               inl=ofsbeta(alpha_a) + iv
               do i=nb_s,nb_e
                  betapsi (iv,i)=bp(inl,i)
                  dbetapsi(iv,i)=dbp(inl,i,ipol)
               end do
            end do
            !
            ! dproj(m,i) = \sum_iv wfcdbeta(m,iv)*betapsi (iv,i) +
            !                      wfcbeta (m,iv)*dbetapsi(iv,i) 
            !
            CALL dgemm( 'N', 'N', nwfcU, nb_e-nb_s+1, nh(alpha_s), 1.0_DP, &
                  wfcdbeta, nwfcU, betapsi(1,nb_s), nh(alpha_s), &
                  1.0_DP, dproj(1,nb_s), nwfcU )
            CALL dgemm( 'N', 'N', nwfcU, nb_e-nb_s+1, nh(alpha_s), 1.0_DP, &
                  wfcbeta, nwfcU, dbetapsi(1,nb_s), nh(alpha_s), &
                  1.0_DP, dproj(1,nb_s), nwfcU )
            !
            deallocate (dbetapsi)
            deallocate (betapsi)
            !
         end if
         ! end band parallelization - only dproj(1,nb_s:nb_e) are calculated
         !
         deallocate (wfcbeta)
         deallocate (wfcdbeta)

      END IF

      return
      end subroutine dprojdtau
!
!-----------------------------------------------------------------------
      SUBROUTINE projwfc_hub( c, nx, eigr, betae, n, nwfcU,  &
     &                        offset, Hubbard_l, wfcU, becwfc, swfc, proj )
!-----------------------------------------------------------------------
      !! Projection on atomic wavefunctions.  
      !! Atomic wavefunctions are not orthogonalized.
      !
      USE kinds,              ONLY: DP
      USE io_global,          ONLY: stdout
      USE mp_global,          ONLY: intra_bgrp_comm
      USE mp,                 ONLY: mp_sum
      USE gvecw,              ONLY: ngw
      USE gvect,              ONLY: gstart
      USE ions_base,          ONLY: nsp, nat
      USE uspp,               ONLY: nkb
      USE cp_interfaces,      only: calbec
!
      IMPLICIT NONE
      INTEGER,     INTENT(IN) :: nx, n, nwfcU, offset(nat), &
                                 Hubbard_l(nsp)
      COMPLEX(DP), INTENT(IN) :: c( ngw, nx ), eigr(ngw,nat)
      COMPLEX(DP), INTENT(INOUT) :: betae(ngw,nkb)
!
      COMPLEX(DP), INTENT(OUT):: wfcU(ngw, nwfcU),    &
     &                           swfc(ngw, nwfcU)
      real(DP), intent(out):: becwfc(nkb,nwfcU), proj(nwfcU,n)
      !
      IF ( nwfcU .EQ. 0 ) RETURN
      !
      CALL start_clock('projwfc_hub')
      !
      ! calculate wfcU = atomic states with associated Hubbard U
      !
      CALL atomic_wfc_hub( offset, Hubbard_l, eigr, nwfcU, wfcU )
      !
      ! calculate bec = <beta|wfc>
      !
      CALL calbec( nwfcU, betae, wfcU, becwfc )
      !
      ! calculate swfc = S|wfc>
      !
      CALL s_wfc( nwfcU, becwfc, betae, wfcU, swfc )
      !
      ! calculate proj = <wfcU|S|c>
      !
      CALL dgemm( 'C', 'N', nwfcU, n, 2*ngw, 2.0_DP, swfc, 2*ngw, c, &
                    2*ngw, 0.0_DP, proj, nwfcU )
      IF ( gstart == 2 ) &
         CALL dger( nwfcU, n,  -1.0_DP, swfc, 2*ngw, c, 2*ngw, proj, nwfcU )
      CALL mp_sum( proj, intra_bgrp_comm )
!
      CALL stop_clock('projwfc_hub')
!
      RETURN
      END SUBROUTINE projwfc_hub
!
!-----------------------------------------------------------------------
      SUBROUTINE atomic_wfc_hub( offset, Hubbard_l, eigr, nwfcU, wfcU )
!-----------------------------------------------------------------------
      !! Compute atomic wavefunctions (not orthogonalized) in G-space.
      !
      USE kinds,              ONLY: DP
      USE gvecw,              ONLY: ngw
      USE gvect,              ONLY: gstart, gg, g
      USE ions_base,          ONLY: nsp, nat, ityp
      USE cell_base,          ONLY: tpiba, omega
      USE atom,               ONLY: rgrid
      USE uspp_param,         ONLY: upf
      USE constants,          ONLY: fpi
!
      IMPLICIT NONE
      INTEGER,     INTENT(in) :: nwfcU, offset(nat), &
                                 Hubbard_l(nsp)
      COMPLEX(DP), INTENT(in) :: eigr( ngw, nat )
      COMPLEX(DP), INTENT(out):: wfcU( ngw, nwfcU )
!
      INTEGER :: natwfc, ndm, is, ia, ir, nb, l, m, lm, i, lmax_wfc
      REAL(DP), ALLOCATABLE ::  ylm(:,:), q(:), jl(:), vchi(:), chiq(:)

      IF( .NOT. ALLOCATED( rgrid ) ) &
         CALL errore( ' atomic_wfc_hub ', ' rgrid not allocated ', 1 )
!
! calculate max angular momentum required in wavefunctions
!
      lmax_wfc=-1
      DO is = 1,nsp
         lmax_wfc = MAX ( lmax_wfc, MAXVAL (upf(is)%lchi(1:upf(is)%nwfc) ) )
      ENDDO
      !
      ALLOCATE(ylm(ngw,(lmax_wfc+1)**2))
      !
      CALL ylmr2 ((lmax_wfc+1)**2, ngw, g, gg, ylm)
      ndm = MAXVAL(rgrid(1:nsp)%mesh)
      !
      ALLOCATE(jl(ndm), vchi(ndm))
      ALLOCATE(q(ngw), chiq(ngw))
!
      DO i=1,ngw
         q(i) = SQRT(gg(i))*tpiba
      END DO
      !
      DO is=1,nsp
         !
         !   radial fourier transform of the chi functions. NOTA BENE:
         !   chi is r times the radial part of the atomic wavefunction
         !
         natwfc=0
         DO nb = 1,upf(is)%nwfc
            l = upf(is)%lchi(nb)
            IF ( l /= Hubbard_l(is) ) GO TO 10 
            DO i=1,ngw
               CALL sph_bes (rgrid(is)%mesh, rgrid(is)%r, q(i), l, jl)
               DO ir=1,rgrid(is)%mesh
                  vchi(ir) = upf(is)%chi(ir,nb)*rgrid(is)%r(ir)*jl(ir)
               ENDDO
               CALL simpson_cp90(rgrid(is)%mesh,vchi,rgrid(is)%rab,chiq(i))
            ENDDO
            !
            !   multiply by angular part and structure factor
            !   NOTA BENE: the factor i^l MUST be present!!!
            !
            DO m = 1,2*l+1
               lm = l**2 + m
               natwfc = natwfc + 1
               DO ia = 1, nat
                  IF( ityp(ia) == is ) THEN
                     wfcU(:,natwfc+offset(ia)) = (0.d0,1.d0)**l * eigr(:,ia) * ylm(:,lm)*chiq(:)
                  END IF
               ENDDO
            ENDDO
  10        CONTINUE
         ENDDO
      ENDDO
!
      IF (natwfc+offset(nat) .NE. nwfcU )  CALL errore('atomic_wfc','unexpected error',natwfc)
!
      do i = 1,nwfcU
        call dscal(2*ngw,fpi/sqrt(omega),wfcU(1,i),1)
      end do
      DEALLOCATE(q, chiq, vchi, jl, ylm)
!
      RETURN
      END SUBROUTINE atomic_wfc_hub
