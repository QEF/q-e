!
! Copyright (C) 2004-2009 Dario Alfe' and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE write_casino_pwfn(gather)

   USE kinds, ONLY: DP
   USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
   USE cell_base, ONLY: omega, alat, tpiba2, at, bg
   USE printout_base, ONLY: title    ! title of the run
   USE constants, ONLY: tpi, e2
   USE ener, ONLY: ewld, ehart, etxc, vtxc, etot, etxcc, demet, ef
   USE gvect, ONLY: ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                    nrxx, g, gg, ecutwfc, gcutm, nl, igtongl
   USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss
   USE lsda_mod, ONLY: lsda, nspin
   USE scf, ONLY: rho, rho_core, rhog_core, vnew
   USE ldaU, ONLY : eth
   USE vlocal, ONLY: vloc, strf
   USE wvfct, ONLY: npw, npwx, nbnd, igk, g2kin, wg, et
   USE control_flags, ONLY : gamma_only
   USE uspp, ONLY: nkb, vkb, dvan
   USE uspp_param, ONLY: nh
   USE io_global, ONLY: stdout, ionode, ionode_id
   USE io_files, ONLY: nd_nmbr, nwordwfc, iunwfc
   USE wavefunctions_module, ONLY : evc
   USE funct, ONLY : dft_is_meta
   USE mp_global, ONLY: inter_pool_comm, intra_pool_comm, nproc_pool
   USE mp, ONLY: mp_sum, mp_gather, mp_bcast
   USE dfunct, ONLY : newd

   IMPLICIT NONE
   LOGICAL, INTENT(in) :: gather
   INTEGER :: ig, ibnd, ik, io, ispin, nbndup, nbnddown, &
              nk, ngtot, ig7, ikk, id, ip
   INTEGER, ALLOCATABLE :: idx(:), igtog(:)
   LOGICAL :: exst
   REAL(DP) :: ek, eloc, enl
   INTEGER, EXTERNAL :: atomic_number
   REAL (DP), EXTERNAL :: ewald, w1gauss

   INTEGER ngtot_g
   INTEGER, ALLOCATABLE :: ngtot_d(:), ngtot_cumsum(:), indx(:)
   REAL(DP), ALLOCATABLE :: g_l(:,:), g_g(:,:), g2(:)
   COMPLEX(DP), ALLOCATABLE :: evc_l(:), evc_g(:)

   call init_us_1
   call newd

   ! four times npwx should be enough
   allocate (idx (4*npwx) )
   allocate (igtog (4*npwx) )
   idx(:) = 0
   igtog(:) = 0

   if( lsda )then
      nbndup = nbnd
      nbnddown = nbnd
      nk = nks/2
      !     nspin = 2
   else
      nbndup = nbnd
      nbnddown = 0
      nk = nks
      !     nspin = 1
   endif

   call calc_energies

   do ispin = 1, nspin
      do ik = 1, nk
         ikk = ik + nk*(ispin-1)
         call gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
         do ig =1, npw
            if( igk(ig) > 4*npwx ) &
               call errore ('pw2casino','increase allocation of index', ig)
            idx( igk(ig) ) = 1
         enddo
      enddo
   enddo

   ngtot = 0
   do ig = 1, 4*npwx
      if( idx(ig) >= 1 )then
         ngtot = ngtot + 1
         igtog(ngtot) = ig
      endif
   enddo

   if(ionode.or..not.gather)then
      io = 77
      write (6,'(/,5x,''Writing file pwfn.data for program CASINO'')')
      call seqopn( 77, 'pwfn.data', 'formatted',exst)
      call write_header
   endif

   allocate ( g_l(3,ngtot), evc_l(ngtot) )
   do ig = 1, ngtot
      g_l(:,ig) = g(:,igtog(ig))
   enddo

   if(gather)then
      allocate ( ngtot_d(nproc_pool), ngtot_cumsum(nproc_pool) )
      call mp_gather( ngtot, ngtot_d, ionode_id, intra_pool_comm )
      call mp_bcast( ngtot_d, ionode_id, intra_pool_comm )
      id = 0
      do ip = 1,nproc_pool
         ngtot_cumsum(ip) = id
         id = id + ngtot_d(ip)
      enddo
      ngtot_g = id

      allocate ( g_g(3,ngtot_g), evc_g(ngtot_g) )
      call mp_gather( g_l, g_g, ngtot_d, ngtot_cumsum, ionode_id, intra_pool_comm)

      if(ionode)then
         allocate ( indx(ngtot_g) )
         call create_index2(g_g,indx)
         call write_gvecs(g_g,indx)
      endif
   else
      allocate ( indx(ngtot) )
      call create_index2(g_l,indx)
      call write_gvecs(g_l,indx)
   endif

   if(ionode.or..not.gather)call write_wfn_head

   do ik = 1, nk
      if(ionode.or..not.gather)call write_kpt_head

      do ispin = 1, nspin
         ikk = ik + nk*(ispin-1)
         if( nks > 1 )then
            call gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
            call davcio(evc,nwordwfc,iunwfc,ikk,-1)
         endif
         do ibnd = 1, nbnd
            if(ionode.or..not.gather)call write_bnd_head
            evc_l(:) = (0.d0, 0d0)
            do ig=1, ngtot
               ! now for all G vectors find the PW coefficient for this k-point
               find_ig: do ig7 = 1, npw
                  if( igk(ig7) == igtog(ig) )then
                     evc_l(ig) = evc(ig7,ibnd)
                     exit find_ig
                  endif
               enddo find_ig
            enddo
            if(gather)then
               call mp_gather( evc_l, evc_g, ngtot_d, ngtot_cumsum, ionode_id, intra_pool_comm)

               if(ionode)call write_wfn_data(evc_g,indx)
            else
               call write_wfn_data(evc_l,indx)
            endif
         enddo
      enddo
   enddo
   if(ionode.or..not.gather)close(io)

   deallocate (igtog)
   deallocate (idx)

   deallocate ( g_l, evc_l )
   if(gather) deallocate ( ngtot_d, ngtot_cumsum, g_g, evc_g )
   if(ionode.or..not.gather) deallocate (indx)

CONTAINS

   SUBROUTINE calc_energies
      USE becmod, ONLY: becp, calbec, allocate_bec_type, deallocate_bec_type

      COMPLEX(DP), ALLOCATABLE :: aux(:)
      INTEGER :: ibnd, j, ig, ik, ikk, ispin, na, nt, ijkb0, ikb, ih, jh, jkb

      REAL(DP) :: charge, etotefield

      allocate (aux(nrxx))
      call allocate_bec_type ( nkb, nbnd, becp )

      ek  = 0.d0
      eloc= 0.d0
      enl = 0.d0
      demet=0.d0
      !
      do ispin = 1, nspin
         !
         !     calculate the local contribution to the total energy
         !
         !      bring rho to G-space
         !
         aux(:) = cmplx( rho%of_r(:,ispin), 0.d0,kind=DP)
         call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
         !
         do nt=1,ntyp
            do ig = 1, ngm
               eloc = eloc + vloc(igtongl(ig),nt) * strf(ig,nt) &
                    * conjg(aux(nl(ig)))
            enddo
         enddo

         do ik = 1, nk
            ikk = ik + nk*(ispin-1)
            call gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
            call davcio (evc, nwordwfc, iunwfc, ikk, - 1)
            call init_us_2 (npw, igk, xk (1, ikk), vkb)
            call calbec ( npw, vkb, evc, becp )
            !
            ! -TS term for metals (ifany)
            !
            if( degauss > 0.0_dp)then
               do ibnd = 1, nbnd
                  demet = demet + wk (ik) * &
                     degauss * w1gauss ( (ef-et(ibnd,ik)) / degauss, ngauss)
               enddo
            endif
            !
            ! calculate the kinetic energy
            !
            do ibnd = 1, nbnd
               do j = 1, npw
                  ek = ek +  conjg(evc(j,ibnd)) * evc(j,ibnd) * &
                                  g2kin(j) * wg(ibnd,ikk)
               enddo

               !
               ! Calculate Non-local energy
               !
               ijkb0 = 0
               do nt = 1, ntyp
                  do na = 1, nat
                     if(ityp (na) .eq. nt)then
                        do ih = 1, nh (nt)
                           ikb = ijkb0 + ih
                           enl=enl+conjg(becp%k(ikb,ibnd))*becp%k(ikb,ibnd) &
                                *wg(ibnd,ikk)* dvan(ih,ih,nt)
                           do jh = ( ih + 1 ), nh(nt)
                              jkb = ijkb0 + jh
                              enl=enl + &
                                   (conjg(becp%k(ikb,ibnd))*becp%k(jkb,ibnd)+&
                                    conjg(becp%k(jkb,ibnd))*becp%k(ikb,ibnd))&
                                   * wg(ibnd,ikk) * dvan(ih,jh,nt)

                           enddo

                        enddo
                        ijkb0 = ijkb0 + nh (nt)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo

#ifdef __PARA
      call mp_sum( eloc,  intra_pool_comm )
      call mp_sum( ek,    intra_pool_comm )
      call mp_sum( ek,    inter_pool_comm )
      call mp_sum( enl,   inter_pool_comm )
      call mp_sum( demet, inter_pool_comm )
#endif
      eloc = eloc * omega
      ek = ek * tpiba2

      !
      ! compute ewald contribution
      !
      ewld = ewald( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, &
           g, gg, ngm, gcutm, gstart, gamma_only, strf )
      !
      ! compute hartree and xc contribution
      !
      call v_of_rho( rho, rho_core, rhog_core, &
                     ehart, etxc, vtxc, eth, etotefield, charge, vnew )
      !
      etot=(ek + (etxc-etxcc)+ehart+eloc+enl+ewld)+demet

      call deallocate_bec_type (becp)
      deallocate (aux)

      write (stdout,*) 'Kinetic energy   ', ek/e2
      write (stdout,*) 'Local energy     ', eloc/e2
      write (stdout,*) 'Non-Local energy ', enl/e2
      write (stdout,*) 'Ewald energy     ', ewld/e2
      write (stdout,*) 'xc contribution  ',(etxc-etxcc)/e2
      write (stdout,*) 'hartree energy   ', ehart/e2
      if( degauss > 0.0_dp ) &
         write (stdout,*) 'Smearing (-TS)   ', demet/e2
      write (stdout,*) 'Total energy     ', etot/e2


   END SUBROUTINE calc_energies

   SUBROUTINE write_header
      INTEGER j, na, nt, at_num

      write(io,'(a)') title
      write(io,'(a)')
      write(io,'(a)') ' BASIC INFO'
      write(io,'(a)') ' ----------'
      write(io,'(a)') ' Generated by:'
      write(io,'(a)') ' PWSCF'
      write(io,'(a)') ' Method:'
      write(io,'(a)') ' DFT'
      write(io,'(a)') ' DFT Functional:'
      write(io,'(a)') ' unknown'
      write(io,'(a)') ' Pseudopotential'
      write(io,'(a)') ' unknown'
      write(io,'(a)') ' Plane wave cutoff (au)'
      write(io,*) ecutwfc/2
      write(io,'(a)') ' Spin polarized:'
      write(io,*)lsda
      if( degauss > 0.0_dp )then
         write(io,'(a)') ' Total energy (au per primitive cell; includes -TS term)'
         write(io,*)etot/e2, demet/e2
      else
         write(io,'(a)') ' Total energy (au per primitive cell)'
         write(io,*)etot/e2
      endif
      write(io,'(a)') ' Kinetic energy (au per primitive cell)'
      write(io,*)ek/e2
      write(io,'(a)') ' Local potential energy (au per primitive cell)'
      write(io,*)eloc/e2
      write(io,'(a)') ' Non local potential energy(au per primitive cel)'
      write(io,*)enl/e2
      write(io,'(a)') ' Electron electron energy (au per primitive cell)'
      write(io,*)ehart/e2
      write(io,'(a)') ' Ion ion energy (au per primitive cell)'
      write(io,*)ewld/e2
      write(io,'(a)') ' Number of electrons per primitive cell'
      write(io,*)nint(nelec)
      ! uncomment the following ifyou want the Fermi energy - KN 2/4/09
      !  write(io,'(a)') ' Fermi energy (au)'
      !  write(io,*) ef/e2
      write(io,'(a)') ' '
      write(io,'(a)') ' GEOMETRY'
      write(io,'(a)') ' -------- '
      write(io,'(a)') ' Number of atoms per primitive cell '
      write(io,*) nat
      write(io,'(a)')' Atomic number and position of the atoms(au) '
      do na = 1, nat
         nt = ityp(na)
         at_num = atomic_number(trim(atm(nt)))
         write(io,'(i6,3f20.14)') at_num, (alat*tau(j,na),j=1,3)
      enddo
      write(io,'(a)') ' Primitive lattice vectors (au) '
      write(io,100) alat*at(1,1), alat*at(2,1), alat*at(3,1)
      write(io,100) alat*at(1,2), alat*at(2,2), alat*at(3,2)
      write(io,100) alat*at(1,3), alat*at(2,3), alat*at(3,3)
      write(io,'(a)') ' '

  100 format (3(1x,f20.15))

   END SUBROUTINE write_header


   SUBROUTINE write_gvecs(g,indx)
      REAL(DP),INTENT(in) :: g(:,:)
      INTEGER,INTENT(in) :: indx(:)
      INTEGER ig

      write(io,'(a)') ' G VECTORS'
      write(io,'(a)') ' ---------'
      write(io,'(a)') ' Number of G-vectors'
      write(io,*) size(g,2)
      write(io,'(a)') ' Gx Gy Gz (au)'
      do ig = 1, size(g,2)
         write(io,100) tpi/alat*g(1,indx(ig)), tpi/alat*g(2,indx(ig)), &
               tpi/alat*g(3,indx(ig))
      enddo

  100 format (3(1x,f20.15))

      write(io,'(a)') ' '
   END SUBROUTINE write_gvecs


   SUBROUTINE write_wfn_head
      write(io,'(a)') ' WAVE FUNCTION'
      write(io,'(a)') ' -------------'
      write(io,'(a)') ' Number of k-points'
      write(io,*) nk
   END SUBROUTINE write_wfn_head


   SUBROUTINE write_kpt_head
      INTEGER j

      write(io,'(a)') ' k-point # ; # of bands (up spin/down spin); &
            &           k-point coords (au)'
      write(io,'(3i4,3f20.16)') ik, nbndup, nbnddown, &
            (tpi/alat*xk(j,ik),j=1,3)
   END SUBROUTINE write_kpt_head


   SUBROUTINE write_bnd_head
      ! KN: if you want to print occupancies, replace these two lines ...
      write(io,'(a)') ' Band, spin, eigenvalue (au)'
      write(io,*) ibnd, ispin, et(ibnd,ikk)/e2
      ! ...with the following two - KN 2/4/09
      ! write(io,'(a)') ' Band, spin, eigenvalue (au), occupation number'
      ! write(io,*) ibnd, ispin, et(ibnd,ikk)/e2, wg(ibnd,ikk)/wk(ikk)
      write(io,'(a)') ' Eigenvectors coefficients'
   END SUBROUTINE write_bnd_head


   SUBROUTINE write_wfn_data(evc,indx)
      COMPLEX(DP),INTENT(in) :: evc(:)
      INTEGER,INTENT(in) :: indx(:)
      INTEGER ig

      do ig=1, size(evc,1)
         write(io,*)evc(indx(ig))
      enddo
   END SUBROUTINE write_wfn_data


   SUBROUTINE create_index2(y,x_index)
      DOUBLE PRECISION,INTENT(in) :: y(:,:)
      INTEGER,INTENT(out) :: x_index(size(y,2))
      DOUBLE PRECISION y2(size(y,2))
      INTEGER i
      do i = 1,size(y,2)
         y2(i) = sum(y(:,i)**2)
      enddo
      call create_index(y2,x_index)
   END SUBROUTINE create_index2


   SUBROUTINE create_index(y,x_index)
 !-----------------------------------------------------------------------------!
 ! This subroutine creates an index array x_index for the n items of data in   !
 ! the array y.  Adapted from Numerical Recipes.                               !
 ! Copied from merge_pwfn.f90, included with CASINO distribution               !
 !-----------------------------------------------------------------------------!
      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(in) :: y(:)
      INTEGER,INTENT(out) :: x_index(:)
      INTEGER,PARAMETER :: ins_sort_thresh=7,stacksize=80
      INTEGER n,i,x_indexj,ir,itemp,j,jstack,k,l,lp1,istack(stacksize)
      DOUBLE PRECISION yj
      n=size(x_index)
      do j=1,n
         x_index(j)=j
      enddo ! j
      if(n<=1)return
      jstack=0
      l=1
      ir=n
      do
         if(ir-l<ins_sort_thresh)then
    jloop : do j=l+1,ir
               x_indexj=x_index(j) ; yj=y(x_indexj)
               do i=j-1,l,-1
                  if(y(x_index(i))<=yj)then
                     x_index(i+1)=x_indexj
                     cycle jloop
                  endif! y(x_index(i))<=yj
                  x_index(i+1)=x_index(i)
               enddo ! i
               x_index(l)=x_indexj
            enddo jloop ! j
            if(jstack==0)return
            ir=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
         else
            k=(l+ir)/2
            lp1=l+1
            itemp=x_index(k)    ; x_index(k)=x_index(lp1)  ; x_index(lp1)=itemp
            if(y(x_index(l))>y(x_index(ir)))then
               itemp=x_index(l)   ; x_index(l)=x_index(ir)   ; x_index(ir)=itemp
            endif
            if(y(x_index(lp1))>y(x_index(ir)))then
               itemp=x_index(lp1) ; x_index(lp1)=x_index(ir) ; x_index(ir)=itemp
            endif
            if(y(x_index(l))>y(x_index(lp1)))then
               itemp=x_index(l)   ; x_index(l)=x_index(lp1)  ; x_index(lp1)=itemp
            endif
            i=lp1
            j=ir
            x_indexj=x_index(lp1)
            yj=y(x_indexj)
            do
               do
                  i=i+1
                  if(y(x_index(i))>=yj)exit
               enddo ! i
               do
                  j=j-1
                  if(y(x_index(j))<=yj)exit
               enddo ! j
               if(j<i)exit
               itemp=x_index(i) ; x_index(i)=x_index(j) ; x_index(j)=itemp
            enddo
            x_index(lp1)=x_index(j)
            x_index(j)=x_indexj
            jstack=jstack+2
            if(jstack>stacksize)then
               write(6,*)'stacksize is too small.'
               stop
            endif! jstack>stacksize
            if(ir-i+1>=j-l)then
               istack(jstack)=ir
               istack(jstack-1)=i
               ir=j-1
            else
               istack(jstack)=j-1
               istack(jstack-1)=l
               l=i
            endif! ir-i+1>=j-l
         endif! ir-l<ins_sort_thresh
      enddo
   END SUBROUTINE create_index

END SUBROUTINE write_casino_pwfn
