!
! Copyright (C) 2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
program nmr
!-----------------------------------------------------------------------
!
! Calculation of the NMR chemical shift tensor using pseudopotentials
! and the GIPAW method to reconstruct all-electrons response
!
 
#include "f_defs.h"
 USE global_version,   ONLY : version_number
 use kinds,            only : dp
 USE io_files,         ONLY : nd_nmbr, prefix, outdir, tmp_dir, nwordwfc, &
                              iunwfc
 USE ions_base,        ONLY : ntyp => nsp, ityp
 use parameters,       only : ntypx
 use gvect,            only : g, ngm, ecutwfc, nrxx
 use klist,            only : nks, xk, wk
 use cell_base,        only : tpiba2
 use wavefunctions_module, only: evc
 use wvfct,            only : npwx, nbnd, npw, igk, g2kin, et
 use paw,              only : read_recon
 use uspp,             only : nkb, vkb
 USE uspp_param,       ONLY : nh
 use control_ph,       ONLY : nbnd_occ, alpha_pv
 use nmr_mod,          ONLY : igk_q,  vkb_q, init_nmr
 use phcom,            only : igkq, evq
 use becmod,           only : becp          
!
 implicit none 
 !
 CHARACTER (LEN=9)   :: code = 'NMR'
 character (len=256) :: filerec(ntypx)
 integer :: ios


 integer :: ik, ik_q, iq, isign, iperm0, b0, p0, ibnd, idir
 integer :: iperm1, p1, b1


 complex (DP), allocatable :: dpsi(:,:), phi(:), phi1(:), j_bare(:,:,:) 
 complex (DP), allocatable :: dev (:,:), dev2(:,:), auxg (:)&
      , tmp1(:), tmp2(:)
 complex (DP) :: kkterm_hh(3,3), kkterm_vv(3,3), dchi, &
      chihh(3,3), chivv(3,3)

 real (DP) :: emin, emax, chi_macro, test , b_length

 complex (DP) :: ZDOTC

!
 namelist / inputnmr / prefix, filerec, outdir
 !
 CALL init_clocks( .TRUE. )
 CALL start_clock( code )
 !
 !
 
 call startup (nd_nmbr, code, version_number)

  !
  ! set default value
  !
  prefix = 'pwscf'
  outdir = './'

  read (5, inputnmr, err=200, iostat=ios)
200 call errore('nmr.x', 'reading inputnmr namelist', abs(ios))

  tmp_dir = trim(outdir)




  call read_file 

  call openfil

  call read_recon(filerec)

  allocate ( becp(nkb,nbnd) )
  allocate ( dpsi(npwx,nbnd) )  
  allocate ( evq(npwx,nbnd) )  
  allocate ( dev(npwx,nbnd) )
  allocate ( dev2(npwx,nbnd) )
  allocate ( igkq(npwx))
  allocate (j_bare(nrxx,3,3))
  allocate (tmp1(nrxx),tmp2(nrxx))

!allocate (phi)
!allocate (phi1)

! The algorithm is written for insulators
  nbnd_occ(:)=nbnd



  call init_us_1
  call newd
  call init_nmr
  do ik = 1, nks,7 
     kkterm_hh=(0.d0,0.d0)
     kkterm_vv=(0.d0,0.d0)
     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     g2kin=g2kin*tpiba2

     call davcio (evc, nwordwfc, iunwfc, ik, - 1)
     call init_us_2(npw, igk, xk (1, ik), vkb)


!     call diamagnetic_correction()

     vkb_q=vkb
     igkq=igk
     evq=evc
     call magnetic_kkterm(ik,kkterm_vv,kkterm_hh)

     print *,'kkterm_vv',kkterm_vv
     print *,'kkterm_hh',kkterm_hh


! compute alpha_pv

     emin = et (1, ik)
     do ibnd = 1, nbnd
        emin = min (emin, et (ibnd, ik) )
     enddo
#ifdef __PARA
     ! find the minimum across pools
     call poolextreme (emin, -1)
#endif
     emax = et (1, 1)
     do ibnd = 1, nbnd
        emax = max (emax, et (ibnd, ik) )
     enddo
#ifdef __PARA
     ! find the maximum across pools
     call poolextreme (emax, + 1)
#endif
     alpha_pv = 2.d0 * (emax - emin)
     ! avoid zero value for alpha_pv
     alpha_pv = max (alpha_pv, 1.0d-2)
  
     print *,'alpha_pv',alpha_pv   

     do iq=1,3            !loop over the 3 cartesion direction
        do isign=0,3,3    !the plus sign are before the minus sign, so
                          !to have minus kpoint, you should add 3 to ik

        ! get the wavefunction u_kq
           ik_q=ik+iq+isign

           call gk_sort (xk (1, ik_q), ngm, g, ecutwfc / tpiba2, npw,&
                igkq, g2kin)

           g2kin=g2kin*tpiba2
           call davcio (evq, nwordwfc, iunwfc, ik_q, - 1)
           call init_us_2(npw, igkq, xk (1, ik_q), vkb_q)



! Loop over the two directions of the applied field perpendicular
! to the q direction
           do iperm0=1,-1,-2

              b0= mod(iq+iperm0+2,3) +1 !give b0 different direction from iq
              p0= mod(iq-iperm0+2,3) +1 !give p0 third orthogonal direction

              dev=(0.d0,0.d0)
              call take_nloc_k_kq(ik, ik_q, evc, p0, dev )


              call solve_cg(dev,ik, dpsi)


           
           !! call paramagnetic_correction()
              
!              do ibnd =1, nbnd_occ(ik)
                 
              do idir = 1,3 

                 dev=(0.d0,0.d0)
                 call grad(ik_q, dpsi, idir, dev )

!                    tmp1=0.d0
!                    tmp2=0.d0
!                    tmp1(:)=evc(:,ibnd)
!                    tmp2(:)=dev(:,ibnd)

                 call add_j_bare (evc, dev, &
                      sign(DBLE(iperm0),DBLE(1-isign)) *wk(ik),&
                      j_bare(1,b0,idir))

              enddo
              do idir = 1,3 

                 dev=(0.d0,0.d0)
                 call grad(ik, evc, idir, dev )

                 call add_j_bare (dpsi, dev, &
                      - sign(DBLE(iperm0),DBLE(1-isign)) &  
                      *wk(ik), j_bare(1,b0,idir))   
              enddo
              
              do iperm1=1,-1,-2
                 
                 b1= mod(iq+iperm1+2,3) +1 
                 p1= mod(iq-iperm1+2,3) +1
                 
                 dev=(0.d0,0.d0)
                 call grad(ik, dpsi, p1, dev )
                 
                 dev2=(0.d0,0.d0)
                 call grad(ik_q, dpsi,  p1, dev2 )

                 dev = 0.5d0*(dev+dev2)

                 dchi=(0.d0,0.d0)
                 do ibnd=1,nbnd
                    print *, ZDOTC(npw,evc(1,ibnd),1,dev(1,ibnd),1)
                    dchi = dchi+  ZDOTC(npw,evc(1,ibnd),1,dev(1,ibnd),1)
                 enddo

                    dchi = DBLE(iperm1)*DBLE(iperm0)*wk(ik)* &
                         (dchi-kkterm_hh(p0,p1)*DBLE(nbnd))
                 if (b0 .eq. b1) then
                    chihh(b0,b1) = chihh(b0,b1) + dchi 
                 else
                    chihh(b0,b1) = chihh(b0,b1) + 2.d0* dchi
                 endif
                 
                 print *,'chihh',chihh

                 dev=(0.d0,0.d0)
                 call take_nloc_k_kq(ik_q, ik, dpsi, p1, dev )
                
                 dchi=(0.d0,0.d0)
                 do ibnd=1,nbnd
                    dchi = dchi+  ZDOTC(npw,evc(1,ibnd),1,dev(1,ibnd),1)
                 enddo
                 print *,'dchi',dchi
                    dchi = DBLE(iperm1)*DBLE(iperm0)*wk(ik)* &
                         (dchi-kkterm_vv(p0,p1)*DBLE(nbnd))
                 if (b0 .eq. b1) then
                    chivv(b0,b1) = chivv(b0,b1) + dchi 
                 else
                    chihh(b0,b1) = chivv(b0,b1) + 2.d0* dchi
                 endif
                 
                 print *,'chivv',chivv 
              enddo
              
           enddo
        enddo
     enddo
  enddo
  b_length =  0.001d0
  chihh = -2.d0*DBLE(chihh(:,:))/(b_length**2)
  chihh = chihh *(.529177d-8)**3*0.6022d24/(137.036**2)*1.d6

  chi_macro= DBLE(chihh(1,1)+chihh(2,2)+chihh(3,3))/3.d0

  print *,'chihh_macro',chi_macro
  

  chivv = -2.d0*DBLE(chivv(:,:))/(b_length**2)
  chivv = chivv *(.529177d-8)**3*0.6022d24/(137.036**2)*1.d6

  chi_macro= DBLE(chivv(1,1)+chivv(2,2)+chivv(3,3))/3.d0

  print *,'chivv_macro',chi_macro


end program nmr
