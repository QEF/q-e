!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------



module exx
#include "f_defs.h"

  USE klist,    ONLY :  nks 
  USE kinds,    ONLY : DP
  implicit none
  integer :: currentk
  real (kind=DP):: exxalfa=1.d0 ! 1 if exx, 0 elsewhere
  logical:: exxstart=.false. !1 if initialited
  integer :: iunexx=49
  CHARACTER(len=80) :: exx_file = 'os.exx'
  integer :: exx_nwordwfc
  real (kind=DP) :: exxdivergency=0d0
  real (kind=DP) :: exxdivfac=0d0

contains

  subroutine exxinit()

    !This subroutine is run before the first H_psi()
    !It saves the wavefunctions for the right density matrix. in real space
    !It saves all the wavefunctions in an only file called prefix.exx
    USE wavefunctions_module, ONLY : evc  
    USE io_files,             ONLY : nwordwfc
    USE io_files,             ONLY : prefix
    USE io_files,             ONLY : tmp_dir, iunwfc, iunigk
    USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
    USE gsmooth,              ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
         nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
    USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et


    integer :: ios, ik,ibnd
    COMPLEX(KIND=DP) :: temppsic(nrxx)   
    complex(kind=dp) :: tempevc( npwx, nbnd )

    call start_clock ('exxinit')
    exx_nwordwfc=2*nrxx
    exx_file = TRIM( tmp_dir ) // TRIM( prefix ) //'.exx'

    if (.not.exxstart) then 
       open (unit=iunexx, file=exx_file, form="unformatted", &
             action="readwrite", status="replace", &
             recl=DIRECT_IO_FACTOR*exx_nwordwfc, access="direct")
    endif

    IF ( nks > 1 ) REWIND( iunigk )
    DO ik = 1, nks
       call davcio (tempevc, nwordwfc, iunwfc, ik, -1)
       IF ( nks > 1 ) READ( iunigk ) npw, igk

       do ibnd =1, nbnd     
          temppsic(:) = ( 0.D0, 0.D0 )
          temppsic(nls(igk(1:npw))) = tempevc(1:npw,ibnd)
          CALL cft3s( temppsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
          CALL davcio( temppsic, exx_nwordwfc, iunexx, (ik-1)*nbnd+ibnd, 1 )
       end do
    end do
    exxstart=.true.
    
    call facdiv ()


    call stop_clock ('exxinit')  

  end subroutine exxinit

  subroutine vexx(lda, n, m, psi, hpsi)
    !This routine calculates V_xx \Psi

    ! ... This routine computes the product of the Hamiltonian
    ! ... matrix with m wavefunctions contained in psi
    !
    ! ... input:
    ! ...    lda   leading dimension of arrays psi, spsi, hpsi
    ! ...    n     true dimension of psi, spsi, hpsi
    ! ...    m     number of states psi
    ! ...    psi
    !
    ! ... output:
    ! ...    hpsi  Vexx*psi
    !
    USE cell_base,  ONLY : alat, omega
    USE gvect,      ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm, gstart
    USE gsmooth,    ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
         nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
    USE wvfct,      ONLY : nbnd, npwx, npw, igk, wg, et
    USE klist,      ONLY : xk,wk
    USE lsda_mod,   ONLY : lsda, current_spin, isk


    USE gvect,      ONLY : g, nl

    INTEGER          :: lda, n, m, kpsi
    COMPLEX(KIND=DP) :: psi(lda,m) 
    COMPLEX(KIND=DP) :: hpsi(lda,m)
    COMPLEX(KIND=DP) :: tempphic(nrxx)   
    COMPLEX(KIND=DP) :: temppsic(nrxx)   
    COMPLEX(KIND=DP) :: result(nrxx)
    COMPLEX(KIND=DP) :: rhoc(nrxx)   
    COMPLEX(KIND=DP) :: vc(nrxx)   
    real (kind=DP)   :: fac(ngm)
    integer          :: ibnd, ik, im , ig, inrxx, ikdiv
    real(kind=DP), parameter  :: pi =  3.14159265358979d0
    real(kind=DP), parameter  :: fpi = 4.d0 * 3.14159265358979d0, &
         e2  = 2.d0
    real(kind=DP) :: tpiba2, qq
    tpiba2 = (fpi / 2.d0 / alat) **2

    call start_clock ('exx')
    ! write (*,*) exx_nwordwfc,lda,n,m, lda*n

    do im=1,m !for each band of psi (the k cycle is outside band)
       temppsic(:) = ( 0.D0, 0.D0 )
       temppsic(nls(igk(1:npw))) = psi(1:npw,im)
       CALL cft3s( temppsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )

       result(:)= ( 0.D0, 0.D0 )
       do ik=1,nks !for each k point of phi
          if (ik /= currentk) cycle
          if (lsda) then
             if ( isk(ik) /= current_spin ) cycle
          end if
             do ig=1,ngm
                qq = (xk(1,ik)-xk(1,currentk)+g(1,ig))**2 + &
                     (xk(2,ik)-xk(2,currentk)+g(2,ig))**2 + &
                     (xk(3,ik)-xk(3,currentk)+g(3,ig))**2
                if (qq.gt.1.d-8) then
                   fac(ig)=e2*fpi/tpiba2/qq
                else
                   fac(ig)=e2*fpi/tpiba2*exxdivfac
                end if
             end do
             do ibnd=1,nbnd !for each band of phi
                if ( abs(wg(ibnd,ik)) < 1.d-6) cycle
                !
                !loads the phi from file
                !
                CALL davcio(tempphic,exx_nwordwfc,iunexx,(ik-1)*nbnd+ibnd,-1)
                !calculate rho in real space
                rhoc(:)=conjg(tempphic(:))*temppsic(:) / omega
                !brings it to G-space
                CALL cft3s( rhoc,nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1 )

                vc(:) = ( 0.D0, 0.D0 )
                do ig=1,ngm
                   vc(nls(ig)) = fac(ig) * rhoc(nls(ig))! v in G-space
                end do
                vc = vc * wg (ibnd, ik) 

                !brings back v in real space
                CALL cft3s( vc, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 1 ) 

                do inrxx=1,nrxx
                   !accumulates over bands and k points
                   result(inrxx)=result(inrxx)+vc(inrxx)*tempphic(inrxx)
                end do
             end do
          end do
          !brings back result in G-space
          CALL cft3s( result, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2 )
          !adds it to hpsi
          hpsi(1:npw,im)=hpsi(1:npw,im) - exxalfa*result(nls(igk(1:npw)))
       end do
       call stop_clock ('exx')

     end subroutine vexx

  function exxenergy ()
    ! This function is called to correct the deband value and have the correct energy 
    USE io_files,             ONLY : nwordwfc

    USE cell_base,  ONLY : alat, omega
    USE gvect,      ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm, gstart
    USE gsmooth,    ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
         nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
    USE wvfct,      ONLY : nbnd, npwx, npw, igk, wg, et
    USE klist,      ONLY : xk,wk
    USE lsda_mod,   ONLY : lsda, current_spin, isk

    USE io_files,             ONLY : tmp_dir, iunwfc, iunigk
    USE parameters, only : ndmx
    USE gvect,      ONLY : g, nl
    implicit none
    integer          :: m
    REAL (KIND=DP)   :: exxenergy
    COMPLEX(KIND=DP) :: psi(ndmx,nbnd) 
    COMPLEX(KIND=DP) :: hpsi(ndmx,nbnd)
    COMPLEX(KIND=DP) :: tempphic(nrxx)   
    COMPLEX(KIND=DP) :: temppsic(nrxx)   
    COMPLEX(KIND=DP) :: result(nrxx)
    COMPLEX(KIND=DP) :: rhoc(nrxx)   
    COMPLEX(KIND=DP) :: vc(nrxx)   
    real (kind=DP)   :: fac(ngm)
    integer          :: ibnd, ik, im , ig, inrxx, ikk
    real(kind=DP), parameter  :: fpi = 4.d0 * 3.14159265358979d0, &
         e2  = 2.d0
    real(kind=DP) :: tpiba2, qq
    complex(kind=dp) :: tempevc( npwx, nbnd )

    tpiba2 = (fpi / 2.d0 / alat) **2

    call start_clock ('exxenergy')
    exxenergy=0
    IF ( nks > 1 ) REWIND( iunigk )

    do ikk=1,nks
       call davcio (tempevc, nwordwfc, iunwfc, ikk, -1)
       IF ( nks > 1 ) READ( iunigk ) npw, igk

       m=nbnd ! We have to compute all bands to have the correct energy
       do im=1,m !for each band of psi (the k cycle )

          result(:)= ( 0.D0, 0.D0 )
          temppsic(:) = ( 0.D0, 0.D0 )
          temppsic(nls(igk(1:npw))) = tempevc(1:npw,im)

          CALL cft3s( temppsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
          do ik=1,nks !for each k point of phi
             if (ik /=ikk ) cycle
             if (lsda) then
                if ( isk(ik) /= isk(ikk) ) cycle
             end if
             do ig=1,ngm
                qq = (xk(1,ik)-xk(1,currentk)+g(1,ig))**2 + &
                     (xk(2,ik)-xk(2,currentk)+g(2,ig))**2 + &
                     (xk(3,ik)-xk(3,currentk)+g(3,ig))**2
                if (qq.gt.1.d-8) then
                   fac(ig)=e2*fpi/tpiba2/qq
                else
                   fac(ig)=exxdivergency
                end if
             end do
             do ibnd=1,nbnd !for each band of psi
                if ( abs(wg(ibnd,ik)) < 1.d-6) cycle
                !
                !loads the phi from file
                !
                CALL davcio(tempphic,exx_nwordwfc,iunexx,(ik-1)*nbnd+ibnd,-1)
                !calculate rho in real space
                rhoc(:)=conjg(tempphic(:))*temppsic(:) / omega
                !brings it to G-space
                CALL cft3s( rhoc,nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1 )

                vc(:) = ( 0.D0, 0.D0 )
                do ig=1,ngm
                   vc(nls(ig)) = fac(ig) * rhoc(nls(ig))! v in G-space
                end do
                vc = vc * wg (ibnd, ik) 

                !brings back v in real space
                CALL cft3s( vc, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 1 ) 

                do inrxx=1,nrxx
                   !accumulates over bands and k points
                   result(inrxx)=result(inrxx)+vc(inrxx)*tempphic(inrxx)
                end do
             end do
          end do
          !brings back result in G-space
          CALL cft3s( result, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2 )

          !hpsi(1:npw,im)=hpsi(1:npw,im) - exxalfa*result(nls(igk(1:npw)))
          do ig=1,npw
             exxenergy=exxenergy - exxalfa*real(conjg ( tempevc(ig,im)) *result(nls(igk(ig))))
          end do

       end do
    end do
    call stop_clock ('exxenergy')
  end function exxenergy

  function exxenergy2 ()
    ! This function is called to correct the deband value and have the correct energy 
    USE io_files,   ONLY : iunigk,iunwfc, nwordwfc
    USE wvfct,      ONLY : nbnd, npwx, npw, igk, wg
    USE wavefunctions_module, ONLY : evc
    USE lsda_mod,   ONLY : lsda, current_spin, isk

    implicit none
    REAL (KIND=DP)   :: exxenergy2,  energy
    INTEGER          :: ibnd, ik
    COMPLEX(KIND=DP) :: vxpsi ( npwx, nbnd ), psi(npwx,nbnd)
    COMPLEX(KIND=DP) :: ZDOTC

    call start_clock ('exxenergy')

    energy=0.d0
    IF ( nks > 1 ) REWIND( iunigk )
    do ik=1,nks
       currentk = ik
       IF ( nks > 1 ) READ( iunigk ) npw, igk
       IF ( lsda ) current_spin = isk(ik)
       call davcio (psi, nwordwfc, iunwfc, ik, -1)
       vxpsi(:,:) = (0.d0, 0.d0)
!!    subroutine vexx(lda, n, m, psi, hpsi)
       call vexx(npwx,npw,nbnd,psi,vxpsi)
       do ibnd=1,nbnd
          energy = energy + &
                   wg(ibnd,ik) * ZDOTC(npw,psi(1,ibnd),1,vxpsi(1,ibnd),1)
       end do
    end do
    exxenergy2 = energy

    call stop_clock ('exxenergy')
  end function exxenergy2

  subroutine facdiv ()
    USE cell_base,  ONLY : alat, omega
    USE gvect,      ONLY : ecutwfc
    USE klist,      ONLY : xk


    real(kind=DP), parameter  :: fpi = 4.d0 * 3.14159265358979d0, &
         e2  = 2.d0
    real(kind=DP), parameter  :: pi =  3.14159265358979d0
    real(kind=DP) ::alfa,tpiba2,qq
    integer ::ikdiv
    tpiba2 = (fpi / 2.d0 / alat) **2
    alfa=(((alat)**2)/tpiba2+1/(ecutwfc))/2  !To check alfa
    ! alfa = omega **(2d0/3d0)
!!!!!!!!!calculate factor of divergency F(q) calculated on grid minus calculated through integration
    do ikdiv=1,nks
       qq = (xk(1,ikdiv)-xk(1,currentk))**2 + &
            (xk(2,ikdiv)-xk(2,currentk))**2 + &
            (xk(3,ikdiv)-xk(3,currentk))**2
       if (qq.gt.1.d-8) then
          exxdivfac=exxdivfac+exp(-alfa*qq)/qq
       else
          exxdivfac=exxdivfac+alfa
       endif
    end do
    exxdivfac=exxdivfac-(omega/((2*pi)**3))*2*((pi/alfa)**(1.5))


  end subroutine facdiv
end module exx
