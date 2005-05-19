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
    integer          :: ibnd, ik, im , ig, inrxx
    real(kind=DP), parameter  :: fpi = 4.d0 * 3.14159265358979d0, &
                               e2  = 2.d0
    real(kind=DP) :: tpiba2, qq
    tpiba2 = (fpi / 2.d0 / alat) **2
  
    call start_clock ('exx')
    ! write (*,*) exx_nwordwfc,lda,n,m, lda*n

    do im=1,m !for each band of psi (the k cycle is outside band)

       result(:)= ( 0.D0, 0.D0 )
       temppsic(:) = ( 0.D0, 0.D0 )
       temppsic(nls(igk(1:npw))) = psi(1:npw,im)
       CALL cft3s( temppsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
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
                fac(ig)=0.d0
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
       !adds it to hpsi
       hpsi(1:npw,im)=hpsi(1:npw,im) - exxalfa*result(nls(igk(1:npw)))
    end do
    call stop_clock ('exx')

  end subroutine vexx

  subroutine brutalexxinit()

    !This subroutine is run before the first H_psi()
    !It saves the wavefunctions for the right density matrix. 
    !It saves all the wavefunctions in an only file called prefix.exx
    USE wavefunctions_module, ONLY : evc  
    USE io_files,             ONLY : nwordwfc
    USE io_files, ONLY : prefix
    use io_files, only : tmp_dir,iunwfc
    USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et
    integer :: ios,ik
    complex (kind=dp) :: tempevc( npwx, nbnd )
    call start_clock ('exxinit')
    exx_nwordwfc=nwordwfc
    exx_file = TRIM( tmp_dir ) // TRIM( prefix ) //'.exx'

    if (.not.exxstart) then 
       open (unit=iunexx, file=exx_file, form="unformatted", &
            action="readwrite", status="scratch", &
            recl=DIRECT_IO_FACTOR*exx_nwordwfc, access="direct")
    endif

    DO ik = 1, nks     
       call davcio (tempevc, nwordwfc, iunwfc, ik, -1)
       CALL davcio (tempevc, nwordwfc, iunexx, ik, 1 )

    end do
    exxstart=.true.
    call stop_clock ('exxinit')  
  end subroutine brutalexxinit

  subroutine brutalvexx(lda, n, m, psi, hpsi)
    !This routine calculates V_xx \Psi
    !It uses the \rho extimated in exxinit 

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
    USE wvfct,          ONLY : npwx, nbnd, nbndx
    USE cell_base,      ONLY : tpiba2, omega
    USE constants,      ONLY : fpi,e2
    USE klist,          ONLY : xk,wk
    USE gvect,          ONLY : g, nl
    INTEGER          :: lda, n, m, kspi
    COMPLEX(KIND=DP) :: psi(lda,m) 
    COMPLEX(KIND=DP) :: hpsi(lda,m)   
    REAL (KIND=DP) :: fac
    INTEGER ::  ibnd, igp, igpp, ikp,ig,ik
    COMPLEX(KIND=DP) :: phi(npwx, nbnd)   
    COMPLEX(KIND=DP) :: temppsi(lda,m)   
    REAL (KIND=DP) :: sqfac

    CALL start_clock ('exx')
    temppsi=0D0
    fac= fpi * e2 *omega / (2*tpiba2)
    write (*,*) 'fac=   ', fac
    do ik =1, nks
       write (*,*) ik
       CALL davcio( phi, exx_nwordwfc, iunexx, ik, -1 )
       do ibnd=1,m
          do igp=1,npwx
             do igpp=1,npwx
                do ig=1,npwx
                   sqfac=(xk(1,ik)-xk(1,currentk)+g(1,igp)-g(1,ig))**2 + &
                         (xk(2,ik)-xk(2,currentk)+g(2,igp)-g(2,ig))**2 + &
                         (xk(3,ik)-xk(3,currentk)+g(3,igp)-g(3,ig))**2
                   if ((sqfac).ne.0) then
                      temppsi(ig,ibnd)=temppsi(ig,ibnd) + &
                                       (phi(igp,ibnd)*conjg(phi(igpp,ibnd))* &
                                        psi(nl(ig+igp-igpp),ibnd)/(((sqfac))))
                   endif
                end do
             end do
          end do
       end do
    end do
    do ibnd=1,m
       do ig=1,npwx
          write (*,*) temppsi(ig,ibnd),'     ',hpsi(ig,ibnd)
       end do
    end do

    do ibnd=1,m
       do ig=1,npwx
          hpsi(ig,ibnd)=fac*temppsi(ig,ibnd)+hpsi(ig,ibnd)
       end do
    end do
    call stop_clock ('exx')

  end subroutine brutalvexx

end module exx
