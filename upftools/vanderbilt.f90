
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE Vanderbilt
  !
  ! All variables read from Vanderbilt's file format
  !
  ! trailing underscore means that a variable with the same name
  ! is used in module 'upf' containing variables to be written
  !
  INTEGER :: nvalps, nang, nbeta_, kkbeta, nchi, ifpcor, keyps, &
       mesh_, iver(3), idmy(3), nnlz, ifqopt, nqf_, irel,  npf, &
       nlc, lloc
  real(8) ::  z_, zp_, exfact, etot, eloc, rcloc_, rpcor, &
       qtryc, ptryc, rinner1_
  real(8), ALLOCATABLE::  wwnlps(:), eeps(:), rinner_(:), rc(:), &
       beta(:,:), ddd0(:,:), ddd(:,:), qqq_(:,:), eee(:), rho_atc_(:), &
       r_(:), rab_(:), rho_at_(:), qfunc_(:,:,:), vloc(:), vloc_(:), &
       wf(:,:), qfcoef_(:,:,:,:)
  INTEGER, ALLOCATABLE :: lll_(:), nnlzps(:), iptype(:)
  CHARACTER(len=20):: title
END MODULE Vanderbilt
!
!     ----------------------------------------------------------
SUBROUTINE read_uspp(iunit)
  !     ----------------------------------------------------------
  !
  USE Vanderbilt
  IMPLICIT NONE
  INTEGER :: iunit
  !
  INTEGER :: i, j, k, lp
  real(8) :: rinner1
  !
  !
  READ (iunit) (iver(i),i=1,3),(idmy(i),i=1,3)
  READ (iunit) title, z_, zp_, exfact, nvalps, mesh_, etot

  ALLOCATE(nnlzps(nvalps), wwnlps(nvalps), eeps(nvalps))
  READ (iunit) (nnlzps(i),wwnlps(i),eeps(i),i=1,nvalps)

  READ (iunit) keyps, ifpcor, rinner1

  IF ( iver(1) == 1 ) THEN
     nang = nvalps
     nqf_ = 3
     nlc = 5
  ELSEIF ( iver(1) == 2 ) THEN
     nang = nvalps
     nqf_ = 3
     nlc = 2 * nvalps - 1
  ELSEIF ( iver(1) >= 3 ) THEN
     READ (iunit) nang, lloc, eloc, ifqopt, nqf_, qtryc
     nlc = 2 * nang - 1
  ENDIF

  ALLOCATE(rinner_(2*nang-1))
  rinner_(1) = rinner1
  rinner1_ = rinner1
  IF (10*iver(1)+iver(2)>=51) &
       READ (iunit) (rinner_(i),i=1,nang*2-1)

  IF ( iver(1) >= 4 ) THEN
     READ (iunit) irel
  ELSE
     irel = 0
  ENDIF

  ALLOCATE(rc(nang))
  READ (iunit) (rc(i),i=1,nang)

  READ (iunit) nbeta_,kkbeta
  !
  ALLOCATE(beta(kkbeta,nbeta_))
  ALLOCATE(qfunc_(kkbeta,nbeta_,nbeta_))
  ALLOCATE(ddd0(nbeta_,nbeta_))
  ALLOCATE(ddd (nbeta_,nbeta_))
  ALLOCATE(qqq_(nbeta_,nbeta_))
  ALLOCATE(lll_(nbeta_))
  ALLOCATE(eee(nbeta_))
  ALLOCATE(qfcoef_(nqf_,nlc,nbeta_,nbeta_))
  !
  DO j=1,nbeta_
     READ (iunit) lll_(j),eee(j),(beta(i,j),i=1,kkbeta)
     DO k=j,nbeta_
        READ (iunit) ddd0(j,k),ddd(j,k),qqq_(j,k), &
             (qfunc_(i,j,k),i=1,kkbeta), &
             ((qfcoef_(i,lp,j,k),i=1,nqf_),lp=1,2*nang-1)
     ENDDO
  ENDDO
  !
  ALLOCATE(iptype(nbeta_))
  IF (10*iver(1)+iver(2)>=72) &
       READ (iunit) (iptype(j),j=1,nbeta_),npf,ptryc
  !
  ALLOCATE(vloc_(mesh_))
  READ (iunit) rcloc_,(vloc_(i),i=1,mesh_)
  !
  ALLOCATE(rho_atc_(mesh_))
  IF (ifpcor>0) THEN
     READ (iunit) rpcor
     READ (iunit) (rho_atc_(i),i=1,mesh_)
  ENDIF
  !
  ALLOCATE(rho_at_(mesh_), vloc(mesh_))
  READ (iunit) (vloc(i),i=1,mesh_)
  READ (iunit) (rho_at_(i),i=1,mesh_)

  ALLOCATE(r_(mesh_), rab_(mesh_))
  READ (iunit) (r_(i),i=1,mesh_)
  READ (iunit) (rab_(i),i=1,mesh_)
  IF (iver(1) >= 6) THEN
     nchi = nvalps
     IF (iver(1) >= 7)  READ (iunit) nchi
     ALLOCATE(wf(mesh_,nchi))
     READ (iunit) ((wf(i,j), i=1,mesh_),j=1,nchi)
  ENDIF
  !
  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential successfully read'
  !     ----------------------------------------------------------
  !
END SUBROUTINE read_uspp
!     ----------------------------------------------------------
!     ----------------------------------------------------------
SUBROUTINE read_vdb(iunit)
  !     ----------------------------------------------------------
  !
  USE Vanderbilt
  IMPLICIT NONE
  INTEGER :: iunit
  !
  INTEGER :: i, j, k, lp
  real(8) :: rinner1
  !
  !
  READ(iunit, *) (iver(i),i=1,3),(idmy(i),i=1,3)
  READ(iunit,'(a20,3f15.9)' ) title, z_, zp_, exfact
  READ(iunit, *) nvalps, mesh_, etot

  ALLOCATE(nnlzps(nvalps), wwnlps(nvalps), eeps(nvalps))
  DO i = 1,nvalps
     READ(iunit, *) nnlzps(i), wwnlps(i), eeps(i)
  ENDDO

  READ(iunit, *)  keyps, ifpcor, rinner1

  IF ( iver(1) == 1 ) THEN
     nang = nvalps
     nqf_ = 3
     nlc = 5
  ELSEIF ( iver(1) == 2 ) THEN
     nang = nvalps
     nqf_ = 3
     nlc = 2 * nvalps - 1
  ELSEIF ( iver(1) >= 3 ) THEN
     READ(iunit, *) nang, lloc, eloc, ifqopt, nqf_, qtryc
     nlc = 2 * nang - 1
  ENDIF

  ALLOCATE(rinner_(2*nang-1))
  rinner_(1) = rinner1
  IF (10*iver(1)+iver(2)>=51) &
       READ (iunit, *) (rinner_(i),i=1,nang*2-1)
  IF ( iver(1) >= 4 ) THEN
     READ (iunit, *) irel
  ELSE
     irel = 0
  ENDIF

  ALLOCATE(rc(nang))
  READ(iunit, *) ( rc(i), i=1,nang)

  READ (iunit,* ) nbeta_, kkbeta

  ALLOCATE(beta(kkbeta,nbeta_))
  ALLOCATE(qfunc_(kkbeta,nbeta_,nbeta_))
  ALLOCATE(ddd0(nbeta_,nbeta_))
  ALLOCATE(ddd (nbeta_,nbeta_))
  ALLOCATE(qqq_(nbeta_,nbeta_))
  ALLOCATE(lll_(nbeta_))
  ALLOCATE(eee (nbeta_))
  ALLOCATE(qfcoef_(nqf_,nlc,nbeta_,nbeta_))

  DO j=1,nbeta_
     READ ( iunit, *) lll_(j)
     READ ( iunit, *) eee(j), ( beta(i,j), i=1,kkbeta )
     DO k=j,nbeta_
        READ( iunit, *) ddd0(j,k), ddd(j,k), qqq_(j,k), &
             (qfunc_(i,j,k),i=1,kkbeta),&
             ((qfcoef_(i,lp,j,k),i=1,nqf_),lp=1,2*nang-1)
     ENDDO
  ENDDO

  ALLOCATE(iptype(nbeta_))
  IF (10*iver(1)+iver(2)>=72) THEN
     READ ( iunit, * ) (iptype(i), i=1,nbeta_)
     READ ( iunit, * )  npf, ptryc
  ENDIF

  ALLOCATE(vloc_(mesh_))
  READ(iunit, *) rcloc_, ( vloc_(i), i=1,mesh_)

  ALLOCATE(rho_atc_(mesh_))
  IF ( ifpcor>0 ) THEN
     READ(iunit, *) rpcor
     READ(iunit, *) ( rho_atc_(i), i=1,mesh_)
  ENDIF

  ALLOCATE(rho_at_(mesh_), vloc(mesh_))
  READ(iunit, *)  (vloc(i), i=1,mesh_)
  READ(iunit, *)  (rho_at_(i), i=1,mesh_)

  ALLOCATE(r_(mesh_),rab_(mesh_))
  READ(iunit, *)  (r_(i), i=1,mesh_)
  READ(iunit, *)  (rab_(i),i=1,mesh_)

  IF (iver(1) >= 6) THEN
     nchi = nvalps
     IF (iver(1) >= 7)  READ (iunit, *) nchi
     ALLOCATE(wf(mesh_,nchi))
     READ (iunit, *) ((wf(i,j), i=1,mesh_),j=1,nchi)
  ENDIF

  RETURN
END SUBROUTINE read_vdb

SUBROUTINE convert_uspp
  !     ----------------------------------------------------------
  !
  USE Vanderbilt
  USE constants, ONLY : fpi
  USE upf
  IMPLICIT NONE
  INTEGER i
  CHARACTER(len=1), DIMENSION(0:3) :: convel=(/'S','P','D','F'/)

  WRITE(generated, '("Generated using Vanderbilt code, version ",3i3)') iver
  WRITE(date_author,'("Author: unknown    Generation date:",3i5)') idmy
  WRITE(comment,'("Automatically converted from original format")')
  IF (irel == 0) THEN
    rel = 0
  ELSEIF (irel == 1) THEN
    rel = 2
  ELSEIF (irel == 2) THEN
    rel = 1
  ENDIF
  rcloc = rcloc_
  nwfs = nvalps
  ALLOCATE( els(nwfs), oc(nwfs), epseu(nwfs))
  ALLOCATE(lchi(nwfs), nns(nwfs) )
  ALLOCATE(rcut (nwfs), rcutus (nwfs))
  DO i=1, nwfs
     nns (i)  = nnlzps(i)/100
     lchi(i)  = mod (nnlzps(i)/10,10)
     rcut(i)  = rinner1_
     rcutus(i)= rc(lchi(i)+1)
     oc (i) = wwnlps(i)
     WRITE(els(i),'(i1,a1)') nns(i), convel(lchi(i))
     epseu(i) = eeps(i)
  ENDDO
  DEALLOCATE (nnlzps, rc, wwnlps, eeps)

  psd = title
  IF (keyps<=2) THEN
     pseudotype = 'NC'
  ELSE
     pseudotype = 'US'
  ENDIF
  nlcc = ifpcor>0
  zp = zp_
  etotps = etot
  ecutrho=0.0d0
  ecutwfc=0.0d0
  lmax = nang - 1
  mesh = mesh_
  nbeta = nbeta_
  IF (nvalps /= nchi) THEN
     PRINT *, 'WARNING: verify info on atomic wavefunctions'
  ENDIF
  ntwfc = nchi
  ALLOCATE( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  DO i=1, min(ntwfc,nwfs)
     elsw(i) = els(i)
     ocw(i)  = oc (i)
     lchiw(i)=lchi(i)
  ENDDO

  print *, "I got this exfact", exfact
  IF ( exfact==0) THEN
     iexch=1; icorr=1; igcx=0; igcc=0 ! Perdew-Zunger
  ELSEIF ( exfact==1) THEN
     iexch=1; icorr=3; igcx=1; igcc=3 ! Becke-Lee-Yang-Parr
  ELSEIF ( exfact==2) THEN
     iexch=1; icorr=1; igcx=1; igcc=0 ! Becke88 exchange
  ELSEIF (exfact==-5.or.exfact==3) THEN
     iexch=1; icorr=1; igcx=1; igcc=1 ! Becke88-Perdew 86
  ELSEIF (exfact==-6.or.exfact==4) THEN
     iexch=1; icorr=4; igcx=2; igcc=2 ! Perdew-Wang 91
  ELSEIF (exfact== 5) THEN
     iexch=1; icorr=4; igcx=3; igcc=4 ! Perdew-Becke-Erkerhof
  ELSEIF (exfact== 6) THEN
     iexch=1; icorr=4; igcx=10; igcc=8 ! PBEsol
     print *, "I got the PBEsol correctly"
  ELSE
     WRITE (6,'("convert: wrong xc in pseudopotential ",f12.6)') exfact
     STOP
  ENDIF

  ALLOCATE (r(mesh), rab(mesh))
  r  =  r_
  rab=rab_
  DEALLOCATE (r_, rab_)

  ALLOCATE (rho_atc(mesh))
  ! Vanderbilt rho_core(r) =  4pi*r^2*rho_core(r) UPF
  rho_atc (1) = 0.d0
  rho_atc (2:mesh) = rho_atc_(2:mesh) / fpi / r(2:mesh)**2
  DEALLOCATE (rho_atc_)

  ALLOCATE (vloc0(mesh))
  vloc0(2:mesh) = vloc_(2:mesh)/r(2:mesh)
  vloc0(1) = vloc0(2)
  DEALLOCATE (vloc_)

  ALLOCATE(ikk2(nbeta), lll(nbeta))
  ikk2 = kkbeta
  lll  = lll_
  DEALLOCATE (lll_)
  ALLOCATE(betar(kkbeta,nbeta))
  betar = beta
  DEALLOCATE (beta)

  ALLOCATE(dion(nbeta,nbeta))
  dion = ddd0
  DEALLOCATE (ddd0)

  ALLOCATE(qqq(nbeta,nbeta))
  qqq = qqq_
  DEALLOCATE (qqq_)

  ALLOCATE(qfunc(mesh,nbeta,nbeta))
  qfunc(1:kkbeta,:,:) = qfunc_(1:kkbeta,:,:)
  qfunc(kkbeta+1:mesh,:,:) = 0.d0
  DEALLOCATE (qfunc_)

  nqf = nqf_
  nqlc= nlc
  ALLOCATE(rinner(nqlc))
  rinner = rinner_
  DEALLOCATE(rinner_)
  ALLOCATE(qfcoef(nqf,nqlc,nbeta,nbeta))
  qfcoef = qfcoef_
  DEALLOCATE (qfcoef_)

  ALLOCATE (rho_at(mesh))
  rho_at = rho_at_
  DEALLOCATE (rho_at_)

  ALLOCATE (chi(mesh,ntwfc))
  chi = wf
  DEALLOCATE (wf)

  RETURN
END SUBROUTINE convert_uspp

