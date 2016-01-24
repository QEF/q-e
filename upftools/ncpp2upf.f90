!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
PROGRAM ncpp2upf
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in PWSCF format
  !     (norm-conserving) to unified pseudopotential format

  IMPLICIT NONE
  CHARACTER(len=256) filein, fileout
  !
  !
  CALL get_file ( filein )
  OPEN(unit=1,file=filein,status='old',form='formatted')
  ! underscore added to prevent conflict with read_ncpp in Modules/
  CALL read_ncpp_(1)
  CLOSE (unit=1)

  ! convert variables read from NCPP format into those needed
  ! by the upf format - add missing quantities

  CALL convert_ncpp

  fileout=trim(filein)//'.UPF'
  PRINT '(''Output PP file in US format :  '',a)', fileout

  OPEN(unit=2,file=fileout,status='unknown',form='formatted')
  CALL write_upf_v1(2)
  CLOSE (unit=2)
  STOP
20 CALL errore ('ncpp2upf', 'Reading pseudo file name ', 1)

END PROGRAM ncpp2upf

MODULE ncpp
  !
  ! All variables read from NCPP file format
  !
  ! trailing underscore means that a variable with the same name
  ! is used in module 'upf' containing variables to be written
  !
  CHARACTER(len=20) :: dft_
  CHARACTER(len=2)  :: psd_
  real(8) :: zp_
  INTEGER nlc, nnl, lmax_, lloc, nchi
  LOGICAL :: numeric, bhstype, nlcc_
  real(8) :: alpc(2), cc(2), alps(3,0:3), aps(6,0:3)
  real(8) :: a_nlcc, b_nlcc, alpha_nlcc

  real(8) :: zmesh, xmin, dx
  real(8), ALLOCATABLE::  r_(:), rab_(:)
  INTEGER :: mesh_

  real(8), ALLOCATABLE::  vnl(:,:), rho_atc_(:), rho_at_(:)
  INTEGER, ALLOCATABLE:: lchi_(:)
  real(8), ALLOCATABLE:: chi_(:,:),  oc_(:)

END MODULE ncpp
!
!     ----------------------------------------------------------
SUBROUTINE read_ncpp_(iunps)
  !     ----------------------------------------------------------
  !
  USE ncpp
  USE upf , ONLY : els
  IMPLICIT NONE
  INTEGER :: iunps
  !
  CHARACTER(len=1), DIMENSION(0:3) :: convel=(/'S','P','D','F'/)
  CHARACTER(len=2) :: label
  real (8) :: x, qe_erf
  INTEGER :: l, i, ir, nb, n
  CHARACTER (len=255) line
  EXTERNAL qe_erf

  READ(iunps, *, end=300, err=300 ) dft_
  IF (dft_(1:2)=='**') dft_ = 'PZ'

  READ (iunps, *, err=300) psd_, zp_, lmax_, nlc, nnl, nlcc_, &
       lloc, bhstype
  IF ( nlc>2 .or. nnl>3) &
       CALL errore( 'read_ncpp','Wrong nlc or nnl',1 )
  IF ( nlc* nnl < 0 ) &
       CALL errore( 'read_ncpp','nlc*nnl < 0 ? ',1 )
  IF ( zp_<=0d0 ) &
      CALL errore( 'read_ncpp','Wrong zp ',1 )
  IF ( lmax_>3.or.lmax_<0 ) &
      CALL errore( 'read_ncpp','Wrong lmax ',1 )

  IF (lloc==-1000) lloc=lmax_
  !
  !   In numeric pseudopotentials both nlc and nnl are zero.
  !
  numeric = nlc<=0 .and. nnl<=0

  IF (.not.numeric) THEN
     !
     !   read pseudopotentials in analytic form
     !
     READ(iunps, *, err=300) &
          ( alpc(i), i=1, 2 ), ( cc(i), i=1,2 )
     IF ( abs(cc(1)+cc(2)-1.d0)>1.0d-6) &
          CALL errore ('read_ncpp','wrong pseudopotential coefficients',1)
     DO l = 0, lmax_
        READ (iunps, *, err=300) &
             ( alps(i,l),i=1,3 ), (aps(i,l),i=1,6)
     ENDDO

     IF (nlcc_) THEN
        READ(iunps, *, err=300) &
             a_nlcc, b_nlcc, alpha_nlcc
        IF (alpha_nlcc<=0.d0) &
             CALL errore('read_ncpp','nlcc but alpha=0',1)
     ENDIF

     IF (bhstype) CALL bachel(alps,aps,1,lmax_)
  ENDIF

  READ(iunps, *, err=300) zmesh, xmin, dx, mesh_, nchi

  IF ( mesh_<=0) CALL errore( 'read_ncpp', 'mesh too small', 1)
  IF ( (nchi<lmax_   .and. lloc==lmax_).or.          &
       (nchi<lmax_+1 .and. lloc/=lmax_)     )        &
       CALL errore( 'read_ncpp', 'wrong no. of wfcts', 1 )
  !
  !    compute the radial mesh
  !
  ALLOCATE(  r_(mesh_))
  ALLOCATE(rab_(mesh_))

  DO ir = 1, mesh_
     x = xmin + dble(ir-1) * dx
     r_  (ir) = exp(x) / zmesh
     rab_(ir) = dx * r_(ir)
  ENDDO

  ALLOCATE(vnl(mesh_,0:lmax_))
  IF (numeric) THEN
     !
     !  read pseudopotentials in numeric form
     !
     DO l = 0, lmax_
        READ(iunps, '(a)', err=300)
        READ(iunps, *, err=300)  (vnl(ir,l),ir=1,mesh_)
     ENDDO

     ALLOCATE(rho_atc_(mesh_))
     IF(nlcc_) THEN
        READ(iunps, *, err=300) ( rho_atc_(ir), ir=1,mesh_ )
     ENDIF

  ELSE
     !
     !  convert analytic to numeric form
     !
     DO l=0,lmax_
     !
     !  DO NOT USE f90 ARRAY SYNTAX: qe_erf IS NOT AN INTRINSIC FUNCTION!!!
     !
        DO ir=1,mesh_
           vnl(ir,l)= - ( cc(1)*qe_erf(sqrt(alpc(1))*r_(ir)) +     &
                          cc(2)*qe_erf(sqrt(alpc(2))*r_(ir)) ) * zp_/r_(ir)
        ENDDO

        DO n=1,nnl
           vnl(:,l)= vnl(:,l)+ (aps(n,l)+ aps(n+3,l)*r_(:)**2 )* &
                exp(-alps(n,l)*r_(:)**2)
        ENDDO
        !
        ! convert to Rydberg
        !
        vnl(:,l) = vnl(:,l)*2.0d0
     ENDDO

     ALLOCATE(rho_atc_(mesh_))
     IF (nlcc_) THEN
        rho_atc_(:) =(a_nlcc+b_nlcc*(r_(:)**2))*exp(-alpha_nlcc*r_(:)**2)
        WHERE(abs(rho_atc_) < 1.0d-15)
           rho_atc_ = 0
        END WHERE
     ENDIF
  ENDIF
  !
  ! subtract the local part
  !
   DO l = 0, lmax_
     IF ( l/=lloc ) vnl(:,l) = vnl(:,l) - vnl(:,lloc)
  ENDDO
  !
  ! read pseudowavefunctions
  !
  ALLOCATE(lchi_(nchi), els(nchi))
  ALLOCATE(oc_(nchi))
  ALLOCATE(chi_(mesh_,nchi))
  DO nb = 1, nchi
     ! read wavefunction label and store for later
     READ(iunps, '(a)', err=300) line
     READ(iunps, *, err=300) lchi_( nb), oc_( nb )
     !
     !     Test lchi and occupation numbers
     !
     IF ( nb<=lmax_.and.lchi_(nb)+1/=nb) &
          CALL errore('read_ncpp','order of wavefunctions',nb)
     IF (lchi_(nb)>lmax_ .or. lchi_(nb)<0) &
          CALL errore('read_ncpp','wrong lchi',nb)
     IF ( oc_(nb)<0.d0 .or.            &
          oc_(nb)>2.d0*(2*lchi_(nb)+1)) &
             CALL errore('read_ncpp','wrong oc',nb)
     !
     ! parse and check wavefunction label
     READ(line,'(14x,a2)', err=222, end=222) label
     IF (label(2:2)/=convel(lchi_(nb))) GOTO 222
     DO l = 0, lmax_
        IF (label(2:2)==convel(l)) THEN
           els(nb) = label(1:2)
           GOTO 223
        ENDIF
     ENDDO
222  CONTINUE
     els(nb)   = '*'//convel(lchi_(nb))
223  CONTINUE
     !
     ! finally read the wavefunction
     READ(iunps, *, err=300) (chi_(ir,nb),ir=1,mesh_)
  ENDDO
  !
  !    compute the atomic charges
  !
  ALLOCATE(rho_at_(mesh_))
  rho_at_(:)=0.d0
  DO nb = 1, nchi
     IF( oc_(nb)/=0.d0) &
          rho_at_(:) = rho_at_(:) + oc_(nb)*chi_(:,nb)**2
  ENDDO
  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential successfully read'
  !     ----------------------------------------------------------
  RETURN

300 CALL errore('read_ncpp','pseudo file is empty or wrong',1)

END SUBROUTINE read_ncpp_

!     ----------------------------------------------------------
SUBROUTINE convert_ncpp
  !     ----------------------------------------------------------
  USE ncpp
  USE upf
  USE funct, ONLY : set_dft_from_name, get_iexch, get_icorr, get_igcx, get_igcc
  IMPLICIT NONE
  real(8), PARAMETER :: rmax = 10.0d0
  real(8), ALLOCATABLE :: aux(:)
  real(8) :: vll
  INTEGER :: kkbeta, l, iv, ir, i

  WRITE(generated, '("Generated using ld1 code (maybe, or maybe not)")')
  WRITE(date_author,'("Author: unknown    Generation date: as well")')
  comment = 'Info: automatically converted from PWSCF format'
  ! reasonable assumption
  IF (zmesh > 18) THEN
     rel = 1
  ELSE
     rel = 0
  ENDIF
  rcloc = 0.0d0
  nwfs  = nchi
  ALLOCATE( oc(nwfs), epseu(nwfs))
  ALLOCATE(lchi(nwfs), nns(nwfs) )
  ALLOCATE(rcut (nwfs), rcutus (nwfs))
  DO i=1, nwfs
     nns (i)  = 0
     lchi(i)  = lchi_(i)
     rcut(i)  = 0.0d0
     rcutus(i)= 0.0d0
     oc (i)   = oc_(i)
     epseu(i) = 0.0d0
  ENDDO
  DEALLOCATE (lchi_, oc_)

  psd = psd_
  pseudotype = 'NC'
  nlcc = nlcc_
  zp = zp_
  etotps = 0.0d0
  ecutrho=0.0d0
  ecutwfc=0.0d0
  IF ( lmax_ == lloc) THEN
     lmax = lmax_-1
  ELSE
     lmax = lmax_
  ENDIF
  nbeta= lmax_
  mesh = mesh_
  ntwfc= nchi
  ALLOCATE( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  DO i=1, nchi
     lchiw(i) = lchi(i)
     ocw(i)   = oc(i)
     elsw(i)  = els(i)
  ENDDO
  CALL set_dft_from_name(dft_)
  iexch = get_iexch()
  icorr = get_icorr()
  igcx  = get_igcx()
  igcc  = get_igcc()

  ALLOCATE(rab(mesh))
  ALLOCATE(  r(mesh))
  rab = rab_
  r = r_

  ALLOCATE (rho_atc(mesh))
  rho_atc = rho_atc_
  DEALLOCATE (rho_atc_)

  ALLOCATE (vloc0(mesh))
  vloc0(:) = vnl(:,lloc)

  IF (nbeta > 0) THEN

     ALLOCATE(ikk2(nbeta), lll(nbeta))
     kkbeta=mesh
     DO ir = 1,mesh
        IF ( r(ir) > rmax ) THEN
           kkbeta=ir
           exit
        ENDIF
     ENDDO

! make sure kkbeta is odd as required for simpson
     IF(mod(kkbeta,2) == 0) kkbeta=kkbeta-1
     ikk2(:) = kkbeta
     ALLOCATE(aux(kkbeta))
     ALLOCATE(betar(mesh,nbeta))
     ALLOCATE(qfunc(mesh,nbeta,nbeta))
     ALLOCATE(dion(nbeta,nbeta))
     ALLOCATE(qqq (nbeta,nbeta))
     qfunc(:,:,:)=0.0d0
     dion(:,:) =0.d0
     qqq(:,:)  =0.d0
     iv=0
     DO i=1,nchi
        l=lchi(i)
        IF (l/=lloc) THEN
           iv=iv+1
           lll(iv)=l
           DO ir=1,kkbeta
              betar(ir,iv)=chi_(ir,i)*vnl(ir,l)
              aux(ir) = chi_(ir,i)**2*vnl(ir,l)
           ENDDO
           CALL simpson(kkbeta,aux,rab,vll)
           dion(iv,iv) = 1.0d0/vll
        ENDIF
        IF(iv >= nbeta) exit  ! skip additional pseudo wfns
     ENDDO
     DEALLOCATE (vnl, aux)
!
!   redetermine ikk2
!
     DO iv=1,nbeta
        ikk2(iv)=kkbeta
        DO ir = kkbeta,1,-1
           IF ( abs(betar(ir,iv)) > 1.d-12 ) THEN
              ikk2(iv)=ir
              exit
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  ALLOCATE (rho_at(mesh))
  rho_at = rho_at_
  DEALLOCATE (rho_at_)

  ALLOCATE (chi(mesh,ntwfc))
  chi = chi_
  DEALLOCATE (chi_)

  RETURN
END SUBROUTINE convert_ncpp
