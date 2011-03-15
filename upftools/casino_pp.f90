
MODULE casino_pp

  !
  ! All variables read from CASINO file format
  !
  ! trailing underscore means that a variable with the same name
  ! is used in module 'upf' containing variables to be written
  !
  USE kinds, ONLY : DP

  CHARACTER(len=20) :: dft_
  CHARACTER(len=2)  :: psd_
  REAL(DP) :: zp_
  INTEGER nlc, nnl, lmax_, lloc, nchi, rel_
  LOGICAL :: numeric, bhstype, nlcc_
  REAL(DP) :: alpc(2), cc(2), alps(3,0:3), aps(6,0:3)
  REAL(DP) :: a_nlcc, b_nlcc, alpha_nlcc

  REAL(DP) :: zmesh, xmin, dx
  REAL(DP), ALLOCATABLE::  r_(:), rab_(:)
  INTEGER :: mesh_

  REAL(DP), ALLOCATABLE::  vnl(:,:), rho_atc_(:), rho_at_(:)
  INTEGER, ALLOCATABLE:: lchi_(:), nns_(:)
  REAL(DP), ALLOCATABLE:: chi_(:,:),  oc_(:)

CONTAINS
  !
  !     ----------------------------------------------------------
  SUBROUTINE read_casino(iunps,nofiles)
    !     ----------------------------------------------------------
    !

    USE upf , ONLY : els
    USE kinds,  ONLY : DP
    IMPLICIT NONE
    TYPE :: wavfun_list
       INTEGER :: occ,eup,edwn, nquant, lquant
       CHARACTER(len=2) :: label
#ifdef __STD_F95
       REAL(DP), POINTER :: wavefunc(:)
#else
       REAL(DP), ALLOCATABLE :: wavefunc(:)
#endif
       TYPE (wavfun_list), POINTER :: p

    END TYPE wavfun_list

    TYPE (wavfun_list), POINTER :: mhead
    TYPE (wavfun_list), POINTER :: mptr
    TYPE (wavfun_list), POINTER :: mtail

    INTEGER :: iunps, nofiles
    !
    LOGICAL :: groundstate, found
    CHARACTER(len=2) :: label, rellab
    REAL(DP), PARAMETER :: r_exp=20._DP/1500._DP
    INTEGER :: l, i, ir, nb, gsorbs, j,k,m,tmp, lquant, orbs, nquant
    INTEGER, ALLOCATABLE :: gs(:,:)

    NULLIFY (  mhead, mptr, mtail )
    dft_ = 'HF'   !Hardcoded at the moment should eventually be HF anyway

    nlc = 0              !These two values are always 0 for numeric pps
    nnl = 0              !so lets just hard code them

    nlcc_ = .false.       !Again these two are alwas false for CASINO pps
    bhstype = .false.



    READ(iunps,'(a2,35x,a2)') rellab, psd_
    READ(iunps,*)
    IF ( rellab == 'DF' ) THEN
       rel_=1
    ELSE
       rel_=0
    ENDIF

    READ(iunps,*) zmesh,zp_  !Here we are reading zmesh (atomic #) and
    DO i=1,3                 !zp_ (pseudo charge)
       READ(iunps,*)
    ENDDO
    READ(iunps,*) lmax_               !reading in lmax
    IF ( zp_<=0d0 ) &
         CALL errore( 'read_casino','Wrong zp ',1 )
    IF ( lmax_>3.or.lmax_<0 ) &
         CALL errore( 'read_casino','Wrong lmax ',1 )

    lloc=lmax_ !think lloc shoudl always = lmax for this case yes/no ??

    !
    !    compute the radial mesh
    !

    DO i=1,3
       READ(iunps,*)
    ENDDO
    READ(iunps,*) mesh_   !Reading in total no. of mesh points


    ALLOCATE(  r_(mesh_))
    ALLOCATE(rab_(mesh_))
    READ(iunps,*)
    DO i=1,mesh_
       READ(iunps,*) r_(i)
    ENDDO
    DO ir = 1, mesh_
       rab_(ir) = r_exp  * r_(ir) !hardcoded at the moment
    ENDDO


    ALLOCATE(vnl(mesh_,0:lmax_))

    DO l = 0, lmax_
       READ(iunps, '(a)', err=300)
       READ(iunps, *, err=300)  (vnl(ir,l),ir=1,mesh_)
    ENDDO

    DO l = 0, lmax_
       DO ir = 1, mesh_
          vnl(ir,l) = vnl(ir,l)/r_(ir) !Removing the factor of r CASINO has
       ENDDO
       vnl(1,l) = 0                   !correcting for the divide by zero
    ENDDO

    ALLOCATE(rho_atc_(mesh_))
    IF(nlcc_) THEN
       READ(iunps, *, err=300) ( rho_atc_(ir), ir=1,mesh_ )
    ENDIF

    !
    ! subtract the local part
    !

    DO l = 0, lmax_
       IF ( l/=lloc ) vnl(:,l) = vnl(:,l) - vnl(:,lloc)
    ENDDO


    ALLOCATE(mhead)
    mtail => mhead

    mptr => mhead

    NULLIFY(mtail%p)
    groundstate=.true.
    DO j=1,nofiles

       DO i=1,4

          READ(j,*)
       ENDDO

       READ(j,*) orbs

       IF ( groundstate ) THEN

          ALLOCATE( gs(orbs,3) )

          gs = 0
          gsorbs = orbs
       ENDIF

       DO i=1,2
          READ(j,*)
       ENDDO

       READ(j,*) mtail%eup, mtail%edwn
       READ(j,*)

       DO i=1,mtail%eup+mtail%edwn
          READ(j,*) tmp, nquant, lquant

          IF ( groundstate ) THEN
             found = .true.

             DO m=1,orbs

                IF ( (nquant==gs(m,1) .and. lquant==gs(m,2)) ) THEN
                   gs(m,3) = gs(m,3) + 1
                   exit
                ENDIF

                found = .false.

             ENDDO

             IF (.not. found ) THEN

                DO m=1,orbs

                   IF ( gs(m,1) == 0 ) THEN
                      gs(m,1) = nquant
                      gs(m,2) = lquant
                      gs(m,3) = 1

                      exit
                   ENDIF

                ENDDO

             ENDIF

          ENDIF

       ENDDO

       READ(j,*)
       READ(j,*)

       DO i=1,mesh_
          READ(j,*)
       ENDDO

       DO k=1,orbs
          READ(j,'(13x,a2)', err=300) label
          READ(j,*) tmp, nquant, lquant

          IF ( .not. groundstate ) THEN
             found = .false.

             DO m = 1,gsorbs

                IF ( nquant == gs(m,1) .and. lquant == gs(m,2) ) THEN
                   found = .true.
                   exit
                ENDIF
             ENDDO
             mptr => mhead
             DO
                IF ( .not. associated(mptr) )exit
                IF ( nquant == mptr%nquant .and. lquant == mptr%lquant ) found = .true.
                mptr =>mptr%p
             ENDDO
             IF ( found ) THEN
                DO i=1,mesh_
                   READ(j,*)
                ENDDO

                CYCLE
             ENDIF
          ENDIF
#ifdef __STD_F95
          IF ( associated(mtail%wavefunc) ) THEN
#else
             IF ( allocated(mtail%wavefunc) ) THEN
#endif
                ALLOCATE(mtail%p)
                mtail=>mtail%p
                NULLIFY(mtail%p)
                ALLOCATE( mtail%wavefunc(mesh_) )
             ELSE
                ALLOCATE( mtail%wavefunc(mesh_) )
             ENDIF
             mtail%label = label
             mtail%nquant = nquant
             mtail%lquant = lquant


             READ(j, *, err=300) (mtail%wavefunc(ir),ir=1,mesh_)
          ENDDO
          groundstate = .false.
       ENDDO

       nchi =0
       mptr => mhead
       DO
          IF ( .not. associated(mptr) )exit
          nchi=nchi+1

          mptr =>mptr%p
       ENDDO

       ALLOCATE(lchi_(nchi), els(nchi), nns_(nchi))
       ALLOCATE(oc_(nchi))
       ALLOCATE(chi_(mesh_,nchi))
       oc_ = 0

       !Sort out the occupation numbers
       DO i=1,gsorbs
          oc_(i)=gs(i,3)
       ENDDO
       DEALLOCATE( gs )

       i=1
       mptr => mhead
       DO
          IF ( .not. associated(mptr) )exit
          nns_(i) = mptr%nquant
          lchi_(i) = mptr%lquant
          els(i) = mptr%label

          DO ir=1,mesh_

             chi_(ir:,i) = mptr%wavefunc(ir)
          ENDDO
          DEALLOCATE( mptr%wavefunc )
          mptr =>mptr%p
          i=i+1
       ENDDO

       !Clean up the linked list (deallocate it)
       DO
          IF ( .not. associated(mhead) )exit
          mptr => mhead
          mhead => mhead%p
          DEALLOCATE( mptr )
       ENDDO


       !
       !    compute the atomic charges
       !
       ALLOCATE(rho_at_(mesh_))
       rho_at_(:)=0.d0
       DO nb = 1, nchi
          IF( oc_(nb)/=0.d0) &
               &  rho_at_(:) = rho_at_(:) + oc_(nb)*chi_(:,nb)**2
       ENDDO
       !     ----------------------------------------------------------
       WRITE (6,'(a)') 'Pseudopotential successfully read'
       !     ----------------------------------------------------------
       RETURN

300    CALL errore('read_casino','pseudo file is empty or wrong',1)

     END SUBROUTINE read_casino

     !     ----------------------------------------------------------
     SUBROUTINE convert_casino
       !     ----------------------------------------------------------
       USE kinds, ONLY : DP

       USE upf
       USE funct, ONLY : set_dft_from_name, get_iexch, get_icorr, get_igcx, get_igcc
       IMPLICIT NONE
       REAL(DP), PARAMETER :: rmax = 10.0d0
       REAL(DP), ALLOCATABLE :: aux(:)
       REAL(DP) :: vll
       INTEGER :: kkbeta, l, iv, ir, i

       WRITE(generated, '("From a Trail & Needs tabulated PP for CASINO")')
       WRITE(date_author,'("Author: unknown    Generation date: as well")')
       comment = 'Info: automatically converted from CASINO Tabulated format'

       rel = rel_


       rcloc = 0.0d0
       nwfs  = nchi
       ALLOCATE( oc(nwfs), epseu(nwfs))
       ALLOCATE(lchi(nwfs), nns(nwfs) )
       ALLOCATE(rcut (nwfs), rcutus (nwfs))
       DO i=1, nwfs
          nns (i)  = nns_(i)
          lchi(i)  = lchi_(i)
          rcut(i)  = 0.0d0
          rcutus(i)= 0.0d0
          oc (i)   = oc_(i)
          epseu(i) = 0.0d0
       ENDDO
       DEALLOCATE (lchi_, oc_, nns_)

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
     END SUBROUTINE convert_casino


     SUBROUTINE write_casino_tab(upf_in, grid)

       USE upf_module
       USE radial_grids, ONLY: radial_grid_type, deallocate_radial_grid 


       IMPLICIT NONE

       TYPE(pseudo_upf), INTENT(IN)       :: upf_in
       TYPE(radial_grid_type), INTENT(IN) :: grid
       INTEGER :: i

       INTEGER, EXTERNAL :: atomic_number

       WRITE(6,*) "Converted Pseudopotential in REAL space for ", upf_in%psd
       WRITE(6,*) "Atomic number and pseudo-charge"
       WRITE(6,"(I3,F5.2)") atomic_number( upf_in%psd ),upf_in%zp  
       WRITE(6,*) "Energy units (rydberg/hartree/ev):"
       WRITE(6,*) "rydberg"
       WRITE(6,*) "Angular momentum of local component (0=s,1=p,2=d..)"
       WRITE(6,"(I2)") upf_in%lloc
       WRITE(6,*) "NLRULE override (1) VMC/DMC (2) config gen (0 ==> &
            &input/default VALUE)"
       WRITE(6,*) "0 0"
       WRITE(6,*) "Number of grid points"
       WRITE(6,*) grid%mesh
       WRITE(6,*) "R(i) in atomic units"
       WRITE(6, "(T4,E22.15)") grid%r(:)

       DO i=1,SIZEOF(vnl)/SIZEOF(vnl(:,1))
          WRITE(6, "(A,I1,A)") 'r*potential (L=',i-1,') in Ry'
          WRITE(6, "(T4,E22.15)") vnl(:,i)
       END DO

     END SUBROUTINE write_casino_tab

     SUBROUTINE conv_upf2casino(upf_in,grid)

       USE upf_module
       USE radial_grids, ONLY: radial_grid_type, deallocate_radial_grid 


       IMPLICIT NONE

       TYPE(pseudo_upf), INTENT(IN)       :: upf_in
       TYPE(radial_grid_type), INTENT(IN) :: grid
       INTEGER :: i, l, channels

       REAL(DP), PARAMETER :: offset=1E-20_DP
       !This is an offset added to the wavefunctions to 
       !eliminate any divide by zeros that may be caused by 
       !zeroed wavefunction terms.

       channels=upf_in%nbeta+1
       ALLOCATE ( vnl(grid%mesh,channels) )

       !Set up the local component of each channel
       DO i=1,channels
          vnl(:,i)=grid%r(:)*upf_in%vloc(:)
       END DO


       DO i=1,upf_in%nbeta
          l=upf_in%lll(i)+1

          !Check if any wfc components have been zeroed 
          !and apply the offset IF they have

          IF ( MINVAL(ABS(upf_in%chi(:,l))) .ne. 0 ) THEN
             vnl(:,l)= (upf_in%beta(:,l)/(upf_in%chi(:,l)) &
                  *grid%r(:)) + vnl(:,l)
          ELSE
             WRITE(0,"(A,ES10.3,A)") 'Applying ',offset , ' offset to &
                  &wavefunction to avoid divide by zero'
             vnl(:,l)= (upf_in%beta(:,l)/(upf_in%chi(:,l)+offset) &
                  *grid%r(:)) + vnl(:,l)
          END IF

       END DO

     END SUBROUTINE conv_upf2casino

END MODULE casino_pp
