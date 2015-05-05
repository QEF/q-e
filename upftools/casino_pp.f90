
MODULE casino_pp

  !
  ! All variables read from CASINO file format
  !
  ! trailing underscore means that a variable with the same name
  ! is used in module 'upf' containing variables to be written
  !
  USE kinds, ONLY : dp

  CHARACTER(len=20) :: dft_
  CHARACTER(len=2)  :: psd_
  REAL(dp) :: zp_
  INTEGER nlc, nnl, lmax_, lloc, nchi, rel_
  LOGICAL :: numeric, bhstype, nlcc_
  CHARACTER(len=2), ALLOCATABLE :: els_(:)
  REAL(dp) :: zmesh
  REAL(dp) :: xmin      = -7.0_dp
  REAL(dp) :: dx        = 20.0_dp/1500.0_dp
  REAL(dp) :: tn_prefac = 0.75E-6_dp
  LOGICAL  :: tn_grid   = .true.


  REAL(dp), ALLOCATABLE::  r_(:)
  INTEGER :: mesh_

  REAL(dp), ALLOCATABLE::  vnl(:,:)
  INTEGER, ALLOCATABLE:: lchi_(:), nns_(:)
  REAL(dp), ALLOCATABLE:: chi_(:,:),  oc_(:)

CONTAINS
  !
  !     ----------------------------------------------------------
  SUBROUTINE read_casino(iunps,nofiles, waveunit)
    !     ----------------------------------------------------------
    !
    !     Reads in a CASINO tabulated pp file and it's associated
    !     awfn files. Some basic processing such as removing the
    !     r factors from the potentials is also performed.


    USE kinds,  ONLY : dp
    IMPLICIT NONE
    TYPE :: wavfun_list
       INTEGER :: occ,eup,edwn, nquant, lquant
       CHARACTER(len=2) :: label
       REAL(dp), ALLOCATABLE :: wavefunc(:)
       TYPE (wavfun_list), POINTER :: p

    END TYPE wavfun_list

    TYPE :: channel_list
       INTEGER :: lquant
       REAL(dp), ALLOCATABLE :: channel(:)
       TYPE (channel_list), POINTER :: p

    END TYPE channel_list


    TYPE (channel_list), POINTER :: phead
    TYPE (channel_list), POINTER :: pptr
    TYPE (channel_list), POINTER :: ptail

    TYPE (wavfun_list), POINTER :: mhead
    TYPE (wavfun_list), POINTER :: mptr
    TYPE (wavfun_list), POINTER :: mtail

    INTEGER :: iunps, nofiles, ios
    !
    LOGICAL :: groundstate, found
    CHARACTER(len=2) :: label, rellab

    INTEGER :: l, i, ir, nb, gsorbs, j,k,m,tmp, lquant, orbs, nquant
    INTEGER, ALLOCATABLE :: gs(:,:)
    INTEGER, INTENT(in) :: waveunit(nofiles)

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
    READ(iunps,*) lloc               !reading in lloc
    IF ( zp_<=0d0 ) &
         CALL errore( 'read_casino','Wrong zp ',1 )
    IF ( lloc>3.or.lloc<0 ) &
         CALL errore( 'read_casino','Wrong lloc ',1 )


    !
    !    compute the radial mesh
    !

    DO i=1,3
       READ(iunps,*)
    ENDDO
    READ(iunps,*) mesh_   !Reading in total no. of mesh points


    ALLOCATE(  r_(mesh_))

    READ(iunps,*)
    DO i=1,mesh_
       READ(iunps,*) r_(i)
    ENDDO


    ! Read in the different channels of V_nl
    ALLOCATE(phead)
    ptail => phead
    pptr  => phead

    ALLOCATE( pptr%channel(mesh_) )
    READ(iunps, '(15x,I1,7x)') l
    pptr%lquant=l
    READ(iunps, *)  (pptr%channel(ir),ir=1,mesh_)


    DO
       READ(iunps, '(15x,I1,7x)', IOSTAT=ios) l

       IF (ios /= 0 ) THEN
          exit
       ENDIF

       ALLOCATE(pptr%p)
       pptr=> pptr%p
       ptail=> pptr
       ALLOCATE( pptr%channel(mesh_) )
       pptr%lquant=l
       READ(iunps, *)  (pptr%channel(ir),ir=1,mesh_)

    ENDDO

    !Compute the number of channels read in.
    lmax_ =-1
    pptr => phead
    DO
       IF ( .not. associated(pptr) )exit
       lmax_=lmax_+1

       pptr =>pptr%p
    ENDDO

    ALLOCATE(vnl(mesh_,0:lmax_))
    i=0
    pptr => phead
    DO
       IF ( .not. associated(pptr) )exit
       !         lchi_(i) = pptr%lquant

       DO ir=1,mesh_
          vnl(ir,i) = pptr%channel(ir)
       ENDDO
       DEALLOCATE( pptr%channel )
       pptr =>pptr%p
       i=i+1
    ENDDO

    !Clean up the linked list (deallocate it)
    DO
       IF ( .not. associated(phead) )exit
       pptr => phead
       phead => phead%p
       DEALLOCATE( pptr )
    ENDDO

    DO l = 0, lmax_
       DO ir = 1, mesh_
          vnl(ir,l) = vnl(ir,l)/r_(ir) !Removing the factor of r CASINO has
       ENDDO
       ! correcting for possible divide by zero
       IF ( r_(1) == 0 ) THEN
          vnl(1,l) = 0
       ENDIF
    ENDDO

    ALLOCATE(mhead)
    mtail => mhead

    mptr => mhead

    NULLIFY(mtail%p)
    groundstate=.true.
    DO j=1,nofiles

       DO i=1,4

          READ(waveunit(j),*)
       ENDDO

       READ(waveunit(j),*) orbs

       IF ( groundstate ) THEN

          ALLOCATE( gs(orbs,3) )

          gs = 0
          gsorbs = orbs
       ENDIF

       DO i=1,2
          READ(waveunit(j),*)
       ENDDO

       READ(waveunit(j),*) mtail%eup, mtail%edwn
       READ(waveunit(j),*)

       DO i=1,mtail%eup+mtail%edwn
          READ(waveunit(j),*) tmp, nquant, lquant

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

       READ(waveunit(j),*)
       READ(waveunit(j),*)

       DO i=1,mesh_
          READ(waveunit(j),*)
       ENDDO

       DO k=1,orbs
          READ(waveunit(j),'(13x,a2)', err=300) label
          READ(waveunit(j),*) tmp, nquant, lquant

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
                   READ(waveunit(j),*)
                ENDDO

                CYCLE
             ENDIF
          ENDIF
            IF ( allocated(mtail%wavefunc) ) THEN
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


             READ(waveunit(j), *, err=300) (mtail%wavefunc(ir),ir=1,mesh_)
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

       ALLOCATE(lchi_(nchi), els_(nchi), nns_(nchi))
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
          els_(i) = mptr%label

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


       !     ----------------------------------------------------------
       WRITE (0,'(a)') 'Pseudopotential successfully read'
       !     ----------------------------------------------------------
       RETURN

300    CALL errore('read_casino','pseudo file is empty or wrong',1)

     END SUBROUTINE read_casino

     !     ----------------------------------------------------------
     SUBROUTINE convert_casino(upf_out)
       !     ----------------------------------------------------------
       USE kinds, ONLY : dp
       USE upf_module
       USE radial_grids, ONLY: radial_grid_type, deallocate_radial_grid
       USE funct, ONLY : set_dft_from_name, get_iexch, get_icorr, &
                         get_igcx, get_igcc

       IMPLICIT NONE

       TYPE(pseudo_upf), INTENT(inout)       :: upf_out

       REAL(dp), ALLOCATABLE :: aux(:)
       REAL(dp) :: vll
       INTEGER :: kkbeta, l, iv, ir, i, nb

       WRITE(upf_out%generated, '("From a Trail & Needs tabulated &
            &PP for CASINO")')
       WRITE(upf_out%author,'("unknown")')
       WRITE(upf_out%date,'("unknown")')
       upf_out%comment = 'Info: automatically converted from CASINO &
            &Tabulated format'

       IF (rel_== 0) THEN
          upf_out%rel = 'no'
       ELSEIF (rel_==1 ) THEN
          upf_out%rel = 'scalar'
       ELSE
          upf_out%rel = 'full'
       ENDIF

       IF (xmin == 0 ) THEN
          xmin= log(zmesh * r_(2) )
       ENDIF

       ! Allocate and assign the raidal grid

       upf_out%mesh  = mesh_
       upf_out%zmesh = zmesh
       upf_out%dx    = dx
       upf_out%xmin  = xmin

       ALLOCATE(upf_out%rab(upf_out%mesh))
       ALLOCATE(  upf_out%r(upf_out%mesh))

       upf_out%r = r_
       DEALLOCATE( r_ )

       upf_out%rmax = maxval(upf_out%r)


       !
       ! subtract out the local part from the different
       ! potential channels
       !

       DO l = 0, lmax_
          IF ( l/=lloc ) vnl(:,l) = vnl(:,l) - vnl(:,lloc)
       ENDDO

       ALLOCATE (upf_out%vloc(upf_out%mesh))
       upf_out%vloc(:) = vnl(:,lloc)


       ! Compute the derivatives of the grid. The Trail and Needs
       ! grids use r(i) = (tn_prefac / zmesh)*( exp(i*dx) - 1 ) so
       ! must be treated differently to standard QE grids.

       IF ( tn_grid ) THEN
          DO ir = 1, upf_out%mesh
             upf_out%rab(ir) = dx * ( upf_out%r(ir) + tn_prefac / zmesh )
          ENDDO
       ELSE
          DO ir = 1, upf_out%mesh
             upf_out%rab(ir) = dx  * upf_out%r(ir)
          ENDDO
       ENDIF


       !
       !    compute the atomic charges
       !
       ALLOCATE (upf_out%rho_at(upf_out%mesh))
       upf_out%rho_at(:) = 0.d0

       DO nb = 1, nchi
          IF( oc_(nb)/=0.d0) THEN
             upf_out%rho_at(:) = upf_out%rho_at(:) +&
               &  oc_(nb)*chi_(:,nb)**2
          ENDIF
       ENDDO

       ! This section deals with the pseudo wavefunctions.
       ! These values are just given directly to the pseudo_upf structure
       upf_out%nwfc  = nchi

       ALLOCATE( upf_out%oc(upf_out%nwfc), upf_out%epseu(upf_out%nwfc) )
       ALLOCATE( upf_out%lchi(upf_out%nwfc), upf_out%nchi(upf_out%nwfc) )
       ALLOCATE( upf_out%els(upf_out%nwfc) )
       ALLOCATE( upf_out%rcut_chi(upf_out%nwfc) )
       ALLOCATE( upf_out%rcutus_chi (upf_out%nwfc) )

       DO i=1, upf_out%nwfc
          upf_out%nchi(i)  = nns_(i)
          upf_out%lchi(i)  = lchi_(i)
          upf_out%rcut_chi(i)  = 0.0d0
          upf_out%rcutus_chi(i)= 0.0d0
          upf_out%oc (i)   = oc_(i)
          upf_out%els(i) = els_(i)
          upf_out%epseu(i) = 0.0d0
       ENDDO
       DEALLOCATE (lchi_, oc_, nns_)

       upf_out%psd = psd_
       upf_out%typ = 'NC'
       upf_out%nlcc = nlcc_
       upf_out%zp = zp_
       upf_out%etotps = 0.0d0
       upf_out%ecutrho=0.0d0
       upf_out%ecutwfc=0.0d0
       upf_out%lloc=lloc

       IF ( lmax_ == lloc) THEN
          upf_out%lmax = lmax_-1
       ELSE
          upf_out%lmax = lmax_
       ENDIF
       upf_out%nbeta = lmax_

       ALLOCATE ( upf_out%els_beta(upf_out%nbeta) )
       ALLOCATE ( upf_out%rcut(upf_out%nbeta) )
       ALLOCATE ( upf_out%rcutus(upf_out%nbeta) )

       upf_out%rcut=0.0d0
       upf_out%rcutus=0.0d0
       upf_out%dft =dft_


       IF (upf_out%nbeta > 0) THEN

          ALLOCATE(upf_out%kbeta(upf_out%nbeta), upf_out%lll(upf_out%nbeta))
          upf_out%kkbeta=upf_out%mesh
          DO ir = 1,upf_out%mesh
             IF ( upf_out%r(ir) > upf_out%rmax ) THEN
                upf_out%kkbeta=ir
                exit
             ENDIF
          ENDDO

          ! make sure kkbeta is odd as required for simpson
          IF(mod(upf_out%kkbeta,2) == 0) upf_out%kkbeta=upf_out%kkbeta-1
          upf_out%kbeta(:) = upf_out%kkbeta
          ALLOCATE(aux(upf_out%kkbeta))
          ALLOCATE(upf_out%beta(upf_out%mesh,upf_out%nbeta))
          ALLOCATE(upf_out%dion(upf_out%nbeta,upf_out%nbeta))

          upf_out%dion(:,:) =0.d0

          iv=0
          DO i=1,upf_out%nwfc
             l=upf_out%lchi(i)
             IF (l/=upf_out%lloc) THEN
                iv=iv+1
                upf_out%els_beta(iv)=upf_out%els(i)
                upf_out%lll(iv)=l
                DO ir=1,upf_out%kkbeta

                   upf_out%beta(ir,iv)=chi_(ir,i)*vnl(ir,l)
                   aux(ir) = chi_(ir,i)**2*vnl(ir,l)

                ENDDO
                CALL simpson(upf_out%kkbeta,aux,upf_out%rab,vll)
                upf_out%dion(iv,iv) = 1.0d0/vll
             ENDIF

             IF(iv >= upf_out%nbeta) exit  ! skip additional pseudo wfns
          ENDDO


          DEALLOCATE (vnl, aux)

          !
          !   redetermine ikk2
          !
          DO iv=1,upf_out%nbeta
             upf_out%kbeta(iv)=upf_out%kkbeta
             DO ir = upf_out%kkbeta,1,-1
                IF ( abs(upf_out%beta(ir,iv)) > 1.d-12 ) THEN
                   upf_out%kbeta(iv)=ir
                   exit
                ENDIF
             ENDDO
          ENDDO
       ENDIF

       ALLOCATE (upf_out%chi(upf_out%mesh,upf_out%nwfc))
       upf_out%chi = chi_
       DEALLOCATE (chi_)

       RETURN
     END SUBROUTINE convert_casino


     SUBROUTINE write_casino_tab(upf_in, grid)

       USE upf_module
       USE radial_grids, ONLY: radial_grid_type, deallocate_radial_grid


       IMPLICIT NONE

       TYPE(pseudo_upf), INTENT(in)       :: upf_in
       TYPE(radial_grid_type), INTENT(in) :: grid
       INTEGER :: i, lp1

       INTEGER, EXTERNAL :: atomic_number

       WRITE(6,*) "Converted Pseudopotential in REAL space for ", upf_in%psd
       WRITE(6,*) "Atomic number and pseudo-charge"
       WRITE(6,"(I3,F8.2)") atomic_number( upf_in%psd ),upf_in%zp
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

       lp1 = size ( vnl, 2 )
       DO i=1,lp1
          WRITE(6, "(A,I1,A)") 'r*potential (L=',i-1,') in Ry'
          WRITE(6, "(T4,E22.15)") vnl(:,i)
       ENDDO

     END SUBROUTINE write_casino_tab

     SUBROUTINE conv_upf2casino(upf_in,grid)

       USE upf_module
       USE radial_grids, ONLY: radial_grid_type, deallocate_radial_grid


       IMPLICIT NONE

       TYPE(pseudo_upf), INTENT(in)       :: upf_in
       TYPE(radial_grid_type), INTENT(in) :: grid
       INTEGER :: i, l, channels

       REAL(dp), PARAMETER :: offset=1E-20_dp
       !This is an offset added to the wavefunctions to
       !eliminate any divide by zeros that may be caused by
       !zeroed wavefunction terms.

       channels=upf_in%nbeta+1
       ALLOCATE ( vnl(grid%mesh,channels) )

       !Set up the local component of each channel
       DO i=1,channels
          vnl(:,i)=grid%r(:)*upf_in%vloc(:)
       ENDDO


       DO i=1,upf_in%nbeta
          l=upf_in%lll(i)+1

          !Check if any wfc components have been zeroed
          !and apply the offset IF they have

          IF ( minval(abs(upf_in%chi(:,l))) /= 0 ) THEN
             vnl(:,l)= (upf_in%beta(:,l)/(upf_in%chi(:,l)) &
                  *grid%r(:)) + vnl(:,l)
          ELSE
             WRITE(0,"(A,ES10.3,A)") 'Applying ',offset , ' offset to &
                  &wavefunction to avoid divide by zero'
             vnl(:,l)= (upf_in%beta(:,l)/(upf_in%chi(:,l)+offset) &
                  *grid%r(:)) + vnl(:,l)
          ENDIF

       ENDDO

     END SUBROUTINE conv_upf2casino

END MODULE casino_pp
