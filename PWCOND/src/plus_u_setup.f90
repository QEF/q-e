subroutine plus_u_setup(natih, lsr)
!
! Add additional +U orbitals (if DFT+U) to the full list of projectors
!
!
  USE kinds,            ONLY : DP
  USE constants,        ONLY : rytoev
  use noncollin_module, ONLY : noncolin
  USE ldaU,             ONLY : lda_plus_U, lda_plus_u_kind, U_projection, &
                               Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_alpha, &
                               Hubbard_J0, Hubbard_beta
  use atom,             ONLY : rgrid
  USE scf,              ONLY : rho
  use radial_grids,     ONLY : ndmx
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp, atm
  USE cell_base,        ONLY : alat
  use uspp_param,       only : nhm, upf 
  USE io_global,        ONLY : stdout
  USE cond,             ONLY : norbs, nocrosl, noinss, nocrosr, tblms, taunews, &
                               nenergy, earr, nrzs, zs, tran_tot, norbf, nbrx,  &
                               cross, zpseus, zpseus_nc, betars, iofspin  
  implicit none

  integer               :: lsr, iorb, iorb1, it, iwfc, iwfc1, mesh, i, ipol, ldim, &
                           norbs_new, nocrosl_new, noinss_new, nocrosr_new, na,    &
                           natih(2,norbs), lll, kkbeta
  integer, allocatable  :: ind(:,:), tblms_new(:,:), cross_new(:,:)   
  real(DP), parameter   :: epswfc=1.d-4, eps=1.d-8
  REAL(DP)              :: r1, beta1, beta2, norm, ledge, redge 
  REAL(DP), ALLOCATABLE :: bphi(:,:), rsph(:), taunews_new(:,:),  &
                           gi(:), zpseus_new(:,:,:)
 
!-- 
! Some checks 

  if (.not.lda_plus_u) return
  if (lda_plus_u_kind.eq.1) call errore('plus_u_setup','Full LDA+U not yet implemented',1)

  WRITE( stdout, '(/,/,5x,"Simplified LDA+U calculation (l_max = ",i1, &
        & ") with parameters (eV):")') Hubbard_lmax
  WRITE( stdout, '(5x,A)') &
         "atomic species    L          U    alpha       J0     beta"
  DO it = 1, ntyp
    IF ( Hubbard_U(it) /= 0.D0 .OR. Hubbard_alpha(it) /= 0.D0 .OR. &
         Hubbard_J0(it) /= 0.D0 .OR. Hubbard_beta(it) /= 0.D0 ) THEN
       WRITE( stdout,'(5x,a6,12x,i1,2x,4f9.4)') atm(it), Hubbard_L(it), &
             Hubbard_U(it)*rytoev, Hubbard_alpha(it)*rytoev, &
             Hubbard_J0(it)*rytoev, Hubbard_beta(it)*rytoev
    END IF
  END DO

  if (U_projection.eq."pseudo") return
  if (U_projection.ne."atomic") &
   call errore('plus_u_setup','+U works only for U_projection=''pseudo'' or ''atomic'' ',1)
  if (noncolin) call errore('plus_u_setup','+U for noncollinear case not yet implemented',1)
  if (lsr.ne.2) call errore('plus_u_setup','+U atoms are allowed only in scatt. region',1)
!--

  allocate ( gi(ndmx) )
  allocate ( bphi(nbrx,ntyp) )
  allocate ( rsph(ntyp) )
  allocate ( ind(2,norbs) )

  bphi(:,:) = 0.d0
  rsph(:) = 0.d0
  ind(:,:) = 0

!--
! Calculate the total number of orbitals (beta + U WF's) 
!

  noinss_new = noinss
  norbs_new = norbs

  iorb = 1
  do while (iorb.le.norbs)
    it = tblms(1,iorb)
    if (Hubbard_U(it).ne.0.d0) then
      iorb1 = iorb
      do while (natih(1,iorb1).eq.natih(1,iorb))
        iorb1 = iorb1 + 1
      enddo
!
!     The last beta for the atom with U is provided with the
!     index of the 1st (iorb) and the last (iorb1) beta   
!
      iorb1 = iorb1 - 1
      ind(1,iorb1) = iorb
      ind(2,iorb1) = iorb1
!--
      ldim = 2*Hubbard_l(it)+1
      noinss_new = noinss_new + ldim
      norbs_new = norbs_new + ldim
      iorb = iorb1
    endif
    iorb = iorb + 1
  enddo
!--

!--
! Determine the radii of atomic U WF's
!
  do it = 1, ntyp
    if (Hubbard_U(it).ne.0.d0) then
     do iwfc = 1, upf(it)%nwfc
      if (upf(it)%lchi(iwfc).eq.Hubbard_l(it)) then
       r1 = 0.d0
       do i = 2, rgrid(it)%mesh
         r1 = max(r1, ABS(upf(it)%chi(i,iwfc)/rgrid(it)%r(i)))
       enddo
       i = rgrid(it)%mesh
       do while (abs(upf(it)%chi(i,iwfc)/rgrid(it)%r(i)).le.epswfc*r1)
         i = i - 1
       enddo
       rsph(it) = rgrid(it)%r(i) / alat
       mesh = i
      endif
     enddo
    endif
  enddo
!--

!--
! Check that all +U orbitals are totally inside the scatt. region 

  i = 0
  write(6,*) 'Scatt. region L =  ', zs(nrzs+1)
  do iorb = 1, norbs
    if (ind(2,iorb).eq.iorb) then
      it = tblms(1,iorb)
      beta1 = taunews(3,iorb)-rsph(it)
      beta2 = taunews(3,iorb)+rsph(it)
      if (beta1.le.1.d-4.or.beta2.gt.zs(nrzs+1)-1.d-4) i = 1
    endif
  enddo
  if (i.eq.1) call errore('plus_u_setup','some +U orbitals cross the boundary (not allowed) ...',1)
!--

!--
!  Calculate the integrals of betas with U atomic orbitals, 
!          bphi(i) = \sum_j q_ij <beta_j|phi> 
!
 do it = 1, ntyp
   if (Hubbard_U(it).ne.0.d0) then

     mesh = upf(it)%grid%mesh
     kkbeta = upf(it)%kkbeta
     do iwfc = 1, upf(it)%nwfc
       if (upf(it)%lchi(iwfc).eq.Hubbard_l(it)) then
         do iorb = 1, upf(it)%nbeta 
           if (upf(it)%lll(iorb).eq.Hubbard_l(it)) then
              gi(1:kkbeta)= upf(it)%beta(1:kkbeta,iorb) * &
                            upf(it)%chi (1:kkbeta,iwfc)
              call simpson (kkbeta, gi, upf(it)%grid%rab,bphi(iorb,it))
           endif  
         enddo
       endif 
     enddo
     gi(:) = 0.d0
     do iorb = 1, upf(it)%nbeta
       do iorb1 = 1, upf(it)%nbeta
         gi(iorb) = gi(iorb) + upf(it)%qqq(iorb,iorb1)*bphi(iorb1,it)
       enddo
     enddo
     bphi(1:upf(it)%nbeta,it) = gi(1:upf(it)%nbeta) 

   endif 
 enddo
!--

!--
! Allocate the arrays with all the orbitals (beta + U WF's) 
!
  allocate( taunews_new(4,norbs_new) )
  allocate( tblms_new(4,norbs_new) )
  allocate( cross_new(norbs_new, nrzs) )
  allocate( zpseus_new(2, norbs_new, norbs_new) )
  zpseus_new(:,:,:) = 0.d0
!--

!--
! Set up new extended arrays (beta + U orbitals)  
!
! iorb  --> old list
! iorb1 --> new list 
!
!                  old list                              new list
!
!           (1st atom beta        )               (1st atom beta        )
! iorb  --> (2nd atom beta        )               (         +U orbitals )
!           ( ...                 )     iorb1 --> (2nd atom beta        )
!                                                 (         +U orbitals )
!                                                 ( ...                 ) 
  iorb1 = 0
  do iorb = 1, norbs
    iorb1 = iorb1 + 1
    it = tblms(1,iorb)
    na = natih(1,iorb)
!--
!   setting up some beta arrays from old ones (just shifting)

    do i = 1, 4
      taunews_new(i,iorb1) = taunews(i,iorb)
    enddo
    do i = 1, 4
      tblms_new(i,iorb1) = tblms(i,iorb)
    enddo
    do i = 1, nrzs
      cross_new(iorb1,i) = cross(iorb,i)
    enddo
!--

!--
!   beta-beta block of zpseu (again just shifting)
!
    do i = 1, norbs
     if(natih(1,i).eq.na) then
       zpseus_new(:,i-iorb+iorb1,iorb1) = zpseus(:,i,iorb)
     endif
    enddo  
!--

!   entering into +U orbitals part (if any)

    if (ind(2,iorb).eq.iorb) then

      lll = Hubbard_l(it)
      ldim = 2*lll + 1 
!--
!     beta-beta additional block of zpseu
!
      do iwfc = ind(1,iorb), iorb
        if (tblms(3,iwfc).eq.lll) then
          do iwfc1 = ind(1,iorb), iorb
            if (tblms(3,iwfc1).eq.lll) then
             r1 = -2.d0*rho%ns(tblms(4,iwfc),tblms(4,iwfc1),iofspin,na) 
             if (tblms(4,iwfc).eq.tblms(4,iwfc1)) r1 = r1 + 1.d0
             zpseus_new(1,iwfc-iorb+iorb1,iwfc1-iorb+iorb1) = &
              zpseus_new(1,iwfc-iorb+iorb1,iwfc1-iorb+iorb1) + &
              0.5d0*Hubbard_U(it)*bphi(tblms(2,iwfc),it)*bphi(tblms(2,iwfc1),it)*r1 
            endif 
          enddo
        endif
      enddo
!--

!--
!     beta-atomic block of zpseu
!
      lll = Hubbard_l(it)
      do iwfc = ind(1,iorb), iorb
        if (tblms(3,iwfc).eq.lll) then
          do iwfc1 = 1, ldim
            r1 = -2.d0*rho%ns(tblms(4,iwfc),iwfc1,iofspin,na) 
            if (tblms(4,iwfc).eq.iwfc1) r1 = r1 + 1.d0
            zpseus_new(1,iwfc-iorb+iorb1,iorb1+iwfc1) =   &
             0.5d0*Hubbard_U(it)*bphi(tblms(2,iwfc),it)*r1    
          enddo
        endif
      enddo
!--

!--
!     atomic-beta block of zpseu 
!
      do iwfc1 = 1, ldim
        do iwfc = ind(1,iorb), iorb
          zpseus_new(1,iorb1+iwfc1,iwfc-iorb+iorb1) = &
           zpseus_new(1,iwfc-iorb+iorb1,iorb1+iwfc1)
        enddo
      enddo
!--

!--
!     atomic-atomic block of zpseu
!
      do iwfc = 1, ldim
        do iwfc1 = 1, ldim
          zpseus_new(1,iorb1+iwfc,iorb1+iwfc1) = &
           - Hubbard_U(it) * rho%ns(iwfc,iwfc1,iofspin,na)   
        enddo
        zpseus_new(1,iorb1+iwfc,iorb1+iwfc) =   &
          zpseus_new(1,iorb1+iwfc,iorb1+iwfc) + 0.5d0*Hubbard_U(it)  
      enddo
!--

!--
!     setting up some +U orbitals arrays from those of beta  
!

      do iwfc = 1, ldim

        iorb1 = iorb1 + 1

        do i = 1, 3
          taunews_new(i,iorb1) = taunews(i,iorb)
        enddo
        taunews_new(4,iorb1) = rsph(it)

        tblms_new(1,iorb1) = tblms(1,iorb)
        tblms_new(2,iorb1) = tblms(2,iorb) + 1
        tblms_new(3,iorb1) = Hubbard_l(it)
        tblms_new(4,iorb1) = iwfc

        ledge = taunews(3,iorb) - rsph(it)
        redge = taunews(3,iorb) + rsph(it)
        do i = 1, nrzs
          if (ledge.gt.zs(i+1).or.redge.lt.zs(i)) then
            cross_new(iorb1,i)=0
          else
            cross_new(iorb1,i)=1
          endif
        enddo

      enddo
!--

    endif

  enddo
!--

!--
! Add the atomic radial WF's with U to the list betars
!
  do it = 1, ntyp
    if (Hubbard_U(it).ne.0.d0) then
      do iwfc = 1, upf(it)%nwfc
        if (upf(it)%lchi(iwfc).eq.Hubbard_l(it)) then
          betars(1:rgrid(it)%mesh,upf(it)%nbeta+1,it) = &
                            upf(it)%chi(1:rgrid(it)%mesh,iwfc)
        endif
      enddo
    endif
  enddo
!--

!--
! Reallocate the orbital arrays with new dimensions and date
!
  deallocate( taunews )
  deallocate( tblms )
  deallocate( cross )
  deallocate( zpseus )

  noinss = noinss_new
  norbs = norbs_new
  norbf = norbs

  allocate( taunews(4,norbs) )
  allocate( tblms(4,norbs) )
  allocate( cross(norbs, nrzs) )

  taunews(:,:) = taunews_new(:,:)
  tblms(:,:)   = tblms_new(:,:)
  cross(:,:)   = cross_new(:,:)

  allocate( zpseus(2, norbs, norbs) )
  zpseus(:,:,:) = zpseus_new(:,:,:)
!--
  deallocate( gi )
  deallocate( bphi )
  deallocate( rsph )
  deallocate( ind )

  deallocate( taunews_new )
  deallocate( tblms_new )
  deallocate( cross_new )
  deallocate( zpseus_new )

  return
end subroutine plus_u_setup

