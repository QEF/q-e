!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#include "machine.h"
#undef DEBUG
!-----------------------------------------------------------------------
subroutine mix_rho (rhout, rhoin, nsout, nsin, alphamix, dr2, iter, &
                    n_iter, filename, conv)
  !-----------------------------------------------------------------------
  !
  ! Modified Broyden's method for charge density mixing
  !             d.d. johnson prb 38, 12807 (1988)
  ! On output: the mixed density is in rhoin, rhout is UNCHANGED
  !
  use parameters, only : DP
  use pwcom
  use gamma
  !
  !   First the I/O variable
  !
  character (len=42) ::  &
                filename     !  (in) I/O filename for mixing history
                             !  if absent everything is kept in memory
  integer ::    &
                iter,       &!  (in)  counter of the number of iterations
                n_iter       !  (in)  numb. of iterations used in mixing

  real (kind=DP) :: &
                rhout(nrxx,nspin), &! (in) the "out" density; (out) rhout-rhoin
                rhoin(nrxx,nspin), &! (in) the "in" density; (out) the new dens.
                nsout(5,5,nspin,nat), &!
                nsin(5,5,nspin,nat),  &!
                alphamix,          &! (in) mixing factor
                dr2                 ! (out) the estimated errr on the energy

  logical ::    &
                conv        ! (out) if true the convergence has been reached
  !
  integer, parameter:: &
                maxmix = 25 ! max number of iterations for charge mixing

  !
  !   Here the local variables
  !
  integer ::    &
                iunmix,    &! I/O unit number of charge density file
                iunmix2,   &! I/O unit number of ns file
                iunit,     &! counter on I/O unit numbers
                iter_used, &! actual number of iterations used
                ipos,      &! index of the present iteration
                inext,     &! index of the next iteration
                i, j,      &! counters on number of iterations
                is,        &! counter on spin component
                ig,        &! counter on G-vectors
                iwork(maxmix),&! dummy array used as output by libr. routines
                info        ! flag saying if the exec. of libr. routines was ok

  complex (kind=DP), allocatable  :: aux(:), rhocin(:,:), rhocout(:,:), &
                rhoinsave(:), rhoutsave(:),  &
                nsinsave(:,:,:,:),  nsoutsave(:,:,:,:)
  complex (kind=DP), allocatable, save :: df(:,:), dv(:,:), &
                                          df_ns(:,:,:,:,:), dv_ns(:,:,:,:,:)
                ! aux(nrxx)            : auxiliary array used for FFT
                ! rhocin(ngm0,nspin)
                ! rhocout(ngm0,nspin)
                ! rhoinsave(ngm0*nspin): work space
                ! rhoutsave(ngm0*nspin): work space
                ! df(ngm0*nspin,n_iter): information from preceding iterations
                ! dv(ngm0*nspin,n_iter):    "  "       "     "        "  "
                ! df_ns(5,5,nspin,nat,n_iter):  idem
                ! dv_ns(5,5,nspin,nat,n_iter):  idem

  real (kind=DP) :: betamix(maxmix,maxmix), gamma0, work(maxmix), &
                    fn_dehar, dehar

  logical ::    &
                saveonfile, &! save intermediate steps on file "filename"
                opnd,       &! if true the file is already opened
                exst         ! if true the file exists

  real (kind=DP) :: rho_dot_product

  external DCOPY, DSYTRF, DSYTRI, DSCAL
  external diropn, davcio, rho_dot_product, fn_dehar

  call start_clock('mix_rho')

  if (iter.lt.1) call errore('mix_rho','iter is wrong',1)
  if (n_iter.gt.maxmix) call errore('mix_rho','n_iter too big',1)

  saveonfile=filename.ne.' '

!  call DAXPY(nrxx*nspin,-1.d0,rhoin,1,rhout,1)

    allocate(aux(nrxx), rhocin(ngm0,nspin), rhocout(ngm0,nspin))

  do is=1,nspin
     aux(:) = DCMPLX(rhoin(:,is),0.d0)
     call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
     do ig=1,ngm0
        rhocin(ig,is) = aux(nl(ig))
     end do
     aux(:) = DCMPLX(rhout(:,is),0.d0)
     call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
     do ig=1,ngm0
        rhocout(ig,is) = aux(nl(ig)) - rhocin(ig,is)
     end do
  end do
  if (lda_plus_u) nsout(:,:,:,:) = nsout(:,:,:,:) - nsin(:,:,:,:)

  dr2=rho_dot_product(rhocout,rhocout)
  conv = dr2.lt.tr2
  dehar = fn_dehar(rhocout)
#ifdef DEBUG
  if (conv) then
     write (6,100) dr2, rho_dot_product(rhocout,rhocout)
     write (6,'(" dehar =",f15.8)') dehar
  end if
#endif

  if (saveonfile) then
     do iunit=99,1,-1
        inquire(unit=iunit,opened=opnd)
        iunmix=iunit
        if (.not.opnd) go to 10
     end do
     call errore('mix_rho','free unit not found?!?',1)
10   continue
     if (lda_plus_u) then
        do iunit=iunmix-1,1,-1
           inquire(unit=iunit,opened=opnd)
           iunmix2=iunit
           if (.not.opnd) go to 20
        end do
        call errore('mix_rho','second free unit not found?!?',1)
20      continue
     end if
     if (conv) then
        call diropn (iunmix, filename, 2*ngm0*nspin, exst)
        close (unit=iunmix, status='delete')
        if (lda_plus_u) then
           call diropn (iunmix2, trim(filename)//'.ns',25*nspin*nat, exst)
           close (unit=iunmix2, status='delete')
        end if
        deallocate (aux, rhocin, rhocout)
        call stop_clock('mix_rho')
        return
     end if

     call diropn(iunmix,filename,2*ngm0*nspin,exst)
     if (lda_plus_u) call diropn (iunmix2, trim(filename)//'.ns',25*nspin*nat, exst)

    if (iter.gt.1 .and. .not.exst) then
        call errore('mix_rho','file not found, restarting',-1)
        iter=1
     end if
     allocate (df(ngm0*nspin,n_iter), dv(ngm0*nspin,n_iter))
     if (lda_plus_u) &
        allocate (df_ns(5,5,nspin,nat,n_iter), dv_ns(5,5,nspin,nat,n_iter))
 else
     if (iter.eq.1) then
        allocate (df(ngm0*nspin,n_iter), dv(ngm0*nspin,n_iter))
        if (lda_plus_u) &
           allocate (df_ns(5,5,nspin,nat,n_iter), dv_ns(5,5,nspin,nat,n_iter))
     end if
     if (conv) then
        if (lda_plus_u) deallocate(df_ns, dv_ns)
        deallocate (df, dv)
        deallocate (aux, rhocin, rhocout)
        call stop_clock('mix_rho')
        return
     end if
     allocate (rhoinsave(ngm0*nspin), rhoutsave(ngm0*nspin))
     if (lda_plus_u) allocate(nsinsave(5,5,nspin,nat),nsoutsave(5,5,nspin,nat))
  end if
  !
  ! copy only the high frequency Fourier component into rhoin
  !                                                (NB: rhout=rhout-rhoin)
  !
  rhoin(:,:) = rhout(:,:)
  do is=1,nspin
     aux(:) = (0.d0, 0.d0)
     do ig=1,ngm0
        aux(nl (ig)) = rhocin(ig,is)+rhocout(ig,is)
        aux(nlm(ig)) = conjg(aux(nl (ig)))
     end do
     call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
     rhoin(:,is) = rhoin(:,is) - DREAL(aux(:))
  end do
  !
  ! iter_used = iter-1  if iter <= n_iter
  ! iter_used = n_iter  if iter >  n_iter
  !
  iter_used=min(iter-1,n_iter)
  !
  ! ipos is the position in which results from the present iteration
  ! are stored. ipos=iter-1 until ipos=n_iter, then back to 1,2,...
  !
  ipos =iter-1-((iter-2)/n_iter)*n_iter
  !
  if (iter.gt.1) then
     if (saveonfile) then
        call davcio(df(1,ipos),2*ngm0*nspin,iunmix,1,-1)
        call davcio(dv(1,ipos),2*ngm0*nspin,iunmix,2,-1)
        if (lda_plus_u) then
           call davcio(df_ns(1,1,1,1,ipos),25*nspin*nat,iunmix2,1,-1)
           call davcio(dv_ns(1,1,1,1,ipos),25*nspin*nat,iunmix2,2,-1)
        end if
     end if
     call DAXPY(2*ngm0*nspin,-1.d0,rhocout,1,df(1,ipos),1)
     call DAXPY(2*ngm0*nspin,-1.d0,rhocin ,1,dv(1,ipos),1)
!        norm = sqrt(rho_dot_product(df(1,ipos),df(1,ipos)))
!        call DSCAL (2*ngm0*nspin,-1.d0/norm,df(1,ipos),1)
!        call DSCAL (2*ngm0*nspin,-1.d0/norm,dv(1,ipos),1)
     if (lda_plus_u) then
        call DAXPY(25*nspin*nat,-1.d0,nsout,1,df_ns(1,1,1,1,ipos),1)
        call DAXPY(25*nspin*nat,-1.d0,nsin ,1,dv_ns(1,1,1,1,ipos),1)
     end if
  end if
  !
  if (saveonfile) then
     do i=1,iter_used
        if (i.ne.ipos) then
           call davcio(df(1,i),2*ngm0*nspin,iunmix,2*i+1,-1)
           call davcio(dv(1,i),2*ngm0*nspin,iunmix,2*i+2,-1)
           if (lda_plus_u) then
              call davcio(df_ns(1,1,1,1,i),25*nspin*nat,iunmix2,2*i+1,-1)
              call davcio(dv_ns(1,1,1,1,i),25*nspin*nat,iunmix2,2*i+2,-1)
           end if
        end if
     end do
     call davcio(rhocout,2*ngm0*nspin,iunmix,1,1)
     call davcio(rhocin ,2*ngm0*nspin,iunmix,2,1)
     if (iter.gt.1) then
        call davcio(df(1,ipos),2*ngm0*nspin,iunmix,2*ipos+1,1)
        call davcio(dv(1,ipos),2*ngm0*nspin,iunmix,2*ipos+2,1)
     end if
     if (lda_plus_u) then
        call davcio(nsout,25*nspin*nat,iunmix2,1,1)
        call davcio(nsin ,25*nspin*nat,iunmix2,2,1)
        if (iter.gt.1) then
           call davcio(df_ns(1,1,1,1,ipos),25*nspin*nat,iunmix2,2*ipos+1,1)
           call davcio(dv_ns(1,1,1,1,ipos),25*nspin*nat,iunmix2,2*ipos+2,1)
        end if
     end if
  else
     call DCOPY(2*ngm0*nspin,rhocin ,1,rhoinsave,1)
     call DCOPY(2*ngm0*nspin,rhocout,1,rhoutsave,1)
     if (lda_plus_u) then
        call DCOPY(25*nspin*nat,nsin ,1,nsinsave ,1)
        call DCOPY(25*nspin*nat,nsout,1,nsoutsave,1)
     end if
end if
  !
  do i=1,iter_used
     do j=i,iter_used
        betamix(i,j) = rho_dot_product(df(1,j),df(1,i))
     end do
  end do
  !
  call DSYTRF ('U',iter_used,betamix,maxmix,iwork,work,maxmix,info)
  call errore('broyden','factorization',info)
  call DSYTRI ('U',iter_used,betamix,maxmix,iwork,work,info)
  call errore('broyden','DSYTRI',info)
  !
  do i=1,iter_used
     do j=i+1,iter_used
        betamix(j,i)=betamix(i,j)
     end do
  end do
  !
  do i=1,iter_used
     work(i) = rho_dot_product(df(1,i),rhocout)
  end do
  !
  do i=1,iter_used
     gamma0=0.d0
     do j=1,iter_used
        gamma0 = gamma0 + betamix(j,i)*work(j)
     end do

     call DAXPY(2*ngm0*nspin,-gamma0,dv(1,i),1,rhocin,1)
     call DAXPY(2*ngm0*nspin,-gamma0,df(1,i),1,rhocout,1)
     if (lda_plus_u) then
        call DAXPY(25*nspin*nat,-gamma0,dv_ns(1,1,1,1,i),1,nsin(1,1,1,1) ,1)
        call DAXPY(25*nspin*nat,-gamma0,df_ns(1,1,1,1,i),1,nsout(1,1,1,1),1)
     end if
  end do
  !
#ifdef DEBUG
  write (6,100) dr2, rho_dot_product(rhocout,rhocout)
  write (6,'(" dehar =",f15.8)') dehar
#endif
100  format (' dr2 =',1pe15.1, ' internal_best_dr2= ', 1pe15.1)

  ! - auxiliary vectors dv and df not needed anymore
  if (saveonfile) then
     if (lda_plus_u) then
        close(iunmix2,status='keep')
        deallocate (df_ns, dv_ns)
     end if
     close(iunmix, status='keep')
     deallocate(dv, df)
  else
     inext=iter-((iter-1)/n_iter)*n_iter
     if (lda_plus_u) then
        call DCOPY(25*nspin*nat,nsoutsave,1,df_ns(1,1,1,1,inext),1)
        call DCOPY(25*nspin*nat,nsinsave ,1,dv_ns(1,1,1,1,inext),1)
        deallocate (nsinsave, nsoutsave)
     end if
     call DCOPY(2*ngm0*nspin,rhoutsave,1,df(1,inext),1)
     call DCOPY(2*ngm0*nspin,rhoinsave ,1,dv(1,inext),1)
     deallocate(rhoutsave, rhoinsave)
  end if

  ! - preconditioning the new search direction (if imix.gt.0)

  if (imix.eq.1) then
     call approx_screening(rhocout)
  else if (imix.eq.2) then
     call approx_screening2(rhocout,rhocin)
  end if

  ! - set new trial density

  call DAXPY(2*ngm0*nspin,alphamix,rhocout,1,rhocin,1)

  do is=1,nspin
     aux(:) =(0.d0, 0.d0)
     do ig=1,ngm0
        aux(nl (ig)) = rhocin(ig,is)
        aux(nlm(ig)) = conjg(aux(nl (ig)))
     end do
     call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
     rhoin(:,is) = rhoin(:,is) + DREAL(aux(:))
  end do
  if (lda_plus_u) call DAXPY(25*nspin*nat,alphamix,nsout,1,nsin,1)

  ! - clean up

  deallocate(rhocout)
  deallocate(rhocin)
  deallocate(aux)
  call stop_clock('mix_rho')

  return
end subroutine mix_rho

!
!--------------------------------------------------------------------
function rho_dot_product (rho1,rho2)
  !--------------------------------------------------------------------
  ! this function evaluates the dot product between two input densities
  !
  use parameters, only : DP
  use pwcom
  !
  ! I/O variables
  !
  real (kind=DP) :: rho_dot_product ! (out) the function value

  complex (kind=DP) :: rho1(ngm0,nspin), rho2(ngm0,nspin) ! (in) the two densities

  !
  ! and the local variables
  !
  real (kind=DP) :: fac   ! a multiplicative factors

  integer  :: is, ig

  rho_dot_product = 0.d0

  if (nspin.eq.1) then
     is=1
     do ig = gstart,ngm0
        fac = 2.d0*e2*fpi / (tpiba2*gg(ig))
        rho_dot_product = rho_dot_product +  fac * &
                          DREAL(conjg(rho1(ig,is))*rho2(ig,is))
     end do
  else
     do ig = gstart,ngm0
        fac = 2.d0*e2*fpi / (tpiba2*gg(ig))
        rho_dot_product = rho_dot_product +  fac * &
                          DREAL(conjg(rho1(ig,1)+rho1(ig,2))* &
                                     (rho2(ig,1)+rho2(ig,2)))
     end do
     ! G=0 term
     if (gstart == 2) then
        rho_dot_product = rho_dot_product + e2*fpi / (tpi**2) * &
                          DREAL(conjg(rho1(1,1)-rho1(1,2))* &
                                     (rho2(1,1)-rho2(1,2)))
     end if
     fac = 2.d0*e2*fpi / (tpi**2)  ! lambda=1 a.u.
     do ig = gstart,ngm0
        rho_dot_product = rho_dot_product +  fac * &
                          DREAL(conjg(rho1(ig,1)-rho1(ig,2))* &
                                     (rho2(ig,1)-rho2(ig,2)))
     end do
  end if

  rho_dot_product = rho_dot_product * omega / 2.d0
#ifdef __PARA
  call reduce(1,rho_dot_product)
#endif

  return
end function rho_dot_product

!--------------------------------------------------------------------
function fn_dehar (drho)
  !--------------------------------------------------------------------
  ! this function evaluates the residual hartree energy of drho
  !
  use parameters, only : DP
  use pwcom
  !
  ! I/O variables
  !
  real (kind=DP) :: fn_dehar ! (out) the function value

  complex (kind=DP) :: drho(ngm0,nspin) ! (in) the density difference

  !
  ! and the local variables
  !
  real (kind=DP) :: fac   ! a multiplicative factors

  integer  :: is, ig

  fn_dehar = 0.d0

  if (nspin.eq.1) then
     is=1
     do ig = gstart,ngm0
        fac = 2.d0*e2*fpi / (tpiba2*gg(ig))
        fn_dehar = fn_dehar +  fac * abs(drho(ig,is))**2
     end do
  else
     do ig = gstart,ngm0
        fac = 2.d0*e2*fpi / (tpiba2*gg(ig))
        fn_dehar = fn_dehar +  fac * abs(drho(ig,1)+drho(ig,2))**2
     end do
  end if

  fn_dehar = fn_dehar * omega / 2.d0

#ifdef __PARA
  call reduce(1,fn_dehar)
#endif

  return
end function fn_dehar

!--------------------------------------------------------------------
subroutine approx_screening (drho)
  !--------------------------------------------------------------------
  ! apply an average TF preconditioning to drho
  !
  use parameters, only : DP
  use pwcom
  !
  ! I/O
  !
  complex (kind=DP) drho(ngm0,nspin) ! (in/out)
  !
  ! and the local variables
  !
  real (kind=DP) :: rrho, rmag, rs, agg0

  integer :: is, ig

  rs = (3.d0*omega/fpi/nelec)**(1.d0/3.d0)
  agg0 = (12.d0/pi)**(2.d0/3.d0)/tpiba2/rs

#ifdef DEBUG
  write (6,'(a,f12.6,a,f12.6)') ' avg rs  =', rs, ' avg rho =', nelec/omega
#endif

  if (nspin.eq.1) then
     is = 1
     do ig = 1,ngm0
        drho(ig,is) =  drho(ig,is) * gg(ig)/(gg(ig)+agg0)
     end do
  else
     do ig = 1,ngm0
        rrho = (drho(ig,1) + drho(ig,2)) * gg(ig)/(gg(ig)+agg0)
        rmag = (drho(ig,1) - drho(ig,2))
        drho(ig,1) =  0.5d0 * (rrho + rmag)
        drho(ig,2) =  0.5d0 * (rrho - rmag)
     end do
  end if

  return
end subroutine approx_screening

!
!--------------------------------------------------------------------
  subroutine approx_screening2 (drho,rhobest)
  !--------------------------------------------------------------------
  ! apply a local-density dependent TF preconditioning to drho
  !
  use parameters, only : DP
  use pwcom
  use gamma
  !
  ! I/O
  !
  !
  complex (kind=DP) ::  drho(ngm0,nspin), rhobest(ngm0,nspin)
  !
  !    and the local variables
  !
  integer :: mmx
  parameter (mmx=12)
  integer :: iwork(mmx),i,j,m,info, nspin_save
  real (kind=DP) :: rs, min_rs, max_rs, avg_rsm1, target, &
                    dr2_best, ccc, cbest, l2smooth
  real (kind=DP) :: aa(mmx,mmx), invaa(mmx,mmx), bb(mmx), work(mmx), &
                    vec(mmx),agg0
  complex (kind=DP) :: rrho, rmag

  complex (kind=DP), allocatable :: aux(:), v(:,:), w(:,:), dv(:), &
                                vbest(:), wbest(:)
  ! aux(nrxx), v(ngm0,mmx), w(ngm0,mmx), dv(ngm0), vbest(ngm0), wbest(ngm0)
  real (kind=DP), allocatable :: alpha(:)
  ! alpha(nrxx)

  integer :: is, ir, ig

  real (kind=DP) rho_dot_product
  external rho_dot_product

  if (nspin.eq.2) then
     do ig=1,ngm0
        rrho       = drho(ig,1) + drho(ig,2)
        rmag       = drho(ig,1) - drho(ig,2)
        drho(ig,1) = rrho
        drho(ig,2) = rmag
     end do
  end if

  nspin_save = nspin
  nspin = 1
  is = 1
  target = 0.d0

!  write (6,*) ' eccoci qua '

  if (gg(1).lt.1.d-8) drho(1,is) = (0.d0,0.d0)

  allocate (alpha(nrxx), aux(nrxx), v(ngm0,mmx), w(ngm0,mmx), &
            dv(ngm0), vbest(ngm0), wbest(ngm0))

  v(:,:) = (0.d0, 0.d0)
  w(:,:) = (0.d0, 0.d0)
  dv(:)  = (0.d0, 0.d0)
  vbest(:) = (0.d0, 0.d0)
  wbest(:) = (0.d0, 0.d0)

  !
  ! - calculate alpha from density smoothed with a lambda=0 a.u.
  !
  l2smooth = 0.d0
  aux(:) = (0.d0, 0.d0)
  if (nspin.eq.1) then
     do ig=1,ngm0
        aux(nl(ig)) = rhobest(ig,1) * exp(-0.5*l2smooth*tpiba2*gg(ig))
        aux(nlm(ig)) = conjg(aux(nl (ig)))
     end do
  else
     do ig=1,ngm0
        aux(nl(ig)) =(rhobest(ig,1) + rhobest(ig,2)) &
                                    * exp(-0.5*l2smooth*tpiba2*gg(ig))
        aux(nlm(ig)) = conjg(aux(nl (ig)))
     end do
  end if
  call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
  alpha(:) = DREAL(aux(:))

  min_rs = (3.d0*omega/fpi/nelec)**(1.d0/3.d0)
  max_rs = min_rs
  avg_rsm1 = 0.d0

  do ir=1,nrxx

     alpha(ir)=abs(alpha(ir))
     rs = (3.d0/fpi/alpha(ir))**(1.d0/3.d0)

     min_rs = min(min_rs,rs)
     avg_rsm1 =avg_rsm1 + 1.d0/rs
     max_rs = max(max_rs,rs)

     alpha(ir) = rs

  end do

#ifdef __PARA
  call reduce  (1, avg_rsm1)
  call extreme (min_rs, -1)
  call extreme (max_rs, +1)
#endif

  call DSCAL(nrxx, 3.d0*(tpi/3.d0)**(5.d0/3.d0), alpha, 1)

  avg_rsm1 = (nr1*nr2*nr3)/avg_rsm1
  rs = (3.d0*omega/fpi/nelec)**(1.d0/3.d0)
  agg0 = (12.d0/pi)**(2.d0/3.d0)/tpiba2/avg_rsm1
#ifdef DEBUG
  write (6,'(a,5f12.6)') ' min/avgm1/max rs  =', min_rs,avg_rsm1,max_rs,rs
#endif

  !
  ! - calculate deltaV and the first correction vector
  !
  aux(:) =(0.d0, 0.d0)
  do ig=1,ngm0
     aux(nl(ig)) = drho(ig,is)
     aux(nlm(ig)) = conjg(aux(nl (ig)))
  end do
  call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
  aux(:) = aux(:) * alpha(:)
  call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
  do ig=1,ngm0
     dv(ig) = aux(nl(ig))*gg(ig)*tpiba2
     v(ig,1)= aux(nl(ig))*gg(ig)/(gg(ig)+agg0)
  end do
  m=1
  ccc = rho_dot_product(dv,dv)
  aa(:,:) = 0.d0
  bb(:) =0.d0

3 continue
  !
  ! - generate the vector w
  !
  do ig=1,ngm0
     w(ig,m) = gg(ig)*tpiba2*v(ig,m)
  end do
  aux(:) =(0.d0, 0.d0)
  do ig=1,ngm0
     aux(nl(ig)) = v(ig,m)
     aux(nlm(ig)) = conjg(aux(nl (ig)))
  end do
  call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
  aux(:) = aux(:)*fpi*e2/alpha(:)
  call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
  do ig=1,ngm0
     w(ig,m) = w(ig,m) + aux(nl(ig))
  end do

  !
  ! - build the linear system
  !
  do i=1,m
     aa(i,m) = rho_dot_product(w(1,i),w(1,m))
     aa(m,i) = aa(i,m)
  end do
  bb(m) = rho_dot_product(w(1,m),dv)

  !
  ! - solve it -> vec
  !
  call DCOPY (mmx*mmx,aa,1,invaa,1)
  call DSYTRF ('U',m,invaa,mmx,iwork,work,mmx,info)
  call errore('BROYDEN','factorization',info)
  call DSYTRI ('U',m,invaa,mmx,iwork,work,info)
  call errore('broyden','DSYTRI',info)
  !
  do i=1,m
     do j=i+1,m
        invaa(j,i)=invaa(i,j)
     end do
  end do
  do i=1,m
     vec(i) = 0.d0
     do j=1,m
        vec(i) = vec(i) + invaa(i,j)*bb(j)
     end do
  end do
  ! -
  vbest(:) = (0.d0,0.d0)
  wbest(:) = dv(:)
  do i=1,m
     call DAXPY(2*ngm0, vec(i), v(1,i),1, vbest,1)
     call DAXPY(2*ngm0,-vec(i), w(1,i),1, wbest,1)
  end do

  cbest = ccc
  do i=1,m
     cbest = cbest - bb(i)*vec(i)
  end do

  dr2_best= rho_dot_product(wbest,wbest)
  if (target.eq.0.d0) target = 1.d-6 * dr2_best
!  write (6,*) m, dr2_best, cbest

  if (dr2_best .lt. target) then
!     write(6,*) ' last', dr2_best/target * 1.d-6
     aux(:) = (0.d0, 0.d0)
     do ig=1,ngm0
        aux(nl(ig)) = vbest(ig)
        aux(nlm(ig)) = conjg(aux(nl (ig)))
     end do
     call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
     aux(:) = aux(:)/alpha(:)
     call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
     do ig=1,ngm0
        drho(ig,is) = aux(nl(ig))
     end do
     nspin = nspin_save
     if (nspin.eq.2) then
        do ig=1,ngm0
           rrho = drho(ig,1)
           rmag = drho(ig,2)
           drho(ig,1) = 0.5d0 * ( rrho + rmag )
           drho(ig,2) = 0.5d0 * ( rrho - rmag )
        end do
     end if
     deallocate (alpha, aux, v, w, dv, vbest, wbest)
     return
  else if (m.ge.mmx) then
!     write (6,*) m, dr2_best, cbest
     m=1
     do ig=1,ngm0
        v(ig,m)=vbest(ig)
     end do
     aa(:,:) = (0.d0, 0.d0)
     bb(:) = 0.d0
     go to 3
  end if

  m = m + 1
  v(:,m)=wbest(:)/(gg(:)+agg0)

  go to 3

end subroutine approx_screening2
