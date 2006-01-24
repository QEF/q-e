!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------------
subroutine descreening
  !--------------------------------------------------------------------------
  !
  !     This routine descreens the local potential and the ddd
  !     coefficients (the latter only in the US case)
  !     The charge density is computed with the test configuration,
  !     not the one used to generate the pseudopotential
  !      
  use ld1inc
  implicit none


  integer ::  &
       ns,    &  ! counter on pseudo functions
       ns1,   &  ! counter on pseudo functions
       ib,jb, &  ! counter on beta functions
       lam,   &  ! the angular momentum
       ind

  real(DP) :: &
       vaux(ndm,2),&   ! work space
       phist(ndm,nwfsx)! auxiliary to save the phi

  real(DP), external :: int_0_inf_dr ! the integral function

  real(DP), parameter :: &
       thresh= 1.e-12_dp          ! threshold for selfconsistency

  integer  :: &
       m, n, l, n1, n2, nwf0, nst, ikl, imax, iwork(nwfsx), &
       is, nbf, nc, ios
  !
  !     descreening the local potential: NB: this descreening is done with
  !     the occupation of the test configuration. This is required
  !     for pseudopotentials with semicore states. In the other cases
  !     a test configuration equal to the one used for pseudopotential
  !     generation is strongly suggested
  !
  nc=1
  nwfts=nwftsc(nc)
  do n=1,nwfts
     nnts(n)=nntsc(n,nc)
     llts(n)=lltsc(n,nc)
     elts(n)=eltsc(n,nc)
     !         rcutts(n)=rcut(n)
     !         rcutusts(n)=rcutus(n)
     jjts(n) = jjtsc(n,nc)
     iswts(n)=iswtsc(n,nc)
     octs(n)=octsc(n,nc)
     nstoae(n)=nstoaec(n,nc)
     enlts(n)=enl(nstoae(n))
     new(n)=.false.
  enddo

  do ns=1,nwfs
     do n=1,mesh
        phist(n,ns)=phis(n,ns)
     enddo
  enddo
  !
  !    compute the pseudowavefunctions in the test configuration
  !
  if (pseudotype.eq.1) then
     nbf=0
  else
     nbf=nbeta
  endif

  do ns=1,nwfts
     if (octs(ns).gt.0.0_dp) then
        is=iswts(ns)
        if (pseudotype ==1) then
           if ( rel < 2 .or. llts(ns) == 0 .or. &
                abs(jjts(ns)-llts(ns)+0.5_dp) < 0.001_dp) then
              ind=1
           else if ( rel == 2 .and. llts(ns) > 0 .and. &
                abs(jjts(ns)-llts(ns)-0.5_dp) < 0.001_dp) then
              ind=2
           endif
           do n=1,mesh
              vaux(n,1)=vpsloc(n)+vnl(n,llts(ns),ind)
           enddo
        else
           do n=1,mesh
              vaux(n,1)=vpsloc(n)
           enddo
        endif
        call ascheqps(nnts(ns),llts(ns),jjts(ns),enlts(ns),    &
             mesh,ndm,dx,r,r2,sqr,vaux,thresh,phis(1,ns), & 
             betas,bmat,qq,nbf,nwfsx,lls,jjs,ikk)
        !            write(6,*) ns, nnts(ns),llts(ns), jjts(ns), enlts(ns)
     endif
  enddo
  !
  !    descreening the D coefficients
  !
  if (pseudotype.eq.3) then
     do ib=1,nbeta
        do jb=1,ib
           if (lls(ib).eq.lls(jb).and.abs(jjs(ib)-jjs(jb)).lt.1.e-7_dp) then
              lam=lls(ns)
              nst=(lam+1)*2
              do n=1,ikk(ib)
                 vaux(n,1)=qvan(n,ib,jb)*vpsloc(n)
              enddo
              bmat(ib,jb)= bmat(ib,jb)  &
                   - int_0_inf_dr(vaux(1,1),r,r2,dx,ikk(ib),nst)
           endif
           bmat(jb,ib)=bmat(ib,jb)
        enddo
     enddo
     write(6,'(/5x,'' The ddd matrix'')')
     do ns1=1,nbeta
        write(6,'(6f12.5)') (bmat(ns1,ns),ns=1,nbeta)
     enddo
  endif
  !
  !    descreening the local pseudopotential
  !
  iwork=1
  call normalize
  call chargeps(nwfts,llts,jjts,octs,iwork)

  call new_potential(ndm,mesh,r,r2,sqr,dx,0.0_dp,vxt,lsd,nlcc,latt,enne, &
       rhoc,rhos,vh,vaux)

  do n=1,mesh
     vpstot(n,1)=vpsloc(n)
     vpsloc(n)=vpsloc(n)-vaux(n,1)
  enddo

  if (file_screen .ne.' ') then
     open(unit=20,file=file_screen, status='unknown', iostat=ios, &
          err=100 )
100  call errore('descreening','opening file'//file_screen,abs(ios))
     do n=1,mesh
        write(20,'(i5,7e12.4)') n,r(n), vpsloc(n)+vaux(n,1), vpsloc(n), &
             vaux(n,1),   rhos(n,1)
     enddo
     close(20)
  endif
  !
  !  copy the phis used to construct the pseudopotential
  !
  do ns=1,nwfs
     do n=1,mesh
        phis(n,ns)=phist(n,ns)
     enddo
  enddo

  return
end subroutine descreening
