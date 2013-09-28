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
  use kinds, only: dp
  use io_global, only : stdout, ionode, ionode_id
  use mp,        only : mp_bcast
  use mp_world,  only : world_comm
  use radial_grids, only: ndmx
  use ld1_parameters, only: nwfsx
  use ld1inc, only: grid, nlcc, vxt, lsd, vpstot, vpsloc, file_screen, &
                    vh, enne, rhoc, latt, rhos, enl, &
                    nbeta, bmat, qvan, qvanl, jjs, lls, ikk, pseudotype, &
                    nwfts, enlts, octs, llts, jjts, phits, nstoaets, lpaw, &
                    which_augfun
  implicit none

  integer ::  &
       ns,    &  ! counter on pseudo functions
       ns1,   &  ! counter on pseudo functions
       ib,jb, &  ! counter on beta functions
       lam       ! the angular momentum

  real(DP) :: &
       vaux(ndmx,2)     ! work space

  real(DP), external :: int_0_inf_dr ! the integral function

  real(DP), parameter :: &
       thresh= 1.e-12_dp          ! threshold for selfconsistency

  integer  :: &
       n, nst, iwork(nwfsx), ios, nerr
  !
  !     descreening the local potential: NB: this descreening is done with
  !     the occupation of the test configuration. This is required
  !     for pseudopotentials with semicore states. In the other cases
  !     a test configuration equal to the one used for pseudopotential
  !     generation is strongly suggested
  !
  do n=1,nwfts
     enlts(n)=enl(nstoaets(n))
  enddo
  !
  !    compute the pseudowavefunctions in the test configuration
  !
  call ascheqps_drv(vpsloc, 1, thresh, .false., nerr)
  !
  !    descreening the D coefficients
  !
  if (pseudotype.eq.3) then
     do ib=1,nbeta
        do jb=1,ib
           if (lls(ib).eq.lls(jb).and.abs(jjs(ib)-jjs(jb)).lt.1.e-7_dp) then
              lam=lls(ib)
              nst=(lam+1)*2
              IF (which_augfun=='PSQ') then
                 do n=1,ikk(ib)
                    vaux(n,1)=qvanl(n,ib,jb,0)*vpsloc(n)
                 enddo
              ELSE
                 do n=1,ikk(ib)
                    vaux(n,1)=qvan(n,ib,jb)*vpsloc(n)
                 enddo
              ENDIF
              bmat(ib,jb)= bmat(ib,jb)  &
                   - int_0_inf_dr(vaux(1,1),grid,ikk(ib),nst)
           endif
           bmat(jb,ib)=bmat(ib,jb)
        enddo
     enddo
     write(stdout,'(/5x,'' The ddd matrix'')')
     do ns1=1,nbeta
        write(stdout,'(6f12.5)') (bmat(ns1,ns),ns=1,nbeta)
     enddo
  endif
  !
  !    descreening the local pseudopotential
  !
  iwork=1
  call chargeps(rhos,phits,nwfts,llts,jjts,octs,iwork)

  call new_potential(ndmx,grid%mesh,grid,0.0_dp,vxt,lsd,nlcc,latt,enne,&
       rhoc,rhos,vh,vaux,1)

  do n=1,grid%mesh
     vpstot(n,1)=vpsloc(n)
     vpsloc(n)=vpsloc(n)-vaux(n,1)
  enddo

  if (file_screen .ne.' ') then
     if (ionode) &
        open(unit=20,file=file_screen, status='unknown', iostat=ios, err=100 )
100  call mp_bcast(ios, ionode_id, world_comm)
     call errore('descreening','opening file'//file_screen,abs(ios))
     if (ionode) then
        do n=1,grid%mesh
           write(20,'(i5,7e12.4)') n,grid%r(n), vpsloc(n)+vaux(n,1), &
                vpsloc(n), vaux(n,1), rhos(n,1)
        enddo
        close(20)
     endif
  endif

  return
end subroutine descreening
