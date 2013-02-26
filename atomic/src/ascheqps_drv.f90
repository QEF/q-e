!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------------
subroutine ascheqps_drv(veff, ncom, thresh, flag_all, nerr)
  !--------------------------------------------------------------------------
  !
  !     This routine is a driver that calculates for the test
  !     configuration the solutions of the Kohn and Sham equation
  !     with a fixed pseudo-potential. The potentials are assumed
  !     to be screened. The effective potential veff is given in input.
  !     The output wavefunctions are written in phits and are normalized.
  !     If flag is .true. compute all wavefunctions, otherwise only
  !     the wavefunctions with positive occupation.
  !      
  use kinds, only: dp
  use ld1_parameters, only: nwfsx
  use radial_grids, only: ndmx
  use ld1inc, only: grid, pseudotype, rel, &
                    lls, jjs, qq, ikk, ddd, betas, nbeta, vnl, &
                    nwfts, iswts, octs, llts, jjts, nnts, enlts, phits 
  implicit none

  integer ::    &
          nerr, &     ! control the errors of the routine ascheqps
          ncom        ! number of components of the pseudopotential

  real(DP) :: &
       veff(ndmx,ncom)    ! work space for writing the potential 

  logical :: flag_all    ! if true calculates all the wavefunctions

  integer ::  &
       ns,    &  ! counter on pseudo functions
       is,    &  ! counter on spin
       nbf,   &  ! auxiliary nbeta
       n,     &  ! index on r point
       nstop, &  ! errors in each wavefunction
       ind

  real(DP) :: &
       vaux(ndmx,2)     ! work space for writing the potential 

  real(DP) :: thresh         ! threshold for selfconsistency
  !
  !    compute the pseudowavefunctions in the test configuration
  !
  if (pseudotype.eq.1) then
     nbf=0
  else
     nbf=nbeta
  endif

  nerr=0
  do ns=1,nwfts
     if ( octs(ns) > 0.0_dp .or. ( octs(ns) > -1.0_dp .and. flag_all ) ) then
        is=iswts(ns)
        if (ncom==1.and.is==2) call errore('ascheqps_drv','incompatible spin',1)
        if (pseudotype == 1) then
           if ( rel < 2 .or. llts(ns) == 0 .or. &
                abs(jjts(ns)-llts(ns)+0.5_dp) < 0.001_dp) then
              ind=1
           else if ( rel == 2 .and. llts(ns) > 0 .and. &
                abs(jjts(ns)-llts(ns)-0.5_dp) < 0.001_dp) then
              ind=2
           else
              call errore('ascheqps_drv','unexpected case',1)
           endif
           do n=1,grid%mesh
              vaux(n,is)=veff(n,is)+vnl(n,llts(ns),ind)
           enddo
        else
           do n=1,grid%mesh
              vaux(n,is)=veff(n,is)
           enddo
        endif
        call ascheqps(nnts(ns),llts(ns),jjts(ns),enlts(ns),grid%mesh,ndmx,&
             grid,vaux(1,is),thresh,phits(1,ns),betas,ddd(1,1,is),qq,nbf, &
             nwfsx,lls,jjs,ikk,nstop)
        !           write(6,*) ns, nnts(ns),llts(ns), jjts(ns), enlts(ns)
        !
        !   normalize the wavefunctions 
        !
        call normalize(phits(1,ns),llts(ns),jjts(ns), ns)
        !
        !   not sure whether the "best" error code should be like this:
        ! IF ( octs(ns) > 0.0_dp ) nerr = nerr + nstop
        !   i.e. only for occupied states, or like this:
        nerr = nerr + nstop
     endif
  enddo

  return
end subroutine ascheqps_drv
