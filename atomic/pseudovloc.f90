!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------------
subroutine pseudovloc
  !--------------------------------------------------------------------------
  !
  !     This routine generate a local pseudopotential 
  !     The output of the routine are:
  !     vpsloc: the local pseudopotential
  !
  !      
  use ld1inc
  implicit none

  real(DP) :: &
       fae,   &  ! the value of the all-electron function
       f1ae,  &  ! its first derivative
       f2ae,  &  ! the second derivative
       faenorm   ! the norm of the function

  integer :: &
       nwf0, &  ! used to specify the all electron function
       nst,  &  ! auxiliary
       iok,  &  ! if 0 there are no problems
       ik       ! the point corresponding to rc

  real(DP) ::             &
       f1aep1,f1aem1,jnor, &  ! auxilairy quantities
       p1aep1, p1aem1,     &  ! derivatives of the bessel functions
       xc(8),              &  ! the coefficients of the fit
       bm(2),              &  ! the derivative of the bessel
       vaux(ndm,2),        &  ! keeps the potential
       j1(ndm,4)              ! the bessel functions

  real(DP) :: &
       deriv_7pts, deriv2_7pts


  integer ::         &
       n,        &  ! counter on mesh points
       indi,rep, &  ! auxiliary
       indns(0:1), & ! auxiliary
       nc           ! counter on bessel

  if (lloc < 0) then
     !
     !   Compute the potential by smoothing the AE potential
     !
     !   Compute the ik which correspond to this cutoff radius
     !
     write(6,"(/,5x,' Generating local potential from pseudized AE ', &
          & 'potential, matching radius rcloc = ',f8.4)") rcloc
     ik=0
     do n=1,mesh
        if (r(n) < rcloc) ik=n
     enddo
     if (mod(ik,2) == 0) ik=ik+1
     if (ik <= 1 .or. ik > mesh) &
          call errore('pseudovloc','wrong matching point',1)
     !
     !    compute first and second derivative
     !
     fae=vpot(ik,1)
     f1ae=deriv_7pts(vpot,ik,r(ik),dx)
     f2ae=deriv2_7pts(vpot,ik,r(ik),dx)
     !
     !    find the q_i of the bessel functions
     !      
     call find_qi(f1ae/fae,xc(3),ik,0,2,0,iok)
     if (iok /= 0) &
          call errore('pseudovloc','problems with find_qi',1)
     !
     !    compute the functions
     !
     do nc=1,2
        call sph_bes(ik+1,r,xc(2+nc),0,j1(1,nc))
        jnor=j1(ik,nc)
        do n=1,ik+1
           j1(n,nc)=j1(n,nc)*vpot(ik,1)/jnor
        enddo
     enddo
     !
     !    compute the second derivative and impose continuity of zero, 
     !    first and second derivative
     !

     do nc=1,2
        p1aep1=(j1(ik+1,nc)-j1(ik,nc))/(r(ik+1)-r(ik))
        p1aem1=(j1(ik,nc)-j1(ik-1,nc))/(r(ik)-r(ik-1))
        bm(nc)=(p1aep1-p1aem1)*2.0_dp/(r(ik+1)-r(ik-1))
     enddo

     xc(2)=(f2ae-bm(1))/(bm(2)-bm(1))
     xc(1)=1.0_dp-xc(2)
     write(6, 110) rcloc,xc(4)**2 
110  format (/5x, ' Local pseudo, rcloc=',f6.3, &
          ' Estimated cut-off energy= ', f8.2,' Ry')
     !
     !    define the local pseudopotential
     !
     do n=1,ik
        vpsloc(n)=xc(1)*j1(n,1)+xc(2)*j1(n,2)
     enddo

     do n=ik+1,mesh
        vpsloc(n)=vpot(n,1)
     enddo

  else
     !
     !    if a given angular momentum gives the local component this is done 
     !    here
     !
     nst=(lloc+1)*2
     if (rel==2 .and. lloc > 0) then
        rep=1
        indns(0)=nsloc
        indns(1)=nsloc+1
        if (jjs(nsloc) > jjs(nsloc+1) ) then
           indns(0)=nsloc+1
           indns(1)=nsloc
        endif
     else
        rep=0
        indns(0)=nsloc
     endif
     vpsloc=0.0_dp
     vaux=0.0_dp
     do indi=0,rep
        nwf0=nstoae(nsloc+indi)
        if (enls(nsloc+indi) == 0.0_dp) &
             enls(nsloc+indi)=enl(nstoae(nsloc+indi))
        !
        !    compute the ik closer to r_cut
        !
        ik=0
        do n=1,mesh
           if (r(n) < rcut(nsloc+indi)) ik=n
        enddo
        if (mod(ik,2).eq.0) ik=ik+1
        if (ik <= 1 .or. ik > mesh) &
           call errore('pseudovloc','wrong matching point',1)
        rcloc=rcut(nsloc+indi)
        if (rep == 0) then
           write(6,"(/,5x,' Generating local pot.: lloc=',i1, &
                  & ', matching radius rcloc = ',f8.4)") lloc, rcloc
        else
           if (rel==2) then
              write(6,"(/,5x,' Generating local pot.: lloc=',i1, &
                  &', j=',f5.2,', matching radius rcloc = ',f8.4)") &
                  lloc, lloc-0.5+indi, rcloc
           else
              write(6,"(/,5x,' Generating local pot.: lloc=',i1, &
                  &', spin=',i1,', matching radius rcloc = ',f8.4)") &
                  lloc, indi+1, rcloc
           endif
        endif
        !
        !   compute the phi functions
        !
        call compute_phipot(lloc,ik,nwf0,indns(indi),xc)
        !
        !     set the local potential equal to the all-electron one at large r
        !
        do n=1,mesh
           if (r(n) > rcloc) then
              vaux(n,indi+1)=vpot(n,1)
           else
              vaux(n,indi+1)=chis(n,indns(indi))/phis(n,indns(indi))
           endif
        enddo
     enddo
     if (rep==0) then
        do n=1,mesh
           vpsloc(n)=vaux(n,1)
        enddo
     else
        do n=1,mesh
           vpsloc(n)=(lloc*vaux(n,1)+(lloc+1.0_dp)*vaux(n,2))/ &
                (2.0_dp*lloc+1.0_dp)
        enddo
     endif
  endif

  return
end subroutine pseudovloc
