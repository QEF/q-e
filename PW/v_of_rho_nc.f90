
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#include "f_defs.h"
!--------------------------------------------------------------------
subroutine v_of_rho_nc (rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, &
     nrx3, nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
     ehart, etxc, vtxc, charge, v,lambda, vtcon, i_cons, mcons, &
          pointlist, pointnum, factlist, nat, ntyp, ityp)
  !--------------------------------------------------------------------
  !
  !     This routine computes the Hartree and Exchange and Correlation
  !     potential and energies which corresponds to a given charge density
  !     The XC potential is computed in real space, while the
  !     Hartree potential is computed in reciprocal space.
  !
  !
  USE kinds, only: DP
  USE noncollin_module, ONLY : magtot_nc, bfield
  implicit none
  !
  !    first the dummy variables
  !

  integer :: nat, ntyp, nspin, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
       ngm, nl (ngm),gstart, i_cons, pointlist(nrxx,nat), pointnum(nat), &
       ityp(nat)
  ! input: number of atoms
  ! input: 1=lda, 2=lsda, 4=noncolin
  ! input: the FFT indices
  ! input: the true array dimensions
  ! input: the total dimension
  ! input: the number of G vectors
  ! input: correspondence G <-> FFT
  ! input: first nonzero G-vector
  ! local integrated magnetization for every atom
  ! local integrated charge for every atom
  ! input: flag indicating if constrains on local magn. is calcul
  ! input: number of sch points
  real(kind=DP) :: rho (nrxx, nspin), rho_core (nrxx), g (3, ngm), &
       gg (ngm), alat, omega, vtxc, etxc, ehart, charge, v (nrxx, nspin), &
       vtcon, mcons(3,ntyp), factlist(nrxx,nat),lambda,fact,ma,xx,xxx,m1(3), &
       m_loc(3,nat), r_loc(nat)
  ! input: the valence charge
  ! input: the core charge
  ! input: the G vectors
  ! input: the norm of G vector
  ! input: the length of the cell
  ! input: the volume of the cell
  ! output: the integral V_xc * rho
  ! output: the E_xc energy
  ! output: the hartree energy
  ! output: the integral of the charge
  ! output: the H+xc_up  potential
  ! output: constrain-contribution to the total energy
  ! input: constrained value of the local moments
  ! input: weightening factors for local integrations
  ! input: prefactor for penalty-constrain-calculations
  ! auxiliary variables, used to calculate constraining fields
  ! (fact,ma,xx,xxx,m1(3))
  ! input: points in the integration volume around atom na
  !
  real(kind=DP), parameter  :: fpi = 4.d0 * 3.14159265358979d0, &
                               e2  = 2.d0
  !
  !    and the local variables
  !

  real(kind=DP) :: tpiba2, fac
  ! the measure unit in reciprocal space
  ! a multiplicative factors
  real(kind=DP), allocatable ::  aux (:,:), aux1 (:,:)
  ! used to do the fft
  ! auxiliary variable for the potential

  integer :: ir, is, ig, ipol ,nt, na
  ! counters

  call start_clock ('v_of_rho')
  !
  !  calculate exchange-correlation potential
  !
  call v_xc_nc (rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
       nl, ngm, g, nspin, alat, omega, etxc, vtxc, v)
! If noncollinear is set, one can calculate constrains on the local
! magnetization, calculated in get_locals
! To this end, a "penalty term" of the form
! E_p = lambda * ( m_loc - m_loc_constr)^2
! is added to the energy. Here we calculate the resulting
! "constraining B-field" and add it to v(ir,2..4)

! vtcon is the contribution of the constraining field to the total
! energy 

  vtcon = 0.d0
  if (i_cons.ne.0) then

! get the actual values of the local integrated quantities
     if (i_cons.lt.3) then
        call get_locals(r_loc, m_loc)

        do na = 1,ntyp
           nt = ityp(na)

          ! i_cons = 1 means that the 3 components of the magnetization
          ! are constrained, they are given in the input-file
           if (i_cons.eq.1) then
              do ipol = 1,3
                 m1(ipol) = m_loc(ipol,nt) - mcons(ipol,nt)
              enddo
           else if (i_cons.eq.2) then
           ! i_cons = 2 means that the angle theta between the local
           ! magn. moment and the z-axis is constrained
           ! mcons (1,nt) is the cos of the constraining angle theta
           ! the penalty functional in this case is
           ! lambda*(acos(m_loc(z)/|m_loc|) - theta )^2
               ma = dsqrt(m_loc(1,nt)**2+m_loc(2,nt)**2+m_loc(3,nt)**2)
               xx = m_loc(3,nt)/ma
               if (abs(xx - 1.d0).gt.1.d-10) then
                  xxx = - (acos(xx) - mcons(1,nt))/sqrt(1.d0 - xx*xx)
               else
                  xxx = -1.d0
               endif

               m1(1) = - xxx * m_loc(1,nt)*m_loc(3,nt) / (ma*ma*ma)
               m1(2) = - xxx * m_loc(2,nt)*m_loc(3,nt) / (ma*ma*ma)
               m1(3) =   xxx * (ma*ma-m_loc(3,nt)*m_loc(3,nt)) / (ma*ma*ma)

            endif

            do ir =1,pointnum(na)
               fact = 2.d0*lambda*factlist(ir,na)*omega/(nr1*nr2*nr3)
               do ipol = 1,3
                  v(pointlist(ir,na),ipol+1) = v(pointlist(ir,na) &
                         ,ipol+1) +fact*m1(ipol)
!                  vtcon = vtcon + fact*m1(ipol)*rho(pointlist(ir,na),ipol+1)
               enddo            ! ipol
            enddo               ! points
        enddo ! na
!        call reduce(1,vtcon)
!        vtcon = omega*vtcon/(nr1*nr2*nr3)
     else if (i_cons==3) then
        fact = 2.d0*lambda
        do ipol=1,3
           bfield(ipol)=-fact*(magtot_nc(ipol)-mcons(ipol,1))
        enddo
        do ipol = 2,4
           fact=bfield(ipol-1)
           do ir =1,nrxx
              v(ir,ipol) = v(ir,ipol)-fact
           enddo            
        enddo              
     else if (i_cons==4) then
        bfield(:)=mcons(:,1)
        do ipol = 2,4
           fact=bfield(ipol-1)
           do ir =1,nrxx
              v(ir,ipol) = v(ir,ipol)-fact
           enddo            
        enddo              
     else
        call errore('v_of_rho_nc','i_cons not programmed',1)
     endif
  endif  ! end penalty functional calculation


  allocate (aux(2,nrxx),aux1(2,ngm) )
  tpiba2 = (fpi / 2.d0 / alat) **2
  !
  !  copy total rho in aux
  !
  aux(2,:) = 0.d0
  aux(1,:) = rho(:,1)
  if (nspin.eq.2) aux(1,:) = aux(1,:) + rho(:,2)
  !
  !  bring rho (aux) to G space
  !

  call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
  charge = 0.d0
  if (gg(1) .lt.1.0d-8) charge = omega * aux(1,nl(1))
#ifdef __PARA
  call reduce (1, charge)
#endif
  !
  !      calculate hartree potential in G-space (NB: V(G=0)=0 )
  !
  ehart = 0.d0
  aux1(:,:) = 0.d0
  do ig = gstart, ngm
     fac = e2 * fpi / (tpiba2 * gg (ig) )
     ehart = ehart + (aux(1,nl(ig))**2 + aux(2,nl(ig))**2) * fac
     aux1(1,ig) = fac * aux(1,nl(ig))
     aux1(2,ig) = fac * aux(2,nl(ig))
  enddo
  ehart = ehart * omega / 2.d0
#ifdef __PARA

  call reduce (1, ehart)
#endif
  aux(:,:) = 0.d0
  do ig = 1, ngm
     aux(1,nl(ig)) = aux1(1,ig)
     aux(2,nl(ig)) = aux1(2,ig)
  enddo
  !
  !      transform hartree potential to real space
  !
  call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
  !
  !      add hartree potential to the xc potential
  !
  if ( nspin == 4 ) then
     !
     ! Noncollinear case: Add the Hartree potential only to the first
     ! component of the potential
     do ir=1,nrxx
        v(ir,1)=v(ir,1)+aux(1,ir)
     end do
!
! LDA and LSDA case

  else

     do is = 1, nspin
        do ir = 1, nrxx
          v(ir,is) = v(ir,is) + aux(1,ir)
        enddo
     enddo
  endif

  deallocate (aux,aux1)

  call stop_clock ('v_of_rho')

  return
end subroutine v_of_rho_nc
!
!--------------------------------------------------------------------
subroutine v_xc_nc (rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
     nrxx, nl, ngm, g, nspin, alat, omega, etxc, vtxc, v)
  !--------------------------------------------------------------------
  !
  !     Exchange-Correlation potential Vxc(r) from n(r)
  !
  USE io_global,  ONLY : stdout
  USE spin_orb,   ONLY : domag
  USE kinds, only : DP
  implicit none
  !
  ! input
  !

  integer :: nspin, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm, &
       nl (ngm), i
  !  nspin=1 :unpolarized, =2 :spin-polarized, nspin=4 :noncollinear
  ! the FFT indices
  ! the true FFT array dimensions
  ! the total dimension
  ! the number of G vectors
  ! correspondence G <-> FFT
  real(kind=DP) :: rho (nrxx, nspin), rho_core (nrxx), g (3, ngm), &
       alat, omega
  ! the valence charge
  ! the core charge
  ! the G vectors
  ! the length of the cell
  ! the volume of the cell
  !
  ! output
  !
  real(kind=DP) :: v (nrxx, nspin), vtxc, etxc
  ! V_xc potential
  ! integral V_xc * rho
  ! E_xc energy
  !
  !    local variables
  !
  ! the square of the e charge
  real(kind=DP) :: e2
  parameter (e2 = 2.d0)

  real(kind=DP) :: rhox, arhox, zeta, ex, ec, vx (2), vc (2), vs, amag
  ! the total charge in each point
  ! the absolute value of the charge
  ! the absolute value of the charge
  ! local exchange energy
  ! local correlation energy
  ! local exchange potential
  ! local correlation potential
  ! ??
  ! absolute magnetization 
  integer :: ir, is, ig, neg (3), ipol
  ! counter on mesh points
  ! counter on spin polarizations
  ! counter on G vectors
  ! number of points with wrong zeta/cha
  ! counter of wfc coordinates
  !
  !
  !      call start_clock('vxc')
  !
  ! initialization
  !
  etxc = 0.d0
  vtxc = 0.d0
  v(:,:) = 0.d0
  !
  ! case of noncollinear magnetism
  !
  neg=0
  if (domag) then
     do ir = 1,nrxx
        amag = sqrt(rho(ir,2)**2+rho(ir,3)**2+rho(ir,4)**2)
        rhox = rho(ir,1) + rho_core(ir)
        arhox = abs(rhox)
        if (arhox.gt.1.d-30) then
           zeta = amag / arhox
           if(abs(zeta).gt.1.d0) then
               neg(3)=neg(3)+1
               zeta=sign(1.d0,zeta)
           endif
           call xc_spin(arhox,zeta,ex,ec,vx(1),vx(2),vc(1),vc(2))
           vs=0.5d0*(vx(1)+vc(1)-vx(2)-vc(2))
           v(ir,1) = e2*(0.5d0*(vx(1)+vc(1)+vx(2)+vc(2)))
           !               WRITE( stdout,'(4f18.8)') arhox,zeta,vs,v(ir,1)
           if (amag.gt.1.d-20) then
              do ipol=2,4
                 v(ir,ipol) = e2*vs*rho(ir,ipol)/amag
              enddo
           endif
           etxc= etxc + e2*(ex+ec)*rhox
           vtxc= vtxc + v(ir,1)*rho(ir,1)
        endif
     enddo
  else
     do ir = 1, nrxx
        rhox = rho (ir, 1) + rho_core (ir)
        arhox = abs (rhox)
        if (arhox > 1.d-30) then
           call xc (arhox, ex, ec, vx, vc)
           v(ir,1) = e2 * (vx(1) + vc(1) )
           etxc = etxc + e2 * (ex + ec) * rhox
           vtxc = vtxc + v(ir,1) * rho(ir,1)
        endif
     enddo
  endif

#ifdef __PARA
  call ireduce (3, neg)
#endif

  if (neg(3).gt.0) WRITE( stdout,'(/,4x," npt with |zeta| > 1: ",i8, &
       &  ", npt tot ",i8, ",",f10.2, " %" )') neg(3), &
       &  nr1*nr2*nr3, float(neg(3)*100) / real(nr1*nr2*nr3)
  if (neg(1).gt.0) WRITE( stdout,'(/,4x," npt with rhoup < 0: ",i8, &
       &  ", npt tot ",i8, ",",f10.2, " %" )') neg(1), &
       &  nr1*nr2*nr3, float(neg(1)*100) / real(nr1*nr2*nr3)
  if (neg(2).gt.0) WRITE( stdout,'(/,4x," npt with rhodw < 0: ",i8, &
       &  ", npt tot ",i8, ",",f10.2, " %" )') neg(2), &
       &  nr1*nr2*nr3, float(neg(2)*100) / real(nr1 * nr2 * nr3)

  !
  ! energy terms, local-density contribution
  !
  vtxc = omega * vtxc / (nr1 * nr2 * nr3)
  etxc = omega * etxc / (nr1 * nr2 * nr3)
  !
  ! add gradient corrections (if any)
  !
  call gradcorr_nc (rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
       nrxx, nl, ngm, g, alat, omega, e2, etxc, vtxc, v, nspin)

#ifdef __PARA
  call reduce (1, vtxc)
  call reduce (1, etxc)
#endif
  !      call stop_clock('vxc')
  !
  return
end subroutine v_xc_nc
