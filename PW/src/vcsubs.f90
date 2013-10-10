!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!*
!*
subroutine vcinit (mxdtyp, mxdatm, ntype, natot, rat, ityp, avec, &
     vcell, force, if_pos, frr, calc, temp, vx2, vy2, vz2, rms, vmean, ekin, &
     avmod, theta, atmass, cmass, press, p, dt, aveci, avecd, avec2d, &
     avec2di, sigma, sig0, avec0, v0, rati, ratd, rat2d, rat2di, enew, &
     uta, eka, eta, ekla, utl, etl, ut, ekint, etot, iforceh)
  !
  ! rmw (18/8/99)
  ! Cesar RS Silva (04/12/2005)
  !
  ! input:
  ! mxdtyp = array dimension for type of atoms
  ! mxdatm = array dimension for atoms (irrespective of type)
  ! ntype = number of types of atoms
  ! atmass(nt) = atomic masses for atoms of type nt (in proton masses)
  ! natot = total number of atoms
  ! rat(j,na) = atomic positions in lattice coordinates
  ! ityp(na) = atomic type of na-th atom
  ! avec(3,3) = lattice vectors
  ! enew = DFT total energy
  ! calc = calculation type
  ! temp = temperature in Kelvin
  !
  ! output:
  ! rat(j,na) = atomic positions in lattice coordinates
  ! rati(j,na) = atomic positions for previous step
  ! ratd(j,na) = atomic velocities     "        "
  ! rat2d(i,na) =   "   acceleration    "        "
  ! rat2di(i,na) =   "   acceleration    "        " (previous step)
  ! avec(3,3) = lattice vectors
  ! aveci(3,3) = lattice vectors for "previous" step
  ! avecd(3,3) = 1st lattice vectors derivatives
  ! avec2d(3,3) = 2nd lattice vectors derivatives
  ! avec2di(3,3) = 2nd lattice vectors derivatives (previous step)
  ! p = internal (virial) pressure
  ! ut = new total potential energy
  ! ekin = new total kinetic energy
  ! etot = total energy
  ! we also obtain the same quantities for atomic and lattice components
  ! uta,eka,eta,utl,ekla,etl
  ! theta(3,3) = angle between lattice vectors
  ! avmod(3) = lattice vectors moduli
  !
  USE kinds
  implicit none
  !
  real(DP) :: zero, um, dois, tres, quatro, seis
  parameter (zero = 0.0d0, um = 1.0d0, dois = 2.0d0, tres = 3.0d0, &
       quatro = 4.0d0, seis = 6.0d0)
  !
  character (len=2) :: calc
  !
  integer :: mxdatm, mxdtyp
  real(DP) :: avec (3, 3), avecd (3, 3), avec2d (3, 3), avec2di (3, &
       3), aveci (3, 3), g (3, 3), gm1 (3, 3), gd (3, 3), sigma (3, 3), &
       sigav (3, 3), gmgd (3, 3), avec0 (3, 3), sig0 (3, 3), avmod (3), &
       theta (3, 3), pim (3, 3), piml (3, 3), frr (3, 3)
  !
  integer :: ityp (mxdatm), natot, if_pos(3,mxdatm), iforceh(3,3)
  real(DP) :: atmass (mxdtyp), rat (3, mxdatm), ratd (3, mxdatm), &
       rati (3, mxdatm), rat2d (3, mxdatm), rat2di (3, mxdatm)
  !
  real(DP) :: force (3, mxdatm), d2 (3, 3)
  !
  real(DP) :: vx2 (mxdtyp), vy2 (mxdtyp), vz2 (mxdtyp)
  real(DP) :: rms (mxdtyp), vmean (mxdtyp), ekin (mxdtyp)

  real(DP) :: ekint, ut, etot, tr, ekk, etl, &
       cmass, uta, enew, v0, eka, utl, ekla, eta, dt, vcell, p, press, &
       temp, ww, pv

  integer :: na, nt, i, j, l, k, m, ntype
  !
  real(DP) :: scaloff=1.0d0
  !
  IF ( COUNT( iforceh == 2 ) > 0 ) scaloff=0.5d0
  !
  ! calculate the metric for the current step
  !
  call setg (avec, g)
  !
  ! initialize cell related quantities
  !
  do j = 1, 3
     do i = 1, 3
        avecd (i, j) = zero
        avec2d (i, j) = zero
        avec2di (i, j) = zero
     enddo
  enddo
  !
  ! update metric related quantities
  !
  call updg (avec, avecd, g, gd, gm1, gmgd, sigma, vcell)
  !
  ! define reference cell
  !
  do j = 1, 3
     do i = 1, 3
        avec0 (i, j) = avec (i, j)
        sig0 (i, j) = sigma (i, j)
     enddo
  enddo
  v0 = vcell
  !
  ! establish maxwellian distribution of velocities
  !
  if (calc (2:2) .eq.'d') then
     !
     ! NB: velocities are generated in cartesian coordinates by ranv
     !     and converted to lattice coordinates immediately after.
     !     In order to avoid the use of an additional array just for
     !     this call, rat2di is used and contains therefore the velocities
     !     in cartesian coordinate. It is set to zero shortly after.
     !
     !                   I apologize, sdg. :-)
     !
     call ranv (ntype, natot, ityp, atmass, mxdtyp, mxdatm, temp, &
          ekint, rat2di, vmean, rms, vx2, vy2, vz2, ekin)
     !
     do na = 1, natot
        do l = 1, 3
           ratd(l,na) = zero
           do k = 1, 3
              IF ( if_pos(l,na) == 1 ) &
                 ratd(l,na) = rat2di(k,na) * sigma(k,l) / vcell + ratd(l,na) 
           enddo
        enddo
     enddo
  else
     do na = 1, natot
        do k = 1, 3
           ratd(k,na) = zero
        enddo
     enddo
  endif
  !
  ! define (uncorrected) accelerations and initialize rat2di
  !
  do na = 1, natot
     nt = ityp(na)
     do l = 1, 3
        rat2d (l, na) = if_pos(l,na) * force (l, na) / atmass (nt)
        rat2di(l, na) = zero
     enddo
  enddo
  !
  ! update cell related quantities
  !
  if (calc (1:1) .ne.'m') then
     !
     ! initialize piml (virial stress in lattice coordinates)
     !
     do j = 1, 3
        do i = 1, 3
           piml (i, j) = zero
        enddo
     enddo
     !
     ! correct forces on atoms
     !
     do na = 1, natot
        nt = ityp (na)
        do k = 1, 3
           do m = 1, 3
              rat2d (k, na) = rat2d (k, na) - gmgd (k, m) * ratd (m, na)
           enddo
        enddo
        !
        ! calculate virial stress in lattice coordinates
        !
        do j = 1, 3
           do i = 1, 3
              piml(i,j) = piml(i,j) + atmass(nt) * ratd(i,na) * ratd(j,na)
           enddo
        enddo
     enddo
     !
     ! calculate virial stress in cartesian coordinates
     !
     do j = 1, 3
        do i = 1, 3
           pim (i, j) = zero
           do l = 1, 3
              do m = 1, 3
                 pim(i,j) = pim(i,j) + avec(i,l) * piml(l,m) * avec(j,m)
              enddo
           enddo
        enddo
     enddo
     !
     ! add potential energy contribution to stress
     !
     do j = 1, 3
        do i = 1, 3
           pim (i, j) = (pim (i, j) + frr (i, j) ) / vcell
           avec2d (i, j) = zero
        enddo
     enddo
     !
     ! subtract external pressure from diagonal term
     !
     pim (1, 1) = pim (1, 1) - press
     pim (2, 2) = pim (2, 2) - press
     pim (3, 3) = pim (3, 3) - press
     !
     do j = 1, 3
        do i = 1, 3
           do k = 1, 3
              avec2d (i, j) = avec2d (i, j) + pim (i, k) * sigma (k, j)
           enddo
           avec2d (i, j) = avec2d (i, j) / cmass
        enddo
     enddo
     !
     ! if new cell dynamics...
     !
     if (calc (1:1) .eq.'n') then
        call sigp (avec, avecd, avec2d, sigma, vcell)
     endif
     !
     ! strain/stress symmetrization
     !
     do i = 1, 3
        do j = 1, 3
           d2 (i, j) = zero
           do k = 1, 3
              d2 (i, j) = d2 (i, j) + avec2d (i, k) * sig0 (j, k)
           enddo
           d2 (i, j) = d2 (i, j) / v0
        enddo
     enddo
     !
     d2 (1, 2) = (d2 (1, 2) + d2 (2, 1) ) / dois
     d2 (1, 3) = (d2 (1, 3) + d2 (3, 1) ) / dois
     d2 (2, 3) = (d2 (2, 3) + d2 (3, 2) ) / dois
     d2 (2, 1) = d2 (1, 2)
     d2 (3, 1) = d2 (1, 3)
     d2 (3, 2) = d2 (2, 3)
     !
     do i = 1, 3
        do j = 1, 3
           avec2d (i, j) = zero
           do k = 1, 3
              avec2d (i, j) = avec2d (i, j) + d2 (i, k) * avec0 (k, j)
           enddo
        enddo
     enddo
  else
     do i = 1, 3
        do j = 1, 3
           avec2d (i, j) = zero
        enddo
     enddo
  endif
  !
  !      WRITE( stdout,*) avec2d(2,1),avec2d(3,1), avec2d(3,2)
  !
  ! compute atomic energies
  !
  eka = zero
  do na = 1, natot
     nt = ityp (na)
     do i = 1, 3
        ekk = zero
        do j = 1, 3
           ekk = ekk + ratd (i, na) * g (i, j) * ratd (j, na)
        enddo
        eka = eka + ekk * atmass (nt) / dois
     enddo
  enddo
  uta = enew
  eta = eka + uta
  !
  !      WRITE( stdout,*) 'eka,ekint', eka, ekint
  !
  !      lattice contribution
  !
  ekla = zero
  if (calc (1:1) .ne.'m') then
     !
     ! new dynamics case
     !
     if (calc (1:1) .eq.'n') then
        do j = 1, 3
           do i = 1, 3
              sigav (i, j) = zero
              do l = 1, 3
                 sigav (i, j) = sigav (i, j) + sigma (l, i) * avecd (l, j)
              enddo
           enddo
        enddo
        do k = 1, 3
           tr = zero
           do m = 1, 3
              tr = tr + sigav (m, k) * sigav (m, k)
           enddo
           ekla = ekla + tr
        enddo
     endif
     !
     ! parrinello rahman case
     !
     if (calc (1:1) .eq.'c') then
        do k = 1, 3
           tr = zero
           do m = 1, 3
              tr = tr + avecd (m, k) * avecd (m, k)
           enddo
           ekla = ekla + tr
        enddo
     endif
  endif
  !
  ekla = ekla * cmass / dois
  utl = + press * vcell
  etl = utl + ekla
  !
  !      total energy
  !
  ekint = eka + ekla
  ut = uta + utl
  etot = ekint + ut
  !
  !      calculate "internal (virial) pressure"
  !
  ww = frr (1, 1) + frr (2, 2) + frr (3, 3)
  p = (dois * eka + ww) / tres / vcell

  pv = p * vcell
  !
  !       WRITE( stdout,1001) ekint,ut,etot
  !
  ! now make the initial move
  !
  !
  !      update atomic positions and calculate intermediate velocities
  !      and accelerations
  !
  do na = 1, natot
     do k = 1, 3
        rati (k, na) = rat (k, na)
        rat (k, na) = rat (k, na) + dt * ratd (k, na) + dt * dt * (quatro &
             * rat2d (k, na) - rat2di (k, na) ) / seis
        rat2di (k, na) = rat2d (k, na)
     enddo
  enddo
  !
  ! update lattice vectors if cell dynamics
  !
  if (calc (1:1) .ne.'m') then
     do j = 1, 3
        do i = 1, 3
           aveci (i, j) = avec (i, j)
           avec (i, j) = avec (i, j) + dt * avecd (i, j) + (dt * dt * &
                (quatro * avec2d (i, j) - avec2di (i, j) ) / seis)  * dble(iforceh(i,j))*scaloff
           avec2di (i, j) = avec2d (i, j)
        enddo
     enddo
     !
     ! update cell quantities just in case forclj need them
     !
     call updg (avec, avecd, g, gd, gm1, gmgd, sigma, vcell)

  endif

  return

end subroutine vcinit
!*
!*
subroutine vcmove (mxdtyp, mxdatm, ntype, ityp, rat, avec, vcell, &
     force, if_pos, frr, calc, avmod, theta, atmass, cmass, press, p, dt, &
     avecd, avec2d, aveci, avec2di, sigma, sig0, avec0, v0, ratd, &
     rat2d, rati, rat2di, enew, uta, eka, eta, ekla, utl, etl, ut, &
     ekint, etot, temp, tolp, ntcheck, ntimes, nst, tnew, nzero, natot, &
     acu, ack, acp, acpv, avu, avk, avp, avpv, iforceh)
  !
  !      rmw (18/8/99)
  !
  !     input:
  !     mxdtyp = array dimension for type of atoms
  !     mxdatm = array dimension for atoms (irrespective of type)
  !     ntype = number of types of atoms
  !     atmass(nt) = atomic masses for atoms of type nt (in proton masses)
  !     ityp(na) = atomic type of na-th atom
  !     rat(j,na) = atomic positions in lattice coordinates
  !     rati(j,na) = atomic positions in lattice coordinates (previous ste
  !     ratd(j,na) = atomic velocities     "        "
  !     rat2di(i,na) =   "   acceleration    "        " (previous step)
  !     avec(3,3) = lattice vectors
  !     aveci(3,3) = lattice vectors (previous step)
  !     avecd(3,3) = 1st lattice vectors derivatives
  !     avec2d(3,3) = 2nd lattice vectors derivatives
  !     avec2di(3,3) = 2nd lattice vectors derivatives (previous step)
  !     avec0(3,3) = initial lattice vectors
  !     sig0(3,3) = initial reciprocal lattice vectors *  vcell / 2 pi
  !     v0 = initial volume
  !     enew = DFT total energy
  !
  !     output:
  !     rat(j,na) = atomic positions in lattice coordinates (updated)
  !     ratd(j,na) = atomic velocities     "        " (updated)
  !     rat2d(i,na) =   "   acceleration    "        " (updated)
  !     rati(j,na) and rat2di(i,na) (updated)
  !     avec(3,3) = lattice vectors
  !     avecd(3,3) = 1st lattice vectors derivatives
  !     avec2d(3,3) = 2nd lattice vectors derivatives
  !     aveci(3,3) and avec2di(3,3) (updated)
  !     p = internal (virial) pressure
  !     ut = new total potential energy
  !     ekin = new total kinetic energy
  !     etot = total energy
  !     we also obtain the same quantities for atomic and lattice componen
  !     uta,eka,eta,utl,ekl,etl
  !     theta(3,3) = angle between lattice vectors
  !     avmod(3) = lattice vectors moduli
  !
  !
  USE kinds,         only : DP
  USE constants,     ONLY : pi, eps16, k_boltzmann_ry
  USE io_global,     ONLY : stdout
  
  implicit none
  !
  real(DP) :: zero, um, dois, tres, quatro, seis
  parameter (zero = 0.0d0, um = 1.0d0, dois = 2.0d0, tres = 3.0d0, &
       quatro = 4.0d0, seis = 6.0d0)
  !
  character (len=2) :: calc
  !
  integer :: mxdatm, mxdtyp

  integer :: ityp (mxdatm), if_pos(3,mxdatm), iforceh(3,3)
  real(DP) :: avec (3, 3), rat (3, mxdatm)
  !
  real(DP) :: atmass (mxdtyp), ratd (3, mxdatm), rat2d (3, mxdatm), &
       avecd (3, 3), avec2d (3, 3), g (3, 3), gm1 (3, 3), gd (3, 3), &
       sigma (3, 3), avec0 (3, 3), sig0 (3, 3), avmod (3), theta (3, 3), &
       pim (3, 3), piml (3, 3), frr (3, 3), rati (3, mxdatm), rat2di (3, &
       mxdatm), sigav (3, 3), gmgd (3, 3), aveci (3, 3), avec2di (3, 3)

  integer :: i, j, k, l, m, na, nt, nst, natot, nzero, ntimes, &
       ntcheck, ntype, i_update, n_update

  real(DP) :: avpv,  pv, ww, ts, xx, alpha, x, &
       tr, tnew, tolp, temp, avk, avu, ekk, avp, ack, acu, acpv, acp, dt, &
       p, enew, v0, vcell, press, ut, etl, etot, ekint, utl, uta, cmass, &
       eka, ekla, eta
  logical :: symmetrize_stress
  !
  real(DP) :: force (3, mxdatm), d2 (3, 3)
  !
  real(DP) :: scaloff=1.0d0
  !
  IF ( COUNT( iforceh == 2 ) > 0 ) scaloff=0.5d0
  !
  !
  !      zero energy components
  !
  ut = zero
  ekint = zero
  etot = zero
  uta = zero
  eka = zero
  eta = zero
  utl = zero
  ekla = zero
  etl = zero
  p = zero
  !
  ! set the metric for the current step
  !
  
  call setg (avec, g)
  !
  ! calculate (uncorrected) rat2d
  !
  
  do na = 1, natot
     nt = ityp (na)
     do i = 1, 3
        rat2d (i, na) = if_pos(i,na) * force (i, na) / atmass (nt)
     enddo
  enddo
  !
  ! if variable cell, estimate velocities and set the number of update to
  ! be performed in order to have them accurate. This is needed only for
  ! variable cell shape dynamics (where accelerations depends on velocities)
  ! and a few, even just one, iteration is usually enough
  !
  if (calc (1:1) .ne.'m') then
     do na = 1, natot
        do k = 1, 3
           ratd (k, na) = ratd (k, na) + dt * rat2di (k, na)
        enddo
     enddo
     do j = 1, 3
        do i = 1, 3
           avecd (i, j) = avecd (i, j) + dt * avec2di (i, j) !* dble(iforceh(i,j))
        enddo
     enddo
     n_update = 19
  else
     n_update = 1
  endif

  do i_update = 1, n_update
     if (calc (1:1) .ne.'m') then
        !
        !      update metric related quantities
        !
        call updg (avec, avecd, g, gd, gm1, gmgd, sigma, vcell)
        !
        !  zero piml (virial stress in lattice coordinates)
        !
        do j = 1, 3
           do i = 1, 3
              piml (i, j) = zero
           enddo
        enddo
        !
        !  correct forces on atoms and set cell forces
        !
        do na = 1, natot
           nt = ityp (na)
           do k = 1, 3
              rat2d (k, na) = if_pos(k,na) * force (k, na) / atmass (nt)
              do m = 1, 3
                 rat2d (k, na) = rat2d (k, na) - gmgd (k, m) * ratd (m, na)
              enddo
           enddo
           !
           !  calculate virial stress in lattice coordinates
           !
           do j = 1, 3
              do i = 1, 3
                 piml(i,j) = piml(i,j) + atmass(nt) * ratd(i,na) * ratd(j,na)
              enddo
           enddo
        enddo
        !
        !  calculate virial stress in cartesian coordinates
        !
        do i = 1, 3
           do j = 1, 3
              pim (i, j) = zero
              do l = 1, 3
                 do m = 1, 3
                    pim(i,j) = pim(i,j) + avec(i,l) * piml(l,m) * avec(j,m)
                 enddo
              enddo
           enddo
        enddo
        !
        !  add potential energy contribution to stress
        !
        do j = 1, 3
           do i = 1, 3
              pim (i, j) = (pim (i, j) + frr (i, j) ) / vcell
           enddo
        enddo
        !
        ! subtract external pressure from diagonal term
        !
        pim (1, 1) = pim (1, 1) - press
        pim (2, 2) = pim (2, 2) - press
        pim (3, 3) = pim (3, 3) - press
        !
        do j = 1, 3
           do i = 1, 3
              avec2d (i, j) = zero
              do k = 1, 3
                 avec2d (i, j) = avec2d (i, j) + pim (i, k) * sigma (k, j)
              enddo
              avec2d (i, j) = avec2d (i, j) / cmass
           enddo
        enddo
        !
        !      if new cell dynamics...
        !
        if (calc (1:1) .eq.'n') call sigp (avec, avecd, avec2d, sigma, vcell)
        !
        !      strain/stress symmetrization
        !
        symmetrize_stress = .true.

        if (.not.symmetrize_stress) goto 666
        do i = 1, 3
           do j = 1, 3
              d2 (i, j) = zero
              do k = 1, 3
                 d2 (i, j) = d2 (i, j) + avec2d (i, k) * sig0 (j, k)
              enddo
              d2 (i, j) = d2 (i, j) / v0
           enddo
        enddo
        !
        d2 (1, 2) = (d2 (1, 2) + d2 (2, 1) ) / dois
        d2 (1, 3) = (d2 (1, 3) + d2 (3, 1) ) / dois
        d2 (2, 3) = (d2 (2, 3) + d2 (3, 2) ) / dois
        d2 (2, 1) = d2 (1, 2)
        d2 (3, 1) = d2 (1, 3)
        d2 (3, 2) = d2 (2, 3)
        !
        do i = 1, 3
           do j = 1, 3
              avec2d (i, j) = zero
              do k = 1, 3
                 avec2d (i, j) = avec2d (i, j) + d2 (i, k) * avec0 (k, j)
              enddo
           enddo

        enddo

666     continue
        !
        ! calculate correct lattice velocities and ...
        !
        do j = 1, 3
           do i = 1, 3
              avecd (i, j) = (avec (i, j) - aveci (i, j) ) / dt + (dt * &
                   (dois * avec2d (i, j) + avec2di (i, j) ) / seis)  * dble(iforceh(i,j))*scaloff
           enddo

        enddo

     endif
     !
     ! calculate correct atomic velocities
     !
     do na = 1, natot
        do k = 1, 3
           ratd (k, na) = (rat (k, na) - rati (k, na) ) / dt + dt * (dois * &
                rat2d (k, na) + rat2di (k, na) ) / seis
        enddo
     enddo
     ! and do-loop on n_update
  enddo
  !
  !      calculate basis vectors' moduli and angles
  !
  if (calc (1:1) .ne.'m') then
     do k = 1, 3
        avmod (k) = zero
        do l = 1, 3
           theta (l, k) = zero
           avmod (k) = avmod (k) + avec (l, k) * avec (l, k)
           do m = 1, 3
              theta (l, k) = theta (l, k) + avec (m, l) * avec (m, k)
           enddo
        enddo
        avmod (k) = dsqrt (avmod (k) )
     enddo
     do k = 1, 3
        do l = 1, 3
           x = theta (l, k) / avmod (l) / avmod (k)
           if (x.ge.0.d0) then
              x = dmin1 (1.d0, x)
           else
              x = dmax1 ( - 1.d0, x)
           endif
           theta (l, k) = dacos (x) * 180.d0 / pi
        enddo
     enddo
  endif
  !
  !     compute atomic energies
  !
  do na = 1, natot
     nt = ityp (na)
     do i = 1, 3
        ekk = zero
        do j = 1, 3
           ekk = ekk + ratd (i, na) * g (i, j) * ratd (j, na)
        enddo
        eka = eka + ekk * atmass (nt) / dois
     enddo
  enddo
  !
  uta = enew
  eta = eka + uta
  !
  !      lattice contribution
  !
  ekla = zero

  if (calc (1:1) .ne.'m') then
     if (calc (1:1) .eq.'n') then
        !
        ! new dynamics or new minimization cases
        !
        do j = 1, 3
           do i = 1, 3
              sigav (i, j) = zero
              do l = 1, 3
                 sigav (i, j) = sigav (i, j) + sigma (l, i) * avecd (l, j)
              enddo
           enddo
        enddo
        do k = 1, 3
           tr = zero
           do m = 1, 3
              tr = tr + sigav (m, k) * sigav (m, k)
           enddo
           ekla = ekla + tr
        enddo
     endif
     !
     if (calc (1:1) .eq.'c') then
        !
        ! cell dynamics or cell minimization cases
        !
        do k = 1, 3
           tr = zero
           do m = 1, 3
              tr = tr + avecd (m, k) * avecd (m, k)
           enddo
           ekla = ekla + tr
        enddo
     endif
  endif
  !
  ekla = ekla * cmass / dois
  utl = + press * vcell
  etl = utl + ekla
  !
  !      total energy
  !
  ekint = eka + ekla
  ut = uta + utl
  etot = ekint + ut
  !
  !      calculate "internal (virial) pressure"
  !
  ww = frr (1, 1) + frr (2, 2) + frr (3, 3)
  p = (dois * eka + ww) / tres / vcell
  pv = p * vcell
  !
  !       update accumulators and set averages
  !
  nzero = nzero + 1
  acu = acu + ut
  ack = ack + ekint
  acp = acp + p
  acpv = acpv + pv
  avu = acu / DBLE (nzero)
  avk = ack / DBLE (nzero)
  avp = acp / DBLE (nzero)
  avpv = acpv / DBLE (nzero)
  !
  !       choose # of degrees of freedom and calculate tnew
  !
  if (calc (1:1) .ne.'m') then
     tnew = dois / tres / DBLE (natot + 1) * avk / k_boltzmann_ry
  else
     tnew = dois / tres / DBLE (natot - 1) * avk / k_boltzmann_ry
  endif
  !
  !       rescale velocities
  !
  if ( mod (nst, ntcheck) == 0 ) then
     !
     !       with the new definition of tolp, this is the test to perform
     !
     if ( ( ABS (tnew - temp ) > tolp) .and. ( abs(ntimes) > 0) )  then
        !
        if ( tnew < 1.0d-12) then
           alpha = 1.0_dp
        else
           alpha = sqrt (temp / tnew)
        endif
        do na = 1, natot
           do k = 1, 3
              ratd (k, na) = alpha * ratd (k, na)
           enddo
        enddo
        if (calc (2:2) .eq.'d') then
           do k = 1, 3
              do l = 1, 3
                 avecd (l, k) = alpha * avecd (l, k)  !* dble(iforceh(i,j))
              enddo
           enddo
        endif
        !
        ! update ntimes and nzero and reset accumulators
        !
        acu = zero
        ack = zero
        acp = zero
        acpv = zero
        if ( ntimes > 0 ) ntimes = ntimes - 1
        nzero = 0
     endif

  endif
  if (calc (2:2) .eq.'m') then
     !         WRITE( stdout,109) alpha,nst
     ! if(.true. ) = original version modified by Cesar Da Silva
     ! if(.false.) = modified algorithm by SdG
     if (.false.) then
        do na = 1, natot
           do k = 1, 3
              xx = rat2di (k, na) * rat2d (k, na)
              if (xx.lt.zero) then
                 ratd (k, na) = zero
                 rat(k,na)=rat2d(k,na)*rati(k,na)-rat2di(k,na)*rat(k,na)
                 rat(k,na)=rat(k,na)/(rat2d(k,na)-rat2di(k,na))
                 rat2d(k,na)=zero
                 rat2di(k,na)=zero
              endif
           enddo
        enddo
     else
        do na = 1, natot
           xx = 0.d0
           do k=1,3
              xx = rat2d(1,na) * g(1,k) * ratd(k,na) + &
                   rat2d(2,na) * g(2,k) * ratd(k,na) + &
                   rat2d(3,na) * g(3,k) * ratd(k,na) + xx
           end do
           if (xx.gt.eps16) then
              ratd (:,na) =  rat2d (:,na) * xx
              xx = 0.d0
              do k=1,3
                 xx = rat2d(1,na) * g(1,k) * rat2d(k,na) + &
                      rat2d(2,na) * g(2,k) * rat2d(k,na) + &
                      rat2d(3,na) * g(3,k) * rat2d(k,na) + xx
              end do
              ratd(:,na) = ratd(:,na) / xx
           else
              ratd(:, na) = zero
           endif
        enddo
     endif

     if (calc (1:1) .ne.'m') then
        do k = 1, 3
           do l = 1, 3
              xx = avec2d (l, k) * avec2di (l, k)
              if (xx.lt.zero) then
                 avecd (l, k) = zero
                 avec(l, k)=avec2d(l,k)*aveci(l,k)-avec2di(l,k)*avec(l,k)
                 avec(l, k)=avec(l,k)/(avec2d(l,k)-avec2di(l,k))
                 avec2d(l,k)=zero
                 avec2di(l,k)=zero
              endif
           enddo
        enddo
     endif
  endif
  !
  !      update atomic positions and calculate intermediate velocities
  !      and accelerations
  !
  do na = 1, natot
     do k = 1, 3
        rati (k, na) = rat (k, na)
        rat (k, na) = rat (k, na) + dt * ratd (k, na) + dt * dt * (quatro &
             * rat2d (k, na) - rat2di (k, na) ) / seis
        rat2di (k, na) = rat2d (k, na)
     enddo
  enddo
  !
  !      update lattice vectors if cell dynamics
  !
  if (calc (1:1) .ne.'m') then
     do j = 1, 3
        do i = 1, 3
           aveci (i, j) = avec (i, j)
           avec (i, j) = avec (i, j) + dt * avecd (i, j) + (dt * dt * &
                (quatro * avec2d (i, j) - avec2di (i, j) ) / seis)  * dble(iforceh(i,j))*scaloff
           avec2di (i, j) = avec2d (i, j)
        enddo
     enddo
  endif
  !
  ! update metric related quantities just in case are needed by forclj
  !

  call updg (avec, avecd, g, gd, gm1, gmgd, sigma, vcell)

  return
end subroutine vcmove
!*
!*
subroutine ranv (ntype, natot, ityp, atmass, mxdtyp, mxdatm, temp, &
     ekint, v, vmean, rms, vx2, vy2, vz2, ekin)
  !
  !     sets up random velocities with maxwellian distribution
  !     at temperature t. total linear momentum components are zero
  !     rewritten on 1/31/90 by rmw
  !     extracted from car & parrinello 's program
  !
  !     input:
  !     mxdtyp = array dimension for type of atoms
  !     mxdatm = array dimension for atoms (irrespective of type)
  !     ntype = number of types of atoms
  !     natot = total number of atoms
  !     ityp(na) = atomic type of na-th atom
  !     atmass(i) = atomic masses for atoms of type i (in proton masses)
  !     temp = temperature in k
  !
  !     output:
  !     v(i,na) = initial velocity of atom na of type nt
  !     vmean(nt), rms(nt),vx2(nt),vy2(nt),vz2(nt)
  !
  USE io_global,   ONLY : stdout
  USE constants,   ONLY : k_boltzmann_ry
  USE kinds , only : DP
  implicit none
  !

  integer :: mxdtyp, mxdatm
  real(DP) :: atmass (mxdtyp)
  real(DP) :: v (3, mxdatm), p (3)
  real(DP) :: vx2 (mxdtyp), vy2 (mxdtyp), vz2 (mxdtyp)
  real(DP) :: rms (mxdtyp), vmean (mxdtyp), ekin (mxdtyp)
  !

  integer :: ityp (mxdatm), natot

  integer :: na, nt, j, k, ntype, iseed, natom


  real(DP) :: ran3, vfac, sig, tfac, vr, atemp, eps, temp, ekint, t
  real(DP) :: b0, b1, c0, c1
  real(DP) :: zero, um, dois, tres
  data b0, b1, c0, c1 / 2.30753d0, 0.27061d0, 0.99229d0, 0.04481d0 /
  data zero, um, dois, tres / 0.d0, 1.d0, 2.d0, 3.0d0 /
  !
  !     example run
  !
  do nt = 1, ntype
     ekin (nt) = zero
  enddo
  ekint = zero
  !
  !
  if (natot.ne.1) then
     !
     !     assign random velocities
     !
     t = temp
     if (temp.lt.1.d-14) t = 1.d-14
     iseed = - 119
     eps = ran3 (iseed)
     !
     !     establish gaussian distribution for each atom kind
     !
     ! natom (the number of atoms of a given type) is calculated when needed
     !
     do nt = 1, ntype
        natom = 0
        vfac = dsqrt (k_boltzmann_ry * t / atmass (nt) )
        !            WRITE( stdout,901)
        !            WRITE( stdout,*) 'vfac = ',vfac
        iseed = iseed+382
        do na = 1, natot
           if (ityp (na) .eq.nt) then
              natom = natom + 1
              do j = 1, 3
                 eps = ran3 (iseed)
                 if (eps.lt.1.d-10) eps = 1.d-10
                 if (eps.le.0.5d0) goto 100
                 eps = eps - um
                 if (eps.gt. - 1.d-10) eps = - 1.d-10
100              sig = dsqrt (log (um / (eps * eps) ) )
                 vr = sig - (b0 + b1 * sig) / (um + c0 * sig + c1 * sig * &
                      sig)
                 vr = vr * vfac
                 if (eps.lt.zero) vr = - vr
                 v (j, na) = vr
              enddo
           endif
        enddo
        !
        p (1) = zero
        p (2) = zero
        p (3) = zero
        ekin (nt) = zero
        if (natom.eq.0) then
           WRITE( stdout,*) 'natom=0 for type',nt,'in sub ranv (1) !!!! '
           go to 111
        end if
        !
        !       calculate linear-momentum.
        !
        do na = 1, natot
           if (ityp (na) .eq.nt) then
              p (1) = p (1) + v (1, na)
              p (2) = p (2) + v (2, na)
              p (3) = p (3) + v (3, na)
           endif
        enddo
        p (1) = p (1) / DBLE (natom)
        p (2) = p (2) / DBLE (natom)
        p (3) = p (3) / DBLE (natom)
        !
        !       zero linear momentum for atom type nt
        !
        do na = 1, natot
           if (ityp (na) .eq.nt) then
              v (1, na) = v (1, na) - p (1)
              v (2, na) = v (2, na) - p (2)
              v (3, na) = v (3, na) - p (3)
           endif
        enddo
        do na = 1, natot
           if (ityp (na) .eq.nt) then
              ekin(nt) = ekin(nt) + ( v(1,na)*v(1,na) + v(2,na)*v(2,na) + &
                                      v(3,na)*v(3,na) ) / dois
           endif
        enddo
        !            WRITE( stdout,*) 'ekin(nt)',ekin(nt)
        ekin (nt) = atmass (nt) * ekin (nt)
        ekint = ekint + ekin (nt)
 111    continue
     enddo
     !
     !     rescale velocities to give correct temperature
     !
     atemp = dois * ekint / tres / DBLE (natot - 1) / k_boltzmann_ry
     tfac = dsqrt (t / atemp)
     if (temp.lt.1d-14) tfac = zero
     !         WRITE( stdout,*) 'atemp = ',atemp,' k'
     !         WRITE( stdout,*) 'tfac = ',tfac
     do nt = 1, ntype
        vmean (nt) = zero
        rms (nt) = zero
        vx2 (nt) = zero
        vy2 (nt) = zero
        vz2 (nt) = zero
     enddo
     do na = 1, natot
        nt = ityp (na)
        v (1, na) = v (1, na) * tfac
        v (2, na) = v (2, na) * tfac
        v (3, na) = v (3, na) * tfac
        vmean(nt) = vmean(nt) + dsqrt (v(1,na)**2 + v(2,na)**2 + v(3,na)**2)
        vx2 (nt) = vx2 (nt) + v (1, na) **2
        vy2 (nt) = vy2 (nt) + v (2, na) **2
        vz2 (nt) = vz2 (nt) + v (3, na) **2
     enddo
     do nt = 1, ntype
        natom = 0
        do na = 1, natot
           if (ityp (na) .eq.nt) natom = natom + 1
        enddo

        if (natom.gt.0) then
           vmean (nt) = vmean (nt) / DBLE (natom)
           rms (nt) = dsqrt ( (vx2 (nt) + vy2 (nt) + vz2 (nt) ) /  &
                      DBLE ( natom) )
           vx2 (nt) = dsqrt (vx2 (nt) / DBLE (natom) )
           vy2 (nt) = dsqrt (vy2 (nt) / DBLE (natom) )
           vz2 (nt) = dsqrt (vz2 (nt) / DBLE (natom) )
        else
           vmean (nt) = zero
           rms (nt) = zero
           vx2 (nt) = zero
           vy2 (nt) = zero
           vz2 (nt) = zero
        end if
     enddo
     ekint = ekint * tfac * tfac
  else
     ekint = zero
     do k = 1, 3
        v (k, 1) = zero
     enddo
     vmean (1) = zero
     rms (1) = zero
     vx2 (1) = zero
     vy2 (1) = zero
     ekin (1) = zero

  endif

  return

end subroutine ranv
!*
!*
subroutine sigp (avec, avecd, avec2d, sigma, vcell)
  !
  !     calculates sigmap matrices and avec2d for
  !     new dynamics(rmw 5/30/90)
  !
  !      input:
  !            avec = lattice vectors
  !            avecd = time derivative of lattice vectors
  !            avec2d = 2nd time derivative of lattice vectors
  !            sigma = volume * rec. latt. vectors / 2 pi
  !            vcell = cell volume
  !
  !      output:
  !            avec2d = new 2nd time derivative of lattice vectors
  !
  USE kinds, only : DP
  implicit none
  !
  real(DP) :: avec (3, 3), avecd (3, 3), avec2d (3, 3), sigmap (3, &
       3, 3, 3), sigmad (3, 3), sigma (3, 3), e (3, 3), fp (3, 3, 3, 3), &
       fd (3, 3), fm1 (3, 3), fm (3, 3), sm (3, 3), avint (3, 3), &
       vcell
  !

  integer :: i, j, k, l, m, n
  real(DP) :: zero, dois
  parameter (zero = 0.d0, dois = 2.d0)
  !
  ! sigmap_ijkl = d sigma_ij / d h_kl
  !             =( sigma_ij * sigma_kl - sigma_kj * sigma_il ) / vcell
  !
  do i = 1, 3
     do j = 1, 3
        do k = 1, 3
           do l = 1, 3
              sigmap(i,j,k,l) = ( sigma(i,j)*sigma(k,l) - &
                                  sigma(k,j)*sigma(i,l) ) / vcell
           enddo
        enddo
     enddo
  enddo
  !                _1  t           2
  !     calculate f = h * h / vcell
  !
  do j = 1, 3
     do i = 1, 3
        fm1 (i, j) = zero
        do l = 1, 3
           fm1 (i, j) = fm1 (i, j) + avec (l, i) * avec (l, j)
        enddo
        fm1 (i, j) = fm1 (i, j) / vcell / vcell
     enddo
  enddo
  !                   .t  .
  !     calculate e = h * h
  !
  do j = 1, 3
     do i = 1, 3
        e (i, j) = zero
        do m = 1, 3
           e (i, j) = e (i, j) + avecd (m, i) * avecd (m, j)
        enddo
     enddo
  enddo
  !                          ij t             t       ij
  !     calculate f' = sigma'  * sigma + sigma * sigma'
  !
  do n = 1, 3
     do m = 1, 3
        do j = 1, 3
           do i = 1, 3
              fp(i,j,m,n) = zero
              do l = 1, 3
                 fp(i,j,m,n) = fp(i,j,m,n) + sigmap(i,j,l,m) * sigma(l,n) + &
                                             sigma(l,m) * sigmap(i,j,l,n)
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  !     calculate sigmad
  !
  do n = 1, 3
     do m = 1, 3
        sigmad(m,n) = zero
        do j = 1, 3
           do i = 1, 3
              sigmad(m,n) = sigmad(m,n) + sigmap(i,j,m,n)*avecd(i,j)
           enddo
        enddo
     enddo
  enddo
  !               .
  !     calculate f
  !
  do j = 1, 3
     do i = 1, 3
        fd(i,j) = zero
        do l = 1, 3
           fd(i,j) = fd(i,j) + sigmad(l,i)*sigma(l,j) + sigma(l,i)*sigmad(l,j)
        enddo
     enddo
  enddo
  !
  !     calculate fm
  !
  do j = 1, 3
     do i = 1, 3
        fm (i, j) = zero
        do l = 1, 3
           do k = 1, 3
              fm (i, j) = fm (i, j) + e (l, k) * fp (i, j, k, l)
           enddo
        enddo
        fm (i, j) = fm (i, j) / dois
     enddo
  enddo
  !
  !     calculate sm
  !
  do j = 1, 3
     do i = 1, 3
        sm (i, j) = zero
        do l = 1, 3
           sm (i, j) = sm (i, j) + avecd (i, l) * fd (l, j)
        enddo
     enddo
  enddo
  !
  !     calculate new avec2d
  !
  do j = 1, 3
     do i = 1, 3
        avint (i, j) = avec2d (i, j) + fm (i, j) - sm (i, j)
     enddo
  enddo
  !
  !
  do j = 1, 3
     do i = 1, 3
        avec2d (i, j) = zero
        do m = 1, 3
           avec2d (i, j) = avec2d (i, j) + avint (i, m) * fm1 (m, j)
        enddo
     enddo
  enddo
  !
  return
end subroutine sigp
!*
!*
subroutine updg (avec, avecd, g, gd, gm1, gmgd, sigma, vcell)
  !
  !
  !     update metric related quantities
  !     (rmw 18/8/99)
  !
  !      input:
  !      avec(3,3) = lattice vectors
  !      avecd(3,3) = derivative of lattice vectors
  !
  !      output:      t
  !      g(3,3) = avec * avec
  !                     t                     t
  !      gd(3,3) = avecd * avec + avecd * avec
  !                  _1
  !      gm1(3,3) = g
  !                   _1
  !      gmgd(3,3) = g * gd
  !      sigma(3,3) = reciprocal lattice vectors / twopi
  !      vcell = cell volume
  !
  USE kinds, only : DP
  implicit none
  !
  real(DP) :: zero, um, dois, tres
  parameter (zero = 0.0d0, um = 1.0d0, dois = 2.0d0, tres = 3.0d0)
  !
  real(DP) :: avec (3, 3), avecd (3, 3), sigma (3, 3)
  real(DP) :: g (3, 3), gd (3, 3), gmgd (3, 3), gm1 (3, 3)

  real(DP) :: vcell
  integer :: i, j, m
  !
  !     compute the lattice wave-vectors/twopi and the cell volume
  !
  !     vcell = abs (det (h_ij))       ! NOTE the abs value !
  !
  !     sigma_ij = d vcell / d h_ij
  !
  sigma (1, 1) = avec (2, 2) * avec (3, 3) - avec (3, 2) * avec (2, 3)
  sigma (2, 1) = avec (3, 2) * avec (1, 3) - avec (1, 2) * avec (3, 3)
  sigma (3, 1) = avec (1, 2) * avec (2, 3) - avec (2, 2) * avec (1, 3)
  sigma (1, 2) = avec (2, 3) * avec (3, 1) - avec (3, 3) * avec (2, 1)
  sigma (2, 2) = avec (3, 3) * avec (1, 1) - avec (1, 3) * avec (3, 1)
  sigma (3, 2) = avec (1, 3) * avec (2, 1) - avec (2, 3) * avec (1, 1)
  sigma (1, 3) = avec (2, 1) * avec (3, 2) - avec (3, 1) * avec (2, 2)
  sigma (2, 3) = avec (3, 1) * avec (1, 2) - avec (1, 1) * avec (3, 2)
  sigma (3, 3) = avec (1, 1) * avec (2, 2) - avec (2, 1) * avec (1, 2)
  !
  !     compute cell volume and modify sigma if needed
  !

  vcell = sigma (1, 1) * avec (1, 1) + sigma (2, 1) * avec (2, 1) &
        + sigma (3, 1) * avec (3, 1)
  if (vcell.lt.0.d0) then
     vcell = - vcell
     do i = 1, 3
        do j = 1, 3
           sigma (i, j) = - sigma (i, j)
        enddo
     enddo
  endif
  !
  !      calculate g, gd, and gm1 matrices
  !
  do j = 1, 3
     do i = 1, 3
        g (i, j) = zero
        gm1 (i, j) = zero
        gd (i, j) = zero
     enddo
  enddo
  do j = 1, 3
     do i = 1, 3
        do m = 1, 3
           g(i, j)   = g(i, j)   + avec(m,i)*avec(m,j)
           gm1(i, j) = gm1(i, j) + sigma(m,i)*sigma(m,j)
           gd(i, j)  = gd(i, j)  + avec(m,i)*avecd(m,j) + avecd(m,i)*avec(m,j)
        enddo
        gm1(i,j) = gm1(i,j) / vcell / vcell
     enddo
  enddo
  !                _1  .
  !     calculate g * g    ( = gmgd)
  !
  do j = 1, 3
     do i = 1, 3
        gmgd (i, j) = zero
        do m = 1, 3
           gmgd (i, j) = gmgd (i, j) + gm1 (i, m) * gd (m, j)
        enddo
     enddo
  enddo

  return
end subroutine updg
!*
!*
subroutine setg (avec, g)
  !
  !
  !     update metric related quantities
  !     (rmw 18/8/99)
  !
  !      input:
  !      avec(3,3) = lattice vectors
  !
  !      output:      t
  !      g(3,3) = avec * avec
  !
  USE kinds, only : DP
  implicit none
  !
  real(DP) :: zero
  parameter (zero = 0.0d0)
  !
  real(DP) :: avec (3, 3), g (3, 3)
  integer :: i, j, m
  !
  !     calculate g
  !
  do j = 1, 3
     do i = 1, 3
        g (i, j) = zero
     enddo
  enddo
  do j = 1, 3
     do i = 1, 3
        do m = 1, 3
           g (i, j) = g (i, j) + avec (m, i) * avec (m, j)
        enddo
     enddo

  enddo
  return
end subroutine setg
!*
!*
real(8) function ran3 (idum)
  USE kinds, only : DP
  implicit none

  save
  !         implicit real*4(m)
  !         parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
  integer :: mbig, mseed, mz
  real(DP) :: fac
  parameter (mbig = 1000000000, mseed = 161803398, mz = 0, fac = 1.d-9)

  integer :: ma (55), iff, k, inext, inextp, ii, mj, idum, i, mk
  !     common /ranz/ ma,inext,inextp
  data iff / 0 /
  if (idum.lt.0.or.iff.eq.0) then
     iff = 1
     mj = mseed-iabs (idum)
     mj = mod (mj, mbig)
     ma (55) = mj
     mk = 1
     do i = 1, 54
        ii = mod (21 * i, 55)
        ma (ii) = mk
        mk = mj - mk
        if (mk.lt.mz) mk = mk + mbig
        mj = ma (ii)
     enddo
     do k = 1, 4
        do i = 1, 55
           ma (i) = ma (i) - ma (1 + mod (i + 30, 55) )
           if (ma (i) .lt.mz) ma (i) = ma (i) + mbig
        enddo
     enddo
     inext = 0
     inextp = 31
     idum = 1
  endif
  inext = inext + 1
  if (inext.eq.56) inext = 1
  inextp = inextp + 1
  if (inextp.eq.56) inextp = 1
  mj = ma (inext) - ma (inextp)
  if (mj.lt.mz) mj = mj + mbig
  ma (inext) = mj
  ran3 = mj * fac
  return
end function ran3

