!-----------------------------------------------------------------------
!
!      EPA (electron-phonon-averaged) approximation
!      Averages electron-phonon matrix elements over k and q wavevectors
!          in the bins of incident and scattered electron energy grids
!      Written by Georgy Samsonidze from 2015-01-28 to 2015-05-11
!
!      Adv. Energy Mater. 2018, 1800246
!      doi:10.1002/aenm.201800246
!      https://doi.org/10.1002/aenm.201800246
!
!-----------------------------------------------------------------------

program epa
  use kinds, only : dp
  use constants, only : ry_to_cmm1, rytoev, degspin, eps12
  implicit none

  integer, parameter :: nwin = 2
  integer, parameter :: uni = 101
  integer, parameter :: uno = 102
  real(dp), parameter :: epsw2 = (20_dp / ry_to_cmm1)**2

  logical :: minus_q
  integer :: nmodes, nqs, nspin, nbnd, nkstot, nksqtot, &
      nat, nsymq, irotmq, iq, ik, ikk, ikq, ibnd, jbnd, &
      nu, mu, vu, ipert, jpert, ii, jj, kk, ll, ijob, &
      nev, bvmin, bvmax, nec, bcmin, bcmax, ndist, nepmax, &
      nemin, nemax, ne(nwin), nepair(nwin), bpair(2, nwin), &
      s(3, 3, 48), invs(48)
  real(dp) :: wtot, weight, factor, gbuf, gsum, &
      gmax, gstep, ev, dev, ec, dec, wspin, &
      ee(nwin), de(nwin), xq(3), at(3, 3), bg(3, 3)
  character(len=256) :: fni, fno, job, fmt
  character(len=32) :: s1, s2, s3, s4, s5

  real(dp), allocatable :: gdist(:,:,:)
  real(dp), allocatable :: gepair(:,:,:)
  real(dp), allocatable :: gavg(:,:,:,:)
  real(dp), allocatable :: wavg(:)
  real(dp), allocatable :: gtot(:,:,:)
  integer, allocatable :: gnum(:,:,:)
  real(dp), allocatable :: egrid(:,:,:)
  real(dp), allocatable :: x_q(:,:)
  real(dp), allocatable :: wq(:)
  real(dp), allocatable :: w2(:)
  complex(dp), allocatable :: dyn(:,:)
  real(dp), allocatable :: wk(:)
  real(dp), allocatable :: et(:,:)
  integer, allocatable :: ikks(:)
  integer, allocatable :: ikqs(:)
  integer, allocatable :: irt(:,:)
  real(dp), allocatable :: rtau(:,:,:)
  complex(dp), pointer :: u(:,:)
  complex(dp), allocatable :: el_ph_mat(:,:,:,:)
  complex(dp), allocatable :: el_ph_sum(:,:)

  write(6, '("Reading standard input")')
  read(5, '(a)') fni
  read(5, '(a)') fno
  read(5, '(a)') job
  read(5, *) ev, dev, nev, bvmin, bvmax
  read(5, *) ec, dec, nec, bcmin, bcmax
  read(5, *) gmax, ndist, nepmax

  write(6, '("    input file name: ", a)') trim(fni)
  write(6, '("    output file name: ", a)') trim(fno)
  write(6, '("    job type: ", a)') trim(job)
  write(s1, '(f16.8)') ev
  write(s2, '(f16.8)') dev
  write(s3, '(i8)') nev
  write(s4, '(i8)') bvmin
  write(s5, '(i8)') bvmax
  write(6, '("    ev = ", a, " eV  dev = ", a, " eV  nev = ", &
      a, "  bvmin = ", a, "  bvmax = ", a)') trim(adjustl(s1)), &
      trim(adjustl(s2)), trim(adjustl(s3)), trim(adjustl(s4)), &
      trim(adjustl(s5))
  write(s1, '(f16.8)') ec
  write(s2, '(f16.8)') dec
  write(s3, '(i8)') nec
  write(s4, '(i8)') bcmin
  write(s5, '(i8)') bcmax
  write(6, '("    ec = ", a, " eV  dec = ", a, " eV  nec = ", &
      a, "  bcmin = ", a, "  bcmax = ", a)') trim(adjustl(s1)), &
      trim(adjustl(s2)), trim(adjustl(s3)), trim(adjustl(s4)), &
      trim(adjustl(s5))
  write(s1, '(f16.8)') gmax
  write(s2, '(i8)') ndist
  write(s3, '(i8)') nepmax
  write(6, '("    gmax = ", a, " eV^2  ndist = ", a, "  nepmax = ", &
      a)') trim(adjustl(s1)), trim(adjustl(s2)), trim(adjustl(s3))

  if (trim(job) .eq. 'bpair') then
    ijob = 1
  elseif (trim(job) .eq. 'egrid') then
    ijob = 2
  elseif (trim(job) .eq. 'gdist') then
    ijob = 3
  else
    write(0, '("Error: wrong job type")')
    stop 1
  endif

  open(uni, file = fni, form = 'unformatted', status = 'old')
  write(6, '("Reading file ", a)') trim(fni)
  read(uni) s1, s2, s3
  write(6, '("    title: ", a, ", date: ", a, ", time: ", a)') &
      trim(s1), trim(s2), trim(s3)
  read(uni) ii, nat, ii, ii, ii, ii, ii, &
      ii, nspin, nbnd, nmodes, nqs
  write(s1, '(i8)') nqs
  write(s2, '(i8)') nbnd
  write(s3, '(i8)') nspin
  write(s4, '(i8)') nmodes
  write(6, '("    nqs = ", a, " nbnd = ", a, " nspin = ", a, &
      " nmodes = ", a)') trim(adjustl(s1)), trim(adjustl(s2)), &
      trim(adjustl(s3)), trim(adjustl(s4))

  if (nspin .eq. 1) then
    wspin = 1.0d0 / degspin
  else
    wspin = 1.0d0
  endif

  ee = (/ev, ec/)
  de = (/dev, dec/)
  ne = (/nev, nec/)
  nemin = minval(ne)
  nemax = maxval(ne)
  bpair = reshape((/bvmin, bvmax, bcmin, bcmax/), shape(bpair))
  do ii = 1, nwin
    if (bpair(1, ii) .lt. 1 .or. bpair(1, ii) .gt. nbnd) bpair(1, ii) = 1
    if (bpair(2, ii) .lt. 1 .or. bpair(2, ii) .gt. nbnd) bpair(2, ii) = nbnd
  enddo
  if (ijob .eq. 1) then
    allocate(gavg(nmodes, nbnd, nbnd, 1))
    allocate(gtot(1, 1, 1))
  elseif (ijob .eq. 2) then
    allocate(gavg(nmodes, nemax, nemax, nwin))
    allocate(gtot(nemax, nemax, nwin))
    allocate(gnum(nemax, nemax, nwin))
    allocate(egrid(2, nwin, nemax))
  elseif (ijob .eq. 3) then
    allocate(gdist(2, ndist, nwin))
    allocate(gepair(4, nepmax, nwin))
    allocate(egrid(2, nwin, 1))
  endif
  if (ijob .eq. 1 .or. ijob .eq. 2) then
    allocate(wavg(nmodes))
  endif
  if (ijob .eq. 2 .or. ijob .eq. 3) then
    do ii = 1, nwin
      egrid(ii, ii, 1) = ee(ii) + de(ii)
      egrid(nwin + 1 - ii, ii, 1) = sum(ee) / nwin
    enddo
    if (ijob .eq. 2 .and. nemax .gt. 1) then
      do jj = 2, nemax
        do ii = 1, nwin
          egrid(ii, ii, jj) = ee(ii) + de(ii) * jj
          egrid(nwin + 1 - ii, ii, jj) = ee(ii) + de(ii) * (jj - 1)
        enddo
      enddo
    endif
    egrid(:,:,:) = egrid(:,:,:) / rytoev
  endif

  allocate(x_q(3, nqs))
  allocate(wq(nqs))
  allocate(w2(nmodes))
  allocate(dyn(nmodes, nmodes))
  allocate(irt(48, nat))
  allocate(rtau(3, 48, nat))
  allocate(u(nmodes, nmodes))
  allocate(el_ph_sum(nmodes, nmodes))
  read(uni)
  read(uni)
  read(uni)
  read(uni)
  read(uni)
  read(uni)
  read(uni) ((at(ii, jj), ii = 1, 3), jj = 1, 3)
  read(uni) ((bg(ii, jj), ii = 1, 3), jj = 1, 3)
  read(uni)
  read(uni)
  read(uni)
  read(uni) ((x_q(ii, jj), ii = 1, 3), jj = 1, nqs)
  read(uni) (wq(ii), ii = 1, nqs)
  read(uni)
  call cryst_to_cart(nqs, x_q, bg, 1)

  if (ijob .eq. 1 .or. ijob .eq. 2) then
    gavg(:,:,:,:) = 0.0d0
    wavg(:) = 0.0d0
    gtot(:,:,:) = 0.0d0
    wtot = 0.0d0
  elseif (ijob .eq. 3) then
    gstep = gmax / ndist
    do ii = 1, nwin
      do jj = 1, ndist
        gdist(1, jj, ii) = (jj - 1) * gstep
        gdist(2, jj, ii) = 0.0d0
      enddo
    enddo
    nepair(:) = 0
  endif
  if (ijob .eq. 2) gnum(:,:,:) = 0
  do iq = 1, nqs
    xq(:) = x_q(:, iq)
    read(uni) nsymq, irotmq, ii, ii, nkstot, nksqtot
    write(s1, '(i8)') iq
    write(s2, '(i8)') nkstot
    write(s3, '(i8)') nksqtot
    write(6, '("    iq = ", a, " nkstot = ", a, " nksqtot = ", a)') &
        trim(adjustl(s1)), trim(adjustl(s2)), trim(adjustl(s3))
    allocate(wk(nkstot))
    allocate(et(nbnd, nkstot))
    allocate(ikks(nksqtot))
    allocate(ikqs(nksqtot))
    allocate(el_ph_mat(nbnd, nbnd, nksqtot, nmodes))
    read(uni) minus_q
    read(uni)
    read(uni)
    read(uni) (((rtau(ii, jj, kk), ii = 1, 3), jj = 1, 48), kk = 1, nat)
    read(uni)
    read(uni)
    read(uni) ((u(ii, jj), ii = 1, nmodes), jj = 1, nmodes)
    read(uni)
    read(uni)
    read(uni)
    read(uni)
    read(uni) (((s(ii, jj, kk), ii = 1, 3), jj = 1, 3), kk = 1, 48)
    read(uni) (invs(ii), ii = 1, 48)
    read(uni)
    read(uni)
    read(uni)
    read(uni)
    read(uni)
    read(uni) ((irt(ii, jj), ii = 1, 48), jj = 1, nat)
    read(uni)
    read(uni) (wk(ii), ii = 1, nkstot)
    read(uni) ((et(ii, jj), ii = 1, nbnd), jj = 1, nkstot)
    read(uni)
    read(uni)
    read(uni)
    read(uni) (ikks(ii), ii = 1, nksqtot)
    read(uni) (ikqs(ii), ii = 1, nksqtot)
    read(uni)
    read(uni) (w2(ii), ii = 1, nmodes)
    read(uni) ((dyn(ii, jj), ii = 1, nmodes), jj = 1, nmodes)
    read(uni) ((((el_ph_mat(ii, jj, kk, ll), ii = 1, nbnd), &
        jj = 1, nbnd), kk = 1, nksqtot), ll = 1, nmodes)

    if (ijob .eq. 1) then
      do ibnd = 1, nbnd
        do jbnd = 1, nbnd
          el_ph_sum(:,:) = (0.0d0, 0.0d0)
          do ik = 1, nksqtot
            ikk = ikks(ik)
            ikq = ikqs(ik)
            weight = wq(iq) * wk(ikk) * wspin
            if (ibnd .eq. 1 .and. jbnd .eq. 1) gtot(1, 1, 1) = &
                gtot(1, 1, 1) + weight
            do jpert = 1, nmodes
              do ipert = 1, nmodes
                el_ph_sum(ipert, jpert) = el_ph_sum(ipert, jpert) &
                    + conjg(el_ph_mat(jbnd, ibnd, ik, ipert)) &
                    * el_ph_mat(jbnd, ibnd, ik, jpert) * weight
              enddo
            enddo
          enddo
          call symdyn_munu_new(el_ph_sum, u, xq, s, invs, rtau, irt, &
              at, bg, nsymq, nat, irotmq, minus_q)
          do nu = 1, nmodes
            if (w2(nu) > epsw2) then
              factor = rytoev**2 / (2.0d0 * sqrt(w2(nu)))
            else
              factor = 0.0d0
            endif
            gbuf = 0.0d0
            do mu = 1, nmodes
              do vu = 1, nmodes
                gbuf = gbuf + dble(conjg(dyn(mu, nu)) * &
                    el_ph_sum(mu, vu) * dyn(vu, nu))
              enddo
            enddo
            gavg(nu, jbnd, ibnd, 1) = gavg(nu, jbnd, ibnd, 1) + gbuf * factor
          enddo
        enddo
      enddo
    elseif (ijob .eq. 2) then
      do ii = 1, nwin
        do jj = 1, nemax
          do kk = 1, nemax
            el_ph_sum(:,:) = (0.0d0, 0.0d0)
            do ik = 1, nksqtot
              ikk = ikks(ik)
              ikq = ikqs(ik)
              weight = wq(iq) * wk(ikk) * wspin
              do ibnd = bpair(1, ii), bpair(2, ii)
                if (et(ibnd, ikk) .gt. egrid(1, ii, jj) .and. &
                    et(ibnd, ikk) .le. egrid(2, ii, jj)) then
                  do jbnd = bpair(1, ii), bpair(2, ii)
                    if (et(jbnd, ikq) .gt. egrid(1, ii, kk) .and. &
                        et(jbnd, ikq) .le. egrid(2, ii, kk)) then
                      gnum(kk, jj, ii) = gnum(kk, jj, ii) + 1
                      gtot(kk, jj, ii) = gtot(kk, jj, ii) + weight
                      do jpert = 1, nmodes
                        do ipert = 1, nmodes
                          el_ph_sum(ipert, jpert) = el_ph_sum(ipert, jpert) &
                              + conjg(el_ph_mat(jbnd, ibnd, ik, ipert)) &
                              * el_ph_mat(jbnd, ibnd, ik, jpert) * weight
                        enddo
                      enddo
                    endif
                  enddo
                endif
              enddo
            enddo
            call symdyn_munu_new(el_ph_sum, u, xq, s, invs, rtau, irt, at, &
                bg, nsymq, nat, irotmq, minus_q)
            do nu = 1, nmodes
              if (w2(nu) > epsw2) then
                factor = rytoev**2 / (2.0d0 * sqrt(w2(nu)))
              else
                factor = 0.0d0
              endif
              gbuf = 0.0d0
              do mu = 1, nmodes
                do vu = 1, nmodes
                  gbuf = gbuf + dble(conjg(dyn(mu, nu)) * &
                      el_ph_sum(mu, vu) * dyn(vu, nu))
                enddo
              enddo
              gavg(nu, kk, jj, ii) = gavg(nu, kk, jj, ii) + gbuf * factor
            enddo
          enddo
        enddo
      enddo
    elseif (ijob .eq. 3) then
      do ii = 1, nwin
        do ik = 1, nksqtot
          ikk = ikks(ik)
          ikq = ikqs(ik)
          weight = wq(iq) * wk(ikk) * wspin
          do ibnd = bpair(1, ii), bpair(2, ii)
            if (et(ibnd, ikk) .gt. egrid(1, ii, 1) .and. &
                et(ibnd, ikk) .le. egrid(2, ii, 1)) then
              do jbnd = bpair(1, ii), bpair(2, ii)
                if (et(jbnd, ikq) .gt. egrid(1, ii, 1) .and. &
                    et(jbnd, ikq) .le. egrid(2, ii, 1)) then
                  do jpert = 1, nmodes
                    do ipert = 1, nmodes
                      el_ph_sum(ipert, jpert) = &
                          conjg(el_ph_mat(jbnd, ibnd, ik, ipert)) &
                          * el_ph_mat(jbnd, ibnd, ik, jpert)
                    enddo
                  enddo
                  call symdyn_munu_new(el_ph_sum, u, xq, s, invs, &
                      rtau, irt, at, bg, nsymq, nat, irotmq, minus_q)
                  gsum = 0.0d0
                  do nu = 1, nmodes
                    if (w2(nu) > epsw2) then
                      factor = rytoev**2 / (2.0d0 * sqrt(w2(nu)))
                    else
                      factor = 0.0d0
                    endif
                    gbuf = 0.0d0
                    do mu = 1, nmodes
                      do vu = 1, nmodes
                        gbuf = gbuf + dble(conjg(dyn(mu, nu)) * &
                            el_ph_sum(mu, vu) * dyn(vu, nu))
                      enddo
                    enddo
                    gsum = gsum + gbuf * factor
                  enddo
                  jj = nint(abs(gsum) / gstep)
                  if (jj .ge. 1 .and. jj .le. ndist) then
                    gdist(2, jj, ii) = gdist(2, jj, ii) + weight / gstep
                  endif
                  nepair(ii) = nepair(ii) + 1
                  gepair(1, nepair(ii), ii) = et(ibnd, ikk) * rytoev - ee(ii)
                  gepair(2, nepair(ii), ii) = et(jbnd, ikq) * rytoev - ee(ii)
                  gepair(3, nepair(ii), ii) = gsum
                  gepair(4, nepair(ii), ii) = weight
                endif
              enddo
            endif
          enddo
        enddo
      enddo
    endif

    if (ijob .eq. 1 .or. ijob .eq. 2) then
      weight = wq(iq)
      wtot = wtot + weight
      do nu = 1, nmodes
        if (w2(nu) > epsw2) then
          wavg(nu) = wavg(nu) + weight * sqrt(w2(nu))
        endif
      enddo
    endif

    deallocate(wk)
    deallocate(et)
    deallocate(ikks)
    deallocate(ikqs)
    deallocate(el_ph_mat)
  enddo
  close(uni, status = 'keep')

  if (ijob .eq. 1) then
    do ibnd = 1, nbnd
      do jbnd = 1, nbnd
        do nu = 1, nmodes
          gavg(nu, jbnd, ibnd, 1) = gavg(nu, jbnd, ibnd, 1) / gtot(1, 1, 1)
        enddo
      enddo
    enddo
  elseif (ijob .eq. 2) then
    do ii = 1, nwin
      do jj = 1, nemax
        do kk = 1, nemax
          if (gtot(kk, jj, ii) .gt. eps12) then
            do nu = 1, nmodes
              gavg(nu, kk, jj, ii) = gavg(nu, kk, jj, ii) / gtot(kk, jj, ii)
            enddo
          endif
        enddo
      enddo
    enddo
  endif

  if (ijob .eq. 1 .or. ijob .eq. 2) then
    do nu = 1, nmodes
      wavg(nu) = wavg(nu) * ry_to_cmm1 / wtot
    enddo
  endif

  if (ijob .eq. 2) then
    write(fmt, '("(", i0, "(5x, ""e"", a, "" (eV)"", 4x, ""e"", a, ""''", &
        " (eV)"", 3x, ""<|g"", a, ""|^2> (eV^2)"", 2x, ""count"", a, 5x,", &
        " ""weight"", a))")') nwin
    write(6, fmt) ('v', ll = 1, 5), ('c', ll = 1, 5)
    write(fmt, '("(", i0, "(2f12.6, e18.8, i8, f12.6))")') nwin
    do jj = 1, nemin
      do kk = 1, nemin
        write(6, fmt) ((jj - 0.5d0) * de(ii), (kk - 0.5d0) * de(ii), &
            sum(gavg(:, kk, jj, ii)), gnum(kk, jj, ii), gtot(kk, jj, ii), &
            ii = 1, nwin)
      enddo
    enddo
  endif

  open(uno, file = fno, form = 'formatted', status = 'replace')
  write(6, '("Writing file ", a)') trim(fno)
  if (ijob .eq. 1) then
    write(uno, '("    mode       <w> (cm^-1)")')
    write(uno, '("--------------------------")')
    do nu = 1, nmodes
      write(uno, '(i8, e18.8)') nu, wavg(nu)
    enddo
    write(uno, *)
    write(uno, '("    band    band    mode    <|g|^2> (eV^2)")')
    write(uno, '("------------------------------------------")')
    do ibnd = 1, nbnd
      do jbnd = 1, nbnd
        do nu = 1, nmodes
          write(uno, '(3i8, e18.8)') ibnd, jbnd, nu, gavg(nu, jbnd, ibnd, 1)
        enddo
      enddo
    enddo
  elseif (ijob .eq. 2) then
    write(uno, '(2i8)') nwin, nmodes
    do ii = 1, nwin
      write(uno, '(2f14.8, i8)') ee(ii), de(ii), ne(ii)
    enddo
    write(fmt, '("(", i0, "e18.8)")') nmodes
    write(uno, fmt) (wavg(nu), nu = 1, nmodes)
    write(fmt, '("(3i8, ", i0, "e18.8)")') nmodes
    do ii = 1, nwin
      do jj = 1, ne(ii)
        do kk = 1, ne(ii)
          write(uno, fmt) ii, jj, kk, (gavg(nu, kk, jj, ii), &
              nu = 1, nmodes)
        enddo
      enddo
    enddo
  elseif (ijob .eq. 3) then
    write(uno, '(i8)') ndist
    do jj = 1, ndist
      write(uno, '(f14.8, 2e18.8)') gdist(1, jj, 1), gdist(2, jj, 1), &
          gdist(2, jj, 2)
    enddo
    write(uno, '(i8)') (nepair(ii), ii = 1, nwin)
    do ii = 1, nwin
      do jj = 1, nepair(ii)
        write(uno, '(2f14.8, 2e18.8)') (gepair(kk, jj, ii), kk = 1, 4)
      enddo
    enddo
  endif
  close(uno, status = 'keep')

  if (ijob .eq. 1 .or. ijob .eq. 2) then
    deallocate(gavg)
    deallocate(gtot)
    deallocate(wavg)
  endif
  if (ijob .eq. 2) deallocate(gnum)
  if (ijob .eq. 3) then
    deallocate(gdist)
    deallocate(gepair)
  endif
  if (ijob .eq. 2 .or. ijob .eq. 3) deallocate(egrid)
  deallocate(x_q)
  deallocate(wq)
  deallocate(w2)
  deallocate(dyn)
  deallocate(irt)
  deallocate(rtau)
  deallocate(u)
  deallocate(el_ph_sum)

end program epa

