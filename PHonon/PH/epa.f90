!
! Copyright (C) 2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
program epa
  !----------------------------------------------------------------------------
  !
  !! epa.x:
  !!    reads electron-phonon coupling matrix elements produced by the
  !!    phonon code with electron_phonon = 'epa', and makes a transformation
  !!    from momentum to energy space according to electron-phonon averaged
  !!    (EPA) approximation as described in G. Samsonidze & B. Kozinsky,
  !!    Adv. Energy Mater. 2018, 1800246 doi:10.1002/aenm.201800246
  !!    arXiv:1511.08115
  !!
  !! Input data:
  !!
  !!    fin
  !!    fout
  !!    job
  !!    edgev, stepv, nbinv, ivmin, ivmax
  !!    edgec, stepc, nbinc, icmin, icmax
  !!    stepg, nbing, ngmax
  !!
  !!    fin        :  input file produced by the phonon code
  !!    fout       :  output file used by the BoltzTraP code
  !!    job        :  job type determines how to average the matrix elements
  !!                  - 'egrid' is to average within bins of regular energy
  !!                            grids (standard procedure)
  !!                  - 'bpair' is to average over individual pairs of bands
  !!                            (obsolete procedure)
  !!                  - 'ghist' is to construct the histogram for plotting
  !!                            a distribution of squared matrix elements
  !!    edgev      :  grid edge of the valence energy grid (in eV), set to
  !!                  the valence band maximum
  !!    stepv      :  grid step of the valence energy grid (in eV), set to
  !!                  a small negative value
  !!    nbinv      :  number of bins in the valence energy grid
  !!    ivmin/max  :  range of valence bands (0/0 for all valence bands)
  !!    edgec      :  grid edge of the conduction energy grid (in eV), set to
  !!                  the conduction band minimum
  !!    stepc      :  grid step of the conduction energy grid (in eV), set to
  !!                  a small positive value
  !!    nbinc      :  number of bins in the conduction energy grid
  !!    icmin/max  :  range of conduction bands (0/0 for all conduction bands)
  !!    stepg      :  grid step for the histogram (in eV^2)
  !!    nbing      :  number of bins for the histogram
  !!    ngmax      :  maximum number of elements in one bin of the histogram
  !!
  !! For more details on how to set up the energy grids please refer to
  !! online documentation at https://github.com/mir-group/EPA
  !
  use kinds, only : dp
  use constants, only : ry_to_cmm1, rytoev, degspin, eps12
  implicit none

  integer, parameter :: ngrid = 2
  integer, parameter :: uni = 101
  integer, parameter :: uno = 102
  real(dp), parameter :: epsw2 = (20_dp / ry_to_cmm1)**2

  logical :: minus_q
  integer :: nmodes, nqs, nspin, nbnd, nkstot, nksqtot, &
      nat, nsymq, irotmq, iq, ik, ikk, ikq, ibnd, jbnd, &
      nu, mu, vu, ipert, jpert, ii, jj, kk, ll, ijob, &
      nbinv, ivmin, ivmax, nbinc, icmin, icmax, nbing, ngmax, &
      nbinmax, nbin(ngrid), nhist(ngrid), bpair(2, ngrid), &
      s(3, 3, 48), invs(48)
  real(dp) :: wtot, weight, factor, gbuf, gsum, &
      stepg, edgev, stepv, edgec, stepc, wspin, &
      edge(ngrid), step(ngrid), xq(3), at(3, 3), bg(3, 3)
  character(len=256) :: fin, fout, job, fmt
  character(len=32) :: s1, s2, s3, s4, s5

  real(dp), allocatable :: ghist(:,:,:)
  real(dp), allocatable :: gdist(:,:,:)
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
  read(5, '(a)') fin
  read(5, '(a)') fout
  read(5, '(a)') job
  read(5, *) edgev, stepv, nbinv, ivmin, ivmax
  read(5, *) edgec, stepc, nbinc, icmin, icmax
  read(5, *) stepg, nbing, ngmax

  write(6, '("    input file name: ", a)') trim(fin)
  write(6, '("    output file name: ", a)') trim(fout)
  write(6, '("    job type: ", a)') trim(job)
  write(s1, '(f16.8)') edgev
  write(s2, '(f16.8)') stepv
  write(s3, '(i8)') nbinv
  write(s4, '(i8)') ivmin
  write(s5, '(i8)') ivmax
  write(6, '("    edgev = ", a, " eV  stepv = ", a, " eV  nbinv = ", &
      a, "  ivmin = ", a, "  ivmax = ", a)') trim(adjustl(s1)), &
      trim(adjustl(s2)), trim(adjustl(s3)), trim(adjustl(s4)), &
      trim(adjustl(s5))
  write(s1, '(f16.8)') edgec
  write(s2, '(f16.8)') stepc
  write(s3, '(i8)') nbinc
  write(s4, '(i8)') icmin
  write(s5, '(i8)') icmax
  write(6, '("    edgec = ", a, " eV  stepc = ", a, " eV  nbinc = ", &
      a, "  icmin = ", a, "  icmax = ", a)') trim(adjustl(s1)), &
      trim(adjustl(s2)), trim(adjustl(s3)), trim(adjustl(s4)), &
      trim(adjustl(s5))
  write(s1, '(f16.8)') stepg
  write(s2, '(i8)') nbing
  write(s3, '(i8)') ngmax
  write(6, '("    stepg = ", a, " eV^2  nbing = ", a, "  ngmax = ", &
      a)') trim(adjustl(s1)), trim(adjustl(s2)), trim(adjustl(s3))

  if (trim(job) .eq. 'bpair') then
    ijob = 1
  elseif (trim(job) .eq. 'egrid') then
    ijob = 2
  elseif (trim(job) .eq. 'ghist') then
    ijob = 3
  else
    write(0, '("Error: wrong job type")')
    stop 1
  endif

  open(uni, file = fin, form = 'unformatted', status = 'old')
  write(6, '("Reading file ", a)') trim(fin)
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

  edge = (/edgev, edgec/)
  step = (/stepv, stepc/)
  nbin = (/nbinv, nbinc/)
  nbinmax = maxval(nbin)
  bpair = reshape((/ivmin, ivmax, icmin, icmax/), shape(bpair))
  do ii = 1, ngrid
    if (bpair(1, ii) .lt. 1 .or. bpair(1, ii) .gt. nbnd) bpair(1, ii) = 1
    if (bpair(2, ii) .lt. 1 .or. bpair(2, ii) .gt. nbnd) bpair(2, ii) = nbnd
  enddo
  if (ijob .eq. 1) then
    allocate(gavg(nmodes, nbnd, nbnd, 1))
    allocate(gtot(1, 1, 1))
  elseif (ijob .eq. 2) then
    allocate(gavg(nmodes, nbinmax, nbinmax, ngrid))
    allocate(gtot(nbinmax, nbinmax, ngrid))
    allocate(gnum(nbinmax, nbinmax, ngrid))
    allocate(egrid(2, ngrid, nbinmax))
  elseif (ijob .eq. 3) then
    allocate(ghist(2, nbing, ngrid))
    allocate(gdist(4, ngmax, ngrid))
    allocate(egrid(2, ngrid, 1))
  endif
  if (ijob .eq. 1 .or. ijob .eq. 2) then
    allocate(wavg(nmodes))
  endif
  if (ijob .eq. 2 .or. ijob .eq. 3) then
    do ii = 1, ngrid
      egrid(ii, ii, 1) = edge(ii) + step(ii)
      egrid(ngrid + 1 - ii, ii, 1) = sum(edge) / ngrid
    enddo
    if (ijob .eq. 2 .and. nbinmax .gt. 1) then
      do jj = 2, nbinmax
        do ii = 1, ngrid
          egrid(ii, ii, jj) = edge(ii) + step(ii) * jj
          egrid(ngrid + 1 - ii, ii, jj) = edge(ii) + step(ii) * (jj - 1)
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
    do ii = 1, ngrid
      do jj = 1, nbing
        ghist(1, jj, ii) = (jj - 1) * stepg
        ghist(2, jj, ii) = 0.0d0
      enddo
    enddo
    nhist(:) = 0
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
      do ii = 1, ngrid
        do jj = 1, nbin(ii)
          do kk = 1, nbin(ii)
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
      do ii = 1, ngrid
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
                  jj = nint(abs(gsum) / stepg)
                  if (jj .ge. 1 .and. jj .le. nbing) then
                    ghist(2, jj, ii) = ghist(2, jj, ii) + weight / stepg
                  endif
                  nhist(ii) = nhist(ii) + 1
                  gdist(1, nhist(ii), ii) = et(ibnd, ikk) * rytoev - edge(ii)
                  gdist(2, nhist(ii), ii) = et(jbnd, ikq) * rytoev - edge(ii)
                  gdist(3, nhist(ii), ii) = gsum
                  gdist(4, nhist(ii), ii) = weight
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
    do ii = 1, ngrid
      do jj = 1, nbin(ii)
        do kk = 1, nbin(ii)
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
        " ""weight"", a))")') ngrid
    write(6, fmt) ('v', ll = 1, 5), ('c', ll = 1, 5)
    write(fmt, '("(", i0, "(2f12.6, e18.8, i8, f12.6))")') ngrid
    do jj = 1, nbinmax
      do kk = 1, nbinmax
        write(6, fmt) ((jj - 0.5d0) * step(ii), (kk - 0.5d0) * step(ii), &
            sum(gavg(:, kk, jj, ii)), gnum(kk, jj, ii), gtot(kk, jj, ii), &
            ii = 1, ngrid)
      enddo
    enddo
  endif

  open(uno, file = fout, form = 'formatted', status = 'replace')
  write(6, '("Writing file ", a)') trim(fout)
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
    write(uno, '(2i8)') ngrid, nmodes
    do ii = 1, ngrid
      write(uno, '(2f14.8, i8)') edge(ii), step(ii), nbin(ii)
    enddo
    write(fmt, '("(", i0, "e18.8)")') nmodes
    write(uno, fmt) (wavg(nu), nu = 1, nmodes)
    write(fmt, '("(3i8, ", i0, "e18.8)")') nmodes
    do ii = 1, ngrid
      do jj = 1, nbin(ii)
        do kk = 1, nbin(ii)
          write(uno, fmt) ii, jj, kk, (gavg(nu, kk, jj, ii), &
              nu = 1, nmodes)
        enddo
      enddo
    enddo
  elseif (ijob .eq. 3) then
    write(uno, '(i8)') nbing
    do jj = 1, nbing
      write(uno, '(f14.8, 2e18.8)') ghist(1, jj, 1), ghist(2, jj, 1), &
          ghist(2, jj, 2)
    enddo
    write(uno, '(i8)') (nhist(ii), ii = 1, ngrid)
    do ii = 1, ngrid
      do jj = 1, nhist(ii)
        write(uno, '(2f14.8, 2e18.8)') (gdist(kk, jj, ii), kk = 1, 4)
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
    deallocate(ghist)
    deallocate(gdist)
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
!
!----------------------------------------------------------------------------
