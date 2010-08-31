!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine transmit(ik, ien, tk_out, left_to_right)

!
! This subroutine constructs the scattering states
! using the CBS of the left and right tips (saved by compbs)
! and the functions and integrals computed by scatter_forw in
! the scattering region.
!
  use io_global,  ONLY :  stdout
  use cond_files,  ONLY :  prefixl, prefixs
  use lsda_mod, only: nspin
  use noncollin_module, ONLY : noncolin, npol
  use spin_orb, only : lspinorb
 use cond
implicit none

  integer, intent(in) :: ik, ien
  real(DP), intent(out) :: tk_out
  !
  integer :: n, iorb, iorb1, iorb2, iorba, ipol, nt, ir, it, &
             ig, ntran, ij, is, js, nchan_in, nchan_out, info
  integer, allocatable :: ipiv(:)
  real(DP) :: tk, tj, rj, tij, rij, eev
  real(DP), allocatable :: zps(:,:), eigen(:)
  complex(DP) :: x1, x2, xi1(2)
  complex(DP), allocatable :: amat(:,:), vec1(:,:), &
                     tmat(:,:), veceig(:,:), zps_nc(:,:), &
                     vec2(:,:), smat(:,:)
  LOGICAL :: left_to_right ! if .t./.f. solves the scatt. problem
  !          for right/left moving states

!--
! Set up the number of incoming and outgoing channels
!
  if (left_to_right) then
    nchan_in = nchanl
    nchan_out = nchanr
  else
    nchan_in = nchanr
    nchan_out = nchanl
  endif
!--

! electron scattering energy in eV with respect to the Fermi level
  eev = earr(ien)

!--
! Goes further only if nchan_in, nchan_out <> 0 or
! nchan_in <> 0 and nchan_out = 0 but lorb = .t.
!
  if (nchan_in*nchan_out.eq.0) then
    tk = 0.d0
    tk_out = tk
    WRITE( stdout,'(a24, 2f12.7)') 'E-Ef(ev), T = ',   &
                           eev, tk
    if (nchan_in.eq.0.or..not.lorb) return
  endif
!--

  if(lorb) then
    write(stdout,*)
    if(left_to_right) then
      write(stdout,*) 'RIGHT MOVING --> scattering states ...'
    else
      write(stdout,*) 'LEFT MOVING <-- scattering states ...'
    endif
    write(stdout,*)
  endif

  ntran=4*n2d+npol*(norbs+nocrosl+nocrosr)

  allocate( ipiv( ntran ) )
  allocate( amat( ntran, ntran ) )
  allocate( vec1( ntran, nchan_in ) )
  allocate( smat(nchan_in+nchan_out, nchan_in) )
  if (nchan_out.ne.0) allocate( tmat( nchan_out, nchan_in ) )
  allocate( veceig( nchan_in, nchan_in ) )
  allocate( eigen( nchan_in ) )

  if (noncolin) then
    allocate( zps_nc( norbs*npol, norbs*npol ) )
  else
    allocate( zps( norbs, norbs ) )
  endif
  amat=(0.d0,0.d0)
!
! To form  zps=zpseu-e*qq
!
  do iorb=1, norbs
    do iorb1=1, norbs
      if (noncolin) then
        ij=0
        do is=1,npol
          do js=1,npol
            ij=ij+1
            zps_nc(npol*(iorb-1)+is, npol*(iorb1-1)+js)=         &
                zpseus_nc(1,iorb,iorb1,ij)
            if (lspinorb) then
              zps_nc(npol*(iorb-1)+is,npol*(iorb1-1)+js)=        &
                    zps_nc(npol*(iorb-1)+is,npol*(iorb1-1)+js) &
                            -eryd*zpseus_nc(2,iorb,iorb1,ij)
            else
              if (is.eq.js)                                      &
                  zps_nc(npol*(iorb-1)+is,npol*(iorb1-1)+js)=    &
                      zps_nc(npol*(iorb-1)+is,npol*(iorb1-1)+js) &
                             -eryd*zpseus_nc(2,iorb,iorb1,ij)
            endif
          enddo
        enddo
      else
        zps(iorb,iorb1)=zpseus(1,iorb,iorb1)-eryd*zpseus(2,iorb,iorb1)
      endif
    enddo
  enddo

!
! Compute the part of amat which comes from the matching of
! the wave function on the boundaries
!
!     1) local functions
  do n=1, n2d
    do ig=1, n2d
      amat(ig,n)=fun0(ig,n)
      amat(ig+n2d,n)=fund0(ig,n)
      amat(ig+2*n2d,n)=fun1(ig,n)
      amat(ig+3*n2d,n)=fund1(ig,n)
      amat(ig,n+n2d)=fun0(ig,n2d+n)
      amat(ig+n2d,n+n2d)=fund0(ig,n2d+n)
      amat(ig+2*n2d,n+n2d)=fun1(ig,n2d+n)
      amat(ig+3*n2d,n+n2d)=fund1(ig,n2d+n)
    enddo
  enddo
!     2) nonlocal functions
  do iorb=1, norbs*npol
    do ig=1, n2d
      amat(ig,2*n2d+iorb)=funl0(ig, iorb)
      amat(ig+n2d,2*n2d+iorb)=fundl0(ig, iorb)
      amat(ig+2*n2d,2*n2d+iorb)=funl1(ig, iorb)
      amat(ig+3*n2d,2*n2d+iorb)=fundl1(ig, iorb)
    enddo
  enddo
!     4) to add reflection and transmission parts
  do ig=1, n2d
    do n=1, n2d+npol*nocrosl
      amat(ig,2*n2d+npol*norbs+n)=-kfunl(ig,n2d+npol*nocrosl+n)
      amat(ig+n2d,2*n2d+npol*norbs+n)=-kfundl(ig,n2d+npol*nocrosl+n)
    enddo
    do n=1, n2d+npol*nocrosr
      amat(ig+2*n2d,3*n2d+npol*(norbs+nocrosl)+n)=-kfunr(ig,n)
      amat(ig+3*n2d,3*n2d+npol*(norbs+nocrosl)+n)=-kfundr(ig,n)
    enddo
  enddo
!
!     3) Part coming from the definion of C_alpha
!
  do iorb=1, norbs*npol
    do n=1, 2*n2d
      do iorb1=1, norbs*npol
        if (noncolin) then
           amat(4*n2d+iorb,n)=amat(4*n2d+iorb,n)-             &
                    zps_nc(iorb,iorb1)*intw1(iorb1,n)
        else
           amat(4*n2d+iorb,n)=amat(4*n2d+iorb,n)-             &
                    zps(iorb,iorb1)*intw1(iorb1,n)
        endif
      enddo
    enddo
    do iorb1=1, norbs*npol
      do iorb2=1, norbs*npol
        if (noncolin) then
           amat(4*n2d+iorb,2*n2d+iorb1)=amat(4*n2d+iorb,2*n2d+iorb1)-     &
                     zps_nc(iorb,iorb2)*intw2(iorb2, iorb1)
        else
           amat(4*n2d+iorb,2*n2d+iorb1)=amat(4*n2d+iorb,2*n2d+iorb1)-     &
                     zps(iorb,iorb2)*intw2(iorb2, iorb1)
        endif
      enddo
    enddo
    amat(4*n2d+iorb,2*n2d+iorb)=amat(4*n2d+iorb,2*n2d+iorb)+(1.d0,0.d0)
  enddo
!     To set up terms depending on kint and kcoef of
!     5) left tip
  do ig=1, nocrosl*npol
    do n=1, n2d+nocrosl*npol
      amat(4*n2d+ig,2*n2d+npol*norbs+n)=-kintl(ig,n2d+npol*nocrosl+n)
      amat(4*n2d+npol*norbs+ig,2*n2d+npol*norbs+n)= &
                               -kcoefl(ig,n2d+npol*nocrosl+n)
    enddo
  enddo
!     6) right tip
  do ig=1, nocrosr*npol
    do n=1, n2d+nocrosr*npol
      amat(4*n2d+npol*(nocrosl+noinss)+ig,3*n2d+npol*(norbs+nocrosl)+n)= &
                         -kintr(ig,n)
      amat(4*n2d+npol*(norbs+nocrosl)+ig,3*n2d+npol*(norbs+nocrosl)+n)= &
                         -kcoefr(ig,n)
    enddo
  enddo
!     7) to match C_alpha for crossing orbitals with the tips ones
  do ig=1, nocrosl*npol
    amat(4*n2d+norbs*npol+ig,2*n2d+ig)=(1.d0,0.d0)
  enddo
  do ig=1, nocrosr*npol
    amat(4*n2d+npol*(norbs+nocrosl)+ig,2*n2d+npol*(nocrosl+noinss)+ig)= &
           (1.d0,0.d0)
  enddo

!--
! Form the vector of free coefficients
!
  vec1 = (0.d0,0.d0)
  if (left_to_right) then
    do n = 1, nchan_in
      do ig = 1, n2d
        vec1(ig,n) = kfunl(ig,n)
        vec1(n2d+ig,n) = kfundl(ig,n)
      enddo
      do ig = 1, nocrosl*npol
        vec1(4*n2d+ig,n) = kintl(ig,n)
        vec1(4*n2d+npol*norbs+ig,n) = kcoefl(ig,n)
      enddo
    enddo
  else
    do n = 1, nchan_in
      do ig = 1, n2d
        vec1(2*n2d+ig,n) = kfunr(ig,n2d+npol*nocrosr+n)
        vec1(3*n2d+ig,n) = kfundr(ig,n2d+npol*nocrosr+n)
      enddo
      do ig = 1, nocrosr*npol
        vec1(4*n2d+npol*(norbs-nocrosr)+ig,n) = kintr(ig,n2d+npol*nocrosr+n)
        vec1(4*n2d+npol*(norbs+nocrosl)+ig,n) = kcoefr(ig,n2d+npol*nocrosr+n)
      enddo
    enddo
  endif
!--

!
! Solve the system on the coefficiens vec1 of scattering states
!
  call ZGESV(ntran, nchan_in, amat, ntran, ipiv, vec1, ntran, info)

  if (info.ne.0) call errore('transmit','problems with the linear system', &
                                                              abs(info))
!--
! transmission matrix tmat n --> ig and
! the S-matrix, (r_{ig,n}, t_{ig,n})
!
  if (left_to_right) then
    ir = 2*n2d + npol*norbs
    it = ir + n2d + npol*nocrosl
  else
    it = 2*n2d + npol*norbs
    ir = it + n2d + npol*nocrosl
  endif

  do n = 1, nchan_in
    do ig = 1, nchan_out
      tmat(ig,n) = vec1(it+ig,n)
    enddo
  enddo
  do n = 1, nchan_in
    do ig = 1, nchan_in
      smat(ig,n) = vec1(ir+ig,n)
    enddo
    do ig = 1, nchan_out
      smat(nchan_in+ig,n) = tmat(ig,n)
    enddo
  enddo
!--

!--
! Check for the S matrix unitarity
!
  if (nchan_out.ne.0) then
    call sunitary(nchan_in, nchan_out, smat, info)
    call errore('transmit','S matrix is not unitary',-abs(info))
  endif
!--

!--
! transmission and reflection coefficients of each band
!
  WRITE( stdout,*) 'Band j to band i transmissions and reflections:'
  WRITE( stdout,'(4x,''j'',9x,''i'',5x,''|T_ij|^2'',4x,''|R_ij|^2'')')
  WRITE( stdout,*)
  ij = MIN(nchan_in, nchan_out)
  do n = 1, nchan_in
    tj = 0.d0
    rj = 0.d0
    do ig = 1, ij
     tij = DBLE(smat(nchan_in+ig,n))**2+AIMAG(smat(nchan_in+ig,n))**2
     tj = tj+tij
     rij = DBLE(smat(ig,n))**2+AIMAG(smat(ig,n))**2
     rj = rj+rij
     WRITE(stdout,'(i5,'' --> '',i5,2f12.5)') n, ig, tij, rij
    enddo
    do ig = ij + 1, nchan_out
     tij = DBLE(smat(nchan_in+ig,n))**2+AIMAG(smat(nchan_in+ig,n))**2
     tj = tj+tij
     WRITE(stdout,'(i5,'' --> '',i5,f12.5)') n, ig, tij
    enddo
    do ig = ij + 1, nchan_in
     rij = DBLE(smat(ig,n))**2+AIMAG(smat(ig,n))**2
     rj = rj+rij
     WRITE(stdout,'(i5,'' --> '',i5,12x,f12.5)') n, ig, rij
    enddo
    WRITE(stdout,'(3x,''Total T_j, R_j =  '',2f9.5)') tj, rj
    WRITE(stdout,*)
  enddo
!--

!--------------
! eigenchannel decomposition
!
  tk = 0
  if (nchan_out.ne.0) then
    call eigenchnl(nchan_in, nchan_out, tmat, veceig, eigen)

!-- calculate and output of T(k) on a general file
    if (left_to_right) then
      do n = 1, nchan_in
        tk = tk + eigen(n)
      enddo
      if (nspin.eq.1) then
        tk = 2.d0*tk
        WRITE(stdout,'(a24, 2f12.7)') 'E-Ef(ev), T(x2 spins) = ', eev, tk
      else
        WRITE(stdout,'(a24, 2f12.7)') 'E-Ef(ev), T = ', eev, tk
      endif
    endif
!--

    if (prefixl.ne.prefixs) then
      WRITE( stdout,*) 'Eigenchannel decomposition:'
      do n = 1, nchan_in
        WRITE( stdout,'(''#'',i5, 2f9.5)') n, eev, eigen(n)
        do ig = 1, nchan_in
          tj = DBLE(veceig(ig,n))**2+AIMAG(veceig(ig,n))**2
          WRITE( stdout,'(20x, f9.5)') tj
        enddo
      enddo
    endif
  else
    veceig(:,:) = 0.d0
    do n = 1, nchan_in
      veceig(n,n) = 1.d0
    enddo
  endif
!--------------

!
! output T(k), to be added to the total T
!
   tk_out = tk

!---------------------------
!   Angular momentum projection of transmission
!
  if(left_to_right.and.(orbj_in*orbj_fin).gt.0) then
    nt = orbj_fin - orbj_in + 1
    allocate( vec2( ntran, nchanl ) )
    x1 = (1.d0,0.d0)
    x2 = (0.d0,0.d0)
    call zgemm('n', 'n', ntran, nchanl, nchanl, x1, vec1, ntran,  &
              veceig, nchanl, x2, vec2, ntran)
    write(stdout,*) 'Nchannel, Norbital, projection'
!---------------------------
!   Angular momentum projection of eigenchannels
!
    do n = 1, nchanl
     if(eigen(n).gt.1.d-5) then
      do iorb = orbj_in, orbj_fin
         do ipol=1, npol
            iorba=(iorb-1)*npol+ipol
            xi1(ipol) = 0.d0
            do ig = 1, 2*n2d
               xi1(ipol) = xi1(ipol)+intw1(iorba, ig)*vec2(ig, n)
            enddo
            do ig = 1, norbs*npol
               xi1(ipol) = xi1(ipol)+intw2(iorba, ig)*vec2(2*n2d+ig, n)
            enddo
         enddo
         write(stdout,'(2i5,2f16.12,2x,2f16.12)') n, iorb-orbj_in+1,   &
                                      (xi1(ipol),ipol=1,npol)

      enddo
     endif
    enddo
!------------------------
    deallocate(vec2)

  endif

!--
! Constructs and writes the transmission eigenchannels
!
  if (lorb) then
     allocate( vec2(ntran,nchan_in) )
     x1 = (1.d0,0.d0)
     x2 = (0.d0,0.d0)
     call zgemm('n', 'n', ntran, nchan_in, nchan_in, x1, vec1, ntran,  &
              veceig, nchan_in, x2, vec2, ntran)

     call scat_states_plot(ik,ien,norbs,nocrosl,nchan_in,vec2,veceig, &
                           left_to_right)
     deallocate( vec2 )
  endif
!--

  deallocate(ipiv)
  deallocate(amat)
  deallocate(smat)
  if (nchan_out.ne.0) deallocate(tmat)
  deallocate(eigen)
  if (noncolin) then
    deallocate(zps_nc)
  else
    deallocate(zps)
  endif


  deallocate(vec1)
  deallocate(veceig)

  return
end subroutine transmit

