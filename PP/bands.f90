!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
program bands
  !-----------------------------------------------------------------------
#include "machine.h"
  use io_files, only: nd_nmbr, prefix, tmp_dir
#ifdef __PARA
  use para, only: me, npool
  use io_global, only: ionode_id
  use mp, only: mp_bcast
#endif
  implicit none
  !
  character (len=80) :: filband
  character(len=256) :: outdir
  integer :: ios
  namelist / inputpp / outdir, prefix, filband
                                                                                
  CHARACTER (LEN=80)  :: input_file
  INTEGER             :: nargs, iiarg, ierr, ilen
  INTEGER, EXTERNAL   :: iargc

  !
  call start_postproc (nd_nmbr)
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  outdir = './'
  filband = 'bands.out'
  !
#ifdef __PARA
  if (npool /= 1) call errore('bands','pools not implemented',npool)
  if (me == 1)  then
#endif
  !
  ! ... Input from file ?
  !
  nargs = iargc()
  !
  DO iiarg = 1, ( nargs - 1 )
     !
     CALL getarg( iiarg, input_file )
     IF ( TRIM( input_file ) == '-input' .OR. &
          TRIM( input_file ) == '-inp'   .OR. &
          TRIM( input_file ) == '-in' ) THEN
        !
        CALL getarg( ( iiarg + 1 ) , input_file )
        OPEN ( UNIT = 5, FILE = input_file, FORM = 'FORMATTED', &
               STATUS = 'OLD', IOSTAT = ierr )
        CALL errore( 'iosys', 'input file ' // TRIM( input_file ) // &
                   & ' not found' , ierr )
        !
     END IF
     !
  END DO

  read (5, inputpp, err = 200, iostat = ios)
200 call errore ('do_bands', 'reading inputpp namelist', abs (ios) )
  !
  tmp_dir = trim(outdir)
  !
#ifdef __PARA
  end if
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( filband, ionode_id )
#endif
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  call read_file
  call openfil_pp
  call init_us_1
  !
  call punch_band (filband)
  !
  call stop_pp
  stop
end program bands
!
!-----------------------------------------------------------------------
subroutine punch_band (filband)
  !-----------------------------------------------------------------------
  !
  !    This routine writes the band energies on a file. The routine orders
  !    the eigenvalues using the overlap of the eigenvectors to give
  !    an estimate crossing and anticrossing of the bands. This simplified
  !    method works in many, but not in all the cases.
  !
  !
#ifdef __PARA
  use para, only: me
#endif
  use atom
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  use cell_base
  use constants, only: rytoev
  use gvect
  use lsda_mod, only: nspin
  use klist
  use io_files, only: iunpun, nwordwfc, iunwfc
  use wvfct
  use uspp, only: nkb, vkb, qq
  use uspp_param, only: tvanp, nh, nhm
  use wavefunctions_module, only: evc

  implicit none
  character (len=*) :: filband
  real(kind=DP) :: proold
  ! the best overlap product
  complex(kind=DP) :: pro
  ! the product of wavefunctions

  complex(kind=DP), allocatable :: psiold (:,:), old (:), new (:)
  ! psiold: eigenfunctions at previous k-point, ordered
  ! old, new: contain one band resp. at previous and current k-point
  complex(kind=DP), allocatable :: becp(:,:), becpold (:,:)
  ! becp   : <psi|beta> at current  k-point
  ! becpold: <psi|beta> at previous k-point

  integer :: ibnd, jbnd, ik, ikb, ig, npwold, ios
  ! counters
  integer, allocatable :: ok (:), igkold (:), il (:)
  ! ok: keeps track of which bands have been already ordered
  ! igkold: indices of k+G at previous k-point
  ! il: band ordering
  integer, parameter :: maxdeg = 4
  ! maxdeg : max allowed degeneracy
  integer :: ndeg, deg, nd
  ! ndeg : number of degenerate states
  integer, allocatable :: degeneracy(:), degbands(:,:), index(:)
  ! degbands keeps track of which states are degenerate
  real(kind=DP), allocatable:: edeg(:)
  real(kind=DP), parameter :: eps = 0.001
  ! threshold (Ry) for degenerate states 
  complex(kind=DP), external :: cgracsc
 ! scalar product with the S matrix

  if (filband == ' ') return
  if (nspin == 2) call errore('bands','LSDA bands not implemented',-nspin)
  iunpun = 18
#ifdef __PARA
  if (me == 1) then
#endif
  open (unit = iunpun, file = filband, status = 'unknown', form = &
       'formatted', err = 100, iostat = ios)
100 call errore ('punch_band', 'Opening filband file', abs (ios) )
  rewind (iunpun)
#ifdef __PARA
  endif
#endif
  !
  allocate (psiold( npwx, nbnd))    
  allocate (old(ngm), new(ngm))    
  allocate (becp(nkb, nbnd), becpold(nkb, nbnd))    
  allocate (igkold (npwx))    
  allocate (ok (nbnd), il (nbnd))    
  allocate (degeneracy(nbnd), edeg(nbnd))
  allocate (index(maxdeg), degbands(nbnd,maxdeg))
  !
  do ik = 1, nks
     !
     !    prepare the indices of this k point
     !
     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, &
          igk, g2kin)
     !
     !   read eigenfunctions
     !
     call davcio (evc, nwordwfc, iunwfc, ik, - 1)
     !
     ! calculate becp = <psi|beta> 
     ! 
     call init_us_2 (npw, igk, xk (1, ik), vkb)
     call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
     !
     if (ik == 1) then
        !
        !  first k-point in the list:
        !  save eigenfunctions in the current order (increasing energy)
        !
        do ibnd = 1, nbnd
           il (ibnd) = ibnd
        end do
     else
        !
        !  following  k-points in the list:
        !  determine eigenfunction order in array il
        !
        do ibnd = 1, nbnd
           ok (ibnd) = 0
        enddo
        do ibnd = 1, nbnd
           old(:) = (0.d0, 0.d0)
           do ig = 1, npwold
              old (igkold (ig) ) = psiold (ig, ibnd)
           enddo
           proold = 0.d0
           do jbnd = 1, nbnd
              if (ok (jbnd) == 0) then
                 new (:) = (0.d0, 0.d0)
                 do ig = 1, npw
                    new (igk (ig) ) = evc (ig, jbnd)
                 enddo
                 pro = cgracsc (nkb, becp (1, jbnd), becpold (1, ibnd), &
                      nhm, ntyp, nh, qq, nat, ityp, ngm, new, old, tvanp)
                 if (abs (pro) > proold) then
                    il (ibnd) = jbnd
                    proold = abs (pro)
                 endif
              endif
           enddo
           ok (il (ibnd) ) = 1
        enddo
        !
        !  if there were bands crossing at degenerate eigenvalues
        !  at previous k-point, re-order those bands so as to keep
        !  lower band indices corresponding to lower bands
        !
        do nd = 1, ndeg
           do deg = 1, degeneracy (nd)
              index(deg) = il(degbands(nd,deg))
              edeg (deg) = et(il(degbands(nd,deg)), ik)
           end do
           call hpsort(degeneracy (nd), edeg, index)
           do deg = 1, degeneracy (nd)
              il(degbands(nd,deg)) = index(deg)
           end do
        end do
     end if
     !
     !   Now the order of eigenfunctions has been established
     !   for this k-point -- prepare data for next k point
     !
     do ibnd = 1, nbnd
        do ig = 1, npw
           psiold (ig, ibnd) = evc (ig, il (ibnd) )
        enddo
        do ikb = 1, nkb
           becpold (ikb, ibnd) = becp (ikb, il (ibnd) )
        enddo
     enddo
     do ig = 1, npw
        igkold (ig) = igk (ig)
     enddo
     npwold = npw
     !
     !  find degenerate eigenvalues
     !
     deg  = 0
     ndeg = 0
     do ibnd = 2, nbnd
        if ( abs (et(ibnd, ik) - et(ibnd-1, ik)) < eps ) then
           if ( deg == 0 ) then
              ndeg = ndeg + 1
              edeg (ndeg) = et(ibnd, ik)
           end if
           deg = 1
        else
           deg = 0
        end if
     end do
     !
     !  locate band crossings at degenerate eigenvalues
     !
     do nd = 1, ndeg
        deg = 0
        do ibnd = 1, nbnd
           if ( abs (et(il(ibnd), ik) - edeg (nd)) < eps ) then
              deg = deg + 1
              if (deg > maxdeg) call errore ('punch_band', &
                   ' increase maxdeg', deg)
              degbands(nd,deg) = ibnd
           end if
        end do
        degeneracy (nd) = deg
     end do
     !
#ifdef __PARA
     if (me == 1) then
#endif
     if (ik == 1) then
        write (iunpun, '(" &plot nbnd=",i4,", nks=",i4," /")') &
             nbnd, nks
     end if
     write (iunpun, '(14x,3f7.4)') xk(1,ik),xk(2,ik),xk(3,ik)
     write (iunpun, '(10f8.3)') (et (il (ibnd) , ik) &
             * rytoev, ibnd = 1, nbnd)
#ifdef __PARA
     endif
#endif
  enddo

  deallocate (index, degbands)
  deallocate (edeg, degeneracy)
  deallocate (il, ok)
  deallocate (igkold)
  deallocate (becpold, becp)
  deallocate (new, old)
  deallocate (psiold)

#ifdef __PARA
  if (me == 1) then
#endif
  close (iunpun)
#ifdef __PARA
  endif
#endif
  return
end subroutine punch_band
