  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------
  subroutine wigner_seitz2 (nk1, nk2, nk3, nq1, nq2, nq3,&
       nrr_k, nrr_q, irvec, wslen, ndegen_k, ndegen_q)
  !-----------------------------------------------------------------
  !!
  !! We have nk1*nk2*nk3 electron points and nq1*nq2*nq3 phonon points 
  !! on the same grid. Assuming nq_i <= nk_i, i=1..3 we sort the corresponding
  !! wigner-seitz points in such a way that a subset 1...nrr_q < nrr_k gives 
  !! the WS points for the phonons, while the full set 1..nrr_k gives the 
  !! WS points for the electrons
  !!
  !! the unsorted electron and phonon grids are obtained by calling
  !! wigner_seitz.f90
  !!
  !! On exit, we have the same irvec for electrons and phonons, but
  !! different ndegen
  !!
  !-----------------------------------------------------------------
  USE kinds, ONLY : DP
  !
  implicit none
  !
  INTEGER, INTENT (in) :: nk1
  !! size of the uniform k mesh
  INTEGER, INTENT (in) :: nk2
  !! size of the uniform k mesh
  INTEGER, INTENT (in) :: nk3
  !! size of the uniform k mesh
  INTEGER, INTENT (in) :: nq1
  !! size of the uniform q mesh
  INTEGER, INTENT (in) :: nq2
  !! size of the uniform q mesh
  INTEGER, INTENT (in) :: nq3
  !! size of the uniform q mesh
  INTEGER, INTENT (out) :: irvec(3,20*nk1*nk2*nk3)
  !! integer components of the ir-th Wigner-Seitz grid point in the basis of the lattice vectors
  INTEGER, INTENT (out) :: nrr_k
  !! number of Wigner-Seitz grid points for electrons
  INTEGER, INTENT (out) :: nrr_q
  !! number of Wigner-Seitz grid points for electrons
  !
  REAL(kind=DP), INTENT (out) :: wslen(2*nk1*nk2*nk3)
  !! real-space length, in units of alat
  !
  ! work variables
  INTEGER :: irvec_k (3,20*nk1*nk2*nk3), ndegen_k (20*nk1*nk2*nk3), &
             irvec_q (3,20*nk1*nk2*nk3), ndegen_q (20*nk1*nk2*nk3), &
             ind2(20*nk1*nk2*nk3), ire, ir, nind
  INTEGER, ALLOCATABLE :: ind(:)
  REAL(kind=DP) :: wslen_k (20*nk1*nk2*nk3), wslen_q (20*nk1*nk2*nk3)
  !
  !  The allocation of the sorting arrays is not very clean.  However,
  !  for the moment it works.
  !
  nind = 20*nk1*nk2*nk3
  IF (nind .lt. 125) then
      nind = 125
  ENDIF
  ALLOCATE (ind(nind))
  !
  ! initialization for ihpsort (not to be removed!)
  !
  ind = 0
  ind2 = 0
  !
  !  check the bounds
  !
  IF ( nq1.gt.nk1 .or. nq2.gt.nk2 .or. nq3.gt.nk3 ) call errore &
     ('wigner_seitz2','phonon grid should be smaller than electron grid',1)
  !
  !  now generated the un-sorted points for both electrons and phonons
  !
  CALL wigner_seitz ( nk1, nk2, nk3, irvec_k, nrr_k, ndegen_k, wslen_k)  
  CALL wigner_seitz ( nq1, nq2, nq3, irvec_q, nrr_q, ndegen_q, wslen_q)  
  !
  ! loop on phonon points and find the match in the corresponding electronic list
  !
  DO ir = 1, nrr_k
    !
    ! fake index which is useful in ihpsort
    ! (I need to split the two subsets 1...nrr_q and nrr_q+1...nrr_k: 
    ! here below ind() will be between 1 and nrr_q. Therefore, if ind()
    ! is larger than nrr_q it must belong to the second subset)
    !
    ind (ir) = nrr_q + ir
  ENDDO
  DO ire = 1, nrr_k
    ir = 1
    DO while ( ( irvec_k(1,ire).ne.irvec_q(1,ir)   .or.  &
                 irvec_k(2,ire).ne.irvec_q(2,ir)   .or.  &
                 irvec_k(3,ire).ne.irvec_q(3,ir) ) .and. &
                 ir.le.nrr_q )
      ir = ir + 1
    ENDDO 
    IF ( ir.le.nrr_q ) ind (ire) = ir - 1
  ENDDO
  !
  !  Sort the electronic points accordingly
  ! 
  CALL ihpsort ( nrr_k, ind, ind2 )
  ! 
  DO ir = 1, nrr_k
    irvec (:, ir) = irvec_k (:, ind2(ir) )
    ndegen_k (ir) = ndegen_k ( ind2(ir) )
    wslen (ir) = wslen_k ( ind2(ir) )
  ENDDO
  !
  end subroutine wigner_seitz2

