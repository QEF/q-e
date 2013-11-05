!
! Copyright (C) 2003-2012 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine summary_tran_tot()
!
! It writes transmission coefficients onto the file tran_file
!
  USE kinds, only : DP
  USE cond, ONLY : nenergy, earr, start_e, last_e, tran_tot, tran_file
  implicit none
  integer ::  i

!
! Output of T onto the file
!
  IF (tran_file.NE.' ') then 
    open (4,file=trim(tran_file),form='formatted', status='unknown')
    write(4,'("# E-Ef, T")')
    do i = start_e, last_e
      write(4,'(F12.5,3X,E14.5)') earr(i), tran_tot(i)
    enddo
    close(unit=4)
  ENDIF
  do i = start_e, last_e
    write(6,'(a8,F12.5,3X,E14.5)') 'T_tot', earr(i), tran_tot(i)
  enddo

  return

end subroutine summary_tran_tot

subroutine summary_tran_k(ien, nk1, nk2, k1, k2)
!
! Writes the k-resolved transmission 
!
! 
  USE kinds, only : DP
  USE cell_base,  ONLY : bg
  USE symm_base, ONLY: nsym, s, t_rev, time_reversal
  USE cond, ONLY : xyk, nkpts, tran_file, tran_k, tk_plot  
  implicit none

  integer ::  ien, nk1, nk2, k1, k2, c_tab, nk_full, i, j, l, k  
  integer ::  isym, ik
  integer, allocatable :: fbz_ibz(:)
  logical :: f
  real(DP) :: ktmp(2), xmin, xmax, ymin, ymax, segno  
  real(DP), parameter :: eps = 1.0e-5
  real(DP), allocatable :: xyk_full(:,:), xyk_full_cart(:,:)  
  CHARACTER(LEN=50) :: filename

  IF(tran_file.eq.' ') return

!--
! Output of T(k) into a file for the IBZ
!
  c_tab = LEN("ibz_cryst_"//trim(tran_file)) + 1
  filename =  "ibz_cryst_"//trim(tran_file)
  WRITE (filename(c_tab:c_tab+1),'(a2)') '_e'
  c_tab = c_tab + 2

  IF (ien>99) THEN
    write(filename( c_tab : c_tab+2 ),'(i3)') ien
  ELSEIF (ien>9) THEN
    write(filename( c_tab : c_tab+1 ),'(i2)') ien
  ELSE
    write(filename( c_tab : c_tab ),'(i1)') ien
  ENDIF

  open (4,file=filename,form='formatted', status='replace')
  write(4,*) "# T(k) in IBZ, [k in cryst. coor.]"
  write(4,*) "# kx ", "     ky ", "     T"
!  write(4,*) nkpts
  do i=1, nkpts
    write(4,'(3E16.5)') xyk(:,i), tran_k(i)
  end do
  close(4)
!--

!--
! Expand in the full BZ

  nk_full = nk1*nk2
  ALLOCATE( xyk_full(2,nk_full) )
  ALLOCATE( xyk_full_cart(2,nk_full) )
  ALLOCATE( fbz_ibz(nk_full) )

  xyk_full(:,:) = 0.d0
  xyk_full_cart(:,:) = 0.d0
  fbz_ibz(:) = 0

! generate k-points in the full BZ, FBZ
  DO i = 1, nk1
    DO j = 1, nk2
      k = (i-1)*nk2 + j
      xyk_full(1,k) = dble(i-1)/nk1 + dble(k1)/2/nk1
      xyk_full(2,k) = dble(j-1)/nk2 + dble(k2)/2/nk2
      xyk_full(:,k) = xyk_full(:,k) - nint( xyk_full(:,k) ) 
    ENDDO
  ENDDO

! for each k in FBZ find equivalent k-point in IBZ

do i=1, nk_full 

! run over ik in IBZ and over symmetries, isym 
  ik = 1
  isym = 1 
  do while (fbz_ibz(i).eq.0)  

    if(ik.gt.nkpts) call errore('summary_tran_k','k in IBZ not found',1)

! shift by a small 1.d-8 to get -0.5 from both +-0.5
    if (t_rev(isym)==1) then
      segno = -1.d0
    else
      segno = 1.d0
    endif 
    ktmp(1) = segno*dot_product(s(1,1:2,isym),xyk(:,ik)) + 1.d-8
    ktmp(2) = segno*dot_product(s(2,1:2,isym),xyk(:,ik)) + 1.d-8 
    ktmp(:) = ktmp(:) - nint( ktmp(:) ) - 1.d-8 

! Check if rotated k-point matches to the point in the FBZ
    f = (abs(ktmp(1)-xyk_full(1,i))<eps.AND.abs(ktmp(2)-xyk_full(2,i))<eps)
    if (f) fbz_ibz(i) = ik 

! applying time reversal
!
    if ( fbz_ibz(i).eq.0.and.time_reversal ) then
      ktmp(1) = -segno*dot_product(s(1,1:2,isym),xyk(:,ik)) + 1.d-8
      ktmp(2) = -segno*dot_product(s(2,1:2,isym),xyk(:,ik)) + 1.d-8
      ktmp(:) = ktmp(:) - nint( ktmp(:) ) - 1.d-8

! Check again
      f = (abs(ktmp(1)-xyk_full(1,i))<eps.AND.abs(ktmp(2)-xyk_full(2,i))<eps)
      if (f) fbz_ibz(i) = ik
    endif

! going either to next symm. operation or to next k in IBZ 
    isym = isym + 1
    if (isym.gt.nsym) then
      ik = ik + 1
      isym = 1
    endif

  enddo

!  write(6,*) 'FBZ_to_IBZ(i) = ', i, fbz_ibz(i) 
  if(fbz_ibz(i).eq.0) call errore('summary_tran_k','equivalent k in IBZ not found',1)

enddo
!--

!--
! cartesian coordinates (in units of 2pi/a_0)
!
  do i=1, nk_full
    do j = 1, 2
      xyk_full_cart(j, i) = bg(j, 1)*xyk_full(1, i) + bg(j, 2) &
                            *xyk_full(2, i)
    enddo
  enddo
!--

!--
! Output of T(k) into file for the whole BZ
!

! in cryst. coordinates

  c_tab = LEN("bz_cryst_"//trim(tran_file)) + 1
  filename = "bz_cryst_"//trim(tran_file)
  WRITE (filename(c_tab:c_tab+1),'(a2)') '_e'
  c_tab = c_tab + 2

  IF (ien>99) THEN
    write(filename( c_tab : c_tab+2 ),'(i3)') ien
  ELSEIF (ien>9) THEN
    write(filename( c_tab : c_tab+1 ),'(i2)') ien
  ELSE
    write(filename( c_tab : c_tab ),'(i1)') ien
  ENDIF

  open (4,file=filename,form='formatted', status='replace')
  write(4,*) "# T(k) in full BZ, [k in cryst. coor.]"
  write(4,*) "# kx ", "     ky ", "     T"
!  write(4,*) nk_full
  do i=1, nk_full               
     write(4,'(3E16.5)') xyk_full(:,i), tran_k( fbz_ibz(i) )
  end do    
  close(4)

! in cartesian coordinates (in 2pi/a_0)

  c_tab = LEN("bz_cart_"//trim(tran_file)) + 1
  filename = "bz_cart_"//trim(tran_file)
  WRITE (filename(c_tab:c_tab+1),'(a2)') '_e'
  c_tab = c_tab + 2

  IF (ien>99) THEN
    write(filename( c_tab : c_tab+2 ),'(i3)') ien
  ELSEIF (ien>9) THEN
    write(filename( c_tab : c_tab+1 ),'(i2)') ien
  ELSE
    write(filename( c_tab : c_tab ),'(i1)') ien
  ENDIF
    
  open (4,file=filename,form='formatted', status='replace')
  write(4,*) "# T(k) in full BZ, [k in cart. coor.]"
  write(4,*) "# kx(2pi/a_0)", "    ky(2pi/a_0)", "    T"
!  write(4,*) nk_full

!-- 
! plot within (tk_plot x FBZ) for better visualization

  xmin = MINVAL( xyk_full_cart(1,1:nk_full) ) 
  xmax = MAXVAL( xyk_full_cart(1,1:nk_full) ) - xmin 
  ymin = MINVAL( xyk_full_cart(2,1:nk_full) ) 
  ymax = MAXVAL( xyk_full_cart(2,1:nk_full) ) - ymin 

  xmax = xmin + tk_plot*xmax 
  ymax = ymin + tk_plot*ymax
!--

  k = 2*tk_plot

  do i=1, nk_full

   do j = -k, k  
    do l = -k, k 

     ktmp(:) = xyk_full_cart(:,i)+bg(1:2,1)*j+bg(1:2,2)*l
     f = ktmp(1).ge.xmin.and.ktmp(1).le.xmax.and.   &
         ktmp(2).ge.ymin.and.ktmp(2).le.ymax

     if (f) write(4,'(3E16.5)') ktmp(1), ktmp(2), tran_k( fbz_ibz(i) )

    enddo
   enddo

  end do
  close(4)
!--

  deallocate( xyk_full )
  deallocate( xyk_full_cart )
  deallocate( fbz_ibz )

  return
end subroutine summary_tran_k

