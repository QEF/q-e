MODULE tetra_ip
!this module provides rouitnes for the tetrahedra method

  USE kinds, ONLY : dp

contains
  !---------------------------------------------------------------------------
  subroutine tetrahedra1 (nk1, nk2, nk3, ntetra, tetra)
    !-------------------------------------------------------------------------
    !
    ! Tetrahedron method according to P. E. Bloechl et al, PRB49, 16223 (1994)
    !
    !-------------------------------------------------------------------------
    implicit none
    
    integer, intent(in):: nk1, nk2, nk3, ntetra
    integer, intent(out) :: tetra(4,ntetra) 
    integer :: nkr, i,j,k, n, ip1,jp1,kp1, &
         n1,n2,n3,n4,n5,n6,n7,n8
    integer, parameter :: stderr = 0
    
    ! Total number of k-points
    nkr = nk1*nk2*nk3

    !  Construct tetrahedra

    do i=1,nk1
       do j=1,nk2
          do k=1,nk3
             !  n1-n8 are the indices of k-point 1-8 forming a cube
             ip1 = mod(i,nk1)+1
             jp1 = mod(j,nk2)+1
             kp1 = mod(k,nk3)+1
             n1 = (  k-1) + (  j-1)*nk3 + (  i-1)*nk2*nk3 + 1
             n2 = (  k-1) + (  j-1)*nk3 + (ip1-1)*nk2*nk3 + 1
             n3 = (  k-1) + (jp1-1)*nk3 + (  i-1)*nk2*nk3 + 1
             n4 = (  k-1) + (jp1-1)*nk3 + (ip1-1)*nk2*nk3 + 1
             n5 = (kp1-1) + (  j-1)*nk3 + (  i-1)*nk2*nk3 + 1
             n6 = (kp1-1) + (  j-1)*nk3 + (ip1-1)*nk2*nk3 + 1
             n7 = (kp1-1) + (jp1-1)*nk3 + (  i-1)*nk2*nk3 + 1
             n8 = (kp1-1) + (jp1-1)*nk3 + (ip1-1)*nk2*nk3 + 1
             !  there are 6 tetrahedra per cube (and nk1*nk2*nk3 cubes)
             n  = 6 * ( (k-1) + (j-1)*nk3 + (i-1)*nk3*nk2 )
             
             tetra (1,n+1) = n1
             tetra (2,n+1) = n2
             tetra (3,n+1) = n3
             tetra (4,n+1) = n6
             
             tetra (1,n+2) = n2
             tetra (2,n+2) = n3
             tetra (3,n+2) = n4
             tetra (4,n+2) = n6
             
             tetra (1,n+3) = n1
             tetra (2,n+3) = n3
             tetra (3,n+3) = n5
             tetra (4,n+3) = n6
             
             tetra (1,n+4) = n3
             tetra (2,n+4) = n4
             tetra (3,n+4) = n6
             tetra (4,n+4) = n8
             
             tetra (1,n+5) = n3
             tetra (2,n+5) = n6
             tetra (3,n+5) = n7
             tetra (4,n+5) = n8
             
             tetra (1,n+6) = n3
             tetra (2,n+6) = n5
             tetra (3,n+6) = n6
             tetra (4,n+6) = n7
          enddo
       enddo
    enddo
    
    !  check
    do n=1, ntetra
       do i=1,4        
          if ( tetra(i,n)<1 .or. tetra(i,n)> (nk1*nk2*nk3) ) then
             write(stderr,*) 'Something wrong with the construction of tetrahedra' 
             call exit(1)
          endif
       enddo
    enddo
  end subroutine tetrahedra1
  
  !-----------------------------------------------------------------------
  subroutine hpsort1 (n, ra, ind)  
    !---------------------------------------------------------------------
    ! sort an array ra(1:n) into ascending order using heapsort algorithm.
    ! n is input, ra is replaced on output by its sorted rearrangement.
    ! create an index table (ind) by making an exchange in the index array
    ! whenever an exchange is made on the sorted data array (ra).
    ! in case of equal values in the data array (ra) the values in the
    ! index array (ind) are used to order the entries.
    ! if on input ind(1)  = 0 then indices are initialized in the routine,
    ! if on input ind(1) != 0 then indices are assumed to have been
    !                initialized before entering the routine and these
    !                indices are carried around during the sorting process
    !
    ! no work space needed !
    ! free us from machine-dependent sorting-routines !
    !
    ! adapted from Numerical Recipes pg. 329 (new edition)
    !
    !---------------------------------------------------------------------
    implicit none  
    !-input/output variables
    integer :: n  
    integer :: ind (*)  
    real(kind=dp) :: ra (*)  
    !-local variables
    integer :: i, ir, j, l, iind  
    real(kind=dp) :: rra  
    ! initialize index array
    if (ind (1) .eq.0) then  
       do i = 1, n  
          ind (i) = i  
       enddo
    endif
    ! nothing to order
    if (n.lt.2) return  
    ! initialize indices for hiring and retirement-promotion phase
    l = n / 2 + 1  
    ir = n  
10  continue  
    ! still in hiring phase
    if (l.gt.1) then  
       l = l - 1  
       rra = ra (l)  
       iind = ind (l)  
       ! in retirement-promotion phase.
    else  
       ! clear a space at the end of the array
       rra = ra (ir)  
       !
       iind = ind (ir)  
       ! retire the top of the heap into it
       ra (ir) = ra (1)  
       !
       ind (ir) = ind (1)  
       ! decrease the size of the corporation
       ir = ir - 1  
       ! done with the last promotion
       if (ir.eq.1) then  
          ! the least competent worker at all !
          ra (1) = rra  
          !
          ind (1) = iind  
          return  
       endif
    endif
    ! wheter in hiring or promotion phase, we
    i = l  
    ! set up to place rra in its proper level
    j = l + l  
    !
    do while (j.le.ir)  
       if (j.lt.ir) then  
          ! compare to better underling
          if (ra (j) .lt.ra (j + 1) ) then  
             j = j + 1  
          elseif (ra (j) .eq.ra (j + 1) ) then  
             if (ind (j) .lt.ind (j + 1) ) j = j + 1  
          endif
       endif
       ! demote rra
       if (rra.lt.ra (j) ) then  
          ra (i) = ra (j)  
          ind (i) = ind (j)  
          i = j  
          j = j + j  
       elseif (rra.eq.ra (j) ) then  
          ! demote rra
          if (iind.lt.ind (j) ) then  
             ra (i) = ra (j)  
             ind (i) = ind (j)  
             i = j  
             j = j + j  
          else  
             ! set j to terminate do-while loop
             j = ir + 1  
          endif
          ! this is the right place for rra
       else  
          ! set j to terminate do-while loop
          j = ir + 1  
       endif
    enddo
    ra (i) = rra  
    ind (i) = iind  
    goto 10  
    !
  end subroutine hpsort1
  
  !--------------------------------------------------------------------
  subroutine weights_delta1 (energy, band,  ntetra, tetra, nkpoints,  weights)
    !--------------------------------------------------------------------
    !
    ! ... calculates weights with the tetrahedron method (P.E.Bloechl)
    !     assuming the integrand is convoluted with a delta function
    !     centered about the band 
    !
    !--------------------------------------------------------------------
    
    implicit none
    ! I/O variables
    integer, intent(in) :: nkpoints, ntetra, tetra(4, ntetra)
    real(kind=dp), intent(in) :: energy, band(nkpoints)
    real(kind=dp), intent(out) :: weights(nkpoints)
    ! local variables
    real(kind=dp) :: e1, e2, e3, e4, fact, etetra (4)
    integer :: ik, ibnd, nt, nk, ns, i, kp1, kp2, kp3, kp4, itetra (4)
    
    ! Initialize to zero the weights
    weights = 0.0d0
    
    do nt = 1, ntetra      
       !
       ! etetra are the energies at the vertexes of the nt-th tetrahedron
       !
       do i = 1, 4
          etetra (i) = band (tetra (i, nt))
       enddo
       itetra (1) = 0
       call hpsort1 (4, etetra, itetra)
       !
       ! ...sort in ascending order: e1 < e2 < e3 < e4
       !
       e1 = etetra (1)
       e2 = etetra (2)
       e3 = etetra (3)
       e4 = etetra (4)
       !
       ! kp1-kp4 are the irreducible k-points corresponding to e1-e4
       !
       kp1 = tetra (itetra (1), nt) 
       kp2 = tetra (itetra (2), nt) 
       kp3 = tetra (itetra (3), nt) 
       kp4 = tetra (itetra (4), nt) 
       !
       ! calculate weights 
       !
       ! Remember that if energy.lt.e1 or energy.ge.e4 
       ! there is no contribution in this case
       if (energy.lt.e4.and.energy.ge.e3) then                                
          fact = 1.0d0/ntetra *  (e4 - energy)**2 / (e4 - e1) / (e4 - e2) &
               / (e4 - e3)
          weights (kp1) = weights (kp1) + fact * (e4 - energy) / (e4 - e1)
          weights (kp2) = weights (kp2) + fact * (e4 - energy) / (e4 - e2)
          weights (kp3) = weights (kp3) + fact * (e4 - energy) / (e4 - e3)
          weights (kp4) = weights (kp4) + fact * ( 3.0d0 - (e4 - energy) * &
               (1.d0 / (e4 - e1) + 1.d0 / (e4 - e2) + 1.d0 / (e4 - e3) ) )
       elseif (energy.lt.e3.and.energy.ge.e2) then
          weights (kp1) = weights (kp1) + 1.0d0/ntetra * ( &
               (energy - e2)**3 * (1.0d0/ (e3 - e1)**2 / (e3 - e2) / (e4 - e2) + &
               1.0d0/ (e3 - e1) / (e4 - e1)**2 / (e4 - e2) + &
               1.0d0/ (e3 - e1)**2 / (e4 - e1) / (e4 - e2))  &
               -3.0d0 * (energy - e2)**2 * (1.0d0/ (e3 - e1)**2 / (e4 - e1) + &
               1.0d0/ (e3 - e1) / (e4 - e1)**2)  &
               -3.0d0 * (energy - e2) * ((e2 - e1)**2 - (e3 - e2)*(e4 - e2))  &
               / (e3 - e1)**2 / (e4 - e1)**2         &                                            
               +(e2 - e1) * ( (e3 - e2)/(e3 - e1) + (e4 - e2)/(e4 - e1)) & 
               / (e3 - e1) / (e4 - e1))
          weights (kp2) = weights (kp2) + 1.0d0/ntetra * ( &
               (energy - e2)**3 * (1.0d0/ (e3 - e1) / (e4 - e1) / (e4 - e2)**2 + &
               1.0d0/ (e3 - e1) / (e3 - e2) / (e4 - e2)**2 + &
               1.0d0/ (e3 - e1) / (e3 - e2)**2 / (e4 - e2))  &
               -3.0d0 * (energy - e2)**2 * (1.0d0/ (e3 - e1) / (e4 - e1) / (e4 - e2) + &
               1.0d0/ (e3 - e1) / (e3 - e2) / (e4 - e2))  &
               +3.0d0 * (energy - e2) * (1.0d0/ (e3 - e1) / (e4 - e1))  &
               +(e2 - e1) / (e3 - e1) / (e4 - e1))
          weights (kp3) = weights (kp3) + 1.0d0/ntetra * ( &
               -(energy - e2)**3 * (1.0d0/ (e3 - e1)**2 / (e4 - e1) / (e4 - e2) + &
               1.0d0/ (e3 - e1) / (e3 - e2)**2 / (e4 - e2) + &
               1.0d0/ (e3 - e1)**2 / (e3 - e2) / (e4 - e2))  &
               +3.0d0 * (energy - e2)**2 * (1.0d0/ (e3 - e1)**2 / (e4 - e1)) &
               +3.0d0 * (energy - e2) * ((e2 - e1)/ (e3 - e1)**2 / (e4 - e1))  &
               +(e2 - e1)**2 / (e3 - e1)**2 / (e4 - e1))
          weights (kp4) = weights (kp4) + 1.0d0/ntetra * ( &
               -(energy - e2)**3 * (1.0d0/ (e3 - e1) / (e3 - e2) / (e4 - e1)**2 + &
               1.0d0/ (e3 - e2) / (e4 - e1) / (e4 - e2)**2 + &
               1.0d0/ (e3 - e2) / (e4 - e1)**2 / (e4 - e2))  &
               +3.0d0 * (energy - e2)**2 * (1.0d0/ (e3 - e1) / (e4 - e1)**2) &
               +3.0d0 * (energy - e2) * ((e2 - e1)/ (e3 - e1) / (e4 - e1)**2)  &
               +(e2 - e1)**2 / (e3 - e1) / (e4 - e1)**2)
       elseif (energy.lt.e2.and.energy.ge.e1) then
          fact = 1.0d0 / ntetra * (energy - e1) **2 / (e2 - e1) / (e3 - e1) &
               / (e4 - e1)
          weights (kp1) = weights (kp1) + fact * (3.d0 - (energy - e1) * &
               (1.d0 / (e2 - e1) + 1.d0 / (e3 - e1) + 1.d0 / (e4 - e1) ) )
          weights (kp2) = weights (kp2) + fact * (energy - e1) / (e2 - e1) 
          weights (kp3) = weights (kp3) + fact * (energy - e1) / (e3 - e1) 
          weights (kp4) = weights (kp4) + fact * (energy - e1) / (e4 - e1)
       endif
    enddo
  end subroutine weights_delta1

END MODULE tetra_ip

