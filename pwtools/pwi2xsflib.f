c
c---------------------------------------------------------------------
      subroutine latgen_conventional
     +     ( ibrav, celldm, p1, p2, p3, c1, c2, c3 )
c---------------------------------------------------------------------
c
c   Conventional crystallographic vectors c1, c2, and c3.
c   See "latgen" for the meaning of variables
c
      implicit none
c
c     First the input variables
c
      real*8 
     +     celldm( 6 ),         ! input : the dimensions of the lattice
     +     p1( 3 ),             ! input : first lattice vector (PRIMITIVE)
     +     p2( 3 ),             ! input : second lattice vector
     +     p3( 3 ),             ! input : third lattice vector
     +     c1( 3 ),             ! output: first lattice vector(CONVENTIONAL)
     +     c2( 3 ),             ! output: second lattice vector
     +     c3( 3 )              ! output: third lattice vector
      integer
     +       ibrav          ! input: the index of the Bravais lattice
c
      integer i
c
c
      do i = 1, 3
         c1(i) =0.d0
         c2(i) =0.d0
         c3(i) =0.d0
      end do
c
      if ( ibrav .eq. 2 .or. ibrav .eq.3 ) then
c
c     fcc and bcc lattice
c
         c1( 1 ) = 1.0d0
         c2( 2 ) = 1.0d0
         c3( 3 ) = 1.0d0
c
      else if ( ibrav .eq. 7 ) then
c
c     body centered tetragonal lattice
c
         if ( celldm( 1 ) .le. 0.d0 .or. celldm( 3 ) .le. 0.d0 ) 
     +      call errore( 'latgen', 'wrong celldm', 7 )
         c1(1) = 1.0d0
         c2(2) = 1.0d0
         c3(3) = celldm(3)
c
      else if ( ibrav .eq. 10 ) then
c
c     All face centered orthorombic lattice
c
         if ( celldm( 1 ) .le. 0.d0 .or. celldm( 2 ) .le. 0.d0
     +        .or. celldm( 3 ) .le. 0.d0 ) 
     +      call errore( 'latgen', 'wrong celldm', 10 )
         c1(1) = 1.0d0
         c2(2) = celldm(2)
         c3(3) = celldm(3)
c
      elseif ( ibrav .eq. 11 ) then
c
c     Body centered orthorombic lattice
c
         if ( celldm( 1 ) .le. 0.d0 .or. celldm( 2 ) .le. 0.d0
     +        .or. celldm( 3 ) .le. 0.d0 ) 
     +      call errore( 'latgen', 'wrong celldm', 11 )
         c1(1) = 1.0d0
         c2(2) = celldm(2)
         c3(3) = celldm(3)
      else
c     **********
c     all other cases : just copy p vectors to c vectors !!!
c     **********
         do  i = 1, 3
            c1( i ) = p1( i )
            c2( i ) = p2( i )
            c3( i ) = p3( i )
         enddo
      end if
c
      return
      end

c     ------------------------------------------------------------------------
      subroutine write_XSF_header (alat, p, c, nat, ounit)
c     writes the header for XSF structure file
c     ------------------------------------------------------------------------
      real*8
     $     alat,                ! lattice parameter
     $     p(3,3), c(3,3),      ! lattive vectors (PRIMITIVE & CONVENTIONAL)
     $     p1(3,3), c1(3,3)     ! lattive vectors in ANGSTROMS unit
      integer
     $     nat,                 ! number of atoms
     $     ounit                ! output unit     
      integer
     $     i, j                 ! dummies

      do i=1,3
         do j=1,3
            p1(i,j) = alat*p(i,j)
            c1(i,j) = alat*c(i,j)
         enddo
      enddo

      write(ounit,'('' CRYSTAL'')')
      write(ounit,'(/,'' PRIMVEC'')')
      write(ounit,'(3(f15.10,2x,f15.10,2x,f15.10,/))')
     $     ((p1(i,j),i=1,3),j=1,3)
      write(ounit,'('' CONVVEC'')')
      write(ounit,'(3(f15.10,2x,f15.10,2x,f15.10,/))')
     $     ((c1(i,j),i=1,3),j=1,3)
      write(ounit,'('' PRIMCOORD'')')
      write(ounit,*) nat, 1
      return
      end
c
c     -------------------------------------------------
      integer function i_trimleft_white_space(word)
c     trim left white spaces out of word
c     -------------------------------------------------
      character word*(*), auxword*80
 
      ilen=len(word)
      auxword=word
      do i=1,ilen
         if ( word(i:i) .eq. ' ' ) then
            auxword=word(i+1:ilen)
         else
            goto 1
         endif
      enddo
 1    continue
      i_trimleft_white_space=len(word)
      word=auxword(1:i_trimleft_white_space)
      return
      END
