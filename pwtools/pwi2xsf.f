! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Tone: File adapted from pwi2xsf.f file of XCRYSDEN distribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     ------------------------------------------------------------------------
      program pwi2xsf
c     Reads pre-procesed (with pwi2xsf.sh) PW-input file
c     and converts to XSF format file
c
c     This program reads the NEWLY formated preprocessed-PW.X input
c     
c     Usage: pwi2xsf.sh < PW-preprocessed file
c     ------------------------------------------------------------------------

      implicit none

c maxtyp : maximum number of types of atoms
c maxatom: maximum number of atoms

      integer
     $     maxtyp,
     $     maxatom,
     $     ALAT_UNIT,
     $     BOHR_UNIT,
     $     ANGSTROM_UNIT,
     $     CRYSTAL_UNIT

      real*8
     $     bohr

      parameter (
     $     maxtyp  = 100,
     $     maxatom = 10000,
     $     bohr    = 0.529177d0,
     $
     $     ALAT_UNIT     = 1,
     $     BOHR_UNIT     = 2,
     $     ANGSTROM_UNIT = 3,
     $     CRYSTAL_UNIT  = 4 )          

      integer
     $     ibrav,               ! label for Bravais lattice
     $     nat,                 ! number of atoms
     $     ntyp,                ! number of pseudopotentials
     $     num_of_images,       ! number of NEB images
     $     atomic_posunit       ! length-unit of atomic positions

      real*8
     $     celldm(6),           ! cell parameters
     $     omega,               ! cell volume (not used)
     $     alat,                ! lattice parameter
     $     a, b, c, cosab, cosac, cosbc ! lattice parameters

      character
     $     calculation*80,      ! type of calculation
     $     line*120             ! line of input
      character*3
     $     atm(maxatom)       ! atomic symbols

      integer
     $     ityp,                ! type of PP
     $     ounit,               ! output unit
     $     i, j,                ! dummies
     $     inat, iim, m,        ! counters
     $     i_trimleft_white_space, ! string whitespace-triming function
     $     len                  ! length of string

      real*8
     $     x,y,z,w1,w2,         ! Cartesian coordinates & weights
     $     tau(3,maxatom),      ! atomic coordinates
     $     tau2(3,maxatom),     ! atomic coordinates (2nd image)
     +     pv( 3,3 ),           ! lattice vectors (PRIMITIVE)
     +     cv( 3,3 )            ! lattice vectors (CONVENTIONAL)

      logical
     $     ltaucry

      namelist/system/
     $     ibrav, nat, celldm, a, b, c, cosab, cosac, cosbc,
     $     calculation, num_of_images
      
      ounit=6

c     set default values
      calculation   = 'scf'
      num_of_images = 1
      nat       = 0
      ibrav     = 0
      celldm(1) = 0.0d0
      a         = 0.0D0 
      b         = 0.0D0
      c         = 0.0D0
      cosab     = 0.0D0
      cosac     = 0.0D0
      cosbc     = 0.0D0
c
c     read namelist system
c
      read (5,system)
      if ( nat.eq.0 .or. celldm(1).eq.0.0d0 ) then
         print *,'ERROR: while reading INPUT !!!'
         STOP
      endif

c     was lattice specified in terms of A,B,C,...
      if ( celldm(1) .eq. 0.0D0 .AND. a .ne. 0.0D0 ) THEN
         if ( ibrav .eq. 0 ) ibrav = 14 
         celldm(1) = a / bohr
         celldm(2) = b / a
         celldm(3) = c / a
         celldm(4) = cosab
         celldm(5) = cosac
         celldm(6) = cosbc 
      else if ( celldm(1) .ne. 0.0D0 .AND. a .ne. 0.0D0 ) THEN
         print *, 'ERROR: do not specify both celldm and a,b,c !!!'
      endif
c
c     read the rest of the input
c
 990  continue
      read(5,'(a120)',end=999) line
      len = i_trimleft_white_space(line)

c     
c     CELL_PARAMETERS
c     
      if (     line(1:15) .eq. 'CELL_PARAMETERS' ) then
         read (5,*) ((pv(i,j),i=1,3),j=1,3)
         do j=1,3
            do i=1,3
               cv(i,j) = pv(i,j)
            end do
         end do
c
c     ATOMIC_POSITIONS
c
      elseif ( line(1:16) .eq. 'ATOMIC_POSITIONS' ) then
c     find out the length-unit
         line = line(17:len)
         len  = i_trimleft_white_space(line)
         atomic_posunit = ALAT_UNIT         
         if (len.gt.0 ) then
            if ( line(1:4) .eq. 'ALAT' ) then
               atomic_posunit = ALAT_UNIT
            elseif ( line(1:4) .eq. 'BOHR' ) then
               atomic_posunit = BOHR_UNIT
            elseif ( line(1:7) .eq. 'CRYSTAL' ) then
               atomic_posunit = CRYSTAL_UNIT
            elseif ( line(1:8) .eq.'ANGSTROM') then
               atomic_posunit = ANGSTROM_UNIT
            endif
         endif

c     
c     read atoms
c     
         if ( calculation(1:3) .ne. 'NEB' )  then
            call read_atoms(nat,atm,tau)
         else
c
c     NEB: read atoms
c
            if (num_of_images.lt.2) num_of_images=2
            read (5,'(a120)') line ! line: first_image
            call read_atoms(nat,atm,tau)
            read (5,'(a120)') line ! line: second_image
            call read_atoms(nat,atm,tau2)            
         endif
      endif
      goto 990
      
 999  continue
      
      if ( ibrav.ne.0 ) then
         call latgen( ibrav, celldm,
     $        pv(1,1), pv(1,2), pv(1,3),
     $        cv(1,1), cv(1,2), cv(1,3), omega )
         do j=1,3
            do i=1,3
               pv(i,j) = pv(i,j)/celldm(1)
            end do
         end do
         call latgen_conventional(ibrav, celldm,
     $        pv(1,1), pv(1,2), pv(1,3),
     $        cv(1,1), cv(1,2), cv(1,3))
      endif

      alat = bohr*celldm(1)
      call write_XSF_header (num_of_images,alat, pv, cv, nat, ounit)
      
      do inat=1,nat
         if     ( atomic_posunit .eq. BOHR_UNIT ) then            
            tau(1,inat) = bohr * tau(1,inat)
            tau(2,inat) = bohr * tau(2,inat)
            tau(3,inat) = bohr * tau(3,inat)
            tau2(1,inat) = bohr * tau2(1,inat)
            tau2(2,inat) = bohr * tau2(2,inat)
            tau2(3,inat) = bohr * tau2(3,inat)
            
         elseif ( atomic_posunit .eq. ALAT_UNIT ) then
            tau(1,inat) = alat * tau(1,inat)
            tau(2,inat) = alat * tau(2,inat)
            tau(3,inat) = alat * tau(3,inat)
            tau2(1,inat) = alat * tau2(1,inat)
            tau2(2,inat) = alat * tau2(2,inat)
            tau2(3,inat) = alat * tau2(3,inat)
            
         elseif ( atomic_posunit .eq. CRYSTAL_UNIT ) then
            call cryst_to_cart(1, tau(1,inat), pv, 1)
            call cryst_to_cart(1, tau2(1,inat),pv, 1)
            
            tau(1,inat) = alat * tau(1,inat)
            tau(2,inat) = alat * tau(2,inat)
            tau(3,inat) = alat * tau(3,inat)
            tau2(1,inat) = alat * tau2(1,inat)
            tau2(2,inat) = alat * tau2(2,inat)
            tau2(3,inat) = alat * tau2(3,inat)
         endif
         if ( num_of_images .lt. 2 ) then
            write(ounit,'(a3,2x,3f15.10)')
     $           atm(inat), tau(1,inat), tau(2,inat), tau(3,inat)
         endif
      enddo

      if ( num_of_images .ge. 2 ) then
         m = num_of_images - 1         
         do iim=1,num_of_images            
            w1 = dble(m-(iim-1))/dble(m)
            w2 = dble(iim-1)/dble(m)
            write(ounit,'('' PRIMCOORD '',i5)') iim
            write(ounit,*) nat, 1
            do inat=1,nat
               x = w1*tau(1,inat) + w2*tau2(1,inat)
               y = w1*tau(2,inat) + w2*tau2(2,inat)
               z = w1*tau(3,inat) + w2*tau2(3,inat)
               write(ounit,'(a3,2x,3f15.10)') atm(inat), x, y, z
            enddo
         enddo
      endif
      END



c---------------------------------------------------------------------
      subroutine latgen_conventional
     +     ( ibrav, celldm, p1, p2, p3, c1, c2, c3 )
c     Generate convetional lattice
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
      subroutine read_atoms(nat,atm,coor)
c     read atomic coordinates
c     ------------------------------------------------------------------------
      implicit none
      integer
     $     nat,                 ! number of atoms
     $     ipol,inat,len,       ! counters
     $     i_trimleft_white_space ! integer-function
      character
     $     line*120             ! line of input
      character*3
     $     atm(*)               ! atomic symbols
      real*8
     $     coor(3,*)
      
      do inat=1,nat
 10      continue
         read (5,'(a120)') line
         len = i_trimleft_white_space(line)
         
         if (len.eq.0) then
c     an empty line, read again
            goto 10
         endif
         atm(inat) = line(1:3)
         line      = line(3:len)
         read (line,*) (coor(ipol,inat),ipol=1,3)
      enddo
      return
      end


c     ------------------------------------------------------------------------
      subroutine write_XSF_header (num_of_images,alat, p, c, nat, ounit)
c     writes the header for XSF structure file
c     ------------------------------------------------------------------------
      real*8
     $     alat,                ! lattice parameter
     $     p(3,3), c(3,3),      ! lattive vectors (PRIMITIVE & CONVETIONAL)
     $     p1(3,3), c1(3,3)     ! lattive vectors in ANGSTROMS unit
      integer
     $     nat,                 ! number of atoms
     $     num_of_images,       ! number of NEB images
     $     ounit                ! output unit     
      integer
     $     i, j                 ! dummies

      do i=1,3
         do j=1,3
            p1(i,j) = alat*p(i,j)
            c1(i,j) = alat*c(i,j)
         enddo
      enddo

      if (num_of_images .gt. 1)
     $     write(ounit,'('' ANIMSTEPS '',i5)') num_of_images

      write(ounit,'('' CRYSTAL'')')
      write(ounit,'(/,'' PRIMVEC'')')
      write(ounit,'(3(f15.10,2x,f15.10,2x,f15.10,/))')
     $     ((p1(i,j),i=1,3),j=1,3)
      write(ounit,'('' CONVVEC'')')
      write(ounit,'(3(f15.10,2x,f15.10,2x,f15.10,/))')
     $     ((c1(i,j),i=1,3),j=1,3)
      if (num_of_images .eq. 1) then
         write(ounit,'('' PRIMCOORD'')')
         write(ounit,*) nat, 1
      endif
      return
      end

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
