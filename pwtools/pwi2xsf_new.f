! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .


c     ------------------------------------------------------------------------
      program pwi2xsf_new
c     reads preprocessed (with pwi2xsf.sh) PWscf v1.2 input-file
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
     $     atomic_posunit       ! length-unit of atomic positions

      real*8
     $     celldm(6),           ! cell parameters
     $     omega,               ! cell volume (not used)
     $     alat                 ! lattice parameter

      character
     $     line*120             ! line of input
      character*3
     $     atm(maxatom)       ! atomic symbols

      integer
     $     ityp,                ! type of PP
     $     ounit,               ! output unit
     $     i, j,                ! dummies
     $     ipol,inat,           ! counters
     $     i_trimleft_white_space, ! string whitespace-triming function
     $     len                  ! length of string

      real*8
     $     tau(3,maxatom),
     +     p( 3,3 ),            ! lattice vectors (PRIMITIVE)
     +     c( 3,3 )             ! lattice vectors (CONVENTIONAL)

      logical
     $     ltaucry

      namelist/system/ ibrav, nat, celldm
      
      ounit=6

c     set default values
      nat       = 0
      ibrav     = 0
      celldm(1) = 0.0d0

c     read namelist system
      read (5,system)
      if ( nat.eq.0 .or. celldm(1).eq.0.0d0 ) then
         print *,'ERROR reading INPUT.    STOPING !!!'
         STOP
      endif
c      print *,'DEBUG: namelist read'
c      print *,'DEBUG: nat = ', nat

c     read the rest of the input
 990  continue
      read(5,'(a120)',end=999) line
c      print *,'DEBUG: read line     :', line(1:40)
      len = i_trimleft_white_space(line)
c      print *,'DEBUG: read line_trim:', line(1:40)
c     
c     CELL_PARAMETERS
c     
      if (     line(1:15) .eq. 'CELL_PARAMETERS' ) then
         read (5,*) ((p(i,j),i=1,3),j=1,3)
         do j=1,3
            do i=1,3
               c(i,j) = p(i,j)
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

c     read atoms
         do inat=1,nat
 10         continue
            read (5,'(a120)') line
c            print *,'DEBUG: line_len: line:',line(1:40)
            len = i_trimleft_white_space(line)
c            print *,'DEBUG: line_len: len: ',len, line

            if (len.eq.0) then
c     an empty line, read again
               goto 10
            endif
            atm(inat) = line(1:3)
c            print *,'DEBUG: atm():',atm(inat)
            line      = line(3:len)
            read (line,*) (tau(ipol,inat),ipol=1,3)
         enddo
      endif
      goto 990
      
 999  continue
      
      if ( ibrav.ne.0 ) then
         call latgen( ibrav, celldm,
     $        p(1,1), p(1,2), p(1,3), c(1,1), c(1,2), c(1,3), omega )
         do j=1,3
            do i=1,3
               p(i,j) = p(i,j)/celldm(1)
            end do
         end do
         call latgen_conventional(ibrav, celldm, p(1,1), p(1,2), p(1,3),
     $        c(1,1), c(1,2), c(1,3))
      endif

      alat = bohr*celldm(1)
      call write_XSF_header (alat, p, c, nat, ounit)

      do inat=1,nat
         if (     atomic_posunit .eq. BOHR_UNIT ) then            
            tau(1,inat) = bohr * tau(1,inat)
            tau(2,inat) = bohr * tau(2,inat)
            tau(3,inat) = bohr * tau(3,inat)

         elseif ( atomic_posunit .eq. ALAT_UNIT ) then
            tau(1,inat) = alat * tau(1,inat)
            tau(2,inat) = alat * tau(2,inat)
            tau(3,inat) = alat * tau(3,inat)

         elseif ( atomic_posunit .eq. CRYSTAL_UNIT ) then
            call cryst_to_cart(1, tau(1,inat), p, 1)

            tau(1,inat) = alat * tau(1,inat)
            tau(2,inat) = alat * tau(2,inat)
            tau(3,inat) = alat * tau(3,inat)
         endif

         write(ounit,'(a3,2x,3f15.10)')
     $        atm(inat), tau(1,inat), tau(2,inat), tau(3,inat)
      enddo
      END
