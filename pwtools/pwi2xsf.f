! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .


c     ------------------------------------------------------------------------
      program pwi2xsf
c     reads preprocessed (with pwi2xsf.sh) PWscf v1.0 input-file
c     and converts to XSF format file
c     
c     Usage: pwi2xsf.sh < PW-preprocessed file
c     ------------------------------------------------------------------------

      implicit none

      integer maxtyp
      real*8
     $     bohr
      parameter (maxtyp = 100, bohr = 0.529177d0)

      integer
     $     ibrav,               ! label for Bravais lattice
     $     nat,                 ! number of atoms
     $     ntyp                 ! number of pseudopotentials

      real*8
     $     celldm(6),           ! cell parameters
     $     omega,               ! cell volume (not used)
     $     alat                 ! lattice parameter

      character
     $     dummy*80

      integer
     $     ityp,                ! type of PP
     $     atn(maxtyp),         ! nuclear charge
     $     ounit,               ! output unit
     $     i, j                 ! dummies

      real*8
     $     tau( 3 ),
     +     p( 3,3 ),            ! lattice vectors (PRIMITIVE)
     +     c( 3,3 )             ! lattice vectors (CONVENTIONAL)

      logical
     $     ltaucry

      namelist/input/ ibrav, nat, celldm, ltaucry
      
      ltaucry = .false.

      ounit=6

c     set default values
      nat       = 0
      ibrav     = 0
      celldm(1) = 0.0d0
      do i=1,maxtyp
         atn(i) = 0
      enddo

      open(1, file='nuclei.charges', status='old')
      read (1,*) ntyp
      do i=1,ntyp
         read (1,*) j,atn(j)
      enddo
      close(1)

      read (5,input)

      if ( nat.eq.0 .or. celldm(1).eq.0.0d0 ) then
         print *,'ERROR reading INPUT.    STOPPING !!!'
         STOP
      endif

      if ( ibrav.eq.0 ) then
c     read custom lattice
         read (5,*) ((p(i,j),i=1,3),j=1,3)
         read (5,*) dummy      
      else
         call latgen( ibrav, celldm, p(1,1), p(1,2), p(1,3), omega)
         do j=1,3
            do i=1,3
               p(i,j) = p(i,j)/celldm(1)
            end do
         end do
      endif
      call latgen_conventional (ibrav, celldm, p(1,1), p(1,2), p(1,3),
     $     c(1,1), c(1,2), c(1,3))

      alat = bohr*celldm(1)
      call write_XSF_header (alat, p, c, nat, ounit)

      do i=1,nat
         read(5,*) tau(1), tau(2), tau(3), ityp
         if (ltaucry) call cryst_to_cart(1, tau, p, 1)
         write(ounit,'(i3,2x,3f15.10)') atn(ityp),
     $        alat*tau(1), alat*tau(2), alat*tau(3)
      enddo
      END
