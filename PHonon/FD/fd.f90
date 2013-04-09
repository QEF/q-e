!-----------------------------------------------------------------------
!By Luis Agapito @ Marco Buongiorno Nardelli's UNT group
!Wed Apr  3 18:48:06 CDT 2013
!
!Based on the ifcdisp.f90 program of
!Marco Buongiorno Nardelli
!-----------------------------------------------------------------------
program fd
  use constants,  ONLY : pi, bohr_radius_angs, amu_au, amu_ry
  use io_files,   ONLY : prefix, tmp_dir, outdir
  use io_files,   ONLY : psfile, pseudo_dir
  use io_global,  ONLY : stdout, ionode, ionode_id
  USE mp_global,  ONLY : mp_startup,mpime,kunit
  USE environment,ONLY : environment_start
  USE mp,         ONLY : mp_bcast
  USE cell_base,  ONLY : tpiba2, alat,omega, at, bg, ibrav, celldm
  USE ions_base,  ONLY : amass, nat, atm, zv, tau, ntyp => nsp, ityp
  USE kinds,      ONLY : dp 
  USE wvfct,      ONLY : ecutwfc

  implicit none
  character(len=9) :: code = 'FD'
  integer :: ios, kunittmp
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  character(len=200) :: pp_file
  logical :: uspp_spsi, ascii, single_file, raw

  INTEGER :: i, ipol, apol, na, nt
  ! counter on the celldm elements
  ! counter on polarizations
  ! counter on direct or reciprocal lattice vect
  ! counter on atoms
  ! counter on type of atoms
  integer           :: nrx1,nrx2,nrx3,nr1,nr2,nr3,nb,nax,natx,inn
  integer           :: j
  real(kind=dp)     :: r1(3),r2(3),r3(3),rr(3,3),de
  REAL(KIND=DP),    ALLOCATABLE     :: taut(:,:) 
  REAL(KIND=DP),    ALLOCATABLE     :: atomx(:,:) 
  REAL(KIND=DP),    ALLOCATABLE     :: taux(:,:)
  CHARACTER(4),     ALLOCATABLE    :: atom_name(:), atom_namex(:)
  CHARACTER(100)    :: file_out, fd_outfile, fd_outfile_dir
  CHARACTER(100)    :: cna,ci,cnx, cnb, cnc,cnn,cj,ck
  character(100)    :: disp_dir,calc,fd_prefix,fd_outdir
  character(300)    :: control,electrons,system2,kpoints
  logical           :: do_000

  NAMELIST /inputfd/ fd_prefix,nrx1,nrx2,nrx3,de,fd_outfile,fd_outdir,fd_outfile_dir
  NAMELIST /verbatim/ control,electrons,system2,kpoints
  CALL mp_startup ( )
  CALL environment_start ( code )

  !  IF ( TRIM( outdir ) == ' ' ) outdir = './'


  IF ( ionode ) THEN
    CALL input_from_file ( )
    READ(5,inputfd,IOSTAT=ios)
    IF (ios /= 0) CALL errore ('FD', 'reading inputfd namelist', ABS(ios) )
    READ(5,verbatim,IOSTAT=ios)
    IF (ios /= 0) CALL errore ('FD', 'reading verbatim namelist', ABS(ios) )
    prefix=trim(fd_prefix)
    outdir=trim(fd_outdir)
    tmp_dir = trimcheck( outdir )
    call system('mkdir '//trim(fd_outfile_dir))
  endif

  CALL mp_bcast( outdir, ionode_id )
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )

  !reading the xml file
  call read_file

  if (ionode) then
    write(6,*) '**************************************************'
    write(6,*) '* Info from the preceding pw.x run:              *'
    write(6,*) '**************************************************'
    write(6,*) ''
    write(6,*) '    prefix=  ',trim(prefix)
    write(6,*) '    outdir=  ',trim(outdir)
    write(6,*) '    ectuwfc= ',ecutwfc, 'Ry'
    natx=nat*nrx1*nrx2*nrx3
    ALLOCATE(atom_name(nat))
    ALLOCATE(atom_namex(natx))
    ALLOCATE(taut(3,nat))
    ALLOCATE(atomx(3,natx)) 
    ALLOCATE(taux(3,natx))

    WRITE( stdout, 100) ibrav, alat, omega, nat, ntyp
    100 FORMAT(5X, &
      &     'bravais-lattice index     = ',I12,/,5X, &
      &     'lattice parameter (alat)  = ',F12.4,'  a.u.',/,5X, &
      &     'unit-cell volume          = ',F12.4,' (a.u.)^3',/,5X, &
      &     'number of atoms/cell      = ',I12,/,5X, &
      &     'number of atomic types    = ',I12)
    !
    WRITE( stdout, '(/2(3X,3(2X,"celldm(",I1,")=",F11.6),/))' ) &
      ( i, celldm(i), i = 1, 6 )
    !
    !lattice vectors in Angs
    rr = at*alat*bohr_radius_angs
    r1(:) = rr(:,1)
    r2(:) = rr(:,2)
    r3(:) = rr(:,3)
    WRITE( stdout, '(5X, &
      &     "Lattice vectors: (cart. coord. in Angs)",/, &
      &       3(15x,"a(",i1,") = (",3f11.6," )  ",/ ) )')  (apol,  &
      (rr (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
    !atoms
    do na=1,nat
      atom_name(na)=atm(ityp(na))
    enddo
    do na=1,nat
      taut(:,na) = tau(:,na)*alat*bohr_radius_angs
    enddo
    WRITE( stdout, '(/,5x,"Atomic coordiantes")')                                                  
    WRITE( stdout, '(5x,"site n.     atom                  positions (Angs)")')
    WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")') &
      (na, atom_name(na), na, (taut(ipol,na), ipol=1,3), na=1,nat)

    write(6,*) '**************************************************'
    write(6,*) '*Printing dynamical matrix header --> header.txt *'
    write(6,*) '1 amu (Ry) = ', amu_ry 
    write(6,*) '**************************************************'

     OPEN(unit=2,file=trim(fd_outfile_dir)//'/header.txt',status='unknown',form='formatted')
     WRITE(2,'(i3,i5,i3,6f11.7)') ntyp,nat,ibrav,celldm
     if (ibrav==0) then
        write (2,'(2x,3f15.9)') ((at(i,j),i=1,3),j=1,3)
     end if
     DO nt = 1,ntyp
        WRITE(2,*) nt," '",atm(nt),"' ",amass(nt)*amu_ry
     END DO
     DO na=1,nat
        WRITE(2,'(2i5,3f18.10)') na,ityp(na),(tau(j,na),j=1,3)
     END DO
     close(2)
              
    write(6,*) '**************************************************'
    write(6,*) '* Info for macrocell calculation                 *'
    write(6,*) '**************************************************'
    write(6,*) ''
    WRITE(6,*) 'a1: ',(r1(j), j=1,3)
    WRITE(6,*) 'a2: ',(r2(j), j=1,3)
    WRITE(6,*) 'a3: ',(r3(j), j=1,3)

    nax=0 
    do nr1=1,nrx1
      do nr2=1,nrx2
        do nr3=1,nrx3
          do na=1,nat
            nax=nax+1
            do i=1,3
              atomx(i,nax)=taut(i,na)+(nr1-1)*r1(i)+(nr2-1)*r2(i)+(nr3-1)*r3(i)
            enddo
            atom_namex(nax)=atom_name(na)
            !        write(2,'(a4,3(f15.9,1x))') atom_namex(nax), atomx(1,nax),atomx(2,nax),atomx(3,nax)
          enddo
        enddo
      enddo
    enddo

    do_000 = .true.
    do inn=1,2 
      do i=1,3
        do na=1,nat
          taux(:,:) = atomx(:,:)
          if (inn==1) then
            taux(i,na) = taux(i,na)+de
          else
            taux(i,na) = taux(i,na)-de
          endif
          write(cnx,*) inn
          !cnx=1/2 positive/negative displacement
          !cna=atomic index within the primitive cell
          !ci =1,2,3 x,y,z cartesian displacement direction 
          cnx=adjustl(cnx)  
          write(cna,*) na
          cna=adjustl(cna)
          write(ci,*) i
          ci=adjustl(ci)
          101   if (do_000) then
            OPEN(2,FILE=TRIM(fd_outfile_dir)//'/'//trim(fd_outfile)//'.0.0.0'//'.in',FORM='formatted')
          else
            file_out=trim(fd_outfile_dir)//'/'//TRIM(fd_outfile)//'.'//TRIM(cnx)//'.'//TRIM(ci)//'.'//TRIM(cna)//'.in'
            OPEN(2,FILE=TRIM(file_out),FORM='formatted')
          endif
          !Writing the input file
          write(2,*) '&CONTROL'
          write(2,*) 'calculation="scf"'
          write(2,*) 'tprnfor    = .true.'
          write(2,*) trim(control)
          write(2,*) '/'
          write(2,*) '&SYSTEM'
          write(2,*) 'ibrav  = 0'
          write(2,*) 'nat    = ', natx
          write(2,*) 'ntyp   = ', ntyp
          write(2,*) 'ecutwfc= ',ecutwfc
          write(2,*) trim(system2)
          write(2,*) '/'
          write(2,*) '&ELECTRONS'
          write(2,*) trim(electrons)
          write(2,*) '/'
          write(2,*) 'ATOMIC_SPECIES'
          DO nt = 1, ntyp
            WRITE(2, '(a6,3x,f10.3,6x,a)') &
              atm(nt), amass(nt),TRIM (psfile(nt))
          ENDDO
          write(2,*) 'ATOMIC_POSITIONS {angstrom}'
          if (do_000) then
            do nb=1,natx
              write(2,'(a4,3(f15.9,1x))') atom_namex(nb), atomx(1,nb),atomx(2,nb),atomx(3,nb) 
            enddo 
          else
            do nb=1,natx
              write(2,'(a4,3(f15.9,1x))') atom_namex(nb), taux(1,nb),taux(2,nb),taux(3,nb) 
            enddo 
          endif
          call replace_cr(kpoints)
          write(2,*) trim(kpoints) 
          write(2,*) 'CELL_PARAMETERS {angstrom}'
          WRITE(2,*) (nrx1*r1(j), j=1,3)
          WRITE(2,*) (nrx2*r2(j), j=1,3)
          WRITE(2,*) (nrx3*r3(j), j=1,3)
          CLOSE(2)
          if (do_000) then
            do_000=.false.
            go to 101
          endif
        enddo
      enddo
    enddo
  END IF


  call stop_pp
  stop
end program fd

subroutine replace_cr(control)
  character(300),intent(inout) :: control 
  character(300) :: control_split=''
  character      :: letter
  integer :: i

  do i=1,len(trim(control))
    letter = control(i:i)
    if (letter.eq.'*') then
      control_split(i:i)=char(10)
    else
      control_split(i:i)=control(i:i)
    endif
  enddo
  control=control_split
end subroutine replace_cr

