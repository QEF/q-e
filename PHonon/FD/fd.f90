!==========================================

!* File Name : fd.f90

!* Creation Date : 15-05-2013

!* Last Modified : Wed Nov  6 08:40:35 2013

!* Created By : Marco Buongiorno Nardelli 

!==========================================
program fd
  use constants,  ONLY : pi, bohr_radius_angs, amu_au, amu_ry
  use io_files,   ONLY : prefix, tmp_dir, psfile, pseudo_dir
  use io_global,  ONLY : stdout, ionode, ionode_id
  USE mp_global,  ONLY : mp_startup
  USE environment,ONLY : environment_start
  USE mp,         ONLY : mp_bcast
  USE mp_world,   ONLY : world_comm
  USE cell_base,  ONLY : tpiba2, alat,omega, at, bg, ibrav, celldm
  USE ions_base,  ONLY : amass, nat, atm, zv, tau, ntyp => nsp, ityp
  USE kinds,      ONLY : dp 
  USE gvecw,      ONLY : ecutwfc
  USE gvect,     ONLY : ecutrho
  USE wrappers,  ONLY : f_mkdir_safe

  USE symm_base
  USE symme
  USE rap_point_group, ONLY : code_group, nclass, nelem, elem, &
       which_irr, char_mat, name_rap, name_class, gname, ir_ram
  USE rap_point_group_so, ONLY : nrap, nelem_so, elem_so, has_e, &
       which_irr_so, char_mat_so, name_rap_so, name_class_so, d_spin, &
       name_class_so1
  USE rap_point_group_is, ONLY : nsym_is, sr_is, ftau_is, d_spin_is, &
       gname_is, sname_is, code_group_is
  USE fft_base, ONLY : dfftp

  implicit none
  character(len=9) :: code = 'FD'
  integer :: ios
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  character(len=200) :: pp_file
  logical :: uspp_spsi, ascii, single_file, raw, disp_only

  INTEGER :: i, ipol, apol, na, nt
  ! counter on the celldm elements
  ! counter on polarizations
  ! counter on direct or reciprocal lattice vect
  ! counter on atoms
  ! counter on type of atoms
  integer           :: nrx1,nrx2,nrx3,nr1,nr2,nr3,nb,nax,natx,inn
  integer           :: j,k,natdp, innx
  real(kind=dp)     :: r1(3),r2(3),r3(3),rr(3,3),de
  REAL(KIND=DP),    ALLOCATABLE     :: taut(:,:) 
  REAL(KIND=DP),    ALLOCATABLE     :: atomx(:,:), ttomx(:,:), atomx_a0(:,:)
  integer,     ALLOCATABLE    :: itypx(:)
  REAL(KIND=DP),    ALLOCATABLE     :: taux(:,:)
  CHARACTER(4),     ALLOCATABLE    :: atom_name(:), atom_namex(:)
  CHARACTER(100)    :: file_out, fd_outfile, fd_outfile_dir
  CHARACTER(100)    :: cna,ci,cnx, cnb, cnc,cnn,cj,ck
  character(100)    :: disp_dir,calc,fd_prefix,fd_outdir
  character(3000)    :: control,electrons,system2,kpoints
  logical           :: do_000, verbose
  LOGICAL :: atom_in_list
  LOGICAL, EXTERNAL :: eqdisp

  INTEGER, ALLOCATABLE :: ieq(:,:),seq(:,:),neq(:), atdp(:)

  INTEGER :: nclass_ref   ! The number of classes of the point group
  INTEGER :: isym
  INTEGER, ALLOCATABLE  ::  sxy(:), sxz(:), syz(:)
  REAL (dp) :: ft1, ft2, ft3
  REAL (dp) :: d(3,3),rd(3,3),dhex(3,3), dcub(3,3)
  REAL (dp) :: accep=1.0d-5
  LOGICAL :: nodispsym, noatsym,hex
  LOGICAL, ALLOCATABLE :: move_sl(:,:)
  real(DP), PARAMETER :: sin3 = 0.866025403784438597d0, cos3 = 0.5d0, &
                        msin3 =-0.866025403784438597d0, mcos3 = -0.5d0

  data dcub/ 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
  data dhex/ 1.d0, 0.d0, 0.d0, mcos3, sin3, 0.d0, 0.d0, 0.d0,  1.d0/

  NAMELIST /inputfd/ fd_prefix,nrx1,nrx2,nrx3,de,fd_outfile,fd_outdir,fd_outfile_dir, disp_only, verbose,innx, &
                     noatsym, nodispsym,hex
  NAMELIST /verbatim/ control,electrons,system2,kpoints
  CALL mp_startup ( )
  CALL environment_start ( code )

  nrx1=1
  nrx2=1
  nrx3=1
  innx=2
  verbose=.false.
  nodispsym=.false.
  noatsym=.false.
  hex=.false.

  IF ( ionode ) THEN
    CALL input_from_file ( )
    READ(5,inputfd,IOSTAT=ios)
    IF (ios /= 0) CALL errore ('FD', 'reading inputfd namelist', ABS(ios) )
    READ(5,verbatim,IOSTAT=ios)
    IF (ios /= 0) CALL errore ('FD', 'reading verbatim namelist', ABS(ios) )
    prefix=trim(fd_prefix)
    tmp_dir = trimcheck( fd_outdir )
    ios = f_mkdir_safe( TRIM( fd_outfile_dir ) )
  endif

  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )

  !reading the xml file
  call read_xml_file

  if (ionode) then
    write(6,*) '**************************************************'
    write(6,*) '* Info from the preceding pw.x run:              *'
    write(6,*) '**************************************************'
    write(6,*) ''
    write(6,*) '    prefix=  ',trim(prefix)
    write(6,*) '    outdir=  ',trim(tmp_dir)
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
    WRITE( stdout, '(/,5x,"Atomic coordinates")')                                                  
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
     write(2,'("F")')
     close(2)

    write(6,*) '**************************************************'
    write(6,*) '*Printing symmetry information and independent atoms'
    write(6,*) '**************************************************'

    IF (verbose) THEN
     WRITE( stdout, '(36x,"s",24x,"frac. trans.")')
     nsym_is=0
     DO isym = 1, nsym
        WRITE( stdout, '(/6x,"isym = ",i2,5x,a45/)') isym, sname(isym)
        IF ( ftau(1,isym).NE.0 .OR. ftau(2,isym).NE.0 .OR. &
             ftau(3,isym).NE.0) THEN
           ft1 = at(1,1)*ftau(1,isym)/dfftp%nr1 + at(1,2)*ftau(2,isym)/dfftp%nr2 + &
                at(1,3)*ftau(3,isym)/dfftp%nr3
           ft2 = at(2,1)*ftau(1,isym)/dfftp%nr1 + at(2,2)*ftau(2,isym)/dfftp%nr2 + &
                at(2,3)*ftau(3,isym)/dfftp%nr3
           ft3 = at(3,1)*ftau(1,isym)/dfftp%nr1 + at(3,2)*ftau(2,isym)/dfftp%nr2 + &
                at(3,3)*ftau(3,isym)/dfftp%nr3
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                &        " )    f =( ",f10.7," )")') &
                isym, (s(1,ipol,isym),ipol=1,3), DBLE(ftau(1,isym))/DBLE(dfftp%nr1)
           WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )")') &
                (s(2,ipol,isym),ipol=1,3), DBLE(ftau(2,isym))/DBLE(dfftp%nr2)
           WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )"/)') &
                (s(3,ipol,isym),ipol=1,3), DBLE(ftau(3,isym))/DBLE(dfftp%nr3)
           WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7, &
                &        " )    f =( ",f10.7," )")') &
                isym, (sr(1,ipol,isym),ipol=1,3), ft1
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )")') &
                (sr(2,ipol,isym),ipol=1,3), ft2
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )"/)') &
                (sr(3,ipol,isym),ipol=1,3), ft3
        ELSE
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), " )")') &
                isym,  (s (1, ipol, isym) , ipol = 1,3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )")')  (s(2,ipol,isym), ipol=1,3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )"/)') (s(3,ipol,isym), ipol=1,3)
           WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7," )")') &
                isym,  (sr (1, ipol,isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )")')  (sr (2, ipol,isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )"/)') (sr (3, ipol,isym) , ipol = 1, 3)
        END IF
     END DO
     END IF

! find equivalent atoms and determine atoms to be displaced

    ALLOCATE (ieq(nat,nsym))
    allocate (neq(nat))
    allocate (seq(nat,nsym))
    allocate (atdp(nat))
  
    ieq=0
    seq=0
    neq=0
    atdp=0

    natdp = 1
    atdp(1)=1
    na=1
    call equiv_atoms (na,nsym,irt(:,na),ieq(na,:),seq(na,:),neq(na))
    outer: do na=2,nat
       do j=1,na-1
          do k=1,neq(j)
             if(na == ieq(j,k)) then
               cycle outer
             end if
          end do
       end do 
       call equiv_atoms (na,nsym,irt(:,na),ieq(na,:),seq(na,:),neq(na))
       natdp=natdp+1
       atdp(natdp)=na
    end do outer

! bypass for testing

    if (noatsym) then

       natdp =nat
       do j=1,nat
          atdp(j)=j
       end do

    end if

    write(*,fmt='(a,i0)') 'Number of independent atoms: ',natdp
    write(*,advance='no',fmt='(a)') 'Atoms to be displaced:'
    write(*,*) atdp(1:natdp)

    write(6,*) '**************************************************'
    write(6,*) '* Info for supercell calculation                 *'
    write(6,*) '**************************************************'
    write(6,*) ''
    WRITE(6,*) 'a1: ',(r1(j), j=1,3)
    WRITE(6,*) 'a2: ',(r2(j), j=1,3)
    WRITE(6,*) 'a3: ',(r3(j), j=1,3)

   allocate(sxy(nat))
   allocate(sxz(nat))
   allocate(syz(nat))
   allocate(move_sl(3,nat))

   ALLOCATE(itypx(natx))
   ALLOCATE(ttomx(3,natx))
   ALLOCATE(atomx_a0(3,natx))

   move_sl = .true.
   sxy = 0
   sxz = 0
   syz = 0

    ! find if the cartesian displacements are independent (only for symmorphic operations)

    nax=0 
    do nr1=1,nrx1
      do nr2=1,nrx2
        do nr3=1,nrx3
          do na=1,nat
            nax=nax+1
            do i=1,3
             atomx(i,nax)=taut(i,na)+(nr1-1)*r1(i)+(nr2-1)*r2(i)+(nr3-1)*r3(i)
             atomx_a0(i,nax)=tau(i,na)+(nr1-1)*at(i,1)+(nr2-1)*at(i,2)+(nr3-1)*at(i,3)
            enddo
            atom_namex(nax)=atom_name(na)
            itypx(nax)=ityp(na)
          enddo
        enddo
      enddo
    enddo

   ! define the supercell vectors

   at(:,1)=nrx1*at(:,1)
   at(:,2)=nrx2*at(:,2)
   at(:,3)=nrx3*at(:,3)
   CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )

   if ( .not.nodispsym ) then

   IF ( ALLOCATED( irt ) ) DEALLOCATE( irt )
   ALLOCATE( irt(48,natx))

   do k=1,natdp
      nb=atdp(k)
 ! center on the nb atom to find the correct symmetry operations
      do na=1,natx
         ttomx(:,na) = atomx_a0(:,na) - atomx_a0(:,nb)
      end do

      call find_sym_ifc( natx, ttomx, itypx )

!if (verbose) then
!   do na=1,natx
!      write(*,'(24i3)') (irt(i,na),i=1,nsym)
!   end do
!   print*, '============='
!end if

if (hex) then
   d=dhex
else
   d=dcub
end if

    ! find if the cartesian displacements are independent also in the supercell (as it should be!) and
    ! identify one symmetry operation to rotate forces

    xy_sl: do isym=1,nsym
       rd(:,1) = sr(:,1,isym) * d(1,1) + &
                 sr(:,2,isym) * d(2,1) + &
                 sr(:,3,isym) * d(3,1)
       if (eqdisp(rd(:,1),d(:,2),accep) .and. &
          ft(1,isym) == 0.d0 .and. &
          ft(2,isym) == 0.d0 .and. &
          ft(3,isym) == 0.d0) then
          if(verbose) write(6,'("on atom",i4," x and y displacements are equivalent for symmetry operation",i4)') nb,isym
          move_sl(2,nb) = .false.
          sxy(nb)=isym
          exit
       end if
    end do xy_sl

    xz_sl: do isym=1,nsym
       rd(:,1) = sr(:,1,isym) * d(1,1) + &
                 sr(:,2,isym) * d(2,1) + &
                 sr(:,3,isym) * d(3,1)
       if (eqdisp(rd(:,1),d(:,3),accep) .and. &
          ft(1,isym) == 0.d0 .and. &
          ft(2,isym) == 0.d0 .and. &
          ft(3,isym) == 0.d0) then
          if(verbose) write(6,'("on atom",i4," x and z displacements are equivalent for symmetry operation",i4)') nb,isym
          move_sl(3,nb) = .false.
          sxz(nb)=isym
          exit
       end if
    end do xz_sl

    yz_sl: do isym=1,nsym
       rd(:,2) = sr(:,1,isym) * d(1,2) + &
                 sr(:,2,isym) * d(2,2) + &
                 sr(:,3,isym) * d(3,2)
       if (eqdisp(rd(:,2),d(:,3),accep) .and. &
          ft(1,isym) == 0.d0 .and. &
          ft(2,isym) == 0.d0 .and. &
          ft(3,isym) == 0.d0) then
          if(verbose) write(6,'("on atom",i4," y and z displacements are equivalent for symmetry operation",i4)') nb,isym
          move_sl(3,nb) = .false.
          syz(nb)=isym
          exit
       end if
    end do yz_sl

  end do

end if

    do_000 = .true.
    do inn=1,innx 
       do k=1,natdp
          na=atdp(k)
          do i=1,3
             if (move_sl(i,na)) then
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

             if(.not.disp_only) then
   
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
             write(2,*) 'ecutrho= ',ecutrho
             write(2,'(300a)') trim(system2)
             write(2,*) '/'
             write(2,*) '&ELECTRONS'
             write(2,*) trim(electrons)
             write(2,*) '/'
             write(2,*) 'ATOMIC_SPECIES'
             DO nt = 1, ntyp
               WRITE(2, '(a6,3x,f10.3,6x,a)') &
                 atm(nt), amass(nt),TRIM (psfile(nt))
             ENDDO

             end if  
   
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
   
             if(.not.disp_only) then
   
             call replace_cr(kpoints)
             write(2,*) trim(kpoints) 
   
             end if
   
             write(2,*) 'CELL_PARAMETERS {angstrom}'
             WRITE(2,*) (nrx1*r1(j), j=1,3)
             WRITE(2,*) (nrx2*r2(j), j=1,3)
             WRITE(2,*) (nrx3*r3(j), j=1,3)

             CLOSE(2)
             if (do_000) then
               do_000=.false.
               go to 101
             endif
            end if
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

subroutine equiv_atoms(na,nsym,example,res,pos,k)
  implicit none
  integer :: nsym
  integer :: example(nsym),na         ! The input
  integer :: res(nsym),pos(nsym)  ! The output
  integer :: k                   ! The number of unique elements
  integer :: i, j

  k = 1
  res(1) = example(1)
  pos(1) = 1
  outer: do i=2,size(example)
     do j=1,k
        if (res(j) == example(i) .or. example(i) == na) then
           ! Found a match so start looking again
           cycle outer
        end if
     end do
     ! No match found so add it to the output
     k = k + 1
     res(k) = example(i)
     pos(k) = i
  end do outer

  write(*,advance='no',fmt='(a,i0,a,i0,a)') 'Atom ',na, ' has ',k,' equivalent(s): '
  write(*,*) res(1:k)
  write(*,advance='no',fmt='(a)') 'for symmetry operation(s):  '
  write(*,*) pos(1:k)
end subroutine equiv_atoms

!-----------------------------------------------------------------------
logical function eqdisp (x, y, accep )
  !-----------------------------------------------------------------------
  !
  !   This function test if the difference x-y is zero
  !   x, y = 3d vectors
  !
  USE kinds
  implicit none
  real(DP), intent(in) :: x (3), y (3), accep
  !
  eqdisp = abs( x(1)-y(1) ) < accep .and. &
           abs( x(2)-y(2) ) < accep .and. &
           abs( x(3)-y(3) ) < accep
  !
  return
end function eqdisp
