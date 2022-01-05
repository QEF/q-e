!
! Copyright (C) 2021-2022 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE wxml
  !
  ! Poor-man FoX_wxml replacement - Paolo Giannozzi, 2022
  !
  USE upf_kinds, ONLY : dp
  use xmltools
  !
  implicit none
  type :: xmlf_t
     integer :: unit = -1
  end type xmlf_t
  character(len=80), save :: opentag = ''
  !
  private
  public :: xmlf_t, xml_openfile, xml_close, xml_addcharacters, &
            xml_addattribute, xml_declarenamespace, xml_newelement, &
            xml_endelement, xml_addnewline, xml_addcomment
  !
  INTERFACE xml_addcharacters
     MODULE PROCEDURE xml_addcharacters_c, xml_addcharacters_l, &
                      xml_addcharacters_r, xml_addcharacters_rv,&
                      xml_addcharacters_i, xml_addcharacters_iv,&
                      xml_addcharacters_rm
  END INTERFACE xml_addcharacters
  !
  INTERFACE xml_addattribute
     MODULE PROCEDURE xml_addattribute_c, xml_addattribute_r, &
             xml_addattribute_i, xml_addattribute_l, xml_addattribute_iv
  END INTERFACE xml_addattribute
  !
CONTAINS
  !
  subroutine xml_openfile( filename, xf, unit, pretty_print, replace, &
       namespace, iostat)
    !
    character(len=*), intent(in) :: filename
    type(xmlf_t), intent(out) :: xf
    ! ignored
    integer, intent(in) :: unit
    integer, intent(out) :: iostat
    ! ignored
    logical, intent(in)  :: pretty_print, replace, namespace 
    integer :: iun
    !
    iun = xml_open_file ( filename )
    if ( iun == -1 ) then
       iostat = 1
    else
       iostat = 0
       ! dirty  trick to have the same format with no changes to qexsd.f90
       if ( replace) then
          call add_attr('version','1.0')
          call add_attr('encoding','UTF-8')
          call xmlw_writetag ( 'xml', '?' )
       end if
    end if
    xf%unit = iun
    !
  end subroutine xml_openfile
  !
  subroutine xml_close ( xf )
    type(xmlf_t), intent(inout) :: xf
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       call xmlw_closetag ( )
       call xml_closefile ( )
       xf%unit = -1 
    end if
    !
  end subroutine xml_close
  !
  subroutine xml_declarenamespace ( xf, prefix, nsURI )
    !
    type(xmlf_t), intent(in) :: xf
    character(len=*), intent(in) :: prefix, nsURI
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       call add_attr('xmlns:'//trim(prefix), nsURI)
    end if
    !
  end subroutine xml_declarenamespace
  !
  subroutine xml_addattribute_c ( xf, name, value )
    !
    type(xmlf_t), intent(in) :: xf
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       call add_attr(name, value)
    end if
    !
  end subroutine xml_addattribute_c
  !
  subroutine xml_addattribute_r ( xf, name, value )
    !
    type(xmlf_t), intent(in) :: xf
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       call add_attr(name, value)
    end if
    !
  end subroutine xml_addattribute_r
  !
  subroutine xml_addattribute_iv ( xf, name, value )
    !
    type(xmlf_t), intent(in) :: xf
    character(len=*), intent(in) :: name
    integer, intent(in) :: value(:)
    character(len=80) :: cvalue
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       write(cvalue,*) value
       call add_attr(name, cvalue)
    end if
    !
  end subroutine xml_addattribute_iv
  !
  subroutine xml_addattribute_i ( xf, name, value )
    !
    type(xmlf_t), intent(in) :: xf
    character(len=*), intent(in) :: name
    integer, intent(in) :: value
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       call add_attr(name, value)
    end if
    !
  end subroutine xml_addattribute_i
  !
  subroutine xml_addattribute_l ( xf, name, value )
    !
    type(xmlf_t), intent(in) :: xf
    character(len=*), intent(in) :: name
    logical, intent(in) :: value
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       ! FoX compatibility
       if ( value ) then
          call add_attr(name, 'true')
       else
          call add_attr(name, 'false')
       end if
    end if
    !
  end subroutine xml_addattribute_l
  !
  subroutine xml_newelement (xf, name)
    !
    type(xmlf_t), intent(in) :: xf
    character(len=*), intent(in) :: name
    !
    integer :: ierr
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       ! workaround for different logic
       if ( opentag /= '' ) call xmlw_opentag ( opentag, ierr )
       opentag = name
    end if
  end subroutine xml_newelement
  !
  subroutine xml_endelement (xf, name)
    !
    type(xmlf_t), intent(in) :: xf
    character(len=*), intent(in) :: name
    integer :: ierr
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       ! workaround for different logic
       if ( opentag /= '' ) call xmlw_opentag ( opentag, ierr )
       call xmlw_closetag ( )
       opentag = ''
    end if
  end subroutine xml_endelement
  !
  subroutine xml_addcharacters_c ( xf, field, fmt )
    !
    type(xmlf_t), intent(in) :: xf
    character(len=*), intent(in) :: field
    character(len=*), intent(in), optional :: fmt
    !
    integer :: ierr
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       ! workaround: write previous tag when contents is added
       ! (FoX adds attributes after tag, xmltools expects them before tag)
       if ( opentag /= '') then
          call xmlw_opentag ( opentag, ierr )
          if ( ierr /= 0 ) print *, 'xml_addcharacter: ierr = ', ierr
          opentag = ''
       end if
       write( xf%unit, '(A)' ) trim(field)
    end if
    !
  end subroutine xml_addcharacters_c
  !
  subroutine xml_addcharacters_l ( xf, field, fmt )
    !
    type(xmlf_t), intent(in) :: xf
    logical, intent(in) :: field
    character(len=*), intent(in), optional :: fmt
    !
    integer :: ierr
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       ! workaround: write previous tag when contents is added
       ! (FoX adds attributes after tag, xmltools expects them before tag)
       if ( opentag /= '') then
          call xmlw_opentag ( opentag, ierr )
          if ( ierr /= 0 ) print *, 'xml_addcharacter: ierr = ', ierr
          opentag = ''
       end if
       ! FoX compatibility
       if ( field ) then
          write( xf%unit, '("true")' )
       else
          write( xf%unit, '("false")' )
       end if
    end if
    !
  end subroutine xml_addcharacters_l
  !
  subroutine xml_addcharacters_r ( xf, field, fmt )
    !
    type(xmlf_t), intent(in) :: xf
    real(dp), intent(in) :: field
    character(len=*), intent(in), optional :: fmt
    !
    integer :: ierr
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       ! workaround: write previous tag when contents is added
       ! (FoX adds attributes after tag, xmltools expects them before tag)
       if ( opentag /= '') then
          call xmlw_opentag ( opentag, ierr )
          if ( ierr /= 0 ) print *, 'xml_addcharacter: ierr = ', ierr
          opentag = ''
       end if
       write( xf%unit, * ) field
    end if
    !
  end subroutine xml_addcharacters_r
  !
  subroutine xml_addcharacters_rv( xf, field, fmt )
    !
    type(xmlf_t), intent(in) :: xf
    real(dp), intent(in) :: field(:)
    character(len=*), intent(in), optional :: fmt
    !
    integer :: ierr
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       ! workaround: write previous tag when contents is added
       ! (FoX adds attributes after tag, xmltools expects them before tag)
       if ( opentag /= '') then
          call xmlw_opentag ( opentag, ierr )
          if ( ierr /= 0 ) print *, 'xml_addcharacter: ierr = ', ierr
          opentag = ''
       end if
       write( xf%unit, '(1p3es24.15)' ) field
    end if
    !
  end subroutine xml_addcharacters_rv
  !
  subroutine xml_addcharacters_rm( xf, field, fmt )
    !
    type(xmlf_t), intent(in) :: xf
    real(dp), intent(in) :: field(:,:)
    character(len=*), intent(in), optional :: fmt
    !
    integer :: ierr
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       ! workaround: write previous tag when contents is added
       ! (FoX adds attributes after tag, xmltools expects them before tag)
       if ( opentag /= '') then
          call xmlw_opentag ( opentag, ierr )
          if ( ierr /= 0 ) print *, 'xml_addcharacter: ierr = ', ierr
          opentag = ''
       end if
       write( xf%unit, '(1p3es24.15)' ) field
    end if
    !
  end subroutine xml_addcharacters_rm
  !
  subroutine xml_addcharacters_i ( xf, field, fmt )
    !
    type(xmlf_t), intent(in) :: xf
    integer, intent(in) :: field
    character(len=*), intent(in), optional :: fmt
    !
    integer :: ierr
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       ! workaround: write previous tag when contents is added
       ! (FoX adds attributes after tag, xmltools expects them before tag)
       if ( opentag /= '') then
          call xmlw_opentag ( opentag, ierr )
          if ( ierr /= 0 ) print *, 'xml_addcharacter: ierr = ', ierr
          opentag = '' 
       end if
       write( xf%unit, * ) field
    end if
    !
  end subroutine xml_addcharacters_i
  !
  subroutine xml_addcharacters_iv( xf, field, fmt )
    !
    type(xmlf_t), intent(in) :: xf
    integer, intent(in) :: field(:)
    character(len=*), intent(in), optional :: fmt
    !
    integer :: ierr
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       ! workaround: write previous tag when contents is added
       ! (FoX adds attributes after tag, xmltools expects them before tag)
       if ( opentag /= '') then
          call xmlw_opentag (opentag, ierr )
          if ( ierr /= 0 ) print *, 'xml_addcharacter: ierr = ', ierr
          opentag = ''
       end if
       write( xf%unit, '(6i12)' ) field
    end if
    !
  end subroutine xml_addcharacters_iv
  !
  subroutine xml_addnewline ( xf )
    !
    type(xmlf_t), intent(in) :: xf
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       continue
       ! write( xf%unit, * )
    end if
    !
  end subroutine xml_addnewline
  !
  subroutine xml_addcomment ( xf, comment )
    !
    type(xmlf_t), intent(in) :: xf
    character(len=*), intent(in) :: comment
    integer :: ierr
    logical, save :: first=.true.
    !
    if ( xf%unit == -1 ) then
       print *, 'xml file not opened'
    else
       ! dirty  trick to have the same format with no changes to qexsd.f90
       if ( first .and. opentag /= '') then
          call xmlw_opentag ( opentag, ierr )
          if ( ierr /= 0 ) print *, 'xml_addcharacter: ierr = ', ierr
          opentag=''
          first = .false.
       end if
       write( xf%unit, '("<!-- ",A," -->")' ) trim(comment)
    end if
    !
  end subroutine xml_addcomment
  !
END MODULE wxml
