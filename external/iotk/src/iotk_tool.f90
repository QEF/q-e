! Input/Output Tool Kit (IOTK)
! Copyright (C) 2006 Giovanni Bussi
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

#include "iotk_auxmacros.h"


subroutine iotk_tool_x(args)
  use iotk_base
  use iotk_error_interf
  use iotk_str_interf
  use iotk_tool_interf
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  character(len=*), intent(in) :: args(:)
  integer :: iarg,ierrl
  character(iotk_linlenx) :: arg
  logical :: print_help_options,print_help_commands,print_help_basic,print_copyright,print_version
  logical :: check
  integer :: linlen,indent,maxindent
  ierrl = 0
  iarg = 1

  print_version = .false.
  print_help_options  = .false.
  print_help_commands = .false.
  print_help_basic = .false.
  print_copyright = .false.

  if(size(args)==0) then
    print_help_basic = .true.
  end if

  do iarg = 1 , size(args)
    arg = args(iarg)
    if(iotk_strcomp(arg(1:1),"-")) then
! options here
      if(iotk_strcomp(arg,"--help") .or. iotk_strcomp(arg,"-H")) then
        print_help_basic = .true.
        exit
      else if(iotk_strcomp(arg,"--version")) then
        print_version = .true.
        exit
      else if(iotk_strcomp(arg,"--do-nothing")) then
        exit
      else if(iotk_strcomp(arg,"--copyright")) then
        print_copyright = .true.
        exit
      else if(iotk_strcomp(arg,"--help-options")) then
        print_help_options = .true.
        exit
      else if(iotk_strcomp(arg,"--help-commands")) then
        print_help_commands = .true.
        exit
      else if(arg(1:13)=="--set-linlen=") then
        call iotk_atoi(linlen,arg(14:iotk_strlen(arg)),check=check)
        if(.not.check) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'')
          goto 1
        end if
        call iotk_set(linlen=linlen,ierr=ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'')
          goto 1
        end if
      else if(arg(1:13)=="--set-indent=") then
        call iotk_atoi(indent,arg(14:iotk_strlen(arg)),check=check)
        if(.not.check) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'')
          goto 1
        end if
        call iotk_set(indent=indent,ierr=ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'')
          goto 1
        end if
      else if(arg(1:16)=="--set-maxindent=") then
        call iotk_atoi(maxindent,arg(17:iotk_strlen(arg)),check=check)
        if(.not.check) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'')
          goto 1
        end if
        call iotk_set(maxindent=maxindent,ierr=ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'')
          goto 1
        end if
      else
        write(iotk_error_unit,"(a)") "unrecognized option `"//arg(1:iotk_strlen(arg))//"'"
        print_help_basic = .true.
        exit
      end if
    else
! commands here
      if(iotk_strcomp(arg,"convert")) then
        call iotk_tool_convert(args(iarg+1:),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'Error converting file')
          goto 1
        end if
      else if(iotk_strcomp(arg,"dump")) then
        call iotk_tool_dump(args(iarg+1:),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'Error dumping file')
          goto 1
        end if
      else if(iotk_strcomp(arg,"info")) then
        call iotk_tool_info(args(iarg+1:),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'Error')
          goto 1
        end if
      else if(iotk_strcomp(arg,"man")) then
        call iotk_tool_man(args(iarg+1:),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'Error')
          goto 1
        end if
      else
        write(iotk_error_unit,"(a)") "Unknown command: `"//arg(1:iotk_strlen(arg))//"'"
        write(iotk_error_unit,"(a)") ""
        print_help_commands = .true.
      end if
      exit
    end if
  end do

  if(print_help_basic) then
    write(iotk_error_unit,"(a)") "Usage: iotk [iotk-options] command [command-options-and-arguments]"
    write(iotk_error_unit,"(a)") "  where iotk-options are ..."
    write(iotk_error_unit,"(a)") "    (specify --help-options for a list of options)"
    write(iotk_error_unit,"(a)") "  where command is convert, dump, etc."
    write(iotk_error_unit,"(a)") "    (specify --help-commands for a list of commands)"
    write(iotk_error_unit,"(a)") "  where command-options-and-arguments depend on the specific command"
    write(iotk_error_unit,"(a)") "    (specify a command followed by --help for command-specific help)"
    write(iotk_error_unit,"(a)") "  Specify --help to receive this message"
  end if

  if(print_help_commands) then
    write(iotk_error_unit,"(a)") "IOTK commands are:"
    write(iotk_error_unit,"(a)") "  convert    to convert a file"
    write(iotk_error_unit,"(a)") "  dump       to dump a file"
    write(iotk_error_unit,"(a)") "  info       to obtain informations about how iotk was compiled"
    write(iotk_error_unit,"(a)") "  man        to print manual pages"
  end if

  if(print_help_options) then
    write(iotk_error_unit,"(a)") "IOTK options are:"
    write(iotk_error_unit,"(a)") "  --iotk-exe EXE     set the full path of iotk.x executable (first option)"
    write(iotk_error_unit,"(a)") "  --copyright        print copyright informations"
    write(iotk_error_unit,"(a)") "  --version          print version informations"
    write(iotk_error_unit,"(a)") "  --help             print a short, generic help"
    write(iotk_error_unit,"(a)") "  --help-options     print a list of options (this list)"
    write(iotk_error_unit,"(a)") "  --help-commands    print a list of commands"
    write(iotk_error_unit,"(a)") "  --set-linlen=N     to set the length of an output line"
    write(iotk_error_unit,"(a)") "  --set-indent=N     to set the number of spaces for an indent level"
    write(iotk_error_unit,"(a)") "  --set-maxindent=N  to set the maximum number of spaces when indenting"
  end if

  if(print_version) then
    write(*,"(a)") "Input/Output Tool Kit (IOTK) version: "//trim(iotk_version)
  end if

  if(print_copyright) then
    write(iotk_error_unit,"(a)") "Input/Output Tool Kit (IOTK)"
    write(iotk_error_unit,"(a)") "Copyright (C) 2004-2006 Giovanni Bussi"
    write(iotk_error_unit,"(a)") ""
    write(iotk_error_unit,"(a)") "This library is free software; you can redistribute it and/or"
    write(iotk_error_unit,"(a)") "modify it under the terms of the GNU Lesser General Public"
    write(iotk_error_unit,"(a)") "License as published by the Free Software Foundation; either"
    write(iotk_error_unit,"(a)") "version 2.1 of the License, or (at your option) any later version."
    write(iotk_error_unit,"(a)") ""
    write(iotk_error_unit,"(a)") "This library is distributed in the hope that it will be useful,"
    write(iotk_error_unit,"(a)") "but WITHOUT ANY WARRANTY; without even the implied warranty of"
    write(iotk_error_unit,"(a)") "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU"
    write(iotk_error_unit,"(a)") "Lesser General Public License for more details."
    write(iotk_error_unit,"(a)") ""
    write(iotk_error_unit,"(a)") "You should have received a copy of the GNU Lesser General Public"
    write(iotk_error_unit,"(a)") "License along with this library; if not, write to the Free Software"
    write(iotk_error_unit,"(a)") "Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA"
  end if

1 continue
  if(ierrl/=0) call iotk_error_handler(ierrl)

end subroutine iotk_tool_x

subroutine iotk_tool_convert_x(args,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_str_interf
  use iotk_misc_interf
  use iotk_files_interf
  implicit none
  character(len=*),           intent(in)  :: args(:)
  integer,          optional, intent(out) :: ierr
  integer :: iarg,ierrl,outfile_len
  character(len=iotk_fillenx) :: infile,outfile
  logical :: binary
  character(len=iotk_attlenx) :: attr
  character(len=iotk_taglenx) :: root
  integer :: maxsize
  logical :: autofmt
  infile=""
  outfile=""
  binary=.true.
  maxsize=-1
  ierrl = 0
  autofmt = .true.
  do iarg = 1 , size(args)
    if(iotk_strcomp(args(iarg)(1:1),"-")) then
      if(iotk_strcomp(args(iarg),"--help")) then
        write(iotk_error_unit,"(a)") "Usage: iotk convert [OPTIONS] infile outfile"
        write(iotk_error_unit,"(a)") "Options:"
        write(iotk_error_unit,"(a)") "  --mode=X  set the output file to be X, where X can be"
        write(iotk_error_unit,"(a)") "            'textual', 'binary' or 'auto'."
        write(iotk_error_unit,"(a)") "  -b        equivalent to --mode=binary"
        write(iotk_error_unit,"(a)") "  -t        equivalent to --mode=textual"
        write(iotk_error_unit,"(a)") "This command converts a iotk data file into another iotk data file."
        write(iotk_error_unit,"(a)") "The infile can be textual or binary, and its format is automatically detected."
        write(iotk_error_unit,"(a)") "The outfile can be textual or binary depending on the --mode option."
        write(iotk_error_unit,"(a)") "If the mode is 'auto', the decision is driven by outfile extension,"
        write(iotk_error_unit,"(a)") "i.e. a file matching *.txt of *.xml will be considered textual, otherwise binary"
        goto 1
      else if(iotk_strcomp(args(iarg),"-b") .or. iotk_strcomp(args(iarg),"--mode=binary")) then
        binary = .true.
        autofmt = .false.
      else if(iotk_strcomp(args(iarg),"-t") .or. iotk_strcomp(args(iarg),"--mode=textual")) then
        binary = .false.
        autofmt = .false.
      else if(iotk_strcomp(args(iarg),"--mode=auto")) then
        binary = .true.
        autofmt = .true.
      else
        call iotk_error_issue(ierrl,"iotk_tool_convert",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'Unknown option')
        goto 1
      end if
    else
      if(infile=="") then
        call iotk_strcpy(infile,args(iarg),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool_convert",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'File name too long')
          goto 1
        end if
      else if(outfile=="") then
        call iotk_strcpy(outfile,args(iarg),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool_convert",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'File name too long')
          goto 1
        end if
      else
        call iotk_error_issue(ierrl,"iotk_tool_convert",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'Three files. What do you mean?')
        goto 1
      end if
    end if
  end do
  if(outfile=="") then
    call iotk_error_issue(ierrl,"iotk_tool_convert",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'Convert: bad usage')
    goto 1
  end if

  outfile_len = iotk_strlen(outfile)
  if(outfile_len>3) then
    select case(outfile(outfile_len-3:outfile_len))
    case(".xml")
      binary = .false.
    case(".txt")
      binary = .false.
    case default
      binary = .true.
    end select
  end if

  call iotk_open_read(60,infile,root=root,attr=attr)
  call iotk_open_write(61,outfile,binary=binary,root=root,attr=attr)
  call iotk_copy_tag(60,61,maxsize=-1)
  call iotk_close_write(61)
  call iotk_close_read(60)

1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_tool_convert_x


subroutine iotk_tool_dump_x(args,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_str_interf
  use iotk_misc_interf
  use iotk_files_interf
  implicit none
  character(len=*),           intent(in)  :: args(:)
  integer,          optional, intent(out) :: ierr
  integer :: iarg,ierrl
  character(len=iotk_fillenx) :: infile
  character(len=iotk_attlenx) :: attr
  character(len=iotk_taglenx) :: root
  integer :: maxsize
  infile=""
  maxsize=-1
  ierrl = 0
  do iarg = 1 , size(args)
    if(iotk_strcomp(args(iarg)(1:1),"-")) then
      if(iotk_strcomp(args(iarg),"--help")) then
        write(iotk_error_unit,"(a)") "Usage: iotk dump file"
        write(iotk_error_unit,"(a)") "This command dumps a iotk data file on standard out."
        write(iotk_error_unit,"(a)") "The file can be textual or binary, and its format is automatically detected."
        goto 1
      else
        call iotk_error_issue(ierrl,"iotk_tool_dump",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'Unknown option')
        goto 1
      end if
    else
      if(infile=="") then
        call iotk_strcpy(infile,args(iarg),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool_dump",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'File name too long')
          goto 1
        end if
      else
        call iotk_error_issue(ierrl,"iotk_tool_dump",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'Two files. What do you mean?')
        goto 1
      end if
    end if
  end do

  call iotk_open_read(60, trim(infile),root=root,attr=attr)
  call iotk_open_write(iotk_output_unit,root=root,attr=attr)
  call iotk_copy_tag(60,iotk_output_unit,maxsize=-1)
  call iotk_close_write(iotk_output_unit)
  call iotk_close_read(60)

1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_tool_dump_x

subroutine iotk_tool_info_x(args,ierr)
  use iotk_base
  use iotk_misc_interf
  use iotk_xtox_interf
  use iotk_error_interf
  implicit none
  character(len=*),           intent(in)  :: args(:)
  integer,          optional, intent(out) :: ierr
  integer :: ierrl
  ierrl = 0
  write(*,"(a)") "IOTK (Input/Output Tool Kit) version: "//trim(iotk_version)
  write(*,"(a)") "Limits:"
  write(*,"(a)") "  maximum rank (soft limit): "//trim(iotk_itoa(iotk_maxrank))
  write(*,"(a)") "  maximum rank (hard limit): "//trim(iotk_itoa(iotk_maxrank_hard))
  write(*,"(a)") "Special kinds:"
  write(*,"(a)") "  headers in binary files are integer(kind="//trim(iotk_itoa(iotk_header_kind))//")"
  write(*,"(a)") "  default integers are integer(kind="//trim(iotk_itoa(iotk_integer_defkind))//")"
  write(*,"(a)") "  default logicals are logical(kind="//trim(iotk_itoa(iotk_logical_defkind))//")"
  write(*,"(a)") "  default characters are character(kind="//trim(iotk_itoa(iotk_character_defkind))//")"
  write(*,"(a)") "Kinds configured for i/o operations:"
#ifdef __IOTK_LOGICAL1
  write(*,"(a)") "  logical(kind="//trim(iotk_itoa(iotk_logical1))//")"
#endif
#ifdef __IOTK_LOGICAL2
  write(*,"(a)") "  logical(kind="//trim(iotk_itoa(iotk_logical2))//")"
#endif
#ifdef __IOTK_LOGICAL3
  write(*,"(a)") "  logical(kind="//trim(iotk_itoa(iotk_logical3))//")"
#endif
#ifdef __IOTK_LOGICAL4
  write(*,"(a)") "  logical(kind="//trim(iotk_itoa(iotk_logical4))//")"
#endif
#ifdef __IOTK_INTEGER1
  write(*,"(a)") "  integer(kind="//trim(iotk_itoa(iotk_integer1))//")"
#endif
#ifdef __IOTK_INTEGER2
  write(*,"(a)") "  integer(kind="//trim(iotk_itoa(iotk_integer2))//")"
#endif
#ifdef __IOTK_INTEGER3
  write(*,"(a)") "  integer(kind="//trim(iotk_itoa(iotk_integer3))//")"
#endif
#ifdef __IOTK_INTEGER4
  write(*,"(a)") "  integer(kind="//trim(iotk_itoa(iotk_integer4))//")"
#endif
#ifdef __IOTK_REAL1
  write(*,"(a)") "  real(kind="//trim(iotk_itoa(iotk_real1))//")"
#endif
#ifdef __IOTK_REAL2
  write(*,"(a)") "  real(kind="//trim(iotk_itoa(iotk_real2))//")"
#endif
#ifdef __IOTK_REAL3
  write(*,"(a)") "  real(kind="//trim(iotk_itoa(iotk_real3))//")"
#endif
#ifdef __IOTK_REAL4
  write(*,"(a)") "  real(kind="//trim(iotk_itoa(iotk_real4))//")"
#endif
#ifdef __IOTK_REAL1
  write(*,"(a)") "  complex(kind="//trim(iotk_itoa(iotk_real1))//")"
#endif
#ifdef __IOTK_REAL2
  write(*,"(a)") "  complex(kind="//trim(iotk_itoa(iotk_real2))//")"
#endif
#ifdef __IOTK_REAL3
  write(*,"(a)") "  complex(kind="//trim(iotk_itoa(iotk_real3))//")"
#endif
#ifdef __IOTK_REAL4
  write(*,"(a)") "  complex(kind="//trim(iotk_itoa(iotk_real4))//")"
#endif
  write(*,"(a)") "  character(kind="//trim(iotk_itoa(iotk_character1))//")"

1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_tool_info_x

subroutine iotk_tool_man_x(args,ierr)
  use iotk_base
  use iotk_misc_interf
  use iotk_xtox_interf
  use iotk_error_interf
  use iotk_str_interf
  implicit none
  character(len=*),           intent(in)  :: args(:)
  integer,          optional, intent(out) :: ierr
  character(len=iotk_linlenx) :: keyword
  integer :: ierrl,iarg
  logical :: printme,printlist

  ierrl = 0
  printme = .false.
  printlist = .false.
  keyword(1:1) = iotk_eos

  do iarg = 1 , size(args)
    if(iotk_strcomp(args(iarg)(1:1),"-")) then
      if(iotk_strcomp(args(iarg),"--help")) then
        write(iotk_error_unit,"(a)") "Usage: iotk man [keyword]"
        write(iotk_error_unit,"(a)") "This command prints on stdout the page of the built-in manual associated with the keyword."
        write(iotk_error_unit,"(a)") "If the keyword is not given a list of all the available keywords will be printed."
        goto 1
      else
        call iotk_error_issue(ierrl,"iotk_tool_dump",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'Unknown option')
        goto 1
      end if
    else
      if(iotk_strcomp(keyword,"")) then
        call iotk_strcpy(keyword,args(iarg),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool_dump",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
          goto 1
        end if
      else
        call iotk_error_issue(ierrl,"iotk_tool_dump",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.18 ")
call iotk_error_msg(ierrl,'Two keywords. What do you mean?')
        goto 1
      end if
    end if
  end do

  if(iotk_strcomp(keyword,"")) then
    write(iotk_output_unit,"(a)") "List of available pages:"
    printlist = .true.
  end if

  if(printlist) write(iotk_output_unit,"(a)") &
" intro"
printme=.false.
if(iotk_strcomp(keyword,"all")) printme=.true.
if(iotk_strcomp(keyword,'intro')) printme=.true.
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"IOTK: INTRODUCTION"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"The input/output tool kit (IOTK) is a FORTRAN90 library intended"
if(printme) write(iotk_output_unit,"(a)") &
"to provide an easy access to tagged files formatted using some specific rule."
if(printme) write(iotk_output_unit,"(a)") &
"In this context, a tagged file is a file containing tags and data."
if(printme) write(iotk_output_unit,"(a)") &
"Tagged files can be textual, in which case a XML-like format is used,"
if(printme) write(iotk_output_unit,"(a)") &
"or binary, in which case a special format is used."
if(printme) write(iotk_output_unit,"(a)") &
"Note that IOTK is not an XML parser, but it can be used as a writer/parser"
if(printme) write(iotk_output_unit,"(a)") &
"for a limited subset of the XML language."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"To use the IOTK library from a FORTRAN90 source, the user should"
if(printme) write(iotk_output_unit,"(a)") &
"use the module 'iotk_module'."
if(printme) write(iotk_output_unit,"(a)") &
"To minimize the possibility of name clashes, all public names exported"
if(printme) write(iotk_output_unit,"(a)") &
"from this module has the "//'"'//&
"iotk_"//'"'//&
" prefix."
if(printme) write(iotk_output_unit,"(a)") &
"Communication between user and library is based on"
if(printme) write(iotk_output_unit,"(a)") &
"integers, characters and logicals of the default kind (notice that"
if(printme) write(iotk_output_unit,"(a)") &
"these kinds can be changed using proper compiler options, so that"
if(printme) write(iotk_output_unit,"(a)") &
"the actual kinds depend on how the library was compiled on your machine)."
if(printme) write(iotk_output_unit,"(a)") &
"However, the library can handle formatted input/output for"
if(printme) write(iotk_output_unit,"(a)") &
"all intrinsic datatypes, kinds and ranks if properly configured."
if(printme) write(iotk_output_unit,"(a)") &
"This is obtained interfacing procedures which acts on all kinds,"
if(printme) write(iotk_output_unit,"(a)") &
"types and (in almost all cases) ranks. Thus, a single generic"
if(printme) write(iotk_output_unit,"(a)") &
"name has to be remembered for each subroutine, and the compiler will"
if(printme) write(iotk_output_unit,"(a)") &
"link the correct one dependening on type, kind and rank of the arguments."
if(printme) write(iotk_output_unit,"(a)") &
"Backward API compatibility will be mantained (as long as it is possible)"
if(printme) write(iotk_output_unit,"(a)") &
"in future versions."
if(printme) write(iotk_output_unit,"(a)") &
"Backward file compatibility will be mantained (as long as it is possible) in"
if(printme) write(iotk_output_unit,"(a)") &
"future versions."
if(printme) write(iotk_output_unit,"(a)") &
"The library writes on files informations about the version of the library."
if(printme) write(iotk_output_unit,"(a)") &
"It also writes informations about the version of the file format (file_version)."
if(printme) write(iotk_output_unit,"(a)") &
"The later has to be older or equal to the format supported in the actual library."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printlist) write(iotk_output_unit,"(a)") &
" error_handling"
printme=.false.
if(iotk_strcomp(keyword,"all")) printme=.true.
if(iotk_strcomp(keyword,'error_handling')) printme=.true.
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"IOTK: ERROR HANDLING"
if(printme) write(iotk_output_unit,"(a)") &
"The way iotk handles error is sophisticated and allows for a trace back"
if(printme) write(iotk_output_unit,"(a)") &
"of the error condition inside the library."
if(printme) write(iotk_output_unit,"(a)") &
"Every iotk routines which possibly leads to an error condition has an optional"
if(printme) write(iotk_output_unit,"(a)") &
"intent(out) integer argument ierr. The returned value is conventionally"
if(printme) write(iotk_output_unit,"(a)") &
"0 when the routine returns correctly, and different from 0 when the routines"
if(printme) write(iotk_output_unit,"(a)") &
"raise an error. The value is effectively a handler for a more complex"
if(printme) write(iotk_output_unit,"(a)") &
"object containing the error message. When an error is raised in a low-level"
if(printme) write(iotk_output_unit,"(a)") &
"iotk routine, a message is written on the error object. Any intermediate routine"
if(printme) write(iotk_output_unit,"(a)") &
"can add other messages to the error object, at least the number of the line in"
if(printme) write(iotk_output_unit,"(a)") &
"the source file. In this way, the error message contains a complete trace of"
if(printme) write(iotk_output_unit,"(a)") &
"the error plus some additional information."
if(printme) write(iotk_output_unit,"(a)") &
"At any point in the chain the messages can be exctracted from the error object."
if(printme) write(iotk_output_unit,"(a)") &
"At some point in the chain the error is really handled, usually by writing the"
if(printme) write(iotk_output_unit,"(a)") &
"message on the appropriate unit and aborting the execution."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"Scanning routines (iotk_scan_*) have an optional logical argument "//'"'//&
"found"//'"'//&
""
if(printme) write(iotk_output_unit,"(a)") &
"which returns true or false. When scanning for data, also a "//'"'//&
"default"//'"'//&
" argument"
if(printme) write(iotk_output_unit,"(a)") &
"can be used. If one of these two argument is present, the searched"
if(printme) write(iotk_output_unit,"(a)") &
"object is considered as an optional object. Otherwise, it is considered as a needed object."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"If the ierr optional argument is absent, the error handling is leaved to the iotk library."
if(printme) write(iotk_output_unit,"(a)") &
"In this case, if a needed object is not present, the library handles the error with a"
if(printme) write(iotk_output_unit,"(a)") &
"forced stop."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"If the ierr optional argument is present, it returns an error code."
if(printme) write(iotk_output_unit,"(a)") &
"ierr = 0 means that no error has occurred"
if(printme) write(iotk_output_unit,"(a)") &
"ierr > 0 means that an error has occurred probably related to file corruption"
if(printme) write(iotk_output_unit,"(a)") &
"ierr < 0 means that the item that was searched for has not been found"
if(printme) write(iotk_output_unit,"(a)") &
"(it is possible only for scanning routines and only if the"
if(printme) write(iotk_output_unit,"(a)") &
"found and the default keywords are both missing, i.e. only for no-optional objetcs)"
if(printme) write(iotk_output_unit,"(a)") &
"In scanning routines, if the argument "//'"'//&
"found"//'"'//&
" is present it returns .true."
if(printme) write(iotk_output_unit,"(a)") &
"if the item has been found, .false. otherwise."
if(printme) write(iotk_output_unit,"(a)") &
"If a library routine returns an ierr /= 0 it is STRONGLY RECOMMENDED to"
if(printme) write(iotk_output_unit,"(a)") &
"clear that error with "//'"'//&
"call iotk_error_clear(ierr)"//'"'//&
" before proceeding."
if(printme) write(iotk_output_unit,"(a)") &
"Thus, the final recipe is:"
if(printme) write(iotk_output_unit,"(a)") &
"* if you want to handle errors, always use the 'ierr' optional argument."
if(printme) write(iotk_output_unit,"(a)") &
"looking at the sign, you will discern between lacking data and file corruption."
if(printme) write(iotk_output_unit,"(a)") &
"with iotk_error_print you can obtain a description of the error."
if(printme) write(iotk_output_unit,"(a)") &
"* if you want to leave the error handling to the library, don't use"
if(printme) write(iotk_output_unit,"(a)") &
"the 'ierr' optional argument."
if(printme) write(iotk_output_unit,"(a)") &
"- if the object you are searching is optional, use 'found' or 'default' optional arguments."
if(printme) write(iotk_output_unit,"(a)") &
"- if the object you are searching is non-optional, don't use 'found' nor 'default' optional arguments."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printlist) write(iotk_output_unit,"(a)") &
" binary_and_textual_files"
printme=.false.
if(iotk_strcomp(keyword,"all")) printme=.true.
if(iotk_strcomp(keyword,'binary_and_textual_files')) printme=.true.
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"IOTK: BINARY AND TEXTUAL FILES"
if(printme) write(iotk_output_unit,"(a)") &
"Units can be opened on textual or binary files."
if(printme) write(iotk_output_unit,"(a)") &
"The word 'binary' is used instead of the fortran 'unformatted' since"
if(printme) write(iotk_output_unit,"(a)") &
"using this libray also binary files have a degree of formattation."
if(printme) write(iotk_output_unit,"(a)") &
"After a unit has been opened, the library automatically detects"
if(printme) write(iotk_output_unit,"(a)") &
"its format through an INQUIRE and acts consequently."
if(printme) write(iotk_output_unit,"(a)") &
"Note that the iotk routines check for necessary properties of an opened unit"
if(printme) write(iotk_output_unit,"(a)") &
"access="//'"'//&
"sequential"//'"'//&
""
if(printme) write(iotk_output_unit,"(a)") &
"blank ="//'"'//&
"null"//'"'//&
" (only textual i/o)"
if(printme) write(iotk_output_unit,"(a)") &
"pad   ="//'"'//&
"yes"//'"'//&
"  (only textual i/o)"
if(printme) write(iotk_output_unit,"(a)") &
"Moreover, a textual or binary unit can be designed as raw."
if(printme) write(iotk_output_unit,"(a)") &
"In that case, no tags are placed on the file and everything"
if(printme) write(iotk_output_unit,"(a)") &
"has to be read and written in the same order."
if(printme) write(iotk_output_unit,"(a)") &
"This feature is provided for compatibility reasons but it should be"
if(printme) write(iotk_output_unit,"(a)") &
"used as few as possible."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printlist) write(iotk_output_unit,"(a)") &
" optional_arguments"
printme=.false.
if(iotk_strcomp(keyword,"all")) printme=.true.
if(iotk_strcomp(keyword,'optional_arguments')) printme=.true.
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"IOTK: OPTIONAL ARGUMENTS"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"Most iotk routines accept optional arguments."
if(printme) write(iotk_output_unit,"(a)") &
"The calling routine will not compile if the names of the"
if(printme) write(iotk_output_unit,"(a)") &
"arguments are not indicated.  For instance, use"
if(printme) write(iotk_output_unit,"(a)") &
"call iotk_scan_dat(10,"//'"'//&
"pippo"//'"'//&
",aa(:),ierr=ii)"
if(printme) write(iotk_output_unit,"(a)") &
"and NOT:"
if(printme) write(iotk_output_unit,"(a)") &
"call iotk_scan_dat(10,"//'"'//&
"pippo"//'"'//&
",aa(:),ii)"
if(printme) write(iotk_output_unit,"(a)") &
"The only exeption is the attr argument, for which the name can be"
if(printme) write(iotk_output_unit,"(a)") &
"omitted if it is placed as the first of the optional arguments."
if(printme) write(iotk_output_unit,"(a)") &
"In any case, it is better to always explicitly label optional arguments."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printlist) write(iotk_output_unit,"(a)") &
" basic_writing_routines iotk_write_begin iotk_write_end iotk_write_empty iotk_write_pi iotk_write_comment"
printme=.false.
if(iotk_strcomp(keyword,"all")) printme=.true.
if(iotk_strcomp(keyword,'basic_writing_routines')) printme=.true.
if(iotk_strcomp(keyword,'iotk_write_begin')) printme=.true.
if(iotk_strcomp(keyword,'iotk_write_end')) printme=.true.
if(iotk_strcomp(keyword,'iotk_write_empty')) printme=.true.
if(iotk_strcomp(keyword,'iotk_write_pi')) printme=.true.
if(iotk_strcomp(keyword,'iotk_write_comment')) printme=.true.
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"IOTK: BASIC WRITING ROUTINES"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_write_begin  (unit,name[,attr][,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_write_end    (unit,name[,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_write_empty  (unit,name[,attr][,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_write_pi     (unit,name[,attr][,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_write_comment(unit,text[,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(in) :: unit"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in) :: name"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in) :: text"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in) :: attr"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(out):: ierr ! see error_handling page"
if(printme) write(iotk_output_unit,"(a)") &
"These routines write a tag named 'name' on fortran unit 'unit'."
if(printme) write(iotk_output_unit,"(a)") &
"The type of the tag is determined from the name of the routine:"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_write_begin   => <name attr>"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_write_end     => </name>"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_write_empty   => <name attr/>"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_write_pi      => <?name attr?>"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_write_comment => <!--text-->"
if(printme) write(iotk_output_unit,"(a)") &
"An optional attribute string can be supplied in 'attr'"
if(printme) write(iotk_output_unit,"(a)") &
"In end tags, no attribute is allowed."
if(printme) write(iotk_output_unit,"(a)") &
"To build the attribute string, use iotk_write_attr."
if(printme) write(iotk_output_unit,"(a)") &
"DON'T TRY TO MANIPULATE THE ATTRIBUTE STRING DIRECTLY!"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printlist) write(iotk_output_unit,"(a)") &
" basic_scanning_routines iotk_scan_begin iotk_scan_end iotk_scan_empty iotk_scan_pi"
printme=.false.
if(iotk_strcomp(keyword,"all")) printme=.true.
if(iotk_strcomp(keyword,'basic_scanning_routines')) printme=.true.
if(iotk_strcomp(keyword,'iotk_scan_begin')) printme=.true.
if(iotk_strcomp(keyword,'iotk_scan_end')) printme=.true.
if(iotk_strcomp(keyword,'iotk_scan_empty')) printme=.true.
if(iotk_strcomp(keyword,'iotk_scan_pi')) printme=.true.
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"IOTK: BASIC SCANNING ROUTINES"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_scan_begin(unit,name[,attr][,found][,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_scan_end  (unit,name[,found][,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_scan_empty(unit,name[,attr][,found][,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_scan_pi   (unit,name[,attr][,found][,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(in) :: unit"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in) :: name  ! len less or equal iotk_namlenx"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(out):: attr  ! len possibily equal iotk_attlenx"
if(printme) write(iotk_output_unit,"(a)") &
"logical,          intent(out):: found ! see error_handling page"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(out):: ierr  ! see error_handling page"
if(printme) write(iotk_output_unit,"(a)") &
"These routines scan for a tag named 'name' on fortran unit 'unit'."
if(printme) write(iotk_output_unit,"(a)") &
"The type of the tag is determined from the name of the routine:"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_scan_begin => <name attr>"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_scan_end   => </name>"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_scan_empty => <name attr/>"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_scan_pi    => <?name attr?>"
if(printme) write(iotk_output_unit,"(a)") &
"These routines (except for iotk_scan_end) also fills the"
if(printme) write(iotk_output_unit,"(a)") &
"attr string, which can be subsequently decoded with iotk_scan_attr."
if(printme) write(iotk_output_unit,"(a)") &
"DON'T TRY TO MANIPULATE THE ATTRIBUTE STRING DIRECTLY!"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printlist) write(iotk_output_unit,"(a)") &
" writing_attributes iotk_write_attr"
printme=.false.
if(iotk_strcomp(keyword,"all")) printme=.true.
if(iotk_strcomp(keyword,'writing_attributes')) printme=.true.
if(iotk_strcomp(keyword,'iotk_write_attr')) printme=.true.
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"IOTK: WRITING ATTRIBUTES"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_write_attr (attr,name,val[,first][,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(out):: attr  ! len less or equal iotk_namlenx"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in) :: name  ! len less or equal iotk_attlenx"
if(printme) write(iotk_output_unit,"(a)") &
"TYPE(KIND),       intent(in) :: val   ! any type, any kind, any rank [but only scalars for character]"
if(printme) write(iotk_output_unit,"(a)") &
"logical,          intent(in) :: first"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(out):: ierr"
if(printme) write(iotk_output_unit,"(a)") &
"This routine adds one attribute to the 'attr' string."
if(printme) write(iotk_output_unit,"(a)") &
"To clean the string (for the first attribute) use first=.true."
if(printme) write(iotk_output_unit,"(a)") &
"Example:"
if(printme) write(iotk_output_unit,"(a)") &
"call iotk_write_attr(attr,"//'"'//&
"pippo"//'"'//&
",1,first=.true.)"
if(printme) write(iotk_output_unit,"(a)") &
"call iotk_write_attr(attr,"//'"'//&
"paperino"//'"'//&
",2)"
if(printme) write(iotk_output_unit,"(a)") &
"call iotk_write_attr(attr,"//'"'//&
"pluto"//'"'//&
",3)"
if(printme) write(iotk_output_unit,"(a)") &
"This is equivalent to attr="//'"'//&
""//'"'//&
" before the call, but more efficient."
if(printme) write(iotk_output_unit,"(a)") &
"The attribute is added in the form name="//'"'//&
"value"//'"'//&
","
if(printme) write(iotk_output_unit,"(a)") &
"where "//'"'//&
"value"//'"'//&
" is a string containing a textual representation"
if(printme) write(iotk_output_unit,"(a)") &
"of the val variable."
if(printme) write(iotk_output_unit,"(a)") &
"If one of <>&"//'"'//&
"' appears in val, it is automatically escaped."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printlist) write(iotk_output_unit,"(a)") &
" scanning_attributes iotk_scan_attr"
printme=.false.
if(iotk_strcomp(keyword,"all")) printme=.true.
if(iotk_strcomp(keyword,'scanning_attributes')) printme=.true.
if(iotk_strcomp(keyword,'iotk_scan_attr')) printme=.true.
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"IOTK: SCANNING ATTRIBUTES"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_scan_attr  (attr,name,val[,found][,default][,eos][,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in) :: attr    ! len possibily equal iotk_attlenx"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in) :: name    ! len less or equal iotk_namlenx"
if(printme) write(iotk_output_unit,"(a)") &
"TYPE(KIND),       intent(out):: val     ! any type, any kind, any rank [but only scalars for character]"
if(printme) write(iotk_output_unit,"(a)") &
"logical,          intent(out):: found   ! see error_handling page"
if(printme) write(iotk_output_unit,"(a)") &
"TYPE(KIND),       intent(in) :: default ! same type, kind and rank as val"
if(printme) write(iotk_output_unit,"(a)") &
"logical,          intent(in) :: eos"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(out):: ierr    ! see error_handling page"
if(printme) write(iotk_output_unit,"(a)") &
"This routine scans for one attribute named 'name' from the 'attr' string."
if(printme) write(iotk_output_unit,"(a)") &
"If the attribute is found, it is read to variable 'val'."
if(printme) write(iotk_output_unit,"(a)") &
"If it is not found and default is present, default is copied onto val."
if(printme) write(iotk_output_unit,"(a)") &
"If TYPE is character and eos is present and true,"
if(printme) write(iotk_output_unit,"(a)") &
"an end-of-string terminator will be attached at the end of the read string,"
if(printme) write(iotk_output_unit,"(a)") &
"and the following bytes will not be touched. This is faster, but requires"
if(printme) write(iotk_output_unit,"(a)") &
"the user to take care directly of the end-of-string. Thus, it is discouraged."
if(printme) write(iotk_output_unit,"(a)") &
"The attribute can be delimited with "//'"'//&
""//'"'//&
" or with ''"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printlist) write(iotk_output_unit,"(a)") &
" writing_data iotk_write_dat"
printme=.false.
if(iotk_strcomp(keyword,"all")) printme=.true.
if(iotk_strcomp(keyword,'writing_data')) printme=.true.
if(iotk_strcomp(keyword,'iotk_write_dat')) printme=.true.
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"IOTK: WRITING DATA"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_write_dat  (unit,name,dat[,fmt][,columns][,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(in) :: unit"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in) :: name    ! len less or equal iotk_namlenx"
if(printme) write(iotk_output_unit,"(a)") &
"TYPE(KIND),       intent(in) :: dat     ! any type, any kind, any rank"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in) :: fmt"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(in) :: columns"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(out):: ierr    ! see error_handling page"
if(printme) write(iotk_output_unit,"(a)") &
"This routines write a data object, that is a self-described"
if(printme) write(iotk_output_unit,"(a)") &
"object containg fortran data."
if(printme) write(iotk_output_unit,"(a)") &
"A single data object has the following form"
if(printme) write(iotk_output_unit,"(a)") &
"<name type="//'"'//&
"TYPE"//'"'//&
" kind="//'"'//&
"KIND"//'"'//&
" size="//'"'//&
"SIZE"//'"'//&
" columns="//'"'//&
"COLUMNS"//'"'//&
" len="//'"'//&
"LEN"//'"'//&
" fmt="//'"'//&
"FMT"//'"'//&
">"
if(printme) write(iotk_output_unit,"(a)") &
".. DATA ..."
if(printme) write(iotk_output_unit,"(a)") &
"</name>"
if(printme) write(iotk_output_unit,"(a)") &
"where"
if(printme) write(iotk_output_unit,"(a)") &
"TYPE    is the intrinsic type (logical,integer,real,complex or character),"
if(printme) write(iotk_output_unit,"(a)") &
"KIND    is the data kind (stored in binary files only)"
if(printme) write(iotk_output_unit,"(a)") &
"SIZE    is the array size (shape informations are not stored)"
if(printme) write(iotk_output_unit,"(a)") &
"COLUMNS is the number of data per line"
if(printme) write(iotk_output_unit,"(a)") &
"LEN     is the string length"
if(printme) write(iotk_output_unit,"(a)") &
"FMT     is a fortran format string used to write data"
if(printme) write(iotk_output_unit,"(a)") &
"If the optional 'fmt' is not passed, default format ('columns' element per line)"
if(printme) write(iotk_output_unit,"(a)") &
"is used and the fmt attribute is not written. Otherwise, the string"
if(printme) write(iotk_output_unit,"(a)") &
"fmt is used as a FORTRAN format specifierfor the write statement. In this"
if(printme) write(iotk_output_unit,"(a)") &
"case it is also written on the file (and used for reading the data back)."
if(printme) write(iotk_output_unit,"(a)") &
"fmt="//'"'//&
"*"//'"'//&
" can be used and correspond to the "//'"'//&
"write(unit,*)"//'"'//&
" statement."
if(printme) write(iotk_output_unit,"(a)") &
"If the optional 'columns' is not passed, it is assumed to be 1 and"
if(printme) write(iotk_output_unit,"(a)") &
"the columns attribute is not written. Note that this attribute is completely"
if(printme) write(iotk_output_unit,"(a)") &
"ininfluent when reading."
if(printme) write(iotk_output_unit,"(a)") &
"columns and fmt arguments are incompatible."
if(printme) write(iotk_output_unit,"(a)") &
"For complex data, one element is a couple of comma separated real numbers."
if(printme) write(iotk_output_unit,"(a)") &
"If one of <>& appears in dat, it is escaped."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printlist) write(iotk_output_unit,"(a)") &
" scanning_data iotk_scan_dat"
printme=.false.
if(iotk_strcomp(keyword,"all")) printme=.true.
if(iotk_strcomp(keyword,'scanning_data')) printme=.true.
if(iotk_strcomp(keyword,'iotk_scan_dat')) printme=.true.
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"IOTK: SCANNING DATA"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_scan_dat  (unit,name,dat[,found][,default][,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(in) :: unit"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in) :: name    ! len less or equal iotk_namlenx"
if(printme) write(iotk_output_unit,"(a)") &
"TYPE(KIND),       intent(out):: dat     ! any type, any kind, any rank"
if(printme) write(iotk_output_unit,"(a)") &
"logical,          intent(out):: found   ! see error_handling page"
if(printme) write(iotk_output_unit,"(a)") &
"TYPE(KIND),       intent(in) :: default ! same type, kind and rank as dat"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(out):: ierr    ! see error_handling page"
if(printme) write(iotk_output_unit,"(a)") &
"A data object written with iotk_write_dat is read."
if(printme) write(iotk_output_unit,"(a)") &
"If it is not found and default is present, default is copied onto dat."
if(printme) write(iotk_output_unit,"(a)") &
"If a keyword is absent in the file, the value is deduced from the"
if(printme) write(iotk_output_unit,"(a)") &
"dat argument and no check is performed. This allows to write"
if(printme) write(iotk_output_unit,"(a)") &
"rapidly by hand data objects. For instance"
if(printme) write(iotk_output_unit,"(a)") &
"<datum> 1.0 </datum>"
if(printme) write(iotk_output_unit,"(a)") &
"can be read with"
if(printme) write(iotk_output_unit,"(a)") &
"real :: val"
if(printme) write(iotk_output_unit,"(a)") &
"call iotk_scan_dat(unit,"//'"'//&
"datum"//'"'//&
",val)"
if(printme) write(iotk_output_unit,"(a)") &
"If fmt is not present on file, the default format is used."
if(printme) write(iotk_output_unit,"(a)") &
"Types and sizes are checked."
if(printme) write(iotk_output_unit,"(a)") &
"Different kinds (for binary i/o) are automatically converted."
if(printme) write(iotk_output_unit,"(a)") &
"Length (for characters) are not checked. If strings on files"
if(printme) write(iotk_output_unit,"(a)") &
"are longer then len(dat), only the first characters are read; if strings"
if(printme) write(iotk_output_unit,"(a)") &
"on files are shorter, dat is padded with blanks."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printlist) write(iotk_output_unit,"(a)") &
" opening_files iotk_open_write iotk_open_read"
printme=.false.
if(iotk_strcomp(keyword,"all")) printme=.true.
if(iotk_strcomp(keyword,'opening_files')) printme=.true.
if(iotk_strcomp(keyword,'iotk_open_write')) printme=.true.
if(iotk_strcomp(keyword,'iotk_open_read')) printme=.true.
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"IOTK: OPENING FILES"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"iotk_open_write(unit[,file][,attr][,binary][,raw][,new][,root][,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(in)  :: unit"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in)  :: file"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in)  :: attr"
if(printme) write(iotk_output_unit,"(a)") &
"logical,          intent(in)  :: binary"
if(printme) write(iotk_output_unit,"(a)") &
"logical,          intent(in)  :: new"
if(printme) write(iotk_output_unit,"(a)") &
"logical,          intent(in)  :: raw"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in)  :: root   ! len less or equal iotk_namlenx"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(out) :: ierr   ! see error_handling page"
if(printme) write(iotk_output_unit,"(a)") &
"If file is present, this routines opens file 'file' on"
if(printme) write(iotk_output_unit,"(a)") &
"unit 'unit' with the proper options."
if(printme) write(iotk_output_unit,"(a)") &
"If binary is present and true, the file is binary."
if(printme) write(iotk_output_unit,"(a)") &
"If new is present and true, the file must not exist already."
if(printme) write(iotk_output_unit,"(a)") &
"If raw is present and true, the file is considered as a raw data file"
if(printme) write(iotk_output_unit,"(a)") &
"(use of raw data files is discouraged)."
if(printme) write(iotk_output_unit,"(a)") &
"If file is not present, unit is assumed to be already connected."
if(printme) write(iotk_output_unit,"(a)") &
"If root is present, it is used as the name of the root begin/end tag."
if(printme) write(iotk_output_unit,"(a)") &
"If it is absent, the default "//'"'//&
"Root"//'"'//&
" is used."
if(printme) write(iotk_output_unit,"(a)") &
"An optional attribute string can be supplied in 'attr', and will be used"
if(printme) write(iotk_output_unit,"(a)") &
"as an attribute list for the begin root tag."
if(printme) write(iotk_output_unit,"(a)") &
"Also informations about iotk version and binary format are written as"
if(printme) write(iotk_output_unit,"(a)") &
"pi informations."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"iotk_open_read(unit[,file][,attr][,binary][,raw][,root][,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(in)  :: unit"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in)  :: file"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(out) :: attr"
if(printme) write(iotk_output_unit,"(a)") &
"logical,          intent(in)  :: binary"
if(printme) write(iotk_output_unit,"(a)") &
"logical,          intent(in)  :: raw"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(out) :: root   ! len possibly equal iotk_namlenx"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(out) :: ierr   ! see error_handling page"
if(printme) write(iotk_output_unit,"(a)") &
"If file is present, this routines opens file 'file' on"
if(printme) write(iotk_output_unit,"(a)") &
"unit 'unit' with the proper options."
if(printme) write(iotk_output_unit,"(a)") &
"If binary is present and true, the file is binary."
if(printme) write(iotk_output_unit,"(a)") &
"If raw is present and true, the file is considered as a raw data file"
if(printme) write(iotk_output_unit,"(a)") &
"(use of raw data files is discouraged)."
if(printme) write(iotk_output_unit,"(a)") &
"If file is not present, unit is assumed to be already connected."
if(printme) write(iotk_output_unit,"(a)") &
"If root is present, the name of root in file is read onto that variable."
if(printme) write(iotk_output_unit,"(a)") &
"If attr is present, the attributes of the root tag are read onto that variable,"
if(printme) write(iotk_output_unit,"(a)") &
"which can be subsequently decoded with iotk_scan_attr."
if(printme) write(iotk_output_unit,"(a)") &
"DON'T TRY TO MANIPULATE THE ATTRIBUTE STRING DIRECTLY!"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printlist) write(iotk_output_unit,"(a)") &
" closing_files iotk_close_write iotk_close_read"
printme=.false.
if(iotk_strcomp(keyword,"all")) printme=.true.
if(iotk_strcomp(keyword,'closing_files')) printme=.true.
if(iotk_strcomp(keyword,'iotk_close_write')) printme=.true.
if(iotk_strcomp(keyword,'iotk_close_read')) printme=.true.
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"IOTK: CLOSING FILES"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"iotk_close_write(unit[,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"iotk_close_read(unit[,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"integer,      intent(in)  :: unit"
if(printme) write(iotk_output_unit,"(a)") &
"integer,      intent(out) :: ierr ! see error_handling page"
if(printme) write(iotk_output_unit,"(a)") &
"This routines close a file opened with iotk_open_*"
if(printme) write(iotk_output_unit,"(a)") &
"Note that if the units were already connected before iotk_open_*, they"
if(printme) write(iotk_output_unit,"(a)") &
"are left connected here."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printlist) write(iotk_output_unit,"(a)") &
" multiple_files iotk_link"
printme=.false.
if(iotk_strcomp(keyword,"all")) printme=.true.
if(iotk_strcomp(keyword,'multiple_files')) printme=.true.
if(iotk_strcomp(keyword,'iotk_link')) printme=.true.
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"IOTK: MULTIPLE FILES"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"When reading, if a begin tag with an attribute iotk_link="//'"'//&
"FILENAME"//'"'//&
" is found,"
if(printme) write(iotk_output_unit,"(a)") &
"file FILENAME is mounted in its place"
if(printme) write(iotk_output_unit,"(a)") &
"If FILENAME begins with a "//'"'//&
"/"//'"'//&
", the path is absolute, otherwise it is relative"
if(printme) write(iotk_output_unit,"(a)") &
"to the original file."
if(printme) write(iotk_output_unit,"(a)") &
"Note that the mounting is completely transparent for users, which can access"
if(printme) write(iotk_output_unit,"(a)") &
"the new file using the old unit. However, if the user wants to access"
if(printme) write(iotk_output_unit,"(a)") &
"directly the new file, iotk_physical_unit should be used."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"When writing, the user can switch a logical unit to a different file using"
if(printme) write(iotk_output_unit,"(a)") &
"the following routine"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"iotk_link(unit,name,file,dummy[,binary][,raw][,create][,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(in)  :: unit"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in)  :: name"
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*), intent(in)  :: file"
if(printme) write(iotk_output_unit,"(a)") &
"logical,          intent(in)  :: binary"
if(printme) write(iotk_output_unit,"(a)") &
"logical,          intent(in)  :: raw"
if(printme) write(iotk_output_unit,"(a)") &
"logical,          intent(in)  :: create"
if(printme) write(iotk_output_unit,"(a)") &
"integer,          intent(out) :: ierr"
if(printme) write(iotk_output_unit,"(a)") &
"name is the name of the tag which represents the link."
if(printme) write(iotk_output_unit,"(a)") &
"file is the name of the new file"
if(printme) write(iotk_output_unit,"(a)") &
"if binary is present and true, the new file will be binary"
if(printme) write(iotk_output_unit,"(a)") &
"if raw is present and true, the new file will be raw"
if(printme) write(iotk_output_unit,"(a)") &
"if create is present and true, the new file is actually created"
if(printme) write(iotk_output_unit,"(a)") &
"and the next write statement will act on this new file automatically."
if(printme) write(iotk_output_unit,"(a)") &
"Otherwise, only the symbolic link is created."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printlist) write(iotk_output_unit,"(a)") &
" utilities"
printme=.false.
if(iotk_strcomp(keyword,"all")) printme=.true.
if(iotk_strcomp(keyword,'utilities')) printme=.true.
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"IOTK: OTHER UTILITIES"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"Here a number of additional routines/parameters available"
if(printme) write(iotk_output_unit,"(a)") &
"from the iotk_module is listed"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*) :: iotk_index (index)"
if(printme) write(iotk_output_unit,"(a)") &
"integer, intent(in) :: index ! scalar or rank 1"
if(printme) write(iotk_output_unit,"(a)") &
"Returns a string representing the index in an array."
if(printme) write(iotk_output_unit,"(a)") &
"Example: index = (/1,2,3/) => iotk_index = "//'"'//&
".1.2.3"//'"'//&
""
if(printme) write(iotk_output_unit,"(a)") &
"The correct way for writing an array of derived types is"
if(printme) write(iotk_output_unit,"(a)") &
"to build the names as follows"
if(printme) write(iotk_output_unit,"(a)") &
"! ONE-DIMENSIONAL ARRAY"
if(printme) write(iotk_output_unit,"(a)") &
"do i = 1 , n"
if(printme) write(iotk_output_unit,"(a)") &
"call iotk_write_begin(unit,"//'"'//&
"dummy"//'"'//&
"//iotk_index(i))"
if(printme) write(iotk_output_unit,"(a)") &
"! WRITE THE OBJECT HERE"
if(printme) write(iotk_output_unit,"(a)") &
"call iotk_write_end  (unit,"//'"'//&
"dummy"//'"'//&
"//iotk_index(i))"
if(printme) write(iotk_output_unit,"(a)") &
"end do"
if(printme) write(iotk_output_unit,"(a)") &
"do i = 1 , n"
if(printme) write(iotk_output_unit,"(a)") &
"do j = 1 , m"
if(printme) write(iotk_output_unit,"(a)") &
"! NOTE THE ORDER OF INDEXES, THE FASTER IS THE LAST"
if(printme) write(iotk_output_unit,"(a)") &
"call iotk_write_begin(unit,"//'"'//&
"dummy"//'"'//&
"//iotk_index((/i,j/)))"
if(printme) write(iotk_output_unit,"(a)") &
"! WRITE THE OBJECT HERE"
if(printme) write(iotk_output_unit,"(a)") &
"call iotk_write_end  (unit,"//'"'//&
"dummy"//'"'//&
"//iotk_index((/i,j/)))"
if(printme) write(iotk_output_unit,"(a)") &
"end do"
if(printme) write(iotk_output_unit,"(a)") &
"end do"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"iotk_free_unit(unit[,ierr])"
if(printme) write(iotk_output_unit,"(a)") &
"integer, intent(out) :: unit"
if(printme) write(iotk_output_unit,"(a)") &
"integer, intent(out) :: ierr"
if(printme) write(iotk_output_unit,"(a)") &
"This routine returns the number of a free FORTRAN unit."
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"character(len=*) :: iotk_version"
if(printme) write(iotk_output_unit,"(a)") &
"version string of iotk"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"character :: iotk_newline"
if(printme) write(iotk_output_unit,"(a)") &
"newline sequence"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"character :: iotk_eos"
if(printme) write(iotk_output_unit,"(a)") &
"end-of-string character"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"integer :: iotk_taglenx"
if(printme) write(iotk_output_unit,"(a)") &
"max length of a tag"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"integer :: iotk_namlenx"
if(printme) write(iotk_output_unit,"(a)") &
"max length of a tag or attribute name"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"integer :: iotk_attlenx"
if(printme) write(iotk_output_unit,"(a)") &
"max length of the attribute string"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"integer :: iotk_vallenx"
if(printme) write(iotk_output_unit,"(a)") &
"max length of the value of an attribute"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"integer :: iotk_linlenx"
if(printme) write(iotk_output_unit,"(a)") &
"max length of a line in textual files"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"integer :: iotk_fillenx"
if(printme) write(iotk_output_unit,"(a)") &
"max length of a file name"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
"integer :: iotk_header_kind"
if(printme) write(iotk_output_unit,"(a)") &
"integer kind of headers in binary files"
if(printme) write(iotk_output_unit,"(a)") &
""
if(printme) write(iotk_output_unit,"(a)") &
""

  1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_tool_man_x


