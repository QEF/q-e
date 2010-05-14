#!/usr/bin/env python

# (C) 2010 Norbert Nemec
#
# USAGE: src-normal.py < input.f90 > output.f90
#
# Script to normalize Fortran source code:
#   a) expand tabs to spaces (tab width 8 characters
#   b) remove trailing space
#   c) normalize multiword keywords
#   d) normalize capitalization of keywords and intrinsics
#   d) replace old relational operators (.eq., .gt., etc.) by new ones (==, >, etc.)
# The script skips comments and strings within the code

import sys,re

dropspace_list = [
 "BLOCK *DATA",
 "CASE *DEFAULT",     # SPLIT NOT OPTIONAL !
 "DOUBLE *PRECISION",
 "DO *WHILE",         # SPLIT NOT OPTIONAL !
 "ELSE *IF",
 "END *BLOCK *DATA",
 "END *DO",
 "END *FILE",
 "END *FORALL",
 "END *FUNCTION",
 "END *IF",
 "END *INTERFACE",
 "END *MODULE",
 "END *PROGRAM",
 "END *SELECT",
 "END *SUBROUTINE",
 "END *TYPE",
 "END *WHERE",
 "GO *TO",
 "IN *OUT",
 "MODULE *PROCEDURE", # SPLIT NOT OPTIONAL !
 "SELECT *CASE",
]

splitword_list = [
 "BLOCK DATA",
 "CASE DEFAULT", # SPLIT NOT OPTIONAL
 "DOUBLE PRECISION",
 "DO WHILE", # SPLIT NOT OPTIONAL
# "ELSEIF", # leave as one word
 "END BLOCK DATA",
# "ENDDO",  # leave as one word
 "END FILE",
 "END FORALL",
 "END FUNCTION",
# "ENDIF",  # leave as one word
 "END INTERFACE",
 "END MODULE",
 "END PROGRAM",
 "END SELECT",
 "END SUBROUTINE",
 "END TYPE",
 "END WHERE",
# "GOTO",   # leave as one word
# "INOUT",  # leave as one word
 "MODULE PROCEDURE", # SPLIT NOT OPTIONAL
 "SELECT CASE",
]


dropspace_re = re.compile(r"\b("+"|".join(dropspace_list)+r")\b",re.I)

def dropspace_fn(s):
    return s.group(0).replace(" ","")

splitword_dict = dict( (a.replace(" ","").lower(),a) for a in splitword_list )
splitword_re = re.compile(r"\b("+"|".join(splitword_list).replace(" ","")+r")\b",re.I)

def splitword_fn(s):
    return splitword_dict[s.group(0).lower()]


uppercase_keywords = r"""
MODULE SUBROUTINE PROGRAM FUNCTION INTERFACE
ENDMODULE ENDSUBROUTINE ENDPROGRAM ENDFUNCTION ENDINTERFACE
BLOCKDATA DOUBLEPRECISION
MODULEPROCEDURE
TYPE ENDTYPE
CONTAINS
USE ONLY
ALLOCATABLE DIMENSION INTENT EXTERNAL INTRINSIC OPTIONAL PARAMETER POINTER
COMMON
FORMAT
IMPLICIT NONE
PRIVATE PUBLIC
CHARACTER COMPLEX INTEGER LOGICAL
ENTRY EQUIVALENCE INCLUDE NAMELIST SAVE SEQUENCE TARGET
ELEMENTAL PURE RECURSIVE RESULT

SELECTCASE CASE CASEDEFAULT ENDSELECT
IF THEN ELSEIF ELSE ENDIF
WHERE ELSEWHERE ENDWHERE
FORALL ENDFORALL
DO DOWHILE ENDDO

ALLOCATE ASSIGN BACKSPACE CALL CLOSE CONTINUE CYCLE DEALLOCATE ENDFILE
EXIT FORMAT GOTO INQUIRE NULLIFY OPEN PAUSE PRINT READ RETURN REWIND STOP WRITE

""".split()

lowercase_keywords = r"""
in inout out
""".split()

intrinsics = r"""
abort abs achar acos acosd acosh adjustl adjustr aimag aint all allocated and anint any asin
asind asinh associated atan atan2 atan2d atand atanh
baddress bit_size btest
ceiling char cmplx conjg cos cosd cosh count cshift
date date_and_time dble dcmplx dfloat digits dim dnum dot_product dprod dreal
eoshift epsilon exit exp exponent
floor flush fnum fraction free fset fstream
getarg getenv gran
hfix huge
iachar iaddr iand iargc ibclr ibits ibset ichar idate idim ieor igetarg ijint imag index int int1
int2 int4 int8 inum iomsg ior iqint irand iranp ishft ishftc isign ixor izext
jnum jzext
kind kzext
lbound len len_trim lge lgt lle llt loc log log10 lshft lshift
malloc matmul max maxexponent maxloc maxval mclock merge min minexponent minloc minval mod modulo mvbits
nearest nint not
or
pack precision present product
qext qfloat qnum qprod
radix ran rand random_number random_seed range repeat reshape rnum rrspacing rshft rshift
scale scan secnds selected_int_kind selected_real_kind set_exponent shape sign sin sind sinh size
sizeof spacing spread sqrt srand sum system system_clock
tan tand tanh time tiny transfer transpose trim
ubound unpack
verify xor zext
""".split()

ignore_for_the_moment = r"""
real REAL isnan
"""

special_keywords = r"""
.and. .or. .not. .true. .false. .eqv. .neqv.
.eq. .ge. .gt. .le. .lt. .ne.
""".replace(".","\\.").split()

def uppercase_fn(s):
    return s.group(0).upper()

def lowercase_fn(s):
    return s.group(0).lower()

def special_fn(s):
    res = s.group(0).lower()
    res = {
    '.eq.': '==',
    '.ge.': '>=',
    '.gt.': '>',
    '.le.': '<=',
    '.lt.': '<',
    '.ne.': '/=',
    }.get(res,res)
    return res

uppercase_re = re.compile(r"\b("+"|".join(uppercase_keywords)+r")\b",re.I)
lowercase_re = re.compile(r"\b("+"|".join(lowercase_keywords+intrinsics)+r")\b",re.I)
special_re = re.compile(r"("+"|".join(special_keywords)+r")",re.I)

def correctcase(line):
    line = dropspace_re.sub(dropspace_fn,line)
    line = uppercase_re.sub(uppercase_fn,line)
    line = lowercase_re.sub(lowercase_fn,line)
    line = special_re.sub(special_fn,line)
    line = splitword_re.sub(splitword_fn,line)
    return line

##############

quote = " "
QUOTES = "'\""

for lin in sys.stdin:
    lin = lin.rstrip().expandtabs()
    pos = 0
    lout = ""
    if lin[:1] == "#":
        lout=lin
        pos=len(lin)
    while pos < len(lin):
        if quote in QUOTES:
            npos = lin.find(quote,pos)
            if npos >= 0:
                assert lin[npos] == quote
                lout += lin[pos:npos+1]
                pos = npos+1
                quote = " "
            elif lin[-1] == "&":
                lout += lin[pos:]
                break
            else:
                raise "unterminated string in line ["+lin+"]"

        cpos = lin.find("!",pos) % (len(lin)+1)
        qpos = lin.find("'",pos) % (len(lin)+1)
        dpos = lin.find('"',pos) % (len(lin)+1)
        npos = min(cpos,qpos,dpos)
        lout += correctcase(lin[pos:npos])
        pos = npos
        if pos == len(lin):
            break
        elif lin[pos] == "!":
            lout += lin[pos:]
            break
        elif lin[pos] in QUOTES:
            quote = lin[pos]
            lout += quote
            pos += 1
            continue
        else:
            raise "Strange internal error"

    sys.stdout.write(lout+"\n")
