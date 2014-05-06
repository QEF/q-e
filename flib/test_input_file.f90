!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
LOGICAL FUNCTION test_input_xml (myunit)
   !
   ! check if file opened as unit "myunit" is a xml file or not
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(in) :: myunit
   !
   CHARACTER(LEN=256) :: dummy
   CHARACTER(LEN=1), EXTERNAL :: capital
   INTEGER :: i, j 
   LOGICAL :: exst
   !
   test_input_xml = .false.
   INQUIRE ( UNIT=myunit, EXIST=exst )
   IF ( .NOT. exst ) GO TO 10
   
   ! read until a non-empty line is found

   dummy = ' '
   DO WHILE ( LEN_TRIM(dummy) < 1 )
      READ ( myunit,'(A)', ERR=10, END=10) dummy
   END DO

   ! remove blanks from line, convert to capital, clean trailing characters

   j=1
   DO i=1, LEN_TRIM(dummy) 
      IF ( dummy(i:i) /= ' ' ) THEN
         dummy(j:j) = capital(dummy(i:i))
         j=j+1
      END IF
   END DO
   DO i=j, LEN_TRIM(dummy) 
      dummy(i:i) = ' '
   END DO

   ! check for string "<?xml" or "<xml" in the beginning, ">" at the end

   j = LEN_TRIM (dummy)
   test_input_xml = ( (dummy(1:5) == "<?XML") .OR. (dummy(1:4) == "<XML") ) &
                      .AND. (dummy(j:j) == ">")
   RETURN

10 WRITE (0,"('from test_input_xml: input file not opened or empty')")

END FUNCTION test_input_xml
