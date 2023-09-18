program test
  implicit none
  complex*16 :: a(2), b(2)
  complex*16 res
  complex*16, external :: ZDOTC

  a(1) = (-9.6224246089610388d-2, 3.2442340359108593d-3)
  a(2) = (-0.9037769140058165, 3.2441868631152768d-3)

  b(1) = (-0.9999890875238262, -5.3357405582201908d-7)
  b(2) = (-0.9999998761616069, 4.3341267956954060d-8)

  res = ZDOTC(2, a, 1, b, 1)
  if (ABS(res - (1.0000000308637980,6.48839729724212839E-003)) >= 1.d-6) then
     write(*,*) "zdotc check failed. Expected (1.0000000308637980,6.48839729724212839E-003) but got ", res
     stop 1
  endif
end program
