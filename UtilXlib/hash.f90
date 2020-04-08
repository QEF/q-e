program hash
    use fletcher32_mod
    implicit none
    integer, parameter :: dat_size=6
    integer :: i, cksum
    integer(2) :: dat(1)
    complex(8) :: mydata(dat_size)
    real(8) :: myreal(dat_size*2)

    do i =1,dat_size
       mydata(i) = cmplx(33.0*i-22.0, 24.9 + i)
       myreal(2*i-1) = dreal(mydata(i))
       myreal(2*i)   = dimag(mydata(i))
    end do
    cksum = fletcher32(transfer(mydata,dat),size(transfer(mydata,dat)))
    write (*,*) 'checksum of a complex array directly casting into the c function ', cksum
    call fletcher32_cksum(mydata,cksum)
    write (*,*) 'checksum of a complex array using fortran interface ', cksum
    call fletcher32_cksum(myreal,cksum)
    write (*,*) 'checksum of the same data seen as a real array using fortran interface ', cksum
    do i =1,dat_size
       myreal(2*i-1) = dimag(mydata(i))
       myreal(2*i)   = dreal(mydata(i))
    end do
    call fletcher32_cksum(myreal,cksum)
    write (*,*) 'checksum of the same data in different order using fortran interface ', cksum

end program
