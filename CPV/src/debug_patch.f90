
do i = 1, nbsp
  do j = 1, dffts%nr1*dffts%nr2*dffts%my_nr3p
    write(20000+100*i+me_image,*) exx_potential(j,i)
  end do ! j
  flush(20000+100*i+me_image)
end do ! i
CALL mp_barrier( intra_image_comm )
stop
