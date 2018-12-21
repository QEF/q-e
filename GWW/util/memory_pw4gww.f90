program memory_pw4gww

  implicit none
  
  integer :: numpw,numt,nsetpola,nsteps,nproc
  real :: mem
  
  write(*,*) 'NUMPW:'
  read(*,*) numpw
  write(*,*) 'NUMT:'
  read(*,*) numt
  write(*,*) 'N_SET_POLA:'
  read(*,*) nsetpola
  write(*,*) 'NSTEPS_LANCZOS_POLA:'
  read(*,*) nsteps
  write(*,*) 'N MPI TASKS:'
  read(*,*) nproc

  mem=4.*numpw*numt*nsetpola
  mem=mem+8.*numpw*numt+8.*numpw*nsteps*nsetpola*numt/nproc
  mem=mem+8.*numpw*numpw
  mem=mem+4.*numt*numt/nproc
  mem=mem+8*numpw*numt/nproc
  mem=mem+8.*2.*numpw*numpw

  mem=mem/(1024)**2.d0

  write(*,*) 'Required memory(MB): ', mem

stop
end program memory_pw4gww
