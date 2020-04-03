program create_hz_input
implicit none
!
character(len=256) ::prefix,current_output,file_dativel,&
  file_datipos,pseudodir,outdir,thermodir
character(len=256),allocatable :: name(:),pseudo(:)
real*8, allocatable ::masses(:),posizioni(:,:,:)
real*8 ::ecutwfc,cell,dt,eta
integer, allocatable ::natoms(:)
integer ::nsp,isp,tot_at,n_max
!
!lettura del file tabella
!
!open(unit=20,file='tabella',status='old')
read(*,*) nsp
allocate(name(nsp))
allocate(natoms(nsp))
allocate(pseudo(nsp))
allocate(masses(nsp))
do isp=1,nsp 
   read(*,*) name(isp)
   read(*,*) natoms(isp)
   read(*,*) pseudo(isp)
   read(*,*) masses(isp) 
end do
read(*,*) ecutwfc
read(*,*) cell
read(*,*) dt
read(*,*) eta
read(*,*) n_max
read(*,*) file_datipos
read(*,*) file_dativel
read(*,*) pseudodir
read(*,*) outdir
read(*,*) thermodir
read(*,*) current_output
!close(20)
!
!calcolo del numero totale di atomi
!
tot_at=0
do isp=1,nsp
   tot_at=tot_at+natoms(isp)
end do
!
!
!creazione del file
!
open(unit=20,file='input_hz_init',status='unknown')
write(20,'(A)') " &inputhartree"
write(20,'(A)') "    prefix_uno= 'a_blocco',"
write(20,'(A)') "    prefix_due= 'b_blocco',"
write(20,'(A,f7.3,A)') "    delta_t= ",dt,","
write(20,'(A)') "    init_linear= 'niente',"
write(20,'(A,A,A)') "    file_output= '",trim(current_output),"',"
write(20,'(A,A,A)') "    file_dativel= '",trim(file_dativel),"',"
write(20,'(A,A,A)') "    outdir= '",trim(outdir),"',"
write(20,'(A,A,A)') "    thermodir= '",trim(thermodir),"',"
write(20,'(A)') " /"
write(20,'(A)') " &inputzero"
write(20,'(A,f7.3,A)') "    eta= ",eta,","
write(20,'(A,I5,A)') "    n_max= ",n_max,","
write(20,'(A)')  "    status= 'initialize',"
write(20,'(A)') " /"
close(20)
 
!
!creazione del file
!
open(unit=20,file='input_hz_run',status='unknown')
write(20,'(A)') " &inputhartree"
write(20,'(A,A,A)') "    prefix_uno= 'a_blocco',"
write(20,'(A,A,A)') "    prefix_due= 'b_blocco',"
write(20,'(A,f7.3,A)') "    delta_t= ",dt,","
write(20,'(A)') "    init_linear= 'niente',"
write(20,'(A,A,A)') "    file_output= '",trim(current_output),"',"
write(20,'(A,A,A)') "    file_dativel= '",trim(file_dativel),"',"
write(20,'(A,A,A)') "    outdir= '",trim(outdir),"',"
write(20,'(A,A,A)') "    thermodir= '",trim(thermodir),"',"
write(20,'(A)') " /"
write(20,'(A)') " &inputzero"
write(20,'(A,f7.3,A)') "    eta= ",eta,","
write(20,'(A,I5,A)') "    n_max= ",n_max,","
write(20,'(A)')  "    status= 'compute',"
write(20,'(A)') " /"
close(20)

end program create_hz_input
