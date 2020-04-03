program create_pw_input
implicit none
!
character(len=256)::prefix,pseudo_dir,out_dir,pseudodir,&
  file_datipos,file_dativel,outdir,thermodir,current_output
character(len=256),allocatable :: name(:),pseudo(:)
real*16, allocatable ::masses(:),posizioni(:,:,:),velocita(:,:,:)
real*16 ::ecutwfc,cell,dt,eta
integer, allocatable ::natoms(:)
integer ::nsp,isp,tot_at,iat,n_max

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

pseudo_dir=trim(pseudodir)

!
!carico le posizioni dal file
!
allocate(posizioni(tot_at,nsp,3))
open(unit=21,file=trim(file_datipos),status='unknown')
read(21,*)
do isp=1,nsp
   do iat=1,natoms(isp)
      read(21,*) posizioni(iat,isp,1:3)
   end do
end do
close(21)

!
!carico le velocità dal file
!
allocate(velocita(tot_at,nsp,3))
open(unit=21,file=trim(file_dativel),status='unknown')
read(21,*)
do isp=1,nsp
   do iat=1,natoms(isp)
      read(21,*) velocita(iat,isp,1:3)
   end do
end do
close(21)


!!!! FIRST STEP
!
!creazione del file PW_FIRST
!
open(unit=20,file='input_pw_1',status='unknown')

write(20,'(A)')     " &CONTROL"
write(20,'(A)')   "    calculation='scf',"
write(20,'(A)')     "    restart_mode='from_scratch',"  
write(20,'(A,A,A)') "    pseudo_dir='", trim(pseudo_dir), "',"
write(20,'(A,A,A)') "    outdir='", trim(outdir), "',"
write(20,'(A)') "    prefix= 'a_blocco',"
write(20,'(A)')     "    tprnfor=.true.,"
write(20,'(A)')     '/'

write(20,'(A)')      " &SYSTEM"
write(20,'(A)')      "    ibrav=1,"
write(20,'(A,f20.15,A)') "    celldm(1)= ", cell, ","
write(20,'(A,I5,A)')   "    nat= ", tot_at, ","
write(20,'(A,I5,A)')   "    ntyp= ", nsp, ","
write(20,'(A,f7.3,A)') "    ecutwfc= ", ecutwfc, ","
write(20,'(A)')      '/'

write(20,'(A)') " &ELECTRONS"
write(20,'(A)') "    conv_thr = 1.D-10,"
write(20,'(A)') "    mixing_beta = 0.7,"
write(20,'(A)') '/'

write(20,'(A)') " ATOMIC_SPECIES"
do isp=1,nsp
   write(20,'(3x,A,1x,f15.8,1x,A)') trim(name(isp)),masses(isp),trim(pseudo(isp)) 
end do

write(20,'(A)') " ATOMIC_POSITIONS {bohr}"
do isp=1,nsp 
   do iat=1,natoms(isp)
      write(20,'(A,3F20.14)') trim(name(isp)),posizioni(iat,isp,1),posizioni(iat,isp,2),posizioni(iat,isp,3)
   end do
end do

write(20,'(A)') " K_POINTS {Gamma}"
close(20)
 

!!!! SECOND STEP
!
!aggiornamento posizioni (notare il fattore due dovuto al cambio di unità tra CP e PW)
!
do isp=1,nsp
   do iat=1,natoms(isp)     
      posizioni(iat,isp,1:3) = posizioni(iat,isp,1:3) - 2.d0 * velocita(iat,isp,1:3) * dt
   end do
end do

!
!creazione del file PW_SECOND
!
open(unit=20,file='input_pw_2',status='unknown')

write(20,'(A)')     " &CONTROL"
write(20,'(A)')   "    calculation='scf',"
write(20,'(A)')     "    restart_mode='from_scratch',"  
write(20,'(A,A,A)') "    pseudo_dir='", trim(pseudo_dir), "',"
write(20,'(A,A,A)') "    outdir='", trim(outdir), "',"
write(20,'(A)') "    prefix= 'b_blocco',"
write(20,'(A)')     "    tprnfor=.true.,"
write(20,'(A)')     '/'

write(20,'(A)')      " &SYSTEM"
write(20,'(A)')      "    ibrav=1,"
write(20,'(A,f20.15,A)') "    celldm(1)= ", cell, ","
write(20,'(A,I5,A)')   "    nat= ", tot_at, ","
write(20,'(A,I5,A)')   "    ntyp= ", nsp, ","
write(20,'(A,f7.3,A)') "    ecutwfc= ", ecutwfc, ","
write(20,'(A)')      '/'

write(20,'(A)') " &ELECTRONS"
write(20,'(A)') "    conv_thr = 1.D-10,"
write(20,'(A)') "    mixing_beta = 0.7,"
write(20,'(A)') '/'

write(20,'(A)') " ATOMIC_SPECIES"
do isp=1,nsp
   write(20,'(3x,A,1x,f15.8,1x,A)') trim(name(isp)),masses(isp),trim(pseudo(isp)) 
end do

write(20,'(A)') " ATOMIC_POSITIONS {bohr}"
do isp=1,nsp 
   do iat=1,natoms(isp)
      write(20,'(A,3F20.14)') trim(name(isp)),posizioni(iat,isp,1),posizioni(iat,isp,2),posizioni(iat,isp,3)
   end do
end do

write(20,'(A)') " K_POINTS {Gamma}"
close(20)

end program create_pw_input
