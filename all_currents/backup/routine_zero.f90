subroutine routine_zero()
use kinds,  only:DP
use wvfct,  only : nbnd,npwx,npw
use io_files,  only :nwordwfc,diropn,iunwfc,prefix, tmp_dir
use fft_base,  only:dffts
use mp,        only :mp_sum,mp_bcast,mp_get
use mp_pools,  only : intra_pool_comm
use io_global, only :stdout,ionode_id,ionode
use ions_base, only :nsp,zv,nat,ityp,amass,tau
use cell_base,  only:tpiba,tpiba2
use uspp_param,only :upf
use atom,      only: rgrid
use mp_world, only:mpime
use fft_base,  only:dffts
use cell_base, only:at,alat,omega
use wavefunctions_module,  only :psic
use fft_interfaces, only :invfft,fwfft
use gvect,  only :ngm,gg,g,nl,nlm,gstart
use gvecs,  only :nls,nlsm
use constants, only : e2, AMU_RY
use uspp, only :nkb
use splines
use zero_mod
use hartree_mod

implicit none


!three dimensional auxiliary real vectors
  real(dp), parameter  :: amconv = AMU_RY
  real(DP) ::n(3),u(3),u_pbc(3),u_x(3),u_y(3),u_z(3),value,x(1:3),ics(3)
  real(DP), external :: qe_erfc
  real(DP) ::modul,erf_value,charge_atom
  real(DP) ::fac_uno,fac_due,fac_tre,fac_quattro
  real(DP), external :: qe_erf
  integer, external  :: find_free_unit


!counters and indexes
integer      ::a,b !cartesian components
integer      ::igm,isp,jgm !reciprocal lattice,types, ions
integer      ::nr1a,nr2a,nr3a
integer      ::im,iqq !pseudo mesh, real_mesh_one_coordinate
integer      ::ii,iv !processors,bands
integer      ::ix,iy,iz !real_mesh_three_coordinates 
integer      ::iun !unit
integer      ::iatom,jatom !atoms
integer      ::n_x,n_y,n_z,n_maxl
integer      ::istep

!variables for calling functions
logical :: exst,l_test,l_scambio,l_scambio_alt

!auxiliary variables 
integer ::err,ir,ieta,enne
real(DP) ::R
character(256) ::filename,pref_box
real(DP), allocatable :: values(:)

write(stdout,*)'INIZIO ROUTINE ZERO'
call start_clock( 'zero_current' )
if (ionode) print*,'ETA',eta
!
l_zero=.true.
!
call allocate_zero()

!
! Loris: commentate queste variabili, nr3p non e' piu' compatibile. 
! Ma nr3_start e nr3_end vengono usati solo in modify_evc, chiamato da 
! add_nl_evc, chiamato solo come test qui sotto
!
!!inizializzazione piani assegnati ad ogni CPU
!nr3s_start=0
!nr3s_end =0
!do ii=1,mpime + 1
!   nr3s_start=nr3s_end+1
!   nr3s_end=nr3s_end+dffts%nr3p(ii)
!end do
!!

select case(status)
    case("initialize") 
!
!inizializzazione delle funzioni da scrivere
!
call init_us_1a()
call init_reciprocal_parts_tab()
!
!scrittura trasformate su file
pref_box=prefix
prefix='thermal'
iun=find_free_unit()
call diropn(iun,'ecur',ngm,exst )
do isp=1,nsp
   call davcio(H_g(:,1,1,isp),ngm,iun,(isp-1)*6+1,1)
   call davcio(H_g(:,2,2,isp),ngm,iun,(isp-1)*6+2,1)
   call davcio(H_g(:,3,3,isp),ngm,iun,(isp-1)*6+3,1)
   call davcio(H_g(:,2,1,isp),ngm,iun,(isp-1)*6+4,1)
   call davcio(H_g(:,3,1,isp),ngm,iun,(isp-1)*6+5,1)
   call davcio(H_g(:,3,2,isp),ngm,iun,(isp-1)*6+6,1)
end do
close(iun)
iun=find_free_unit()
call diropn(iun,'i_uno',ngm,exst )
call davcio(I_uno_g(:,1,1),ngm,iun,1,1)
call davcio(I_uno_g(:,2,2),ngm,iun,2,1)
call davcio(I_uno_g(:,3,3),ngm,iun,3,1)
call davcio(I_uno_g(:,2,1),ngm,iun,4,1)
call davcio(I_uno_g(:,3,1),ngm,iun,5,1)
call davcio(I_uno_g(:,3,2),ngm,iun,6,1)
close(iun)
iun=find_free_unit()
call diropn(iun,'i_due+i_primo',ngm,exst )
call davcio(I_due_g(:),ngm,iun,1,1)
call davcio(I_primo,1,iun,2,1)
close(iun)

prefix=pref_box
!
     case ("compute")
!
if (nkb>0) then
   l_non_loc=.true.
else
   l_non_loc=.false.
end if
!!!!!!! nuovo metodo di lettura
!goto 30
!leggi da file (unformatted) velocità e posizioni degli ioni al tempo istep.Legge solo ionode as usual, check fortran labelling.

!if (ionode) then
!   iun=find_free_unit()
!   call diropn_rect(iun,'ion_positions_unf',3*nat,exst)
!   call davcio (ion_pos(:,:),3*nat,iun,2*istep-1,-1)    
!   call davcio (ion_vel(:,:),3*nat,iun,2*istep  ,-1)  
!   close(iun)    
!end if
!call mp_bcast(ion_pos(:,:),ionode_id)
!call mp_bcast(ion_vel(:,:),ionode_id)
!call check_positions(ion_pos)
!30  do iatom=1,nat 
!   ion_vel(:,iatom)=[0.005,0.008660254,0.0]
!    ion_vel(:,iatom)=[alat*0.003,alat*0.d0,alat*0.d0]
!end do

!lettura velocita'
if (ionode) then
    iun=find_free_unit()
    open(unit=iun,file=trim(tmp_dir)//trim(file_dativel),access='sequential',status='old')
!per raggiungere il passo corretto:
    do istep=1,passo-1
       do iatom=1,nat 
           read(iun,*) 
       end do 
    end do
    do iatom=1,nat
       read(iun,*) ion_vel(1:3,iatom)
    end do
    close(iun)
end if
call mp_bcast(ion_vel(:,:),ionode_id, intra_pool_comm)

!cambio unità di misure da velocità CP a velocità PW
!
!
ion_vel(1:3,1:nat) = 2.d0 * ion_vel(1:3,1:nat)
!
!
!!!!!!!!!!!!!!!!!!!!!! 
!lettura funzione d'onda
close(iunwfc)
call start_clock( 'lett_car' )

if (ionode) print*,'uguali? ',npwx*nbnd,2*nwordwfc
call diropn(iunwfc, 'wfc', 2*nwordwfc, exst )
call davcio (evc_uno,2*nwordwfc,iunwfc,1,-1)
!
!calcolo della carica a partire dalle funzioni d'onda
 charge=0.d0
do iv=1,nbnd,2
   psic=0.d0
   if (iv==nbnd) then
       psic(nls(1:npw))=evc_uno(1:npw,iv)
       psic(nlsm(1:npw))=CONJG(evc_uno(1:npw,iv))
   else
       psic(nls(1:npw))=evc_uno(1:npw,iv)+ ( 0.D0, 1.D0 ) *evc_uno(1:npw,iv+1)
       psic(nlsm(1:npw))=CONJG(evc_uno(1:npw,iv)- ( 0.D0, 1.D0 ) *evc_uno(1:npw,iv+1))
   end if
   call invfft ('Wave', psic, dffts)  
   charge(1:dffts%nnr)=charge(1:dffts%nnr)+dble(psic(1:dffts%nnr))**2.0 
   if (iv /=nbnd) then
      charge(1:dffts%nnr)=charge(1:dffts%nnr)+dimag(psic(1:dffts%nnr))**2.0
   end if
end do
!
!moltiplico per due causa degenerazione di spin
 charge(1:dffts%nnr)=charge(1:dffts%nnr)*2.d0
!carica in spazio reciproco
psic=0.d0
psic(1:dffts%nnr)=dcmplx(charge(1:dffts%nnr),0.d0)
call fwfft ('Smooth', psic, dffts) 
 charge_g(1:ngm)=psic(nls(1:ngm))
call stop_clock( 'lett_car' )

!lettura di H_g e simmetrizzazione
pref_box=prefix
prefix='thermal'
!
call start_clock( 'lett_H' )
iun=find_free_unit()
call diropn(iun,'ecur',ngm,exst )
do isp=1,nsp 
   call davcio (H_g(:,1,1,isp),ngm,iun,(isp-1)*6+1,-1)    
   call davcio (H_g(:,2,2,isp),ngm,iun,(isp-1)*6+2,-1)    
   call davcio (H_g(:,3,3,isp),ngm,iun,(isp-1)*6+3,-1)    
   call davcio (H_g(:,2,1,isp),ngm,iun,(isp-1)*6+4,-1)    
   call davcio (H_g(:,3,1,isp),ngm,iun,(isp-1)*6+5,-1)    
   call davcio (H_g(:,3,2,isp),ngm,iun,(isp-1)*6+6,-1)    
end do
close(iun)
!
iun=find_free_unit()
call diropn(iun,'i_uno',ngm,exst )
call davcio (I_uno_g(:,1,1),ngm,iun,1,-1)    
call davcio (I_uno_g(:,2,2),ngm,iun,2,-1)    
call davcio (I_uno_g(:,3,3),ngm,iun,3,-1)    
call davcio (I_uno_g(:,2,1),ngm,iun,4,-1)    
call davcio (I_uno_g(:,3,1),ngm,iun,5,-1)    
call davcio (I_uno_g(:,3,2),ngm,iun,6,-1)    
close(iun)
!
iun=find_free_unit()
call diropn(iun,'i_due+i_primo',ngm,exst )
call davcio (I_due_g(:),ngm,iun,1,-1)
call davcio(I_primo,1,iun,2,-1)    
close(iun)
!
prefix=pref_box
!
!questo è necessario?
do a=1,3 
   do b=1,3
      if (a>b) then
          do isp=1,nsp
             H_g(:,b,a,isp)=H_g(:,a,b,isp)  
          end do  
          I_uno_g(:,b,a)=I_uno_g(:,a,b)
      end if
   end do
end do
call stop_clock( 'lett_H' )
call start_clock( 'init_u' )
!inizializzazione di u_g
u_g=0.d0
do a=1,3
   do b=1,3
       do igm=1,ngm
          do iatom=1,nat
               u_g(igm,a)=u_g(igm,a)-ion_vel(b,iatom)*H_g(igm,a,b,ityp(iatom))*&
&exp(-(0.d0,1.d0)*DOT_PRODUCT(tpiba*g(1:3,igm), alat*tau(1:3,iatom)))
          end do
       end do
   end do
end do 
!
call stop_clock( 'init_u' )
iun=find_free_unit()
!
call start_clock( 'calcolo_z' )
!calcolo della corrente
z_current=0.d0
do a=1,3
   do igm=gstart,ngm
       z_current(a)= z_current(a)+2.d0*dble(charge_g(igm)*conjg(u_g(igm,a)))
   end do
   if (gstart==2) then
       z_current(a)= z_current(a)+dble(charge_g(1)*conjg(u_g(1,a)))
   end if
end do
call mp_sum(z_current, intra_pool_comm)
call stop_clock( 'calcolo_z' )
call start_clock( 'non_locale' )
if (l_non_loc) then
  ! Loris: commentati questi test
  !  l_test=.false.
  !  if (l_test) then
  !      call add_nl_evc(z_current)
  !  else 
        call add_nc_curr(z_current)
  !  end if
end if
call stop_clock( 'non_locale' )


!
!------ IONIC CURRENT ---------------------------------------
!
i_current=0.d0
i_current_a=0.d0
i_current_b=0.d0
i_current_c=0.d0
i_current_d=0.d0
i_current_e=0.d0

call start_clock( 'calcolo_i' )
do iatom=1,nat
!
   i_current(:)=i_current(:)+ion_vel(:,iatom)*(1./2.*amconv*amass(ityp(iatom))*(ion_vel(1,iatom)**2+&
&ion_vel(2,iatom)**2+ion_vel(3,iatom)**2))
   i_current_a(:)=i_current_a(:)+ion_vel(:,iatom)*(1./2.*amconv*amass(ityp(iatom))*(ion_vel(1,iatom)**2+&
&ion_vel(2,iatom)**2+ion_vel(3,iatom)**2))
!
   i_current(:)=i_current(:)+2./3.*e2*zv(ityp(iatom))**2*ion_vel(:,iatom)*I_primo
   i_current_b(:)=i_current_b(:)+2./3.*e2*zv(ityp(iatom))**2*ion_vel(:,iatom)*I_primo

end do

l_scambio=.true.
if (l_scambio) then
    l_scambio_alt=.true.
    if (l_scambio_alt) then
       do iatom=1,nat
         do jatom=1,nat
            if (iatom > jatom) then
               u(1:3) = (tau(:,iatom)-tau(:,jatom))*alat
               call pbc(u(1:3),u_pbc(1:3))
               call I_due_value(value,u_pbc,1)
               i_current(:)=i_current(:)+1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))&
*(ion_vel(:,iatom)+ion_vel(:,jatom))*value
               i_current_c(:)=i_current_c(:)+1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))*&
(ion_vel(:,iatom)+ion_vel(:,jatom))*value

               do a=1,3
                  do b=1,3
                     if (a>b) then
                         call I_uno_value(value,u_pbc,a,b,1)
                          i_current_e(a)=i_current_e(a)-1./2. *e2*zv(ityp(iatom))*zv(ityp(jatom))*&
&(ion_vel(b,jatom)+ion_vel(b,iatom))*value
                          i_current_e(b)=i_current_e(b)-1./2. *e2*zv(ityp(iatom))*zv(ityp(jatom))*&
&(ion_vel(a,jatom)+ion_vel(a,iatom))*value
                          i_current(a)=i_current(a)-1./2. *e2*zv(ityp(iatom))*zv(ityp(jatom))*&
&(ion_vel(b,jatom)+ion_vel(b,iatom))*value
                           i_current(b)=i_current(b)-1./2. *e2*zv(ityp(iatom))*zv(ityp(jatom))*&
&(ion_vel(a,jatom)+ion_vel(a,iatom))*value
                     end if
                     if (a==b) then
                         call I_uno_value(value,u_pbc,a,b,1)
                         i_current_d(a)=i_current_d(a)-1./2. *e2*zv(ityp(iatom))*zv(ityp(jatom))*&
&(ion_vel(b,jatom)+ion_vel(b,iatom))*value
                          i_current(a)=i_current(a)-1./2. *e2*zv(ityp(iatom))*zv(ityp(jatom))*&
&(ion_vel(b,jatom)+ion_vel(b,iatom))*value
                     end if
                  end do
               end do
            end if
         end do
       end do 
    else
       do iatom=1,nat
          do jatom=1,nat
            if (iatom > jatom) then
               u(1:3) = (tau(:,iatom)-tau(:,jatom))*alat
               call pbc(u(1:3),u_pbc(1:3))
               call I_due_value(value,u_pbc,1)
               i_current(:)=i_current(:)+1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))&
*(ion_vel(:,iatom)+ion_vel(:,jatom))*value
               i_current_c(:)=i_current_c(:)+1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))*&
(ion_vel(:,iatom)+ion_vel(:,jatom))*value
               do a=1,3
                  do b=1,3
                     call I_uno_value(value,u_pbc,a,b,1)
                     if (a==b) then
                         i_current_d(a)=i_current_d(a)-1./2. *e2*zv(ityp(iatom))*zv(ityp(jatom))*&
&(ion_vel(b,jatom)+ion_vel(b,iatom))*value
                     else
                         i_current_e(a)=i_current_e(a)-1./2. *e2*zv(ityp(iatom))*zv(ityp(jatom))*&
&(ion_vel(b,jatom)+ion_vel(b,iatom))*value
                     end if
                  i_current(a)=i_current(a)-1./2. *e2*zv(ityp(iatom))*zv(ityp(jatom))*&
&(ion_vel(b,jatom)+ion_vel(b,iatom))*value
                  end do
                end do
            end if
         end do
       end do
    end if
else
   do iatom=1,nat
      do jatom=1,nat
        if (iatom.ne.jatom) then
            u(1:3) = (tau(:,iatom)-tau(:,jatom))*alat
            call pbc(u(1:3),u_pbc(1:3))
            call I_due_value(value,u_pbc,1)
            i_current(:)=i_current(:)+1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))*ion_vel(:,iatom)*value
            i_current_c(:)=i_current_c(:)+1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))*ion_vel(:,iatom)*value          
            do a=1,3  
               do b=1,3            
                  call I_uno_value(value,u_pbc,a,b,1)
                  if (a==b) then
                      i_current_d(a)=i_current_d(a)-1./2. *e2*zv(ityp(iatom))*zv(ityp(jatom))*&
&ion_vel(b,jatom)*value
                  else
                     i_current_e(a)=i_current_e(a)-1./2. *e2*zv(ityp(iatom))*zv(ityp(jatom))*&
&ion_vel(b,jatom)*value
                  end if
               i_current(a)=i_current(a)-1./2. *e2*zv(ityp(iatom))*zv(ityp(jatom))*&
&ion_vel(b,jatom)*value
               end do
             end do 
         end if
      end do
   end do
end if
call stop_clock( 'calcolo_i' )


!!!!!!!commentato temporaneamente
if (ionode) then   
    iun=find_free_unit()
    open (iun, file = trim(tmp_dir)//trim(file_output),position='append')
    write(iun,'(A,3E20.12)') 'ionic:', i_current(:)
    write(iun,'(A,3E20.12)') 'ionic_a:', i_current_a(:)
    write(iun,'(A,3E20.12)') 'ionic_b:', i_current_b(:)
    write(iun,'(A,3E20.12)') 'ionic_c:', i_current_c(:)
    write(iun,'(A,3E20.12)') 'ionic_d:', i_current_d(:)
    write(iun,'(A,3E20.12)') 'ionic_e:', i_current_e(:)
    write(iun,'(A,3E20.12)') 'zero:',  z_current(:)
    close(iun)
end if

!if (ionode) then
!    iun=find_free_unit()
!    open (iun, file = trim(tmp_dir)//trim(file_output),position='append')
!    write(iun,'(3f17.13)')  i_current_b(:)
!    write(iun,'(3f17.13)')  i_current_c(:)
!    write(iun,'(3f17.13)')  i_current_d(:)
!    write(iun,'(3f17.13)')  i_current_e(:)
!    write(iun,'(3f17.13)')  z_current(:)
!    close(iun)
!end if




! 
     case default
! 
 write(stdout,*) "Unknown keyword for zero_current"


end select
!
call deallocate_zero()
!
300 call stop_clock( 'zero_current' )
call print_clock( 'zero_current' )
call print_clock( 'real' )
call print_clock( 'lett_car' )
call print_clock( 'lett_H' )
call print_clock( 'init_u' )
call print_clock( 'calcolo_z' )
call print_clock( 'non_locale' )
call print_clock( 'calcolo_i' )
call print_clock( 'realrec' )


end subroutine routine_zero


