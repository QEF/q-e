
subroutine add_nc_curr(current)
!
use kinds, only :DP
use becmod
use us, only: spline_ps
use uspp,       ONLY : nkb, vkb, deeq
USE uspp_param, ONLY : upf,nh
use zero_mod, only: becpr,becpd,becprd,xvkb,xdvkb,dvkb, ion_vel,evc_uno
use wvfct,                ONLY : nbnd, npw, npwx
use gvect,                ONLY : g
use cell_base,            ONLY : tpiba
use ions_base,  ONLY : nat, ityp, nsp
use klist, only:xk
use io_files, ONLY : nwordwfc,diropn
use io_global, only : ionode
use wavefunctions_module, only:evc
use klist, only : igk_k
!
implicit none 
!
integer ::iun,a,b
logical ::l_test,exst
real(DP), intent(inout) :: current(3)
integer :: ipol,jpol,ipwd,ijk,ibnd
integer ::ijkb,ikb,ih,na,nt,ipw
real(DP) ::J_nl(3),J_1(3),J_2(3)
complex(DP),allocatable ::write(:,:)
complex(DP),allocatable ::vkb1(:,:)
integer, external       :: find_free_unit
!
allocate (vkb1(npwx,nkb))
allocate(write(npwx,nbnd))
!
!non fare nulla se non vi è potenziale non locale.
if (nkb <=0) return
!check se vi sono pseudo ultrasoffici o PAW ed in tal caso blocca il calcolo.
do nt=1,nsp
   if ((upf(nt)%typ.eq."US") .or. (upf(nt)%typ.eq."PAW")) then
      CALL errore ('add_nc_curr', 'US and PAW not implemented', 1) 
    end if
end do
!
!memory initialization
CALL allocate_bec_type( nkb, nbnd, becp ) 
do ipol=1,3
   CALL allocate_bec_type( nkb, nbnd,becpr(ipol) )    
   CALL allocate_bec_type( nkb, nbnd,becpd(ipol) ) 
   do jpol=1,3
       CALL allocate_bec_type(  nkb, nbnd,becprd(ipol,jpol)) 
   end do
end do
!
!allocazione delle vkb (vkb già allocato?)
allocate(xvkb(npwx,nkb,3))
allocate(dvkb(npwx,nkb,3))
allocate(xdvkb(npwx,nkb,3,3))

!
!inizializzazione di tab(serve?),tabr ed indici
call init_us_1a()
!inizializzazione di vkb (per essere sicuri che lo sia) e xvkb
CALL init_us_2( npw, igk_k(1,1), xk(1,1), vkb )
call init_us_3(npw,xvkb)
!
!!test
!l_test=.false.
!if (l_test) then
!    iun=find_free_unit()
!    do ibnd=1,nbnd
!       write(1:npwx,ibnd)=vkb(1:npwx,1)
!    end do
!    call diropn(iun,'onda',2*nwordwfc,exst )
!    call davcio(write,2*nwordwfc,iun,+1,+1)
!    close(iun)
!end if
!!
!
!inizializzazione di dvkb e xdvkb (nb il ciclo su ipw va dentro per essere
!ottimizzato)
!servono altre inizializzazioni?
dvkb=0.d0
do ipol=1,3
   do ikb=1,nkb
      do ipw=1,npw
         dvkb(ipw,ikb,ipol)=(0.d0,-1.d0)*tpiba*g(ipol,ipw)*vkb(ipw,ikb)
!         dvkb(ipw,ikb,ipol)=(0.d0,-1.d0)*g(ipol,ipw)*vkb(ipw,ikb)
      end do
   end do
end do
do jpol=1,3
   do ipol=1,3
      do ikb=1,nkb
         do ipw=1,npw

            xdvkb(ipw,ikb,ipol,jpol)=(0.d0,-1.d0)*tpiba*g(jpol,ipw)*xvkb(ipw,ikb,ipol)
!            xdvkb(ipw,ikb,ipol,jpol)=(0.d0,-1.d0)*g(jpol,ipw)*xvkb(ipw,ikb,ipol)
            if (ipol==jpol) then
                 xdvkb(ipw,ikb,ipol,jpol)=xdvkb(ipw,ikb,ipol,jpol)+vkb(ipw,ikb)
                 
            end if
         end do
      end do
   end do
end do
!
!
!prodotti scalari (si possono evitare di fare tutti?), qui si esegue comunicazione MPI.
call calbec(npw,vkb,evc,becp)
do ipol=1,3
   CALL calbec(npw,xvkb(1:npwx,1:nkb,ipol),evc_uno,becpr(ipol))  
   do ipw=1,npw
      do ikb=1,nkb
         vkb1(ipw,ikb)=dvkb(ipw,ikb,ipol)
      end do
   end do 
   CALL calbec(npw,vkb1,evc,becpd(ipol))
   do jpol=1,3
       call calbec(npw,xdvkb(1:npwx,1:nkb,ipol,jpol),evc_uno,becprd(ipol,jpol))
   end do
end do
!!------------------------------------------------------------
!      l_test=.false.
!      if (l_test) then   
!          if (ionode) then
!              open(unit=16,file='becprd(nc)->rdbecp_mod(nl)',status='unknown')
!              open(unit=17,file='becpr(nc)->becp_mod(nl)',status='unknown')
!              open(unit=18,file='becpd(nc)->rdbecp(nl)',status='unknown')
!              open(unit=19,file='becp(nc)->becp(nl)',status='unknown')
!              open(unit=20,file='vkb1(nc)',status='unknown')
!              open(unit=21,file='vkb(nc)',status='unknown')
!              open(unit=22,file='evc(nc)',status='unknown')
!!          
!              do a=1,3
!                 do ibnd=1,nbnd
!                    do ikb=1,nkb
!                       do b=1,3
!                          write(16,"('BECPRD: A-B-IKB-IBND..value',4I4,F12.7)" ) &
!a,b,ikb,ibnd,becprd(a,b)%r(ikb,ibnd)
!                       end do  
!                    end do
!                 end do
!              end do
!!          
!!confronto con becpr 
!              do a=1,3
!                 do ibnd=1,nbnd
!                    do ikb=1,nkb
!                       write(17, "('BECPR: A-IKB-IBND..value',3I4,F12.7)")&
!a,ikb,ibnd,becpr(a)%r(ikb,ibnd)
!                    end do
!                 end do
!              end do
!!confrondo per becpd
!              do ibnd=1,nbnd
!                 do ikb=1,nkb
!                    do ipol=1,3
!                        write(18, "('BECPR: A-IKB-IBND..value',3I4,F12.7)")&
!ipol,ikb,ibnd,becpd(ipol)%r(ikb,ibnd)
!                    end do  
!                 end do
!              end do
!!confronto per becp 
!              do ibnd=1,nbnd
!                 do ikb=1,nkb
!                     write(19, "('BECPR: IKB-IBND..value',2I4,F12.7)")&
!ikb,ibnd,becp%r(ikb,ibnd) 
!                 end do
!              end do
!              do ipol=1,3
!                 write(20,*)'IPOL: ',ipol
!                 do ikb=1,nkb
!                    write(20,*) 'JKB: ',ikb
!                    do ipw=1,npw
!                       write(20,*) ipw,dvkb(ipw,ikb,ipol)
!                    end do
!                 end do
!              end do 
!              do ikb=1,nkb
!                 write(21,*) 'JKB: ',ikb
!                 do ipw=1,npw
!                    write(21,*) ipw,vkb(ipw,ikb)
!                 end do
!              end do
!              do ibnd=1,nbnd
!                 write(22,*) 'IBND: ',ibnd
!                 do ipw=1,npw
!                    write(22,*) ipw,evc(ipw,ibnd)
!                 end do
!              end do
!
!              close(16)
!              close(17)
!              close(18)
!              close(19)
!              close(20)
!              close(21)
!              close(22) 
!          end if
!      end if 
!!---------------------------------------------------------
!!--------------------------------------------------------
!l_test=.false.
!if (l_test) then
!   open(unit=15,file='becp',status='unknown')
!   do ikb=1,nkb
!      do ibnd=1,nbnd
!         write(15,*) 'RDxP',ikb,ibnd,becprd(2,1)%r(ikb,ibnd)*becp%r(ikb,ibnd),&
!becprd(2,1)%r(ikb,ibnd),becp%r(ikb,ibnd) 
!         write(15,*) 'RxD',ikb,ibnd,becpr(2)%r(ikb,ibnd)*becpd(1)%r(ikb,ibnd),&
!becpr(2)%r(ikb,ibnd),becpd(1)%r(ikb,ibnd)
!      end do
!   end do
!   close(15)
!end if
!!--------------------------------------------------------------------

!
!a questo punto possiamo usare le quantità calcolate per calcolare la corrente (OpenMP do?)
J_nl=0.d0
J_1=0.d0
J_2=0.d0
ijkb=0
do nt=1,nsp
   do na=1,nat
        if (ityp (na) .eq. nt) then
            do ih = 1, nh (nt)
               ikb=ijkb+ih 
               do ipol=1,3
                  do jpol=1,3
                     do ibnd=1,nbnd
                        J_nl(ipol)=J_nl(ipol)+ion_vel(jpol,na)*&
&becprd(ipol,jpol)%r(ikb,ibnd)*becp%r(ikb,ibnd)*deeq(ih,ih,na,1)

                         J_1(ipol)=J_1(ipol)+ion_vel(jpol,na)*&
&becprd(ipol,jpol)%r(ikb,ibnd)*becp%r(ikb,ibnd)*deeq(ih,ih,na,1)

!                         print*,'corrente non locale: ', J_nl(:)
!                          if (ionode) then
!                              print*,'becpr-ikb-ipol-ibnd-ityp',becpr(ipol)%r(ikb,ibnd),ikb,ipol,ibnd,ityp(na)
!                          end if
!                         print*,'VEL',ion_vel(:,na)
                        J_nl(ipol)=J_nl(ipol)+ion_vel(jpol,na)*&
&becpr(ipol)%r(ikb,ibnd)*becpd(jpol)%r(ikb,ibnd)*deeq(ih,ih,na,1)

                        J_2(ipol)=J_2(ipol)+ion_vel(jpol,na)*&
&becpr(ipol)%r(ikb,ibnd)*becpd(jpol)%r(ikb,ibnd)*deeq(ih,ih,na,1)
!                       print*,'corrente',J_nl(ipol)   
                     end do
                  end do
               end do 
            end do
            ijkb=ijkb+nh(nt)
        end if
 
   end do
end do
!
!il fattore due è fatto per la degenrazione di spin
 current(:)=current(:)+2.d0*J_nl(:)
!
!free memory
CALL deallocate_bec_type(  becp ) 
do ipol=1,3
   CALL deallocate_bec_type( becpr(ipol) )    
   CALL deallocate_bec_type( becpd(ipol) ) 
   do jpol=1,3
       CALL deallocate_bec_type(  becprd(ipol,jpol) ) 
   end do
end do
deallocate(xvkb)
deallocate(dvkb)
deallocate(xdvkb)
!
!
end subroutine add_nc_curr

