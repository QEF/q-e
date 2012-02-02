!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE open_dvscf_star_q( q_index )
  !-----------------------------------------------------------------------
  !
  ! This routine reads the dvscf at a point q calculated by ph.x in
  ! the modes basis and
  !    (i)  transform it in cartesian coordinates, 
  !    (ii) obtain dvscf at a point q' in the star of q by using
  !         symmetry operations 
  !    (iii) dumps it in a special directory
  !
  !    Original routine written by Matteo Calandra and Gianni Profeta,
  !          adapted to ldisp case, fractional translations and imq=.true.
  !          by Matteo. Calandra.
  !
  !-----------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE cell_base,     ONLY : omega,at, bg, celldm, ibrav, tpiba2
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau, atm, amass
  USE wvfct,     ONLY : npwx,npw,igk
  USE symm_base,  ONLY : s, ftau,nsym,irt, invs
  USE lsda_mod, ONLY: nspin
  USE phcom
  USE el_phon
  USE io_global, ONLY : stdout , ionode
  use io_files, only: prefix,  tmp_dir, nd_nmbr, diropn
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE mp_global, ONLY : my_pool_id, npool, kunit,  me_pool, root_pool, intra_image_comm
  USE mp,  ONLY : mp_barrier, mp_bcast
  USE constants, ONLY: tpi
  USE control_flags, ONLY : modenum, noinv 
  USE modes,  ONLY : npert, nirr
  USE units_ph, ONLY : iudyn, lrdrho, iudvscf
  USE fft_base,  ONLY : cgather_sym, dfftp
  USE disp, ONLY : x_q
  USE control_ph, ONLY : ldisp, tmp_dir_ph
  USE xml_io_base,     ONLY : create_directory
!  USE noncollin_module, ONLY : nspin
  implicit none
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  integer :: q_index
!

  logical  :: exst, ascii_dvscf
  integer                   :: iq, ichosen_sym(48), lstr_iq(2*48), isym,i,j,ipol
  integer                   :: nsym_rot, q_rot, imode0, irr, is, ipert, na, nar,isym2
  integer                   :: k,n,nn, ri,rj,rk,iudvrot_asc,iudvrot, iudvrot_asc_imq,iudvrot_imq
  integer :: nq, isq (48), imq, index0,index_start_q, index_inequiv
  integer :: iu_qp, ierr
  ! nq :  degeneracy of the star of q
  ! isq: index of q in the star of a given sym.op.
  ! imq: index of -q in the star of q (0 if not present)
  real(DP) :: sxq (3, 48)
  ! list of vectors in the star of q
  real(DP) :: xq_tau,sxq_tau
  complex(dp) :: phase_xq
  COMPLEX(DP), ALLOCATABLE  :: phase_sxq(:) 
  COMPLEX(DP), ALLOCATABLE :: dvscfin(:,:,:), dvscf_at(:,:,:)
  COMPLEX(DP), ALLOCATABLE ::dvrot (:,:,:),dvrot_scr(:,:,:)
  CHARACTER(len=10)         :: label_q(2*48)
  CHARACTER(len=10)     :: star_i, star_minus_i
  CHARACTER(LEN=256) :: fildvscfrot,fildvscfrot_asc
  CHARACTER(LEN=256) :: fildvscfrot_imq,fildvscfrot_asc_imq
  CHARACTER(LEN=256) :: rdvscf_dir

  
  CALL start_clock ('open_dvscf_star_q')

  ascii_dvscf=.false.

  !
  ! I create (if it does not exist already) the 
  ! directory where the rotated dvscf are to be set
  !

  rdvscf_dir=trim(tmp_dir_ph)//'Rotated_DVSCF/'

  IF (ionode) inquire (file =TRIM(rdvscf_dir)//'.', exist = exst)
  CALL mp_bcast( exst, ionode_id, intra_image_comm )
  if(.not.exst) CALL create_directory(rdvscf_dir)
   

  IF ( me_pool.eq.root_pool ) THEN 
     iu_qp=find_free_unit()
     OPEN( UNIT = iu_qp, FILE = trim(rdvscf_dir)//'Q_POINTS.D', &
          STATUS = 'unknown', IOSTAT = ierr)
     IF( ierr /= 0 ) then 
        CALL errore( ' open in open_dvscf_star_q', &
             ' opening file '//trim(rdvscf_dir)//'Q_POINTS.D', 1 )
     END IF
  END IF
  
  !
  ! I identify the label of the q' in the open zone.
  ! If ldisp is not set this label simply goes from 
  ! 1 to nq in the star 
  !
  
  if(ldisp) then
     
    index_start_q=0 
    if(q_index.gt.1) then

       do iq=1,q_index-1
          call star_q(x_q(1,iq), at, bg, nsym , s , invs , nq, sxq, &
               isq, imq, .FALSE. )

          IF ( me_pool.eq.root_pool ) THEN           
             do q_rot=1,nq
                write(iu_qp,'(3f16.10,2i4)') (sxq(j,q_rot),j=1,3),index_start_q+q_rot,iq
             enddo
             if(imq == 0) then
                do q_rot=nq+1,2*nq
                   write(iu_qp,'(3f16.10,2i4)') (-sxq(j,q_rot-nq),j=1,3),index_start_q+q_rot,iq
                enddo
             endif
          endif
          index_inequiv=iq

          if (imq == 0) then
             index_start_q=index_start_q+2*nq
          else
             index_start_q=index_start_q+nq
          endif

       enddo
    endif
  else
    index_start_q=0
  endif
 

!  write(stdout,*) 'index_start_q=',index_start_q

  allocate(phase_sxq(nat))
  phase_sxq=CMPLX(1.d0,0.d0)
  !
  ! Here I calculate the q' in the star of q of the q that I need.
  !

  call star_q(xq, at, bg, nsym , s , invs , nq, sxq, isq, imq, .TRUE. )

  IF ( me_pool.eq.root_pool ) THEN           
     do iq=1,nq
        write(iu_qp,'(3f16.10,2i4)') (sxq(j,iq),j=1,3),index_start_q+iq,index_inequiv+1
     enddo
     if(imq == 0) then
        do iq=nq+1,2*nq
           write(iu_qp,'(3f16.10,2i4, a4)') &
                (-sxq(j,iq-nq),j=1,3),index_start_q+iq, index_inequiv+1, '  -q'
        enddo
     endif
     close(iu_qp)
  endif
  
  !
  ! Define labels for output files
  !

  
  do iq=1,nq
     k=iq+index_start_q
!     write(stdout,*) 'k=',k
     write(label_q(iq),'(1I4)') k
     label_q(iq)=trim(adjustl(label_q(iq)))
     call define_star_string(k,label_q(iq),lstr_iq(iq))
     write(stdout,*) 'label_q(',iq,')=',label_q(iq),'lstr_iq=',lstr_iq(iq)
  enddo

  if (imq == 0) then
     do iq=nq+1,2*nq
        k=iq+index_start_q
        write(label_q(iq),'(1I4)') k
        label_q(iq)=trim(adjustl(label_q(iq)))
        call define_star_string(k,label_q(iq),lstr_iq(iq))
        write(stdout,*) 'label_q(',iq,')=',label_q(iq),'lstr_iq=',lstr_iq(iq)
     enddo
  endif



  !
  !  Between all the possible symmetries I chose the first one
  !        (all of them lead to the same rotated dvscf)
  !
  
  do iq=1,nq
     nsym_rot=0
     do isym=1,nsym
        if (isq(isym) == iq) then
           nsym_rot=nsym_rot+1
           if(nsym_rot.eq.1) ichosen_sym(iq)=isym
        endif
     enddo
     if(nsym_rot.eq.0) then
        call errore('elphonstar','no symmetry relates q at star(q)',iq)
     endif
  enddo

!  do iq=1,nq
!     write(stdout,*) iq,ichosen_sym(iq)
!  enddo

  write(stdout,*) 'Chosen symmetries'
  write(stdout,*) '================='

  do iq=1,nq
     isym2=invs(ichosen_sym(iq))
     write(stdout,*) '-----------------------------------------'
     write(stdout,'(1a3,3f16.8)') 'xq=',(xq(j),j=1,3)
     write(stdout,'(1a4,3f16.8)') 'sxq=',(sxq(j,iq),j=1,3)
     write(stdout,*) '-----------------------------------------'
     write(stdout,*) 'S matrix is'
     do i=1,3
        write(stdout,'(3I14)') (s(i,j,isym2) ,j=1,3)
     enddo
     write(stdout,*) '-----------------------------------------'
     write(stdout,*) 'fractional translation is'
     write(stdout,'(3I4)') (ftau(j,isym2),j=1,3)
  enddo
  
  write(stdout,*) '==============================================='

  ALLOCATE (dvscf_at (dfftp%nr1x * dfftp%nr2x * dfftp%nr3x , nspin , 3*nat ))
  ALLOCATE (dvrot ( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x , nspin , 3*nat) )
  ALLOCATE (dvrot_scr ( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x , nspin , 3*nat) )

  dvscf_at=CMPLX(0.d0,0.d0)
  dvrot=CMPLX(0.d0,0.d0)
  dvrot_scr=CMPLX(0.d0,0.d0)
  
  imode0 = 0

  DO irr = 1, nirr
     
     
!     ALLOCATE (dvscfin ( nr1x * nr2x * nr3x , nspin , npert(irr)) )
     ALLOCATE (dvscfin ( dfftp%nnr , nspin , npert(irr)) )
     dvscfin=CMPLX(0.d0,0.d0)
     DO is = 1, nspin
        DO ipert = 1, npert (irr)
           CALL davcio_drho ( dvscfin(1,is,ipert),  lrdrho, iudvscf, &
                imode0 + ipert, -1 )
        
        END DO
     ENDDO

! for several pe cancel line below
!     dvrot(:,:,1:npert(irr))=dvscfin(:,:,1:npert(irr))

#if defined (__MPI)


     DO ipert = 1, npert(irr)
        DO is = 1, nspin
           CALL cgather_sym (dvscfin (:, is, ipert), dvrot (:, is, ipert) )
!           call cgather_sym( dvscfin(1,is,ipert), dvrot (1,is,ipert) )
        ENDDO
     ENDDO
!     call mp_barrier(intra_image_comm)
#else
     dvrot(:,:,1:npert(irr))=dvscfin(:,:,1:npert(irr))
#endif



     do ipert=1, npert(irr)
        dvscf_at(:,:,imode0+ipert)=dvrot(:, :, ipert)
     enddo
     
     imode0 = imode0 + npert(irr)
     
     DEALLOCATE (dvscfin)

  enddo


  ! 
  !Note that this dvscf is in the pw patterns.
  !   The pw< patterns are in the vector u(3*nat,3*nat) =u_{j}^{PW}
  !   defined in pwcom.
  ! We need the pwscf in the x_qAalpha modes, thus we
  !   need the following transformation
  !
  ! \partial_{x_{qAalpha}} = \sum_{j} C_{Aalpha}^{j PW}(q)*
  !        \partial_{u_{j}^{PW}}
  !
  ! where
  !
  ! C_{Aalpha}^{j PW}(q)=<u_j^{PW}|x_{qA alpha}>
  !
  

  dvrot=CMPLX(0.d0,0.d0) 


  do i=1,3*nat     
     do j=1,3*nat
        dvrot(:,:,i)=dvrot(:,:,i)+ &
             dvscf_at(:,:,j)*conjg(u(i,j))
     enddo
  enddo

  
  do i=1,3*nat
     dvscf_at(:,:,i)=(0.d0,0.d0)
  enddo


!################## FOR DEBUGGING ################################
  !
  ! Write out the unsymmetrized dvscf (this is for check)
  ! 
!
!  if(me_pool.eq.root_pool) then
  
!     open(unit=7878,file='dvscf_3natcart.dat',status='unknown')  
     
!     do na=1,nat
!        do ipol=1,3
!           irr=(na-1)*3+ipol
!           do  k = 1, nr3
!              do j = 1, nr2
!                 do i = 1, nr1
!                    n=(i-1) + (j-1)*nr1 + (k-1)*nr2*nr1 + 1 
                    
!                    write(7878,'(5i6,3f14.8)')  na,ipol , i,j,k, dvrot(n,1,irr), &
!                         sqrt(real(dvrot(n, 1, irr))**2+ &
!                         aimag(dvrot(n, 1, irr))**2)
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo
    
!     close(7878)

!  endif
!
!##########################################################################


  

  ! I transform dvscf in crystalline coordinates
  ! (necessary in order to apply s)

 
  do i=1,nat
     na=(i-1)*3
     do j=1,3
        dvscf_at(:,:,na+j)=dvrot(:,:,na+1)*at(1,j) +&
             dvrot(:,:,na+2)*at(2,j) + dvrot(:,:,na+3)*at(3,j) 
     enddo
  enddo

  !
  ! I take away the phase due to the q-point
  !

  dvrot=CMPLX(0.d0,0.d0)

  do i=1,nat
    xq_tau=(xq(1)*tau(1,i)+xq(2)*tau(2,i)+xq(3)*tau(3,i))*tpi
    phase_xq= CMPLX (cos(xq_tau),sin(xq_tau))
    do ipol=1,3
       imode0=(i-1)*3+ipol
       dvrot(:,:,imode0)=dvscf_at(:,:,imode0)*phase_xq
    enddo
  enddo



  dvscf_at=dvrot

  do q_rot=1,nq
     dvrot=CMPLX(0.d0,0.d0)  
     
     !
     !  Opening files f
     !
     star_i=label_q(q_rot)
     write(stdout,*) 'star_i=',star_i
     fildvscfrot='dvscf_sym_q'//star_i(1:lstr_iq(q_rot))//'_'
     fildvscfrot_asc=trim(rdvscf_dir)//trim(prefix)//"."//'dvscf_asc_sym_q'//star_i(1:lstr_iq(q_rot))//'_'//nd_nmbr 
     if (imq == 0) then
        star_minus_i=label_q(q_rot+nq)
        fildvscfrot_imq='dvscf_sym_q'//star_minus_i(1:lstr_iq(q_rot+nq))//'_'
        fildvscfrot_asc_imq=trim(rdvscf_dir)//trim(prefix)//"."//'dvscf_asc_sym_q'//star_minus_i(1:lstr_iq(q_rot+nq))//'_'//nd_nmbr 
     endif
     
     IF ( me_pool.eq.root_pool.and.ascii_dvscf ) THEN
        iudvrot_asc=find_free_unit()
        open(unit=iudvrot_asc,file=fildvscfrot_asc,status='unknown')
        if (imq == 0) then
           iudvrot_asc_imq=find_free_unit()
           open(unit=iudvrot_asc_imq,file=fildvscfrot_asc_imq,status='unknown')
        endif
     endif
     
     IF ( me_pool.eq.root_pool )  then
        iudvrot=find_free_unit()
        CALL diropn (iudvrot, fildvscfrot, lrdrho, exst, rdvscf_dir)
        if (imq == 0) then
          iudvrot_imq=find_free_unit()
          CALL diropn (iudvrot_imq, fildvscfrot_imq, lrdrho, exst, rdvscf_dir)
        endif
     endif
  
     
     !
     ! not that below isym is S and isym2 refers to S^-1
     !
     
     
     isym=ichosen_sym(q_rot)
     isym2=invs(ichosen_sym(q_rot))
     
     do k=1,nat  
        sxq_tau=(sxq(1,q_rot)*tau(1,k)+ &
             sxq(2,q_rot)*tau(2,k)+ &
             sxq(3,q_rot)*tau(3,k))*tpi
        
        phase_sxq(k)=CMPLX (cos(sxq_tau),sin(sxq_tau))
     enddo
     
     
     do is=1,nspin
        do  k = 1, dfftp%nr3
           do j = 1, dfftp%nr2
              do i = 1, dfftp%nr1
                 
                 
                 !
                 ! Here I rotate r
                 !
                 
                 ri = s(1, 1, isym2) * (i - 1) + s(2, 1, isym2) * (j - 1) &
                      + s(3, 1, isym2) * (k - 1) - ftau (1, isym2)
                 
                 
                 rj = s(1, 2, isym2) * (i - 1) + s(2, 2, isym2) * (j - 1) &
                      + s(3,2, isym2) * (k - 1) - ftau (2, isym2)
                 
                 
                 
                 rk = s(1, 3, isym2) * (i - 1) + s(2,3, isym2) * (j - 1) &
                      + s(3, 3, isym2) * (k - 1) - ftau (3, isym2)
                 
                 ri = mod (ri, dfftp%nr1) + 1
                 
                 rj = mod (rj, dfftp%nr2) + 1                 
                 
                 rk = mod (rk, dfftp%nr3) + 1
                 
                 
                 if (ri < 1) then 
                    ri = ri + dfftp%nr1
                 endif
                 
                 if (rj < 1) then 
                    rj = rj + dfftp%nr2
                 endif
                 
                 if (rk < 1) then 
                    rk = rk + dfftp%nr3
                 endif
                 
                 n=(i-1) + (j-1)*dfftp%nr1 + (k-1)*dfftp%nr2*dfftp%nr1 + 1
                 nn=(ri-1) + (rj-1)*dfftp%nr1 + (rk-1)*dfftp%nr2*dfftp%nr1 + 1
                 
                 
                 do na=1,nat
                    nar=irt(isym2,na)

                    index0=(nar-1)*3
                    
                    do ipol=1,3
                       imode0=(na-1)*3+ipol
                       
                       dvrot(n,is,imode0)=dvrot(n,is,imode0)+ &
                            ( s (ipol, 1, isym2) * dvscf_at (nn,is,index0+1) + &
                            s (ipol, 2, isym2) * dvscf_at (nn,is,index0+2) + &
                            s (ipol, 3, isym2) * dvscf_at (nn,is,index0+3) )
                       
                       
                    enddo
                 enddo
                 
              enddo
           enddo
        enddo
        
     enddo
     

     do na=1,nat
        do ipol=1,3
           imode0=(na-1)*3+ipol
           dvrot_scr(:,:,imode0 )=dvrot(:,:,imode0)/phase_sxq(na)
        enddo
     enddo
     
     
  !
  ! Back to cartesian coordinates
  !

     dvrot=CMPLX(0.d0,0.d0)
     do i=1,nat
        imode0=(i-1)*3
        do j=1,3
           dvrot(:,:,imode0+j)=dvrot_scr(:,:,imode0+1)*bg(j,1) +&
                dvrot_scr(:,:,imode0+2)*bg(j,2) + dvrot_scr(:,:,imode0+3)*bg(j,3)
        enddo
     enddo
     

        

     do na=1,nat
        do ipol=1,3
           imode0=(na-1)*3+ipol
           if(me_pool.eq.root_pool) CALL davcio( dvrot(1,1,imode0), lrdrho, iudvrot, imode0, + 1 )
        enddo
     enddo
        

     if(imq == 0) then
        dvrot_scr=conjg(dvrot)
        do na=1,nat
           do ipol=1,3
              imode0=(na-1)*3+ipol
              if(me_pool.eq.root_pool) CALL davcio( dvrot_scr(1,1,imode0), lrdrho, iudvrot_imq, imode0, + 1 )
           enddo
        enddo
     endif

           
!############ FOR DEBUGGING
!


     !  <MCB> Ascii writing of dvscf
     if (me_pool.eq.root_pool.and.ascii_dvscf) then
        do na=1,nat
           do ipol=1,3
              irr=(na-1)*3+ipol
              do  k = 1, dfftp%nr3
                 do j = 1, dfftp%nr2
                    do i = 1, dfftp%nr1
                       
                       n=(i-1) + (j-1)*dfftp%nr1 + (k-1)*dfftp%nr2*dfftp%nr1 + 1
                       
                       write(iudvrot_asc,'(1i10,2f16.10)')   n, dvrot(n,1,irr)
                    enddo
                 enddo
              enddo
           enddo
        enddo
        if(imq == 0) then
           do na=1,nat
              do ipol=1,3
                 irr=(na-1)*3+ipol
                 do  k = 1, dfftp%nr3
                    do j = 1, dfftp%nr2
                       do i = 1, dfftp%nr1
                          
                          n=(i-1) + (j-1)*dfftp%nr1 + (k-1)*dfftp%nr2*dfftp%nr1 + 1
                          
                          write(iudvrot_asc_imq,'(1i10,2f16.10)')   n, conjg(dvrot(n,1,irr))
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        endif
        !           
        close(iudvrot_asc)
        if(imq == 0) close(iudvrot_asc_imq)
     endif

! ##################################################################
     IF ( me_pool.eq.root_pool )   then
        CLOSE( UNIT = iudvrot, STATUS = 'KEEP' )
        if(imq == 0)  CLOSE( UNIT = iudvrot_imq, STATUS = 'KEEP' )
     endif
           
  enddo  !end loop over the nq
  

  


  deallocate(phase_sxq)
  DEALLOCATE(dvrot_scr)
  DEALLOCATE(dvrot)
  DEALLOCATE(dvscf_at)


  CALL stop_clock ('open_dvscf_star_q')
  RETURN

END SUBROUTINE open_dvscf_star_q


subroutine define_star_string(index,string,lstr)
  implicit none
  integer index,lstr
  character(len=10) string

! here put a check on the string length


  if(index < 10) then
     WRITE( string(1:1), '(I1)' ) index
     lstr=1
  elseif(index < 100) then
     WRITE( string(1:2), '(I2)' ) index
     lstr=2
  elseif(index < 1000) then
     WRITE( string(1:3), '(I3)' ) index
     lstr=3
  endif

  end subroutine define_star_string
!
