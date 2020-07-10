subroutine pimd_get_force_from_pw(forcetmp)

   use pimd_variables, only : nbeadMD,natMD,ndimMD
   use path_variables, only : grad_pes
   implicit none
   integer k,iat,i,cc
   real(8) :: forcetmp(ndimMD,natMD,nbeadMD)
   
   forcetmp=0.d0
   DO k=1,nbeadMD
     cc=0
     DO iat=1,natMD
       DO i=1,ndimMD
         cc=cc+1
         forcetmp(i,iat,k)=-grad_pes(cc,k)
       END DO
     END DO
   END DO
   
   return

end subroutine pimd_get_force_from_pw

subroutine pimd_get_pot_from_pw(epMD)
   
   use path_variables, only : pes  
   use pimd_variables, only : nbeadMD
   implicit none
   real(8) :: epMD
     
   if (nbeadMD.gt.1) then
     epMD=sum(pes)/nbeadMD
   else
     epMD=pes(1)
   end if
   
   return
end subroutine pimd_get_pot_from_pw
   
subroutine pimd_pw_convert_pos(abc)
 
 ! use path_input_parameters_module, only : pos
  use path_variables, only : pos
  use pimd_variables, only : rpos,rcentroid,nbeadMD,natMD,ndimMD
  USE io_global, ONLY : meta_ionode, meta_ionode_id
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  
  implicit none
  integer cc,k,iat,i
  character*8 abc
  real(8), allocatable :: rpostmp(:,:,:)
  
  
  if (abc .eq. 'pw_to_md') then
    
    do k=1,nbeadMD
      cc=0
      DO iat=1,natMD
        DO i=1,ndimMD
          cc=cc+1
          rpos(i,iat,k)=pos(cc,k)
        END DO
      END DO
    end do 

    allocate(rpostmp(ndimMD,natMD,nbeadMD))
    rpostmp=0.0
    
    if (nbeadMD.gt.1) then
      do k=1,nbeadMD
        call refold(k,rpostmp(:,:,k))
      end do
    else
      rpostmp=rpos      
    end if
    
    if (nbeadMD.gt.1) then
      rcentroid=0.0
      do k=1,nbeadMD
        DO iat=1,natMD
          DO i=1,ndimMD
            rcentroid(i,iat)=rcentroid(i,iat)+rpostmp(i,iat,k)
          END DO
        END DO
      end do
      rcentroid=rcentroid/nbeadMD
    end if  
    deallocate(rpostmp)
     
  else
    
    do k=1,nbeadMD
      cc=0
      DO iat=1,natMD
        DO i=1,ndimMD
          cc=cc+1
          pos(cc,k)=rpos(i,iat,k)
        END DO
      END DO
    end do  
  !  write(*,*) 'pos1',pos
  ! CALL mp_bcast( pos,  meta_ionode_id, world_comm )
      
  end if 
  
  return
end subroutine pimd_pw_convert_pos


subroutine pimd_restart_traj()
 
 ! use path_input_parameters_module, only : pos
  use path_variables, only : pos
  use pimd_variables, only : rpos,vel,nbeadMD,natMD,ndimMD
  use pimd_variables, only : unit_dot_positions, unit_dot_velocities
  USE path_io_units_module,         ONLY : iunpath
  
  implicit none
  integer cc,k,iat,i,iflagerr,ngen
  real(8), allocatable :: vec_tmp(:)
  
  allocate(vec_tmp(natMD*ndimMD))

  !count the number of total available snapshots
  ngen=0
  rewind(unit_dot_positions)
  rewind(unit_dot_velocities)
  do while(.true.)
    read(unit_dot_positions,*,iostat=iflagerr)
    if(iflagerr/=0) exit
    ngen=ngen+1
  enddo
  rewind(unit_dot_positions)
  
  do i=1,ngen-nbeadMD
    READ(unit_dot_positions,*)
    READ(unit_dot_velocities,*)
  end do
  
  rpos=0.0
  do k=1,nbeadMD
    vec_tmp=0.0
    read(unit_dot_positions,*) vec_tmp(:)
    cc=0
    DO iat=1,natMD
      DO i=1,ndimMD
        cc=cc+1
        rpos(i,iat,k)=vec_tmp(cc)
      END DO
    END DO
  end do
  
    do k=1,nbeadMD
      cc=0
      DO iat=1,natMD
        DO i=1,ndimMD
          cc=cc+1
          pos(cc,k)=rpos(i,iat,k)
        END DO
      END DO
    end do 

  vel=0.0
  do k=1,nbeadMD
    vec_tmp=0.0
    read(unit_dot_velocities,*) vec_tmp(:)
    cc=0
    DO iat=1,natMD
      DO i=1,ndimMD
        cc=cc+1
        vel(i,iat,k)=vec_tmp(cc)
      END DO
    END DO
  end do
  
  close(unit_dot_positions)
  close(unit_dot_velocities)
  open(unit_dot_positions,file='positions.dat',position='APPEND',form='formatted')
  open(unit_dot_velocities,file='velocities.dat',position='APPEND',form='formatted')
  

  deallocate (vec_tmp)

  
  return
end subroutine pimd_restart_traj


SUBROUTINE pimd_get_amas_and_nat
  !
  USE kinds,         ONLY : DP
  !
  USE pimd_variables, ONLY : nbeadMD,ndimMD,mtot,&
                             natMD,amas,ion_name,fixcm,gMD
  use path_input_parameters_module, only : nat
  !
  USE ions_base, ONLY : amass, ityp, atm
  !
  IMPLICIT NONE
  !
  INTEGER :: iat

  natMD=nat
  if(fixcm) then 
    gMD=ndimMD*(natMD-1)      !   # degrees of freedom
  else 
    gMD=ndimMD*natMD
  endif
  ! is this really necessary? (GS)
  if(.not.allocated(amas)) allocate(amas(natMD))
  if(.not.allocated(ion_name)) allocate(ion_name(natMD))
  amas(:) = 0.0
  mtot=0.0
  DO iat=1,natMD
     amas(iat)=amass(ityp(iat))*10000.d0/5.48579909065d0
     mtot=mtot+amas(iat)
     ion_name(iat)=atm(ityp(iat))
  END DO
  
  RETURN
END SUBROUTINE pimd_get_amas_and_nat


SUBROUTINE match_neb_and_pimd_var
  use path_input_parameters_module, only : nstep_path, num_of_images,&
                                           first_last_opt
  use pimd_variables, only : nbeadMD,nblocks,nstep_block
  USE io_global, ONLY : meta_ionode,meta_ionode_id
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm

  implicit none
  
  nstep_path = nblocks*nstep_block
  num_of_images = nbeadMD
  if (nbeadMD.eq.1) num_of_images=2
  first_last_opt=.true.
  CALL mp_bcast( nstep_path,  meta_ionode_id, world_comm )
  CALL mp_bcast( num_of_images,  meta_ionode_id, world_comm )
  CALL mp_bcast( first_last_opt,  meta_ionode_id, world_comm )
  
  return
  
END SUBROUTINE match_neb_and_pimd_var


SUBROUTINE refold(idx,rpostmp)

  USE kinds,         ONLY : DP

  USE path_input_parameters_module, ONLY : alat 
  USE path_io_units_module,         ONLY : iunpath
  USE pimd_variables,               ONLY : natMD, ndimMD, nbeadMD, rpos
  !
  USE ions_base, ONLY : tau, ityp
  USE cell_base, ONLY : bg, at
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: idx
  !
  INTEGER :: iat
  REAL(DP), ALLOCATABLE :: pos0(:,:), pos1(:,:)
  REAL(8) :: rpostmp(ndimMD,natMD)
  
  ! atomic positions (in crystal units) of the previous and current image
  !
  ! tau is already in the internal QE units
  rpostmp=0.0
  rpostmp(:,:)=rpos(:,:,idx)
  !
  ALLOCATE( pos0(ndimMD,natMD), pos1(ndimMD,natMD) )
  !
  ! atomic positions in current image
  pos1(:,:) = rpos(:,:,idx) / alat
  CALL cryst_to_cart( natMD, pos1(1,1), bg, -1 )
  ! refold them within the unit cell around the origin
  pos1 = pos1(:,:) - anint(pos1(:,:))
  !
  IF ( idx > 1 ) THEN
     !
     ! atomic positions in previous image
     pos0(:,:) = rpos(:,:,idx-1) / alat
     CALL cryst_to_cart( natMD, pos0(1,1), bg, -1 )
     !
     DO iat = 1,natMD
        !
        ! translate atom by a lattice vector if needed
        ! N.B.: this solves the problem only when |p1-p0|<1.0
        !
           !IF ( ANY(ABS(pos1(:,iat) - pos0(:,iat)) > 0.5_DP) ) THEN
           !   WRITE ( iunpath, '(/,5x,A,I5,A,I3,A,I3,/,5x,A)' ) "WARNING: atom", iat, &
           !      " moved more than 1/2 alat from image", idx-1, " to image", idx, &
           !      "You can set minimum_image to true to avoid jumps in the path"
           !ENDIF
           WHERE( (pos1(:,iat) - pos0(:,iat)) > 0.5_DP )
              pos1(:,iat) = pos1(:,iat) - 1.0_DP
           ENDWHERE
           !
           WHERE( (pos1(:,iat) - pos0(:,iat)) < -0.5_DP ) 
              pos1(:,iat) = pos1(:,iat) + 1.0_DP
           ENDWHERE

     ENDDO
  ENDIF
  
  
  IF ( idx .eq. 1 ) THEN
     !
     ! atomic positions in previous image
     pos0(:,:) = rpos(:,:,nbeadMD) / alat
     CALL cryst_to_cart( natMD, pos0(1,1), bg, -1 )
     !
     DO iat = 1,natMD
        !
        ! translate atom by a lattice vector if needed
        ! N.B.: this solves the problem only when |p1-p0|<1.0
        !
           !IF ( ANY(ABS(pos1(:,iat) - pos0(:,iat)) > 0.5_DP) ) THEN
           !   WRITE ( iunpath, '(/,5x,A,I5,A,I3,A,I3,/,5x,A)' ) "WARNING: atom", iat, &
           !      " moved more than 1/2 alat from image", nbeadMD, " to image", idx, &
           !      "You can set minimum_image to true to avoid jumps in the path"
           !ENDIF
           WHERE( (pos1(:,iat) - pos0(:,iat)) > 0.5_DP )
              pos1(:,iat) = pos1(:,iat) - 1.0_DP
           ENDWHERE
           !
           WHERE( (pos1(:,iat) - pos0(:,iat)) < -0.5_DP ) 
              pos1(:,iat) = pos1(:,iat) + 1.0_DP
           ENDWHERE

     ENDDO
  ENDIF
  
     ! update positions for the temp. vector rpostmp
  CALL cryst_to_cart( natMD, pos1(1,1), at, 1 )
  rpostmp(:,:) =  pos1(:,:)*alat  
  
  DEALLOCATE( pos0, pos1 )

  RETURN
  !
END SUBROUTINE refold
!
