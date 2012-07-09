!
!-----------------------------------------------------------------------
subroutine write_eigenvectors (nat,ntyp,amass,ityp,q,w2,z,iout)
  !-----------------------------------------------------------------------
  !
  !   write modes on output file in a readable way
  !
  use kinds, only: dp
  use constants, only: amu_ry
  implicit none
  ! input
  integer nat, iout,ntyp
  integer ityp(nat)
  real(DP) q(3), w2(3*nat),amass(ntyp)
  complex(DP) z(3*nat,3*nat)
  ! local
  integer nat3, na, nta, ipol, i, j
  real(DP):: freq(3*nat)
  real(DP):: rydthz,rydcm1,cm1thz,znorm
  !
  nat3=3*nat
  !
  !  conversion factors RYD=>THZ, RYD=>1/CM e 1/CM=>THZ
  !
  rydthz = 13.6058*241.796
  rydcm1 = 13.6058*8065.5
  cm1thz = 241.796/8065.5
  !
  !  write frequencies and normalised displacements
  !
  write(iout,'(5x,''diagonalizing the dynamical matrix ...''/)')
  write(iout,'(1x,''q = '',3f12.4)') q
  write(iout,'(1x,74(''*''))')

 do i = 1,nat3
    do na = 1,nat
       nta = ityp(na)
       do ipol = 1,3
          z((na-1)*3+ipol,i) = z((na-1)*3+ipol,i)* sqrt(amu_ry*amass(nta))
       end do
    end do
 end do

  do i = 1,nat3
     !
     freq(i)= sqrt(abs(w2(i)))*rydcm1
     if (w2(i).lt.0.0) freq(i) = -freq(i)
     write (iout,9010) i, freq(i)*cm1thz, freq(i)
     do na = 1,nat
        write (iout,9020) (z((na-1)*3+ipol,i),ipol=1,3)
     end do
     !
  end do
  write(iout,'(1x,74(''*''))')
  !
  !      if (flvec.ne.' ') then
  !         open (unit=iout,file=flvec,status='unknown',form='unformatted')
  !         write(iout) nat, nat3, (ityp(i),i=1,nat), (q(i),i=1,3)
  !         write(iout) (freq(i),i=1,nat3), ((z(i,j),i=1,nat3),j=1,nat3)
  !         close(iout)
  !      end if
  !
  return
  !
9010 format(5x,'omega(',i2,') =',f15.6,' [THz] =',f15.6,' [cm-1]')
9020 format (1x,'(',3 (f10.6,1x,f10.6,3x),')')
  !
end subroutine write_eigenvectors
!
