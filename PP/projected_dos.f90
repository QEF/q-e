!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

program doss
  !
  !     calcola la densita' degli stati e la proiezione di questa
  !     sugli orbitali atomici
  !

  use parameters, only : DP
  implicit none
  integer :: maxb, maxk, maxstep
  character (len=30) :: filename, nd * 2
  parameter (maxb = 500, maxk = 64, maxstep = 1001)
  real(kind=DP) :: ei, ef, de, efermi, en, en1, dummy, et (maxb, maxk), &
       peso (maxk), work (0:maxstep), dos (maxstep), dosp (maxstep, maxb) &
       , xdos (maxstep), xdosp (maxstep, maxb), alpha, total, inte, sup, &
       kx, ky, kz, proj (maxb, maxb, maxk)
  ! energie delle bande
  ! pesi dei punti k
  ! la densita` totale
  ! la densita` totale
  ! la densita` proiettata
  ! idem come sopra
  ! ma dopo lo smooth
  ! smooth in eV
  ! contiene le proiezioni
  integer :: e, nchi (4), lchi (4, 4), ityp (40), nstep, na, i, &
       kpoint, nb, npr, n, nat, npseu, np, nks, nbnd, nproj, nfermi, &
       index, npdummy
  namelist / input / ei, ef, alpha, nstep, filename, index
  print *, "InizioProgrammaprojdos"
  ! lower limit for the energy (below efer
  ei = 8
  ! upper limit (above efermi)
  ef = 5
  ! smoothing width
  alpha = 0.2
  ! number of data point between -ei and e
  nstep = 500
  !

  index = 1
  read (5, input)
  !      print*, "Hello, world"
  filename = 'wo3proj'
  ! file contenente le proiezioni
  open (4, file = filename)
  read (4, * ) nproj, nbnd, nks, efermi
  efermi = efermi * 13.6058
  read (4, * ) npseu, nat
  do np = 1, npseu
     read (4, * ) npdummy, nchi (np), (lchi (n, np), n = 1, nchi (np) )
     if (npdummy.ne.np) call errore ('projected_dos', &
          'wrong pseudopotential order', np)
  enddo
  read (4, * ) (ityp (na), na = 1, nat)
  ei = efermi - ei
  ef = efermi + ef
  !      print*, 'Tutto ok'
  do kpoint = 1, nks
     read (4, '(4f14.10)') kx, ky, kz, peso (kpoint)
     do nb = 1, nbnd
        read (4, '(f14.10)') et (nb, kpoint)
        read (4, '(8f10.8)') (proj (npr, nb, kpoint) , npr = 1, nproj + 1)
        !            print*, "Nb ", nb, ", Kpoint ",kpoint
     enddo
  enddo
  close (4)
  !      print*, "Letto tutto"
  ! passo
  de = (ef - ei) / nstep
  do e = 1, nstep
     en = ei + (e-1) * de
     en1 = ei + e * de
     if (en.le.efermi.and.en1.ge.efermi) nfermi = e
     dos (e) = 0.d0
     do npr = 1, nproj + 1
        dosp (e, npr) = 0.d0
     enddo
     do kpoint = 1, nks
        do nb = 1, nbnd
           if (et (nb, kpoint) .ge.en.and.et (nb, kpoint) .lt.en1) then
              do npr = 1, nproj + 1
                 dosp (e, npr) = dosp (e, npr) + proj (npr, nb, kpoint) * peso ( &
                      kpoint) / de
              enddo
              dos (e) = dos (e) + peso (kpoint) / de
           endif
        enddo
     enddo
  enddo
  call smooth (dos, dosp, nproj, xdos, xdosp, nstep, alpha, de, &
       work, maxstep)
  !
  !     write total dos
  !
  open (1, file = 'tot')
  open (2, file = 'projected')
  do e = 1, nstep
     en = ei + de * (e-0.5)
     write (1, * ) en - efermi, xdos (e)

  enddo
  !     do index = 1, nproj
  do e = 1, nstep
     en = ei + de * (e-0.5)
     write (2, * ) en - efermi, (xdosp (e, index), index = 1, nproj)
  enddo
  write (2, * )
  !     end do
  close (1)

  close (2)


end program doss
!***********************************************************************
subroutine smooth (dos, dosp, nproj, xdos, xdosp, nstep, alpha, &
     de, dummy, maxstep)
  !
  use parameters, only : DP
  implicit none
  integer :: maxstep, nstep, i, j, k, npr, nproj

  real(kind=DP) :: xdos (maxstep), dos (maxstep), dosp (maxstep, nproj), &
       xdosp (maxstep, nproj), alpha, de, x, y, somma, dummy (0:maxstep)
  dummy (0) = de
  do i = 1, nstep
     x = i * de
     dummy (i) = de * exp ( - x**2 / alpha**2)
  enddo
  do i = 1, nstep
     somma = 0.d0
     xdos (i) = 0.d0
     do j = 1, nstep
        k = abs (i - j)
        xdos (i) = xdos (i) + dos (j) * dummy (k)
        somma = somma + dummy (k)
     enddo
     xdos (i) = xdos (i) / somma

  enddo
  do i = 1, nstep
     do npr = 1, nproj + 1
        somma = 0.d0
        xdosp (i, npr) = 0.d0
        do j = 1, nstep
           k = abs (i - j)
           xdosp (i, npr) = xdosp (i, npr) + dosp (j, npr) * dummy (k)
           somma = somma + dummy (k)
        enddo
        xdosp (i, npr) = xdosp (i, npr) / somma
     enddo

  enddo
  return
end subroutine smooth
