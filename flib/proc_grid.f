!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


      subroutine proc_grid_init(nproc,npx,npy,npz)

      implicit none
      integer nproc,npx,npy,npz

      npx = 1
      npy = 1
      npz = nproc

      if((npy*npz*npx).ne.nproc) then
        call error("proc_grid_init"," npx*npy*npz .ne. nproc ",
     c       npx*npy*npz)
      end if

      return
      end


      subroutine grid_dist_init(idist,mpime,npx,npy,npz,pex,pey,pez)

      implicit none
      integer idist,mpime,npx,npy,npz,pex,pey,pez

      idist = 0
      pex = 0
      pey = 0
      pez = mpime

      return
      end

      subroutine mesh_desc_init(desca,mpime,nproc,nr1,nr2,nr3,
     c           nr1_l,nr2_l,nr3_l, idist,NFFTX1,NFFTX2)

      implicit none
      integer mpime,nproc
      integer desca(12),nr1,nr2,nr3,nr1_l,nr2_l,nr3_l
      integer idist,NFFTX1,NFFTX2
      integer err

      DESCA(1) = nr1
      DESCA(2) = nr2
      DESCA(3) = nr3
      DESCA(5) = mpime
      DESCA(6) = nproc
      DESCA(8) = nr1
      DESCA(9) = nr2

      if(err.ne.0) then
        call error(' mesh_desc_init ',' descinit3d ',err)
      end if                                          

      return
      end


