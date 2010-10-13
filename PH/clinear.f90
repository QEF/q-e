
       SUBROUTINE clinear(nk1,nk2,nk3,nti,ntj,ntk,point,noint)
       USE kinds, ONLY : DP
       implicit none
       integer :: ll,iold,jold,kold,jnew,knew,istep,jstep,kstep
       integer :: ik1,ik2,ij1,ij2,ii1,ii2,nk1,nk2,nk3,ntk,ntj,nti
       integer :: npr,nkr
       integer :: np1, np2, np3
       complex(DP) :: point(*), noint(*)
!      nk1=nk2=nk3=32  -->  32768
!      np1=np2=np3=96  --> 884736

       nkr = nk1*nk2*nk3
       IF (nti==1.AND.ntj==1.AND.ntk==1) THEN
          noint(1:nkr)=point(1:nkr)
          RETURN
       ENDIF

       np1=nk1*nti
       np2=nk2*ntj
       np3=nk3*ntk
       npr = np1*np2*np3
       nkr = nk1*nk2*nk3

       ll = 0
       iold = 1
          jold = 1
             do kold = 1,nk3-1
                ik1 = kold
                ik2 = kold+1
                do kstep = 0,ntk-1
                   ll = ll+1
                   noint(ll) = point(ik1) + (point(ik2)-point(ik1))/ntk*kstep
                enddo !kstep
             enddo !kold
             ik1 = nk3
             ik2 = 1
             do kstep = 0,ntk-1
                ll=ll+1
                noint(ll) = point(ik1) + (point(ik2)-point(ik1))/ntk*kstep
             enddo

          do jold = 2,nk2
             ll=ll+np3*(ntj-1)
             do kold = 1,nk3-1
                ik1 = nk3*(jold-1) + kold
                ik2 = nk3*(jold-1) + kold+1
                do kstep = 0,ntk-1
                   ll = ll+1
                   noint(ll) = point(ik1) + (point(ik2)-point(ik1))/ntk*kstep
                enddo !kstep
             enddo !kold
             ik1 = jold*nk3
             ik2 = (jold-1)*nk3 + 1
             do kstep = 0,ntk-1
                ll = ll+1
                noint(ll) = point(ik1) + (point(ik2)-point(ik1))/ntk*kstep
             enddo
             ll=ll-np3*ntj
             do jstep=1,ntj-1
                do knew=1,np3
                   ll = ll+1
                   ij1 = (jold-2)*np3*ntj + knew
                   ij2 = (jold-1)*np3*ntj + knew
                   noint(ll) = noint(ij1) + (noint(ij2)-noint(ij1))/ntj*jstep
                enddo !knew
             enddo !jstep
             ll=ll+np3
          enddo !jold

          do jstep=1,ntj-1
             do knew=1,np3
                ll = ll+1
                ij1 = (nk2-1)*np3*ntj + knew
                ij2 = knew
                noint(ll) = noint(ij1) + (noint(ij2)-noint(ij1))/ntj*jstep
             enddo !knew
          enddo !jstep

          ll=ll+(nti-1)*np2*np3

       do iold = 2,nk1
          jold = 1
             do kold = 1,nk3-1
                ik1 = (iold-1)*nk2*nk3 + kold
                ik2 = (iold-1)*nk2*nk3 + kold+1
                do kstep = 0,ntk-1
                   ll = ll+1
                   noint(ll) = point(ik1) + (point(ik2)-point(ik1))/ntk*kstep
                enddo !kstep
             enddo !kold
             ik1 = (iold-1)*nk2*nk3 + nk3
             ik2 = (iold-1)*nk2*nk3 + 1
             do kstep = 0,ntk-1
                ll=ll+1
                noint(ll) = point(ik1) + (point(ik2)-point(ik1))/ntk*kstep
             enddo

          do jold = 2,nk2
             ll=ll+np3*(ntj-1)
             do kold = 1,nk3-1
                ik1 = (iold-1)*nk3*nk2 + nk3*(jold-1) + kold
                ik2 = (iold-1)*nk3*nk2 + nk3*(jold-1) + kold+1
                do kstep = 0,ntk-1
                   ll = ll+1
                   noint(ll) = point(ik1) + (point(ik2)-point(ik1))/ntk*kstep
                enddo !kstep
             enddo !kold
             ik1 = (iold-1)*nk2*nk3 + jold*nk3
             ik2 = (iold-1)*nk2*nk3 + (jold-1)*nk3 + 1
             do kstep = 0,ntk-1
                ll = ll+1
                noint(ll) = point(ik1) + (point(ik2)-point(ik1))/ntk*kstep
             enddo

             ll=ll-np3*ntj
             do jstep=1,ntj-1
                do knew=1,np3
                   ll = ll+1
                   ij1 = (iold-1)*np2*np3*nti + (jold-2)*np3*ntj + knew
                   ij2 = (iold-1)*np2*np3*nti + (jold-1)*np3*ntj + knew
                   noint(ll) = noint(ij1) + (noint(ij2)-noint(ij1))/ntj*jstep
                enddo !knew
             enddo !jstep
             ll=ll+np3
          enddo !jold

          do jstep=1,ntj-1
             do knew=1,np3
                ll = ll+1
                ij1 = (iold-1)*np2*np3*nti + (nk2-1)*np3*ntj + knew
                ij2 = (iold-1)*np2*np3*nti + knew
                noint(ll) = noint(ij1) + (noint(ij2)-noint(ij1))/ntj*jstep
             enddo !knew
          enddo !jstep

       ll=ll-nti*np2*np3

          do istep=1,nti-1
             do jnew=1,np2
                do knew=1,np3
                   ll = ll+1
                   ii1 = (iold-2)*np2*np3*nti + (jnew-1)*np3 + knew
                   ii2 = (iold-1)*np2*np3*nti + (jnew-1)*np3 + knew
                   noint(ll) = noint(ii1) + (noint(ii2)-noint(ii1))/nti*istep
                enddo  !knew
             enddo  !jnew
          enddo !istep
          ll=ll+nti*np2*np3
       enddo !iold

          ll = ll - (nti-1)*np2*np3
          do istep=1,nti-1
             do jnew=1,np2
                do knew=1,np3
                   ll = ll+1
                   ii1 = (nk1-1)*np2*np3*nti + (jnew-1)*np3 + knew
                   ii2 = (jnew-1)*np3 + knew
                   noint(ll) = noint(ii1) + (noint(ii2)-noint(ii1))/nti*istep
                enddo  !knew
             enddo  !jnew
          enddo !istep

       RETURN
       END SUBROUTINE clinear
