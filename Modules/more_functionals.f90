!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!     ==================================================================
      SUBROUTINE LSD_LYP(RHO,ETA,ELYP,VALYP,VBLYP)
!     ==--------------------------------------------------------------==
!     ==  C. LEE, W. YANG, AND R.G. PARR, PRB 37, 785 (1988)          ==
!     ==  THIS IS ONLY THE LDA PART                                   ==
!     ==--------------------------------------------------------------==
      USE kinds, ONLY: DP
!
      IMPLICIT NONE
!     arguments
      REAL(DP) :: RHO,ETA,ELYP,VALYP,VBLYP
!     locals
      REAL(DP) :: RA,RB,RM3,DR,E1,OR,DOR,E2,DE1A,DE1B,DE2A,DE2B
      REAL(DP), PARAMETER :: SMALL=1.D-24, A=0.04918D0, B=0.132D0, &
                         C=0.2533D0, D=0.349D0, CF=2.87123400018819108D0
!     ==--------------------------------------------------------------==
      RA=RHO*0.5D0*(1.D0+ETA)
      RA=MAX(RA,SMALL)
      RB=RHO*0.5D0*(1.D0-ETA)
      RB=MAX(RB,SMALL)
      RM3=RHO**(-1.D0/3.D0)
      DR=(1.D0+D*RM3)
      E1=4.D0*A*RA*RB/RHO/DR
      OR=EXP(-C*RM3)/DR*RM3**11.D0
      DOR=-1.D0/3.D0*RM3**4*OR*(11.D0/RM3-C-D/DR)
      E2=2.D0**(11.D0/3.D0)*CF*A*B*OR*RA*RB*(RA**(8.d0/3.d0)+ RB**(8.d0/3.d0))
      ELYP=(-E1-E2)/RHO
      DE1A=-E1*(1.D0/3.D0*D*RM3**4/DR+1./RA-1./RHO)
      DE1B=-E1*(1.D0/3.D0*D*RM3**4/DR+1./RB-1./RHO)
      DE2A=-2.D0**(11.D0/3.D0)*CF*A*B*(DOR*RA*RB*(RA**(8.d0/3.d0)+ &
            RB**(8.d0/3.d0))+OR*RB*(11.d0/3.d0*RA**(8.d0/3.d0)+ &
            RB**(8.d0/3.d0)))
      DE2B=-2.D0**(11.D0/3.D0)*CF*A*B*(DOR*RA*RB*(RA**(8.d0/3.d0)+ &
      RB**(8.d0/3.d0))+OR*RA*(11.d0/3.d0*RB**(8.d0/3.d0)+ &
      RA**(8.d0/3.d0)))
      VALYP=DE1A+DE2A
      VBLYP=DE1B+DE2B
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_LYP



!     ==================================================================
      SUBROUTINE LSD_PADE(RHO,ETA,EC,VCA,VCB)
!     ==--------------------------------------------------------------==
!     ==  PADE APPROXIMATION                                          ==
!     ==--------------------------------------------------------------==
      USE kinds, ONLY: DP
      IMPLICIT NONE
!     arguments
      REAL(DP) :: RHO,ETA,EC,VCA,VCB
!     locals
      REAL(DP) :: RS,FS,DFS,DFSA,DFSB,A0P,A1P,A2P,A3P,B1P,B2P,B3P,B4P
      REAL(DP) :: TOP,DTOP,TOPX,BOT,DBOT,BOTX,VC,DX
      REAL(DP), PARAMETER :: A0=.4581652932831429d0, A1=2.217058676663745d0, &
                A2=0.7405551735357053d0, A3=0.01968227878617998d0
      REAL(DP), PARAMETER :: B1=1.0D0, B2=4.504130959426697d0, &
                 B3=1.110667363742916d0, B4=0.02359291751427506d0
      REAL(DP), PARAMETER :: DA0=.119086804055547D0, DA1=.6157402568883345d0, &
                 DA2=.1574201515892867d0, DA3=.003532336663397157d0
      REAL(DP), PARAMETER :: DB1=0.0d0, DB2=.2673612973836267d0,  &
                 DB3=.2052004607777787d0, DB4=.004200005045691381d0
      REAL(DP), PARAMETER :: RSFAC=.6203504908994000d0, FSFAC=1.92366105093153617d0
!     ==--------------------------------------------------------------==
      RS=RSFAC*RHO**(-1.d0/3.d0)
      FS=FSFAC*((1.d0+ETA)**(4.d0/3.d0)+(1.d0-ETA)**(4.d0/3.d0)-2.d0)
      DFS=FSFAC*4.d0/3.d0* ((1.d0+ETA)**(1.d0/3.d0)-(1.d0-ETA)**(1.d0/3.d0))
      DFSA=DFS*(1.d0-ETA)
      DFSB=DFS*(-1.d0-ETA)
      A0P=A0+FS*DA0
      A1P=A1+FS*DA1
      A2P=A2+FS*DA2
      A3P=A3+FS*DA3
      B1P=B1+FS*DB1
      B2P=B2+FS*DB2
      B3P=B3+FS*DB3
      B4P=B4+FS*DB4
      TOP=A0P+RS*(A1P+RS*(A2P+RS*A3P))
      DTOP=A1P+RS*(2.d0*A2P+RS*3.d0*A3P)
      TOPX=DA0+RS*(DA1+RS*(DA2+RS*DA3))
      BOT=RS*(B1P+RS*(B2P+RS*(B3P+RS*B4P)))
      DBOT=B1P+RS*(2.d0*B2P+RS*(3.d0*B3P+RS*4.d0*B4P))
      BOTX=RS*(DB1+RS*(DB2+RS*(DB3+RS*DB4)))
      EC=-TOP/BOT
      VC=EC+RS*(DTOP/BOT-TOP*DBOT/(BOT*BOT))/3.d0
      DX=-(TOPX/BOT-TOP*BOTX/(BOT*BOT))
      VCA=VC+DX*DFSA
      VCB=VC+DX*DFSB
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_PADE



!     ==================================================================
      SUBROUTINE LSD_GLYP(RA,RB,GRHOAA,GRHOAB,GRHOBB,SC,  &
                                            V1CA,V2CA,V1CB,V2CB,V2CAB)
!     ==--------------------------------------------------------------==
      USE kinds, ONLY: DP
! LEE, YANG PARR: GRADIENT CORRECTION PART
      IMPLICIT NONE ! REAL(DP) (A-H,O-Z), INTEGER (I-N)
!     arguments
      REAL(DP) :: RA,RB,GRHOAA,GRHOAB,GRHOBB,SC, &
                  V1CA,V2CA,V1CB,V2CB,V2CAB
!     locals
      REAL(DP) :: RHO,RM3,DR,OR,DOR,DER,DDER
      REAL(DP) :: DLAA,DLAB,DLBB,DLAAA,DLAAB,DLABA,DLABB,DLBBA,DLBBB
      REAL(DP), PARAMETER :: A=0.04918D0,B=0.132D0,C=0.2533D0,D=0.349D0
!     ==--------------------------------------------------------------==
      RHO=RA+RB
      RM3=RHO**(-1.D0/3.D0)
      DR=(1.D0+D*RM3)
      OR=EXP(-C*RM3)/DR*RM3**11.D0
      DOR=-1.D0/3.D0*RM3**4*OR*(11.D0/RM3-C-D/DR)
      DER=C*RM3+D*RM3/DR
      DDER=1.d0/3.d0*(D*D*RM3**5/DR/DR-DER/RHO)
      DLAA=-A*B*OR*(RA*RB/9.d0*(1.d0-3*DER-(DER-11.d0)*RA/RHO)-RB*RB)
      DLAB=-A*B*OR*(RA*RB/9.d0*(47.d0-7.d0*DER)-4.d0/3.d0*RHO*RHO)
      DLBB=-A*B*OR*(RA*RB/9.d0*(1.d0-3*DER-(DER-11.d0)*RB/RHO)-RA*RA)
      DLAAA=DOR/OR*DLAA-A*B*OR*(RB/9.d0*(1.d0-3*DER-(DER-11.d0)*RA/RHO)- &
            RA*RB/9.d0*((3.d0+RA/RHO)*DDER+(DER-11.d0)*RB/RHO/RHO))
      DLAAB=DOR/OR*DLAA-A*B*OR*(RA/9.d0*(1.d0-3.d0*DER-(DER-11.d0)*RA/RHO)- &
            RA*RB/9.d0*((3.d0+RA/RHO)*DDER-(DER-11.d0)*RA/RHO/RHO)-2.d0*RB)
      DLABA=DOR/OR*DLAB-A*B*OR*(RB/9.d0*(47.d0-7.d0*DER)-7.d0/9.d0*RA*RB*DDER- &
            8.d0/3.d0*RHO)
      DLABB=DOR/OR*DLAB-A*B*OR*(RA/9.d0*(47.d0-7.d0*DER)-7.d0/9.d0*RA*RB*DDER- &
            8.d0/3.d0*RHO)
      DLBBA=DOR/OR*DLBB-A*B*OR*(RB/9.d0*(1.d0-3.d0*DER-(DER-11.d0)*RB/RHO)- &
            RA*RB/9.d0*((3.d0+RB/RHO)*DDER-(DER-11.d0)*RB/RHO/RHO)-2.d0*RA)
      DLBBB=DOR/OR*DLBB-A*B*OR*(RA/9.d0*(1.d0-3*DER-(DER-11.d0)*RB/RHO)- &
            RA*RB/9.d0*((3.d0+RB/RHO)*DDER+(DER-11.d0)*RA/RHO/RHO))
      SC=DLAA*GRHOAA+DLAB*GRHOAB+DLBB*GRHOBB
      V1CA=DLAAA*GRHOAA+DLABA*GRHOAB+DLBBA*GRHOBB
      V1CB=DLAAB*GRHOAA+DLABB*GRHOAB+DLBBB*GRHOBB
      V2CA=2.d0*DLAA
      V2CB=2.d0*DLBB
      V2CAB=DLAB
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_GLYP


!______________________________________________________________________
      subroutine ggablyp4(nnr,nspin,gradr,rhor,exc)
!     _________________________________________________________________
!     becke-lee-yang-parr gga
!     
!     exchange: becke, pra 38, 3098 (1988) but derived from
!        pw91 exchange formula given in prb 48, 14944 (1993)
!        by setting "b3" and "b4" to 0.0
!     correlation: miehlich et al., cpl 157, 200 (1989) 
!     method by ja white & dm bird, prb 50, 4954 (1994) 
!
!     spin-polarized version by andras stirling 10/1998,
!     using original gga program of alfredo pasquarello 22/09/1994
!     and spin-unpolarized blyp routine of olivier parisel and 
!     alfredo pasquarello (02/1997)  
!
      USE kinds, ONLY : DP
      USE constants, ONLY: pi, fpi
!
      implicit none
! input
      integer nspin, nnr
      real(DP) gradr(nnr,3,nspin), rhor(nnr,nspin)
! output
! on output: rhor contains the exchange-correlation potential
      real(DP)  exc
! local
      integer isdw, isup, isign, ir
!
      real(DP) abo, agdr, agdr2, agr, agr2, agur, agur2, arodw,     &
           arodw2, aroe, aroe2, aroup, aroup2, ax
      real(DP) byagdr, byagr, byagur, cden, cf, cl1, cl11, cl2,     &
           cl21, cl22, cl23, cl24, cl25, cl26, cl27, clyp, csum
      real(DP) dddn, dexcdg, dexcdgd, dexcdgu, df1d, df1u, df2d,    &
          df2u, dfd, dfnum1d, dfnum1u, dfnum2d, dfnum2u, dfs, dfu,     &
          dfxdd, dfxdg, dfxdgd, dfxdgu, dfxdu, dilta, dilta119, dl1dn, &
          dl1dnd, dl1dnu, dl2dd, dl2dg, dl2dgd, dl2dgu, dl2dn,         &
          dl2dnd, dl2dnd1, dl2dnu, dl2dnu1, dl2do, dlt, dodn,          &
          disign, dwsign, dys, dysd, dysu
      real(DP) ex, excupdt, exd, exu, fac1, fac2, factor1, factor2, &
          fx, fxd, fxden, fxdend, fxdenu, fxnum, fxnumd, fxnumu, fxu
      real(DP) gkf, gkfd, gkfu, grdx, grdy, grdz, grux, gruy, gruz, &
          grx, gry, grz
      real(DP) omiga, pd, pi2, pider2, piexch, pu
      real(DP) rhodw, rhoup, roe, roedth, roeth, roeuth, rometh
      real(DP) s, s2, sd, sd2, sddw, sdup, su, su2, sysl, sysld, syslu
      real(DP) t113, upsign, usign
      real(DP) x1124, x113, x118, x13, x143, x19, x23, x43,         &
          x4718, x53, x672, x718, x772, x83
      real(DP) ys, ysd, ysl, ysld, yslu, ysr, ysrd, ysru, ysu   
!===========================================================================
      real(DP) bb1, bb2, bb5, aa, bb, cc, dd, delt, eps
      parameter(bb1=0.19644797d0,bb2=0.2742931d0,bb5=7.79555418d0,            &
          aa=0.04918d0,                                                  &
          bb=0.132d0,cc=0.2533d0,dd=0.349d0,delt=1.0d-12,eps=1.0d-14)
!
!
      x13=1.0d0/3.0d0
      x19=1.0d0/9.0d0
      x23=2.0d0/3.0d0
      x43=4.0d0/3.0d0
      x53=5.0d0/3.0d0
      x83=8.0d0/3.0d0
      x113=11.0d0/3.0d0
      x4718=47.0d0/18.0d0
      x718=7.0d0/18.0d0
      x118=1.0d0/18.0d0
      x1124=11.0d0/24.0d0
      x143=14.0d0/3.0d0
      x772=7.0d0/72.0d0
      x672=6.0d0/72.0d0
!     
!     _________________________________________________________________
!     derived parameters from pi
!
      pi2=pi*pi
      ax=-0.75d0*(3.0d0/pi)**x13
      piexch=-0.75d0/pi 
      pider2=(3.0d0*pi2)**x13
      cf=0.3d0*pider2*pider2
!     _________________________________________________________________
!     other parameters
!
      t113=2.0d0**x113
!
      rhodw=0.0d0
      grdx=0.0d0
      grdy=0.0d0
      grdz=0.0d0
!
      fac1=1.0d0
!     _________________________________________________________________
!     main loop
!
      isup=1
      isdw=2
      do ir=1,nnr
         rhoup=rhor(ir,isup)
         grux=gradr(ir,1,isup)
         gruy=gradr(ir,2,isup)
         gruz=gradr(ir,3,isup)
         if(nspin.eq.2) then
            rhodw=rhor(ir,isdw)
            grdx=gradr(ir,1,isdw)
            grdy=gradr(ir,2,isdw)
            grdz=gradr(ir,3,isdw)
         else
            rhodw=0.0d0
            grdx =0.0d0
            grdy =0.0d0
            grdz =0.0d0
         endif
         roe=rhoup+rhodw
         if(roe.eq.0.0) goto 100
         aroup=abs(rhoup)
         arodw=abs(rhodw)
         aroe=abs(roe)
         grx=grux + grdx
         gry=gruy + grdy
         grz=gruz + grdz
         agur2=grux*grux+gruy*gruy+gruz*gruz
         agur=sqrt(agur2)
         agdr2=grdx*grdx+grdy*grdy+grdz*grdz
         agdr=sqrt(agdr2)
         agr2=grx*grx+gry*gry+grz*grz
         agr=sqrt(agr2)
         roeth=aroe**x13
         rometh=1.0d0/roeth
         gkf=pider2*roeth
         sd=1.0d0/(2.0d0*gkf*aroe)
         s=agr*sd
         s2=s*s
!     _________________________________________________________________
!     exchange 
!
         if(nspin.eq.1) then
!
!
            ysr=sqrt(1.0d0+bb5*bb5*s2)
            ys=bb5*s+ysr
            ysl=log(ys)*bb1
            sysl=s*ysl
            fxnum=1.0d0+sysl+bb2*s2
            fxden=1.0d0/(1.0d0+sysl)
            fx=fxnum*fxden
!
            ex=ax*fx*roeth*aroe
!
!     ### potential contribution ###
!
            dys=bb5*(1.0d0+bb5*s/ysr)/ys
            dfs=-fxnum*(ysl+bb1*s*dys)*fxden*fxden                      &
     &           +(ysl+bb1*s*dys+2.0d0*s*bb2)*fxden
            dfxdu=(ax*roeth*x43)*(fx-dfs*s)
            dfxdg=ax*roeth*dfs*sd
!
!     ### end of potential contribution ###
!     
         else
!
            roeuth=(2.0d0*aroup)**x13
            roedth=(2.0d0*arodw)**x13
            gkfu=pider2*roeuth*aroup
            gkfd=pider2*roedth*arodw
            upsign=sign(1.d0,gkfu-eps)
            dwsign=sign(1.d0,gkfd-eps)
            factor1=0.5d0*(1+upsign)/(gkfu+(1-upsign)*eps)
            fac1=gkfu*factor1
            factor2=0.5d0*(1+dwsign)/(gkfd+(1-dwsign)*eps)
            fac2=gkfd*factor2
            sdup=1.0d0/2.0d0*factor1
            sddw=1.0d0/2.0d0*factor2
            su=agur*sdup
            su2=su*su
            sd=agdr*sddw
            sd2=sd*sd
!
            ysru=sqrt(1.0d0+bb5*bb5*su2)
            ysu=bb5*su+ysru
            yslu=log(ysu)*bb1
            syslu=su*yslu
            fxnumu=1.0d0+syslu+bb2*su2
            fxdenu=1.0d0/(1.0d0+syslu)
            fxu=fxnumu*fxdenu
            exu=piexch*2.0d0*gkfu*fxu*fac1
!
            ysrd=sqrt(1.0d0+bb5*bb5*sd2)
            ysd=bb5*sd+ysrd
            ysld=log(ysd)*bb1
            sysld=sd*ysld
            fxnumd=1.0d0+sysld+bb2*sd2
            fxdend=1.0d0/(1.0d0+sysld)
            fxd=fxnumd*fxdend
            exd=piexch*2.0d0*gkfd*fxd*fac2
!
            ex=0.5d0*(exu+exd)
!
!     ### potential contribution ###
!
            dysu=bb5*(1.0d0+bb5*su/ysru)/ysu
            pu=2.0d0*su*bb2
            dfnum1u=yslu+bb1*su*dysu+pu
            df1u=dfnum1u*fxdenu
            dfnum2u=fxnumu*(yslu+bb1*su*dysu)
            df2u=dfnum2u*fxdenu*fxdenu
            dfu=df1u-df2u
            dfxdu=ax*roeuth*x43*1.0d0*(fxu-dfu*su)*fac1
            dfxdgu=ax*aroup*roeuth*dfu*sdup*fac1
!
            dysd=bb5*(1.0d0+bb5*sd/ysrd)/ysd
            pd=2.0d0*sd*bb2
            dfnum1d=ysld+bb1*sd*dysd+pd
            df1d=dfnum1d*fxdend
            dfnum2d=fxnumd*(ysld+bb1*sd*dysd)
            df2d=dfnum2d*fxdend*fxdend
            dfd=df1d-df2d
            dfxdd=ax*roedth*x43*1.0d0*(fxd-dfd*sd)*fac2
            dfxdgd=ax*arodw*roedth*dfd*sddw*fac2
!     
!     ### end of potential contribution ###
!
         endif
!     _________________________________________________________________
!     correlation lyp(aroe,aroup,arodw,agr,agur,agdr)
!
         cden=1.0d0+dd*rometh
         cl1=-aa/cden
!
         omiga=exp(-cc*rometh)/cden/aroe**x113
         dilta=rometh*(cc+dd/cden)
         aroe2=aroe*aroe
         abo=aa*bb*omiga
!
         dodn=x13*omiga/aroe*(dilta-11.0d0)
         dddn=x13*(dd*dd*aroe**(-x53)/cden/cden-dilta/aroe)
!
         if(nspin.eq.1) then
!
            cl1=cl1*aroe
!
            cl21=4.0d0*cf*aroe**x83
            cl22=(x4718-x718*dilta)*agr2
            cl23=(2.5d0-x118*dilta)*agr2/2.0d0
            cl24=(dilta-11.0d0)/9.0d0*agr2/4.0d0
            cl25=x1124*agr2
!
            cl2=-abo*aroe2*(0.25d0*(cl21+cl22-cl23-cl24)-cl25)
!
!     ### potential contribution ###
!
            dl1dnu=-aa*(1/cden+x13*dd*rometh/cden/cden)
!
            dlt=x672+2.0d0*x772*dilta
            dl2dn=-abo*aroe*(cf*x143*aroe**x83-dlt*agr2)
            dl2do=cl2/omiga
            dl2dd=abo*aroe2*x772*agr2
            dl2dnu=dl2dn+dl2do*dodn+dl2dd*dddn
!     
            dl2dg=abo*aroe2*agr*dlt
!
!     ### end of potential contribution ###
!
         else
!
            cl11=cl1*4.0d0/aroe
            cl1=cl11*aroup*arodw
!
            aroup2=aroup*aroup
            arodw2=arodw*arodw
!
            cl21=t113*cf*(aroup**x83+arodw**x83)
            cl22=(x4718-x718*dilta)*agr2
            cl23=(2.5d0-x118*dilta)*(agur2+agdr2)
            dilta119=(dilta-11.0d0)/9.0d0
            cl24=dilta119/aroe*(aroup*agur2+arodw*agdr2)
            cl25=x23*aroe2*agr2
            cl26=(x23*aroe2-aroup2)*agdr2
            cl27=(x23*aroe2-arodw2)*agur2
!
            csum=cl21+cl22-cl23-cl24
            cl2=-abo*(aroup*arodw*csum-cl25+cl26+cl27)
!
!     ### potential contribution ###
!
!     *** cl1 has changed its form! ***
!
            dl1dn=cl1/aroe*(x13*dd/cden*rometh-1.0d0)
            dl1dnu=dl1dn+cl11*arodw
            dl1dnd=dl1dn+cl11*aroup
!     
            dl2dnu1=arodw*csum+                                         &
     &           arodw*aroup*(t113*cf*x83*aroup**x53-                   &
     &           dilta119*arodw/aroe2*(agur2-agdr2))-x43*aroe*agr2+     &
     &           x23*agdr2*(2.0d0*arodw-aroup)+x43*aroe*agur2
            dl2dnd1=aroup*csum+                                         &
     &           aroup*arodw*(t113*cf*x83*arodw**x53+                   &
     &           dilta119*aroup/aroe2*(agur2-agdr2))-x43*aroe*agr2+     &
     &           x23*agur2*(2.0d0*aroup-arodw)+x43*aroe*agdr2
!
            dl2do=cl2/omiga
            dl2dd=-abo*aroup*arodw*                                     &
     &           (-x718*agr2+x118*(agur2+agdr2)-                        &
     &           x19*(aroup*agur2+arodw*agdr2)/aroe)
!
            dl2dnu=-abo*dl2dnu1+dl2do*dodn+dl2dd*dddn
            dl2dnd=-abo*dl2dnd1+dl2do*dodn+dl2dd*dddn
!
            dl2dg=-abo*                                                 &
     &           (aroup*arodw*2.0d0*(x4718-x718*dilta)*agr-               &
     &           x43*aroe2*agr) 
            dl2dgu=-2.0d0*abo*agur*((x118*dilta-2.5d0-                      &
     &           dilta119*aroup/aroe)*aroup*arodw                       &
     &           +x23*aroe2-arodw2)
            dl2dgd=-2.0d0*abo*agdr*((x118*dilta-2.5d0-                      &
     &           dilta119*arodw/aroe)*aroup*arodw                       &
     &           +x23*aroe2-aroup2)
!
         endif
!
         clyp=cl1+cl2
!     _________________________________________________________________
!     updating of xc-energy
!
         excupdt=ex+clyp
!
         exc=exc+excupdt
!
!     _________________________________________________________________
!     first part xc-potential construction
!
!
         rhor(ir,isup)=dfxdu+(dl1dnu+dl2dnu)*fac1
         isign=sign(1.d0,agr-delt)
         byagr=0.5d0*(1+isign)/(agr+(1-isign)*delt)
!
         if(nspin.eq.1) then
!
            dexcdg=(dfxdg*aroe+dl2dg)*byagr
            gradr(ir,1,isup)=grx*dexcdg
            gradr(ir,2,isup)=gry*dexcdg
            gradr(ir,3,isup)=grz*dexcdg
!
         else
!
            rhor(ir,isdw)=dfxdd+(dl1dnd+dl2dnd)*fac2
!
            usign =sign(1.d0,agur-delt)
            disign=sign(1.d0,agdr-delt)
            byagur=0.5d0*(1+ usign)/(agur+(1- usign)*delt)
            byagdr=0.5d0*(1+disign)/(agdr+(1-disign)*delt)
!
            dexcdgu=(dfxdgu+dl2dgu)*byagur
            dexcdgd=(dfxdgd+dl2dgd)*byagdr
            dexcdg=dl2dg*byagr
!
            gradr(ir,1,isup)=(dexcdgu*grux+dexcdg*grx)*fac1
            gradr(ir,2,isup)=(dexcdgu*gruy+dexcdg*gry)*fac1
            gradr(ir,3,isup)=(dexcdgu*gruz+dexcdg*grz)*fac1
            gradr(ir,1,isdw)=(dexcdgd*grdx+dexcdg*grx)*fac2
            gradr(ir,2,isdw)=(dexcdgd*grdy+dexcdg*gry)*fac2
            gradr(ir,3,isdw)=(dexcdgd*grdz+dexcdg*grz)*fac2
!
         endif
!          
 100     continue
      end do
!
      return
      end subroutine ggablyp4
!
!______________________________________________________________________
      subroutine ggapbe(nnr,nspin,gradr,rhor,excrho)
!     _________________________________________________________________
!     Perdew-Burke-Ernzerhof gga
!     Perdew, et al. PRL 77, 3865, 1996
!
      USE kinds, ONLY: DP
      use constants, only: pi, fpi
!
      implicit none
! input
      integer nspin, nnr
      real(DP)  gradr(nnr,3,nspin), rhor(nnr,nspin)
! output: excrho: exc * rho ;  E_xc = \int excrho(r) d_r
! output: rhor:   contains the exchange-correlation potential
      real(DP)  excrho
! local
      integer ir, icar, iss, isup, isdw, nspinx
      real(DP) lim1, lim2
      parameter ( lim1=1.d-8, lim2=1.d-8, nspinx=2 )
      real(DP) zet, arho(nspinx), grad(3,nspinx), agrad(nspinx),    &
          arhotot, gradtot(3), agradtot,                               &
          scl, scl1, wrkup, wrkdw,                                     &
          exrho(nspinx), dexdrho(nspinx), dexdg(nspinx),               &
          ecrho, decdrho(nspinx), decdg
!
!     main loop
!
      isup=1
      isdw=2
      do ir=1,nnr
!
         arho(isup) = abs(rhor(ir,isup))
         arhotot = arho(isup)
         zet = 0.d0
         do icar = 1, 3
            grad(icar,isup) = gradr(ir,icar,isup)
            gradtot(icar) = gradr(ir,icar,isup)
         enddo
!
         if (nspin.eq.2) then
            arho(isdw) = abs(rhor(ir,isdw))
            arhotot = abs(rhor(ir,isup)+rhor(ir,isdw))
            do icar = 1, 3
               grad(icar,isdw) = gradr(ir,icar,isdw)
               gradtot(icar) = gradr(ir,icar,isup)+gradr(ir,icar,isdw)
            enddo
            zet = (rhor(ir,isup) - rhor(ir,isdw)) / arhotot
            if (zet.ge. 1.d0) zet =  1.d0
            if (zet.le.-1.d0) zet = -1.d0
         endif
!
         do iss = 1, nspin
            agrad(iss) = sqrt( grad(1,iss)*grad(1,iss) +                &
     &                         grad(2,iss)*grad(2,iss) +                &
     &                         grad(3,iss)*grad(3,iss) )
            agradtot = sqrt( gradtot(1)*gradtot(1) +                    &
     &                       gradtot(2)*gradtot(2) +                    &
     &                       gradtot(3)*gradtot(3) )
         enddo
!
!     _________________________________________________________________
!     First it calculates the energy density excrho
!     exrho:  exchange term
!     ecrho:  correlation term
!
         if ( nspin.eq.2 ) then
            scl = 2.d0
            scl1 = 0.5d0
         else
            scl = 1.d0
            scl1 = 1.d0
         endif
         do iss = 1, nspin
            if ( arho(iss).gt.lim1) then
               call exchpbe( scl*arho(iss), scl*agrad(iss),             &
     &                       exrho(iss),dexdrho(iss),dexdg(iss))
               excrho = excrho + scl1*exrho(iss)
            else
               dexdrho(iss) = 0.d0
               dexdg(iss) = 0.d0
            endif
         enddo
         if ( arhotot.gt.lim1) then
            call ecorpbe( arhotot, agradtot, zet, ecrho,                &
     &                    decdrho(1), decdrho(2), decdg, nspin )
            excrho = excrho + ecrho
         else
            decdrho(isup) = 0.d0
            decdrho(isdw) = 0.d0
            decdg = 0.d0
         endif
!     _________________________________________________________________
!     Now it calculates the potential and writes it in rhor
!     it uses the following variables:
!     dexdrho = d ( ex*rho ) / d (rho)
!     decdrho = d ( ec*rho ) / d (rho)
!     dexdg   = (d ( ex*rho ) / d (grad(rho)_i)) * agrad / grad_i
!     decdg   = (d ( ec*rho ) / d (grad(rho)_i)) * agrad / grad_i
!     gradr  here is used as a working array
!
!     _________________________________________________________________
!    first part of the xc-potential : D(rho*exc)/D(rho)
!
         do iss = 1, nspin
            rhor(ir,iss) = dexdrho(iss) + decdrho(iss)
         enddo
!
!     gradr = D(rho*exc)/D(|grad rho|) * (grad rho) / |grad rho|
!
         do iss = 1, nspin
            do icar = 1,3
               wrkup =0.d0
               wrkdw =0.d0
               if (agrad(iss).gt.lim2)                                  &
     &                    wrkup = dexdg(iss)*grad(icar,iss)/agrad(iss)
               if (agradtot.gt.lim2)                                    &
     &                    wrkdw = decdg*gradtot(icar)/agradtot
               gradr(ir,icar,iss) = wrkup + wrkdw
            enddo
         enddo
!
      end do
!
      return
      end subroutine ggapbe
!
!______________________________________________________________________
      subroutine exchpbe(rho,agrad,ex,dexdrho,dexdg)
!     _________________________________________________________________
!
! Perdew-Burke-Ernzerhof gga, Exchange term:
! Calculates the exchange energy density and the two functional derivative
! that will be used to calculate the potential
!
      USE kinds, ONLY: DP
      implicit none
! input
! input rho:     charge density
! input agrad:   abs(grad rho)
      real(DP) rho, agrad
! ouput
! output ex: Ex[rho,grad_rho] = \int ex dr
! output dexdrho: d ex / d rho
! output dexdg:   d ex / d grad_rho(i) = dexdg*grad_rho(i)/abs(grad_rho)
      real(DP) ex, dexdrho, dexdg
! local
      real(DP) thrd, thrd4, pi32td, ax, al, um, uk, ul
      parameter(thrd=.33333333333333333333d0,thrd4=4.d0/3.d0)
      parameter(pi32td=3.09366772628014d0) ! pi32td=(3.d0*pi*pi)**0.333d0
      parameter(al=0.161620459673995d0)    ! al=1.0/(2.0*(pi32)**0.333d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
!
      real(DP) rhothrd, exunif, dexunif, kf, s, s2, p0, fxpbe, fs
!----------------------------------------------------------------------
! construct LDA exchange energy density
!
      rhothrd = rho**thrd
      dexunif = ax*rhothrd
      exunif  = rho*dexunif
!----------------------------------------------------------------------
! construct PBE enhancement factor
!
      kf = pi32td*rhothrd
      s = agrad/(2.d0*kf*rho)
      s2 = s*s
      p0 = 1.d0 + ul*s2
      fxpbe = 1.d0 + uk - uk/p0
      ex = exunif*fxpbe
!----------------------------------------------------------------------
! now calculates the potential terms
!
!  fs=(1/s)*d fxPBE/ ds
!
      fs=2.d0*uk*ul/(p0*p0)
      dexdrho = dexunif*thrd4*(fxpbe-s2*fs)
      dexdg = ax*al*s*fs
!
      return
      end subroutine exchpbe

!----------------------------------------------------------------------
      subroutine ecorpbe(rho,agrad,zet,ectot,decup,decdn,decdg,nspin)
!     -----------------------------------------------------------------
!
!  Adapted from the Official PBE correlation code. K. Burke, May 14, 1996.
!
!   input: rho   = rho_up + rho_down; total  charge density
!   input: agrad = abs( grad(rho) )
!   input: zet   = (rho_up-rho_down)/rho
!   input: nspin
!  output: ectot = ec*rho       ---correlation energy density---
!  output: decup = d ( ec*rho ) / d (rho_up)
!  output: decdn = d ( ec*rho ) / d (rho_down)
!  output: decdg = (d ( ec*rho ) / d (grad(rho)_i)) * agrad / grad_i
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
!     {\sl Generalized gradient approximation made simple}, sub.
!     to Phys. Rev.Lett. May 1996.
! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
!     construction of a generalized gradient approximation:  The PW91
!     density functional}, submitted to Phys. Rev. B, Feb. 1996.
! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE kinds, ONLY: DP
      USE constants, ONLY: pi
      implicit none
      real(DP) rho, agrad, zet, ectot, decup, decdn, decdg
      integer nspin
      real(DP) pi32, alpha, thrd, thrdm, thrd2, sixthm, thrd4,  &
          gam, fzz, gamma, bet, delt, eta
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      gam= 2^(4/3)-2
!      fzz=f''(0)= 8/(9*gam)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      parameter(pi32=29.608813203268075856503472999628d0)
      parameter(alpha=1.91915829267751300662482032624669d0)
      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(thrd4=4.d0*thrd)
      parameter(gam=0.5198420997897463295344212145565d0)
      parameter(fzz=8.d0/(9.d0*gam))
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(bet=0.06672455060314922d0,delt=bet/gamma)
      parameter(eta=1.d-12)
      real(DP) g, fk, rs, sk, twoksg, t
      real(DP) rtrs, eu, eurs, ep, eprs, alfm, alfrsm, z4, f, ec
      real(DP) ecrs, fz, eczet, comm, vcup, vcdn, g3, pon, b, b2, t2, t4
      real(DP) q4, q5, h, g4, t6, rsthrd, gz, fac
      real(DP) bg, bec, q8, q9, hb, hrs, hz, ht, pref
!----------------------------------------------------------------------
      if (nspin.eq.1) then
         g=1.d0
      else
         g=((1.d0+zet)**thrd2+(1.d0-zet)**thrd2)*0.5d0
      endif
      fk=(pi32*rho)**thrd
      rs=alpha/fk
      sk=sqrt(4.d0*fk/pi)
      twoksg=2.d0*sk*g
      t=agrad/(twoksg*rho)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! eu=unpolarized LSD correlation energy
! eurs=deu/drs
! ep=fully polarized LSD correlation energy
! eprs=dep/drs
! alfm=-spin stiffness, [c](3).
! alfrsm=-dalpha/drs
! f=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      call gcor2(0.0310907d0,0.21370d0,7.5957d0,3.5876d0,1.6382d0,      &
     &    0.49294d0,rtrs,eu,eurs)
      if (nspin.eq.2) then
         call gcor2(0.01554535d0,0.20548d0,14.1189d0,6.1977d0,3.3662d0, &
     &       0.62517d0,rtrs,ep,eprs)
         call gcor2(0.0168869d0,0.11125d0,10.357d0,3.6231d0,0.88026d0,  &
     &       0.49671d0,rtrs,alfm,alfrsm)
         z4 = zet**4
         f=((1.d0+zet)**thrd4+(1.d0-zet)**thrd4-2.d0)/gam
         ec = eu*(1.d0-f*z4)+ep*f*z4-alfm*f*(1.d0-z4)/fzz
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ecrs = dec/drs [c](A2)
! eczet=dec/dzeta [c](A3)
! fz = df/dzeta [c](A4)
         ecrs = eurs*(1.d0-f*z4)+eprs*f*z4-alfrsm*f*(1.d0-z4)/fzz
         fz = thrd4*((1.d0+zet)**thrd-(1.d0-zet)**thrd)/gam
         eczet = 4.d0*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu       &
     &           -(1.d0-z4)*alfm/fzz)
         comm = ec -rs*ecrs/3.d0-zet*eczet
         vcup = comm + eczet
         vcdn = comm - eczet
      else
         ecrs = eurs
         ec = eu
         vcup = ec -rs*ecrs/3.d0
      endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! PBE correlation energy
! g=phi(zeta), given after [a](3)
! delt=bet/gamma
! b=a of [a](8)
!      g=((1.d0+zet)**thrd2+(1.d0-zet)**thrd2)/2.d0
      g3 = g**3
      pon=-ec/(g3*gamma)
      b = delt/(dexp(pon)-1.d0)
      b2 = b*b
      t2 = t*t
      t4 = t2*t2
      q4 = 1.d0+b*t2
      q5 = 1.d0+b*t2+b2*t4
      h = g3*(bet/delt)*dlog(1.d0+delt*Q4*t2/Q5)
      ectot = rho*(ec + h)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! energy done. Now the potential, using appendix e of [b].
      t6 = t4*t2
      rsthrd = rs/3.d0
      fac = delt/b+1.d0
      bec = b2*fac/(bet*g3)
      q8 = q5*q5+delt*q4*q5*t2
      q9 = 1.d0+2.d0*b*t2
      hb = -bet*g3*b*t6*(2.d0+b*t2)/q8
      hrs = -rsthrd*hb*bec*ecrs
      ht = 2.d0*bet*g3*q9/q8
      comm = h+hrs-7.d0*t2*ht/6.d0
      if (nspin.eq.2) then
         g4 = g3*g
         bg = -3.d0*b2*ec*fac/(bet*g4)
         gz=(((1.d0+zet)**2+eta)**sixthm-                               &
     &   ((1.d0-zet)**2+eta)**sixthm)/3.d0
         hz = 3.d0*gz*h/g + hb*(bg*gz+bec*eczet)
         pref = hz-gz*t2*ht/g
         decup = vcup + comm + pref*(  1.d0 - zet)
         decdn = vcdn + comm + pref*( -1.d0 - zet)
      else
         decup = vcup + comm
      endif
      decdg = t*ht/twoksg
!
      return
      end subroutine ecorpbe
!______________________________________________________________________
      subroutine gcor2(a,a1,b1,b2,b3,b4,rtrs,gg,ggrs)
!     _________________________________________________________________
! slimmed down version of GCOR used in PW91 routines, to interpolate
! LSD correlation energy, as given by (10) of
! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! K. Burke, May 11, 1996.
!
      USE kinds, ONLY : DP
      implicit none
      real(DP) a, a1, b1, b2, b3, b4, rtrs, gg, ggrs
      real(DP) q0, q1, q2, q3
!
      q0 = -2.d0*a*(1.d0+a1*rtrs*rtrs)
      q1 = 2.d0*a*rtrs*(b1+rtrs*(b2+rtrs*(b3+b4*rtrs)))
      q2 = dlog(1.d0+1.d0/q1)
      gg = q0*q2
      q3 = a*(b1/rtrs+2.d0*b2+rtrs*(3.d0*b3+4.d0*b4*rtrs))
      ggrs = -2.d0*a*a1*q2-q0*q3/(q1*(1.d0+q1))
!
      return
      end subroutine gcor2
!
!______________________________________________________________________
      subroutine ggapw(nnr,nspin,gradr,rhor,exc)
!     _________________________________________________________________
!     perdew-wang gga (PW91)
!
      USE kinds, ONLY: DP
      use constants, only: pi, fpi
!
      implicit none
! input
      integer nspin, nnr
      real(DP) gradr(nnr,3,nspin), rhor(nnr,nspin)
! output
      real(DP) exc
! local
      integer isup, isdw, ir
      real(DP) rhoup, rhodw, roe, aroe, rs, zeta
      real(DP) grxu, gryu, grzu, grhou, grxd, gryd, grzd, grhod, grho
      real(DP) ex, ec,vc, sc, v1x, v2x, v1c, v2c
      real(DP) ecrs, eczeta
      real(DP) exup, vcup, v1xup, v2xup, v1cup
      real(DP) exdw, vcdw, v1xdw, v2xdw, v1cdw
      real(DP), parameter:: pi34 = 0.75d0/pi, third = 1.d0/3.d0, &
                small = 1.d-10
!
!     _________________________________________________________________
!     main loop
!
      isup=1
      isdw=2
      exc=0.0d0
      do ir=1,nnr
         rhoup=rhor(ir,isup)
         if(nspin.eq.2) then
            rhodw=rhor(ir,isdw)
         else
            rhodw=0.0d0
         end if
         roe=rhoup+rhodw
         aroe=abs(roe)
         if (aroe.lt.small) then
            rhor(ir,isup)   =0.0d0
            gradr(ir,1,isup)=0.0d0
            gradr(ir,2,isup)=0.0d0
            gradr(ir,3,isup)=0.0d0
            if(nspin.eq.2) then
               rhor(ir,isdw)   =0.0d0
               gradr(ir,1,isdw)=0.0d0
               gradr(ir,2,isdw)=0.0d0
               gradr(ir,3,isdw)=0.0d0
            end if
            go to 100
         end if
         grxu =gradr(ir,1,isup)
         gryu =gradr(ir,2,isup)
         grzu =gradr(ir,3,isup)
         grhou=sqrt(grxu**2+gryu**2+grzu**2)
         if(nspin.eq.2) then
            grxd =gradr(ir,1,isdw)
            gryd =gradr(ir,2,isdw)
            grzd =gradr(ir,3,isdw)
            grhod=sqrt(grxd**2+gryd**2+grzd**2)
         else
            grxd =0.0d0
            gryd =0.0d0
            grzd =0.0d0
            grhod=0.0d0
         endif
         grho=sqrt((grxu+grxd)**2+(gryu+gryd)**2+(grzu+grzd)**2)
!
         rs=(pi34/aroe)**third
         if (nspin.eq.1) then
            call exchpw91(aroe,grho,ex,v1x,v2x)
            call pwlda(rs,ec,vc,ecrs)
            call corpw91ns(rs,grho,ec,ecrs,sc,v1c,v2c)
            exc = exc + roe*(ex+ec) + sc
            rhor(ir,isup) = vc + v1x + v1c
!
!     gradr = D(rho*exc)/D(|grad rho|) * (grad rho) / |grad rho|
!
            gradr(ir,1,isup)=grxu*(v2x+v2c)
            gradr(ir,2,isup)=gryu*(v2x+v2c)
            gradr(ir,3,isup)=grzu*(v2x+v2c)
         else
            zeta=(rhoup-rhodw)/aroe        
            zeta=min(zeta, 1.d0)
            zeta=max(zeta,-1.d0)
            call exchpw91(2.d0*abs(rhoup),2.0d0*grhou,exup,v1xup,v2xup)
            call exchpw91(2.d0*abs(rhodw),2.0d0*grhod,exdw,v1xdw,v2xdw)
            call pwlsd(rs,zeta,ec,vcup,vcdw,ecrs,eczeta)
            call corpw91(rs,zeta,grho,ec,ecrs,eczeta,sc,v1cup,v1cdw,v2c)
            rhor(ir,isup) = vcup + v1xup + v1cup
            rhor(ir,isdw) = vcdw + v1xdw + v1cdw
            exc = exc+roe*(0.5d0*((1.d0+zeta)*exup+(1.d0-zeta)*exdw)+ec) &
                 + sc
!
!     gradr = D(rho*exc)/D(|grad rho|) * (grad rho) / |grad rho|
!
            gradr(ir,1,isup)=grxu*(2.0d0*v2xup+v2c)+grxd*v2c
            gradr(ir,2,isup)=gryu*(2.0d0*v2xup+v2c)+gryd*v2c
            gradr(ir,3,isup)=grzu*(2.0d0*v2xup+v2c)+grzd*v2c
            gradr(ir,1,isdw)=grxd*(2.0d0*v2xdw+v2c)+grxu*v2c
            gradr(ir,2,isdw)=gryd*(2.0d0*v2xdw+v2c)+gryu*v2c
            gradr(ir,3,isdw)=grzd*(2.0d0*v2xdw+v2c)+grzu*v2c
         end if
 100     continue
      end do
!
      return
      end subroutine ggapw
!
!----------------------------------------------------------------------
      subroutine exchpw91(rho,grho,ex,v1x,v2x)
!----------------------------------------------------------------------
!
!  PW91 exchange for a spin-unpolarized electronic system
!  Modified from the "official" PBE code of Perdew, Burke et al.
!  input rho   : density
!  input grho:  abs(grad rho)
!  output:  exchange energy per electron (ex) and potentials
!          v1x = d(rho*exc)/drho
!          v2x = d(rho*exc)/d|grho| * (1/|grho|)
!
      USE kinds, ONLY : DP
      USE constants, ONLY : pi
      implicit none
!  input
      real(DP) rho, grho
!  output
      real(DP) ex, v1x, v2x
! local
      real(DP) ex0, kf, s, s2, s4, f, fs, p0,p1,p2,p3,p4,p5,p6,p7
! parameters
      real(DP) a1, a2, a3, a4, a, b1, bx, pi34, thrd, thrd4
      parameter(a1=0.19645d0,a2=0.27430d0,a=7.7956d0,a4=100.d0)
! for becke exchange, set a3=b1=0
      parameter(a3=0.15084d0,b1=0.004d0)
! pi34=3/(4pi) ,  bx=(3pi^2)^(1/3)
      parameter(pi34=0.75d0/pi, bx=3.093667726d0, thrd=0.333333333333d0, &
                thrd4=4.d0*thrd)
!
      if (rho.lt.1.d-10) then
         ex =0.0d0
         v1x=0.0d0
         v2x=0.0d0
      end if
!
!  kf=k_Fermi, ex0=Slater exchange energy
!
      kf = bx*(rho**thrd)
      ex0=-pi34*kf
      if (grho.lt.1.d-10) then
         ex =ex0
         v1x=ex0*thrd4
         v2x=0.0d0
      end if
      s  = grho/(2.d0*kf*rho)
      s2 = s*s
      s4 = s2*s2
      p0 = 1.d0/sqrt(1.d0+a*a*s2)
      p1 = log(a*s+1.d0/p0)
      p2 = exp(-a4*s2)
      p3 = 1.d0/(1.d0+a1*s*p1+b1*s4)
      p4 = 1.d0+a1*s*p1+(a2-a3*p2)*s2
!  f is the enhancement factor
      f = p3*p4
      ex = ex0*f
!  energy done. now the potential:
      p5 = b1*s2-(a2-a3*p2)
      p6 = a1*s*(p1+a*s*p0)
      p7 = 2.d0*(a2-a3*p2)+2.d0*a3*a4*s2*p2-4.d0*b1*s2*f
! fs = (1/s) dF(s)/ds
      fs = p3*(p3*p5*p6+p7)
      v1x = ex0*thrd4*(f-s2*fs)
      v2x = 0.5d0*ex0/kf*s*fs/grho
!
      return
      end subroutine exchpw91
!
!----------------------------------------------------------------------
      subroutine corpw91ns(rs,grho,ec,ecrs,h,v1c,v2c)
!----------------------------------------------------------------------
!
!  PW91 correlation (gradient correction term) - no spin case
!  Modified from the "official" PBE code of Perdew, Burke et al.
!
!  input rs:   seitz radius
!  input zeta:  relative spin polarization
!  input grho: abs(grad rho)
!  input ec:   Perdew-Wang correlation energy
!  input ecrs:  d(rho*ec)/d r_s

!  output h  :  nonlocal part of correlation energy per electron
!  output v1c:  nonlocal parts of correlation potential
!         v1c = d(rho*exc)/drho
!         v2c = d(rho*exc)/d|grho|*(1/|grho|)
!
      USE kinds, ONLY : DP
      USE constants, ONLY : pi
      implicit none
! input
      real(DP) rs, grho, ec, ecrs
! output
      real(DP) h, v1c, v2c
! local
      real(DP) rho, t, ks,  bet, delt, pon, b, b2, t2, t4, t6
      real(DP) q4, q5, q6, q7, q8, q9, r0, r1, r2, r3, r4, rs2, rs3
      real(DP) ccrs, rsthrd, fac, bec, coeff, cc
      real(DP) h0, h0b, h0rs, h0t, h1, h1t, h1rs, hrs, ht
! parameters
      real(DP) nu, cc0, cx, alf, c1, c2, c3, c4, c5, c6, a4, ax, pi34
      parameter(nu=15.75592d0,cc0=0.004235d0,cx=-0.001667212d0)
      parameter(c1=0.002568d0,c2=0.023266d0,c3=7.389d-6,c4=8.723d0)
      parameter(c5=0.472d0,c6=7.389d-2,a4=100.d0, alf=0.09d0)
! ax=(4*1.9191583/pi)^(1/2), where k_F=1.9191583/r_s, k_s=boh*r_s^(1/2)
      parameter(ax=1.5631853d0, pi34 = 0.75d0/pi)
!
!
      rs2 = rs*rs
      rs3 = rs2*rs
      rho=pi34/rs3
!  k_s=(4k_F/pi)^(1/2)
      ks=ax/sqrt(rs)
!  t=abs(grad rho)/(rho*2.*ks)
      t=grho/(2.d0*rho*ks)
      bet = nu*cc0
      delt = 2.d0*alf/bet
      pon = -delt*ec/bet
      b = delt/(exp(pon)-1.d0)
      b2 = b*b
      t2 = t*t
      t4 = t2*t2
      t6 = t4*t2
      q4 = 1.d0+b*t2
      q5 = 1.d0+b*t2+b2*t4
      q6 = c1+c2*rs+c3*rs2
      q7 = 1.d0+c4*rs+c5*rs2+c6*rs3
      cc = -cx + q6/q7
      r0 = 0.663436444d0*rs
      r1 = a4*r0
      coeff = cc-cc0-3.d0*cx/7.d0
      r2 = nu*coeff
      r3 = exp(-r1*t2)
      h0 = (bet/delt)*log(1.d0+delt*q4*t2/q5)
      h1 = r3*r2*t2
      h = (h0+h1)*rho
!  energy done. now the potential:
      ccrs = (c2+2.d0*c3*rs)/q7 - q6*(c4+2.d0*c5*rs+3.d0*c6*rs2)/q7**2
      rsthrd = rs/3.d0
      r4 = rsthrd*ccrs/coeff
      fac = delt/b+1.d0
      bec = b2*fac/bet
      q8 = q5*q5+delt*q4*q5*t2
      q9 = 1.d0+2.d0*b*t2
      h0b = -bet*b*t6*(2.d0+b*t2)/q8
      h0rs = -rsthrd*h0b*bec*ecrs
      h0t = 2.d0*bet*q9/q8
      h1rs = r3*r2*t2*(-r4+r1*t2/3.d0)
      h1t = 2.d0*r3*r2*(1.d0-r1*t2)
      hrs = h0rs+h1rs
      ht = h0t+h1t
      v1c = h0+h1+hrs-7.d0*t2*ht/6.d0
      v2c = t*ht/(2.d0*ks*grho)
!
      return
      end subroutine corpw91ns
!
!----------------------------------------------------------------------
      subroutine corpw91(rs,zeta,grho,ec,ecrs,eczeta,h,v1cup,v1cdn,v2c)
!----------------------------------------------------------------------
!
!  PW91 correlation (gradient correction term)
!  Modified from the "official" PBE code of Perdew, Burke et al.
!
!  input rs:   seitz radius
!  input zeta:  relative spin polarization
!  input grho: abs(grad rho)
!  input ec:   Perdew-Wang correlation energy
!  input ecrs:  d(rho*ec)/d r_s ?
!  input eczeta: d(rho*ec)/d zeta ?

!  output h: nonlocal part of correlation energy per electron
!  output v1cup,v1cdn:  nonlocal parts of correlation potentials
!         v1c** = d(rho*exc)/drho               (up and down components)
!         v2c   = d(rho*exc)/d|grho|*(1/|grho|) (same for up and down)
!
      USE kinds, ONLY : DP
      USE constants, ONLY : pi
      implicit none
! input
      real(DP) rs, zeta, grho, ec, ecrs, eczeta
! output
      real(DP) h, v1cup, v1cdn, v2c
! local
      real(DP) rho, g, t, ks, gz, bet, delt, g3, g4, pon, b, b2, t2, t4, t6
      real(DP) q4, q5, q6, q7, q8, q9, r0, r1, r2, r3, r4, rs2, rs3
      real(DP) ccrs, rsthrd, fac, bg, bec, coeff, cc
      real(DP) h0, h0b, h0rs, h0z, h0t, h1, h1t, h1rs, h1z
      real(DP) hz, hrs, ht, comm, pref
! parameters
      real(DP) nu, cc0, cx, alf, c1, c2, c3, c4, c5, c6, a4
      real(DP) thrdm, thrd2, ax, eta, pi34
      parameter(nu=15.75592d0,cc0=0.004235d0,cx=-0.001667212d0)
      parameter(c1=0.002568d0,c2=0.023266d0,c3=7.389d-6,c4=8.723d0)
      parameter(c5=0.472d0,c6=7.389d-2,a4=100.d0, alf=0.09d0)
      parameter(thrdm=-0.333333333333d0,thrd2=0.666666666667d0)
! ax=(4*1.9191583/pi)^(1/2), where k_F=1.9191583/r_s, k_s=boh*r_s^(1/2)
      parameter(ax=1.5631853d0, eta=1.d-12, pi34 = 0.75d0/pi )
!
!
      if (grho.lt.1.d-10) then
         h=0.0d0
         v1cup=0.0d0
         v1cdn=0.0d0
         v2c=0.0d0
      end if
      rs2 = rs*rs
      rs3 = rs2*rs
      rho=pi34/rs3
      g=((1.d0+zeta)**thrd2+(1.d0-zeta)**thrd2)/2.d0
!  k_s=(4k_F/pi)^(1/2)
      ks=ax/sqrt(rs)
!  t=abs(grad rho)/(rho*2.*ks*g)
      t=grho/(2.d0*rho*g*ks)
      bet = nu*cc0
      delt = 2.d0*alf/bet
      g3 = g**3
      g4 = g3*g
      pon = -delt*ec/(g3*bet)
      b = delt/(exp(pon)-1.d0)
      b2 = b*b
      t2 = t*t
      t4 = t2*t2
      t6 = t4*t2
      q4 = 1.d0+b*t2
      q5 = 1.d0+b*t2+b2*t4
      q6 = c1+c2*rs+c3*rs2
      q7 = 1.d0+c4*rs+c5*rs2+c6*rs3
      cc = -cx + q6/q7
      r0 = 0.663436444d0*rs
      r1 = a4*r0*g4
      coeff = cc-cc0-3.d0*cx/7.d0
      r2 = nu*coeff*g3
      r3 = dexp(-r1*t2)
      h0 = g3*(bet/delt)*log(1.d0+delt*q4*t2/q5)
      h1 = r3*r2*t2
      h = (h0+h1)*rho
!  energy done. now the potential:
      ccrs = (c2+2.d0*c3*rs)/q7 - q6*(c4+2.d0*c5*rs+3.d0*c6*rs2)/q7**2
      rsthrd = rs/3.d0
      r4 = rsthrd*ccrs/coeff
!  eta is a small quantity that avoids trouble if zeta=+1 or -1
      gz = ((1.d0+zeta+eta)**thrdm - (1.d0-zeta+eta)**thrdm)/3.d0
      fac = delt/b+1.d0
      bg = -3.d0*b2*ec*fac/(bet*g4)
      bec = b2*fac/(bet*g3)
      q8 = q5*q5+delt*q4*q5*t2
      q9 = 1.d0+2.d0*b*t2
      h0b = -bet*g3*b*t6*(2.d0+b*t2)/q8
      h0rs = -rsthrd*h0b*bec*ecrs
      h0z = 3.d0*gz*h0/g + h0b*(bg*gz+bec*eczeta)
      h0t = 2.d0*bet*g3*q9/q8
      h1rs = r3*r2*t2*(-r4+r1*t2/3.d0)
      h1z = gz*r3*r2*t2*(3.d0-4.d0*r1*t2)/g
      h1t = 2.d0*r3*r2*(1.d0-r1*t2)
      hrs = h0rs+h1rs
      ht = h0t+h1t
      hz = h0z+h1z
      comm = h0+h1+hrs-7.d0*t2*ht/6.d0
      pref = hz-gz*t2*ht/g
      comm = comm-pref*zeta
      v1cup = comm + pref 
      v1cdn = comm - pref
      v2c   = t*ht/(2.d0*ks*g*grho)
!
      return
      end subroutine corpw91
!----------------------------------------------------------------------
      subroutine pwlda(rs,ec,vc,ecrs)
!----------------------------------------------------------------------
!
!  uniform-gas, spin-unpolarised correlation of perdew and wang 1991
!  input:  rs   seitz radius
!  output: ec   correlation energy per electron
!          vc   potential
!          ecrs derivatives of ec wrt rs
!
      USE kinds, ONLY : DP
      implicit none
! input
      real(DP) rs
! output
      real(DP) ec, vc, ecrs
!  local
      real(DP) q0, rs12, q1, q2, q3
! parameters
      real(DP) a, a1, b1, b2, b3, b4
      parameter(a =0.0310907d0, a1=0.21370d0, b1=7.5957d0,              &
                b2=3.5876d0,    b3=1.6382d0,  b4=0.49294d0)
!
      q0 = -2.d0*a*(1.d0+a1*rs)
      rs12 = sqrt(rs)
      q1 = 2.d0*a*rs12*(b1+rs12*(b2+rs12*(b3+b4*rs12)))
      q2 = log(1.d0+1.d0/q1)
      ec = q0*q2
      q3 = a*(b1/rs12+2.d0*b2+3.d0*b3*rs12+2.d0*b4*2.d0*rs)
      ecrs = -2.d0*a*a1*q2-q0*q3/(q1**2+q1)
      vc = ec - rs*ecrs/3.d0
!
      return
      end subroutine pwlda
!----------------------------------------------------------------------
      subroutine pwlsd(rs,zeta,ec,vcup,vcdn,ecrs,eczeta)
!----------------------------------------------------------------------
!
!  uniform-gas correlation of perdew and wang 1991
!  Modified from the "official" PBE code of Perdew, Burke et al.
!  input: seitz radius (rs), relative spin polarization (zeta)
!  output: correlation energy per electron (ec)
!          up- and down-spin potentials (vcup,vcdn)
!          derivatives of ec wrt rs (ecrs) & zeta (eczeta)
!
      USE kinds, ONLY : DP
      implicit none
! input
      real(DP) rs, zeta
! output
      real(DP) ec, vcup, vcdn, ecrs, eczeta
! local
      real(DP) f, eu, ep, eurs, eprs, alfm, alfrsm, z4, fz, comm
      real(DP) rs12, q0, q1, q2, q3
! parameters
      real(DP) gam, fzz, thrd, thrd4
      parameter(gam=0.5198421d0,fzz=1.709921d0)
      parameter(thrd=0.333333333333d0,thrd4=1.333333333333d0)
!
      real(DP) au, au1, bu1, bu2, bu3, bu4
      parameter(au =0.0310907d0, au1=0.21370d0, bu1=7.5957d0,           &
                bu2=3.5876d0,    bu3=1.6382d0,  bu4=0.49294d0)
      real(DP) ap, ap1, bp1, bp2, bp3, bp4
      parameter(ap =0.01554535d0,ap1=0.20548d0, bp1=14.1189d0,          &
                bp2=6.1977d0,    bp3=3.3662d0,  bp4=0.62517d0 )
      real(DP) am, am1, bm1, bm2, bm3, bm4
      parameter(am =0.0168869d0, am1=0.11125d0, bm1=10.357d0,           &
                bm2=3.6231d0,    bm3=0.88026d0, bm4=0.49671d0 )
!
      rs12 = sqrt(rs)
!
      q0 = -2.d0*au*(1.d0+au1*rs)
      q1 = 2.d0*au*rs12*(bu1+rs12*(bu2+rs12*(bu3+bu4*rs12)))
      q2 = log(1.d0+1.d0/q1)
      eu = q0*q2
      q3 = au*(bu1/rs12+2.d0*bu2+3.d0*bu3*rs12+2.d0*bu4*2.d0*rs)
      eurs = -2.d0*au*au1*q2-q0*q3/(q1**2+q1)
!
      q0 = -2.d0*ap*(1.d0+ap1*rs)
      q1 = 2.d0*ap*rs12*(bp1+rs12*(bp2+rs12*(bp3+bp4*rs12)))
      q2 = log(1.d0+1.d0/q1)
      ep = q0*q2
      q3 = ap*(bp1/rs12+2.d0*bp2+3.d0*bp3*rs12+2.d0*bp4*2.d0*rs)
      eprs = -2.d0*ap*ap1*q2-q0*q3/(q1**2+q1)
!
      q0 = -2.d0*am*(1.d0+am1*rs)
      q1 = 2.d0*am*rs12*(bm1+rs12*(bm2+rs12*(bm3+bm4*rs12)))
      q2 = log(1.d0+1.d0/q1)
!  alfm is minus the spin stiffness alfc
      alfm=q0*q2
      q3 = am*(bm1/rs12+2.d0*bm2+3.d0*bm3*rs12+2.d0*bm4*2.d0*rs)
      alfrsm=-2.d0*am*am1*q2-q0*q3/(q1**2+q1)
!
      f = ((1.d0+zeta)**thrd4+(1.d0-zeta)**thrd4-2.d0)/gam
      z4 = zeta**4
      ec = eu*(1.d0-f*z4)+ep*f*z4-alfm*f*(1.d0-z4)/fzz
!  energy done. now the potential:
      ecrs = eurs*(1.d0-f*z4)+eprs*f*z4-alfrsm*f*(1.d0-z4)/fzz
      fz = thrd4*((1.d0+zeta)**thrd-(1.d0-zeta)**thrd)/gam
      eczeta = 4.d0*(zeta**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu        &
     &        -(1.d0-z4)*alfm/fzz)
      comm = ec -rs*ecrs/3.d0-zeta*eczeta
      vcup = comm + eczeta
      vcdn = comm - eczeta
!
      return
      end subroutine pwlsd
!
!______________________________________________________________________
      subroutine ggapwold(nnr,nspin,gradr,rhor,exc)
!     _________________________________________________________________
!     perdew-wang gga 
!     as given in y-m juan & e kaxiras, prb 48, 14944 (1993) 
!     method by ja white & dm bird, prb 50, 4954 (1994) 
!     non-spin polarized case only
!     _________________________________________________________________
!     by alfredo pasquarello 22/09/1994
!
      USE kinds, ONLY: DP
      use constants, only: pi, fpi
!
      implicit none
!
      integer nspin, nnr
      real(DP) gradr(nnr,3), rhor(nnr), exc
!
      real(DP) bb1, bb2, bb3, bb4, bb5, alfa, beta, cc0, cc1, delt, &
           c1, c2, c3, c4, c5, c6, c7, a, alfa1, bt1, bt2, bt3, bt4
      parameter(bb1=0.19645d0,bb2=0.27430d0,bb3=-0.15084d0,bb4=0.004d0,         &
       bb5=7.7956d0,alfa=0.09d0,beta=0.0667263212d0,cc0=15.75592d0,             &
       cc1=0.003521d0,c1=0.001667d0,c2=0.002568d0,c3=0.023266d0,c4=7.389d-6,    &
       c5=8.723d0,c6=0.472d0,c7=7.389d-2,a=0.0621814d0,alfa1=0.2137d0,          &
       bt1=7.5957d0,bt2=3.5876d0,bt3=1.6382d0,bt4=0.49294d0,delt=1.0d-12) 
      real(DP) x13, x43, x76, pi2, ax, pider1, pider2, pider3,      &
           abder1, abder2, abder3
      integer isign, ir
      real(DP)                                                      &
           aexp, abig, abig2, agr, aroe, byagr, ccr, ccrnum, ccrden,    &
           dfxd, dfxdg, dys, dfs, dh1ds, dh1dg, dh1d, dh1dt, dexcdg,    &
           dexcd, dh1drs, dh0da, dadec, decdrs, decd, dh0dg, dcdrs,     &
           dh0d, dh0dt, eclog, ecr, ecden, fx, fxnum, fxden, fxexp,     &
           gkf, grx, gry, grz, h0, h1, h0den, h0arg, h0num,             &
           roeth, roe, rs, rs12, rs2, rs3, rs32, s, sd, s2, s3, s4,     &
           sysl, t, td, t2, t3, t4, xchge, ys, ysl, ysr
!
!
      if (nspin.ne.1) call errore('ggapw','spin not implemented',nspin)
!
      x13=1.0d0/3.0d0
      x43=4.0d0/3.0d0
      x76=7.0d0/6.0d0
!     _________________________________________________________________
!     derived parameters from pi
!
      pi2=pi*pi
      ax=-0.75d0*(3.0d0/pi)**x13
      pider1=(0.75d0/pi)**x13
      pider2=(3.0d0*pi2)**x13
      pider3=(3.0d0*pi2/16.0d0)**x13
!     _________________________________________________________________
!     derived parameters from alfa and beta 
!
      abder1=beta*beta/(2.0d0*alfa)
      abder2=1.0d0/abder1
      abder3=2.0d0*alfa/beta
!     _________________________________________________________________
!     main loop
!
      do ir=1,nnr
         roe=rhor(ir)
         if(roe.eq.0.0) goto 100
         aroe=abs(roe)
         grx=gradr(ir,1)
         gry=gradr(ir,2)
         grz=gradr(ir,3)
         agr=sqrt(grx*grx+gry*gry+grz*grz)
         roeth=aroe**x13 
         rs= pider1/roeth
         gkf=pider2*roeth
         sd=1.0d0/(2.0d0*gkf*aroe)
         s=agr*sd
         s2=s*s
         s3=s*s2
         s4=s2*s2
!     _________________________________________________________________
!     exchange 
!
         ysr=sqrt(1.0d0+bb5*bb5*s2)
         ys=bb5*s+ysr
         ysl=log(ys)*bb1
         sysl=s*ysl
         fxexp=exp(-100.0d0*s2)
         fxnum=1.0d0+sysl+(bb2+bb3*fxexp)*s2
         fxden=1.0d0/(1.0d0+sysl+bb4*s4)
         fx=fxnum*fxden
         xchge=ax*fx*roeth
!     _________________________________________________________________
!     correlation ecr=ec(rho) 
!
         rs12=sqrt(rs)
         rs32=rs12*rs
         rs2=rs*rs
         rs3=rs*rs2
         ecden=a*(bt1*rs12+bt2*rs+bt3*rs32+bt4*rs2) 
         eclog=log(1.0d0+(1.0d0/ecden))
         ecr=-a*(1.0d0+alfa1*rs)*eclog
!     _________________________________________________________________
!     correlation h0(t,ecr)
!
         td=pider3*sd/rs12
         t=agr*td
         t2=t*t
         t3=t*t2
         t4=t2*t2
         aexp=exp(-abder2*ecr)-1.0d0
         abig=abder3/aexp
         abig2=abig*abig
         h0num=t2+abig*t4
         h0den=1.0d0/(1.0d0+abig*t2+abig2*t4)
         h0arg=1.0d0+abder3*h0num*h0den
         h0=abder1*log(h0arg)
!     _________________________________________________________________
!     correlation h1(t,s,aroe)
!
         ccrnum=c2+c3*rs+c4*rs2
         ccrden=1.0d0/(1.0d0+c5*rs+c6*rs2+c7*rs3)
         ccr=c1+ccrnum*ccrden
         h1=cc0*(ccr-cc1)*t2*fxexp
!     _________________________________________________________________
!     updating of xc-energy
!
         exc=exc+(xchge+ecr+h0+h1)*aroe
!     _________________________________________________________________
!     first part xc-potential from exchange  
!     
         dys=bb5*(1.0d0+bb5*s/ysr)/ys
         dfs=-fxnum*(ysl+bb1*s*dys+4.0d0*bb4*s3)*fxden*fxden              &
     &        +(ysl+bb1*s*dys+2.0d0*s*(bb2+bb3*fxexp)                     &
     &        -200.0d0*s3*bb3*fxexp)*fxden
         dfxd=(ax*roeth*x43)*(fx-dfs*s)
         dfxdg=ax*roeth*dfs*sd
!     _________________________________________________________________
!     first part xc-potential from ecr 
!
         decdrs=-a*alfa1*eclog*rs + a*(1+alfa1*rs)                      &
     &        *a*(0.5d0*bt1*rs12+bt2*rs+1.5d0*bt3*rs32+2.0d0*bt4*rs2)         &
     &        /(ecden*ecden+ecden)
         decd=-x13*decdrs
!     _________________________________________________________________
!     first part xc-potential from h0 
!     
         dh0da=abder1/h0arg*abder3*h0den*                               &
     &        (t4-h0num*h0den*(t2+2.0d0*abig*t4))
         dadec=abder3*abder2*(aexp+1.0d0)/(aexp*aexp)
         dh0d=dh0da*dadec*decd 
         dh0dt=abder1/h0arg*abder3*h0den                                &
     &        *(2.0d0*t+4.0d0*abig*t3-h0num*h0den*(2.0d0*abig*t+4.0d0*abig2*t3))
         dh0d=dh0d-x76*t*dh0dt
         dh0dg=dh0dt*td
!     _________________________________________________________________
!     first part xc-potential from h1 
!
         dcdrs=(c3+2.0d0*c4*rs-ccrnum*ccrden*(c5+2.0d0*c6*rs+3.0d0*c7*rs2))   &
     &        *ccrden
         dh1drs=cc0*t2*fxexp*dcdrs
         dh1d=-x13*rs*dh1drs
         dh1dt=2.0d0*t*cc0*(ccr-cc1)*fxexp
         dh1d=dh1d-x76*t*dh1dt
         dh1ds=-200.0d0*s*cc0*(ccr-cc1)*t2*fxexp
         dh1d=dh1d-x43*s*dh1ds
         dh1dg=dh1dt*td+dh1ds*sd
!     _________________________________________________________________
!     first part of xc-potential: D(rho*exc)/D(rho)
!
         dexcd=dfxd+decd+dh0d+dh1d+ecr+h0+h1
         isign=sign(1.d0,agr-delt)
         byagr=0.5d0*(1+isign)/(agr+(1-isign)*delt)
         rhor(ir)=dexcd
!
!     gradr = D(rho*exc)/D(|grad rho|) * (grad rho) / |grad rho|
!
         dexcdg=(dfxdg+dh0dg+dh1dg)*aroe*byagr
         gradr(ir,1)=gradr(ir,1)*dexcdg
         gradr(ir,2)=gradr(ir,2)*dexcdg
         gradr(ir,3)=gradr(ir,3)*dexcdg
 100     continue
      end do
!
      return
      end subroutine ggapwold

!-----------------------------------------------------------------------
      subroutine dftname_cp (exfact, dft)
!-----------------------------------------------------------------------
!
      implicit none
      integer :: exfact
      character(len=25) dft
!
      if (exfact == 0) then
         dft = 'PZ'
      elseif (exfact == 1) then
         dft = 'BLYP'
      elseif (exfact == 2) then
         dft = 'B88'
      elseif (exfact ==  - 5 .or. exfact == 3) then
         dft = 'BP'
      elseif (exfact ==  - 6 .or. exfact == 4) then
         dft = 'PW91'
      elseif (exfact == 5) then
         dft = 'PBE'
      elseif (exfact ==-1) then
         dft = 'WIG'
      elseif (exfact ==-2) then
         dft = 'HL'
      elseif (exfact ==-3) then
         dft = 'GL'
      elseif (exfact == 6) then
         dft = 'TPSS'
      else
         call errore ('dftname','unknown exch-corr functional',exfact)
      end if

      return
      end subroutine dftname_cp


!-------------------------------------------------------------------------
      subroutine expxc(nnr,nspin,rhor,exc)
!----------------------------------------------------------------------
!
!       ceperley & alder's correlation energy
!       after j.p. perdew & a. zunger prb 23, 5048 (1981)
!
!       rhor contains rho(r) on input, vxc(r) on output
!
      USE kinds, ONLY : DP
      use constants, only: pi, fpi
!
      implicit none
!
      integer nnr, nspin
      real(DP) rhor(nnr,nspin), exc
! local variables
      integer ir, iflg, isup, isdw
      real(DP) roe, aroe, rs, rsl, rsq, ecca, vcca, eccp, vccp,    &
          zeta, onemz, zp, zm, fz, dfzdz, exc1, vxc1, vxc2
! constants
      real(DP) x76, x43, x13
      parameter(x76=7.d0/6.d0, x43=4.d0/3.d0, x13=1.d0/3.d0)
      real(DP) ax
      parameter (ax = -0.916330586d0)
! Perdew and Zunger parameters
      real(DP) ap, bp, cp, dp0, af, bf, cf, df,                      &
          bp1, cp1, dp1, bf1, cf1, df1
      parameter                                                         &
      ( ap=0.03110d0*2.0d0, bp=-0.0480d0*2.0d0, cp=0.0020d0*2.0d0, dp0=-0.0116d0*2.0d0   &
      , af=0.01555d0*2.0d0, bf=-0.0269d0*2.0d0, cf=0.0007d0*2.0d0, df=-0.0048d0*2.0d0   &
      , bp1=bp-ap/3.0d0, cp1=2.0d0*cp/3.0d0, dp1=(2.0d0*dp0-cp)/3.0d0              &
      , bf1=bf-af/3.0d0, cf1=2.0d0*cf/3.0d0, df1=(2.0d0*df-cf)/3.0d0 )
      real(DP) va(2), vb(2), vc(2), vd(2), vbt1(2), vbt2(2)
      real(DP)  a(2), b(2), c(2), d(2), g(2), b1(2), b2(2)
      data va/ap ,af /, vb/bp1,bf1/, vc/cp1,cf1/, vd/dp1,df1/,          &
          vbt1/1.0529d0,1.3981d0/, vbt2/0.3334d0,0.2611d0/
      data a/0.0622d0,0.0311d0/, b/-0.096d0,-0.0538d0/, c/0.0040d0,0.0014d0/,       &
          d/-0.0232d0,-0.0096d0/, b1/1.0529d0,1.3981d0/, b2/0.3334d0,0.2611d0/,    &
          g/-0.2846d0,-0.1686d0/
!
      if (nspin.eq.1) then
!
! iflg=1: paramagnetic (unpolarised) results
!
         iflg=1
         do ir=1,nnr
            roe=rhor(ir,1)
            if(roe.lt.1.0d-30) goto 10
            aroe=abs(roe)
            rs= (3.d0/aroe/fpi)**x13
            if(rs.le.1.d0) then
               rsl=log(rs)
               ecca= a(iflg)*rsl+ b(iflg)+ c(iflg)*rs*rsl+ d(iflg)*rs
               vcca=va(iflg)*rsl+vb(iflg)+vc(iflg)*rs*rsl+vd(iflg)*rs
            else
               rsq=sqrt(rs)
               ecca=g(iflg)/(1.d0+b1(iflg)*rsq+b2(iflg)*rs)
               vcca=ecca*(1.d0+x76*vbt1(iflg)*rsq+x43*vbt2(iflg)*rs)/   &
     &                   (1.d0+    vbt1(iflg)*rsq+    vbt2(iflg)*rs)
            end if
            exc1 = ( ax/rs + ecca )/2.d0
            exc = exc + exc1*roe
            rhor(ir,1)= ( x43*ax/rs + vcca )/2.d0
 10         continue
         end do
      else
         isup=1
         isdw=2
         do ir=1,nnr
            roe=rhor(ir,isup)+rhor(ir,isdw)
            if(roe.lt.1.0d-30) goto 20
            aroe=abs(roe)
            rs= (3.d0/aroe/fpi)**x13
            zeta=abs(rhor(ir,isup)-rhor(ir,isdw))/aroe
            zp = (1.d0+zeta)**x13
            onemz=max(0.d0,1.d0-zeta)
            zm = onemz**x13
            fz= ((1.d0+zeta)*zp + onemz*zm - 2.d0)/                     &
     &           (2.d0**x43 -2.d0)
            dfzdz= x43*(zp - zm)/(2.d0**x43-2.d0)
!
! iflg=1:  paramagnetic (unpolarised) results
! iflg=2: ferromagnetic (  polarised) results
!
            if(rs.le.1.d0) then
               rsl=log(rs)
               ecca= a(1)*rsl+ b(1)+ c(1)*rs*rsl+ d(1)*rs
               vcca=va(1)*rsl+vb(1)+vc(1)*rs*rsl+vd(1)*rs
               eccp= a(2)*rsl+ b(2)+ c(2)*rs*rsl+ d(2)*rs
               vccp=va(2)*rsl+vb(2)+vc(2)*rs*rsl+vd(2)*rs
            else
               rsq=sqrt(rs)
               ecca=g(1)/(1.d0+b1(1)*rsq+b2(1)*rs)
               vcca=ecca*(1.d0+x76*vbt1(1)*rsq+x43*vbt2(1)*rs)/         &
     &                   (1.d0+    vbt1(1)*rsq+    vbt2(1)*rs)
               eccp=g(2)/(1.d0+b1(2)*rsq+b2(2)*rs)
               vccp=eccp*(1.d0+x76*vbt1(2)*rsq+x43*vbt2(2)*rs)/         &
     &                   (1.d0+    vbt1(2)*rsq+    vbt2(2)*rs)
            end if
! exchange part
            exc1 = ax/rs*((1.d0+zeta)*zp+(1.d0-zeta)*zm)/2.d0
            vxc1 = x43*ax/rs*zp
            vxc2 = x43*ax/rs*zm
! correlation part
            vxc1 = vxc1 + vcca + fz*(vccp-vcca)                         &
     &           + dfzdz*(eccp-ecca)*( 1.d0-zeta)
            vxc2 = vxc2 + vcca + fz*(vccp-vcca)                         &
     &           + dfzdz*(eccp-ecca)*(-1.d0-zeta)
            exc  = exc + (exc1 + ecca+fz*(eccp-ecca))*roe/2.d0
            rhor(ir,isup)=vxc1/2.d0
            rhor(ir,isdw)=vxc2/2.d0
 20         continue
         end do
      end if

      return
      end subroutine expxc

      SUBROUTINE wrap_b88( rho, grho, sx, v1x, v2x )
        USE kinds, ONLY: DP
        IMPLICIT NONE
        REAL(DP) ::  rho, grho, sx, v1x, v2x 
        REAL(DP) :: b1 = 0.0042d0
        REAL(DP) :: RHOA,RHOB,GRHOA,GRHOB, V1XA,V2XA,V1XB,V2XB
        rhoa = 0.5d0 * rho
        rhob = 0.5d0 * rho
        grhoa = 0.25d0 * grho
        grhob = 0.25d0 * grho
        CALL LSD_B88(B1,RHOA,RHOB,GRHOA,GRHOB,sx,V1XA,V2XA,V1XB,V2XB)
        v1x = V1XA
        v2x = V2XA
      END SUBROUTINE wrap_b88

      SUBROUTINE wrap_glyp( rho, grho, sc, v1c, v2c )
        USE kinds, ONLY: DP
        IMPLICIT NONE
        REAL(DP) :: rho, grho, sc, v1c, v2c
        REAL(DP) :: RA,RB,GRHOAA,GRHOAB,GRHOBB
        REAL(DP) :: V1CA,V2CA,V1CB,V2CB,V2CAB
        ra = rho * 0.5d0
        rb = rho * 0.5d0
        grhoaa = 0.25d0 * grho
        grhobb = 0.25d0 * grho
        grhoab = 0.25d0 * grho
        CALL LSD_GLYP(RA,RB,GRHOAA,GRHOAB,GRHOBB,SC,  &
                      V1CA,V2CA,V1CB,V2CB,V2CAB)
        v1c = V1CA
        v2c = 2.0d0*(v2ca+v2cb+v2cab*2.d0)*0.25d0
      END SUBROUTINE wrap_glyp

!     ==================================================================
      SUBROUTINE LSD_B88(B1,RHOA,RHOB,GRHOA,GRHOB,sx,V1XA,V2XA,V1XB,V2XB)
!     ==--------------------------------------------------------------==
! BECKE EXCHANGE: PRA 38, 3098 (1988)
      USE kinds, ONLY: DP
      IMPLICIT NONE
      REAL(DP),PARAMETER :: OB3=1.D0/3.D0, SMALL=1.D-20
      REAL(DP) :: xs, xs2, sa2b8, br1, br2, br4, ddd, gf, dgf, shm1, dd
      REAL(DP) :: dd2, grhoa, grhob, sx, b1, rhoa, rhob, v2xb, aa, a 
      REAL(DP) :: v1xa, v2xa, v1xb

!     ==--------------------------------------------------------------==
      sx=0.0D0
      V1XA=0.0D0
      V2XA=0.0D0
      V1XB=0.0D0
      V2XB=0.0D0
      IF(ABS(RHOA).GT.SMALL) THEN
        AA    = GRHOA
        A     = SQRT(AA)
        BR1   = RHOA**OB3
        BR2   = BR1*BR1
        BR4   = BR2*BR2
        XS    = A/BR4
        XS2   = XS*XS
        SA2B8 = SQRT(1.0D0+XS2)
        SHM1  = LOG(XS+SA2B8)
        DD    = 1.0D0 + 6.0D0*B1*XS*SHM1
        DD2   = DD*DD
        DDD   = 6.0D0*B1*(SHM1+XS/SA2B8)
        GF    = -B1*XS2/DD
        DGF   = (-2.0D0*B1*XS*DD + B1*XS2*DDD)/DD2
        sx    = GF*BR4
        V1XA  = 4.d0/3.d0*BR1*(GF-XS*DGF)
        V2XA  = DGF/A
      ENDIF
      IF(ABS(RHOB).GT.SMALL) THEN
        AA    = GRHOB
        A     = SQRT(AA)
        BR1   = RHOB**OB3
        BR2   = BR1*BR1
        BR4   = BR2*BR2
        XS    = A/BR4
        XS2   = XS*XS
        SA2B8 = SQRT(1.0D0+XS2)
        SHM1  = LOG(XS+SA2B8)
        DD    = 1.0D0 + 6.0D0*B1*XS*SHM1
        DD2   = DD*DD
        DDD   = 6.0D0*B1*(SHM1+XS/SA2B8)
        GF    = -B1*XS2/DD
        DGF   = (-2.0D0*B1*XS*DD + B1*XS2*DDD)/DD2
        sx    = sx+GF*BR4
        V1XB  = 4.d0/3.d0*BR1*(GF-XS*DGF)
        V2XB  = DGF/A
      ENDIF
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_B88
