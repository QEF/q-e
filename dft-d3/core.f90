! dftd3 program for computing the dispersion energy and forces from cart
! and atomic numbers as described in
!
! S. Grimme, J. Antony, S. Ehrlich and H. Krieg
! J. Chem. Phys, 132 (2010), 154104
!
! S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem, 32 (2011), 1456
! (for BJ-damping)
!
! Copyright (C) 2009 - 2011 Stefan Grimme, University of Muenster, Germany
!
! Repackaging of the original code without any change in the functionality:
!
! Copyright (C) 2016, Bálint Aradi
!
! MPI parallelization  added by Paolo Giannozzi, June 2021
! OpenACC acceleration added by Ivan Carnimeo,   June 2021
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 1, or (at your option)
! any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! For the GNU General Public License, see <http://www.gnu.org/licenses/>
!

module dftd3_core
  use dftd3_sizes
  use dftd3_common
  use dftd3_pars
  implicit none


  ! atomic <r^2>/<r^4> values
  real(wp) r2r4(max_elem)

  ! r2r4 =sqrt(0.5*r2r4(i)*dfloat(i)**0.5 ) with i=elementnumber
  ! the large number of digits is just to keep the results consistent
  ! with older versions. They should not imply any higher accuracy than
  ! the old values
  data r2r4 / &
    &2.00734898_wp, 1.56637132_wp, 5.01986934_wp, 3.85379032_wp, 3.64446594_wp,&
    &3.10492822_wp, 2.71175247_wp, 2.59361680_wp, 2.38825250_wp, 2.21522516_wp,&
    &6.58585536_wp, 5.46295967_wp, 5.65216669_wp, 4.88284902_wp, 4.29727576_wp,&
    &4.04108902_wp, 3.72932356_wp, 3.44677275_wp, 7.97762753_wp, 7.07623947_wp,&
    &6.60844053_wp, 6.28791364_wp, 6.07728703_wp, 5.54643096_wp, 5.80491167_wp,&
    &5.58415602_wp, 5.41374528_wp, 5.28497229_wp, 5.22592821_wp, 5.09817141_wp,&
    &6.12149689_wp, 5.54083734_wp, 5.06696878_wp, 4.87005108_wp, 4.59089647_wp,&
    &4.31176304_wp, 9.55461698_wp, 8.67396077_wp, 7.97210197_wp, 7.43439917_wp,&
    &6.58711862_wp, 6.19536215_wp, 6.01517290_wp, 5.81623410_wp, 5.65710424_wp,&
    &5.52640661_wp, 5.44263305_wp, 5.58285373_wp, 7.02081898_wp, 6.46815523_wp,&
    &5.98089120_wp, 5.81686657_wp, 5.53321815_wp, 5.25477007_wp,11.02204549_wp,&
    &0.15679528_wp, 9.35167836_wp, 9.06926079_wp, 8.97241155_wp, 8.90092807_wp,&
    &8.85984840_wp, 8.81736827_wp, 8.79317710_wp, 7.89969626_wp, 8.80588454_wp,&
    &8.42439218_wp, 8.54289262_wp, 8.47583370_wp, 8.45090888_wp, 8.47339339_wp,&
    &7.83525634_wp, 8.20702843_wp, 7.70559063_wp, 7.32755997_wp, 7.03887381_wp,&
    &6.68978720_wp, 6.05450052_wp, 5.88752022_wp, 5.70661499_wp, 5.78450695_wp,&
    &7.79780729_wp, 7.26443867_wp, 6.78151984_wp, 6.67883169_wp, 6.39024318_wp,&
    &6.09527958_wp,11.79156076_wp,11.10997644_wp, 9.51377795_wp, 8.67197068_wp,&
    &8.77140725_wp, 8.65402716_wp, 8.53923501_wp, 8.85024712_wp /

  ! PBE0/def2-QZVP atomic values
  ! data r2r4 /
  ! . 8.0589, 3.4698, 29.0974, 14.8517, 11.8799, 7.8715, 5.5588,
  ! . 4.7566, 3.8025, 3.1036, 26.1552, 17.2304, 17.7210, 12.7442,
  ! . 9.5361, 8.1652, 6.7463, 5.6004, 29.2012, 22.3934, 19.0598,
  ! . 16.8590, 15.4023, 12.5589, 13.4788, 12.2309, 11.2809, 10.5569,
  ! . 10.1428, 9.4907, 13.4606, 10.8544, 8.9386, 8.1350, 7.1251,
  ! . 6.1971, 30.0162, 24.4103, 20.3537, 17.4780, 13.5528, 11.8451,
  ! . 11.0355, 10.1997, 9.5414, 9.0061, 8.6417, 8.9975, 14.0834,
  ! . 11.8333, 10.0179, 9.3844, 8.4110, 7.5152, 32.7622, 27.5708,
  ! . 23.1671, 21.6003, 20.9615, 20.4562, 20.1010, 19.7475, 19.4828,
  ! . 15.6013, 19.2362, 17.4717, 17.8321, 17.4237, 17.1954, 17.1631,
  ! . 14.5716, 15.8758, 13.8989, 12.4834, 11.4421, 10.2671, 8.3549,
  ! . 7.8496, 7.3278, 7.4820, 13.5124, 11.6554, 10.0959, 9.7340,
  ! . 8.8584, 8.0125, 29.8135, 26.3157, 19.1885, 15.8542, 16.1305,
  ! . 15.6161, 15.1226, 16.1576 /


  ! scale r4/r2 values of the atoms by sqrt(Z)
  ! sqrt is also globally close to optimum
  ! together with the factor 1/2 this yield reasonable
  ! c8 for he, ne and ar. for larger Z, C8 becomes too large
  ! which effectively mimics higher R^n terms neglected due
  ! to stability reasons

  ! r2r4 =sqrt(0.5*r2r4(i)*dfloat(i)**0.5 ) with i=elementnumber
  ! the large number of digits is just to keep the results consistent
  ! with older versions. They should not imply any higher accuracy than
  ! the old values
  !!data r2r4 / &
  !!    & 2.00734898, 1.56637132, 5.01986934, 3.85379032, 3.64446594, &
  !!    & 3.10492822, 2.71175247, 2.59361680, 2.38825250, 2.21522516, &
  !!    & 6.58585536, 5.46295967, 5.65216669, 4.88284902, 4.29727576, &
  !!    & 4.04108902, 3.72932356, 3.44677275, 7.97762753, 7.07623947, &
  !!    & 6.60844053, 6.28791364, 6.07728703, 5.54643096, 5.80491167, &
  !!    & 5.58415602, 5.41374528, 5.28497229, 5.22592821, 5.09817141, &
  !!    & 6.12149689, 5.54083734, 5.06696878, 4.87005108, 4.59089647, &
  !!    & 4.31176304, 9.55461698, 8.67396077, 7.97210197, 7.43439917, &
  !!    & 6.58711862, 6.19536215, 6.01517290, 5.81623410, 5.65710424, &
  !!    & 5.52640661, 5.44263305, 5.58285373, 7.02081898, 6.46815523, &
  !!    & 5.98089120, 5.81686657, 5.53321815, 5.25477007, 11.02204549, &
  !!    &10.15679528, 9.35167836, 9.06926079, 8.97241155, 8.90092807, &
  !!    & 8.85984840, 8.81736827, 8.79317710, 7.89969626, 8.80588454, &
  !!    & 8.42439218, 8.54289262, 8.47583370, 8.45090888, 8.47339339, &
  !!    & 7.83525634, 8.20702843, 7.70559063, 7.32755997, 7.03887381, &
  !!    & 6.68978720, 6.05450052, 5.88752022, 5.70661499, 5.78450695, &
  !!    & 7.79780729, 7.26443867, 6.78151984, 6.67883169, 6.39024318, &
  !!    & 6.09527958, 11.79156076, 11.10997644, 9.51377795, 8.67197068, &
  !!    & 8.77140725, 8.65402716, 8.53923501, 8.85024712 /


  real(wp) rcov(max_elem)

  ! covalent radii
  ! covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
  ! values for metals decreased by 10 %
  ! data rcov/
  ! . 0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67
  ! ., 1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54
  ! ., 1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09
  ! ., 1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39
  ! ., 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26
  ! ., 1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57
  ! ., 1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53
  ! ., 1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32
  ! ., 1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58
  ! ., 1.52, 1.53, 1.54, 1.55 /

  ! these new data are scaled with k2=4./3. and converted a_0 via
  ! autoang=0.52917726d0
  !!data rcov/ &
  !!    & 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865, &
  !!    & 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527, &
  !!    & 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820, &
  !!    & 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730, &
  !!    & 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923, &
  !!    & 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188, &
  !!    & 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349, &
  !!    & 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216, &
  !!    & 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717, &
  !!    & 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967, &
  !!    & 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625, &
  !!    & 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657, &
  !!    & 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833, &
  !!    & 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098, &
  !!    & 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878, &
  !!    & 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790, &
  !!    & 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584, &
  !!    & 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289, &
  !!    & 3.82984466, 3.85504098, 3.88023730, 3.90543362 /

  ! these new data are scaled with k2=4./3. and converted a_0 via
  ! autoang=0.52917726d0
  data rcov/ &
  & 0.80628308_wp, 1.15903197_wp, 3.02356173_wp, 2.36845659_wp, 1.94011865_wp, &
  & 1.88972601_wp, 1.78894056_wp, 1.58736983_wp, 1.61256616_wp, 1.68815527_wp, &
  & 3.52748848_wp, 3.14954334_wp, 2.84718717_wp, 2.62041997_wp, 2.77159820_wp, &
  & 2.57002732_wp, 2.49443835_wp, 2.41884923_wp, 4.43455700_wp, 3.88023730_wp, &
  & 3.35111422_wp, 3.07395437_wp, 3.04875805_wp, 2.77159820_wp, 2.69600923_wp, &
  & 2.62041997_wp, 2.51963467_wp, 2.49443835_wp, 2.54483100_wp, 2.74640188_wp, &
  & 2.82199085_wp, 2.74640188_wp, 2.89757982_wp, 2.77159820_wp, 2.87238349_wp, &
  & 2.94797246_wp, 4.76210950_wp, 4.20778980_wp, 3.70386304_wp, 3.50229216_wp, &
  & 3.32591790_wp, 3.12434702_wp, 2.89757982_wp, 2.84718717_wp, 2.84718717_wp, &
  & 2.72120556_wp, 2.89757982_wp, 3.09915070_wp, 3.22513231_wp, 3.17473967_wp, &
  & 3.17473967_wp, 3.09915070_wp, 3.32591790_wp, 3.30072128_wp, 5.26603625_wp, &
  & 4.43455700_wp, 4.08180818_wp, 3.70386304_wp, 3.98102289_wp, 3.95582657_wp, &
  & 3.93062995_wp, 3.90543362_wp, 3.80464833_wp, 3.82984466_wp, 3.80464833_wp, &
  & 3.77945201_wp, 3.75425569_wp, 3.75425569_wp, 3.72905937_wp, 3.85504098_wp, &
  & 3.67866672_wp, 3.45189952_wp, 3.30072128_wp, 3.09915070_wp, 2.97316878_wp, &
  & 2.92277614_wp, 2.79679452_wp, 2.82199085_wp, 2.84718717_wp, 3.32591790_wp, &
  & 3.27552496_wp, 3.27552496_wp, 3.42670319_wp, 3.30072128_wp, 3.47709584_wp, &
  & 3.57788113_wp, 5.06446567_wp, 4.56053862_wp, 4.20778980_wp, 3.98102289_wp, &
  & 3.82984466_wp, 3.85504098_wp, 3.88023730_wp, 3.90543362_wp /


contains


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! set parameters
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine setfuncpar(func,version,TZ,s6,rs6,s18,rs18,alp)
    integer version
    real(wp) s6,rs6,s18,alp,rs18
    character*(*) func
    logical TZ
    ! double hybrid values revised according to procedure in the GMTKN30 pap


    if(version.eq.6)then
    s6  =1.0d0
    alp =14.0d0
! BJ damping with parameters from ...
    select case (func)
      case ("b2-plyp")
           rs6 =0.486434
           s18 =0.672820
           rs18=3.656466 
     case ("b3-lyp") 
           rs6 =0.278672 
           s18 =1.466677
           rs18=4.606311
      case ("b97-d")
           rs6 =0.240184
           s18 =1.206988
           rs18=3.864426
      case ("b-lyp")
           rs6 =0.448486
           s18 =1.875007
           rs18=3.610679
      case ("b-p")
           rs6 =0.821850
           s18 =3.140281
           rs18=2.728151
      case ("pbe")
           rs6 =0.012092
           s18 =0.358940
           rs18=5.938951
      case ("pbe0")
           rs6 =0.007912
           s18 =0.528823
           rs18=6.162326
      case ("lc-wpbe")
           rs6 =0.563761
           s18 =0.906564
           rs18=3.593680
      case DEFAULT
            call stoprun( 'functional name unknown' )
    end select
    endif

    if(version.eq.5)then
    s6  =1.0d0
    alp =14.0d0
! zero damping with parameters from ...
    select case (func)
      case ("b2-plyp")
           rs6 =1.313134
           s18 =0.717543
           rs18=0.016035
           s6  =0.640000
      case ("b3-lyp")
           rs6 =1.338153
           s18 =1.532981
           rs18=0.013988
      case ("b97-d")
           rs6 =1.151808
           s18 =1.020078
           rs18=0.035964
      case ("b-lyp")
           rs6 =1.279637
           s18 =1.841686
           rs18=0.014370
      case ("b-p")
           rs6 =1.233460
           s18 =1.945174
           rs18=0.000000
      case ("pbe")
           rs6 =2.340218
           s18 =0.000000
           rs18=0.129434
      case ("pbe0")
           rs6 =2.077949
           s18 =0.000081
           rs18=0.116755
      case ("lc-wpbe")
           rs6 =1.366361
           s18 =1.280619
           rs18=0.003160
      case DEFAULT
            call stoprun( 'functional name unknown' )
    end select
    endif

    ! DFT-D3 with Becke-Johnson finite-damping, variant 2 with their radii
    ! SE: Alp is only used in 3-body calculations
    if (version.eq.4)then
      s6=1.0d0
      alp =14.0d0

      select case (func)
      case ("b-p")
        rs6 =0.3946
        s18 =3.2822
        rs18=4.8516
      case ("b-lyp")
        rs6 =0.4298
        s18 =2.6996
        rs18=4.2359
      case ("revpbe")
        rs6 =0.5238
        s18 =2.3550
        rs18=3.5016
      case ("rpbe")
        rs6 =0.1820
        s18 =0.8318
        rs18=4.0094
      case ("b97-d")
        rs6 =0.5545
        s18 =2.2609
        rs18=3.2297
      case ("pbe")
        rs6 =0.4289
        s18 =0.7875
        rs18=4.4407
      case ("rpw86-pbe")
        rs6 =0.4613
        s18 =1.3845
        rs18=4.5062
      case ("b3-lyp")
        rs6 =0.3981
        s18 =1.9889
        rs18=4.4211
      case ("tpss")
        rs6 =0.4535
        s18 =1.9435
        rs18=4.4752
      case ("hf")
        rs6 =0.3385
        s18 =0.9171
        rs18=2.8830
      case ("tpss0")
        rs6 =0.3768
        s18 =1.2576
        rs18=4.5865
      case ("pbe0")
        rs6 =0.4145
        s18 =1.2177
        rs18=4.8593
      case ("hse06")
        rs6 =0.383
        s18 =2.310
        rs18=5.685
      case ("revpbe38")
        rs6 =0.4309
        s18 =1.4760
        rs18=3.9446
      case ("pw6b95")
        rs6 =0.2076
        s18 =0.7257
        rs18=6.3750
      case ("b2-plyp")
        rs6 =0.3065
        s18 =0.9147
        rs18=5.0570
        s6=0.64d0
      case ("dsd-blyp")
        rs6 =0.0000
        s18 =0.2130
        rs18=6.0519
        s6=0.50d0
      case ("dsd-blyp-fc")
        rs6 =0.0009
        s18 =0.2112
        rs18=5.9807
        s6=0.50d0
      case ("bop")
        rs6 =0.4870
        s18 =3.2950
        rs18=3.5043
      case ("mpwlyp")
        rs6 =0.4831
        s18 =2.0077
        rs18=4.5323
      case ("o-lyp")
        rs6 =0.5299
        s18 =2.6205
        rs18=2.8065
      case ("pbesol")
        rs6 =0.4466
        s18 =2.9491
        rs18=6.1742
      case ("bpbe")
        rs6 =0.4567
        s18 =4.0728
        rs18=4.3908
      case ("opbe")
        rs6 =0.5512
        s18 =3.3816
        rs18=2.9444
      case ("ssb")
        rs6 =-0.0952
        s18 =-0.1744
        rs18=5.2170
      case ("revssb")
        rs6 =0.4720
        s18 =0.4389
        rs18=4.0986
      case ("otpss")
        rs6 =0.4634
        s18 =2.7495
        rs18=4.3153
      case ("b3pw91")
        rs6 =0.4312
        s18 =2.8524
        rs18=4.4693
      case ("bh-lyp")
        rs6 =0.2793
        s18 =1.0354
        rs18=4.9615
      case ("revpbe0")
        rs6 =0.4679
        s18 =1.7588
        rs18=3.7619
      case ("tpssh")
        rs6 =0.4529
        s18 =2.2382
        rs18=4.6550
      case ("mpw1b95")
        rs6 =0.1955
        s18 =1.0508
        rs18=6.4177
      case ("pwb6k")
        rs6 =0.1805
        s18 =0.9383
        rs18=7.7627
      case ("b1b95")
        rs6 =0.2092
        s18 =1.4507
        rs18=5.5545
      case ("bmk")
        rs6 =0.1940
        s18 =2.0860
        rs18=5.9197
      case ("cam-b3lyp")
        rs6 =0.3708
        s18 =2.0674
        rs18=5.4743
      case ("lc-wpbe")
        rs6 =0.3919
        s18 =1.8541
        rs18=5.0897
      case ("b2gp-plyp")
        rs6 =0.0000
        s18 =0.2597
        rs18=6.3332
        s6=0.560
      case ("ptpss")
        rs6 =0.0000
        s18 =0.2804
        rs18=6.5745
        s6=0.750
      case ("pwpb95")
        rs6 =0.0000
        s18 =0.2904
        rs18=7.3141
        s6=0.820
        ! special HF/DFT with eBSSE correction
      case ("hf/mixed")
        rs6 =0.5607
        s18 =3.9027
        rs18=4.5622
      case ("hf/sv")
        rs6 =0.4249
        s18 =2.1849
        rs18=4.2783
      case ("hf/minis")
        rs6 =0.1702
        s18 =0.9841
        rs18=3.8506
      case ("b3-lyp/6-31gd")
        rs6 =0.5014
        s18 =4.0672
        rs18=4.8409
      case ("hcth120")
        rs6=0.3563
        s18=1.0821
        rs18=4.3359
        ! DFTB3 old, deprecated parameters:
        ! case ("dftb3")
        ! rs6=0.7461
        ! s18=3.209
        ! rs18=4.1906
        ! special SCC-DFTB parametrization
        ! full third order DFTB, self consistent charges, hydrogen pair damping
        ! exponent 4.2
      case("dftb3")
        rs6=0.5719d0
        s18=0.5883d0
        rs18=3.6017d0
      case ("pw1pw")
        rs6 =0.3807d0
        s18 =2.3363d0
        rs18=5.8844d0
      case ("pwgga")
        rs6 =0.2211d0
        s18 =2.6910d0
        rs18=6.7278d0
      case ("hsesol")
        rs6 =0.4650d0
        s18 =2.9215d0
        rs18=6.2003d0
        ! special HF-D3-gCP-SRB/MINIX parametrization
      case ("hf3c")
        rs6=0.4171d0
        s18=0.8777d0
        rs18=2.9149d0
        ! special HF-D3-gCP-SRB2/ECP-2G parametrization
      case ("hf3cv")
        rs6=0.3063d0
        s18=0.5022d0
        rs18=3.9856d0
        ! special PBEh-D3-gCP/def2-mSVP parametrization
      case ("pbeh3c", "pbeh-3c")
        rs6=0.4860d0
        s18=0.0000d0
        rs18=4.5000d0
      case ("scan")
        ! Parameters from PRB 94, 115144 (2016); doi: 10.1103/PhysRevB.94.115144
        ! Table 1
        rs6 =0.538
        s18 =0.0
        rs18=5.4200
      case ("r2scan")
        ! Parameters from JCP, 154, 061101 (2021); doi: 10.1063/5.0041008
        ! Table 1
        rs6 =0.4948
        s18 =0.7898
        rs18=5.7308

      case DEFAULT
        call stoprun( 'functional name unknown' )
      end select
    end if

    ! DFT-D3
    if (version.eq.3)then
      s6 =1.0d0
      alp =14.0d0
      rs18=1.0d0
      ! default def2-QZVP (almost basis set limit)
      if (.not.TZ) then
        select case (func)
        case ("slater-dirac-exchange")
          rs6 =0.999
          s18 =-1.957
          rs18=0.697
        case ("b-lyp")
          rs6=1.094
          s18=1.682
        case ("b-p")
          rs6=1.139
          s18=1.683
        case ("b97-d")
          rs6=0.892
          s18=0.909
        case ("revpbe")
          rs6=0.923
          s18=1.010
        case ("pbe")
          rs6=1.217
          s18=0.722
        case ("pbesol")
          rs6=1.345
          s18=0.612
        case ("rpw86-pbe")
          rs6=1.224
          s18=0.901
        case ("rpbe")
          rs6=0.872
          s18=0.514
        case ("tpss")
          rs6=1.166
          s18=1.105
        case ("b3-lyp")
          rs6=1.261
          s18=1.703
        case ("pbe0")
          rs6=1.287
          s18=0.928

        case ("hse06")
          rs6=1.129
          s18=0.109
        case ("revpbe38")
          rs6=1.021
          s18=0.862
        case ("pw6b95")
          rs6=1.532
          s18=0.862
        case ("tpss0")
          rs6=1.252
          s18=1.242
        case ("b2-plyp")
          rs6=1.427
          s18=1.022
          s6=0.64
        case ("pwpb95")
          rs6=1.557
          s18=0.705
          s6=0.82
        case ("b2gp-plyp")
          rs6=1.586
          s18=0.760
          s6=0.56
        case ("ptpss")
          rs6=1.541
          s18=0.879
          s6=0.75
        case ("hf")
          rs6=1.158
          s18=1.746
        case ("mpwlyp")
          rs6=1.239
          s18=1.098
        case ("bpbe")
          rs6=1.087
          s18=2.033
        case ("bh-lyp")
          rs6=1.370
          s18=1.442
        case ("tpssh")
          rs6=1.223
          s18=1.219
        case ("pwb6k")
          rs6=1.660
          s18=0.550
        case ("b1b95")
          rs6=1.613
          s18=1.868
        case ("bop")
          rs6=0.929
          s18=1.975
        case ("o-lyp")
          rs6=0.806
          s18=1.764
        case ("o-pbe")
          rs6=0.837
          s18=2.055
        case ("ssb")
          rs6=1.215
          s18=0.663
        case ("revssb")
          rs6=1.221
          s18=0.560
        case ("otpss")
          rs6=1.128
          s18=1.494
        case ("b3pw91")
          rs6=1.176
          s18=1.775
        case ("revpbe0")
          rs6=0.949
          s18=0.792
        case ("pbe38")
          rs6=1.333
          s18=0.998
        case ("mpw1b95")
          rs6=1.605
          s18=1.118
        case ("mpwb1k")
          rs6=1.671
          s18=1.061
        case ("bmk")
          rs6=1.931
          s18=2.168
        case ("cam-b3lyp")
          rs6=1.378
          s18=1.217
        case ("lc-wpbe")
          rs6=1.355
          s18=1.279
        case ("m05")
          rs6=1.373
          s18=0.595
        case ("m052x")
          rs6=1.417
          s18=0.000
        case ("m06l")
          rs6=1.581
          s18=0.000
        case ("m06")
          rs6=1.325
          s18=0.000
        case ("m062x")
          rs6=1.619
          s18=0.000
        case ("m06hf")
          rs6=1.446
          s18=0.000
          ! DFTB3 (zeta=4.0), old deprecated parameters
          ! case ("dftb3")
          ! rs6=1.235
          ! s18=0.673
        case ("hcth120")
          rs6=1.221
          s18=1.206
        case DEFAULT
          call stoprun( 'functional name unknown' )
        end select
      else
        ! special TZVPP parameter
        select case (func)
        case ("b-lyp")
          rs6=1.243
          s18=2.022
        case ("b-p")
          rs6=1.221
          s18=1.838
        case ("b97-d")
          rs6=0.921
          s18=0.894
        case ("revpbe")
          rs6=0.953
          s18=0.989
        case ("pbe")
          rs6=1.277
          s18=0.777
        case ("tpss")
          rs6=1.213
          s18=1.176
        case ("b3-lyp")
          rs6=1.314
          s18=1.706
        case ("pbe0")
          rs6=1.328
          s18=0.926
        case ("pw6b95")
          rs6=1.562
          s18=0.821
        case ("tpss0")
          rs6=1.282
          s18=1.250
        case ("b2-plyp")
          rs6=1.551
          s18=1.109
          s6=0.5
        case DEFAULT
          call stoprun( 'functional name unknown (TZ case)' )
        end select
      end if
    end if
    ! DFT-D2
    if (version.eq.2)then
      rs6=1.1d0
      s18=0.0d0
      alp=20.0d0
      select case (func)
      case ("b-lyp")
        s6=1.2
      case ("b-p")
        s6=1.05
      case ("b97-d")
        s6=1.25
      case ("revpbe")
        s6=1.25
      case ("pbe")
        s6=0.75
      case ("tpss")
        s6=1.0
      case ("b3-lyp")
        s6=1.05
      case ("pbe0")
        s6=0.6
      case ("pw6b95")
        s6=0.5
      case ("tpss0")
        s6=0.85
      case ("b2-plyp")
        s6=0.55
      case ("b2gp-plyp")
        s6=0.4
      case ("dsd-blyp")
        s6=0.41
        alp=60.0d0
      case DEFAULT
        call stoprun( 'functional name unknown' )
      end select

    end if

  end subroutine setfuncpar


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! compute energy
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
      & rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr,&
      & e6,e8,e10,e12,e63)
    integer n,iz(*),max_elem,maxc,version,mxc(max_elem)
    real(wp) xyz(3,*),r0ab(max_elem,max_elem),r2r4(*)
    real(wp) rs6,rs8,rs10,alp6,alp8,alp10,rcov(max_elem)
    real(wp) c6ab(max_elem,max_elem,maxc,maxc,3)
    real(wp) e6, e8, e10, e12, e63
    logical noabc

    integer iat,jat,kat
    real(wp) r,r2,r6,r8,tmp,alp,dx,dy,dz,c6,c8,c10,ang,rav
    real(wp) damp6,damp8,damp10,rr,thr,c9,r42,c12,r10,c14,rthr,cn_thr
    real(wp) cn(n)
    real(wp) r2ab(n*n),cc6ab(n*n),dmp(n*n),d2(3),t1,t2,t3,a1,a2,tmp2
    real(wp) abcthr
    integer(int64) icomp(n*n)
    integer ij,ik,jk,k

    e6 =0
    e8 =0
    e10=0
    e12=0
    e63=0
    ! the threebody term uses the same threshold as the
    abcthr=cn_thr

    ! Becke-Johnson parameters
    a1=rs6
    a2=rs8


    ! DFT-D2
    if (version.eq.2)then

      do iat=1,n-1
        do jat=iat+1,n
          dx=xyz(1,iat)-xyz(1,jat)
          dy=xyz(2,iat)-xyz(2,jat)
          dz=xyz(3,iat)-xyz(3,jat)
          r2=dx*dx+dy*dy+dz*dz
          ! if (r2.gt.rthr) cycle
          r=sqrt(r2)
          c6=c6ab(iz(jat),iz(iat),1,1,1)
          damp6=1./(1.+exp(-alp6*(r/(rs6*r0ab(iz(jat),iz(iat)))-1.)))
          r6=r2**3
          e6 =e6+c6*damp6/r6
        end do
      end do

    else
      ! DFT-D3
      call ncoord(n,rcov,iz,xyz,cn,cn_thr)

      icomp=0
      do iat=1,n-1
        do jat=iat+1,n
          dx=xyz(1,iat)-xyz(1,jat)
          dy=xyz(2,iat)-xyz(2,jat)
          dz=xyz(3,iat)-xyz(3,jat)
          r2=dx*dx+dy*dy+dz*dz
          !THR
          if (r2.gt.rthr) cycle
          r =sqrt(r2)
          tmp2=r0ab(iz(jat),iz(iat))
          rr=tmp2/r
          ! damping
          if(version.eq.3)then
            ! DFT-D3 zero-damp
            tmp=rs6*rr
            damp6 =1.d0/( 1.d0+6.d0*tmp**alp6 )
            tmp=rs8*rr
            damp8 =1.d0/( 1.d0+6.d0*tmp**alp8 )
          else
            ! DFT-D3M zero-damp
            tmp=(r/(rs6*tmp2))+rs8*tmp2
            damp6 =1.d0/( 1.d0+6.d0*tmp**(-alp6) )
            tmp=(r/tmp2)+rs8*tmp2
            damp8 =1.d0/( 1.d0+6.d0*tmp**(-alp8) )
          endif
          ! get C6
          call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat), &
              & cn(iat),cn(jat),c6)

          r6=r2**3
          r8=r6*r2
          ! r2r4 stored in main as sqrt
          c8 =3.0d0*c6*r2r4(iz(iat))*r2r4(iz(jat))

          ! DFT-D3 zero-damp or DFT-D3 M(zero)
          if((version.eq.3).or.(version.eq.5))then
            e6=e6+c6*damp6/r6
            e8=e8+c8*damp8/r8
          end if
          ! DFT-D3(BJ) or DFT-D3M(BJ)
          if((version.eq.4).or.(version.eq.6))then
            ! use BJ radius
            tmp=sqrt(c8/c6)
            e6=e6+ c6/(r6+(a1*tmp+a2)**6)
            e8=e8+ c8/(r8+(a1*tmp+a2)**8)
          end if

          ! if (.not.noabc) then
          if ((.not.noabc).and.(r2.lt.abcthr)) then
            ij=lin(jat,iat)
            icomp(ij)=1
            ! store C6 for C9, calc as sqrt
            cc6ab(ij)=sqrt(c6)
            ! store R^2 for abc
            r2ab(ij)=r2
            ! store for abc damping
            dmp(ij)=(1./rr)**(1./3.)
            !noabc
          end if
        end do
      end do

      if (.not.noabc) then 

      ! compute non-additive third-order energy using averaged C6
      do iat=1,n
        do jat=1,iat-1
          ij=lin(jat,iat)
          if (icomp(ij).eq.1)then
            do kat=1,jat-1

              ik=lin(kat,iat)
              jk=lin(kat,jat)
              if ((icomp(ik).eq.0).or.(icomp(jk).eq.0)) cycle
              ! damping func product
              ! tmp=dmp(ik)*dmp(jk)*dmp(ij)
              rav=(4./3.)/(dmp(ik)*dmp(jk)*dmp(ij))
              tmp=1.d0/( 1.d0+6.d0*rav**alp8 )
              ! triple C6 coefficient (stored as sqrt)
              c9=cc6ab(ij)*cc6ab(ik)*cc6ab(jk)
              ! write(*,*) 'C9', c9
              ! angular terms with law of cosines

              t1 = (r2ab(ij)+r2ab(jk)-r2ab(ik))
              t2 = (r2ab(ij)+r2ab(ik)-r2ab(jk))
              t3 = (r2ab(ik)+r2ab(jk)-r2ab(ij))
              tmp2=r2ab(ij)*r2ab(jk)*r2ab(ik)
              ang=(0.375d0*t1*t2*t3/tmp2+1.0d0)/tmp2**1.5d0

              e63=e63-tmp*c9*ang

            end do
          end if
        end do
      end do

    end if

    end if
 
  end subroutine edisp


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! compute gradient
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine gdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
      & s6,s18,rs6,rs8,rs10,alp6,alp8,alp10,noabc,rthr, &
      & num,version,echo,g,disp,gnorm,cn_thr,fix)

    integer n,iz(*),max_elem,maxc,version,mxc(max_elem)
    real(wp) xyz(3,*),r0ab(max_elem,max_elem),r2r4(*)
    real(wp) c6ab(max_elem,max_elem,maxc,maxc,3)
    real(wp) g(3,*),s6,s18,rcov(max_elem)
    real(wp) rs6,rs8,rs10,alp10,alp8,alp6,a1,a2
    logical noabc,num,echo,fix(n)

    integer iat,jat,i,j,kat,k
    real(wp) R0,C6,alp,R42,disp,x1,y1,z1,x2,y2,z2,rr,e6abc
    real(wp) dx,dy,dz,r2,r,r4,r6,r8,r10,r12,t6,t8,t10,damp1
    real(wp) damp6,damp8,damp10,e6,e8,e10,e12,gnorm,tmp1
    real(wp) s10,s8,gC6(3),term,step,dispr,displ,r235,tmp2
    real(wp) cn(n),gx1,gy1,gz1,gx2,gy2,gz2,rthr,c8,cn_thr
    real(wp) rthr3

    real(wp) rij(3),rik(3),rjk(3),r7,r9
    real(wp) rik_dist,rjk_dist
    !d(E)/d(r_ij) derivative wrt. dist. iat-ja
    real(wp) drij(n*(n+1)/2)
    real(wp) drik,drjk
    real(wp) rcovij
    !d(C6ij)/d(r_ij)
    real(wp) dc6,c6chk
    real(wp) expterm,dcni
    !dCN(iat)/d(r_ij) is equal to
    real(wp) dcn
    !dCN(jat)/d(r_ij)
    ! saves (1/r^6*f_dmp + 3*r4r2/r^8*f_dmp) for kat l
    real(wp) dc6_rest
    integer linij,linik,linjk
    real(wp) vec(3),vec2(3)
    ! dE_disp/dCN(iat) in dc6i(iat)
    real(wp) dc6i(n)
    ! saves dC6(ij)/dCN(iat)
    real(wp) dc6ij(n,n)
    real(wp) dc6iji,dc6ijj
    logical abccalc(n*(n+1)/2)
    real(wp) abcthr
    ! temporary container to store iat gradient
    real(wp) gtmp(3)
    real(wp) labc,rabc
    real(wp) c6abc(n*(n+1)/2)
    real(wp) r2abc(n*(n+1)/2)
    real(wp) r3abc(n*(n+1)/2)
    real(wp) c9,rav,rav3,fdmp,ang,angr9,eabc,dc9,dfdmp,dang
    real(wp) r2ij,r2jk,r2ik,mijk,imjk,ijmk,rijk3
    integer mat,linim,linjm,linkm,kk


    dc6i=0.0d0
    abccalc=.FALSE.
    abcthr=cn_thr




    !NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    if (num) then
      if (echo)write(*,*) 'doing numerical gradient O(N^3) ...'

      call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
          & rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr, &
          & e6,e8,e10,e12,e6abc)
      disp=-s6*e6-s18*e8-e6abc

      step=2.d-5

      do i=1,n
        do j=1,3
          xyz(j,i)=xyz(j,i)+step
          call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
              & rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr, &
              & e6,e8,e10,e12,e6abc)
          dispr=-s6*e6-s18*e8-e6abc
          rabc=e6abc
          xyz(j,i)=xyz(j,i)-2*step
          call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
              & rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr, &
              & e6,e8,e10,e12,e6abc)
          displ=-s6*e6-s18*e8-e6abc
          labc=e6abc
          g(j,i)=0.5*(dispr-displ)/step
          xyz(j,i)=xyz(j,i)+step
        end do
      end do
      goto 999
    end if

    ! this is the crucial threshold to reduce the N^3 to an
    ! effective N^2.

    ! write(*,*)'rthr=',rthr,'rthr2=',rthr2,'rthr3=',rthr3

    if (echo)write(*,*)
    ! 2222222222222222222222222222222222222222222222222222222222222222222222
    if (version.eq.2)then
      if (echo)write(*,*) 'doing analytical gradient O(N^2) ...'
      disp=0
      do iat=1,n-1
        do jat=iat+1,n
          R0=r0ab(iz(jat),iz(iat))*rs6
          dx=(xyz(1,iat)-xyz(1,jat))
          dy=(xyz(2,iat)-xyz(2,jat))
          dz=(xyz(3,iat)-xyz(3,jat))
          r2 =dx*dx+dy*dy+dz*dz
          ! if (r2.gt.rthr) cycle
          r235=r2**3.5
          r =dsqrt(r2)
          damp6=exp(-alp6*(r/R0-1.0d0))
          damp1=1.+damp6
          c6=c6ab(iz(jat),iz(iat),1,1,1)*s6
          tmp1=damp6/(damp1*damp1*r235*R0)
          tmp2=6./(damp1*r*r235)
          gx1=alp6* dx*tmp1-tmp2*dx
          gx2=alp6*(-dx)*tmp1+tmp2*dx
          gy1=alp6* dy*tmp1-tmp2*dy
          gy2=alp6*(-dy)*tmp1+tmp2*dy
          gz1=alp6* dz*tmp1-tmp2*dz
          gz2=alp6*(-dz)*tmp1+tmp2*dz
          g(1,iat)=g(1,iat)-gx1*c6
          g(2,iat)=g(2,iat)-gy1*c6
          g(3,iat)=g(3,iat)-gz1*c6
          g(1,jat)=g(1,jat)-gx2*c6
          g(2,jat)=g(2,jat)-gy2*c6
          g(3,jat)=g(3,jat)-gz2*c6
          disp=disp+c6*(1./damp1)/r2**3
        end do
      end do
      disp=-disp


      goto 999

      ! 3333333333333333333333333333333333333333333333333333333333333333333333
      ! zero damping
    elseif ((version.eq.3).or.(version.eq.5)) then

      if (echo)write(*,*) 'doing analytical gradient O(N^2) ...'
      call ncoord(n,rcov,iz,xyz,cn,cn_thr)
      s8 =s18
      s10=s18

      disp=0

      drij=0.0d0
      dc6_rest=0.0d0
      kat=0


      kk=0
      do iat=1,n
        do jat=1,iat-1
          rij=xyz(:,jat)-xyz(:,iat)
          r2=sum(rij*rij)
          if (r2.gt.rthr) cycle
          linij=lin(iat,jat)

          r0=r0ab(iz(jat),iz(iat))
          r42=r2r4(iz(iat))*r2r4(iz(jat))
          ! rcovij=rcov(iz(iat))+rcov(iz(jat))
          !
          ! get_dC6_dCNij calculates the derivative dC6(iat,jat)/dCN(iat) and
          ! dC6(iat,jat)/dCN(jat).
          !
          call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)), &
              & mxc(iz(jat)),cn(iat),cn(jat),iz(iat),iz(jat),iat,jat, &
              & c6,dc6iji,dc6ijj)


          r=dsqrt(r2)
          r6=r2*r2*r2
          r7=r6*r
          r8=r6*r2
          r9=r8*r

          !
          ! Calculates damping functions:
          if (version.eq.3) then
            t6 = (r/(rs6*R0))**(-alp6)
            damp6 =1.d0/( 1.d0+6.d0*t6 )
            t8 = (r/(rs8*R0))**(-alp8)
            damp8 =1.d0/( 1.d0+6.d0*t8 )
  
            tmp1=s6*6.d0*damp6*C6/r7
            tmp2=s8*6.d0*C6*r42*damp8/r9
            ! d(r^(-6))/d(r_ij)
            drij(linij)=drij(linij)-tmp1 &
                & -4.d0*tmp2
  
  
            drij(linij)=drij(linij) &
                & +tmp1*alp6*t6*damp6 &
                & +3.d0*tmp2*alp8*t8*damp8
            !d(f_dmp)/d(r_ij)

          else ! version.eq.5
            t6 = (r/(rs6*R0)+R0*rs8)**(-alp6)
            damp6 =1.d0/( 1.d0+6.d0*t6 )
            t8 = (r/R0+R0*rs8)**(-alp8)
            damp8 =1.d0/( 1.d0+6.d0*t8 )
  
            tmp1=s6*6.d0*damp6*C6/r7
            tmp2=s8*6.d0*C6*r42*damp8/r9
            ! d(r^(-6))/d(r_ij)
            drij(linij)=drij(linij)-tmp1 &  
                &  -4.d0*tmp2
  
  
            drij(linij)=drij(linij) &
                & +tmp1*alp6*t6*damp6*r/(r+rs6*R0*R0*rs8) & 
                & +3.d0*tmp2*alp8*t8*damp8*r/(r+R0*R0*rs8)
            !d(f_dmp)/d(r_ij)

          end if


          if ((.not.noabc).and.(r2.lt.abcthr)) then
            ! if (.not.noabc) then
            abccalc(linij)=.TRUE.
            dc6ij(iat,jat)=dc6iji
            dc6ij(jat,iat)=dc6ijj
            c6abc(linij)=c6
            r2abc(linij)=r2
            r3abc(linij)=(r/R0)**(1.0/3.0)
            !noabc
          end if

          dc6_rest=s6/r6*damp6+3.d0*s8*r42/r8*damp8

          ! calculate E_disp for sanity check
          disp=disp-dc6_rest*c6
          !
          !
          ! saving all f_dmp/r6*dC6(ij)/dCN(i) for each atom for later
          dc6i(iat)=dc6i(iat)+dc6_rest*dc6iji
          dc6i(jat)=dc6i(jat)+dc6_rest*dc6ijj


          !jat
        end do
        !iat
      end do






      ! end if !version 3+5

      ! BJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJ
      ! Becke-Johnson finite damping
    elseif ((version.eq.4).or.(version.eq.6)) then
      a1 =rs6
      a2 =rs8
      s8 =s18

      if (echo)write(*,*) 'doing analytical gradient O(N^2) ...'
      disp=0
      call ncoord(n,rcov,iz,xyz,cn,cn_thr)

      drij=0.0d0
      dc6_rest=0.0d0
      kat=0


      do iat=1,n
        do jat=1,iat-1
          rij=xyz(:,jat)-xyz(:,iat)
          r2=sum(rij*rij)
          if (r2.gt.rthr) cycle

          linij=lin(iat,jat)
          r0=r0ab(iz(jat),iz(iat))
          r42=r2r4(iz(iat))*r2r4(iz(jat))
          !
          ! get_dC6_dCNij calculates the derivative dC6(iat,jat)/dCN(iat) and
          ! dC6(iat,jat)/dCN(jat).
          !
          call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)), &
              & mxc(iz(jat)),cn(iat),cn(jat),iz(iat),iz(jat),iat,jat, &
              & c6,dc6iji,dc6ijj)


          r=dsqrt(r2)
          r4=r2*r2
          r6=r4*r2
          r8=r6*r2

          if ((.not.noabc).and.(r2.lt.abcthr)) then
            ! if (.not.noabc) then
            abccalc(linij)=.TRUE.
            dc6ij(iat,jat)=dc6iji
            dc6ij(jat,iat)=dc6ijj
            c6abc(linij)=c6
            r2abc(linij)=r2
            r3abc(linij)=(r/r0)**(1.0/3.0)
            !noabc
          end if
          ! use BJ radius
          R0=a1*dsqrt(3.0d0*r42)+a2

          t6=(r6+R0**6)
          t8=(r8+R0**8)

          drij(linij)=drij(linij) &
              & -s6*C6*6.0d0*r4*r/(t6*t6) &
              & -s8*C6*24.0d0*r42*r6*r/(t8*t8)

          dc6_rest=s6/t6+3.d0*s8*r42/t8
          ! calculate E_disp for sanity check
          disp=disp-dc6_rest*c6

          ! saving all (1/r^6...)* dC6/dCN(i) for each atom
          dc6i(iat)=dc6i(iat)+dc6_rest*dc6iji
          dc6i(jat)=dc6i(jat)+dc6_rest*dc6ijj

          !jat
        end do
        !iat
      end do

      !version=4 (BJ)
    end if


    if (.not.noabc)then

      if (echo)write(*,*) 'doing analytical gradient O(N^3) ...'
      do iat=1,n
        do jat=1,iat-1
          linij=lin(iat,jat)
          if (.NOT.abccalc(linij))cycle
          r2ij=r2abc(linij)
          do kat=1,jat-1

            linik=lin(iat,kat)
            linjk=lin(jat,kat)
            !cuto
            if (.NOT.(abccalc(linjk).AND.abccalc(linik)))cycle
            ! calculating the 3body energy:
            r2jk=r2abc(linjk)
            r2ik=r2abc(linik)
            c9=c6abc(linij)*c6abc(linjk)*c6abc(linik)
            c9=dsqrt(c9)
            rav=r3abc(linij)*r3abc(linjk)*r3abc(linik)
            fdmp=1.0d0/(1+6.0d0*(0.75d0*rav)**(-alp8))
            mijk=-r2ij+r2jk+r2ik
            imjk= r2ij-r2jk+r2ik
            ijmk= r2ij+r2jk-r2ik
            rijk3=r2ij*r2jk*r2ik
            rav3=rijk3**1.5
            ang=0.375d0*ijmk*imjk*mijk/rijk3
            angr9=(ang +1.0d0) &
                & /rav3

            eabc=eabc+c9*angr9*fdmp
            !end of 3body energy calculation

            !start calculating the derivatives of each part w.r.t. r_ij

            r=dsqrt(r2ij)
            dfdmp=-2.d0*alp8*(0.75d0*rav)**(-alp8)*fdmp*fdmp

            dang=-0.375d0*(r2ij**3+r2ij**2*(r2jk+r2ik) &
                & +r2ij*(3.0d0*r2jk**2+2.0*r2jk*r2ik+3.0*r2ik**2) &
                & -5.0*(r2jk-r2ik)**2*(r2jk+r2ik)) &
                & /(r*rijk3*rav3)


            tmp1=dfdmp/r*c9*angr9-dang*c9*fdmp
            drij(linij)=drij(linij)+tmp1

            !start calculating the derivatives of each part w.r.t. r_jk
            r=dsqrt(r2jk)

            dang=-0.375d0*(r2jk**3+r2jk**2*(r2ik+r2ij) &
                & +r2jk*(3.0d0*r2ik**2+2.0*r2ik*r2ij+3.0*r2ij**2) &
                & -5.0*(r2ik-r2ij)**2*(r2ik+r2ij)) &
                & /(r*rijk3*rav3)

            drij(linjk)=drij(linjk) &
                & +dfdmp/r*c9*angr9-dang*c9*fdmp


            !start calculating the derivatives of each part w.r.t. r_ik
            r=dsqrt(r2abc(linik))

            dang=-0.375d0*(r2ik**3+r2ik**2*(r2jk+r2ij) &
                & +r2ik*(3.0d0*r2jk**2+2.0*r2jk*r2ij+3.0*r2ij**2) &
                & -5.0*(r2jk-r2ij)**2*(r2jk+r2ij)) &
                & /(r*rijk3*rav3)

            drij(linik)=drij(linik) &
                & +dfdmp/r*c9*angr9-dang*c9*fdmp



            ! calculate rest* dc9/dcn(iat) and sum it up for every atom ijk
            dc6_rest=angr9*fdmp

            dc9=dc6ij(iat,jat)/c6abc(linij)+dc6ij(iat,kat)/c6abc(linik)
            dc9=-0.5d0*c9*dc9
            dc6i(iat)=dc6i(iat)+dc6_rest*dc9

            dc9=dc6ij(jat,iat)/c6abc(linij)+dc6ij(jat,kat)/c6abc(linjk)
            dc9=-0.5d0*c9*dc9
            dc6i(jat)=dc6i(jat)+dc6_rest*dc9

            dc9=dc6ij(kat,iat)/c6abc(linik)+dc6ij(kat,jat)/c6abc(linjk)
            dc9=-0.5d0*c9*dc9
            dc6i(kat)=dc6i(kat)+dc6_rest*dc9

            !kat
          end do
          !jat
        end do
        !iat
      end do

      disp=disp+eabc
      !noabc
    end if



    ! After calculating all derivatives dE/dr_ij w.r.t. distances,
    ! the grad w.r.t. the coordinates is calculated dE/dr_ij * dr_ij/dxyz_i
    do iat=2,n
      gtmp=0.0
      do jat=1,iat-1
        linij=lin(iat,jat)
        rij=xyz(:,jat)-xyz(:,iat)


        r2=sum(rij*rij)
        r=dsqrt(r2)
        if (r2.lt.cn_thr) then
          rcovij=rcov(iz(iat))+rcov(iz(jat))
          expterm=exp(-k1*(rcovij/r-1.d0))
          dcn=-k1*rcovij*expterm/ &
              & (r*r*(expterm+1.d0)*(expterm+1.d0))
        else
          dcn=0.d0
        end if
        x1=drij(linij)+dcn*(dc6i(iat)+dc6i(jat))

        gtmp=gtmp+x1*rij/r
        !g(:,iat)=g(:,iat)+x1*rij/r
        g(:,jat)=g(:,jat)-x1*rij/r

        !iat
      end do
      g(:,iat)=g(:,iat)+gtmp
      !jat
    end do



999 continue
    gnorm=sum(abs(g(1:3,1:n)))
    if (echo)then
      write(*,*)
      write(*,*)'|G|=',gnorm
    end if

    !fixed atoms have no gradient
    do i=1,n
      if (fix(i))g(:,i)=0.0
    end do



  end subroutine gdisp


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! The N E W gradC6 routine C
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !
  subroutine get_dC6_dCNij(maxc,max_elem,c6ab,mxci,mxcj,cni,cnj, &
      & izi,izj,iat,jat,c6check,dc6i,dc6j)

    integer maxc,max_elem
    real(wp) c6ab(max_elem,max_elem,maxc,maxc,3)
    !mxc(iz(iat))
    integer mxci,mxcj
    real(wp) cni,cnj,term
    integer iat,jat,izi,izj
    real(wp) dc6i,dc6j,c6check


    integer i,j,a,b
    real(wp) zaehler,nenner,dzaehler_i,dnenner_i,dzaehler_j,dnenner_j
    real(wp) expterm,cn_refi,cn_refj,c6ref,r
    real(wp) c6mem,r_save



    c6mem=-1.d99
    r_save=9999.0
    zaehler=0.0d0
    nenner=0.0d0

    dzaehler_i=0.d0
    dnenner_i=0.d0
    dzaehler_j=0.d0
    dnenner_j=0.d0


    do a=1,mxci
      do b=1,mxcj
        c6ref=c6ab(izi,izj,a,b,1)
        if (c6ref.gt.0) then
          ! c6mem=c6ref
          cn_refi=c6ab(izi,izj,a,b,2)
          cn_refj=c6ab(izi,izj,a,b,3)
          r=(cn_refi-cni)*(cn_refi-cni)+(cn_refj-cnj)*(cn_refj-cnj)
          if (r.lt.r_save) then
            r_save=r
            c6mem=c6ref
          end if
          expterm=exp(k3*r)
          zaehler=zaehler+c6ref*expterm
          nenner=nenner+expterm
          expterm=expterm*2.d0*k3
          term=expterm*(cni-cn_refi)
          dzaehler_i=dzaehler_i+c6ref*term
          dnenner_i =dnenner_i + term

          term=expterm*(cnj-cn_refj)
          dzaehler_j=dzaehler_j+c6ref*term
          dnenner_j =dnenner_j + term
        end if
        !b
      end do
      !a
    end do

    if (nenner.gt.1.0d-99) then
      c6check=zaehler/nenner
      dc6i=((dzaehler_i*nenner)-(dnenner_i*zaehler)) &
          & /(nenner*nenner)
      dc6j=((dzaehler_j*nenner)-(dnenner_j*zaehler)) &
          & /(nenner*nenner)
    else
      c6check=c6mem
      dc6i=0.0d0
      dc6j=0.0d0
    end if
  end subroutine get_dC6_dCNij



  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! interpolate c6
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine getc6(maxc,max_elem,c6ab,mxc,iat,jat,nci,ncj,c6)
    integer maxc,max_elem
    integer iat,jat,i,j,mxc(max_elem)
    real(wp) nci,ncj,c6,c6mem
    real(wp) c6ab(max_elem,max_elem,maxc,maxc,3)
    ! the exponential is sensitive to numerics
    ! when nci or ncj is much larger than cn1/cn2
    real(wp) cn1,cn2,r,rsum,csum,tmp,tmp1
    real(wp) r_save

    c6mem=-1.d+99
    rsum=0.0d0
    csum=0.0d0
    c6 =0.0d0
    r_save=1.0d99
    do i=1,mxc(iat)
      do j=1,mxc(jat)
        c6=c6ab(iat,jat,i,j,1)
        if (c6.gt.0)then
          ! c6mem=c6
          cn1=c6ab(iat,jat,i,j,2)
          cn2=c6ab(iat,jat,i,j,3)
          ! distance
          r=(cn1-nci)**2+(cn2-ncj)**2
          if (r.lt.r_save) then
            r_save=r
            c6mem=c6
          end if
          tmp1=exp(k3*r)
          rsum=rsum+tmp1
          csum=csum+tmp1*c6
        end if
      end do
    end do

    if (rsum.gt.1.0d-99)then
      c6=csum/rsum
    else
      c6=c6mem
    end if

  end subroutine getc6

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! compute coordination numbers by adding an inverse damping function
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine ncoord(natoms,rcov,iz,xyz,cn,cn_thr)
    integer iz(*),natoms,i,max_elem
    real(wp) xyz(3,*),cn(*),rcov(94),input
    real(wp) cn_thr

    integer iat
    real(wp) dx,dy,dz,r,damp,xn,rr,rco,r2

    do i=1,natoms
      xn=0.0d0
      do iat=1,natoms
        if (iat.ne.i)then
          dx=xyz(1,iat)-xyz(1,i)
          dy=xyz(2,iat)-xyz(2,i)
          dz=xyz(3,iat)-xyz(3,i)
          r2=dx*dx+dy*dy+dz*dz
          if (r2.gt.cn_thr) cycle
          r=sqrt(r2)
          ! covalent distance in Bohr
          rco=rcov(iz(i))+rcov(iz(iat))
          rr=rco/r
          ! counting function exponential has a better long-range behavior than MH
          damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
          xn=xn+damp
        end if
      end do
      ! if (iz(i).eq.19) then
      ! write(*,*) "Input CN of Kalium"
      ! read(*,*),input
      ! cn(i)=input
      ! else
      cn(i)=xn
      ! end if
    end do

  end subroutine ncoord


  integer function lin(i1,i2)
!$acc routine  seq
    integer i1,i2,idum1,idum2
    idum1=max(i1,i2)
    idum2=min(i1,i2)
    lin=idum2+idum1*(idum1-1)/2
    return
  end function lin


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! set cut-off radii
  ! in parts due to INTEL compiler bug
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine setr0ab(max_elem,autoang,r)
    integer max_elem,i,j,k
    real(wp) r(max_elem,max_elem),autoang
    real(wp) r0ab(4465)
    r0ab( 1: 70)=(/ &
        & 2.1823, 1.8547, 1.7347, 2.9086, 2.5732, 3.4956, 2.3550 &
        &, 2.5095, 2.9802, 3.0982, 2.5141, 2.3917, 2.9977, 2.9484 &
        &, 3.2160, 2.4492, 2.2527, 3.1933, 3.0214, 2.9531, 2.9103 &
        &, 2.3667, 2.1328, 2.8784, 2.7660, 2.7776, 2.7063, 2.6225 &
        &, 2.1768, 2.0625, 2.6395, 2.6648, 2.6482, 2.5697, 2.4846 &
        &, 2.4817, 2.0646, 1.9891, 2.5086, 2.6908, 2.6233, 2.4770 &
        &, 2.3885, 2.3511, 2.2996, 1.9892, 1.9251, 2.4190, 2.5473 &
        &, 2.4994, 2.4091, 2.3176, 2.2571, 2.1946, 2.1374, 2.9898 &
        &, 2.6397, 3.6031, 3.1219, 3.7620, 3.2485, 2.9357, 2.7093 &
        &, 2.5781, 2.4839, 3.7082, 2.5129, 2.7321, 3.1052, 3.2962 &
        &/)
    r0ab( 71: 140)=(/ &
        & 3.1331, 3.2000, 2.9586, 3.0822, 2.8582, 2.7120, 3.2570 &
        &, 3.4839, 2.8766, 2.7427, 3.2776, 3.2363, 3.5929, 3.2826 &
        &, 3.0911, 2.9369, 2.9030, 2.7789, 3.3921, 3.3970, 4.0106 &
        &, 2.8884, 2.6605, 3.7513, 3.1613, 3.3605, 3.3325, 3.0991 &
        &, 2.9297, 2.8674, 2.7571, 3.8129, 3.3266, 3.7105, 3.7917 &
        &, 2.8304, 2.5538, 3.3932, 3.1193, 3.1866, 3.1245, 3.0465 &
        &, 2.8727, 2.7664, 2.6926, 3.4608, 3.2984, 3.5142, 3.5418 &
        &, 3.5017, 2.6190, 2.4797, 3.1331, 3.0540, 3.0651, 2.9879 &
        &, 2.9054, 2.8805, 2.7330, 2.6331, 3.2096, 3.5668, 3.3684 &
        &, 3.3686, 3.3180, 3.3107, 2.4757, 2.4019, 2.9789, 3.1468 &
        &/)
    r0ab( 141: 210)=(/ &
        & 2.9768, 2.8848, 2.7952, 2.7457, 2.6881, 2.5728, 3.0574 &
        &, 3.3264, 3.3562, 3.2529, 3.1916, 3.1523, 3.1046, 2.3725 &
        &, 2.3289, 2.8760, 2.9804, 2.9093, 2.8040, 2.7071, 2.6386 &
        &, 2.5720, 2.5139, 2.9517, 3.1606, 3.2085, 3.1692, 3.0982 &
        &, 3.0352, 2.9730, 2.9148, 3.2147, 2.8315, 3.8724, 3.4621 &
        &, 3.8823, 3.3760, 3.0746, 2.8817, 2.7552, 2.6605, 3.9740 &
        &, 3.6192, 3.6569, 3.9586, 3.6188, 3.3917, 3.2479, 3.1434 &
        &, 4.2411, 2.7597, 3.0588, 3.3474, 3.6214, 3.4353, 3.4729 &
        &, 3.2487, 3.3200, 3.0914, 2.9403, 3.4972, 3.7993, 3.6773 &
        &, 3.8678, 3.5808, 3.8243, 3.5826, 3.4156, 3.8765, 4.1035 &
        &/)
    r0ab( 211: 280)=(/ &
        & 2.7361, 2.9765, 3.2475, 3.5004, 3.4185, 3.4378, 3.2084 &
        &, 3.2787, 3.0604, 2.9187, 3.4037, 3.6759, 3.6586, 3.8327 &
        &, 3.5372, 3.7665, 3.5310, 3.3700, 3.7788, 3.9804, 3.8903 &
        &, 2.6832, 2.9060, 3.2613, 3.4359, 3.3538, 3.3860, 3.1550 &
        &, 3.2300, 3.0133, 2.8736, 3.4024, 3.6142, 3.5979, 3.5295 &
        &, 3.4834, 3.7140, 3.4782, 3.3170, 3.7434, 3.9623, 3.8181 &
        &, 3.7642, 2.6379, 2.8494, 3.1840, 3.4225, 3.2771, 3.3401 &
        &, 3.1072, 3.1885, 2.9714, 2.8319, 3.3315, 3.5979, 3.5256 &
        &, 3.4980, 3.4376, 3.6714, 3.4346, 3.2723, 3.6859, 3.8985 &
        &, 3.7918, 3.7372, 3.7211, 2.9230, 2.6223, 3.4161, 2.8999 &
        &/)
    r0ab( 281: 350)=(/ &
        & 3.0557, 3.3308, 3.0555, 2.8508, 2.7385, 2.6640, 3.5263 &
        &, 3.0277, 3.2990, 3.7721, 3.5017, 3.2751, 3.1368, 3.0435 &
        &, 3.7873, 3.2858, 3.2140, 3.1727, 3.2178, 3.4414, 2.5490 &
        &, 2.7623, 3.0991, 3.3252, 3.1836, 3.2428, 3.0259, 3.1225 &
        &, 2.9032, 2.7621, 3.2490, 3.5110, 3.4429, 3.3845, 3.3574 &
        &, 3.6045, 3.3658, 3.2013, 3.6110, 3.8241, 3.7090, 3.6496 &
        &, 3.6333, 3.0896, 3.5462, 2.4926, 2.7136, 3.0693, 3.2699 &
        &, 3.1272, 3.1893, 2.9658, 3.0972, 2.8778, 2.7358, 3.2206 &
        &, 3.4566, 3.3896, 3.3257, 3.2946, 3.5693, 3.3312, 3.1670 &
        &, 3.5805, 3.7711, 3.6536, 3.5927, 3.5775, 3.0411, 3.4885 &
        &/)
    r0ab( 351: 420)=(/ &
        & 3.4421, 2.4667, 2.6709, 3.0575, 3.2357, 3.0908, 3.1537 &
        &, 2.9235, 3.0669, 2.8476, 2.7054, 3.2064, 3.4519, 3.3593 &
        &, 3.2921, 3.2577, 3.2161, 3.2982, 3.1339, 3.5606, 3.7582 &
        &, 3.6432, 3.5833, 3.5691, 3.0161, 3.4812, 3.4339, 3.4327 &
        &, 2.4515, 2.6338, 3.0511, 3.2229, 3.0630, 3.1265, 2.8909 &
        &, 3.0253, 2.8184, 2.6764, 3.1968, 3.4114, 3.3492, 3.2691 &
        &, 3.2320, 3.1786, 3.2680, 3.1036, 3.5453, 3.7259, 3.6090 &
        &, 3.5473, 3.5327, 3.0018, 3.4413, 3.3907, 3.3593, 3.3462 &
        &, 2.4413, 2.6006, 3.0540, 3.1987, 3.0490, 3.1058, 2.8643 &
        &, 2.9948, 2.7908, 2.6491, 3.1950, 3.3922, 3.3316, 3.2585 &
        &/)
    r0ab( 421: 490)=(/ &
        & 3.2136, 3.1516, 3.2364, 3.0752, 3.5368, 3.7117, 3.5941 &
        &, 3.5313, 3.5164, 2.9962, 3.4225, 3.3699, 3.3370, 3.3234 &
        &, 3.3008, 2.4318, 2.5729, 3.0416, 3.1639, 3.0196, 3.0843 &
        &, 2.8413, 2.7436, 2.7608, 2.6271, 3.1811, 3.3591, 3.3045 &
        &, 3.2349, 3.1942, 3.1291, 3.2111, 3.0534, 3.5189, 3.6809 &
        &, 3.5635, 3.5001, 3.4854, 2.9857, 3.3897, 3.3363, 3.3027 &
        &, 3.2890, 3.2655, 3.2309, 2.8502, 2.6934, 3.2467, 3.1921 &
        &, 3.5663, 3.2541, 3.0571, 2.9048, 2.8657, 2.7438, 3.3547 &
        &, 3.3510, 3.9837, 3.6871, 3.4862, 3.3389, 3.2413, 3.1708 &
        &, 3.6096, 3.6280, 3.6860, 3.5568, 3.4836, 3.2868, 3.3994 &
        &/)
    r0ab( 491: 560)=(/ &
        & 3.3476, 3.3170, 3.2950, 3.2874, 3.2606, 3.9579, 2.9226 &
        &, 2.6838, 3.7867, 3.1732, 3.3872, 3.3643, 3.1267, 2.9541 &
        &, 2.8505, 2.7781, 3.8475, 3.3336, 3.7359, 3.8266, 3.5733 &
        &, 3.3959, 3.2775, 3.1915, 3.9878, 3.8816, 3.5810, 3.5364 &
        &, 3.5060, 3.8097, 3.3925, 3.3348, 3.3019, 3.2796, 3.2662 &
        &, 3.2464, 3.7136, 3.8619, 2.9140, 2.6271, 3.4771, 3.1774 &
        &, 3.2560, 3.1970, 3.1207, 2.9406, 2.8322, 2.7571, 3.5455 &
        &, 3.3514, 3.5837, 3.6177, 3.5816, 3.3902, 3.2604, 3.1652 &
        &, 3.7037, 3.6283, 3.5858, 3.5330, 3.4884, 3.5789, 3.4094 &
        &, 3.3473, 3.3118, 3.2876, 3.2707, 3.2521, 3.5570, 3.6496 &
        &/)
    r0ab( 561: 630)=(/ &
        & 3.6625, 2.7300, 2.5870, 3.2471, 3.1487, 3.1667, 3.0914 &
        &, 3.0107, 2.9812, 2.8300, 2.7284, 3.3259, 3.3182, 3.4707 &
        &, 3.4748, 3.4279, 3.4182, 3.2547, 3.1353, 3.5116, 3.9432 &
        &, 3.8828, 3.8303, 3.7880, 3.3760, 3.7218, 3.3408, 3.3059 &
        &, 3.2698, 3.2446, 3.2229, 3.4422, 3.5023, 3.5009, 3.5268 &
        &, 2.6026, 2.5355, 3.1129, 3.2863, 3.1029, 3.0108, 2.9227 &
        &, 2.8694, 2.8109, 2.6929, 3.1958, 3.4670, 3.4018, 3.3805 &
        &, 3.3218, 3.2815, 3.2346, 3.0994, 3.3937, 3.7266, 3.6697 &
        &, 3.6164, 3.5730, 3.2522, 3.5051, 3.4686, 3.4355, 3.4084 &
        &, 3.3748, 3.3496, 3.3692, 3.4052, 3.3910, 3.3849, 3.3662 &
        &/)
    r0ab( 631: 700)=(/ &
        & 2.5087, 2.4814, 3.0239, 3.1312, 3.0535, 2.9457, 2.8496 &
        &, 2.7780, 2.7828, 2.6532, 3.1063, 3.3143, 3.3549, 3.3120 &
        &, 3.2421, 3.1787, 3.1176, 3.0613, 3.3082, 3.5755, 3.5222 &
        &, 3.4678, 3.4231, 3.1684, 3.3528, 3.3162, 3.2827, 3.2527 &
        &, 3.2308, 3.2029, 3.3173, 3.3343, 3.3092, 3.2795, 3.2452 &
        &, 3.2096, 3.2893, 2.8991, 4.0388, 3.6100, 3.9388, 3.4475 &
        &, 3.1590, 2.9812, 2.8586, 2.7683, 4.1428, 3.7911, 3.8225 &
        &, 4.0372, 3.7059, 3.4935, 3.3529, 3.2492, 4.4352, 4.0826 &
        &, 3.9733, 3.9254, 3.8646, 3.9315, 3.7837, 3.7465, 3.7211 &
        &, 3.7012, 3.6893, 3.6676, 3.7736, 4.0660, 3.7926, 3.6158 &
        &/)
    r0ab( 701: 770)=(/ &
        & 3.5017, 3.4166, 4.6176, 2.8786, 3.1658, 3.5823, 3.7689 &
        &, 3.5762, 3.5789, 3.3552, 3.4004, 3.1722, 3.0212, 3.7241 &
        &, 3.9604, 3.8500, 3.9844, 3.7035, 3.9161, 3.6751, 3.5075 &
        &, 4.1151, 4.2877, 4.1579, 4.1247, 4.0617, 3.4874, 3.9848 &
        &, 3.9280, 3.9079, 3.8751, 3.8604, 3.8277, 3.8002, 3.9981 &
        &, 3.7544, 4.0371, 3.8225, 3.6718, 4.3092, 4.4764, 2.8997 &
        &, 3.0953, 3.4524, 3.6107, 3.6062, 3.5783, 3.3463, 3.3855 &
        &, 3.1746, 3.0381, 3.6019, 3.7938, 3.8697, 3.9781, 3.6877 &
        &, 3.8736, 3.6451, 3.4890, 3.9858, 4.1179, 4.0430, 3.9563 &
        &, 3.9182, 3.4002, 3.8310, 3.7716, 3.7543, 3.7203, 3.7053 &
        &/)
    r0ab( 771: 840)=(/ &
        & 3.6742, 3.8318, 3.7631, 3.7392, 3.9892, 3.7832, 3.6406 &
        &, 4.1701, 4.3016, 4.2196, 2.8535, 3.0167, 3.3978, 3.5363 &
        &, 3.5393, 3.5301, 3.2960, 3.3352, 3.1287, 2.9967, 3.6659 &
        &, 3.7239, 3.8070, 3.7165, 3.6368, 3.8162, 3.5885, 3.4336 &
        &, 3.9829, 4.0529, 3.9584, 3.9025, 3.8607, 3.3673, 3.7658 &
        &, 3.7035, 3.6866, 3.6504, 3.6339, 3.6024, 3.7708, 3.7283 &
        &, 3.6896, 3.9315, 3.7250, 3.5819, 4.1457, 4.2280, 4.1130 &
        &, 4.0597, 3.0905, 2.7998, 3.6448, 3.0739, 3.2996, 3.5262 &
        &, 3.2559, 3.0518, 2.9394, 2.8658, 3.7514, 3.2295, 3.5643 &
        &, 3.7808, 3.6931, 3.4723, 3.3357, 3.2429, 4.0280, 3.5589 &
        &/)
    r0ab( 841: 910)=(/ &
        & 3.4636, 3.4994, 3.4309, 3.6177, 3.2946, 3.2376, 3.2050 &
        &, 3.1847, 3.1715, 3.1599, 3.5555, 3.8111, 3.7693, 3.5718 &
        &, 3.4498, 3.3662, 4.1608, 3.7417, 3.6536, 3.6154, 3.8596 &
        &, 3.0301, 2.7312, 3.5821, 3.0473, 3.2137, 3.4679, 3.1975 &
        &, 2.9969, 2.8847, 2.8110, 3.6931, 3.2076, 3.4943, 3.5956 &
        &, 3.6379, 3.4190, 3.2808, 3.1860, 3.9850, 3.5105, 3.4330 &
        &, 3.3797, 3.4155, 3.6033, 3.2737, 3.2145, 3.1807, 3.1596 &
        &, 3.1461, 3.1337, 3.4812, 3.6251, 3.7152, 3.5201, 3.3966 &
        &, 3.3107, 4.1128, 3.6899, 3.6082, 3.5604, 3.7834, 3.7543 &
        &, 2.9189, 2.6777, 3.4925, 2.9648, 3.1216, 3.2940, 3.0975 &
        &/)
    r0ab( 911: 980)=(/ &
        & 2.9757, 2.8493, 2.7638, 3.6085, 3.1214, 3.4006, 3.4793 &
        &, 3.5147, 3.3806, 3.2356, 3.1335, 3.9144, 3.4183, 3.3369 &
        &, 3.2803, 3.2679, 3.4871, 3.1714, 3.1521, 3.1101, 3.0843 &
        &, 3.0670, 3.0539, 3.3890, 3.5086, 3.5895, 3.4783, 3.3484 &
        &, 3.2559, 4.0422, 3.5967, 3.5113, 3.4576, 3.6594, 3.6313 &
        &, 3.5690, 2.8578, 2.6334, 3.4673, 2.9245, 3.0732, 3.2435 &
        &, 3.0338, 2.9462, 2.8143, 2.7240, 3.5832, 3.0789, 3.3617 &
        &, 3.4246, 3.4505, 3.3443, 3.1964, 3.0913, 3.8921, 3.3713 &
        &, 3.2873, 3.2281, 3.2165, 3.4386, 3.1164, 3.1220, 3.0761 &
        &, 3.0480, 3.0295, 3.0155, 3.3495, 3.4543, 3.5260, 3.4413 &
        &/)
    r0ab( 981:1050)=(/ &
        & 3.3085, 3.2134, 4.0170, 3.5464, 3.4587, 3.4006, 3.6027 &
        &, 3.5730, 3.4945, 3.4623, 2.8240, 2.5960, 3.4635, 2.9032 &
        &, 3.0431, 3.2115, 2.9892, 2.9148, 2.7801, 2.6873, 3.5776 &
        &, 3.0568, 3.3433, 3.3949, 3.4132, 3.3116, 3.1616, 3.0548 &
        &, 3.8859, 3.3719, 3.2917, 3.2345, 3.2274, 3.4171, 3.1293 &
        &, 3.0567, 3.0565, 3.0274, 3.0087, 2.9939, 3.3293, 3.4249 &
        &, 3.4902, 3.4091, 3.2744, 3.1776, 4.0078, 3.5374, 3.4537 &
        &, 3.3956, 3.5747, 3.5430, 3.4522, 3.4160, 3.3975, 2.8004 &
        &, 2.5621, 3.4617, 2.9154, 3.0203, 3.1875, 2.9548, 2.8038 &
        &, 2.7472, 2.6530, 3.5736, 3.0584, 3.3304, 3.3748, 3.3871 &
        &/)
    r0ab(1051:1120)=(/ &
        & 3.2028, 3.1296, 3.0214, 3.8796, 3.3337, 3.2492, 3.1883 &
        &, 3.1802, 3.4050, 3.0756, 3.0478, 3.0322, 3.0323, 3.0163 &
        &, 3.0019, 3.3145, 3.4050, 3.4656, 3.3021, 3.2433, 3.1453 &
        &, 3.9991, 3.5017, 3.4141, 3.3520, 3.5583, 3.5251, 3.4243 &
        &, 3.3851, 3.3662, 3.3525, 2.7846, 2.5324, 3.4652, 2.8759 &
        &, 3.0051, 3.1692, 2.9273, 2.7615, 2.7164, 2.6212, 3.5744 &
        &, 3.0275, 3.3249, 3.3627, 3.3686, 3.1669, 3.0584, 2.9915 &
        &, 3.8773, 3.3099, 3.2231, 3.1600, 3.1520, 3.4023, 3.0426 &
        &, 3.0099, 2.9920, 2.9809, 2.9800, 2.9646, 3.3068, 3.3930 &
        &, 3.4486, 3.2682, 3.1729, 3.1168, 3.9952, 3.4796, 3.3901 &
        &/)
    r0ab(1121:1190)=(/ &
        & 3.3255, 3.5530, 3.5183, 3.4097, 3.3683, 3.3492, 3.3360 &
        &, 3.3308, 2.5424, 2.6601, 3.2555, 3.2807, 3.1384, 3.1737 &
        &, 2.9397, 2.8429, 2.8492, 2.7225, 3.3875, 3.4910, 3.4520 &
        &, 3.3608, 3.3036, 3.2345, 3.2999, 3.1487, 3.7409, 3.8392 &
        &, 3.7148, 3.6439, 3.6182, 3.1753, 3.5210, 3.4639, 3.4265 &
        &, 3.4075, 3.3828, 3.3474, 3.4071, 3.3754, 3.3646, 3.3308 &
        &, 3.4393, 3.2993, 3.8768, 3.9891, 3.8310, 3.7483, 3.3417 &
        &, 3.3019, 3.2250, 3.1832, 3.1578, 3.1564, 3.1224, 3.4620 &
        &, 2.9743, 2.8058, 3.4830, 3.3474, 3.6863, 3.3617, 3.1608 &
        &, 3.0069, 2.9640, 2.8427, 3.5885, 3.5219, 4.1314, 3.8120 &
        &/)
    r0ab(1191:1260)=(/ &
        & 3.6015, 3.4502, 3.3498, 3.2777, 3.8635, 3.8232, 3.8486 &
        &, 3.7215, 3.6487, 3.4724, 3.5627, 3.5087, 3.4757, 3.4517 &
        &, 3.4423, 3.4139, 4.1028, 3.8388, 3.6745, 3.5562, 3.4806 &
        &, 3.4272, 4.0182, 3.9991, 4.0007, 3.9282, 3.7238, 3.6498 &
        &, 3.5605, 3.5211, 3.5009, 3.4859, 3.4785, 3.5621, 4.2623 &
        &, 3.0775, 2.8275, 4.0181, 3.3385, 3.5379, 3.5036, 3.2589 &
        &, 3.0804, 3.0094, 2.9003, 4.0869, 3.5088, 3.9105, 3.9833 &
        &, 3.7176, 3.5323, 3.4102, 3.3227, 4.2702, 4.0888, 3.7560 &
        &, 3.7687, 3.6681, 3.6405, 3.5569, 3.4990, 3.4659, 3.4433 &
        &, 3.4330, 3.4092, 3.8867, 4.0190, 3.7961, 3.6412, 3.5405 &
        &/)
    r0ab(1261:1330)=(/ &
        & 3.4681, 4.3538, 4.2136, 3.9381, 3.8912, 3.9681, 3.7909 &
        &, 3.6774, 3.6262, 3.5999, 3.5823, 3.5727, 3.5419, 4.0245 &
        &, 4.1874, 3.0893, 2.7917, 3.7262, 3.3518, 3.4241, 3.5433 &
        &, 3.2773, 3.0890, 2.9775, 2.9010, 3.8048, 3.5362, 3.7746 &
        &, 3.7911, 3.7511, 3.5495, 3.4149, 3.3177, 4.0129, 3.8370 &
        &, 3.7739, 3.7125, 3.7152, 3.7701, 3.5813, 3.5187, 3.4835 &
        &, 3.4595, 3.4439, 3.4242, 3.7476, 3.8239, 3.8346, 3.6627 &
        &, 3.5479, 3.4639, 4.1026, 3.9733, 3.9292, 3.8667, 3.9513 &
        &, 3.8959, 3.7698, 3.7089, 3.6765, 3.6548, 3.6409, 3.5398 &
        &, 3.8759, 3.9804, 4.0150, 2.9091, 2.7638, 3.5066, 3.3377 &
        &/)
    r0ab(1331:1400)=(/ &
        & 3.3481, 3.2633, 3.1810, 3.1428, 2.9872, 2.8837, 3.5929 &
        &, 3.5183, 3.6729, 3.6596, 3.6082, 3.5927, 3.4224, 3.2997 &
        &, 3.8190, 4.1865, 4.1114, 4.0540, 3.6325, 3.5697, 3.5561 &
        &, 3.5259, 3.4901, 3.4552, 3.4315, 3.4091, 3.6438, 3.6879 &
        &, 3.6832, 3.7043, 3.5557, 3.4466, 3.9203, 4.2919, 4.2196 &
        &, 4.1542, 3.7573, 3.7039, 3.6546, 3.6151, 3.5293, 3.4849 &
        &, 3.4552, 3.5192, 3.7673, 3.8359, 3.8525, 3.8901, 2.7806 &
        &, 2.7209, 3.3812, 3.4958, 3.2913, 3.1888, 3.0990, 3.0394 &
        &, 2.9789, 2.8582, 3.4716, 3.6883, 3.6105, 3.5704, 3.5059 &
        &, 3.4619, 3.4138, 3.2742, 3.7080, 3.9773, 3.9010, 3.8409 &
        &/)
    r0ab(1401:1470)=(/ &
        & 3.7944, 3.4465, 3.7235, 3.6808, 3.6453, 3.6168, 3.5844 &
        &, 3.5576, 3.5772, 3.5959, 3.5768, 3.5678, 3.5486, 3.4228 &
        &, 3.8107, 4.0866, 4.0169, 3.9476, 3.6358, 3.5800, 3.5260 &
        &, 3.4838, 3.4501, 3.4204, 3.3553, 3.6487, 3.6973, 3.7398 &
        &, 3.7405, 3.7459, 3.7380, 2.6848, 2.6740, 3.2925, 3.3386 &
        &, 3.2473, 3.1284, 3.0301, 2.9531, 2.9602, 2.8272, 3.3830 &
        &, 3.5358, 3.5672, 3.5049, 3.4284, 3.3621, 3.3001, 3.2451 &
        &, 3.6209, 3.8299, 3.7543, 3.6920, 3.6436, 3.3598, 3.5701 &
        &, 3.5266, 3.4904, 3.4590, 3.4364, 3.4077, 3.5287, 3.5280 &
        &, 3.4969, 3.4650, 3.4304, 3.3963, 3.7229, 3.9402, 3.8753 &
        &/)
    r0ab(1471:1540)=(/ &
        & 3.8035, 3.5499, 3.4913, 3.4319, 3.3873, 3.3520, 3.3209 &
        &, 3.2948, 3.5052, 3.6465, 3.6696, 3.6577, 3.6388, 3.6142 &
        &, 3.5889, 3.3968, 3.0122, 4.2241, 3.7887, 4.0049, 3.5384 &
        &, 3.2698, 3.1083, 2.9917, 2.9057, 4.3340, 3.9900, 4.6588 &
        &, 4.1278, 3.8125, 3.6189, 3.4851, 3.3859, 4.6531, 4.3134 &
        &, 4.2258, 4.1309, 4.0692, 4.0944, 3.9850, 3.9416, 3.9112 &
        &, 3.8873, 3.8736, 3.8473, 4.6027, 4.1538, 3.8994, 3.7419 &
        &, 3.6356, 3.5548, 4.8353, 4.5413, 4.3891, 4.3416, 4.3243 &
        &, 4.2753, 4.2053, 4.1790, 4.1685, 4.1585, 4.1536, 4.0579 &
        &, 4.1980, 4.4564, 4.2192, 4.0528, 3.9489, 3.8642, 5.0567 &
        &/)
    r0ab(1541:1610)=(/ &
        & 3.0630, 3.3271, 4.0432, 4.0046, 4.1555, 3.7426, 3.5130 &
        &, 3.5174, 3.2884, 3.1378, 4.1894, 4.2321, 4.1725, 4.1833 &
        &, 3.8929, 4.0544, 3.8118, 3.6414, 4.6373, 4.6268, 4.4750 &
        &, 4.4134, 4.3458, 3.8582, 4.2583, 4.1898, 4.1562, 4.1191 &
        &, 4.1069, 4.0639, 4.1257, 4.1974, 3.9532, 4.1794, 3.9660 &
        &, 3.8130, 4.8160, 4.8272, 4.6294, 4.5840, 4.0770, 4.0088 &
        &, 3.9103, 3.8536, 3.8324, 3.7995, 3.7826, 4.2294, 4.3380 &
        &, 4.4352, 4.1933, 4.4580, 4.2554, 4.1072, 5.0454, 5.1814 &
        &, 3.0632, 3.2662, 3.6432, 3.8088, 3.7910, 3.7381, 3.5093 &
        &, 3.5155, 3.3047, 3.1681, 3.7871, 3.9924, 4.0637, 4.1382 &
        &/)
    r0ab(1611:1680)=(/ &
        & 3.8591, 4.0164, 3.7878, 3.6316, 4.1741, 4.3166, 4.2395 &
        &, 4.1831, 4.1107, 3.5857, 4.0270, 3.9676, 3.9463, 3.9150 &
        &, 3.9021, 3.8708, 4.0240, 4.1551, 3.9108, 4.1337, 3.9289 &
        &, 3.7873, 4.3666, 4.5080, 4.4232, 4.3155, 3.8461, 3.8007 &
        &, 3.6991, 3.6447, 3.6308, 3.5959, 3.5749, 4.0359, 4.3124 &
        &, 4.3539, 4.1122, 4.3772, 4.1785, 4.0386, 4.7004, 4.8604 &
        &, 4.6261, 2.9455, 3.2470, 3.6108, 3.8522, 3.6625, 3.6598 &
        &, 3.4411, 3.4660, 3.2415, 3.0944, 3.7514, 4.0397, 3.9231 &
        &, 4.0561, 3.7860, 3.9845, 3.7454, 3.5802, 4.1366, 4.3581 &
        &, 4.2351, 4.2011, 4.1402, 3.5381, 4.0653, 4.0093, 3.9883 &
        &/)
    r0ab(1681:1750)=(/ &
        & 3.9570, 3.9429, 3.9112, 3.8728, 4.0682, 3.8351, 4.1054 &
        &, 3.8928, 3.7445, 4.3415, 4.5497, 4.3833, 4.3122, 3.8051 &
        &, 3.7583, 3.6622, 3.6108, 3.5971, 3.5628, 3.5408, 4.0780 &
        &, 4.0727, 4.2836, 4.0553, 4.3647, 4.1622, 4.0178, 4.5802 &
        &, 4.9125, 4.5861, 4.6201, 2.9244, 3.2241, 3.5848, 3.8293 &
        &, 3.6395, 3.6400, 3.4204, 3.4499, 3.2253, 3.0779, 3.7257 &
        &, 4.0170, 3.9003, 4.0372, 3.7653, 3.9672, 3.7283, 3.5630 &
        &, 4.1092, 4.3347, 4.2117, 4.1793, 4.1179, 3.5139, 4.0426 &
        &, 3.9867, 3.9661, 3.9345, 3.9200, 3.8883, 3.8498, 4.0496 &
        &, 3.8145, 4.0881, 3.8756, 3.7271, 4.3128, 4.5242, 4.3578 &
        &/)
    r0ab(1751:1820)=(/ &
        & 4.2870, 3.7796, 3.7318, 3.6364, 3.5854, 3.5726, 3.5378 &
        &, 3.5155, 4.0527, 4.0478, 4.2630, 4.0322, 4.3449, 4.1421 &
        &, 3.9975, 4.5499, 4.8825, 4.5601, 4.5950, 4.5702, 2.9046 &
        &, 3.2044, 3.5621, 3.8078, 3.6185, 3.6220, 3.4019, 3.4359 &
        &, 3.2110, 3.0635, 3.7037, 3.9958, 3.8792, 4.0194, 3.7460 &
        &, 3.9517, 3.7128, 3.5474, 4.0872, 4.3138, 4.1906, 4.1593 &
        &, 4.0973, 3.4919, 4.0216, 3.9657, 3.9454, 3.9134, 3.8986 &
        &, 3.8669, 3.8289, 4.0323, 3.7954, 4.0725, 3.8598, 3.7113 &
        &, 4.2896, 4.5021, 4.3325, 4.2645, 3.7571, 3.7083, 3.6136 &
        &, 3.5628, 3.5507, 3.5155, 3.4929, 4.0297, 4.0234, 4.2442 &
        &/)
    r0ab(1821:1890)=(/ &
        & 4.0112, 4.3274, 4.1240, 3.9793, 4.5257, 4.8568, 4.5353 &
        &, 4.5733, 4.5485, 4.5271, 2.8878, 3.1890, 3.5412, 3.7908 &
        &, 3.5974, 3.6078, 3.3871, 3.4243, 3.1992, 3.0513, 3.6831 &
        &, 3.9784, 3.8579, 4.0049, 3.7304, 3.9392, 3.7002, 3.5347 &
        &, 4.0657, 4.2955, 4.1705, 4.1424, 4.0800, 3.4717, 4.0043 &
        &, 3.9485, 3.9286, 3.8965, 3.8815, 3.8500, 3.8073, 4.0180 &
        &, 3.7796, 4.0598, 3.8470, 3.6983, 4.2678, 4.4830, 4.3132 &
        &, 4.2444, 3.7370, 3.6876, 3.5935, 3.5428, 3.5314, 3.4958 &
        &, 3.4730, 4.0117, 4.0043, 4.2287, 3.9939, 4.3134, 4.1096 &
        &, 3.9646, 4.5032, 4.8356, 4.5156, 4.5544, 4.5297, 4.5083 &
        &/)
    r0ab(1891:1960)=(/ &
        & 4.4896, 2.8709, 3.1737, 3.5199, 3.7734, 3.5802, 3.5934 &
        &, 3.3724, 3.4128, 3.1877, 3.0396, 3.6624, 3.9608, 3.8397 &
        &, 3.9893, 3.7145, 3.9266, 3.6877, 3.5222, 4.0448, 4.2771 &
        &, 4.1523, 4.1247, 4.0626, 3.4530, 3.9866, 3.9310, 3.9115 &
        &, 3.8792, 3.8641, 3.8326, 3.7892, 4.0025, 3.7636, 4.0471 &
        &, 3.8343, 3.6854, 4.2464, 4.4635, 4.2939, 4.2252, 3.7169 &
        &, 3.6675, 3.5739, 3.5235, 3.5126, 3.4768, 3.4537, 3.9932 &
        &, 3.9854, 4.2123, 3.9765, 4.2992, 4.0951, 3.9500, 4.4811 &
        &, 4.8135, 4.4959, 4.5351, 4.5105, 4.4891, 4.4705, 4.4515 &
        &, 2.8568, 3.1608, 3.5050, 3.7598, 3.5665, 3.5803, 3.3601 &
        &/)
    r0ab(1961:2030)=(/ &
        & 3.4031, 3.1779, 3.0296, 3.6479, 3.9471, 3.8262, 3.9773 &
        &, 3.7015, 3.9162, 3.6771, 3.5115, 4.0306, 4.2634, 4.1385 &
        &, 4.1116, 4.0489, 3.4366, 3.9732, 3.9176, 3.8983, 3.8659 &
        &, 3.8507, 3.8191, 3.7757, 3.9907, 3.7506, 4.0365, 3.8235 &
        &, 3.6745, 4.2314, 4.4490, 4.2792, 4.2105, 3.7003, 3.6510 &
        &, 3.5578, 3.5075, 3.4971, 3.4609, 3.4377, 3.9788, 3.9712 &
        &, 4.1997, 3.9624, 4.2877, 4.0831, 3.9378, 4.4655, 4.7974 &
        &, 4.4813, 4.5209, 4.4964, 4.4750, 4.4565, 4.4375, 4.4234 &
        &, 2.6798, 3.0151, 3.2586, 3.5292, 3.5391, 3.4902, 3.2887 &
        &, 3.3322, 3.1228, 2.9888, 3.4012, 3.7145, 3.7830, 3.6665 &
        &/)
    r0ab(2031:2100)=(/ &
        & 3.5898, 3.8077, 3.5810, 3.4265, 3.7726, 4.0307, 3.9763 &
        &, 3.8890, 3.8489, 3.2706, 3.7595, 3.6984, 3.6772, 3.6428 &
        &, 3.6243, 3.5951, 3.7497, 3.6775, 3.6364, 3.9203, 3.7157 &
        &, 3.5746, 3.9494, 4.2076, 4.1563, 4.0508, 3.5329, 3.4780 &
        &, 3.3731, 3.3126, 3.2846, 3.2426, 3.2135, 3.7491, 3.9006 &
        &, 3.8332, 3.8029, 4.1436, 3.9407, 3.7998, 4.1663, 4.5309 &
        &, 4.3481, 4.2911, 4.2671, 4.2415, 4.2230, 4.2047, 4.1908 &
        &, 4.1243, 2.5189, 2.9703, 3.3063, 3.6235, 3.4517, 3.3989 &
        &, 3.2107, 3.2434, 3.0094, 2.8580, 3.4253, 3.8157, 3.7258 &
        &, 3.6132, 3.5297, 3.7566, 3.5095, 3.3368, 3.7890, 4.1298 &
        &/)
    r0ab(2101:2170)=(/ &
        & 4.0190, 3.9573, 3.9237, 3.2677, 3.8480, 3.8157, 3.7656 &
        &, 3.7317, 3.7126, 3.6814, 3.6793, 3.6218, 3.5788, 3.8763 &
        &, 3.6572, 3.5022, 3.9737, 4.3255, 4.1828, 4.1158, 3.5078 &
        &, 3.4595, 3.3600, 3.3088, 3.2575, 3.2164, 3.1856, 3.8522 &
        &, 3.8665, 3.8075, 3.7772, 4.1391, 3.9296, 3.7772, 4.2134 &
        &, 4.7308, 4.3787, 4.3894, 4.3649, 4.3441, 4.3257, 4.3073 &
        &, 4.2941, 4.1252, 4.2427, 3.0481, 2.9584, 3.6919, 3.5990 &
        &, 3.8881, 3.4209, 3.1606, 3.1938, 2.9975, 2.8646, 3.8138 &
        &, 3.7935, 3.7081, 3.9155, 3.5910, 3.4808, 3.4886, 3.3397 &
        &, 4.1336, 4.1122, 3.9888, 3.9543, 3.8917, 3.5894, 3.8131 &
        &/)
    r0ab(2171:2240)=(/ &
        & 3.7635, 3.7419, 3.7071, 3.6880, 3.6574, 3.6546, 3.9375 &
        &, 3.6579, 3.5870, 3.6361, 3.5039, 4.3149, 4.2978, 4.1321 &
        &, 4.1298, 3.8164, 3.7680, 3.7154, 3.6858, 3.6709, 3.6666 &
        &, 3.6517, 3.8174, 3.8608, 4.1805, 3.9102, 3.8394, 3.8968 &
        &, 3.7673, 4.5274, 4.6682, 4.3344, 4.3639, 4.3384, 4.3162 &
        &, 4.2972, 4.2779, 4.2636, 4.0253, 4.1168, 4.1541, 2.8136 &
        &, 3.0951, 3.4635, 3.6875, 3.4987, 3.5183, 3.2937, 3.3580 &
        &, 3.1325, 2.9832, 3.6078, 3.8757, 3.7616, 3.9222, 3.6370 &
        &, 3.8647, 3.6256, 3.4595, 3.9874, 4.1938, 4.0679, 4.0430 &
        &, 3.9781, 3.3886, 3.9008, 3.8463, 3.8288, 3.7950, 3.7790 &
        &/)
    r0ab(2241:2310)=(/ &
        & 3.7472, 3.7117, 3.9371, 3.6873, 3.9846, 3.7709, 3.6210 &
        &, 4.1812, 4.3750, 4.2044, 4.1340, 3.6459, 3.5929, 3.5036 &
        &, 3.4577, 3.4528, 3.4146, 3.3904, 3.9014, 3.9031, 4.1443 &
        &, 3.8961, 4.2295, 4.0227, 3.8763, 4.4086, 4.7097, 4.4064 &
        &, 4.4488, 4.4243, 4.4029, 4.3842, 4.3655, 4.3514, 4.1162 &
        &, 4.2205, 4.1953, 4.2794, 2.8032, 3.0805, 3.4519, 3.6700 &
        &, 3.4827, 3.5050, 3.2799, 3.3482, 3.1233, 2.9747, 3.5971 &
        &, 3.8586, 3.7461, 3.9100, 3.6228, 3.8535, 3.6147, 3.4490 &
        &, 3.9764, 4.1773, 4.0511, 4.0270, 3.9614, 3.3754, 3.8836 &
        &, 3.8291, 3.8121, 3.7780, 3.7619, 3.7300, 3.6965, 3.9253 &
        &/)
    r0ab(2311:2380)=(/ &
        & 3.6734, 3.9733, 3.7597, 3.6099, 4.1683, 4.3572, 4.1862 &
        &, 4.1153, 3.6312, 3.5772, 3.4881, 3.4429, 3.4395, 3.4009 &
        &, 3.3766, 3.8827, 3.8868, 4.1316, 3.8807, 4.2164, 4.0092 &
        &, 3.8627, 4.3936, 4.6871, 4.3882, 4.4316, 4.4073, 4.3858 &
        &, 4.3672, 4.3485, 4.3344, 4.0984, 4.2036, 4.1791, 4.2622 &
        &, 4.2450, 2.7967, 3.0689, 3.4445, 3.6581, 3.4717, 3.4951 &
        &, 3.2694, 3.3397, 3.1147, 2.9661, 3.5898, 3.8468, 3.7358 &
        &, 3.9014, 3.6129, 3.8443, 3.6054, 3.4396, 3.9683, 4.1656 &
        &, 4.0394, 4.0158, 3.9498, 3.3677, 3.8718, 3.8164, 3.8005 &
        &, 3.7662, 3.7500, 3.7181, 3.6863, 3.9170, 3.6637, 3.9641 &
        &/)
    r0ab(2381:2450)=(/ &
        & 3.7503, 3.6004, 4.1590, 4.3448, 4.1739, 4.1029, 3.6224 &
        &, 3.5677, 3.4785, 3.4314, 3.4313, 3.3923, 3.3680, 3.8698 &
        &, 3.8758, 4.1229, 3.8704, 4.2063, 3.9987, 3.8519, 4.3832 &
        &, 4.6728, 4.3759, 4.4195, 4.3952, 4.3737, 4.3551, 4.3364 &
        &, 4.3223, 4.0861, 4.1911, 4.1676, 4.2501, 4.2329, 4.2208 &
        &, 2.7897, 3.0636, 3.4344, 3.6480, 3.4626, 3.4892, 3.2626 &
        &, 3.3344, 3.1088, 2.9597, 3.5804, 3.8359, 3.7251, 3.8940 &
        &, 3.6047, 3.8375, 3.5990, 3.4329, 3.9597, 4.1542, 4.0278 &
        &, 4.0048, 3.9390, 3.3571, 3.8608, 3.8056, 3.7899, 3.7560 &
        &, 3.7400, 3.7081, 3.6758, 3.9095, 3.6552, 3.9572, 3.7436 &
        &/)
    r0ab(2451:2520)=(/ &
        & 3.5933, 4.1508, 4.3337, 4.1624, 4.0916, 3.6126, 3.5582 &
        &, 3.4684, 3.4212, 3.4207, 3.3829, 3.3586, 3.8604, 3.8658 &
        &, 4.1156, 3.8620, 4.1994, 3.9917, 3.8446, 4.3750, 4.6617 &
        &, 4.3644, 4.4083, 4.3840, 4.3625, 4.3439, 4.3253, 4.3112 &
        &, 4.0745, 4.1807, 4.1578, 4.2390, 4.2218, 4.2097, 4.1986 &
        &, 2.8395, 3.0081, 3.3171, 3.4878, 3.5360, 3.5145, 3.2809 &
        &, 3.3307, 3.1260, 2.9940, 3.4741, 3.6675, 3.7832, 3.6787 &
        &, 3.6156, 3.8041, 3.5813, 3.4301, 3.8480, 3.9849, 3.9314 &
        &, 3.8405, 3.8029, 3.2962, 3.7104, 3.6515, 3.6378, 3.6020 &
        &, 3.5849, 3.5550, 3.7494, 3.6893, 3.6666, 3.9170, 3.7150 &
        &/)
    r0ab(2521:2590)=(/ &
        & 3.5760, 4.0268, 4.1596, 4.1107, 3.9995, 3.5574, 3.5103 &
        &, 3.4163, 3.3655, 3.3677, 3.3243, 3.2975, 3.7071, 3.9047 &
        &, 3.8514, 3.8422, 3.8022, 3.9323, 3.7932, 4.2343, 4.4583 &
        &, 4.3115, 4.2457, 4.2213, 4.1945, 4.1756, 4.1569, 4.1424 &
        &, 4.0620, 4.0494, 3.9953, 4.0694, 4.0516, 4.0396, 4.0280 &
        &, 4.0130, 2.9007, 2.9674, 3.8174, 3.5856, 3.6486, 3.5339 &
        &, 3.2832, 3.3154, 3.1144, 2.9866, 3.9618, 3.8430, 3.9980 &
        &, 3.8134, 3.6652, 3.7985, 3.5756, 3.4207, 4.4061, 4.2817 &
        &, 4.1477, 4.0616, 3.9979, 3.6492, 3.8833, 3.8027, 3.7660 &
        &, 3.7183, 3.6954, 3.6525, 3.9669, 3.8371, 3.7325, 3.9160 &
        &/)
    r0ab(2591:2660)=(/ &
        & 3.7156, 3.5714, 4.6036, 4.4620, 4.3092, 4.2122, 3.8478 &
        &, 3.7572, 3.6597, 3.5969, 3.5575, 3.5386, 3.5153, 3.7818 &
        &, 4.1335, 4.0153, 3.9177, 3.8603, 3.9365, 3.7906, 4.7936 &
        &, 4.7410, 4.5461, 4.5662, 4.5340, 4.5059, 4.4832, 4.4604 &
        &, 4.4429, 4.2346, 4.4204, 4.3119, 4.3450, 4.3193, 4.3035 &
        &, 4.2933, 4.1582, 4.2450, 2.8559, 2.9050, 3.8325, 3.5442 &
        &, 3.5077, 3.4905, 3.2396, 3.2720, 3.0726, 2.9467, 3.9644 &
        &, 3.8050, 3.8981, 3.7762, 3.6216, 3.7531, 3.5297, 3.3742 &
        &, 4.3814, 4.2818, 4.1026, 4.0294, 3.9640, 3.6208, 3.8464 &
        &, 3.7648, 3.7281, 3.6790, 3.6542, 3.6117, 3.8650, 3.8010 &
        &/)
    r0ab(2661:2730)=(/ &
        & 3.6894, 3.8713, 3.6699, 3.5244, 4.5151, 4.4517, 4.2538 &
        &, 4.1483, 3.8641, 3.7244, 3.6243, 3.5589, 3.5172, 3.4973 &
        &, 3.4715, 3.7340, 4.0316, 3.9958, 3.8687, 3.8115, 3.8862 &
        &, 3.7379, 4.7091, 4.7156, 4.5199, 4.5542, 4.5230, 4.4959 &
        &, 4.4750, 4.4529, 4.4361, 4.1774, 4.3774, 4.2963, 4.3406 &
        &, 4.3159, 4.3006, 4.2910, 4.1008, 4.1568, 4.0980, 2.8110 &
        &, 2.8520, 3.7480, 3.5105, 3.4346, 3.3461, 3.1971, 3.2326 &
        &, 3.0329, 2.9070, 3.8823, 3.7928, 3.8264, 3.7006, 3.5797 &
        &, 3.7141, 3.4894, 3.3326, 4.3048, 4.2217, 4.0786, 3.9900 &
        &, 3.9357, 3.6331, 3.8333, 3.7317, 3.6957, 3.6460, 3.6197 &
        &/)
    r0ab(2731:2800)=(/ &
        & 3.5779, 3.7909, 3.7257, 3.6476, 3.5729, 3.6304, 3.4834 &
        &, 4.4368, 4.3921, 4.2207, 4.1133, 3.8067, 3.7421, 3.6140 &
        &, 3.5491, 3.5077, 3.4887, 3.4623, 3.6956, 3.9568, 3.8976 &
        &, 3.8240, 3.7684, 3.8451, 3.6949, 4.6318, 4.6559, 4.4533 &
        &, 4.4956, 4.4641, 4.4366, 4.4155, 4.3936, 4.3764, 4.1302 &
        &, 4.3398, 4.2283, 4.2796, 4.2547, 4.2391, 4.2296, 4.0699 &
        &, 4.1083, 4.0319, 3.9855, 2.7676, 2.8078, 3.6725, 3.4804 &
        &, 3.3775, 3.2411, 3.1581, 3.1983, 2.9973, 2.8705, 3.8070 &
        &, 3.7392, 3.7668, 3.6263, 3.5402, 3.6807, 3.4545, 3.2962 &
        &, 4.2283, 4.1698, 4.0240, 3.9341, 3.8711, 3.5489, 3.7798 &
        &/)
    r0ab(2801:2870)=(/ &
        & 3.7000, 3.6654, 3.6154, 3.5882, 3.5472, 3.7289, 3.6510 &
        &, 3.6078, 3.5355, 3.5963, 3.4480, 4.3587, 4.3390, 4.1635 &
        &, 4.0536, 3.7193, 3.6529, 3.5512, 3.4837, 3.4400, 3.4191 &
        &, 3.3891, 3.6622, 3.8934, 3.8235, 3.7823, 3.7292, 3.8106 &
        &, 3.6589, 4.5535, 4.6013, 4.3961, 4.4423, 4.4109, 4.3835 &
        &, 4.3625, 4.3407, 4.3237, 4.0863, 4.2835, 4.1675, 4.2272 &
        &, 4.2025, 4.1869, 4.1774, 4.0126, 4.0460, 3.9815, 3.9340 &
        &, 3.8955, 2.6912, 2.7604, 3.6037, 3.4194, 3.3094, 3.1710 &
        &, 3.0862, 3.1789, 2.9738, 2.8427, 3.7378, 3.6742, 3.6928 &
        &, 3.5512, 3.4614, 3.4087, 3.4201, 3.2607, 4.1527, 4.0977 &
        &/)
    r0ab(2871:2940)=(/ &
        & 3.9523, 3.8628, 3.8002, 3.4759, 3.7102, 3.6466, 3.6106 &
        &, 3.5580, 3.5282, 3.4878, 3.6547, 3.5763, 3.5289, 3.5086 &
        &, 3.5593, 3.4099, 4.2788, 4.2624, 4.0873, 3.9770, 3.6407 &
        &, 3.5743, 3.5178, 3.4753, 3.3931, 3.3694, 3.3339, 3.6002 &
        &, 3.8164, 3.7478, 3.7028, 3.6952, 3.7669, 3.6137, 4.4698 &
        &, 4.5488, 4.3168, 4.3646, 4.3338, 4.3067, 4.2860, 4.2645 &
        &, 4.2478, 4.0067, 4.2349, 4.0958, 4.1543, 4.1302, 4.1141 &
        &, 4.1048, 3.9410, 3.9595, 3.8941, 3.8465, 3.8089, 3.7490 &
        &, 2.7895, 2.5849, 3.6484, 3.0162, 3.1267, 3.2125, 3.0043 &
        &, 2.9572, 2.8197, 2.7261, 3.7701, 3.2446, 3.5239, 3.4696 &
        &/)
    r0ab(2941:3010)=(/ &
        & 3.4261, 3.3508, 3.1968, 3.0848, 4.1496, 3.6598, 3.5111 &
        &, 3.4199, 3.3809, 3.5382, 3.2572, 3.2100, 3.1917, 3.1519 &
        &, 3.1198, 3.1005, 3.5071, 3.5086, 3.5073, 3.4509, 3.3120 &
        &, 3.2082, 4.2611, 3.8117, 3.6988, 3.5646, 3.6925, 3.6295 &
        &, 3.5383, 3.4910, 3.4625, 3.4233, 3.4007, 3.2329, 3.6723 &
        &, 3.6845, 3.6876, 3.6197, 3.4799, 3.3737, 4.4341, 4.0525 &
        &, 3.9011, 3.8945, 3.8635, 3.8368, 3.8153, 3.7936, 3.7758 &
        &, 3.4944, 3.4873, 3.9040, 3.7110, 3.6922, 3.6799, 3.6724 &
        &, 3.5622, 3.6081, 3.5426, 3.4922, 3.4498, 3.3984, 3.4456 &
        &, 2.7522, 2.5524, 3.5742, 2.9508, 3.0751, 3.0158, 2.9644 &
        &/)
    r0ab(3011:3080)=(/ &
        & 2.8338, 2.7891, 2.6933, 3.6926, 3.1814, 3.4528, 3.4186 &
        &, 3.3836, 3.2213, 3.1626, 3.0507, 4.0548, 3.5312, 3.4244 &
        &, 3.3409, 3.2810, 3.4782, 3.1905, 3.1494, 3.1221, 3.1128 &
        &, 3.0853, 3.0384, 3.4366, 3.4562, 3.4638, 3.3211, 3.2762 &
        &, 3.1730, 4.1632, 3.6825, 3.5822, 3.4870, 3.6325, 3.5740 &
        &, 3.4733, 3.4247, 3.3969, 3.3764, 3.3525, 3.1984, 3.5989 &
        &, 3.6299, 3.6433, 3.4937, 3.4417, 3.3365, 4.3304, 3.9242 &
        &, 3.7793, 3.7623, 3.7327, 3.7071, 3.6860, 3.6650, 3.6476 &
        &, 3.3849, 3.3534, 3.8216, 3.5870, 3.5695, 3.5584, 3.5508 &
        &, 3.4856, 3.5523, 3.4934, 3.4464, 3.4055, 3.3551, 3.3888 &
        &/)
    r0ab(3081:3150)=(/ &
        & 3.3525, 2.7202, 2.5183, 3.4947, 2.8731, 3.0198, 3.1457 &
        &, 2.9276, 2.7826, 2.7574, 2.6606, 3.6090, 3.0581, 3.3747 &
        &, 3.3677, 3.3450, 3.1651, 3.1259, 3.0147, 3.9498, 3.3857 &
        &, 3.2917, 3.2154, 3.1604, 3.4174, 3.0735, 3.0342, 3.0096 &
        &, 3.0136, 2.9855, 2.9680, 3.3604, 3.4037, 3.4243, 3.2633 &
        &, 3.1810, 3.1351, 4.0557, 3.5368, 3.4526, 3.3699, 3.5707 &
        &, 3.5184, 3.4085, 3.3595, 3.3333, 3.3143, 3.3041, 3.1094 &
        &, 3.5193, 3.5745, 3.6025, 3.4338, 3.3448, 3.2952, 4.2158 &
        &, 3.7802, 3.6431, 3.6129, 3.5853, 3.5610, 3.5406, 3.5204 &
        &, 3.5036, 3.2679, 3.2162, 3.7068, 3.4483, 3.4323, 3.4221 &
        &/)
    r0ab(3151:3220)=(/ &
        & 3.4138, 3.3652, 3.4576, 3.4053, 3.3618, 3.3224, 3.2711 &
        &, 3.3326, 3.2950, 3.2564, 2.5315, 2.6104, 3.2734, 3.2299 &
        &, 3.1090, 2.9942, 2.9159, 2.8324, 2.8350, 2.7216, 3.3994 &
        &, 3.4475, 3.4354, 3.3438, 3.2807, 3.2169, 3.2677, 3.1296 &
        &, 3.7493, 3.8075, 3.6846, 3.6104, 3.5577, 3.2052, 3.4803 &
        &, 3.4236, 3.3845, 3.3640, 3.3365, 3.3010, 3.3938, 3.3624 &
        &, 3.3440, 3.3132, 3.4035, 3.2754, 3.8701, 3.9523, 3.8018 &
        &, 3.7149, 3.3673, 3.3199, 3.2483, 3.2069, 3.1793, 3.1558 &
        &, 3.1395, 3.4097, 3.5410, 3.5228, 3.5116, 3.4921, 3.4781 &
        &, 3.4690, 4.0420, 4.1759, 4.0078, 4.0450, 4.0189, 3.9952 &
        &/)
    r0ab(3221:3290)=(/ &
        & 3.9770, 3.9583, 3.9434, 3.7217, 3.8228, 3.7826, 3.8640 &
        &, 3.8446, 3.8314, 3.8225, 3.6817, 3.7068, 3.6555, 3.6159 &
        &, 3.5831, 3.5257, 3.2133, 3.1689, 3.1196, 3.3599, 2.9852 &
        &, 2.7881, 3.5284, 3.3493, 3.6958, 3.3642, 3.1568, 3.0055 &
        &, 2.9558, 2.8393, 3.6287, 3.5283, 4.1511, 3.8259, 3.6066 &
        &, 3.4527, 3.3480, 3.2713, 3.9037, 3.8361, 3.8579, 3.7311 &
        &, 3.6575, 3.5176, 3.5693, 3.5157, 3.4814, 3.4559, 3.4445 &
        &, 3.4160, 4.1231, 3.8543, 3.6816, 3.5602, 3.4798, 3.4208 &
        &, 4.0542, 4.0139, 4.0165, 3.9412, 3.7698, 3.6915, 3.6043 &
        &, 3.5639, 3.5416, 3.5247, 3.5153, 3.5654, 4.2862, 4.0437 &
        &/)
    r0ab(3291:3360)=(/ &
        & 3.8871, 3.7741, 3.6985, 3.6413, 4.2345, 4.3663, 4.3257 &
        &, 4.0869, 4.0612, 4.0364, 4.0170, 3.9978, 3.9834, 3.9137 &
        &, 3.8825, 3.8758, 3.9143, 3.8976, 3.8864, 3.8768, 3.9190 &
        &, 4.1613, 4.0566, 3.9784, 3.9116, 3.8326, 3.7122, 3.6378 &
        &, 3.5576, 3.5457, 4.3127, 3.1160, 2.8482, 4.0739, 3.3599 &
        &, 3.5698, 3.5366, 3.2854, 3.1039, 2.9953, 2.9192, 4.1432 &
        &, 3.5320, 3.9478, 4.0231, 3.7509, 3.5604, 3.4340, 3.3426 &
        &, 4.3328, 3.8288, 3.7822, 3.7909, 3.6907, 3.6864, 3.5793 &
        &, 3.5221, 3.4883, 3.4649, 3.4514, 3.4301, 3.9256, 4.0596 &
        &, 3.8307, 3.6702, 3.5651, 3.4884, 4.4182, 4.2516, 3.9687 &
        &/)
    r0ab(3361:3430)=(/ &
        & 3.9186, 3.9485, 3.8370, 3.7255, 3.6744, 3.6476, 3.6295 &
        &, 3.6193, 3.5659, 4.0663, 4.2309, 4.0183, 3.8680, 3.7672 &
        &, 3.6923, 4.5240, 4.4834, 4.1570, 4.3204, 4.2993, 4.2804 &
        &, 4.2647, 4.2481, 4.2354, 3.8626, 3.8448, 4.2267, 4.1799 &
        &, 4.1670, 3.8738, 3.8643, 3.8796, 4.0575, 4.0354, 3.9365 &
        &, 3.8611, 3.7847, 3.7388, 3.6826, 3.6251, 3.5492, 4.0889 &
        &, 4.2764, 3.1416, 2.8325, 3.7735, 3.3787, 3.4632, 3.5923 &
        &, 3.3214, 3.1285, 3.0147, 2.9366, 3.8527, 3.5602, 3.8131 &
        &, 3.8349, 3.7995, 3.5919, 3.4539, 3.3540, 4.0654, 3.8603 &
        &, 3.7972, 3.7358, 3.7392, 3.8157, 3.6055, 3.5438, 3.5089 &
        &/)
    r0ab(3431:3500)=(/ &
        & 3.4853, 3.4698, 3.4508, 3.7882, 3.8682, 3.8837, 3.7055 &
        &, 3.5870, 3.5000, 4.1573, 4.0005, 3.9568, 3.8936, 3.9990 &
        &, 3.9433, 3.8172, 3.7566, 3.7246, 3.7033, 3.6900, 3.5697 &
        &, 3.9183, 4.0262, 4.0659, 3.8969, 3.7809, 3.6949, 4.2765 &
        &, 4.2312, 4.1401, 4.0815, 4.0580, 4.0369, 4.0194, 4.0017 &
        &, 3.9874, 3.8312, 3.8120, 3.9454, 3.9210, 3.9055, 3.8951 &
        &, 3.8866, 3.8689, 3.9603, 3.9109, 3.9122, 3.8233, 3.7438 &
        &, 3.7436, 3.6981, 3.6555, 3.5452, 3.9327, 4.0658, 4.1175 &
        &, 2.9664, 2.8209, 3.5547, 3.3796, 3.3985, 3.3164, 3.2364 &
        &, 3.1956, 3.0370, 2.9313, 3.6425, 3.5565, 3.7209, 3.7108 &
        &/)
    r0ab(3501:3570)=(/ &
        & 3.6639, 3.6484, 3.4745, 3.3492, 3.8755, 4.2457, 3.7758 &
        &, 3.7161, 3.6693, 3.6155, 3.5941, 3.5643, 3.5292, 3.4950 &
        &, 3.4720, 3.4503, 3.6936, 3.7392, 3.7388, 3.7602, 3.6078 &
        &, 3.4960, 3.9800, 4.3518, 4.2802, 3.8580, 3.8056, 3.7527 &
        &, 3.7019, 3.6615, 3.5768, 3.5330, 3.5038, 3.5639, 3.8192 &
        &, 3.8883, 3.9092, 3.9478, 3.7995, 3.6896, 4.1165, 4.5232 &
        &, 4.4357, 4.4226, 4.4031, 4.3860, 4.3721, 4.3580, 4.3466 &
        &, 4.2036, 4.2037, 3.8867, 4.2895, 4.2766, 4.2662, 4.2598 &
        &, 3.8408, 3.9169, 3.8681, 3.8250, 3.7855, 3.7501, 3.6753 &
        &, 3.5499, 3.4872, 3.5401, 3.8288, 3.9217, 3.9538, 4.0054 &
        &/)
    r0ab(3571:3640)=(/ &
        & 2.8388, 2.7890, 3.4329, 3.5593, 3.3488, 3.2486, 3.1615 &
        &, 3.1000, 3.0394, 2.9165, 3.5267, 3.7479, 3.6650, 3.6263 &
        &, 3.5658, 3.5224, 3.4762, 3.3342, 3.7738, 4.0333, 3.9568 &
        &, 3.8975, 3.8521, 3.4929, 3.7830, 3.7409, 3.7062, 3.6786 &
        &, 3.6471, 3.6208, 3.6337, 3.6519, 3.6363, 3.6278, 3.6110 &
        &, 3.4825, 3.8795, 4.1448, 4.0736, 4.0045, 3.6843, 3.6291 &
        &, 3.5741, 3.5312, 3.4974, 3.4472, 3.4034, 3.7131, 3.7557 &
        &, 3.7966, 3.8005, 3.8068, 3.8015, 3.6747, 4.0222, 4.3207 &
        &, 4.2347, 4.2191, 4.1990, 4.1811, 4.1666, 4.1521, 4.1401 &
        &, 3.9970, 3.9943, 3.9592, 4.0800, 4.0664, 4.0559, 4.0488 &
        &/)
    r0ab(3641:3710)=(/ &
        & 3.9882, 4.0035, 3.9539, 3.9138, 3.8798, 3.8355, 3.5359 &
        &, 3.4954, 3.3962, 3.5339, 3.7595, 3.8250, 3.8408, 3.8600 &
        &, 3.8644, 2.7412, 2.7489, 3.3374, 3.3950, 3.3076, 3.1910 &
        &, 3.0961, 3.0175, 3.0280, 2.8929, 3.4328, 3.5883, 3.6227 &
        &, 3.5616, 3.4894, 3.4241, 3.3641, 3.3120, 3.6815, 3.8789 &
        &, 3.8031, 3.7413, 3.6939, 3.4010, 3.6225, 3.5797, 3.5443 &
        &, 3.5139, 3.4923, 3.4642, 3.5860, 3.5849, 3.5570, 3.5257 &
        &, 3.4936, 3.4628, 3.7874, 3.9916, 3.9249, 3.8530, 3.5932 &
        &, 3.5355, 3.4757, 3.4306, 3.3953, 3.3646, 3.3390, 3.5637 &
        &, 3.7053, 3.7266, 3.7177, 3.6996, 3.6775, 3.6558, 3.9331 &
        &/)
    r0ab(3711:3780)=(/ &
        & 4.1655, 4.0879, 4.0681, 4.0479, 4.0299, 4.0152, 4.0006 &
        &, 3.9883, 3.8500, 3.8359, 3.8249, 3.9269, 3.9133, 3.9025 &
        &, 3.8948, 3.8422, 3.8509, 3.7990, 3.7570, 3.7219, 3.6762 &
        &, 3.4260, 3.3866, 3.3425, 3.5294, 3.7022, 3.7497, 3.7542 &
        &, 3.7494, 3.7370, 3.7216, 3.4155, 3.0522, 4.2541, 3.8218 &
        &, 4.0438, 3.5875, 3.3286, 3.1682, 3.0566, 2.9746, 4.3627 &
        &, 4.0249, 4.6947, 4.1718, 3.8639, 3.6735, 3.5435, 3.4479 &
        &, 4.6806, 4.3485, 4.2668, 4.1690, 4.1061, 4.1245, 4.0206 &
        &, 3.9765, 3.9458, 3.9217, 3.9075, 3.8813, 3.9947, 4.1989 &
        &, 3.9507, 3.7960, 3.6925, 3.6150, 4.8535, 4.5642, 4.4134 &
        &/)
    r0ab(3781:3850)=(/ &
        & 4.3688, 4.3396, 4.2879, 4.2166, 4.1888, 4.1768, 4.1660 &
        &, 4.1608, 4.0745, 4.2289, 4.4863, 4.2513, 4.0897, 3.9876 &
        &, 3.9061, 5.0690, 5.0446, 4.6186, 4.6078, 4.5780, 4.5538 &
        &, 4.5319, 4.5101, 4.4945, 4.1912, 4.2315, 4.5534, 4.4373 &
        &, 4.4224, 4.4120, 4.4040, 4.2634, 4.7770, 4.6890, 4.6107 &
        &, 4.5331, 4.4496, 4.4082, 4.3095, 4.2023, 4.0501, 4.2595 &
        &, 4.5497, 4.3056, 4.1506, 4.0574, 3.9725, 5.0796, 3.0548 &
        &, 3.3206, 3.8132, 3.9720, 3.7675, 3.7351, 3.5167, 3.5274 &
        &, 3.3085, 3.1653, 3.9500, 4.1730, 4.0613, 4.1493, 3.8823 &
        &, 4.0537, 3.8200, 3.6582, 4.3422, 4.5111, 4.3795, 4.3362 &
        &/)
    r0ab(3851:3920)=(/ &
        & 4.2751, 3.7103, 4.1973, 4.1385, 4.1129, 4.0800, 4.0647 &
        &, 4.0308, 4.0096, 4.1619, 3.9360, 4.1766, 3.9705, 3.8262 &
        &, 4.5348, 4.7025, 4.5268, 4.5076, 3.9562, 3.9065, 3.8119 &
        &, 3.7605, 3.7447, 3.7119, 3.6916, 4.1950, 4.2110, 4.3843 &
        &, 4.1631, 4.4427, 4.2463, 4.1054, 4.7693, 5.0649, 4.7365 &
        &, 4.7761, 4.7498, 4.7272, 4.7076, 4.6877, 4.6730, 4.4274 &
        &, 4.5473, 4.5169, 4.5975, 4.5793, 4.5667, 4.5559, 4.3804 &
        &, 4.6920, 4.6731, 4.6142, 4.5600, 4.4801, 4.0149, 3.8856 &
        &, 3.7407, 4.1545, 4.2253, 4.4229, 4.1923, 4.5022, 4.3059 &
        &, 4.1591, 4.7883, 4.9294, 3.3850, 3.4208, 3.7004, 3.8800 &
        &/)
    r0ab(3921:3990)=(/ &
        & 3.9886, 3.9040, 3.6719, 3.6547, 3.4625, 3.3370, 3.8394 &
        &, 4.0335, 4.2373, 4.3023, 4.0306, 4.1408, 3.9297, 3.7857 &
        &, 4.1907, 4.3230, 4.2664, 4.2173, 4.1482, 3.6823, 4.0711 &
        &, 4.0180, 4.0017, 3.9747, 3.9634, 3.9383, 4.1993, 4.3205 &
        &, 4.0821, 4.2547, 4.0659, 3.9359, 4.3952, 4.5176, 4.3888 &
        &, 4.3607, 3.9583, 3.9280, 3.8390, 3.7971, 3.7955, 3.7674 &
        &, 3.7521, 4.1062, 4.3633, 4.2991, 4.2767, 4.4857, 4.3039 &
        &, 4.1762, 4.6197, 4.8654, 4.6633, 4.5878, 4.5640, 4.5422 &
        &, 4.5231, 4.5042, 4.4901, 4.3282, 4.3978, 4.3483, 4.4202 &
        &, 4.4039, 4.3926, 4.3807, 4.2649, 4.6135, 4.5605, 4.5232 &
        &/)
    r0ab(3991:4060)=(/ &
        & 4.4676, 4.3948, 4.0989, 3.9864, 3.8596, 4.0942, 4.2720 &
        &, 4.3270, 4.3022, 4.5410, 4.3576, 4.2235, 4.6545, 4.7447 &
        &, 4.7043, 3.0942, 3.2075, 3.5152, 3.6659, 3.8289, 3.7459 &
        &, 3.5156, 3.5197, 3.3290, 3.2069, 3.6702, 3.8448, 4.0340 &
        &, 3.9509, 3.8585, 3.9894, 3.7787, 3.6365, 4.1425, 4.1618 &
        &, 4.0940, 4.0466, 3.9941, 3.5426, 3.8952, 3.8327, 3.8126 &
        &, 3.7796, 3.7635, 3.7356, 4.0047, 3.9655, 3.9116, 4.1010 &
        &, 3.9102, 3.7800, 4.2964, 4.3330, 4.2622, 4.2254, 3.8195 &
        &, 3.7560, 3.6513, 3.5941, 3.5810, 3.5420, 3.5178, 3.8861 &
        &, 4.1459, 4.1147, 4.0772, 4.3120, 4.1207, 3.9900, 4.4733 &
        &/)
    r0ab(4061:4130)=(/ &
        & 4.6157, 4.4580, 4.4194, 4.3954, 4.3739, 4.3531, 4.3343 &
        &, 4.3196, 4.2140, 4.2339, 4.1738, 4.2458, 4.2278, 4.2158 &
        &, 4.2039, 4.1658, 4.3595, 4.2857, 4.2444, 4.1855, 4.1122 &
        &, 3.7839, 3.6879, 3.5816, 3.8633, 4.1585, 4.1402, 4.1036 &
        &, 4.3694, 4.1735, 4.0368, 4.5095, 4.5538, 4.5240, 4.4252 &
        &, 3.0187, 3.1918, 3.5127, 3.6875, 3.7404, 3.6943, 3.4702 &
        &, 3.4888, 3.2914, 3.1643, 3.6669, 3.8724, 3.9940, 4.0816 &
        &, 3.8054, 3.9661, 3.7492, 3.6024, 4.0428, 4.1951, 4.1466 &
        &, 4.0515, 4.0075, 3.5020, 3.9158, 3.8546, 3.8342, 3.8008 &
        &, 3.7845, 3.7549, 3.9602, 3.8872, 3.8564, 4.0793, 3.8835 &
        &/)
    r0ab(4131:4200)=(/ &
        & 3.7495, 4.2213, 4.3704, 4.3300, 4.2121, 3.7643, 3.7130 &
        &, 3.6144, 3.5599, 3.5474, 3.5093, 3.4853, 3.9075, 4.1115 &
        &, 4.0473, 4.0318, 4.2999, 4.1050, 3.9710, 4.4320, 4.6706 &
        &, 4.5273, 4.4581, 4.4332, 4.4064, 4.3873, 4.3684, 4.3537 &
        &, 4.2728, 4.2549, 4.2032, 4.2794, 4.2613, 4.2491, 4.2375 &
        &, 4.2322, 4.3665, 4.3061, 4.2714, 4.2155, 4.1416, 3.7660 &
        &, 3.6628, 3.5476, 3.8790, 4.1233, 4.0738, 4.0575, 4.3575 &
        &, 4.1586, 4.0183, 4.4593, 4.5927, 4.4865, 4.3813, 4.4594 &
        &, 2.9875, 3.1674, 3.4971, 3.6715, 3.7114, 3.6692, 3.4446 &
        &, 3.4676, 3.2685, 3.1405, 3.6546, 3.8579, 3.9637, 4.0581 &
        &/)
    r0ab(4201:4270)=(/ &
        & 3.7796, 3.9463, 3.7275, 3.5792, 4.0295, 4.1824, 4.1247 &
        &, 4.0357, 3.9926, 3.4827, 3.9007, 3.8392, 3.8191, 3.7851 &
        &, 3.7687, 3.7387, 3.9290, 3.8606, 3.8306, 4.0601, 3.8625 &
        &, 3.7269, 4.2062, 4.3566, 4.3022, 4.1929, 3.7401, 3.6888 &
        &, 3.5900, 3.5350, 3.5226, 3.4838, 3.4594, 3.8888, 4.0813 &
        &, 4.0209, 4.0059, 4.2810, 4.0843, 3.9486, 4.4162, 4.6542 &
        &, 4.5005, 4.4444, 4.4196, 4.3933, 4.3741, 4.3552, 4.3406 &
        &, 4.2484, 4.2413, 4.1907, 4.2656, 4.2474, 4.2352, 4.2236 &
        &, 4.2068, 4.3410, 4.2817, 4.2479, 4.1921, 4.1182, 3.7346 &
        &, 3.6314, 3.5168, 3.8582, 4.0927, 4.0469, 4.0313, 4.3391 &
        &/)
    r0ab(4271:4340)=(/ &
        & 4.1381, 3.9962, 4.4429, 4.5787, 4.4731, 4.3588, 4.4270 &
        &, 4.3957, 2.9659, 3.1442, 3.4795, 3.6503, 3.6814, 3.6476 &
        &, 3.4222, 3.4491, 3.2494, 3.1209, 3.6324, 3.8375, 3.9397 &
        &, 3.8311, 3.7581, 3.9274, 3.7085, 3.5598, 4.0080, 4.1641 &
        &, 4.1057, 4.0158, 3.9726, 3.4667, 3.8802, 3.8188, 3.7989 &
        &, 3.7644, 3.7474, 3.7173, 3.9049, 3.8424, 3.8095, 4.0412 &
        &, 3.8436, 3.7077, 4.1837, 4.3366, 4.2816, 4.1686, 3.7293 &
        &, 3.6709, 3.5700, 3.5153, 3.5039, 3.4684, 3.4437, 3.8663 &
        &, 4.0575, 4.0020, 3.9842, 4.2612, 4.0643, 3.9285, 4.3928 &
        &, 4.6308, 4.4799, 4.4244, 4.3996, 4.3737, 4.3547, 4.3358 &
        &/)
    r0ab(4341:4410)=(/ &
        & 4.3212, 4.2275, 4.2216, 4.1676, 4.2465, 4.2283, 4.2161 &
        &, 4.2045, 4.1841, 4.3135, 4.2562, 4.2226, 4.1667, 4.0932 &
        &, 3.7134, 3.6109, 3.4962, 3.8352, 4.0688, 4.0281, 4.0099 &
        &, 4.3199, 4.1188, 3.9768, 4.4192, 4.5577, 4.4516, 4.3365 &
        &, 4.4058, 4.3745, 4.3539, 2.8763, 3.1294, 3.5598, 3.7465 &
        &, 3.5659, 3.5816, 3.3599, 3.4024, 3.1877, 3.0484, 3.7009 &
        &, 3.9451, 3.8465, 3.9873, 3.7079, 3.9083, 3.6756, 3.5150 &
        &, 4.0829, 4.2780, 4.1511, 4.1260, 4.0571, 3.4865, 3.9744 &
        &, 3.9150, 3.8930, 3.8578, 3.8402, 3.8073, 3.7977, 4.0036 &
        &, 3.7604, 4.0288, 3.8210, 3.6757, 4.2646, 4.4558, 4.2862 &
        &/)
    r0ab(4411:4465)=(/ &
        & 4.2122, 3.7088, 3.6729, 3.5800, 3.5276, 3.5165, 3.4783 &
        &, 3.4539, 3.9553, 3.9818, 4.2040, 3.9604, 4.2718, 4.0689 &
        &, 3.9253, 4.4869, 4.7792, 4.4918, 4.5342, 4.5090, 4.4868 &
        &, 4.4680, 4.4486, 4.4341, 4.2023, 4.3122, 4.2710, 4.3587 &
        &, 4.3407, 4.3281, 4.3174, 4.1499, 4.3940, 4.3895, 4.3260 &
        &, 4.2725, 4.1961, 3.7361, 3.6193, 3.4916, 3.9115, 3.9914 &
        &, 3.9809, 3.9866, 4.3329, 4.1276, 3.9782, 4.5097, 4.6769 &
        &, 4.5158, 4.3291, 4.3609, 4.3462, 4.3265, 4.4341 &
        &/)

    k=0
    do i=1,max_elem
      do j=1,i
        k=k+1
        r(i,j)=r0ab(k)/autoang
        r(j,i)=r0ab(k)/autoang
      end do
    end do

  end subroutine setr0ab


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine stoprun(s)
    character*(*) s
    write(*,*)'program stopped due to: ',s
    open(99, file='dscf_problem', action='write', status='replace')
    close(99)
    !call system('touch dscf_problem')
    stop 'must stop!'
  end subroutine stoprun

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! Returns the number of a given element string (h-pu, 1-94)
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  elemental subroutine ELEM(KEY1, NAT)
    CHARACTER(*), intent(in) :: KEY1
    INTEGER, intent(out) :: NAT

    character(2), parameter :: ELEMNT(94) = [&
        & 'h ','he', &
        & 'li','be','b ','c ','n ','o ','f ','ne', &
        & 'na','mg','al','si','p ','s ','cl','ar', &
        & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu', &
        & 'zn','ga','ge','as','se','br','kr', &
        & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag', &
        & 'cd','in','sn','sb','te','i ','xe', &
        & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy', &
        & 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt', &
        & 'au','hg','tl','pb','bi','po','at','rn', &
        & 'fr','ra','ac','th','pa','u ','np','pu' ]

    CHARACTER(2) :: E
    integer :: k, j, n, i

    nat=0
    e=' '
    k=1
    do J=1,len(key1)
      if (k.gt.2)exit
      N=ICHAR(key1(J:J))
      if (n.ge.ichar('A') .and. n.le.ichar('Z') )then
        e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
        k=k+1
      end if
      if (n.ge.ichar('a') .and. n.le.ichar('z') )then
        e(k:k)=key1(j:j)
        k=k+1
      end if
    end do

    do I=1,94
      if (e.eq.elemnt(i))then
        NAT=I
        RETURN
      end if
    end do


  end subroutine ELEM



  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! B E G I N O F P B C P A R T
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! compute coordination numbers by adding an inverse damping function
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine pbcncoord(natoms,rcov,iz,xyz,cn,lat,rep_cn,crit_cn)
    integer,intent(in) :: natoms,iz(*)
    real(wp),intent(in) :: rcov(94)
    real(wp), intent(in):: xyz(3,*),lat(3,3)
    real(wp), intent(in) :: crit_cn
    real(wp), intent(out):: cn(*)

    integer i,max_elem,rep_cn(3)

    integer iat,taux,tauy,tauz
    real(wp) dx,dy,dz,r,damp,xn,rr,rco,tau(3)

    do i=1,natoms
      xn=0.0d0
      do iat=1,natoms
        do taux=-rep_cn(1),rep_cn(1)
          do tauy=-rep_cn(2),rep_cn(2)
            do tauz=-rep_cn(3),rep_cn(3)
              if (iat.eq.i .and. taux.eq.0 .and. tauy.eq.0 .and.&
                  & tauz.eq.0) cycle
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
              dx=xyz(1,iat)-xyz(1,i)+tau(1)
              dy=xyz(2,iat)-xyz(2,i)+tau(2)
              dz=xyz(3,iat)-xyz(3,i)+tau(3)
              r=(dx*dx+dy*dy+dz*dz)
              if (r.gt.crit_cn) cycle
              r=sqrt(r)
              ! covalent distance in Bohr
              rco=rcov(iz(i))+rcov(iz(iat))
              rr=rco/r
              ! counting function exponential has a better long-range behavior than MH
              damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
              xn=xn+damp
              ! print '("cn(",I2,I2,"): ",E14.8)',i,iat,damp

            end do
          end do
        end do
      end do
      cn(i)=xn
    end do

  end subroutine pbcncoord



  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! compute energy
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
      & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,&
      & e6,e8,e10,e12,e63,lat,rthr,rep_vdw,cn_thr,rep_cn)

    USE mp_images,    ONLY : me_image , nproc_image, intra_image_comm
    USE mp,           ONLY : mp_sum
    integer :: mykey, na_s, na_e

    integer, intent(in) :: max_elem,maxc
    real(wp), intent(in):: r2r4(max_elem),rcov(max_elem)
    real(wp), intent(in)::  rs6,rs8,rs10,alp6,alp8,alp10
    real(wp), intent(in):: rthr,cn_thr
    integer, intent(in):: rep_vdw(3),rep_cn(3)
    integer, intent(in):: n,iz(*),version,mxc(max_elem)
    ! integer rep_v(3)=rep_vdw!,rep_cn(3)
    real(wp), intent(in):: xyz(3,*),r0ab(max_elem,max_elem),lat(3,3)
    ! real(wp) rs6,rs8,rs10,alp6,alp8,alp10,rcov(max_elem)
    real(wp), intent(in):: c6ab(max_elem,max_elem,maxc,maxc,3)
    real(wp), intent(out) :: e6, e8, e10, e12, e63
    logical, intent(in):: noabc

    integer iat,jat,kat
    real(wp) :: crit_cn
    real(wp) r,r2,r6,r8,tmp,dx,dy,dz,c6,c8,c10,ang,rav,R0
    real(wp) damp6,damp8,damp10,rr,thr,c9,r42,c12,r10,c14
    real(wp) cn(n),rxyz(3),dxyz(3)
    real(wp) cc6ab(n*n),dmp(n*n),d2(3),t1,t2,t3,tau(3)
    integer ij,ik,jk
    integer taux,tauy,tauz,counter
    real(wp) a1,a2
    real(wp) bj_dmp6,bj_dmp8
    real(wp) tmp1,tmp2

    e6 =0
    e8 =0
    e10=0
    e12=0
    e63=0
    tau=(/0.0,0.0,0.0/)
    counter=0
    crit_cn=cn_thr
    cc6ab(:) = 0.0_wp
    ! Becke-Johnson parameters
    a1=rs6
    a2=rs8


    CALL block_distribute( n, me_image, nproc_image, na_s, na_e, mykey )
    IF ( mykey == 0 ) THEN


    ! DFT-D2
    if (version.eq.2)then


      do iat=na_s,min(na_e,n-1)
        do jat=iat+1,n
          c6=c6ab(iz(jat),iz(iat),1,1,1)
          do taux=-rep_vdw(1),rep_vdw(1)
            do tauy=-rep_vdw(2),rep_vdw(2)
              do tauz=-rep_vdw(3),rep_vdw(3)
                tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
                dx=xyz(1,iat)-xyz(1,jat)+tau(1)
                dy=xyz(2,iat)-xyz(2,jat)+tau(2)
                dz=xyz(3,iat)-xyz(3,jat)+tau(3)
                r2=dx*dx+dy*dy+dz*dz
                if (r2.gt.rthr) cycle
                r=sqrt(r2)
                damp6=1./(1.+exp(-alp6*(r/(rs6*r0ab(iz(jat),iz(iat)))-1.)))
                r6=r2**3
                e6 =e6+c6*damp6/r6
              end do
            end do
          end do
        end do
      end do

      do iat=na_s,na_e
        jat=iat
        c6=c6ab(iz(jat),iz(iat),1,1,1)
        do taux=-rep_vdw(1),rep_vdw(1)
          do tauy=-rep_vdw(2),rep_vdw(2)
            do tauz=-rep_vdw(3),rep_vdw(3)
              if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
              dx=tau(1)
              dy=tau(2)
              dz=tau(3)
              r2=dx*dx+dy*dy+dz*dz
              if (r2.gt.rthr) cycle
              r=sqrt(r2)
              damp6=1./(1.+exp(-alp6*(r/(rs6*r0ab(iz(jat),iz(iat)))-1.)))
              r6=r2**3
              e6 =e6+c6*damp6/r6*0.50d0
            end do
          end do
        end do
      end do



    else if ((version.eq.3).or.(version.eq.5)) then
      ! DFT-D3(zero-damping)

      call pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,crit_cn)

      do iat=na_s,min(na_e,n-1)
        do jat=iat+1,n
          ! get C6
          call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),&
              & cn(iat),cn(jat),c6)

          if (.not.noabc)then
            ij=lin(jat,iat)
            ! store C6 for C9, calc as sqrt
            cc6ab(ij)=sqrt(c6)
          end if
          do taux=-rep_vdw(1),rep_vdw(1)
            do tauy=-rep_vdw(2),rep_vdw(2)
              do tauz=-rep_vdw(3),rep_vdw(3)
                tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

                dx=xyz(1,iat)-xyz(1,jat)+tau(1)
                dy=xyz(2,iat)-xyz(2,jat)+tau(2)
                dz=xyz(3,iat)-xyz(3,jat)+tau(3)
                r2=dx*dx+dy*dy+dz*dz
                ! cutoff

                if (r2.gt.rthr) cycle
                r =sqrt(r2)
                R0=r0ab(iz(jat),iz(iat))
                rr=R0/r
                ! damping
                if(version.eq.3)then
                  ! DFT-D3 zero-damp
                  tmp=rs6*rr
                  damp6 =1.d0/( 1.d0+6.d0*tmp**alp6 )
                  tmp=rs8*rr
                  damp8 =1.d0/( 1.d0+6.d0*tmp**alp8 )
                else
                  ! DFT-D3M zero-damp
                  tmp=(r/(rs6*R0))+rs8*R0
                  damp6 =1.d0/( 1.d0+6.d0*tmp**(-alp6) )
                  tmp=(r/R0)+rs8*R0
                  damp8 =1.d0/( 1.d0+6.d0*tmp**(-alp8) )
                endif

                r6=r2**3
                e6 =e6+damp6/r6* c6
                ! write(*,*)'e6: ',c6*damp6/r6*autokcal

                ! stored in main as sqrt
                c8 =3.0d0*r2r4(iz(iat))*r2r4(iz(jat))*c6
                r8 =r6*r2

                e8 =e8+c8*damp8/r8

              end do
            end do
          end do
        end do
      end do

      do iat=na_s,na_e
        jat=iat
        ! get C6
        call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),&
            & cn(iat),cn(jat),c6)

        if (.not.noabc)then
          ij=lin(jat,iat)
          ! store C6 for C9, calc as sqrt
          cc6ab(ij)=sqrt(c6)
        end if
        do taux=-rep_vdw(1),rep_vdw(1)
          do tauy=-rep_vdw(2),rep_vdw(2)
            do tauz=-rep_vdw(3),rep_vdw(3)
              if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

              dx=tau(1)
              dy=tau(2)
              dz=tau(3)
              r2=dx*dx+dy*dy+dz*dz
              ! cutoff
              if (r2.gt.rthr) cycle
              r =sqrt(r2)
              R0=r0ab(iz(jat),iz(iat))
              rr=R0/r

              ! damping 
              if(version.eq.3)then
                ! DFT-D3 zero-damp
                tmp=rs6*rr
                damp6 =1.d0/( 1.d0+6.d0*tmp**alp6 )
                tmp=rs8*rr
                damp8 =1.d0/( 1.d0+6.d0*tmp**alp8 )
              else
                ! DFT-D3M zero-damp
                tmp=(r/(rs6*R0))+rs8*R0
                damp6 =1.d0/( 1.d0+6.d0*tmp**(-alp6) )
                tmp=(r/R0)+rs8*R0
                damp8 =1.d0/( 1.d0+6.d0*tmp**(-alp8) )
              endif


              r6=r2**3

              e6 =e6+damp6/r6*0.50d0 *C6

              ! stored in main as sqrt
              c8 =3.0d0*r2r4(iz(iat))*r2r4(iz(jat)) *C6
              r8 =r6*r2

              e8 =e8+c8*damp8/r8*0.50d0
              counter=counter+1

            end do
          end do
        end do
      end do
      ! write(*,*)'counter(edisp): ',counter
    else if((version.eq.4).or.(version.eq.6)) then


      ! DFT-D3(BJ-damping)
      call pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,crit_cn)

      do iat=na_s,na_e
        do jat=iat+1,n
          ! get C6
          call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),&
              & cn(iat),cn(jat),c6)

          rxyz=xyz(:,iat)-xyz(:,jat)
          r42=r2r4(iz(iat))*r2r4(iz(jat))
          bj_dmp6=(a1*dsqrt(3.0d0*r42)+a2)**6
          bj_dmp8=(a1*dsqrt(3.0d0*r42)+a2)**8

          if (.not.noabc)then
            ij=lin(jat,iat)
            ! store C6 for C9, calc as sqrt
            cc6ab(ij)=sqrt(c6)
          end if
          do taux=-rep_vdw(1),rep_vdw(1)
            do tauy=-rep_vdw(2),rep_vdw(2)
              do tauz=-rep_vdw(3),rep_vdw(3)
                tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

                dxyz=rxyz+tau

                r2=sum(dxyz*dxyz)
                ! cutoff
                if (r2.gt.rthr) cycle
                r =sqrt(r2)
                rr=r0ab(iz(jat),iz(iat))/r


                r6=r2**3

                e6 =e6+c6/(r6+bj_dmp6)

                ! stored in main as sqrt
                c8 =3.0d0*c6*r42
                r8 =r6*r2

                e8 =e8+c8/(r8+bj_dmp8)

                counter=counter+1

              end do
            end do
          end do
        end do

        ! Now the self interaction
        jat=iat
        ! get C6
        call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),&
            & cn(iat),cn(jat),c6)
        r42=r2r4(iz(iat))*r2r4(iz(iat))
        bj_dmp6=(a1*dsqrt(3.0d0*r42)+a2)**6
        bj_dmp8=(a1*dsqrt(3.0d0*r42)+a2)**8

        if (.not.noabc)then
          ij=lin(jat,iat)
          ! store C6 for C9, calc as sqrt
          cc6ab(ij)=dsqrt(c6)
        end if

        do taux=-rep_vdw(1),rep_vdw(1)
          do tauy=-rep_vdw(2),rep_vdw(2)
            do tauz=-rep_vdw(3),rep_vdw(3)
              if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

              r2=sum(tau*tau)
              ! cutoff
              if (r2.gt.rthr) cycle
              r =sqrt(r2)
              rr=r0ab(iz(jat),iz(iat))/r


              r6=r2**3

              e6 =e6+c6/(r6+bj_dmp6)*0.50d0

              ! stored in main as sqrt
              c8 =3.0d0*c6*r42
              r8 =r6*r2

              e8 =e8+c8/(r8+bj_dmp8)*0.50d0
              counter=counter+1

            end do
          end do
        end do
      end do


    end if


    ENDIF
    CALL mp_sum ( e6 , intra_image_comm )
    CALL mp_sum ( e8 , intra_image_comm )

    if (.not.noabc) then 

    ! compute non-additive third-order energy using averaged C6
    CALL mp_sum ( cc6ab , intra_image_comm )
    call pbcthreebody(max_elem,xyz,lat,n,iz,rep_cn,crit_cn,&
        & cc6ab,r0ab,e63)

    end if 

  end subroutine pbcedisp


  subroutine pbcthreebody(max_elem,xyz,lat,n,iz,repv,cnthr,cc6ab,&
      & r0ab,eabc)

    USE mp_images,    ONLY : me_image , nproc_image, intra_image_comm
    USE mp,           ONLY : mp_sum
    integer :: mykey, na_s, na_smax, na_e

    integer max_elem
    INTEGER :: n,i,j,k,jtaux,jtauy,jtauz,iat,jat,kat
    INTEGER :: ktaux,ktauy,ktauz,counter,ij,ik,jk,idum
    REAL(WP) :: dx,dy,dz,rij2,rik2,rjk2,c9,rr0ij,rr0ik
    REAL(WP) :: rr0jk,geomean,fdamp,rik,rjk,rij
    REAL(WP) :: r0ij,r0ik,r0jk
    REAL(WP),INTENT(OUT)::eabc
    REAL(WP) :: tmp,tmp1,tmp2,tmp3,tmp4,ang

    REAL(WP) ,DIMENSION(3,3),INTENT(IN)::lat
    REAL(WP) ,DIMENSION(3,*),INTENT(IN) :: xyz
    INTEGER,DIMENSION(*),INTENT(IN)::iz
    REAL(WP),DIMENSION(3):: jtau,ktau,jxyz,kxyz,ijvec,ikvec,jkvec,dumvec
    INTEGER,DIMENSION(3):: repv
    REAL(WP),INTENT(IN) ::cnthr
    REAL(WP),DIMENSION(n*n),INTENT(IN)::cc6ab
    REAL(WP),DIMENSION(max_elem,max_elem),INTENT(IN):: r0ab
    REAL(WP),PARAMETER::sr9=0.75d0
    REAL(WP),PARAMETER::alp9=-16.0d0
    REAL(WP) :: abcthr
    INTEGER :: repmin1, repmin2, repmin3, repmax1, repmax2, repmax3    
    INTEGER :: repv1, repv2, repv3 
    REAL(WP)  :: ijvec1, ijvec2, ijvec3
    REAL(WP)  :: ikvec1, ikvec2, ikvec3
    REAL(WP)  :: jkvec1, jkvec2, jkvec3
    REAL(WP)  :: jtau1, jtau2, jtau3
    REAL(WP)  :: ktau1, ktau2, ktau3
    REAL(WP)  :: dumvec11, dumvec12, dumvec13  
    REAL(WP)  :: dumvec21, dumvec22, dumvec23  
    ! REAL(WP) :: time1,time2

    counter=0
    eabc=0.0d0
    abcthr=cnthr
    ! abcthr=1.0d99
    ! write(*,*)'thr:',(abcthr)

    ! call cpu_time(time1)

    repv1 = repv(1)
    repv2 = repv(2)
    repv3 = repv(3)

    CALL block_distribute( n, me_image, nproc_image, na_s, na_e, mykey )
    IF ( mykey == 0 ) THEN
    na_smax = max(3,na_s)

!$acc data copyin(xyz(1:3,1:n),iz(1:n),cc6ab(1:n*n),lat(1:3,1:3),r0ab(1:max_elem,1:max_elem)) 
!$acc kernels  vector_length(32)
!$acc loop collapse(2) gang private(ijvec1,ijvec2,ijvec3, ikvec1,ikvec2,ikvec3, jkvec1,jkvec2,jkvec3, c9, r0ij,r0ik,r0jk, &
!$acc&                              repmin1,repmin2,repmin3, repmax1,repmax2,repmax3, jtau1,jtau2,jtau3, &
!$acc&                              dumvec11,dumvec12,dumvec13,rij2,rr0ij)  reduction(+:eabc)
    do iat=na_smax,na_e
     do jat = 2, n
        if(jat.ge.iat) cycle
        ijvec1=xyz(1,jat)-xyz(1,iat)
        ijvec2=xyz(2,jat)-xyz(2,iat)
        ijvec3=xyz(3,jat)-xyz(3,iat)
        ij=lin(iat,jat)
        r0ij=r0ab(iz(iat),iz(jat))
!$acc loop seq
        do kat=1,jat-1
          ik=lin(iat,kat)
          jk=lin(jat,kat)
          ikvec1=xyz(1,kat)-xyz(1,iat)
          ikvec2=xyz(2,kat)-xyz(2,iat)
          ikvec3=xyz(3,kat)-xyz(3,iat)
          jkvec1=xyz(1,kat)-xyz(1,jat)
          jkvec2=xyz(2,kat)-xyz(2,jat)
          jkvec3=xyz(3,kat)-xyz(3,jat)
          c9=-1.0d0*(cc6ab(ij)*cc6ab(ik)*cc6ab(jk))

          r0ik=r0ab(iz(iat),iz(kat))
          r0jk=r0ab(iz(jat),iz(kat))

          do jtaux=-repv1,repv1
            do jtauy=-repv2,repv2
              do jtauz=-repv3,repv3
                repmin1=max(-repv1,jtaux-repv1)
                repmax1=min(repv1,jtaux+repv1)
                repmin2=max(-repv2,jtauy-repv2)
                repmax2=min(repv2,jtauy+repv2)
                repmin3=max(-repv3,jtauz-repv3)
                repmax3=min(repv3,jtauz+repv3)
                jtau1=jtaux*lat(1,1)+jtauy*lat(1,2)+jtauz*lat(1,3)
                jtau2=jtaux*lat(2,1)+jtauy*lat(2,2)+jtauz*lat(2,3)
                jtau3=jtaux*lat(3,1)+jtauy*lat(3,2)+jtauz*lat(3,3)
                dumvec11=(ijvec1+jtau1)*(ijvec1+jtau1)
                dumvec12=(ijvec2+jtau2)*(ijvec2+jtau2)
                dumvec13=(ijvec3+jtau3)*(ijvec3+jtau3)
                rij2=dumvec11+dumvec12+dumvec13
                if (rij2.gt.abcthr)cycle

                rr0ij=DSQRT(rij2)/r0ij

!$acc loop vector collapse(3) private(ktau1,ktau2,ktau3,dumvec21,dumvec22,dumvec23,rik2,rr0ik,rjk2,rr0jk,geomean,fdamp,tmp1,tmp2,tmp3,tmp4,ang) reduction(+:eabc)
                do ktaux=repmin1,repmax1
                  do ktauy=repmin2,repmax2
                    do ktauz=repmin3,repmax3
                      ktau1=ktaux*lat(1,1)+ktauy*lat(1,2)+ktauz*lat(1,3)
                      ktau2=ktaux*lat(2,1)+ktauy*lat(2,2)+ktauz*lat(2,3)
                      ktau3=ktaux*lat(3,1)+ktauy*lat(3,2)+ktauz*lat(3,3)
                      dumvec21=(ikvec1+ktau1)*(ikvec1+ktau1)
                      dumvec22=(ikvec2+ktau2)*(ikvec2+ktau2)
                      dumvec23=(ikvec3+ktau3)*(ikvec3+ktau3)
                      rik2=dumvec21+dumvec22+dumvec23
                      if (rik2.gt.abcthr)cycle
                      rr0ik=DSQRT(rik2)/r0ik
                      dumvec21=jkvec1+ktau1-jtau1
                      dumvec22=jkvec2+ktau2-jtau2
                      dumvec23=jkvec3+ktau3-jtau3
                      rjk2=dumvec21*dumvec21+dumvec22*dumvec22+dumvec23*dumvec23
                      if (rjk2.gt.abcthr)cycle
                      rr0jk=DSQRT(rjk2)/r0jk


                      geomean=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
                      ! write(*,*)'geomean:',geomean
                      fdamp=1./(1.+6.*(sr9*geomean)**alp9)
                      tmp1 = (rij2+rjk2-rik2)
                      tmp2 = (rij2+rik2-rjk2)
                      tmp3 = (rik2+rjk2-rij2)
                      tmp4=rij2*rjk2*rik2
                      ang=(0.375d0*tmp1*tmp2*tmp3/tmp4+1.0d0)/tmp4**1.5d0

                      eabc=eabc+ang*c9*fdamp

                    end do
                  end do
                end do

              end do
            end do
          end do

        end do
      end do
    end do
!$acc end kernels  
!$acc end data 

    do iat=max(2,na_s),na_e
      jat=iat
      ij=lin(iat,jat)
      ijvec=0.0d0
      r0ij=r0ab(iz(iat),iz(jat))
      do kat=1,iat-1
        jk=lin(jat,kat)
        ik=jk
        ikvec=xyz(:,kat)-xyz(:,iat)
        jkvec=ikvec
        c9=-(cc6ab(ij)*cc6ab(ik)*cc6ab(jk))

        r0ik=r0ab(iz(iat),iz(kat))
        r0jk=r0ab(iz(jat),iz(kat))
        do jtaux=-repv1,repv1
          do jtauy=-repv2,repv2
            do jtauz=-repv3,repv3
              repmin1=max(-repv1,jtaux-repv1)
              repmax1=min(repv1,jtaux+repv1)
              repmin2=max(-repv2,jtauy-repv2)
              repmax2=min(repv2,jtauy+repv2)
              repmin3=max(-repv3,jtauz-repv3)
              repmax3=min(repv3,jtauz+repv3)
              if (jtaux.eq.0 .and. jtauy.eq.0 .and. jtauz.eq.0) cycle
              jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
              dumvec=ijvec+jtau
              dumvec=dumvec*dumvec
              rij2=SUM(dumvec)
              if (rij2.gt.abcthr)cycle

              rr0ij=DSQRT(rij2)/r0ij

              do ktaux=repmin1,repmax1
                do ktauy=repmin2,repmax2
                  do ktauz=repmin3,repmax3
                    ! every result * 0.5
                    ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
                    dumvec=ikvec+ktau
                    dumvec=dumvec*dumvec
                    rik2=SUM(dumvec)
                    if (rik2.gt.abcthr)cycle
                    rr0ik=DSQRT(rik2)/r0ik

                    dumvec=jkvec+ktau-jtau
                    dumvec=dumvec*dumvec
                    rjk2=SUM(dumvec)
                    if (rjk2.gt.abcthr)cycle
                    rr0jk=DSQRT(rjk2)/r0jk


                    geomean=(rr0ij*rr0ik*rr0jk)**(1./3.)
                    fdamp=1./(1.+6.*(sr9*geomean)**alp9)
                    tmp1 = (rij2+rjk2-rik2)
                    tmp2 = (rij2+rik2-rjk2)
                    tmp3 = (rik2+rjk2-rij2)
                    tmp4=rij2*rjk2*rik2
                    ang=(0.375d0*tmp1*tmp2*tmp3/tmp4+1.0d0)/tmp4**1.5d0

                    eabc=eabc+c9*fdamp*ang/2.0
                  end do
                end do
              end do

            end do
          end do
        end do
      end do
    end do

    do iat=max(2,na_s),na_e
      do jat=1,iat-1
        kat=jat
        ij=lin(iat,jat)
        jk=lin(jat,kat)
        ik=ij
        ikvec=xyz(:,kat)-xyz(:,iat)
        ijvec=ikvec
        jkvec=0.0d0
        c9=-(cc6ab(ij)*cc6ab(ik)*cc6ab(jk))

        r0ij=r0ab(iz(iat),iz(jat))
        r0ik=r0ij
        r0jk=r0ab(iz(jat),iz(kat))

        do jtaux=-repv1,repv1
          do jtauy=-repv2,repv2
            do jtauz=-repv3,repv3
              repmin1=max(-repv1,jtaux-repv1)
              repmax1=min(repv1,jtaux+repv1)
              repmin2=max(-repv2,jtauy-repv2)
              repmax2=min(repv2,jtauy+repv2)
              repmin3=max(-repv3,jtauz-repv3)
              repmax3=min(repv3,jtauz+repv3)
              jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
              dumvec=ijvec+jtau
              dumvec=dumvec*dumvec
              rij2=SUM(dumvec)
              if (rij2.gt.abcthr)cycle

              rr0ij=DSQRT(rij2)/r0ij

              do ktaux=repmin1,repmax1
                do ktauy=repmin2,repmax2
                  do ktauz=repmin3,repmax3
                    ! every result * 0.5
                    if (jtaux.eq.ktaux .and. jtauy.eq.ktauy&
                        & .and. jtauz.eq.ktauz) cycle
                    ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
                    dumvec=ikvec+ktau
                    dumvec=dumvec*dumvec
                    rik2=SUM(dumvec)
                    if (rik2.gt.abcthr)cycle
                    rr0ik=DSQRT(rik2)/r0ik

                    dumvec=jkvec+ktau-jtau
                    dumvec=dumvec*dumvec
                    rjk2=SUM(dumvec)
                    if (rjk2.gt.abcthr)cycle
                    rr0jk=DSQRT(rjk2)/r0jk


                    geomean=(rr0ij*rr0ik*rr0jk)**(1./3.)
                    fdamp=1./(1.+6.*(sr9*geomean)**alp9)
                    tmp1 = (rij2+rjk2-rik2)
                    tmp2 = (rij2+rik2-rjk2)
                    tmp3 = (rik2+rjk2-rij2)
                    tmp4=rij2*rjk2*rik2
                    ang=(0.375d0*tmp1*tmp2*tmp3/tmp4+1.0d0)/tmp4**1.5d0

                    eabc=eabc+c9*fdamp*ang/2.0
                  end do
                end do
              end do

            end do
          end do
        end do
      end do
    end do


    ! And finally the self interaction iat=jat=kat all

    idum=0
    do iat=na_s,na_e
      jat=iat
      kat=iat
      ijvec=0.0d0
      ij=lin(iat,iat)
      ik=ij
      jk=ij
      ikvec=ijvec
      jkvec=ikvec
      c9=-(cc6ab(ij)*cc6ab(ik)*cc6ab(jk))

      r0ij=r0ab(iz(iat),iz(iat))
      r0ik=r0ij
      r0jk=r0ij
      do jtaux=-repv1,repv1
        do jtauy=-repv2,repv2
          do jtauz=-repv3,repv3
            repmin1=max(-repv1,jtaux-repv1)
            repmax1=min(repv1,jtaux+repv1)
            repmin2=max(-repv2,jtauy-repv2)
            repmax2=min(repv2,jtauy+repv2)
            repmin3=max(-repv3,jtauz-repv3)
            repmax3=min(repv3,jtauz+repv3)
            if (jtaux.eq.0 .and. jtauy.eq.0 .and. jtauz.eq.0) cycle
            jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
            dumvec=jtau
            dumvec=dumvec*dumvec
            rij2=SUM(dumvec)
            if (rij2.gt.abcthr)cycle
            rr0ij=DSQRT(rij2)/r0ij

            do ktaux=repmin1,repmax1
              do ktauy=repmin2,repmax2
                do ktauz=repmin3,repmax3
                  if ((ktaux.eq.0) .and.( ktauy.eq.0) .and.( ktauz.eq.0))cycle
                  if ((ktaux.eq.jtaux) .and. (ktauy.eq.jtauy)&
                      & .and. (ktauz.eq.jtauz)) cycle

                  ! every result * 1/6 becaues every triple is counted twice due to the tw
                  !
                  !plus 1/3 becaues every triple is three times in each unitcell
                  ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
                  dumvec=ktau
                  dumvec=dumvec*dumvec
                  rik2=SUM(dumvec)
                  if (rik2.gt.abcthr)cycle
                  rr0ik=DSQRT(rik2)/r0ik

                  dumvec=jkvec+ktau-jtau
                  dumvec=dumvec*dumvec
                  rjk2=SUM(dumvec)
                  if (rjk2.gt.abcthr)cycle
                  rr0jk=DSQRT(rjk2)/r0jk


                  geomean=(rr0ij*rr0ik*rr0jk)**(1./3.)
                  fdamp=1./(1.+6.*(sr9*geomean)**alp9)
                  tmp1 = (rij2+rjk2-rik2)
                  tmp2 = (rij2+rik2-rjk2)
                  tmp3 = (rik2+rjk2-rij2)
                  tmp4=rij2*rjk2*rik2
                  ang=(0.375d0*tmp1*tmp2*tmp3/tmp4+1.0d0)/tmp4**1.5d0

                  eabc=eabc+c9*fdamp*ang/6.0d0

                end do
              end do
            end do
          end do
        end do
      end do

    end do

   ENDIF
   CALL mp_sum ( eabc , intra_image_comm )

  end subroutine pbcthreebody


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! compute gradient
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine pbcgdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
      & rcov,s6,s18,rs6,rs8,rs10,alp6,alp8,alp10,noabc,num,&
      & version,g,disp,gnorm,stress,lat,rep_v,rep_cn,&
      & crit_vdw,echo,crit_cn)


    USE mp_images,    ONLY : me_image , nproc_image, intra_image_comm
    USE mp,           ONLY : mp_sum
    integer :: mykey, na_s, na_smax, na_e

    integer n,iz(*),max_elem,maxc,version,mxc(max_elem)
    real(wp) xyz(3,*),r0ab(max_elem,max_elem),r2r4(*)
    real(wp) c6ab(max_elem,max_elem,maxc,maxc,3)
    real(wp) g(3,*),s6,s18,rcov(max_elem)
    real(wp) rs6,rs8,rs10,alp10,alp8,alp6
    real(wp) a1,a2
    real(wp) bj_dmp6,bj_dmp8
    logical noabc,num,echo
    ! conversion factors

    integer iat,jat,i,j,kat,my,ny,a,b,idum,tau2
    real(wp) R0,C6,alp,R42,disp,x1,y1,z1,x2,y2,z2,rr,e6abc,fdum
    real(wp) dx,dy,dz,r2,r,r4,r6,r8,r10,r12,t6,t8,t10,damp1
    real(wp) damp6,damp8,damp9,e6,e8,e10,e12,gnorm,tmp1
    real(wp) s10,s8,gC6(3),term,step,dispr,displ,r235,tmp2
    real(wp) cn(n),gx1,gy1,gz1,gx2,gy2,gz2,rthr,testsum
    real(wp), DIMENSION(3,3) :: lat,stress,sigma,virialstress,lat_1
    real(wp), DIMENSION(3,3) :: gC6_stress
    integer, DIMENSION(3) :: rep_v,rep_cn
    real(wp) crit_vdw,crit_cn
    integer taux,tauy,tauz
    real(wp), DIMENSION(3) :: tau,vec12,dxyz,dxyz0
    real(wp) ::outpr(3,3)
    real(wp), DIMENSION(3,3):: outerprod

    real(wp) rij(3),rik(3),rjk(3),r7,r9
    real(wp) rik_dist,rjk_dist
    real(wp) drik,drjk
    real(wp) rcovij
    real(wp) dc6,c6chk
    real(wp) expterm,dcni
    real(wp), allocatable,dimension(:,:,:,:) :: drij
    real(wp), allocatable,dimension(:,:,:,:) :: dcn
    real(wp) dcnn
    real(wp) :: dc6_rest
    real(wp) vec(3),vec2(3),dummy
    real(wp) dc6i(n)
    real(wp) dc6ij(n,n)
    real(wp) dc6_rest_sum(n*(n+1)/2)
    integer linij,linik,linjk
    real(wp) abc(3,n)

    real(wp) eabc
    real(wp) gabc(3,n),glatabc(3,3)
    real(wp) sigma_abc(3,3)
    real(wp) labc,rabc
    real(wp) ,dimension(3) ::ijvec,ikvec,jkvec,jtau,ktau,dumvec
    integer jtaux,jtauy,jtauz,ktaux,ktauy,ktauz,mtaux,mtauy,mtauz
    integer,dimension(3) :: taumin,taumax
    integer mat,linim,linjm,linkm
    real(wp) rij2,rik2,rjk2,c9,c6ij,c6ik,c6jk,geomean,geomean3
    real(wp) rr0ij,rr0jk,rr0ik,dc6iji,dc6ijj
    real(wp) :: sr9=0.75d0
    real(wp), parameter :: alp9=-16.0d0
    real(wp),DIMENSION(n*(n+1)) ::c6save
    real(wp) abcthr,time1,time2,geomean2,r0av,dc9,dfdmp,dang,ang
    integer,dimension(3) ::repv,repmin,repmax
    integer :: rep_v1, rep_v2, rep_v3  
    integer :: rep_cn1, rep_cn2, rep_cn3  
    integer :: repmin1, repmin2, repmin3   
    integer :: repmax1, repmax2, repmax3   
    real(wp) :: dumvec1, dumvec2, dumvec3  
    real(wp) :: ijvec1, ijvec2, ijvec3   
    real(wp) :: ikvec1, ikvec2, ikvec3   
    real(wp) :: jkvec1, jkvec2, jkvec3   
    real(wp) :: jtau1, jtau2, jtau3  
    real(wp) :: ktau1, ktau2, ktau3  

    ! R^2 cut-off
    rthr=crit_vdw
    abcthr=crit_cn
    ! write(*,*)'abcthr:', abcthr**(1./1.)
    sigma=0.0d0
    virialstress=0.0d0
    stress=0.0d0
    gabc=0.0d0
    glatabc=0.0d0

    ! testsum=0.0d0

    if (echo)write(*,*)

    if (num) then
      if (echo) &
          & write(*,*) 'doing numerical gradient O(N^3) ...'

      call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
          & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,&
          & e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)


      disp=-s6*e6-s18*e8-e6abc

      step=2.d-5

      do i=1,n
        do j=1,3
          xyz(j,i)=xyz(j,i)+step
          call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
              & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,&
              & e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)

          dispr=-s6*e6-s18*e8-e6abc
          rabc=e6abc
          xyz(j,i)=xyz(j,i)-2*step
          call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
              & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,&
              & e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)

          displ=-s6*e6-s18*e8-e6abc
          labc=e6abc
          gabc(j,i)=0.5*(rabc-labc)/step
          g(j,i)=0.5*(dispr-displ)/step
          xyz(j,i)=xyz(j,i)+step
        end do
      end do
      if (echo) write(*,*)'Doing numerical stresstensor...'

      call xyz_to_abc(xyz,abc,lat,n)
      step=2.d-5
      if (echo) write(*,*)'step: ',step
      do i=1,3
        do j=1,3
          lat(j,i)=lat(j,i)+step
          call abc_to_xyz(abc,xyz,lat,n)
          call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
              & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,&
              & e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)

          dispr=-s6*e6-s18*e8-e6abc
          labc=e6abc


          lat(j,i)=lat(j,i)-2*step
          call abc_to_xyz(abc,xyz,lat,n)
          call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
              & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,&
              & e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)

          displ=-s6*e6-s18*e8-e6abc
          rabc=e6abc
          stress(j,i)=(dispr-displ)/(step*2.0)
          glatabc(j,i)=(rabc-labc)/(step*2.0d0)

          lat(j,i)=lat(j,i)+step
          call abc_to_xyz(abc,xyz,lat,n)

        end do
      end do

      sigma=0.0d0
      call inv_cell(lat,lat_1)
      do a=1,3
        do b=1,3
          do my=1,3
            sigma(a,b)=sigma(a,b)-stress(a,my)*lat(b,my)
          end do
        end do
      end do

      goto 999

    end if


    if (version.eq.2)then
      if (echo)write(*,*) 'doing analytical gradient D-old O(N^2) ...'
      disp=0
      stress=0.0d0
      do iat=1,n-1
        do jat=iat+1,n
          R0=r0ab(iz(jat),iz(iat))*rs6
          c6=c6ab(iz(jat),iz(iat),1,1,1)*s6
          do taux=-rep_v(1),rep_v(1)
            do tauy=-rep_v(2),rep_v(2)
              do tauz=-rep_v(3),rep_v(3)
                tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
                dxyz=xyz(:,iat)-xyz(:,jat)+tau
                r2 =sum(dxyz*dxyz)
                if (r2.gt.rthr) cycle
                r235=r2**3.5
                r =dsqrt(r2)
                damp6=exp(-alp6*(r/R0-1.0d0))
                damp1=1.+damp6
                tmp1=damp6/(damp1*damp1*r235*R0)
                tmp2=6./(damp1*r*r235)

                term=alp6*tmp1-tmp2
                g(:,iat)=g(:,iat)-term*dxyz*c6
                g(:,jat)=g(:,jat)+term*dxyz*c6
                disp=disp+c6*(1./damp1)/r2**3

                do ny=1,3
                  do my=1,3
                    sigma(my,ny)=sigma(my,ny)+term*dxyz(ny)*dxyz(my)*c6
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      ! and now the self interaction, only for convenient energy in dispersion
      do iat=1,n
        jat=iat
        R0=r0ab(iz(jat),iz(iat))*rs6
        c6=c6ab(iz(jat),iz(iat),1,1,1)*s6
        do taux=-rep_v(1),rep_v(1)
          do tauy=-rep_v(2),rep_v(2)
            do tauz=-rep_v(3),rep_v(3)
              if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

              dxyz=tau
              ! vec12=(/ dx,dy,dz /)
              r2 =sum(dxyz*dxyz)
              if (r2.gt.rthr) cycle
              r235=r2**3.5
              r =dsqrt(r2)
              damp6=exp(-alp6*(r/R0-1.0d0))
              damp1=1.+damp6
              tmp1=damp6/(damp1*damp1*r235*R0)
              tmp2=6./(damp1*r*r235)
              disp=disp+(c6*(1./damp1)/r2**3)*0.50d0
              term=alp6*tmp1-tmp2
              do ny=1,3
                do my=1,3
                  sigma(my,ny)=sigma(my,ny)+term*dxyz(ny)*dxyz(my)*c6*0.5d0
                end do
              end do


            end do
          end do
        end do
      end do

      call inv_cell(lat,lat_1)
      do a=1,3
        do b=1,3
          do my=1,3
            stress(a,b)=stress(a,b)-sigma(a,my)*lat_1(b,my)
          end do
        end do
      end do

      disp=-disp
      ! sigma=virialstress
      goto 999
    end if

    CALL block_distribute( n, me_image, nproc_image, na_s, na_e, mykey )

    if ((version.eq.3).or.(version.eq.5)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! begin ZERO DAMPING GRADIENT
      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (echo)&
          & write(*,*) 'doing analytical gradient O(N^2) ...'
      ! precompute for analytical part
      call pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,crit_cn)


      s8 =s18
      s10=s18
      allocate(drij(-rep_v(3):rep_v(3),-rep_v(2):rep_v(2),&
          & -rep_v(1):rep_v(1),n*(n+1)/2))

      disp=0

      drij=0.0d0
      dc6_rest=0.0d0
      dc6_rest_sum=0.0d0
      c6save=0.0d0
      kat=0
      dc6i=0.0d0
      dc6ij=0.0d0

      IF ( mykey == 0 ) THEN

      do iat=na_s, na_e
        call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),&
            & mxc(iz(iat)),cn(iat),cn(iat),iz(iat),iz(iat),iat,iat,&
            & c6,dc6iji,dc6ijj)

        c6save(lin(iat,iat))=c6
        dc6ij(iat,iat)=dc6iji
        r0=r0ab(iz(iat),iz(iat))
        r42=r2r4(iz(iat))*r2r4(iz(iat))
        rcovij=rcov(iz(iat))+rcov(iz(iat))


        do taux=-rep_v(1),rep_v(1)
          do tauy=-rep_v(2),rep_v(2)
            do tauz=-rep_v(3),rep_v(3)
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)


              !first dE/d(tau) saved in drij(i,i,counter)
              rij=tau
              r2=sum(rij*rij)
              ! if (r2.gt.rthr) cycle

              if (r2.gt.0.1.and.r2.lt.rthr) then


                r=dsqrt(r2)
                r6=r2*r2*r2
                r7=r6*r
                r8=r6*r2
                r9=r8*r

                !
                ! Calculates damping functions:
                if (version.eq.3) then
                  t6 = (r/(rs6*R0))**(-alp6)
                  damp6 =1.d0/( 1.d0+6.d0*t6 )
                  t8 = (r/(rs8*R0))**(-alp8)
                  damp8 =1.d0/( 1.d0+6.d0*t8 )

                  drij(tauz,tauy,taux,lin(iat,iat))=drij(tauz,tauy,taux,lin(iat,&
                      & iat))&
                      & +(-s6*(6.0/(r7)*C6*damp6)&
                      & -s8*(24.0/(r9)*C6*r42*damp8))*0.5d0


                  drij(tauz,tauy,taux,lin(iat,iat))=drij(tauz,tauy,taux,lin(iat,&
                      & iat))&
                      & +(s6*C6/r7*6.d0*alp6*t6*damp6*damp6&
                      & +s8*C6*r42/r9*18.d0*alp8*t8*damp8*damp8)*0.5d0
                else !version.eq.5
                  t6 = (r/(rs6*R0)+R0*rs8)**(-alp6)
                  damp6 =1.d0/( 1.d0+6.d0*t6 )
                  t8 = (r/(R0)+R0*rs8)**(-alp8)
                  damp8 =1.d0/( 1.d0+6.d0*t8 )
  
                  tmp1=s6*6.d0*damp6*C6/r7
                  tmp2=s8*6.d0*C6*r42*damp8/r9
                  drij(tauz,tauy,taux,lin(iat,iat))=drij(tauz,tauy,taux,lin(iat, &
                      & iat)) - (tmp1 +4.d0*tmp2)*0.5d0               ! d(r^(-6))/d(r_ij)
  
  
                  drij(tauz,tauy,taux,lin(iat,iat))=drij(tauz,tauy,taux,lin(iat, &
                      & iat)) +(tmp1*alp6*t6*damp6*r/(r+rs6*R0*R0*rs8) &
                      & +3.d0*tmp2*alp8*t8*damp8*r/(r+R0*R0*rs8))*0.5d0  !d(f_dmp)/d(r_ij)
                endif
                !
                ! in dC6_rest all terms BUT C6-term is saved for the kat-loop
                !
                dc6_rest=&
                    & (s6/r6*damp6+3.d0*s8*r42/r8*damp8)*0.50d0


                disp=disp-dc6_rest*c6

                dc6i(iat)=dc6i(iat)+dc6_rest*(dc6iji+dc6ijj)
                ! if (r2.lt.crit_cn)
                dc6_rest_sum(lin(iat,iat))=dc6_rest_sum(lin(iat,iat))+dc6_rest


              else
                drij(tauz,tauy,taux,lin(iat,iat))=0.0d0
              end if


            end do
          end do
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! B E G I N jat L O O P
!!!!!!!!!!!!!!!!!!!!!!!!!!
        do jat=1,iat-1
          !
          ! get_dC6_dCNij calculates the derivative dC6(iat,jat)/dCN(iat) and
          ! dC6(iat,jat)/dCN(jat). these are saved in dC6ij for the kat loop
          !
          call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),&
              & mxc(iz(jat)),cn(iat),cn(jat),iz(iat),iz(jat),iat,jat,&
              & c6,dc6iji,dc6ijj)

          r0=r0ab(iz(jat),iz(iat))
          r42=r2r4(iz(iat))*r2r4(iz(jat))
          rcovij=rcov(iz(iat))+rcov(iz(jat))
          linij=lin(iat,jat)

          dc6ij(iat,jat)=dc6iji
          dc6ij(jat,iat)=dc6ijj
          c6save(linij)=c6
          do taux=-rep_v(1),rep_v(1)
            do tauy=-rep_v(2),rep_v(2)
              do tauz=-rep_v(3),rep_v(3)
                tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)


                rij=xyz(:,jat)-xyz(:,iat)+tau
                r2=sum(rij*rij)
                if (r2.gt.rthr) cycle


                r=dsqrt(r2)
                r6=r2*r2*r2
                r7=r6*r
                r8=r6*r2
                r9=r8*r

                !
                ! Calculates damping functions:
                if (version.eq.3) then
                  t6 = (r/(rs6*R0))**(-alp6)
                  damp6 =1.d0/( 1.d0+6.d0*t6 )
                  t8 = (r/(rs8*R0))**(-alp8)
                  damp8 =1.d0/( 1.d0+6.d0*t8 )

                  drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux,&
                      & linij)&
                      & -s6*(6.0/(r7)*C6*damp6)&
                      & -s8*(24.0/(r9)*C6*r42*damp8)

                  drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux,&
                      & linij)&
                      & +s6*C6/r7*6.d0*alp6*t6*damp6*damp6&
                      & +s8*C6*r42/r9*18.d0*alp8*t8*damp8*damp8
                else !version.eq.5
                  t6 = (r/(rs6*R0)+R0*rs8)**(-alp6)
                  damp6 =1.d0/( 1.d0+6.d0*t6 )
                  t8 = (r/(R0)+R0*rs8)**(-alp8)
                  damp8 =1.d0/( 1.d0+6.d0*t8 )
  
                  tmp1=s6*6.d0*damp6*C6/r7
                  tmp2=s8*6.d0*C6*r42*damp8/r9
                  drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux, &
                      & linij) - (tmp1 +4.d0*tmp2)  ! d(r^(-6))/d(r_ij)
  
  
                  drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux,linij) &
                      & +(tmp1*alp6*t6*damp6*r/(r+rs6*R0*R0*rs8) & 
                      & +3.d0*tmp2*alp8*t8*damp8*r/(r+R0*R0*rs8)) !d(f_dmp)/d(r_ij)
                endif
                !
                ! in dC6_rest all terms BUT C6-term is saved for the kat-loop
                !
                dc6_rest=&
                    & (s6/r6*damp6+3.d0*s8*r42/r8*damp8)


                disp=disp-dc6_rest*c6

                dc6i(iat)=dc6i(iat)+dc6_rest*dc6iji
                dc6i(jat)=dc6i(jat)+dc6_rest*dc6ijj
                ! if (r2.lt.crit_cn)
                dc6_rest_sum(linij)=dc6_rest_sum(linij)&
                    & +dc6_rest


              end do
            end do
          end do

        end do

      end do
      END IF

    elseif ((version.eq.4).or.(version.eq.6)) then



!!!!!!!!!!!!!!!!!!!!!!!
      ! NOW THE BJ Gradient !
!!!!!!!!!!!!!!!!!!!!!!!


      if (echo) write(*,*) 'doing analytical gradient O(N^2) ...'
      call pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,crit_cn)

      a1 =rs6
      a2 =rs8
      s8 =s18

      allocate(drij(-rep_v(3):rep_v(3),-rep_v(2):rep_v(2),&
          & -rep_v(1):rep_v(1),n*(n+1)/2))
      disp=0
      drij=0.0d0
      dc6_rest=0.0d0
      dc6_rest_sum=0.0d0
      c6save=0.0d0
      kat=0
      dc6i=0.0d0
      dc6ij=0.0d0

      IF ( mykey == 0 ) THEN

      do iat=na_s, na_e
        call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),&
            & mxc(iz(iat)),cn(iat),cn(iat),iz(iat),iz(iat),iat,iat,&
            & c6,dc6iji,dc6ijj)

        dc6ij(iat,iat)=dc6iji
        c6save(lin(iat,iat))=c6
        r42=r2r4(iz(iat))*r2r4(iz(iat))
        rcovij=rcov(iz(iat))+rcov(iz(iat))

        R0=a1*sqrt(3.0d0*r42)+a2

        do taux=-rep_v(1),rep_v(1)
          do tauy=-rep_v(2),rep_v(2)
            do tauz=-rep_v(3),rep_v(3)
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

              !first dE/d(tau) saved in drij(i,i,counter)
              rij=tau
              r2=sum(rij*rij)
              ! if (r2.gt.rthr) cycle

              ! if (r2.gt.0.1) then
              if (r2.gt.0.1.and.r2.lt.rthr) then
                !
                ! get_dC6_dCNij calculates the derivative dC6(iat,jat)/dCN(iat) and
                ! dC6(iat,jat)/dCN(jat). these are saved in dC6ij for the kat loop
                !
                r=dsqrt(r2)
                r4=r2*r2
                r6=r4*r2
                r7=r6*r
                r8=r6*r2
                r9=r8*r

                !
                ! Calculates damping functions:

                t6=(r6+R0**6)
                t8=(r8+R0**8)

                drij(tauz,tauy,taux,lin(iat,iat))=drij(tauz,tauy,taux,lin(iat, &
                    & iat))&
                    & -s6*C6*6.0d0*r4*r/(t6*t6)*0.5d0&
                    & -s8*C6*24.0d0*r42*r7/(t8*t8)*0.5d0


                !
                ! in dC6_rest all terms BUT C6-term is saved for the kat-loop
                !
                dc6_rest=&
                    & (s6/t6+3.d0*s8*r42/t8)*0.50d0


                disp=disp-dc6_rest*c6

                dc6i(iat)=dc6i(iat)+dc6_rest*(dc6iji+dc6ijj)
                ! if (r2.lt.crit_cn)
                dc6_rest_sum(lin(iat,iat))=dc6_rest_sum(lin(iat,iat))+&
                    & dc6_rest


              else
                drij(tauz,tauy,taux,lin(iat,iat))=0.0d0
              end if


            end do
          end do
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! B E G I N jat L O O P
!!!!!!!!!!!!!!!!!!!!!!!!!!
        do jat=1,iat-1
          !
          ! get_dC6_dCNij calculates the derivative dC6(iat,jat)/dCN(iat) and
          ! dC6(iat,jat)/dCN(jat). these are saved in dC6ij for the kat loop
          !
          call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),&
              & mxc(iz(jat)),cn(iat),cn(jat),iz(iat),iz(jat),iat,jat,&
              & c6,dc6iji,dc6ijj)

          r42=r2r4(iz(iat))*r2r4(iz(jat))
          rcovij=rcov(iz(iat))+rcov(iz(jat))

          R0=a1*dsqrt(3.0d0*r42)+a2

          linij=lin(iat,jat)
          dc6ij(iat,jat)=dc6iji
          dc6ij(jat,iat)=dc6ijj
          c6save(linij)=c6
          do taux=-rep_v(1),rep_v(1)
            do tauy=-rep_v(2),rep_v(2)
              do tauz=-rep_v(3),rep_v(3)
                tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)


                rij=xyz(:,jat)-xyz(:,iat)+tau
                r2=sum(rij*rij)
                if (r2.gt.rthr) cycle


                r=dsqrt(r2)
                r4=r2*r2
                r6=r4*r2
                r7=r6*r
                r8=r6*r2
                r9=r8*r

                !
                ! Calculates damping functions:
                t6=(r6+R0**6)
                t8=(r8+R0**8)


                drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux,&
                    & linij)&
                    & -s6*C6*6.0d0*r4*r/(t6*t6)&
                    & -s8*C6*24.0d0*r42*r7/(t8*t8)

                !
                ! in dC6_rest all terms BUT C6-term is saved for the kat-loop
                !
                dc6_rest=&
                    & (s6/t6+3.d0*s8*r42/t8)


                disp=disp-dc6_rest*c6

                dc6i(iat)=dc6i(iat)+dc6_rest*dc6iji
                dc6i(jat)=dc6i(jat)+dc6_rest*dc6ijj
                ! if (r2.lt.crit_cn)
                dc6_rest_sum(lin(iat,jat))=dc6_rest_sum(linij)&
                    & +dc6_rest


              end do
            end do
          end do

        end do

      end do
      END IF

    end if

!!!!!!!!!!!!!!!!!!!!!!!
    !! BEGIN Threebody gradient
!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.noabc) then

      ! write(*,*)'!!!!!!!!!! THREEBODY GRADIENT !!!!!!!!!!'
      sr9=0.75d0
      eabc=0.0d0
      abcthr=crit_cn
      repv=rep_cn
      ! write(*,*)'thr:',sqrt(abcthr)

      call cpu_time(time1)

      rep_cn1 = rep_cn(1) 
      rep_cn2 = rep_cn(2) 
      rep_cn3 = rep_cn(3) 
      rep_v1 = rep_v(1)
      rep_v2 = rep_v(2)
      rep_v3 = rep_v(3)

      CALL mp_sum ( c6save , intra_image_comm )
      CALL mp_sum ( dc6ij  , intra_image_comm )
      IF ( mykey == 0 ) THEN
      na_smax = max(3,na_s)

!$acc data copyin(xyz(1:3,1:n),iz(1:n),lat(1:3,1:3),r0ab(1:max_elem,1:max_elem),c6save(1:n*(n+1)),dc6ij(1:n,1:n)) &
!$acc&            copy(dc6i(1:n),drij(-rep_v3:rep_v3,-rep_v2:rep_v2,-rep_v1:rep_v1,1:n*(n+1)/2)) 
!$acc parallel vector_length(32) 
!$acc loop collapse(3) gang  private(ijvec1,ijvec2,ijvec3, ikvec1,ikvec2,ikvec3, jkvec1,jkvec2,jkvec3, c6ij,c6ik,c6jk,c9, linij,linik,linjk, &
!$acc&                              jtau1,jtau2,jtau3, rij2,rr0ij, repmin1,repmin2,repmin3,repmax1,repmax2,repmax3 ) &
!$acc&                       reduction(+:eabc) 
      do iat=na_smax,na_e
        do jat=2, n
          do kat=1, n 
            if((jat.ge.iat).or.(kat.ge.jat)) cycle    
            linij=lin(iat,jat)
            ijvec1=xyz(1,jat)-xyz(1,iat)
            ijvec2=xyz(2,jat)-xyz(2,iat)
            ijvec3=xyz(3,jat)-xyz(3,iat)

            c6ij=c6save(linij)

            linik=lin(iat,kat)
            linjk=lin(jat,kat)
            ikvec1=xyz(1,kat)-xyz(1,iat)
            ikvec2=xyz(2,kat)-xyz(2,iat)
            ikvec3=xyz(3,kat)-xyz(3,iat)
            jkvec1=xyz(1,kat)-xyz(1,jat)
            jkvec2=xyz(2,kat)-xyz(2,jat)
            jkvec3=xyz(3,kat)-xyz(3,jat)

            c6ik=c6save(linik)
            c6jk=c6save(linjk)
            c9=-1.0d0*dsqrt(c6ij*c6ik*c6jk)
!$acc loop seq independent
            do jtaux=-rep_cn1,rep_cn1
              repmin1=max(-rep_cn1,jtaux-rep_cn1)
              repmax1=min(rep_cn1,jtaux+rep_cn1)
!$acc loop seq independent   
              do jtauy=-rep_cn2,rep_cn2
                repmin2=max(-rep_cn2,jtauy-rep_cn2)
                repmax2=min(rep_cn2,jtauy+rep_cn2)
!$acc loop seq independent 
                do jtauz=-rep_cn3,rep_cn3
                  repmin3=max(-rep_cn3,jtauz-rep_cn3)
                  repmax3=min(rep_cn3,jtauz+rep_cn3)
                  jtau1=jtaux*lat(1,1)+jtauy*lat(1,2)+jtauz*lat(1,3)
                  jtau2=jtaux*lat(2,1)+jtauy*lat(2,2)+jtauz*lat(2,3)
                  jtau3=jtaux*lat(3,1)+jtauy*lat(3,2)+jtauz*lat(3,3)
                  rij2= (ijvec1+jtau1)*(ijvec1+jtau1) + (ijvec2+jtau2)*(ijvec2+jtau2) + (ijvec3+jtau3)*(ijvec3+jtau3) 
                  if (rij2.gt.abcthr)cycle

                  rr0ij=DSQRT(rij2)/r0ab(iz(iat),iz(jat))

!$acc loop vector collapse(3) private(ktau1,ktau2,ktau3, dumvec1,dumvec2,dumvec3, rik2,rjk2,rr0ik,rr0jk, &
!$acc&                                geomean,geomean2,geomean3,r0av,damp9,ang,dc6_rest,dfdmp,r,dang,tmp1,dc9) &
!$acc&                        reduction(+:eabc) 
                  do ktaux=repmin1,repmax1
                    do ktauy=repmin2,repmax2
                      do ktauz=repmin3,repmax3
                        ktau1=ktaux*lat(1,1)+ktauy*lat(1,2)+ktauz*lat(1,3)
                        ktau2=ktaux*lat(2,1)+ktauy*lat(2,2)+ktauz*lat(2,3)
                        ktau3=ktaux*lat(3,1)+ktauy*lat(3,2)+ktauz*lat(3,3)
                        rik2=(ikvec1+ktau1)*(ikvec1+ktau1)+(ikvec2+ktau2)*(ikvec2+ktau2)+(ikvec3+ktau3)*(ikvec3+ktau3)
                        if (rik2.gt.abcthr)cycle

                        dumvec1=jkvec1+ktau1-jtau1
                        dumvec2=jkvec2+ktau2-jtau2
                        dumvec3=jkvec3+ktau3-jtau3
                        rjk2=dumvec1*dumvec1+dumvec2*dumvec2+dumvec3*dumvec3
                        if (rjk2.gt.abcthr)cycle
                        rr0ik=dsqrt(rik2)/r0ab(iz(iat),iz(kat))
                        rr0jk=dsqrt(rjk2)/r0ab(iz(jat),iz(kat))
                        geomean2=(rij2*rjk2*rik2)
                        ! first calculate the three components for the energy calculation fdmp
                        ! and ang
                        r0av=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
                        damp9=1./(1.+6.*(sr9*r0av)**alp9)

                        geomean=dsqrt(geomean2)
                        geomean3=geomean*geomean2
                        ang=0.375d0*(rij2+rjk2-rik2)*(rij2-rjk2+rik2)&
                            & *(-rij2+rjk2+rik2)/(geomean3*geomean2)&
                            & +1.0d0/(geomean3)

                        dc6_rest=ang*damp9
                        eabc=eabc+dc6_rest*c9
                        !
                        !start calculating the gradient components dfdmp, dang and dc9

                        !dfdmp is the same for all three distances
                        dfdmp=2.d0*alp9*(0.75d0*r0av)**(alp9)*damp9*damp9

                        !start calculating the derivatives of each part w.r.t. r_ij
                        r=dsqrt(rij2)


                        dang=-0.375d0*(rij2**3+rij2**2*(rjk2+rik2)&
                            & +rij2*(3.0d0*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2)&
                            & -5.0*(rjk2-rik2)**2*(rjk2+rik2))&
                            & /(r*geomean3*geomean2)

                        tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
!$acc atomic update
                        drij(jtauz,jtauy,jtaux,linij)= drij(jtauz,jtauy,jtaux,linij)-tmp1
!$acc end atomic 
                        !start calculating the derivatives of each part w.r.t. r_ik

                        r=dsqrt(rik2)


                        dang=-0.375d0*(rik2**3+rik2**2*(rjk2+rij2)&
                            & +rik2*(3.0d0*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2)&
                            & -5.0*(rjk2-rij2)**2*(rjk2+rij2))&
                            & /(r*geomean3*geomean2)

                        tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                        ! tmp1=-dc9
!$acc atomic update
                        drij(ktauz,ktauy,ktaux,linik)= drij(ktauz,ktauy,ktaux,linik)-tmp1
!$acc end atomic

                        !
                        !start calculating the derivatives of each part w.r.t. r_jk

                        r=dsqrt(rjk2)

                        dang=-0.375d0*(rjk2**3+rjk2**2*(rik2+rij2)&
                            & +rjk2*(3.0d0*rik2**2+2.0*rik2*rij2+3.0*rij2**2)&
                            & -5.0*(rik2-rij2)**2*(rik2+rij2))&
                            & /(r*geomean3*geomean2)

                        tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
!$acc atomic update
                        drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)= drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)-tmp1
!$acc end atomic

                        !calculating the CN derivative dE_disp(ijk)/dCN(i)

                        dc9=dc6ij(iat,jat)/c6ij+dc6ij(iat,kat)/c6ik
                        dc9=0.5d0*c9*dc9
!$acc atomic update 
                        dc6i(iat) = dc6i(iat) + dc6_rest*dc9   
!$acc end atomic

                        dc9=dc6ij(jat,iat)/c6ij+dc6ij(jat,kat)/c6jk
                        dc9=0.5d0*c9*dc9
!$acc atomic update 
                        dc6i(jat) = dc6i(jat) + dc6_rest*dc9  
!$acc end atomic

                        dc9=dc6ij(kat,iat)/c6ik+dc6ij(kat,jat)/c6jk
                        dc9=0.5d0*c9*dc9
!$acc atomic update 
                        dc6i(kat) = dc6i(kat) + dc6_rest*dc9   
!$acc end atomic
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
!$acc end parallel 
!$acc end data 

      ! Now the interaction with jat=iat of the triples iat,iat,kat
      do iat=max(2,na_s),na_e
        jat=iat
        linij=lin(iat,jat)
        ijvec=0.0d0

        c6ij=c6save(linij)
        do kat=1,iat-1
          linjk=lin(jat,kat)
          linik=linjk

          c6ik=c6save(linik)
          c6jk=c6ik
          ikvec=xyz(:,kat)-xyz(:,iat)
          jkvec=ikvec
          c9=-dsqrt(c6ij*c6ik*c6jk)
          do jtaux=-repv(1),repv(1)
            repmin(1)=max(-repv(1),jtaux-repv(1))
            repmax(1)=min(repv(1),jtaux+repv(1))
            do jtauy=-repv(2),repv(2)
              repmin(2)=max(-repv(2),jtauy-repv(2))
              repmax(2)=min(repv(2),jtauy+repv(2))
              do jtauz=-repv(3),repv(3)
                repmin(3)=max(-repv(3),jtauz-repv(3))
                repmax(3)=min(repv(3),jtauz+repv(3))
                if (jtaux.eq.0 .and. jtauy.eq.0 .and. jtauz.eq.0) cycle
                jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
                dumvec=jtau
                rij2=SUM(dumvec*dumvec)
                if (rij2.gt.abcthr)cycle

                rr0ij=DSQRT(rij2)/r0ab(iz(iat),iz(jat))

                do ktaux=repmin(1),repmax(1)
                  do ktauy=repmin(2),repmax(2)
                    do ktauz=repmin(3),repmax(3)
                      ! every result * 0.5

                      ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
                      dumvec=ikvec+ktau
                      dumvec=dumvec*dumvec
                      rik2=SUM(dumvec)
                      if (rik2.gt.abcthr)cycle

                      dumvec=jkvec+ktau-jtau
                      dumvec=dumvec*dumvec
                      rjk2=SUM(dumvec)
                      if (rjk2.gt.abcthr)cycle
                      rr0ik=DSQRT(rik2)/r0ab(iz(iat),iz(kat))
                      rr0jk=DSQRT(rjk2)/r0ab(iz(jat),iz(kat))


                      geomean2=(rij2*rjk2*rik2)
                      r0av=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
                      damp9=1./(1.+6.*(sr9*r0av)**alp9)

                      geomean=dsqrt(geomean2)
                      geomean3=geomean*geomean2
                      ang=0.375d0*(rij2+rjk2-rik2)*(rij2-rjk2+rik2)&
                          & *(-rij2+rjk2+rik2)/(geomean3*geomean2)&
                          & +1.0d0/(geomean3)


                      dc6_rest=ang*damp9/2.0d0
                      eabc=eabc+dc6_rest*c9

                      ! iat=jat
                      dfdmp=2.d0*alp9*(0.75d0*r0av)**(alp9)*damp9*damp9

                      !start calculating the derivatives of each part w.r.t. r_ij
                      r=dsqrt(rij2)

                      dang=-0.375d0*(rij2**3+rij2**2*(rjk2+rik2) &
                          & +rij2*(3.0d0*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2)&
                          & -5.0*(rjk2-rik2)**2*(rjk2+rik2))&
                          & /(r*geomean3*geomean2)

                      tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                      drij(jtauz,jtauy,jtaux,linij)=&
                          & drij(jtauz,jtauy,jtaux,linij)-tmp1/2.0

                      !start calculating the derivatives of each part w.r.t. r_ik
                      r=dsqrt(rik2)


                      dang=-0.375d0*(rik2**3+rik2**2*(rjk2+rij2)&
                          & +rik2*(3.0d0*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2)&
                          & -5.0*(rjk2-rij2)**2*(rjk2+rij2))&
                          & /(r*geomean3*geomean2)

                      tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                      drij(ktauz,ktauy,ktaux,linik)=&
                          & drij(ktauz,ktauy,ktaux,linik)-tmp1/2.0
                      !
                      !start calculating the derivatives of each part w.r.t. r_ik
                      r=dsqrt(rjk2)

                      dang=-0.375d0*(rjk2**3+rjk2**2*(rik2+rij2)&
                          & +rjk2*(3.0d0*rik2**2+2.0*rik2*rij2+3.0*rij2**2)&
                          & -5.0*(rik2-rij2)**2*(rik2+rij2))&
                          & /(r*geomean3*geomean2)

                      tmp1=-dang*c9*damp9+dfdmp/r*c9*ang

                      drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)=&
                          & drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)-tmp1/2.0

                      dc9=dc6ij(iat,jat)/c6ij+dc6ij(iat,kat)/c6ik
                      dc9=0.5d0*c9*dc9
                      dc6i(iat)=dc6i(iat)+dc6_rest*dc9

                      dc9=dc6ij(jat,iat)/c6ij+dc6ij(jat,kat)/c6jk
                      dc9=0.5d0*c9*dc9
                      dc6i(jat)=dc6i(jat)+dc6_rest*dc9

                      dc9=dc6ij(kat,iat)/c6ik+dc6ij(kat,jat)/c6jk
                      dc9=0.5d0*c9*dc9
                      dc6i(kat)=dc6i(kat)+dc6_rest*dc9




                    end do
                  end do
                end do

              end do
            end do
          end do
        end do
      end do

      do iat=max(2,na_s),na_e
        do jat=1,iat-1
          kat=jat
          linij=lin(iat,jat)
          linjk=lin(jat,kat)
          linik=linij

          c6ij=c6save(linij)
          c6ik=c6ij

          c6jk=c6save(linjk)
          ikvec=xyz(:,kat)-xyz(:,iat)
          ijvec=ikvec
          jkvec=0.0d0

          c9=-1.0d0*dsqrt(c6ij*c6ik*c6jk)
          do jtaux=-repv(1),repv(1)
            repmin(1)=max(-repv(1),jtaux-repv(1))
            repmax(1)=min(repv(1),jtaux+repv(1))
            do jtauy=-repv(2),repv(2)
              repmin(2)=max(-repv(2),jtauy-repv(2))
              repmax(2)=min(repv(2),jtauy+repv(2))
              do jtauz=-repv(3),repv(3)
                repmin(3)=max(-repv(3),jtauz-repv(3))
                repmax(3)=min(repv(3),jtauz+repv(3))

                jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
                dumvec=ijvec+jtau
                dumvec=dumvec*dumvec
                rij2=SUM(dumvec)
                if (rij2.gt.abcthr)cycle

                rr0ij=SQRT(rij2)/r0ab(iz(iat),iz(jat))

                do ktaux=repmin(1),repmax(1)
                  do ktauy=repmin(2),repmax(2)
                    do ktauz=repmin(3),repmax(3)
                      ! every result * 0.5
                      if (jtaux.eq.ktaux .and. jtauy.eq.ktauy&
                          & .and. jtauz.eq.ktauz) cycle
                      ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
                      dumvec=ikvec+ktau
                      dumvec=dumvec*dumvec
                      rik2=SUM(dumvec)
                      if (rik2.gt.abcthr)cycle
                      rr0ik=SQRT(rik2)/r0ab(iz(iat),iz(kat))

                      dumvec=jkvec+ktau-jtau
                      dumvec=dumvec*dumvec
                      rjk2=SUM(dumvec)
                      if (rjk2.gt.abcthr)cycle
                      rr0jk=SQRT(rjk2)/r0ab(iz(jat),iz(kat))

                      ! if (rij*rjk*rik.gt.abcthr)cycle

                      geomean2=(rij2*rjk2*rik2)
                      r0av=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
                      damp9=1./(1.+6.d0*(sr9*r0av)**alp9)

                      geomean=dsqrt(geomean2)
                      geomean3=geomean*geomean2
                      ang=0.375d0*(rij2+rjk2-rik2)*(rij2-rjk2+rik2)&
                          & *(-rij2+rjk2+rik2)/(geomean2*geomean3)&
                          & +1.0d0/(geomean3)
                      dc6_rest=ang*damp9/2.0d0
                      eabc=eabc+dc6_rest*c9


                      ! jat=kat
                      dfdmp=2.d0*alp9*(0.75d0*r0av)**(alp9)*damp9*damp9
                      !start calculating the derivatives of each part w.r.t. r_ij
                      r=dsqrt(rij2)

                      dang=-0.375d0*(rij2**3+rij2**2*(rjk2+rik2)&
                          & +rij2*(3.0d0*rjk2**2+2.0d0*rjk2*rik2+3.0d0*rik2**2)&
                          & -5.0d0*(rjk2-rik2)**2*(rjk2+rik2))&
                          & /(r*geomean3*geomean2)

                      tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                      drij(jtauz,jtauy,jtaux,linij)=&
                          & drij(jtauz,jtauy,jtaux,linij)-tmp1/2.0d0

                      !start calculating the derivatives of each part w.r.t. r_ik
                      r=dsqrt(rik2)


                      dang=-0.375d0*(rik2**3+rik2**2*(rjk2+rij2)&
                          & +rik2*(3.0d0*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2)&
                          & -5.0*(rjk2-rij2)**2*(rjk2+rij2))&
                          & /(r*geomean3*geomean2)

                      tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                      ! tmp1=-dc9
                      drij(ktauz,ktauy,ktaux,linik)=&
                          & drij(ktauz,ktauy,ktaux,linik)-tmp1/2.0d0
                      !
                      !start calculating the derivatives of each part w.r.t. r_jk
                      r=dsqrt(rjk2)

                      dang=-0.375d0*(rjk2**3+rjk2**2*(rik2+rij2)&
                          & +rjk2*(3.0d0*rik2**2+2.0*rik2*rij2+3.0*rij2**2)&
                          & -5.0d0*(rik2-rij2)**2*(rik2+rij2))&
                          & /(r*geomean3*geomean2)

                      tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                      drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)=&
                          & drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)-tmp1/2.0d0

                      !calculating the CN derivative dE_disp(ijk)/dCN(i)

                      dc9=dc6ij(iat,jat)/c6ij+dc6ij(iat,kat)/c6ik
                      dc9=0.5d0*c9*dc9
                      dc6i(iat)=dc6i(iat)+dc6_rest*dc9

                      dc9=dc6ij(jat,iat)/c6ij+dc6ij(jat,kat)/c6jk
                      dc9=0.5d0*c9*dc9
                      dc6i(jat)=dc6i(jat)+dc6_rest*dc9

                      dc9=dc6ij(kat,iat)/c6ik+dc6ij(kat,jat)/c6jk
                      dc9=0.5d0*c9*dc9
                      dc6i(kat)=dc6i(kat)+dc6_rest*dc9




                    end do
                  end do
                end do

              end do
            end do
          end do
        end do
      end do


      ! And finally the self interaction iat=jat=kat all

      idum=0
      do iat=na_s,na_e
        jat=iat
        kat=iat
        ijvec=0.0d0
        linij=lin(iat,jat)
        linik=lin(iat,kat)
        linjk=lin(jat,kat)
        ikvec=ijvec
        jkvec=ikvec
        c6ij=c6save(linij)
        c6ik=c6ij
        c6jk=c6ij
        c9=-(DSQRT(c6ij*c6ij*c6ij))

        do jtaux=-repv(1),repv(1)
          repmin(1)=max(-repv(1),jtaux-repv(1))
          repmax(1)=min(repv(1),jtaux+repv(1))
          do jtauy=-repv(2),repv(2)
            repmin(2)=max(-repv(2),jtauy-repv(2))
            repmax(2)=min(repv(2),jtauy+repv(2))
            do jtauz=-repv(3),repv(3)
              repmin(3)=max(-repv(3),jtauz-repv(3))
              repmax(3)=min(repv(3),jtauz+repv(3))
              if ((jtaux.eq.0) .and.(jtauy.eq.0) .and.(jtauz.eq.0))cycle
              jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
              dumvec=jtau
              dumvec=dumvec*dumvec
              rij2=SUM(dumvec)
              if (rij2.gt.abcthr)cycle
              rr0ij=SQRT(rij2)/r0ab(iz(iat),iz(jat))

              do ktaux=repmin(1),repmax(1)
                do ktauy=repmin(2),repmax(2)
                  do ktauz=repmin(3),repmax(3)
                    if ((ktaux.eq.0) .and.( ktauy.eq.0) .and.( ktauz.eq.0))cycle
                    if ((ktaux.eq.jtaux) .and. (ktauy.eq.jtauy)&
                        & .and. (ktauz.eq.jtauz)) cycle

                    ! every result * 1/6 becaues every triple is counted twice due to the tw
                    !
                    !plus 1/3 becaues every triple is three times in each unitcell
                    ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
                    dumvec=ktau
                    dumvec=dumvec*dumvec
                    rik2=SUM(dumvec)
                    if (rik2.gt.abcthr)cycle
                    rr0ik=SQRT(rik2)/r0ab(iz(iat),iz(kat))

                    dumvec=jkvec+ktau-jtau
                    dumvec=dumvec*dumvec
                    rjk2=SUM(dumvec)
                    if (rjk2.gt.abcthr)cycle
                    rr0jk=SQRT(rjk2)/r0ab(iz(jat),iz(kat))

                    geomean2=(rij2*rjk2*rik2)
                    r0av=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
                    damp9=1./(1.+6.*(sr9*r0av)**alp9)

                    geomean=dsqrt(geomean2)
                    geomean3=geomean*geomean2
                    ang=0.375d0*(rij2+rjk2-rik2)*(rij2-rjk2+rik2)&
                        & *(-rij2+rjk2+rik2)/(geomean2*geomean3)&
                        & +1.0d0/(geomean3)
                    dc6_rest=ang*damp9/6.0d0
                    eabc=eabc+c9*dc6_rest

                    ! iat=jat=kat
                    dfdmp=2.d0*alp9*(0.75d0*r0av)**(alp9)*damp9*damp9
                    !start calculating the derivatives of each part w.r.t. r_ij

                    r=dsqrt(rij2)
                    dang=-0.375d0*(rij2**3+rij2**2*(rjk2+rik2)&
                        & +rij2*(3.0d0*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2)&
                        & -5.0*(rjk2-rik2)**2*(rjk2+rik2))&
                        & /(r*geomean3*geomean2)


                    tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                    drij(jtauz,jtauy,jtaux,linij)=&
                        & drij(jtauz,jtauy,jtaux,linij)-tmp1/6.0d0

                    !start calculating the derivatives of each part w.r.t. r_ik

                    r=dsqrt(rik2)

                    dang=-0.375d0*(rik2**3+rik2**2*(rjk2+rij2)&
                        & +rik2*(3.0d0*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2)&
                        & -5.0*(rjk2-rij2)**2*(rjk2+rij2))&
                        & /(r*geomean3*geomean2)

                    tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                    drij(ktauz,ktauy,ktaux,linik)=&
                        & drij(ktauz,ktauy,ktaux,linik)-tmp1/6.0d0
                    !
                    !start calculating the derivatives of each part w.r.t. r_jk

                    r=dsqrt(rjk2)
                    dang=-0.375d0*(rjk2**3+rjk2**2*(rik2+rij2)&
                        & +rjk2*(3.0d0*rik2**2+2.0*rik2*rij2+3.0*rij2**2)&
                        & -5.0*(rik2-rij2)**2*(rik2+rij2))&
                        & /(r*geomean3*geomean2)

                    tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                    drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)=&
                        & drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)-tmp1/6.0d0


                    !calculating the CN derivative dE_disp(ijk)/dCN(i)

                    dc9=dc6ij(iat,jat)/c6ij+dc6ij(iat,kat)/c6ik
                    dc9=0.5d0*c9*dc9
                    dc6i(iat)=dc6i(iat)+dc6_rest*dc9

                    dc9=dc6ij(jat,iat)/c6ij+dc6ij(jat,kat)/c6jk
                    dc9=0.5d0*c9*dc9
                    dc6i(jat)=dc6i(jat)+dc6_rest*dc9

                    dc9=dc6ij(kat,iat)/c6ik+dc6ij(kat,jat)/c6jk
                    dc9=0.5d0*c9*dc9
                    dc6i(kat)=dc6i(kat)+dc6_rest*dc9





                  end do
                end do
              end do
            end do
          end do
          !jtaux
        end do

      end do

      END IF
      CALL mp_sum ( eabc , intra_image_comm )

      call cpu_time(time2)

      ! write(*,*)' eabc(gdisp): ',eabc
      ! write(*,'('' time(abc) '',f6.1)')time2-time1
      disp=disp-eabc
      ! write(*,*)'gdisp:',disp
    end if

    CALL mp_sum ( drij , intra_image_comm )
    CALL mp_sum ( dc6i , intra_image_comm )

    sigma_abc=0.0d0
    sigma=0.0d0

    ! After calculating all derivatives dE/dr_ij w.r.t. distances,
    ! the grad w.r.t. the coordinates is calculated dE/dr_ij * dr_ij/dxyz_i
    do iat=2,n
      do jat=1,iat-1
        linij=lin(iat,jat)
        rcovij=rcov(iz(iat))+rcov(iz(jat))
        do taux=-rep_v(1),rep_v(1)
          do tauy=-rep_v(2),rep_v(2)
            do tauz=-rep_v(3),rep_v(3)
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

              rij=xyz(:,jat)-xyz(:,iat)+tau
              r2=sum(rij*rij)
              if (r2.gt.rthr.or.r2.lt.0.5) cycle
              r=dsqrt(r2)

              if (r2.lt.crit_cn) then
                expterm=exp(-k1*(rcovij/r-1.d0))
                dcnn=-k1*rcovij*expterm/&
                    & (r2*(expterm+1.d0)*(expterm+1.d0))
              else
                dcnn=0.0d0
              end if
              x1=drij(tauz,tauy,taux,linij)+dcnn*(dc6i(iat)+dc6i(jat))
              vec=x1*rij/r
              g(:,iat)=g(:,iat)+vec
              g(:,jat)=g(:,jat)-vec
              do i=1,3
                do j=1,3
                  sigma(j,i)=sigma(j,i)+vec(j)*rij(i)
                end do
              end do



            end do
          end do
        end do
      end do
    end do

    do iat=1,n
      rcovij=rcov(iz(iat))+rcov(iz(iat))
      do taux=-rep_v(1),rep_v(1)
        do tauy=-rep_v(2),rep_v(2)
          do tauz=-rep_v(3),rep_v(3)
            if (taux.eq.0.and.tauy.eq.0.and.tauz.eq.0) cycle

            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
            r2=(sum(tau*tau))
            r=dsqrt(r2)
            if (r2.lt.crit_cn) then
              expterm=exp(-k1*(rcovij/r-1.d0))
              dcnn=-k1*rcovij*expterm/&
                  & (r2*(expterm+1.d0)*(expterm+1.d0))
            else
              dcnn=0.0d0
            end if
            x1=drij(tauz,tauy,taux,lin(iat,iat))+dcnn*dc6i(iat)
            vec=x1*tau/r
            vec2(1)=taux
            vec2(2)=tauy
            vec2(3)=tauz
            do i=1,3
              do j=1,3
                sigma(j,i)=sigma(j,i)+vec(j)*tau(i)
              end do
            end do


          end do
        end do
      end do



    end do



    stress=0.0d0
    glatabc=0.0d0
    call inv_cell(lat,lat_1)
    do a=1,3
      do b=1,3
        do my=1,3
          stress(a,b)=stress(a,b)-sigma(a,my)*lat_1(b,my)
        end do
      end do
    end do



    ! write(*,*)'drij:',drij(lin(iat,jat),:)
    ! write(*,*)'g:',g(1,1:3)
    ! write(*,*)'dcn:',sum(dcn(lin(2,1),:))



    deallocate(drij)




999 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !This is where the D2 gradient and the numerical gradient jump.
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do i=1,n
    ! write(*,'(83F17.12)') g(1:3,i)
    ! end do
    gnorm=sum(abs(g(1:3,1:n)))
    if (echo)then
      ! write(*,*)'testsum:',testsum*autoev/autoang
      write(*,*)'|G(force)| =',gnorm
      gnorm=sum(abs(stress(1:3,1:3)))
      write(*,*)'|G(stress)|=',gnorm
    end if

  end subroutine pbcgdisp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine copyc6(fname,maxc,max_elem,c6ab,maxci, minc6,minc6list,maxc6,&
      & maxc6list)
    integer maxc,max_elem,maxci(max_elem),mima
    real(wp)  c6ab(max_elem,max_elem,maxc,maxc,3)
    character*(*) fname
    logical minc6,maxc6,minc6list(max_elem),maxc6list(max_elem)

    character*1  atmp
    character*80 btmp
    real(wp)  x,y,f,cn1,cn2,cmax,xx(10)
    integer iat,jat,i,n,l,j,k,il,iadr,jadr,nn,kk
    logical special

    call init_pars()
    c6ab=-1
    maxci=0
    ! read file
    kk=1
    !only use values for cn=minimum
    if (minc6.or.maxc6) then
      do i=1,94
        if (minc6list(i))then
          c6ab(i,:,1,:,2)=10000000.0
          c6ab(:,i,:,1,3)=10000000.0
        end if
      end do


      do nn=1,nlines
        special=.false.
        iat=int(pars(kk+1))
        jat=int(pars(kk+2))
        call limit(iat,jat,iadr,jadr)

        !only CN=minimum for iat
        if (minc6list(iat)) then
          special=.true.
          maxci(iat)=1
          maxci(jat)=max(maxci(jat),jadr)

          if (pars(kk+3).le.c6ab(iat,jat,1,jadr,2)) then

            c6ab(iat,jat,1,jadr,1)=pars(kk)
            c6ab(iat,jat,1,jadr,2)=pars(kk+3)
            c6ab(iat,jat,1,jadr,3)=pars(kk+4)

            c6ab(jat,iat,jadr,1,1)=pars(kk)
            c6ab(jat,iat,jadr,1,2)=pars(kk+4)
            c6ab(jat,iat,jadr,1,3)=pars(kk+3)
          end if
        end if

        !only CN=minimum for jat
        if (minc6list(jat)) then
          special=.true.
          maxci(iat)=max(maxci(iat),iadr)
          maxci(jat)=1

          if (pars(kk+4).le.c6ab(iat,jat,iadr,1,3)) then

            c6ab(iat,jat,iadr,1,1)=pars(kk)
            c6ab(iat,jat,iadr,1,2)=pars(kk+3)
            c6ab(iat,jat,iadr,1,3)=pars(kk+4)

            c6ab(jat,iat,1,iadr,1)=pars(kk)
            c6ab(jat,iat,1,iadr,2)=pars(kk+4)
            c6ab(jat,iat,1,iadr,3)=pars(kk+3)
          end if
        end if


        !only CN=minimum for
        if (minc6list(iat).and.minc6list(jat)) then
          special=.true.
          maxci(jat)=1
          maxci(iat)=1

          if (pars(kk+4).le.c6ab(iat,jat,1,1,3).and.&
              &      pars(kk+3).le.c6ab(iat,jat,1,1,2)) then

            c6ab(iat,jat,1,1,1)=pars(kk)
            c6ab(iat,jat,1,1,2)=pars(kk+3)
            c6ab(iat,jat,1,1,3)=pars(kk+4)

            c6ab(jat,iat,1,1,1)=pars(kk)
            c6ab(jat,iat,1,1,2)=pars(kk+4)
            c6ab(jat,iat,1,1,3)=pars(kk+3)
          end if
        end if



        !only CN=maximum for iat
        if (maxc6list(iat)) then
          special=.true.

          maxci(iat)=1
          maxci(jat)=max(maxci(jat),jadr)

          if (pars(kk+3).ge.c6ab(iat,jat,1,jadr,2)) then

            c6ab(iat,jat,1,jadr,1)=pars(kk)
            c6ab(iat,jat,1,jadr,2)=pars(kk+3)
            c6ab(iat,jat,1,jadr,3)=pars(kk+4)

            c6ab(jat,iat,jadr,1,1)=pars(kk)
            c6ab(jat,iat,jadr,1,2)=pars(kk+4)
            c6ab(jat,iat,jadr,1,3)=pars(kk+3)
          end if
        end if
        !only CN=maximum for jat
        if (maxc6list(jat)) then
          special=.true.

          maxci(jat)=1
          maxci(iat)=max(maxci(iat),iadr)

          if (pars(kk+4).ge.c6ab(iat,jat,iadr,1,3)) then

            c6ab(iat,jat,iadr,1,1)=pars(kk)
            c6ab(iat,jat,iadr,1,2)=pars(kk+3)
            c6ab(iat,jat,iadr,1,3)=pars(kk+4)

            c6ab(jat,iat,1,iadr,1)=pars(kk)
            c6ab(jat,iat,1,iadr,2)=pars(kk+4)
            c6ab(jat,iat,1,iadr,3)=pars(kk+3)
          end if
        end if

        !only CN=maximum for
        if (maxc6list(iat).and.maxc6list(jat)) then
          special=.true.
          maxci(jat)=1
          maxci(iat)=1

          if (pars(kk+4).ge.c6ab(iat,jat,1,1,3).and.&
              &      pars(kk+3).ge.c6ab(iat,jat,1,1,2)) then

            c6ab(iat,jat,1,1,1)=pars(kk)
            c6ab(iat,jat,1,1,2)=pars(kk+3)
            c6ab(iat,jat,1,1,3)=pars(kk+4)

            c6ab(jat,iat,1,1,1)=pars(kk)
            c6ab(jat,iat,1,1,2)=pars(kk+4)
            c6ab(jat,iat,1,1,3)=pars(kk+3)
          end if
        end if

        !only CN=minimum for
        if (minc6list(iat).and.maxc6list(jat)) then
          !and CN=maximum jat
          special=.true.
          maxci(jat)=1
          maxci(iat)=1

          if (pars(kk+4).ge.c6ab(iat,jat,1,1,3).and.&
              &      pars(kk+3).le.c6ab(iat,jat,1,1,2)) then

            c6ab(iat,jat,1,1,1)=pars(kk)
            c6ab(iat,jat,1,1,2)=pars(kk+3)
            c6ab(iat,jat,1,1,3)=pars(kk+4)

            c6ab(jat,iat,1,1,1)=pars(kk)
            c6ab(jat,iat,1,1,2)=pars(kk+4)
            c6ab(jat,iat,1,1,3)=pars(kk+3)
          end if
        end if

        !only CN=maximum for
        if (maxc6list(iat).and.minc6list(jat)) then
          !  and CN=minimum fo
          special=.true.
          maxci(jat)=1
          maxci(iat)=1

          if (pars(kk+4).le.c6ab(iat,jat,1,1,3).and.&
              &      pars(kk+3).ge.c6ab(iat,jat,1,1,2)) then

            c6ab(iat,jat,1,1,1)=pars(kk)
            c6ab(iat,jat,1,1,2)=pars(kk+3)
            c6ab(iat,jat,1,1,3)=pars(kk+4)

            c6ab(jat,iat,1,1,1)=pars(kk)
            c6ab(jat,iat,1,1,2)=pars(kk+4)
            c6ab(jat,iat,1,1,3)=pars(kk+3)
          end if
        end if

        if (.not.special) then

          maxci(iat)=max(maxci(iat),iadr)
          maxci(jat)=max(maxci(jat),jadr)

          c6ab(iat,jat,iadr,jadr,1)=pars(kk)
          c6ab(iat,jat,iadr,jadr,2)=pars(kk+3)
          c6ab(iat,jat,iadr,jadr,3)=pars(kk+4)

          c6ab(jat,iat,jadr,iadr,1)=pars(kk)
          c6ab(jat,iat,jadr,iadr,2)=pars(kk+4)
          c6ab(jat,iat,jadr,iadr,3)=pars(kk+3)
        end if
        kk=(nn*5)+1
      end do



      !no min/max at all
    else
      do nn=1,nlines
        iat=int(pars(kk+1))
        jat=int(pars(kk+2))
        call limit(iat,jat,iadr,jadr)
        maxci(iat)=max(maxci(iat),iadr)
        maxci(jat)=max(maxci(jat),jadr)

        c6ab(iat,jat,iadr,jadr,1)=pars(kk)
        c6ab(iat,jat,iadr,jadr,2)=pars(kk+3)
        c6ab(iat,jat,iadr,jadr,3)=pars(kk+4)

        c6ab(jat,iat,jadr,iadr,1)=pars(kk)
        c6ab(jat,iat,jadr,iadr,2)=pars(kk+4)
        c6ab(jat,iat,jadr,iadr,3)=pars(kk+3)
        kk=(nn*5)+1
      end do
    end if

  end subroutine copyc6



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SET_CRITERIA(rthr,lat,tau_max)

    REAL(WP) :: r_cutoff,rthr
    REAL(WP) :: lat(3,3)
    REAL(WP) :: tau_max(3)
    REAL(WP) :: norm1(3),norm2(3),norm3(3)
    REAL(WP) :: cos10,cos21,cos32

    r_cutoff=sqrt(rthr)
    ! write(*,*) 'lat',lat
    call kreuzprodukt(lat(:,2),lat(:,3),norm1)
    call kreuzprodukt(lat(:,3),lat(:,1),norm2)
    call kreuzprodukt(lat(:,1),lat(:,2),norm3)
    ! write(*,*) 'norm2',norm2
    norm1=norm1/VECTORSIZE(norm1)
    norm2=norm2/VECTORSIZE(norm2)
    norm3=norm3/VECTORSIZE(norm3)
    ! write(*,*) 'norm2_',norm2
    cos10=SUM(norm1*lat(:,1))
    cos21=SUM(norm2*lat(:,2))
    cos32=SUM(norm3*lat(:,3))
    tau_max(1)=abs(r_cutoff/cos10)
    tau_max(2)=abs(r_cutoff/cos21)
    tau_max(3)=abs(r_cutoff/cos32)
    ! write(*,'(3f8.4)')tau_max(1),tau_max(2),tau_max(3)
  end subroutine SET_CRITERIA


  subroutine kreuzprodukt(A,B,C)

    REAL(WP) :: A(3),B(3)
    REAL(WP) :: X,Y,Z
    REAL(WP) :: C(3)

    X=A(2)*B(3)-B(2)*A(3)
    Y=A(3)*B(1)-B(3)*A(1)
    Z=A(1)*B(2)-B(1)*A(2)
    C=(/X,Y,Z/)
  end subroutine kreuzprodukt

  function VECTORSIZE(VECT)

    REAL(WP) :: VECT(3)
    REAL(WP) :: SVECT(3)
    REAL(WP) :: VECTORSIZE

    SVECT=VECT*VECT
    VECTORSIZE=SUM(SVECT)
    VECTORSIZE=VECTORSIZE**(0.5)
  end function VECTORSIZE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function determinant(aa) result(det)
    real(wp), intent(in) :: aa(:,:)
    real(wp) :: det

    det = aa(1,1) * (aa(2,2) * aa(3,3) - aa(3,2) * aa(2,3))&
        & - aa(1,2) * (aa(2,1) * aa(3,3) - aa(3,1) * aa(2,3))&
        & + aa(1,3) * (aa(2,1) * aa(3,2) - aa(3,1) * aa(2,2))
    
  end function determinant


  subroutine inv_cell(x,a)
    real(wp), intent(in) :: x(3,3)
    real(wp), intent(out) :: a(3,3)
    integer i
    real(wp) det

    a = 0.0
    det = determinant(x)
    ! write(*,*)'Det:',det
    a(1,1)=x(2,2)*x(3,3)-x(2,3)*x(3,2)
    a(2,1)=x(2,3)*x(3,1)-x(2,1)*x(3,3)
    a(3,1)=x(2,1)*x(3,2)-x(2,2)*x(3,1)
    a(1,2)=x(1,3)*x(3,2)-x(1,2)*x(3,3)
    a(2,2)=x(1,1)*x(3,3)-x(1,3)*x(3,1)
    a(3,2)=x(1,2)*x(3,1)-x(1,1)*x(3,2)
    a(1,3)=x(1,2)*x(2,3)-x(1,3)*x(2,2)
    a(2,3)=x(1,3)*x(2,1)-x(1,1)*x(2,3)
    a(3,3)=x(1,1)*x(2,2)-x(1,2)*x(2,1)
    a=a/det
  end subroutine inv_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xyz_to_abc(xyz,abc,lat,n)
    integer,intent(in) :: n
    real(wp), INTENT(in) :: xyz(3,n)
    real(wp), intent(in) :: lat(3,3)
    real(wp), intent(out) :: abc(3,n)

    real(wp) lat_1(3,3)
    integer i,j,k

    call inv_cell(lat,lat_1)

    abc(:,:n)=0.0d0
    do i=1,n
      do j=1,3
        do k=1,3
          abc(j,i)=abc(j,i)+lat_1(j,k)*xyz(k,i)
        end do
        abc(j,i)=dmod(abc(j,i),1.0d0)
      end do
    end do

  end subroutine xyz_to_abc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine abc_to_xyz(abc,xyz,lat,n)
    real(wp), INTENT(in) :: abc(3,*)
    real(wp), intent(in) :: lat(3,3)
    real(wp), intent(out) :: xyz(3,*)
    integer,intent(in) :: n

    integer i,j,k

    xyz(:,:n)=0.0d0
    do i=1,n
      do j=1,3
        do k=1,3
          xyz(j,i)=xyz(j,i)+lat(j,k)*abc(k,i)
        end do
      end do
    end do

  end subroutine abc_to_xyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module dftd3_core
