#!/bin/bash
touch tmp
touch tmp2
i=1
m=200
while [ $i -le $m ]; do
     cd kpt_$i
     awk '{print $1,($2+$3+$4)/3}' epsi_si_333_ZG.dat > ../epsi_si_333_ZG_av_$i.dat
     cd ../
     paste tmp epsi_si_333_ZG_av_$i.dat > aa2
     mv aa2 tmp
     i=$((i+1))
done

awk '{ 
  TR=0
  NUM = 0
  for( I = 2; I <= NF; I+=2 ) {
    TR += $I
    NUM +=1
  }
  printf "%16.8f %16.8f \n",  $1, TR/NUM
}' tmp > epsi_si_333_ZG_$m.dat

rm tmp tmp2 epsi_si_333_ZG_av_*
