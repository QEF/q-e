#!/bin/sh

Convert() {
    geom=$1
    input=$2
    output=${input%.eps}.png
    convert -antialias -geometry $1 $input $output
    #display $output
}

Convert 50% Guib.eps      
Convert 60% GUI-open.eps
Convert 50% GUI-new.eps
Convert 60% myGUI.eps
Convert 50% TkPWscf.eps
Convert 60% parsing.eps    
