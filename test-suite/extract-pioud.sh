# fname=benchmark.out.git.inp=diamond_pimd.in #Is it the output file wile running or reference utput

# temperature=`awk '/^[[:space:]]*[0-9]/'  'benchmark.out.git.inp=diamond_pimd.in' | awk '{print $8}' | awk 'NR > 1'`

 
# if test "$temperature" != ""; then
# 	echo temprature
# 	echo $temperature
# fi
fname=$1
# echo "$fname"
line=$(grep -E "^[[:space:]]*[2-3]+[[:space:]]+[2-3]" $fname)
block=$(echo "$line" | awk '{print $1}')
nmove=$(echo "$line" | awk '{print $2}')
H=$(echo "$line" | awk '{print $3}')
H_Nose=$(echo "$line" | awk '{print $4}')
E_pot=$(echo "$line" | awk '{print $5}')
E_kin=$(echo "$line" | awk '{print $6}')
Temp_Ha=$(echo "$line" | awk '{print $7}')
Temp_K=$(echo "$line" | awk '{print $8}')
quantum_kin_virial=$(echo "$line" | awk '{print $9}')
quantum_kin_primitive=$(echo "$line" | awk '{print $10}')

# echo "line= $line"
# echo "block = $block"
# echo "nmove = $nmove"
# echo "H = $H"
# echo "H_Nose = $H_Nose"
# echo "E_pot = $E_pot"
# echo "E_kin = $E_kin"
# echo "Temp_Ha = $Temp_Ha"
# echo "Temp_K = $Temp_K"
# echo "quantum_kin_virial = $quantum_kin_virial"
# echo "quantum_kin_primitive = $quantum_kin_primitive"


if test "$block" != ""; then
	echo block
	echo $block
fi

if test "$nmove" != ""; then
	echo nmove
	echo $nmove
fi

if test "$H" != ""; then
	echo H
	echo $H
fi

if test "$H_Nose" != ""; then
	echo H_Nose
	echo $H_Nose
fi

if test "$E_pot" != ""; then
	echo E_pot
	echo $E_pot
fi

if test "$E_kin" != ""; then
	echo E_kin
	echo $E_kin
fi

if test "$Temp_K" != ""; then
	echo Temp_K
	echo $Temp_K
fi

if test "$Temp_Ha" != ""; then
	echo Temp_Ha
	echo $Temp_Ha
fi