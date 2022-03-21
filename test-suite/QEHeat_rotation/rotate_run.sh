#!/bin/bash

ALL_CURRENTS=../../../bin/all_currents.x
./generate_all_currents.in.sh
$ALL_CURRENTS -in all_currents.in > /dev/null
$ALL_CURRENTS -in all_currents_not_rotated.in > /dev/null
python3 rotate.py <(awk '{print $3,$4,$5; print $6,$7,$8; print $9,$10,$11;}' current_hz_not_rotated.dat ) current_hz_not_rotated.dat.rot
cat current_hz_not_rotated.dat.rot
awk '{print $3,$4,$5; print $6,$7,$8; print $9,$10,$11;}' current_hz.dat
