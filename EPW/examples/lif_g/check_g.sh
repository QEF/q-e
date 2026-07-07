echo "====================================="
echo "Note: the phase of g-matrix elements"
echo "      between PH and EPW are different" 
echo "      due to the wave functions and "
echo "      acoustic sum rules."
echo "====================================="
m=2
n=5
nu=6
echo 'Degeneracy case (band 3/4/5):'
grep -r " ${m}        ${n}        ${nu}" *

n=2
echo 'Non-degeneracy case:'
grep -r " ${m}        ${n}        ${nu}"

nu=3
echo 'Non-degeneracy case (acoustic mode):'
grep -r " ${m}        ${n}        ${nu}"
