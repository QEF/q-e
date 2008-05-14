#!/bin/sh

html=$1
tex=${1%.html}.tex

gnuhtml2latex $html

# the wiki produced html is so heavy that gnuhtml2latex fails; 
# fix the latex source

cp $tex /tmp/$tex.$$
sed 's/\\par \\\\/\\par~\\\\/g' /tmp/$tex.$$ | \
    sed 's/\\begin{description}]/\\begin{description}\\item/g' | \
    sed 's/\\begin{document}/\\begin{document}\\title{Quantum-Espresso Documentation}\\maketitle/' > $tex