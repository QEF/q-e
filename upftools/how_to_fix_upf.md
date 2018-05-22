# How to fix UPF files v2.01

Many upf files written with by previous version of  `ld1.x` contain
 `&` characters.  This  is a reserved symbol in XML syntax and needs to be 
    escaped.

To fix your UPF  files you can:
 * use the `fix_upf.x` executable present in this directory ;

 *   replace using an editor the 2 or 3 occurrences of `&` with
 the corresponding escape code `&amp;`;

 * Mark  all the text contained in the  `<PP_INPUTFILE>`  element
 as non parsable text,   to do this make these replacements:
 ```
 <PP_INPUTFILE>  ---> <PP_INPUTFILE><![CDATA[
 ```
 and
 ```
 </PP_INPUTFILE>  ---> ]]></PP_INPUTFILE>
 ```
