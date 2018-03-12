# How to fix UPF files v2.01

Many upf files written with the atomic codes contain
 `&` charactets  which are reserved in XML syntax and need to be fixed.

To fix these files you can:
 * use the `fix_upf.x` executable present in this directory ;

 *   replace with an editor the 2 or 3 occurrences of `&` with
 the corresponding escape code `&amp;`;

 * Mark  all the text contained in the  `<PP_INPUT>`  element
 as non _parsable_;  to do this replace
 ```
 <PP_INPUT>  ----> <PP_INPUT><![CDATA[
 ```

 and
 ```
 </PP_INPUT>  ---> ]]></PP_INPUT>
 ```
