 depends ../doc/manpages

 function mangen () {
 local dslash="//"
 local amp="&"
 local LINE
 while read LINE
 do
   if [[ $LINE == @* ]]
   then
     echo "if(printlist) write(iotk_output_unit,\"(a)\") &"
     echo "\"${LINE//@/}\""
     echo "printme=.false."
     echo "if(iotk_strcomp(keyword,\"all\")) printme=.true."
     for name in $LINE
     do
       [[ $name == @ ]] && continue
       echo "if(iotk_strcomp(keyword,'$name')) printme=.true."
     done
   else
     echo "if(printme) write(iotk_output_unit,\"(a)\") &"
     echo "\"${LINE//\"/\"$dslash'\"'$dslash$amp$newline\"}\""
   fi
 done < ../doc/manpages
 }

