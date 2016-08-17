#
# XML
#
proc ::helpdoc::xml_escape_chr {content} {
    # replace xml special characters by escape-characters
    foreach {chr escChr} {
	& 	{\&amp;}  
	< 	{\&lt;}   
	> 	{\&gt;}   
    } {
	regsub -all -- $chr $content $escChr content
    }
    regsub -all -- '  $content {\&apos;} content
    regsub -all -- \" $content {\&quot;} content

    return $content
}
proc ::helpdoc::xml_attr_escape_chr {content} {
    # replace xml special characters by escape-characters
    foreach {chr escChr} {
	& 	{\&amp;}  
	< 	{\&lt;}   
	> 	{\&gt;}   
    } {
	regsub -all -- $chr $content $escChr content
    }

    return $content
}

proc ::helpdoc::xml_http {content} {
    # PURPOSE: transform all instances of http://**** into links

    #set re {http(s)*://([[:alnum:]-]+\.)+[[:alnum:]-]+[-#/\w]*}
    set re {http(s)*://([^/?#\s]*)?([^?#\s]*)(\?([^#\s]*))?(#([^\s]*))?[^\s\.,;:]}
    return [regsub -all $re $content {<link>\0</link>}]
}

proc ::helpdoc::xml_doi {content} {
    # PURPOSE: transform all instances of doi:10.**** into links
    set re {([dD][oO][iI]:?\s*)(10.[0-9]{4})(\.[0-9]+)?(/[-#%/=&@!,:\.\?\w]+[^\s\.,;:])}    
    #       1                  2            3          4

    return [regsub -all $re $content  {<a href="http://dx.doi.org/\2\3\4">\0</a>}]
}

proc ::helpdoc::xml_prb {content} {
    # PURPOSE: transform "PRB vol, page" into link
     set re {(PRB)[\s,]+([0-9]+)[\s,]+([A-Z]?[0-9]+)[\s,]+\(?[12][0-9]{3}\)?}
    return [regsub -all $re $content {<a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.\2.\3">\0</a>}]
    
}

proc ::helpdoc::xml_prl {content} {
    # PURPOSE: transform "PRL vol, page" into link
     set re {(PRL)[\s,]+([0-9]+)[\s,]+([A-Z]?[0-9]+)[\s,]+\(?[12][0-9]{3}\)?}
    return [regsub -all $re $content {<a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.\2.\3">\0</a>}]
    
}

proc ::helpdoc::xml_arxiv {content} {
    # PURPOSE: transform "arxiv:identifier" into link
    set re {(ar[xX]iv:)([0-9]{4}.[0-9]+(v[0-9]+)?|[a-zA-Z\.\-]+/[0-9]+)}
    #       1          2               3
    return [regsub -all $re $content {<a href="https://arxiv.org/abs/\2">\0</a>}]
    
}

proc ::helpdoc::xml_ref {content} {
    # PURPOSE: transform all "@ref var" into <ref>var</ref>    
    #          Note that Fortran structure names struct%var are supported
    set re {(@ref)\s+(\w+([%]\w)*)}
    return [regsub -all $re $content {<ref>\2</ref>}]
}

proc ::helpdoc::xml_link {content} {
    # PURPOSE: transform all "@link document" into <link>document</link>    

    set re {(@link)\s+([.,;:]*[\w\+-]+([.,;:][\w\+-]+)*)}
    return [regsub -all $re $content {<link>\2</link>}]
}

proc ::helpdoc::xml_tag_enter {tag attr content depth} {
    variable fid
    
    set indent [indent $depth]

    set sep ""
    if { $content != "" } {   
	if { [llength [split $content \n]] > 1 } {
	    set content [trimEmpty $content]
	    set sep \n
	} else {
	    set sep " "
	}
    }
    
    set attr    [xml_attr_escape_chr $attr]
    #set content [formatString [xml_ref [xml_link [xml_prb [xml_prl [xml_doi [xml_http [xml_escape_chr $content]]]]]]]]
    set content [formatString [xml_atTags [xml_arxiv [xml_prb [xml_prl [xml_doi [xml_http [xml_escape_chr $content]]]]]]]]

    if { $attr != "" } {
	puts $fid(xml) "${indent}<$tag ${attr}>${sep}${content}"
    } else {
	puts $fid(xml) "${indent}<$tag>${sep}${content}"
    }
}

proc ::helpdoc::xml_tag_leave {tag attr content depth} {
    variable fid
    puts $fid(xml) "[indent $depth]</$tag>"
}


