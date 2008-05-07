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
    set content [formatString [xml_escape_chr $content]]

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


