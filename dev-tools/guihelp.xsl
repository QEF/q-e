<?xml version="1.0" encoding="ISO-8859-1"?>

<!--
	 ***
	 *** THIS FILE IS a XSL STYLESHEET FOR TRANSFORMING INPUT_*.xml into PWgui's help file, *-help.tcl
	 ***
      -->

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">    

  <xsl:preserve-space elements="*"/>
  <xsl:output method="html"/>

  <!-- *** ROOT *** -->

  <xsl:template match="/input_description">
    # FILE AUTOMATICALLY CREATED: DO NOT EDIT, CHANGES WILL BE LOST
    <xsl:apply-templates match="namelist | card | linecard"/>
  </xsl:template>
  
  <!-- igrnore the following elements -->

  <xsl:template match="intro | toc | label | message | section |
		       subsection | subsubsection | paragraph | text">   
  </xsl:template>
  
  <!--    *** NAMELIST ***  -->

  <xsl:template match="namelist">
    <xsl:apply-templates match="descendant::var | descendant::vargroup | descendant::dimension | descendant::dimensiongroup | descendant::table"/>
  </xsl:template>


  <!--    *** CARD *** -->

  <xsl:template match="card">
    
    <!-- info for card's flags -->    
    <xsl:if test="boolean(ancestor::card/@nameless) = false()">
      help <xsl:value-of select="@name"/>_flags {
      <h2>Description of <xsl:value-of select="@name"/> card's flags</h2>
      
      <pre><xsl:value-of select="flag/info"/></pre>
      }
    </xsl:if>
    
    <xsl:apply-templates select="descendant::var | descendant::vargroup | descendant::dimension |
				 descendant::dimensiongroup | descendant::list | descendant::table"
				 mode="card_description"/>
  </xsl:template>

  <xsl:template match="flag" mode="card_description">
  </xsl:template>

  <!--    *** LINECARD *** -->

  <xsl:template match="linecard">
    <xsl:apply-templates match="descendant::var | descendant::vargroup | descendant::list"/>
  </xsl:template>
  

  <!--    *** VAR | DIMENSION | LIST ***  -->

  <xsl:template match="var | dimension | vargroup | dimensiongroup |
		       list | list/format" mode="card_description">
    <xsl:if test="info != '' or options != '' or status != '' or see != ''">
      <xsl:apply-templates select="."/>
    </xsl:if>
  </xsl:template>

  <xsl:template match="var | dimension">
    <xsl:if test="name(..) != 'vargroup' and name(..) != 'dimensiongroup'">
      help <xsl:value-of select="@name"/> {
      <ul>      
	<xsl:choose>
	  <xsl:when test="name(.)='var'">
	    <li><xsl:text>&#160;</xsl:text> <em>Variable: </em> <big><b><xsl:value-of select="@name"/></b></big></li><br/>
	  </xsl:when>
	  <xsl:when test="name(.)='dimension'">
	    <li><xsl:text>&#160;</xsl:text> <em>Variables: </em> <big><b><xsl:value-of select="@name"/>(i), i=<xsl:value-of select="@start"/>,<xsl:value-of select="@end"/></b></big></li><br/>
	  </xsl:when>
	</xsl:choose>

	<li><xsl:text>&#160;</xsl:text> <em>Type: </em> <xsl:value-of select="@type"/></li><br/>
	
	<xsl:apply-templates select="default"/> 
	<xsl:apply-templates select="status"/>
	<xsl:apply-templates select="see"/>
	<xsl:apply-templates select="info | options"/>
      </ul>      
      }
    </xsl:if>
  </xsl:template>

  <xsl:template match="list">
    help <xsl:value-of select="@name"/> {
    <ul>      
      <li><xsl:text>&#160;</xsl:text> <em>Variables: </em> <big><b><xsl:value-of select="format"/></b></big></li><br/>
      
      <li><xsl:text>&#160;</xsl:text> <em>Type: </em> <xsl:value-of select="@type"/></li><br/>
      
      <xsl:apply-templates select="default"/> 
      <xsl:apply-templates select="status"/>
      <xsl:apply-templates select="see"/>
      <xsl:apply-templates select="info | options"/>
    </ul>      
    }

    grouphelp { <xsl:value-of select="format"/> } {
    <ul>      
      <li><xsl:text>&#160;</xsl:text> <em>Variables: </em> <big><b><xsl:value-of select="format"/></b></big></li><br/>
      
      <li><xsl:text>&#160;</xsl:text> <em>Type: </em> <xsl:value-of select="@type"/></li><br/>
      
      <xsl:apply-templates select="default"/> 
      <xsl:apply-templates select="status"/>
      <xsl:apply-templates select="see"/>
      <xsl:apply-templates select="info | options"/>
    </ul>  
    }
  </xsl:template>

  <!-- *** VARGROUP | DIMENSIONGROUP *** -->

  <xsl:template match="vargroup | dimensiongroup">
    grouphelp {
    <xsl:for-each select="var | dimension">
      <xsl:value-of select="@name"/><xsl:text> </xsl:text> 
    </xsl:for-each>
    } {
    <ul>
      <xsl:if test="name(.)='vargroup'">	
	<li><xsl:text>&#160;</xsl:text> <em>Variables: </em> 
	  <big><b>
	      <xsl:for-each select="var">
		<xsl:value-of select="@name"/><xsl:if test="not(position()=last())">, </xsl:if>
	      </xsl:for-each>	
	  </b></big>
	</li><br/>
      </xsl:if>
      
      <xsl:if test="name(.)='dimensiongroup'">
	<li><xsl:text>&#160;</xsl:text> <em>Variables: </em>
	  <big><b>
	      <xsl:for-each select="dimension">
		<xsl:value-of select="@name"/>(i), 
		<xsl:if test="position()=last()"> 
		  i=<xsl:value-of select="../@start"/>,<xsl:value-of select="../@end"/>
		</xsl:if>
	      </xsl:for-each>
	  </b></big>
	</li><br/>
      </xsl:if>

      <li><xsl:text>&#160;</xsl:text> <em>Type: </em> <xsl:value-of select="@type"/></li><br/>

      <xsl:apply-templates select="default"/> 
      <xsl:apply-templates select="status"/>
      <xsl:apply-templates select="see"/>
      <xsl:apply-templates select="info | options"/>
    </ul>
    }
  </xsl:template>


  <!--    *** VAR's elements ***  -->

  <xsl:template match="default">
    <li><xsl:text>&#160;</xsl:text> <em>Default: </em> <xsl:apply-templates/></li><br/>
  </xsl:template>
  
  <xsl:template match="status">
    <li><xsl:text>&#160;</xsl:text> <em>Status: </em> <xsl:apply-templates/></li><br/>
  </xsl:template>
  
  <xsl:template match="see">
    <li><xsl:text>&#160;</xsl:text> <em>See: </em> <xsl:apply-templates/></li><br/>
  </xsl:template>
  
  <xsl:template match="info">
    <li><xsl:text>&#160;</xsl:text> <em>Description:</em></li>
    <blockquote><pre><xsl:apply-templates/></pre></blockquote>
  </xsl:template>
  
  <xsl:template match="options">
    <li><xsl:text>&#160;</xsl:text> <em>Description:</em></li>
    <blockquote>
      <xsl:apply-templates select="info | opt" mode="options"/>
    </blockquote>
  </xsl:template>

   <xsl:template match="info" mode="options">
    <pre><xsl:apply-templates/></pre>
  </xsl:template>

  <xsl:template match="opt" mode="options">
    <!--<span class="flag"><xsl:value-of select="@val"/></span><xsl:text> :</xsl:text>-->
    <dl style="margin-left: 1.5em;">
      <dt><tt><xsl:call-template name="tokenize_clist"><xsl:with-param name="clist" select="@val"/></xsl:call-template><xsl:if test="normalize-space(.) != ''"> :</xsl:if></tt></dt>
      <dd><pre style="margin-top: 0em; margin-bottom: -1em;"><xsl:apply-templates/></pre></dd>
    </dl>
  </xsl:template>

  <xsl:template name="tokenize_clist">
    <xsl:param name="clist"/>
    <xsl:variable name="first-elem" select="normalize-space(substring-before(concat($clist, ','), ','))"/> 
    <xsl:if test="$first-elem">
      <b><xsl:value-of select="$first-elem"/></b><xsl:if test="not($first-elem = normalize-space($clist))"><xsl:text>, </xsl:text></xsl:if>
      <xsl:call-template name="tokenize_clist"> 
	<xsl:with-param name="clist" select="substring-after($clist,',')" /> 
      </xsl:call-template>    
    </xsl:if>  
  </xsl:template>

  <!--    *** TABLE ***  -->

  <xsl:template match="table" mode="card_description">
    help <xsl:value-of select="@name"/> {
    <xsl:apply-templates select="rows | cols" mode="table"/>   
    }
  </xsl:template>
  
  <xsl:template match="rows | cols" mode="table">
    <xsl:apply-templates select="col | colgroup | row | rowgroup | optional | conditional" mode="table"/>
  </xsl:template>

  <xsl:template match="optional | conditional" mode="table">
    <xsl:apply-templates select="col | colgroup | row | rowgroup | optional | conditional" mode="table"/>
  </xsl:template>
  
  <xsl:template match="colgroup | rowgroup" mode="table">
    <xsl:if test="info != '' or status != '' or see != ''">
      <ul>
	<li><xsl:text>&#160;</xsl:text> <em>Variables: </em>
	  <big><b>
	      <xsl:for-each select=".//col | .//row">
		<xsl:value-of select="@name"/><xsl:if test="not(position()=last())"><xsl:text>, </xsl:text></xsl:if>
	      </xsl:for-each>
	  </b></big>
	</li><br/>

	<li><xsl:text>&#160;</xsl:text> <em>Type: </em> <xsl:value-of select="@type"/></li><br/>
	
	<xsl:apply-templates select="default"/> 
	<xsl:apply-templates select="status"/>
	<xsl:apply-templates select="see"/>
	<xsl:apply-templates select="info | options"/>
      </ul>
    </xsl:if>
  </xsl:template>
  
  <xsl:template match="col | row" mode="table">
    <xsl:if test="info != '' or status != '' or see != ''">
      <ul>
	<li><xsl:text>&#160;</xsl:text> <em>Variable: </em>
	  <big><b><xsl:value-of select="@name"/></b></big>
	</li><br/>

	<li><xsl:text>&#160;</xsl:text> <em>Type: </em> <xsl:value-of select="@type"/></li><br/>

	<xsl:apply-templates select="default"/> 
	<xsl:apply-templates select="status"/>
	<xsl:apply-templates select="see"/>
	<xsl:apply-templates select="info | options"/>
      </ul>
    </xsl:if>
  </xsl:template>

  <xsl:template match="link | ref | a">"<xsl:value-of select="."/>"</xsl:template>
  
  <!-- support a subset of HTML tags (subset of those specified in helpdoc.schema) -->
  
  <xsl:template match="b | i | tt | u | sub | sup | code | pre | br | hr | p | ol | ul | li | dl | dt | dd">
    <xsl:copy>
      <xsl:apply-templates/>
    </xsl:copy>
  </xsl:template>
  
</xsl:stylesheet>
