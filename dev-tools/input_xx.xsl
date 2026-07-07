<?xml version="1.0" encoding="ISO-8859-1"?>

<!--
***
*** THIS FILE IS a XSL STYLESHEET FOR TRANSFORMING INPUT_*.xml to INPUT_*.html
***
-->

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <!--<xsl:strip-space elements="*"/>-->
  <xsl:output method="html"/>  

  <!-- params that are passed on command line -->
  <xsl:param name="version"/>  
  <xsl:param name="current-date"/>
  
  <!-- *** ROOT *** -->

  <xsl:template match="/input_description">
    <html>
      <head>
	<xsl:comment> *** FILE AUTOMATICALLY CREATED: DO NOT EDIT, CHANGES WILL BE LOST *** </xsl:comment>
	<meta http-equiv="Content-Style-Type" CONTENT="text/css" />
	<style>
	  body {
	  background-color:#ffffff;
	  font:normal 14px/1.8em arial, helvetica, sans-serif;
	  width:900px;
	  text-align:justify;
	  margin: 30 10 10 30;
	  }

	  h1 {
	  font-size:24px;
	  }
	  
	  h2 {
	  font-size:18px;
	  }
	  
	  h3 {
	  font-size:16px;
	  }
	  pre, tt, code {
	  font-size:14px;
	  }	  
	  .syntax, .syntax table {
	  font-size:14px;
	  }
	  span.namelist {
	  color: #214478;
	  }
	  span.card {
	  color: #782167;
	  }
	  span.flag {
	  color: #008000;
	  font-weight: bold;
	  }
	</style>
	<title><xsl:value-of select="@program"/>: input description</title>
      </head>
      <body>
	<a name="__top__"></a>
	<table style="border-width: 0; table-layout: auto; width: 100%; text-align: left; vertical-align: top; background: #00395a;">
	  <tr>
	    <th style="margin: 3 3 3 10; background: #005789; background: linear-gradient(rgba(0,87,137,1),rgba(0,119,189,1)); color: #ffffee; ">
	      <h1 style="margin: 10 10 10 15; text-align: left;"> Input File Description </h1>
	      <h2 style="margin: 10 10 10 15; text-align: left;"> Program:
	      <xsl:value-of select="@program"/> / <xsl:value-of select="@package"/> / <xsl:value-of select="@distribution"/>
	      <xsl:if test="$version">
		<span style="font-weight: normal;"> (version: <xsl:value-of select="string($version)"/>)</span>
	      </xsl:if>
	      </h2>
	    </th>
	  </tr>
	  <tr><td style="padding: 10 3 3 3; background: #ffffff; color: #222222; ">
	    <xsl:apply-templates/>
	  </td></tr>
	</table>
	<small>
	  This file has been created by helpdoc utility<xsl:if test="$version"> on <xsl:value-of select="$current-date" /></xsl:if>.
	</small>
      </body>
    </html>
  </xsl:template>

  
  <!--  *** TOC ***  -->

  <xsl:template match="toc">
    <blockquote style="margin-bottom: 2em;">
      <h3>TABLE OF CONTENTS</h3>
      <blockquote>	
	<xsl:apply-templates select=".." mode="toc"/>
      </blockquote>
    </blockquote>    
  </xsl:template>

  
  <!-- TOC templates, i.e. mode = "toc" -->

  <xsl:template match="intro" mode="toc">
    <p><a href="#{generate-id(.)}">INTRODUCTION</a></p>
  </xsl:template>

  <xsl:template match="supercard" mode="toc">
    <p><a href="#{generate-id(.)}"><xsl:call-template name="supercard_name"/></a></p>    
    <blockquote>
      <xsl:apply-templates select="./supercard | ./namelist | ./card | ./linecard |
				   ./optional | ./conditional | ./group | ./if | ./choose " mode="toc"/>
    </blockquote>    
  </xsl:template>

  <xsl:template match="optional | conditional | group | if" mode="toc">
    <xsl:apply-templates select="./supercard | ./namelist | ./card | ./linecard |
				 ./optional | ./conditional | ./group | ./if | ./choose " mode="toc"/>
  </xsl:template>
  <xsl:template match="choose" mode="toc">
    <xsl:for-each select="./when | ./elsewhen | ./ otherwise">
      <xsl:apply-templates select="./supercard | ./namelist | ./card | ./linecard |
				   ./optional | ./conditional | ./group | ./if | ./choose " mode="toc"/>
    </xsl:for-each>
  </xsl:template>
  
  <xsl:template match="linecard" mode="toc">
    <p><a href="#{generate-id(.)}">Line-of-input:</a><xsl:text> </xsl:text>
	    <xsl:apply-templates select=".//var | .//dimension | .//multidimension | .//list" mode="toc"/></p>
  </xsl:template>

  <xsl:template match="namelist | card" mode="toc">
    <p><a href="#{generate-id(.)}">
      <xsl:if test="name(.)='namelist'">&#38;</xsl:if>
      <xsl:value-of select="@name"/>
    </a></p>
    <xsl:if test=".//var != '' or
		  .//dimension != '' or
		  .//multidimension != '' or
		  .//list != '' or
		  .//col != '' or
		  .//row != ''">
      <blockquote>
	      <xsl:apply-templates select=".//var | .//dimension | .//multidimension | .//list | .//col | .//row" mode="toc"/>
      </blockquote>
    </xsl:if>
  </xsl:template>

  <xsl:template match="var | dimension | multidimension" mode="toc">
    <xsl:if test="info != '' or 
		  status != '' or 
		  see    != '' or
		  options != '' or
		  ../../vargroup/info != '' or 
		  ../../dimensiongroup/info != '' or
		  ../../dimensiongroup/options != '' or 
		  ../../multidimensiongroup/info != '' or 
		  ../../multidimensiongroup/options != ''">
      <a href="#{generate-id(.)}"><xsl:value-of select="@name"/></a> 
      <xsl:if test="not(position()=last())">
	<xsl:text> | </xsl:text>
      </xsl:if>
    </xsl:if>
  </xsl:template>

  <xsl:template match="list" mode="toc">
    <xsl:if test="info != '' or status != '' or see != '' or options != ''">
      <a href="#{generate-id(.)}"><xsl:value-of select="format"/></a> 
      <xsl:if test="not(position()=last())">
	<xsl:text> | </xsl:text>
      </xsl:if>
    </xsl:if>
  </xsl:template>

  <xsl:template match="col" mode="toc">
    <xsl:if test="info != '' or status != '' or see != '' or options != ''">
      <a href="#{generate-id(.)}"><xsl:value-of select="@name"/></a>
      <xsl:if test="not(position()=last())">
	<xsl:text> | </xsl:text>
      </xsl:if>
    </xsl:if>
    <xsl:if test="ancestor::colgroup/info != '' or ancestor::colgroup/options != ''">
      <a href="#{generate-id(.)}"><xsl:value-of select="@name"/></a>
      <xsl:if test="not(position()=last())">
	<xsl:text> | </xsl:text>
      </xsl:if>
    </xsl:if>
  </xsl:template>
  
  <xsl:template match="row" mode="toc">
    <xsl:if test="info != '' or status != '' or see != '' or options != ''">
      <a href="#{generate-id(.)}"><xsl:value-of select="@name"/></a>
      <xsl:if test="not(position()=last())">
	<xsl:text> | </xsl:text>
      </xsl:if>
    </xsl:if>
    <xsl:if test="ancestor::rowgroup/info != '' or ancestor::colgroup/options != ''">
      <a href="#{generate-id(.)}"><xsl:value-of select="@name"/></a>
      <xsl:if test="not(position()=last())">
	<xsl:text> | </xsl:text>
      </xsl:if>
    </xsl:if>
  </xsl:template>

  <!-- aaaaa -->
  <xsl:template match="section" mode="toc">
    <p><a href="#{generate-id(.)}"><xsl:value-of select="@title"/></a></p>
      <xsl:apply-templates select="subsection" mode="toc"/>
  </xsl:template>
  
  <!--new: END-->
  
  <xsl:template match="subsection" mode="toc">
    <blockquote>
      <a href="#{generate-id(.)}"><xsl:value-of select="@title"/></a>
      <xsl:apply-templates select="subsubsection" mode="toc"/>
    </blockquote>
  </xsl:template>
  <xsl:template match="subsubsection" mode="toc">
    <blockquote>
      <a href="#{generate-id(.)}"><xsl:value-of select="@title"/></a>
    </blockquote>
  </xsl:template>
  <!-- END of TOC templates -->
  
  <xsl:template name="supercard_name">
    <xsl:choose>
      <xsl:when test="@starttag != ''">
	<xsl:value-of select="@starttag"/>
      </xsl:when>
      <xsl:otherwise>
	<xsl:value-of select="@name"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!-- back-to-top -->

  <xsl:template name="back_to_top">
    <div align="right" style="margin-bottom: 5;">[<a href="#__top__">Back to Top</a>]</div>
  </xsl:template>

  
  <!--    *** INTRO ***  -->

  <xsl:template match="intro">
    <blockquote>
      <a name="{generate-id()}"></a>
      <h3>INTRODUCTION</h3>      
      <blockquote>
	<pre>
	  <xsl:apply-templates/>
	</pre>
      </blockquote>
    </blockquote>
  </xsl:template>


  <!--    *** GROUP ***  -->

  <xsl:template match="group">
    <table style="border-color: #bb9977; border-style: solid; border-width: 3; margin-bottom: 10; table-layout: auto; background-color: #FFddbb; width: 100%; padding: 5 5 0 30">
      <tr><td>
	<xsl:apply-templates/>
      </td></tr>
    </table>
  </xsl:template>

  
  <!--    *** NAMELIST ***  -->

  <xsl:template match="namelist">
    <a name="{generate-id(.)}"></a>
    <a name="{@name}"></a>
    <table border="0" width="100%" style="margin-bottom: 20;">
      <tr>
	<th bgcolor="#ddcba6">
	  <h2 style="margin: 10 10 10 15; text-align: left;"> Namelist: <span class="namelist"><span style="font-weight:normal">&#38;</span><xsl:value-of select="@name"/></span> </h2>
	</th>
      </tr>
      <tr><td style="text-align: left; background: #ffebc6; padding: 5 5 5 30; ">	    
	<table style="border-color: #505087; border-style: solid; border-width: 0; margin-bottom: 10; table-layout: auto; width: 800;">
	  <tbody>
	    <tr><td>
	      <xsl:apply-templates/>
	    </td></tr>
	  </tbody>
	</table>
      </td></tr>
    </table>
  </xsl:template>

  <!--    *** SUPERCARD *** -->
  
  <xsl:template match="supercard">
    <a name="{generate-id(.)}"></a>
    <a name="{@name}"></a>    
    <table style="table-layout: auto; width: 100%; border: 3px solid #1b587b; border-collapse: collapse; margin: 10 5 20 5; padding-right: 5px;">
      <tr>
	<th bgcolor="#c8c4b7">
	  <h2 style="margin: 10 10 10 15; text-align: left; white-space: nowrap;"> 
	    <xsl:call-template name="supercard_name"/>
	  </h2>
	</th>
      </tr>
      <tr>
	<td bgcolor="#eeeeee" style="padding: 5 10 5 15;">
	  <i>Syntax of this supercard is the following:</i>
	  <br/>
	  <pre><xsl:call-template name="supercard_name"/><br/><i>&#160;&#160;... content of the supercard here ...</i><br/><xsl:value-of select="@endtag"/></pre>
	  <i>and the content is:</i>
	</td>
      </tr>
      <xsl:if test="@remark != ''">
	<tr>
	  <td style="padding: 10 10 10 15; background: #ffffff; text-align: left;">
	    <i>( <b>Remark:</b> <xsl:value-of select="@remark"/> )</i>
	  </td>
	</tr>
      </xsl:if>
      
      <tr><td style="text-align: left; background: #ffffff; padding: 5 5 5 30; ">	    
	<xsl:apply-templates/>
      </td></tr>
      <xsl:if test="@endtag != ''">
	<tr>
	  <th bgcolor="#c8c4b7">
	    <h2 style="margin: 10 10 10 15; text-align: left; white-space: nowrap;"> 
	      <xsl:value-of select="@endtag"/>
	    </h2>
	  </th>
	</tr>
      </xsl:if>
    </table>
  </xsl:template>
  
  <!--    *** CARD *** -->

  <xsl:template match="card">
    <a name="{generate-id(.)}"></a>
    <a name="{@name}"></a>
    <table border="0" style="margin-bottom: 20; table-layout: auto; width: 100%;">
      <tr>
	<th bgcolor="#ddcba6">
	  <h2 style="margin: 10 10 10 15; text-align: left; white-space: nowrap;"> 
	    Card: <span class="card"><xsl:value-of select="@name"/></span>
	    <xsl:choose>
	      <xsl:when test="flag/@use = 'optional'">
		<xsl:text> { </xsl:text>
		<xsl:call-template name="tokenize_enum"> 
		  <xsl:with-param name="enums" select="flag/enum"/> 
		</xsl:call-template>
		<!--<span class="flag"><xsl:value-of select="flag/enum"/></span>-->
		<xsl:text> } </xsl:text>
	      </xsl:when>
	      <xsl:otherwise>
		<xsl:text> </xsl:text>
		<xsl:call-template name="tokenize_enum"> 
		  <xsl:with-param name="enums" select="flag/enum"/> 
		</xsl:call-template>
		<!--<span class="flag"><xsl:value-of select="flag/enum"/></span>-->
	      </xsl:otherwise>
	    </xsl:choose>
	  </h2>
	</th>
      </tr>
      
      <tr><td style="text-align: left; background: #ffebc6; padding: 5 5 5 30; ">	    
	<table style="border-color: #505087; border-style: solid; border-width: 0; margin-bottom: 10; table-layout: auto; width: 100%;">
	  <tbody>	      	      
	    <tr><td>
	      <xsl:apply-templates select="syntax | if | choose | label | message"/>
	    </td></tr>
	    <xsl:if test=".//info  != '' or
			  ./status != '' or
			  .//see   != '' or
			  .//options != ''">
	      <tr><td>
		<h3>Description of items:</h3>
		<blockquote>
		  <!-- here: apply templates ... -->
		  <xsl:apply-templates select="flag | descendant::vargroup | descendant::var | descendant::list | descendant::table" mode="card_description"/>
		</blockquote>
	      </td></tr>
	    </xsl:if>
	  </tbody>
	</table>
      </td></tr>
    </table>
  </xsl:template>

  <xsl:template name="tokenize_enum">
    <xsl:param name="enums"/>    
    <xsl:variable name="first-enum" select="normalize-space(substring-before(concat($enums, '|'), '|'))"/> 
    <xsl:if test="$first-enum">
      <span class="flag"><xsl:value-of select="$first-enum"/></span>
      <xsl:if test="not($first-enum = normalize-space($enums))"><xsl:text> | </xsl:text></xsl:if>
      <xsl:call-template name="tokenize_enum"> 
	<xsl:with-param name="enums" select="substring-after($enums,'|')" /> 
      </xsl:call-template>    
    </xsl:if>  
  </xsl:template>
    
  <!-- card/syntax -->

  <xsl:template match="syntax">
    <h3>Syntax:</h3>
    <blockquote>
      <xsl:if test="boolean(ancestor::card/@nameless) = false()">
	<b style="white-space: nowrap;">
	  <xsl:value-of select="ancestor::card/@name"/>
	  <xsl:choose>
	    <xsl:when test="normalize-space(@flag) = ''">
	      <xsl:choose>
		<xsl:when test="ancestor::card/flag/@use = 'optional'">
		  <xsl:text> { </xsl:text> 
		  <xsl:value-of select="ancestor::card/flag/enum" /> 
		  <xsl:text> } </xsl:text>
		</xsl:when>
		<xsl:when test="ancestor::card/flag/@use = 'conditional'">
		  <xsl:text> [ </xsl:text> 
		  <xsl:value-of select="ancestor::card/flag/enum" /> 
		  <xsl:text> ] </xsl:text>
		</xsl:when>
		<xsl:otherwise>
		  <xsl:text> </xsl:text><xsl:value-of select="ancestor::card/flag/enum" />
		</xsl:otherwise>
	      </xsl:choose>
	    </xsl:when>
	    <xsl:otherwise>
	      <xsl:text> </xsl:text><xsl:value-of select="@flag" /> 
	    </xsl:otherwise>
	  </xsl:choose>
	</b>
	<br/>
      </xsl:if>
      <div class="syntax">
	<xsl:apply-templates select="table | line | optional | conditional | list" mode="syntax"/>
      </div>
    </blockquote>
  </xsl:template>

  <!-- card//syntax//line -->

  <xsl:template match="line" mode="syntax">
    <xsl:apply-templates select="optional | conditional | var | keyword | constant | etc | vargroup | list" mode="syntax"/>
    <br/>
  </xsl:template>	


  <!-- card//syntax//optional -->      

  <xsl:template match="optional" mode="syntax">
    <!--<div style="background: #eeeeee; color: #555555;">-->
    <xsl:text> { </xsl:text>
    <xsl:apply-templates select="line | var | vargroup | list | keyword | constant | etc | table" mode="syntax"/>
    <xsl:text> } </xsl:text>	    
    <!--</div>-->
  </xsl:template>	      

  <!-- card//syntax//conditional -->      

  <xsl:template match="conditional" mode="syntax">
    <xsl:text> [ </xsl:text>
    <xsl:apply-templates select="line | var | vargroup | list | keyword | constant | etc | table" mode="syntax"/>
    <xsl:text> ] </xsl:text>	    
  </xsl:template>	      

  <!-- card//syntax//keyword -->      

  <xsl:template match="keyword" mode="syntax">
    <b><xsl:value-of select="@name"/></b><xsl:text>&#160;&#160;</xsl:text>
  </xsl:template>

  <!-- card//syntax//constant -->

  <xsl:template match="constant" mode="syntax">
    <xsl:apply-templates/><xsl:text>&#160;&#160;</xsl:text>
  </xsl:template>

  <!-- card//syntax//etc --> 

  <xsl:template match="etc" mode="syntax">
    <xsl:text>. . .&#160;&#160;</xsl:text>
  </xsl:template>

  <!-- card//syntax//list -->

  <xsl:template match="list" mode="syntax">    
    <i>
      <xsl:choose>	
	<xsl:when test="info != '' or options != ''">
	  <a href="#{generate-id(.)}"><xsl:value-of select="format"/></a>
	</xsl:when>
	<xsl:otherwise>
	  <xsl:value-of select="format"/>
	</xsl:otherwise>
      </xsl:choose>
    </i>
    <xsl:text>&#160;&#160;</xsl:text>
  </xsl:template>

  <!-- card//syntax//var -->      

  <xsl:template match="var" mode="syntax">
    <i>
      <xsl:choose>	
	<xsl:when test="info != '' or options != '' or status != '' or see != '' or
			../../vargroup/info != '' or ../../vargroup/options != ''">
	  <a href="#{generate-id(.)}"><xsl:value-of select="@name"/></a>
	</xsl:when>
	<xsl:otherwise>
	  <xsl:value-of select="@name"/>
	</xsl:otherwise>
      </xsl:choose>
    </i>
    <xsl:text>&#160;&#160;</xsl:text>
  </xsl:template>

  <!-- card//syntax//vargroup -->      

  <xsl:template match="vargroup" mode="syntax">
    <xsl:apply-templates select="var | constant | keyword" mode="syntax"/>
  </xsl:template>

  <!-- card//syntax//table -->      

  <xsl:template match="table" mode="syntax">
    <a name="{generate-id(.)}"></a>
    <table>
      <xsl:apply-templates select="rows | cols" mode="syntaxTableMode"/>
    </table>    
  </xsl:template>

  <!-- syntax//table/rows -->

  <xsl:template match="rows" mode="syntaxTableMode">
    <tr>
      <xsl:call-template name="row">
	<xsl:with-param name="rowID"><xsl:value-of select="@start"/></xsl:with-param>
      </xsl:call-template>
    </tr> 
    <tr>
      <xsl:call-template name="row">
	<xsl:with-param name="rowID"><xsl:value-of select="number(@start+1)"/></xsl:with-param>
      </xsl:call-template>
    </tr> 
    <xsl:choose>
      <xsl:when test="number(@end) != @end">
	<tr><td colspan="2"><xsl:text>&#160;. . .</xsl:text></td></tr>
	<tr>
	  <xsl:call-template name="row">
	    <xsl:with-param name="rowID"><xsl:value-of select="@end"/></xsl:with-param>
	  </xsl:call-template>
	</tr> 
      </xsl:when>
      <xsl:otherwise>
	<xsl:if test="number(@end) > number(@start+2)">
	  <tr><td colspan="2"><xsl:text>&#160;. . .</xsl:text></td></tr>
	</xsl:if>
	<xsl:if test="number(@end) > number(@start+1)">
	  <tr>
	    <xsl:call-template name="row">
	      <xsl:with-param name="rowID"><xsl:value-of select="@end"/></xsl:with-param>
	    </xsl:call-template>
	  </tr>
	</xsl:if>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="row">
    <xsl:param name="rowID" select="1"/>   
    <xsl:apply-templates select="col | colgroup | optional | conditional" mode="rowMode">
      <xsl:with-param name="rowID" select="$rowID"/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="colgroup" mode="rowMode">
    <xsl:param name="rowID"/>
    <xsl:apply-templates select="col | optional | conditional" mode="rowMode">
      <xsl:with-param name="rowID" select="$rowID"/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="optional" mode="rowMode">
    <xsl:param name="rowID"/>
    <td><xsl:text> { </xsl:text></td>
    <xsl:apply-templates select="col | colgroup | conditional | optional" mode="rowMode">
      <xsl:with-param name="rowID" select="$rowID"/>
    </xsl:apply-templates>
    <td><xsl:text> } </xsl:text></td>
  </xsl:template>

  <xsl:template match="conditional" mode="rowMode">
    <xsl:param name="rowID"/>
    <td><xsl:text> [ </xsl:text></td>
    <xsl:apply-templates select="col | colgroup | conditional | optional" mode="rowMode">
      <xsl:with-param name="rowID" select="$rowID"/>
    </xsl:apply-templates>
    <td><xsl:text> ] </xsl:text></td>
  </xsl:template>

  <xsl:template match="col" mode="rowMode">
    <xsl:param name="rowID"/>
    <td style="white-space:nowrap">  
      <xsl:text>&#160;</xsl:text>
      <i>
	<xsl:choose>	
	  <xsl:when test="info != '' or options != '' or status != '' or see != ''
			  or ../../colgroup/info != '' or ../../colgroup/options != ''">
	    <a href="#{generate-id(.)}"><xsl:value-of select="@name"/>(<xsl:value-of select="$rowID"/>)</a>
	  </xsl:when>
	  <xsl:otherwise>
	    <xsl:value-of select="@name"/>(<xsl:value-of select="$rowID"/>)
	  </xsl:otherwise>
	</xsl:choose>
      </i>
      <xsl:text>&#160;</xsl:text>
    </td>  
  </xsl:template>
  
  <!-- syntax//table/cols -->

  <xsl:template match="cols" mode="syntaxTableMode">
    <xsl:apply-templates select="row | rowgroup | optional | conditional" mode="colsMode">
      <xsl:with-param name="colsOptional"  select="false()"/>
      <xsl:with-param name="colsConditional" select="false()"/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="row" mode="colsMode">
    <xsl:param name="colsOptional" select="false()"/>
    <xsl:param name="colsConditional" select="false()"/>
    <tr>
      <td align="right" style="white-space:nowrap">
	<xsl:if test="$colsOptional    = true()"><xsl:text>{ &#160;</xsl:text></xsl:if>
        <xsl:if test="$colsConditional = true()"><xsl:text>[ &#160;</xsl:text></xsl:if>
      </td>
      <xsl:call-template name="insertColumns"/>
      <td align="left" style="white-space:nowrap">
	<xsl:if test="$colsConditional = true()"><xsl:text>&#160; ]</xsl:text></xsl:if>
	<xsl:if test="$colsOptional    = true()"><xsl:text>&#160; }</xsl:text></xsl:if>
      </td>
    </tr>
  </xsl:template>

  <xsl:template match="rowgroup" mode="colsMode">
    <xsl:param name="colsOptional"/>
    <xsl:param name="colsConditional"/>
    <xsl:apply-templates select="row | optional | conditional" mode="colsMode">
      <xsl:with-param name="colsOptional" select="$colsOptional"/>
      <xsl:with-param name="colsConditional" select="$colsConditional"/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="optional" mode="colsMode">
    <xsl:param name="colsOptional"/>
    <xsl:param name="colsConditional"/>
    <xsl:apply-templates select="row | rowgroup | conditional" mode="colsMode">
      <xsl:with-param name="colsOptional" select="true()"/>
      <xsl:with-param name="colsConditional" select="$colsConditional"/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="conditional" mode="colsMode">
    <xsl:param name="colsOptional"/>
    <xsl:param name="colsConditional"/>
    <xsl:apply-templates select="row | rowgroup | optional" mode="colsMode">
      <xsl:with-param name="colsOptional" select="$colsOptional"/>
      <xsl:with-param name="colsConditional" select="true()"/>
    </xsl:apply-templates>
  </xsl:template>
  
  <xsl:template name="insertColumns">
    <xsl:call-template name="insertCol">
      <xsl:with-param name="colID" select="ancestor::cols/@start"/>
    </xsl:call-template>
    <xsl:call-template name="insertCol">
      <xsl:with-param name="colID" select="number(ancestor::cols/@start+1)"/>
    </xsl:call-template>
    <xsl:choose>
      <xsl:when test="number(ancestor::cols/@end) != ancestor::cols/@end">
	<td><xsl:text>&#160;. . .</xsl:text></td>
	<xsl:call-template name="insertCol">
	  <xsl:with-param name="colID" select="ancestor::cols/@end"/>
	</xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
	<xsl:if test="number(ancestor::cols/@end) > number(ancestor::cols/@start+2)">
	  <td><xsl:text>&#160;. . .</xsl:text></td>
	</xsl:if>
	<xsl:if test="number(ancestor::cols/@end) > number(ancestor::cols/@start+1)">
	  <xsl:call-template name="insertCol">
	    <xsl:with-param name="colID" select="ancestor::cols/@end"/>
	  </xsl:call-template>
	</xsl:if>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="insertCol">
    <xsl:param name="colID"/>
    <td>  
      <xsl:text>&#160;</xsl:text>
      <i>
	<xsl:choose>	
	  <xsl:when test="info != '' or options != '' or status != '' or see != ''
			  or ../../rowgroup/info != '' or ../../rowgroup/options != ''">
	    <a href="#{generate-id(.)}"><xsl:value-of select="@name"/>(<xsl:value-of select="$colID"/>)</a>
	  </xsl:when>
	  <xsl:otherwise>
	    <xsl:value-of select="@name"/>(<xsl:value-of select="$colID"/>)
	  </xsl:otherwise>
	</xsl:choose>
      </i>
      <xsl:text>&#160;</xsl:text>
    </td>    
  </xsl:template>
  
  
  <!--    *** LINECARD *** -->

  <xsl:template match="linecard">
    <a name="{generate-id(.)}"></a>
    <table border="0" width="100%" style="margin-bottom: 20; ">
      <tr>
	<th bgcolor="#ddcba6">
	  <h3 style="margin: 10 10 10 15; text-align: left;"> 
	    Line of input 
	  </h3>
	</th>
      </tr>
      
      <tr><td style="text-align: left; background: #ffebc6; padding: 5 5 5 30; ">	    
	<table style="border-color: #505087; border-style: solid; border-width: 0; margin-bottom: 10; table-layout: auto; width: 100%; ">
	  <tbody>	      	      
	    <tr><td>
	      <h3>Syntax:</h3>
	      <blockquote>
		<xsl:apply-templates select="keyword | var | vargroup | list | optional | conditional" mode="syntax"/>
	      </blockquote>
	    </td></tr>
	    <tr><td>
	      <h3>Description of items:</h3>
	      <blockquote>
		<xsl:apply-templates select="descendant::vargroup | descendant::var | descendant::list"/>
	      </blockquote>
	    </td></tr>
	  </tbody>
	</table>
      </td></tr>
    </table>
  </xsl:template>
  
  <!--    *** LABEL ***  -->

  <xsl:template match="label">
    <p><b><xsl:apply-templates/></b></p>
  </xsl:template>


  <!--    *** MESSAGE ***  -->

  <xsl:template match="message">
    <p><pre>
      <xsl:apply-templates/>
    </pre></p>
  </xsl:template>


  <!--    *** IF ***  -->

  <xsl:template match="if">
    <table style="border-color: #bb9977; border-style: solid; border-width: 3; margin-bottom: 10; table-layout: auto; background-color: #FFddbb; width: 100%; padding: 5 5 0 5">
      <tr><td>
	<b>IF </b>  <tt><em><xsl:value-of select="@test"/></em> :</tt> 
	<blockquote>
	  <xsl:apply-templates/>
	</blockquote>
      </td></tr>
    </table>
  </xsl:template>


  
  <!--    *** CHOOSE ... ***  -->

  <xsl:template match="choose">
    <table style="border-color: #bb9977; border-style: solid; border-width: 3; margin-bottom: 10; table-layout: auto; width: 100%; padding: 5 5 0 5">
      <tr><td>
	<xsl:apply-templates select="when"/>	
	<xsl:apply-templates select="elsewhen"/>	
	<xsl:apply-templates select="otherwise"/>
      </td></tr>
      <xsl:apply-templates select="message | label" mode="choose"/>
    </table>
  </xsl:template>

  <xsl:template match="when">  
    <b>IF </b> <tt><em><xsl:value-of select="@test"/></em> :</tt>
    <blockquote>
      <table style="border-color: #bb9977; border-style: solid; border-width: 3; margin-bottom: 10; table-layout: auto; background-color: #FFddbb; width: 100%; padding: 5 5 0 30">
	<tr><td>
	  <xsl:apply-templates/>
	</td></tr>
      </table>
    </blockquote>
  </xsl:template>

  <xsl:template match="elsewhen">  
    <b>ELSEIF </b> <tt><em><xsl:value-of select="@test"/></em> :</tt>
    <blockquote>
      <table style="border-color: #bb9977; border-style: solid; border-width: 3; margin-bottom: 10; table-layout: auto; background-color: #FFddbb; width: 100%; padding: 5 5 0 30">
	<tr><td>
	  <xsl:apply-templates/>	  
	</td></tr>
      </table>
    </blockquote>
  </xsl:template>

  <xsl:template match="otherwise">  
    <b>ELSE </b>
    <blockquote>
      <table style="border-color: #bb9977; border-style: solid; border-width: 3; margin-bottom: 10; table-layout: auto; background-color: #FFddbb; width: 100%; padding: 5 5 0 30">
	<tr><td>	
	  <xsl:apply-templates/>
	</td></tr>
      </table>
    </blockquote>
  </xsl:template>
  
  <!-- *** VARGROUP | DIMENSIONGROUP | MULTIDIMENSIONGROUP *** -->

  <xsl:template match="vargroup | dimensiongroup | multidimensiongroup" mode="card_description">
    <!--<xsl:if test="child::node() != ''">-->
    <xsl:if test="info != '' or options != '' or status != '' or see != ''">
      <xsl:apply-templates select="."/>
    </xsl:if>
  </xsl:template>

  <!--    *** VAR | DIMENSION | LIST | FLAG ***  -->

  <xsl:template match="var | list | dimension | multidimension | flag" mode="card_description">
    <!--<xsl:if test="child::node() != ''">-->
    <xsl:if test="info != '' or options != '' or status != '' or see != ''">
      <xsl:apply-templates select="."/>
    </xsl:if>
  </xsl:template>

  <xsl:template match="var | list | dimension | multidimension | flag">
    <xsl:if test="name(..) != 'vargroup' and name(..) != 'dimensiongroup' and name(..) != 'multidimensiongroup' ">
      <a name="{generate-id(.)}"></a>
      <a name="{substring-before(concat(@name,'('),'(')}"></a>  <!-- to take care of cases if varname is specified as
                                                                     dimension, i.e., varname(i,j,k) -->
      <table width="100%" style="border-color:   #b5b500; border-style: solid; border-width: 2; margin-bottom: 10; table-layout: auto; background-color: #FFFFFF;">
	<tr>
	  <xsl:choose>
	    <xsl:when test="name(.)='var'">
	      <th align="left" valign="top" width="20%" style="background: #ffff99; padding: 2 2 2 10; ">
		<xsl:value-of select="@name"/>
	      </th>
	    </xsl:when>
	    <xsl:when test="name(.)='dimension'">
	      <th width="20%" style="white-space: nowrap; text-align: left; vertical-align: top; background: #ffff99; padding: 2 2 2 10; ">
		<xsl:value-of select="@name"/>(i), i=<xsl:value-of select="@start"/>,<xsl:value-of select="@end"/>
	      </th>
	    </xsl:when>
	    <xsl:when test="name(.)='multidimension'">
		    <th width="20%" style="white-space: nowrap; text-align: left; vertical-align: top; background: #ffff99; padding: 2 2 2 10; ">
	            <xsl:value-of select="@name"/>(<xsl:value-of select="@indexes"/>), (<xsl:value-of select="@indexes"/>) = (<xsl:value-of select="@start"/>) . . . (<xsl:value-of select="@end"/>)
                    </th>
            </xsl:when> 
	    <xsl:when test="name(.)='flag'">
	      <th width="20%" style="white-space: nowrap; text-align: left; vertical-align: top; background: #ffff99; padding: 2 2 2 10; ">
		<i>Card's options:</i>
	      </th>
	    </xsl:when>
	    <xsl:otherwise>
	      <th width="20%" style="white-space: nowrap; text-align: left; vertical-align: top; background: #ffff99; padding: 2 2 2 10; ">
		<xsl:value-of select="format"/>
	      </th>
	    </xsl:otherwise>
	  </xsl:choose>
	  <xsl:choose>
	  <xsl:when test="name(.) != 'flag'">
	    <td style="text-align: left; vertical-align: top; background: #ffffc3; padding: 2 2 2 5; ">
	      <xsl:value-of select="@type"/>
	    </td>
	  </xsl:when>
	  <xsl:otherwise>
	    <th style="text-align: left; vertical-align: top; background: #ffffc3; padding: 2 2 2 5; ">
	      <xsl:call-template name="tokenize_enum"> 
		<xsl:with-param name="enums" select="enum"/> 
	      </xsl:call-template>
	    </th>
	  </xsl:otherwise>
	  </xsl:choose>
	</tr>
	<xsl:apply-templates select="default"/> 
	<xsl:apply-templates select="dimensionality"/>
	<xsl:apply-templates select="units"/>
	<xsl:apply-templates select="status"/>
	<xsl:apply-templates select="see"/>
	<xsl:apply-templates select="info | options"/>
      </table>
      <xsl:call-template name="back_to_top"/>
    </xsl:if>
  </xsl:template>

  <!-- *** VARGROUP | DIMENSIONGROUP *** -->

  <xsl:template match="vargroup | dimensiongroup | multidimensiongroup">
    <table width="100%" style="border-color:   #b5b500; border-style: solid; border-width: 2; margin-bottom: 10; table-layout: auto; background-color: #FFFFFF;">
      <tr>
	<th align="left" valign="top" width="20%" style="white-space: nowrap; background: #ffff99; padding: 2 2 2 10; ">
	  <xsl:if test="name(.)='vargroup'">
	    <xsl:for-each select="var">
	      <a name="{generate-id(.)}"></a>
	      <a name="{substring-before(concat(@name,'('),'(')}"></a>  <!-- to take care of cases if varname is specified as
                                                                             dimension, i.e., varname(i,j,k) -->
	      <xsl:value-of select="@name"/><xsl:if test="not(position()=last())">, </xsl:if>
	    </xsl:for-each>	
	  </xsl:if>
	  <xsl:if test="name(.)='dimensiongroup'">
	    <xsl:for-each select="dimension">
	      <a name="{generate-id(.)}"></a>
		<a name="{substring-before(concat(@name,'('),'(')}"></a>
		<xsl:value-of select="@name"/>(i), 
		<xsl:if test="position()=last()"> 
		  i=<xsl:value-of select="../@start"/>,<xsl:value-of select="../@end"/>
		</xsl:if>
	    </xsl:for-each>
	  </xsl:if>
	  <xsl:if test="name(.)='multidimensiongroup'">
	    <xsl:for-each select="multdimension">
	      <a name="{generate-id(.)}"></a>
		<a name="{substring-before(concat(@name,'('),'(')}"></a>
		<xsl:value-of select="@name"/>(i), 
		<xsl:if test="position()=last()"> 
		(<xsl:value-of select="../@indexes"/>)= (<xsl:value-of select="../@start"/>) . . . (<xsl:value-of select="../@end"/>)
		</xsl:if>
	    </xsl:for-each>
	  </xsl:if>

	</th>
	<td style="text-align: left; vertical-align: top; background: #ffffc3; padding: 2 2 2 5; ">
	  <xsl:value-of select="@type"/>
	</td>
      </tr>
      <xsl:apply-templates select="default"/> 
      <xsl:apply-templates select="dimensionality"/>
      <xsl:apply-templates select="units"/>
      <xsl:apply-templates select="status"/>
      <xsl:apply-templates select="see"/>
      <xsl:apply-templates select="info | options"/>
    </table>
    <xsl:call-template name="back_to_top"/>
  </xsl:template>

  <!--    *** VAR's elements ***  -->

  <!-- render a computed-sentinel value in square brackets, others unchanged.
       Keep in sync with helpdoc.d/txt.tcl. -->
  <xsl:template name="computed-sentinel">
    <xsl:param name="value"/>
    <xsl:variable name="v" select="normalize-space($value)"/>
    <xsl:choose>
      <xsl:when test="$v = 'from_pseudopotential' or $v = 'from_xml' or $v = 'from_environment' or $v = 'internal'">[<xsl:value-of select="$v"/>]</xsl:when>
      <xsl:otherwise><xsl:value-of select="$value"/></xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template match="default">
    <tr>
      <td style="text-align: right; vertical-align: top; background: #ffffc3; padding: 2 10 2 10; "> <i>Default:</i> </td>
      <td style="text-align: left;  vertical-align: top; background: #fff3d9; padding: 2 2 2 5; ">
	<xsl:choose>
	  <!-- conditional: render each "case" as header line + indented value -->
	  <xsl:when test="case">
	    <xsl:for-each select="case">
	      <xsl:choose>
		<xsl:when test="@test">
		  <i>if </i><xsl:call-template name="linkify-refs"><xsl:with-param name="text" select="@test"/></xsl:call-template><i>:</i>
		</xsl:when>
		<xsl:otherwise>
		  <i>otherwise:</i>
		</xsl:otherwise>
	      </xsl:choose>
	      <br/>
	      <xsl:text>&#160;&#160;&#160;</xsl:text>
	      <xsl:variable name="caseval">
		<xsl:apply-templates/>
	      </xsl:variable>
	      <xsl:call-template name="computed-sentinel">
		<xsl:with-param name="value" select="$caseval"/>
	      </xsl:call-template>
	      <br/>
	    </xsl:for-each>
	  </xsl:when>
	  <!-- plain default: render the text body -->
	  <xsl:otherwise>
	    <xsl:apply-templates/>
	  </xsl:otherwise>
	</xsl:choose>
      </td>
    </tr>
  </xsl:template>

  <xsl:template match="dimensionality">
    <tr>
      <td style="text-align: right; vertical-align: top; background: #ffffc3; padding: 2 10 2 10; "> <i>Dimensionality:</i> </td>
      <td style="text-align: left;  vertical-align: top; background: #fff3d9; padding: 2 2 2 5; ">
	<xsl:choose>
	  <!-- conditional: render each "case" as header line + indented value -->
	  <xsl:when test="case">
	    <xsl:for-each select="case">
	      <xsl:choose>
		<xsl:when test="@test">
		  <i>if </i><xsl:call-template name="linkify-refs"><xsl:with-param name="text" select="@test"/></xsl:call-template><i>:</i>
		</xsl:when>
		<xsl:otherwise>
		  <i>otherwise:</i>
		</xsl:otherwise>
	      </xsl:choose>
	      <br/>
	      <xsl:text>&#160;&#160;&#160;</xsl:text>
	      <xsl:call-template name="units-gloss">
		<xsl:with-param name="expr" select="."/>
	      </xsl:call-template>
	      <br/>
	    </xsl:for-each>
	  </xsl:when>
	  <!-- plain dimensionality: render the text body (may be a keyed clist) -->
	  <xsl:otherwise>
	    <xsl:call-template name="units-gloss-clist">
	      <xsl:with-param name="expr" select="."/>
	    </xsl:call-template>
	  </xsl:otherwise>
	</xsl:choose>
      </td>
    </tr>
  </xsl:template>

  <!-- render a units product-of-powers expression in conventional notation,
       e.g. "Ry bohr^-1" -> "Ry/bohr". Keep in sync with
       ::helpdoc::unitsNotation_ in dev-tools/helpdoc.d/txt.tcl. -->
  <xsl:template name="units-notation">
    <xsl:param name="expr"/>
    <xsl:variable name="u" select="normalize-space($expr)"/>
    <xsl:variable name="num">
      <xsl:call-template name="units-notation-rec">
        <xsl:with-param name="rest" select="$u"/>
        <xsl:with-param name="mode" select="'num'"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="den">
      <xsl:call-template name="units-notation-rec">
        <xsl:with-param name="rest" select="$u"/>
        <xsl:with-param name="mode" select="'den'"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:value-of select="$num"/>
    <xsl:value-of select="$den"/>
  </xsl:template>

  <!-- recursive helper for units-notation: mode='num' emits positive-power
       tokens, mode='den' emits negative-power tokens as "/base[^exp]". -->
  <xsl:template name="units-notation-rec">
    <xsl:param name="rest"/>
    <xsl:param name="mode"/>
    <xsl:param name="first" select="'yes'"/>
    <xsl:variable name="tok">
      <xsl:choose>
        <xsl:when test="contains($rest, ' ')">
          <xsl:value-of select="substring-before($rest, ' ')"/>
        </xsl:when>
        <xsl:otherwise><xsl:value-of select="$rest"/></xsl:otherwise>
      </xsl:choose>
    </xsl:variable>
    <xsl:variable name="is-neg" select="contains($tok, '^-')"/>
    <!-- 'yes' if this token belongs in the current pass -->
    <xsl:variable name="emit">
      <xsl:choose>
        <xsl:when test="$tok = ''">no</xsl:when>
        <xsl:when test="$mode = 'num' and not($is-neg)">yes</xsl:when>
        <xsl:when test="$mode = 'den' and $is-neg">yes</xsl:when>
        <xsl:otherwise>no</xsl:otherwise>
      </xsl:choose>
    </xsl:variable>
    <xsl:if test="$emit = 'yes'">
      <xsl:choose>
        <xsl:when test="$mode = 'num'">
          <xsl:if test="$first != 'yes'"><xsl:text> </xsl:text></xsl:if>
          <xsl:value-of select="$tok"/>
        </xsl:when>
        <xsl:otherwise>
          <xsl:variable name="base" select="substring-before($tok, '^-')"/>
          <xsl:variable name="exp" select="substring-after($tok, '^-')"/>
          <xsl:text>/</xsl:text>
          <xsl:value-of select="$base"/>
          <xsl:if test="$exp != '1'">^<xsl:value-of select="$exp"/></xsl:if>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:if>
    <xsl:if test="contains($rest, ' ')">
      <xsl:call-template name="units-notation-rec">
        <xsl:with-param name="rest" select="substring-after($rest, ' ')"/>
        <xsl:with-param name="mode" select="$mode"/>
        <xsl:with-param name="first">
          <xsl:choose>
            <xsl:when test="$first = 'yes' and $emit != 'yes'">yes</xsl:when>
            <xsl:otherwise>no</xsl:otherwise>
          </xsl:choose>
        </xsl:with-param>
      </xsl:call-template>
    </xsl:if>
  </xsl:template>

  <!-- render a <units>/<dimensionality> value as a human-readable phrase:
       $kind='dim' names the physical quantity, $kind='units' renders the
       unit notation. Keep in sync with ::helpdoc::unitsGloss in txt.tcl. -->
  <xsl:template name="units-gloss">
    <xsl:param name="expr"/>
    <xsl:param name="kind" select="'dim'"/>
    <!-- collapse surrounding/internal whitespace for a robust lookup -->
    <xsl:variable name="u" select="normalize-space($expr)"/>
    <xsl:choose>
      <xsl:when test="$kind = 'units'">
        <xsl:choose>
          <!-- atomic-unit-system composites: name the unit system only -->
          <xsl:when test="$u = 'bohr electron_mass^1/2 Ry^-1/2'">Rydberg atomic units</xsl:when>
          <xsl:when test="$u = 'bohr electron_mass^1/2 Hartree^-1/2'">Hartree atomic units</xsl:when>
          <xsl:when test="$u = 'Ry e^-1 bohr^-1'">Rydberg atomic units</xsl:when>
          <xsl:when test="$u = 'Hartree e^-1 bohr^-1'">Hartree atomic units</xsl:when>
          <!-- everything else: conventional unit notation -->
          <xsl:otherwise>
            <xsl:call-template name="units-notation">
              <xsl:with-param name="expr" select="$u"/>
            </xsl:call-template>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <!-- $kind='dim': name the physical quantity -->
      <xsl:otherwise>
        <xsl:choose>
          <xsl:when test="$u = 'Ry bohr^-1'">force (Ry/bohr)</xsl:when>
          <xsl:when test="$u = 'Hartree bohr^-1'">force (Hartree/bohr)</xsl:when>
          <xsl:when test="$u = 'Ry bohr^-3'">pressure (Ry/bohr^3)</xsl:when>
          <xsl:when test="$u = 'states eV^-1'">states/eV</xsl:when>
          <xsl:when test="$u = 'energy length^-1'">force</xsl:when>
          <xsl:when test="$u = 'energy length^-3'">pressure</xsl:when>
          <xsl:when test="$u = 'time^-1'">frequency</xsl:when>
          <xsl:when test="$u = 'length time^-1'">velocity</xsl:when>
          <xsl:when test="$u = 'energy charge^-1'">electric potential</xsl:when>
          <xsl:when test="$u = 'energy charge^-1 length^-1'">electric field</xsl:when>
          <xsl:when test="$u = 'charge length^-3'">charge density</xsl:when>
          <!-- single-token value or unrecognized composite: raw expression -->
          <xsl:otherwise><xsl:value-of select="$u"/></xsl:otherwise>
        </xsl:choose>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!-- render a <units>/<dimensionality> value that may be a keyed comma-list;
       keyed entries gloss to "<value> for index <i>". Keep in sync with
       ::helpdoc::unitsGlossClist in txt.tcl. -->
  <xsl:template name="units-gloss-clist">
    <xsl:param name="expr"/>
    <xsl:param name="kind" select="'dim'"/>
    <xsl:variable name="c" select="normalize-space($expr)"/>
    <xsl:choose>
      <!-- not a comma-list and not keyed: a single plain expression -->
      <xsl:when test="not(contains($c, ',')) and not(contains($c, ':'))">
        <xsl:call-template name="units-gloss">
          <xsl:with-param name="expr" select="$c"/>
          <xsl:with-param name="kind" select="$kind"/>
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <!-- Two passes: keyed overrides first (in listed order), then the
             single unkeyed default rendered as "... otherwise" last. -->
        <xsl:variable name="keyed">
          <xsl:call-template name="units-gloss-clist-rec">
            <xsl:with-param name="rest" select="$c"/>
            <xsl:with-param name="mode" select="'keyed'"/>
            <xsl:with-param name="kind" select="$kind"/>
            <xsl:with-param name="first" select="'yes'"/>
          </xsl:call-template>
        </xsl:variable>
        <xsl:variable name="unkeyed">
          <xsl:call-template name="units-gloss-clist-rec">
            <xsl:with-param name="rest" select="$c"/>
            <xsl:with-param name="mode" select="'unkeyed'"/>
            <xsl:with-param name="kind" select="$kind"/>
            <xsl:with-param name="first" select="'yes'"/>
          </xsl:call-template>
        </xsl:variable>
        <xsl:value-of select="$keyed"/>
        <xsl:if test="$keyed != '' and $unkeyed != ''">, </xsl:if>
        <xsl:value-of select="$unkeyed"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!-- recursive helper: consume the comma-list one token at a time, emitting
       only the entries matching $mode ('keyed' or 'unkeyed') -->
  <xsl:template name="units-gloss-clist-rec">
    <xsl:param name="rest"/>
    <xsl:param name="mode"/>
    <xsl:param name="kind" select="'dim'"/>
    <xsl:param name="first"/>
    <xsl:variable name="tok">
      <xsl:choose>
        <xsl:when test="contains($rest, ',')">
          <xsl:value-of select="normalize-space(substring-before($rest, ','))"/>
        </xsl:when>
        <xsl:otherwise><xsl:value-of select="normalize-space($rest)"/></xsl:otherwise>
      </xsl:choose>
    </xsl:variable>
    <xsl:variable name="tok-is-keyed" select="contains($tok, ':')"/>
    <!-- 'yes' if this token should be emitted in the current pass -->
    <xsl:variable name="emit">
      <xsl:choose>
        <xsl:when test="$tok = ''">no</xsl:when>
        <xsl:when test="$mode = 'keyed' and $tok-is-keyed">yes</xsl:when>
        <xsl:when test="$mode = 'unkeyed' and not($tok-is-keyed)">yes</xsl:when>
        <xsl:otherwise>no</xsl:otherwise>
      </xsl:choose>
    </xsl:variable>
    <xsl:if test="$emit = 'yes'">
      <xsl:if test="$first != 'yes'">, </xsl:if>
      <xsl:choose>
        <!-- keyed entry: "<i>:term" or "<lo>-<hi>:term" -->
        <xsl:when test="$tok-is-keyed">
          <xsl:variable name="key" select="normalize-space(substring-before($tok, ':'))"/>
          <xsl:variable name="val" select="normalize-space(substring-after($tok, ':'))"/>
          <xsl:call-template name="units-gloss">
            <xsl:with-param name="expr" select="$val"/>
            <xsl:with-param name="kind" select="$kind"/>
          </xsl:call-template>
          <xsl:choose>
            <xsl:when test="contains($key, '-')"> for indices <xsl:value-of select="$key"/></xsl:when>
            <xsl:otherwise> for index <xsl:value-of select="$key"/></xsl:otherwise>
          </xsl:choose>
        </xsl:when>
        <!-- unkeyed token: the default for all other indices -->
        <xsl:otherwise>
          <xsl:call-template name="units-gloss">
            <xsl:with-param name="expr" select="$tok"/>
            <xsl:with-param name="kind" select="$kind"/>
          </xsl:call-template>
          <xsl:text> otherwise</xsl:text>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:if>
    <xsl:if test="contains($rest, ',')">
      <xsl:call-template name="units-gloss-clist-rec">
        <xsl:with-param name="rest" select="substring-after($rest, ',')"/>
        <xsl:with-param name="mode" select="$mode"/>
        <xsl:with-param name="kind" select="$kind"/>
        <xsl:with-param name="first">
          <xsl:choose>
            <xsl:when test="$first = 'yes' and $emit != 'yes'">yes</xsl:when>
            <xsl:otherwise>no</xsl:otherwise>
          </xsl:choose>
        </xsl:with-param>
      </xsl:call-template>
    </xsl:if>
  </xsl:template>

  <xsl:template match="units">
    <tr>
      <td style="text-align: right; vertical-align: top; background: #ffffc3; padding: 2 10 2 10; "> <i>Units:</i> </td>
      <td style="text-align: left;  vertical-align: top; background: #fff3d9; padding: 2 2 2 5; ">
	<xsl:choose>
	  <!-- conditional: render each "case" as header line + indented value -->
	  <xsl:when test="case">
	    <xsl:for-each select="case">
	      <xsl:choose>
		<xsl:when test="@test">
		  <i>if </i><xsl:call-template name="linkify-refs"><xsl:with-param name="text" select="@test"/></xsl:call-template><i>:</i>
		</xsl:when>
		<xsl:otherwise>
		  <i>otherwise:</i>
		</xsl:otherwise>
	      </xsl:choose>
	      <br/>
	      <xsl:text>&#160;&#160;&#160;</xsl:text>
	      <xsl:call-template name="units-gloss">
		<xsl:with-param name="expr" select="."/>
		<xsl:with-param name="kind" select="'units'"/>
	      </xsl:call-template>
	      <br/>
	    </xsl:for-each>
	  </xsl:when>
	  <!-- plain units: render the text body (may be a keyed clist) -->
	  <xsl:otherwise>
	    <xsl:call-template name="units-gloss-clist">
	      <xsl:with-param name="expr" select="."/>
	      <xsl:with-param name="kind" select="'units'"/>
	    </xsl:call-template>
	  </xsl:otherwise>
	</xsl:choose>
      </td>
    </tr>
  </xsl:template>

  <xsl:template match="status">
    <tr>
      <td style="text-align: right; vertical-align: top; background: #ffffc3; padding: 2 10 2 10; "> <i>Status:</i> </td>
      <td style="text-align: left;  vertical-align: top; background: #fff3d9; padding: 2 2 2 5; ">
	<xsl:apply-templates/>
      </td>
    </tr>
  </xsl:template>

  <xsl:template match="see">
    <tr>
      <td style="text-align: right; vertical-align: top; background: #ffffc3; padding: 2 10 2 10; "> <i>See:</i> </td>
      <td style="text-align: left;  vertical-align: top; background: #fff3d9; padding: 2 2 2 5; ">
	<xsl:call-template name="tokenize_see"> 
	  <xsl:with-param name="refs" select="."/> 
	</xsl:call-template>  
      </td>
    </tr>
  </xsl:template>

  <xsl:template name="tokenize_see">
    <xsl:param name="refs"/>
    <xsl:variable name="first-ref" select="normalize-space(substring-before(concat($refs, ','), ','))" /> 
    <xsl:if test="$first-ref">
      <a href="#{normalize-space($first-ref)}"><xsl:value-of select="$first-ref"/></a>
      <xsl:if test="not($first-ref = normalize-space($refs))"><xsl:text>, </xsl:text></xsl:if>
      <xsl:call-template name="tokenize_see"> 
	<xsl:with-param name="refs" select="substring-after($refs,',')" /> 
      </xsl:call-template>    
    </xsl:if>  
  </xsl:template>
  
  <xsl:template match="info">
    <tr><td align="left" valign="top" colspan="2">
      <blockquote>
	<pre style="margin-bottom: -1em;">
	  <!--<xsl:apply-templates/>-->
	  <xsl:apply-templates/>
	</pre>
      </blockquote>
    </td></tr>
  </xsl:template>
  

  <!-- Options -->
  
  <xsl:template match="options">
    <tr><td align="left" valign="top" colspan="2">
      <blockquote>
	<xsl:apply-templates select="info | opt" mode="options"/>
      </blockquote>
    </td></tr>
  </xsl:template>

  <xsl:template match="info" mode="options">
    <pre style="margin-bottom: -1em;"><xsl:apply-templates/></pre>
  </xsl:template>

  <xsl:template match="opt" mode="options">
    <!--<span class="flag"><xsl:value-of select="@val"/></span><xsl:text> :</xsl:text>-->
    <dl style="margin-left: 1.5em;">
      <dt><tt><xsl:call-template name="tokenize_clist"><xsl:with-param name="clist" select="@val"/></xsl:call-template><xsl:if test="@alias != ''"><xsl:text>  (synonyms: </xsl:text><xsl:call-template name="tokenize_clist"><xsl:with-param name="clist" select="@alias"/></xsl:call-template><xsl:text>)</xsl:text></xsl:if><xsl:if test="normalize-space(.) != ''"> :</xsl:if></tt></dt>
      <dd><pre style="margin-top: 0em; margin-bottom: -1em;"><xsl:apply-templates/></pre></dd>
    </dl>
  </xsl:template>

  <xsl:template name="tokenize_clist">
    <xsl:param name="clist"/>
    <xsl:variable name="first-elem" select="normalize-space(substring-before(concat($clist, ','), ','))"/> 
    <xsl:if test="$first-elem">
      <span class="flag"><xsl:value-of select="$first-elem"/></span><xsl:if test="not($first-elem = normalize-space($clist))"><xsl:text>, </xsl:text></xsl:if>
      <xsl:call-template name="tokenize_clist"> 
	<xsl:with-param name="clist" select="substring-after($clist,',')" /> 
      </xsl:call-template>    
    </xsl:if>  
  </xsl:template>

  <!--    *** TABLE ***  -->

  <xsl:template match="table" mode="card_description">
    <xsl:apply-templates select="rows | cols" mode="table"/>   
  </xsl:template>

  <xsl:template match="rows | cols" mode="table">
    <xsl:apply-templates select="col | colgroup | row | rowgroup | optional | conditional" mode="table"/>
  </xsl:template>

  <xsl:template match="optional | conditional" mode="table">
    <xsl:apply-templates select="col | colgroup | row | rowgroup | optional | conditional" mode="table"/>
  </xsl:template>
  
  <xsl:template match="colgroup | rowgroup" mode="table">
    <xsl:if test="info != '' or options != '' or status != '' or see != ''">
      <table width="100%" style="border-color:   #b5b500; border-style: solid; border-width: 2; margin-bottom: 10; table-layout: auto; background-color: #FFFFFF;">
	<tr>
	  <th width="20%" align="left" valign="top" style="background: #ffff99; padding: 2 2 2 10; ">
	    <xsl:for-each select=".//col | .//row">
	      <a name="{@name}"><a name="{generate-id(.)}">
		<xsl:value-of select="@name"/>
	      </a></a>
	      <xsl:if test="not(position()=last())">
		<xsl:text>, </xsl:text>
	      </xsl:if>
	    </xsl:for-each>
	  </th>
	  <td style="text-align: left; vertical-align: top; background: #ffffc3; padding: 2 2 2 5; ">
	    <xsl:value-of select="@type"/>
	  </td>
	</tr>
	<xsl:apply-templates select="default"/> 
	<xsl:apply-templates select="dimensionality"/>
	<xsl:apply-templates select="units"/>
	<xsl:apply-templates select="status"/>
	<xsl:apply-templates select="see"/>
	<xsl:apply-templates select="info | options"/>
      </table>
      <xsl:call-template name="back_to_top"/>
    </xsl:if>
  </xsl:template>
  
  <xsl:template match="col | row" mode="table">
    <!--<xsl:if test="child::node() != ''">-->
    <xsl:if test="info != '' or options != '' or status != '' or see != ''">
      <table width="100%" style="border-color:   #b5b500; border-style: solid; border-width: 2; margin-bottom: 10; table-layout: auto; background-color: #FFFFFF;">
	<tr>
	  <th align="left" valign="top" width="20%" style="background: #ffff99; padding: 2 2 2 10; ">
	    <a name="{@name}"><a name="{generate-id(.)}">
	      <xsl:value-of select="@name"/>
	    </a></a>
	  </th>
	  <td style="text-align: left; vertical-align: top; background: #ffffc3; padding: 2 2 2 5; ">
	    <xsl:value-of select="@type"/>
	  </td>
	</tr>
	<xsl:apply-templates select="default"/> 
	<xsl:apply-templates select="dimensionality"/>
	<xsl:apply-templates select="units"/>
	<xsl:apply-templates select="status"/>
	<xsl:apply-templates select="see"/>
	<xsl:apply-templates select="info | options"/>
      </table>
      <xsl:call-template name="back_to_top"/>
    </xsl:if>
  </xsl:template>
  
  <!-- *** SECTION *** -->
  <xsl:template match="section">
    <blockquote>
      <a name="{generate-id(.)}">
	<h3><xsl:value-of select="@title"/></h3>
      </a>
      <xsl:apply-templates/>
    </blockquote>
  </xsl:template>

  <xsl:template match="subsection">
    <blockquote>
      <a name="{generate-id(.)}">
	<h4><xsl:value-of select="@title"/></h4>
      </a>
      <xsl:apply-templates/>
    </blockquote>
  </xsl:template>

  <xsl:template match="subsubsection">
    <blockquote>
      <a name="{generate-id(.)}">
	<h5><xsl:value-of select="@title"/></h5>
      </a>
      <xsl:apply-templates/>     
    </blockquote>
  </xsl:template>

  <xsl:template match="paragraph">
    <blockquote>
      <a name="{generate-id(.)}">
	<h6><xsl:value-of select="@title"/></h6>
	<xsl:apply-templates/>
      </a>
    </blockquote>
  </xsl:template>
  
  <xsl:template match="text">
    <blockquote>
      <pre><xsl:apply-templates/></pre>
    </blockquote>
  </xsl:template>

  <xsl:template match="link">
    <a href="{.}">
      <xsl:value-of select="."/>
    </a>
  </xsl:template>
  
  <xsl:template match="ref">
    <a href="#{.}">
      <xsl:value-of select="."/>
    </a>
  </xsl:template>

  <!-- linkify-refs: expand "@ref name" markup inside a plain string (e.g. a
       conditional <case> "test" attribute) into hyperlinks, leaving the rest
       of the string verbatim. -->
  <xsl:template name="linkify-refs">
    <xsl:param name="text"/>
    <xsl:choose>
      <xsl:when test="contains($text, '@ref ')">
	<xsl:value-of select="substring-before($text, '@ref ')"/>
	<xsl:variable name="after" select="substring-after($text, '@ref ')"/>
	<xsl:variable name="word">
	  <xsl:call-template name="leading-word">
	    <xsl:with-param name="s" select="$after"/>
	  </xsl:call-template>
	</xsl:variable>
	<a href="#{$word}"><xsl:value-of select="$word"/></a>
	<xsl:call-template name="linkify-refs">
	  <xsl:with-param name="text" select="substring($after, string-length($word) + 1)"/>
	</xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
	<xsl:value-of select="$text"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!-- leading-word: maximal prefix of $s made of identifier characters
       (letters, digits, "_", "%" for Fortran struct%member references). -->
  <xsl:template name="leading-word">
    <xsl:param name="s"/>
    <xsl:if test="string-length($s) &gt; 0">
      <xsl:variable name="c" select="substring($s, 1, 1)"/>
      <xsl:if test="contains('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_%', $c)">
	<xsl:value-of select="$c"/>
	<xsl:call-template name="leading-word">
	  <xsl:with-param name="s" select="substring($s, 2)"/>
	</xsl:call-template>
      </xsl:if>
    </xsl:if>
  </xsl:template>
  
  <xsl:template match="a">
    <!-- this is used by ref & link -->
    <xsl:copy-of select="."/>
  </xsl:template>
  
  <!-- support a subset of HTML tags (those specified in helpdoc.schema) -->
  
  <xsl:template match="b | i | tt | u | sub | sup | code | pre | br | hr | p | ol | ul | li | dl | dt | dd">
    <xsl:copy>
      <xsl:apply-templates/>
    </xsl:copy>
  </xsl:template>
  
</xsl:stylesheet>
