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
    <xsl:apply-templates select=".//var | .//dimension | .//list" mode="toc"/></p>
  </xsl:template>

  <xsl:template match="namelist | card" mode="toc">
    <p><a href="#{generate-id(.)}">
      <xsl:if test="name(.)='namelist'">&#38;</xsl:if>
      <xsl:value-of select="@name"/>
    </a></p>
    <xsl:if test=".//var != '' or
		  .//dimension != '' or
		  .//list != '' or
		  .//col != '' or
		  .//row != ''">
      <blockquote>
	<xsl:apply-templates select=".//var | .//dimension | .//list | .//col | .//row" mode="toc"/>
      </blockquote>
    </xsl:if>
  </xsl:template>

  <xsl:template match="var | dimension" mode="toc">
    <xsl:if test="info != '' or 
		  status != '' or 
		  see    != '' or
		  options != '' or
		  ../../vargroup/info != '' or 
		  ../../dimensiongroup/info != '' or
		  ../../dimensiongroup/options != ''">
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
    <xsl:apply-templates select="optional | conditional | var | keyword | vargroup | list" mode="syntax"/>
    <br/>
  </xsl:template>	


  <!-- card//syntax//optional -->      

  <xsl:template match="optional" mode="syntax">
    <!--<div style="background: #eeeeee; color: #555555;">-->
    <xsl:text> { </xsl:text>
    <xsl:apply-templates select="line | var | list | keyword | table" mode="syntax"/>
    <xsl:text> } </xsl:text>	    
    <!--</div>-->
  </xsl:template>	      

  <!-- card//syntax//conditional -->      

  <xsl:template match="conditional" mode="syntax">
    <xsl:text> [ </xsl:text>
    <xsl:apply-templates select="line | var | list | keyword | table" mode="syntax"/>
    <xsl:text> ] </xsl:text>	    
  </xsl:template>	      

  <!-- card//syntax//keyword -->      

  <xsl:template match="keyword" mode="syntax">
    <b><xsl:value-of select="@name"/></b><xsl:text>&#160;&#160;</xsl:text>
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
    <xsl:apply-templates select="var" mode="syntax"/>
  </xsl:template>

  <!-- card//syntax//table -->      

  <xsl:template match="table" mode="syntax">
    <a name="{generate-id(.)}"></a>
    <table>
      <xsl:apply-templates select="rows | cols" mode="syntaxTableMode"/>
    </table>    
  </xsl:template>

  <!-- sytntax//table/rows -->

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
  
  <!-- *** VARGROUP | DIMENSIONGROUP *** -->

  <xsl:template match="vargroup | dimensiongroup" mode="card_description">
    <!--<xsl:if test="child::node() != ''">-->
    <xsl:if test="info != '' or options != '' or status != '' or see != ''">
      <xsl:apply-templates select="."/>
    </xsl:if>
  </xsl:template>

  <!--    *** VAR | DIMENSION | LIST | FLAG ***  -->

  <xsl:template match="var | list | dimension | flag" mode="card_description">
    <!--<xsl:if test="child::node() != ''">-->
    <xsl:if test="info != '' or options != '' or status != '' or see != ''">
      <xsl:apply-templates select="."/>
    </xsl:if>
  </xsl:template>

  <xsl:template match="var | list | dimension | flag">
    <xsl:if test="name(..) != 'vargroup' and name(..) != 'dimensiongroup'">
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
	<xsl:apply-templates select="status"/>
	<xsl:apply-templates select="see"/>
	<xsl:apply-templates select="info | options"/>
      </table>
      <xsl:call-template name="back_to_top"/>
    </xsl:if>
  </xsl:template>

  <!-- *** VARGROUP | DIMENSIONGROUP *** -->

  <xsl:template match="vargroup | dimensiongroup">
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
	</th>
	<td style="text-align: left; vertical-align: top; background: #ffffc3; padding: 2 2 2 5; ">
	  <xsl:value-of select="@type"/>
	</td>
      </tr>
      <xsl:apply-templates select="default"/> 
      <xsl:apply-templates select="status"/>
      <xsl:apply-templates select="see"/>
      <xsl:apply-templates select="info | options"/>
    </table>
    <xsl:call-template name="back_to_top"/>
  </xsl:template>

  <!--    *** VAR's elements ***  -->

  <xsl:template match="default">
    <tr>
      <td style="text-align: right; vertical-align: top; background: #ffffc3; padding: 2 10 2 10; "> <i>Default:</i> </td>
      <td style="text-align: left;  vertical-align: top; background: #fff3d9; padding: 2 2 2 5; ">
	<xsl:apply-templates/>
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
      <dt><tt><xsl:call-template name="tokenize_clist"><xsl:with-param name="clist" select="@val"/></xsl:call-template><xsl:if test="normalize-space(.) != ''"> :</xsl:if></tt></dt>
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
