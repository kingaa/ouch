<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>
<! --- R-Forge Logo --- >
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<table>
<tr>
<td>
<table cellpadding="0" cellspacing="20" border="0">
<tr>
<td colspan="0" rowspan="0" align="left">
<img alt="ouch logo" src="ouch_ctr.cgi" height="79" width="52"></td>
<td colspan="0" rowspan="0" align="right">
<div align="right">
<h2><font color="#E32644">OUCH</font>: <font color="#E32644">O</font>rnstein-<font color="#E32644">U</font>hlenbeck models<br> 
for phylogenetic <font color="#E32644">C</font>omparative hypot<font color="#E32644">H</font>eses</h2>
</div>
</td>
</tr>
</table>
</td>
<td>
<ul>
<li><a href="./index.php">About OUCH</a></li>
<li><a target="_blank" href="http://cran.at.r-project.org/web/packages/ouch/ouch.pdf">OUCH Manual (PDF)</a></li>
<li><a href="./index.php?nav=bib">References to the literature</a></li>
<li><a href="http:<?php echo $domain; ?>/projects/<?php echo $group_name;?>">Development version of <tt>ouch</tt> (on r-forge)</a></li>
<li><a href="http://cran.at.r-project.org/web/packages/ouch/">Release version of <tt>ouch</tt> (on CRAN)</a></li>
<li><a href="./links.html">Related software</a></li>
</ul>
</td>
</tr>
</table>

<br>
<br>
<p>OUCH is software for phylogenetic comparative analysis. It is implemented in the <a href="http://www.r-project.org" target="_new">R language (www.r-project.org)</a> and can be downloaded from either the Comprehensive R Archive Network (<a href="http://cran.r-project.org/web/packages/ouch">cran.r-project.org</a>, stable version) or R-forge (<a href="http://ouch.r-forge.r-project.org">ouch.r-forge.r-project.org</a>, development version). </p>
<p>The method is based on the ideas of Thomas F. Hansen (see
T.&nbsp;F.&nbsp;Hansen, 1997. Stabilizing selection and the
comparative analysis of adaptation. <em>Evolution</em>,
<b>51</b>:1341-1351). It is explained fully in</p>
<dl>
<dd>Butler, M.A. and A.A. King, 2004. Phylogenetic comparative analysis: a modeling approach for adaptive evolution. <em>American Naturalist</em> <b>164</b>:683-695.</dd>
</dl>
<br>

<p> You can find the <strong>project summary page</strong> <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>
<p> The release version of the package is available on <a href="http://cran.at.r-project.org/web/packages/ouch/index.html">CRAN</a>.</p>

</body>
</html>
