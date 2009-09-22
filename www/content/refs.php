<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<head>
<meta name="generator" content="Bluefish 1.0.7">
<title>Colleagues</title>
<link rel="stylesheet" href="./mystyle.css" type="text/css">
</head>
<body>
<?php
$dbfile = "content/refs.csv";
$abst = $_REQUEST["abst"];

function read_csv_file ($fname) {
    $recs = array();
    $row = 1;
    $handle = fopen($fname, "r");
    $keys = fgetcsv($handle, 1000, ",");

    while (($data = fgetcsv($handle, 1000, ",")) !== FALSE) {
        $recs[$row] = array_fill_keys($keys,$data[1]);
        foreach ($data as $col => $entry) {
            $recs[$row][$keys[$col]] = $entry;
        }
        foreach (array("title","pages","booktitle","journal") as $key) {
            $recs[$row][$key] = preg_replace("/\.$/","",$recs[$row][$key]);
        }
        $recs[$row]["abstract"] = preg_replace("/\\\\\"/","\"",$recs[$row]["abstract"]);
        $row++;
    }
    fclose($handle);
    return($recs);
}

function earlier ($a, $b) {
    if ($a["year"] == $b["year"]) {
        if ($a["author"]==$b["author"])
            return ($a["title"]<$b["title"]) ? -1 : 1;
        else
            return ($a["author"]<$b["author"]) ? -1 : 1;
    }
    return ($a["year"]>$b["year"]) ? -1 : 1;
}

$recs = read_csv_file($dbfile);
uasort($recs,"earlier");

if (!$abst) {
    echo "<h2>Selected Publications.</h2>\n\n";
    echo "<table border=\"0\" summary=\"List of publications\">\n";
    foreach ($recs as $row => $rec) {
        echo "<tr>\n";
        if ($rec["abstract"]) {
            echo "<td><a method=\"post\" href=\"index.php?nav=bib&abst=",$row,"\"><img border=\"0\" src=\"graphics/abstract.gif\" alt=\"Abstract\" height=\"48\" width=\"48\" align=\"left\"></a></td>\n";
        } else {
            echo "<td></td>\n";
        }
        echo "<td>";
        echo $rec["author"]," (",$rec["year"],") ",$rec["title"],".&nbsp;";
        if ($rec["journal"]) {
            echo "<i>",$rec["journal"],"</i>,&nbsp;";
        } 
        if ($rec["volume"]) {
            echo "<b>",$rec["volume"],"</b>:&nbsp;";
        } 
        if ($rec["booktitle"]) {
            echo "<i>",$rec["booktitle"],"</i>&nbsp;";
        }
        if ($rec["editor"]) {
            echo "(edited by ",$rec["editor"],")&nbsp;";
        }
        if ($rec["pages"]) {
            echo $rec["pages"],".&nbsp;";
        } 
        if ($rec["school"]) {
            echo "Ph.D. thesis, ",$rec["school"],".&nbsp;";
        }
        if ($rec["doi"]) {
            echo "<a href=\"",$rec["doi"],"\" target=\"_blank\">DOI</a>&nbsp;";
        }
        echo "</td></tr>\n";
    }
    echo "</table>";
} else {
    $rec = $recs[$abst];
    echo "<hr>\n";
    echo "<center><h3>",$rec["title"],"</h3></center>\n";
    echo "<center>",$rec["author"],"</center>\n";
    echo "<center>\n";
    if ($rec["journal"]) {
        echo "<i>",$rec["journal"],"</i>,&nbsp;";
    } 
    if ($rec["volume"]) {
        echo "<b>",$rec["volume"],"</b>:&nbsp;";
    } 
    if ($rec["pages"]) {
        echo $rec["pages"],",&nbsp;";
    } 
    if ($rec["booktitle"]) {
        echo "<i>",$rec["booktitle"],"</i>&nbsp;";
    }
    if ($rec["school"]) {
        echo "Ph.D. thesis, ",$rec["school"],".&nbsp;";
    }
    echo $rec["year"],".\n";
    echo "</center>\n";
    echo "<p><b>Abstract</b></p>\n";
    echo "<p>",$rec["abstract"],"</p>\n";
    echo "<hr><br>\n";
    if ($rec["doi"]) {
        echo "<b>The official version of the paper is <a href=\"",$rec["doi"],"\" target=\"_blank\">here</a>.</b>&nbsp;&nbsp;\n";
    }
}
?>
<hr width="100%">
</body>
</html>
