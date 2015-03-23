<?php

/* Usage
 Grab some XML data, either from a file, URL, etc. however you want. Assume storage in $strYourXML;

 $objXML = new xml2Array();
 $arrOutput = $objXML->parse($strYourXML);
 print_r($arrOutput); //print it out, or do whatever!
  
*/
class xml2Array {
    
    var $arrOutput = array();
	var $pathName = array();
	var $pathUrl = array();
	var $categories = array();
	var $alphabetical = array();
	var $pageDescriptions = array();
	var $indiceXML;
	var $fileIndex;
	var $htmlTemplate;
	var $resParser;
    var $strXmlData;
	var $dirTarget;
	var $dirBase;
    
    function parse($strInputXML) {
    
            $this->resParser = xml_parser_create ();
            xml_set_object($this->resParser,$this);
            xml_set_element_handler($this->resParser, "tagOpen", "tagClosed");
            
            xml_set_character_data_handler($this->resParser, "tagData");
        
            $this->strXmlData = xml_parse($this->resParser,$strInputXML );
            if(!$this->strXmlData) {
               die(sprintf("XML error: %s at line %d",
            xml_error_string(xml_get_error_code($this->resParser)),
            xml_get_current_line_number($this->resParser)));
            }
                            
            xml_parser_free($this->resParser);
            
            return $this->arrOutput;
    }
    function tagOpen($parser, $name, $attrs) {
       $tag=array("name"=>$name,"attrs"=>$attrs); 
       array_push($this->arrOutput,$tag);
    }
    
    function tagData($parser, $tagData) {
       if(trim($tagData)) {
            if(isset($this->arrOutput[count($this->arrOutput)-1]['tagData'])) {
                $this->arrOutput[count($this->arrOutput)-1]['tagData'] .= str_replace("\t","",$tagData);
            } 
            else {
                $this->arrOutput[count($this->arrOutput)-1]['tagData'] = str_replace("\t","",$tagData);
            }
       }
    }
    
    function tagClosed($parser, $name) {
       $this->arrOutput[count($this->arrOutput)-2]['children'][] = $this->arrOutput[count($this->arrOutput)-1];
       array_pop($this->arrOutput);
    }
	
	function getNavigation()
	{
		$result = '<ul><li itemscope="" itemtype="http://www.data-vocabulary.org/Breadcrumb" itemprop="breadcrumb" class="breadcrumb_first"><a href="#" onclick="window.location.href=home;" itemprop="url"><span itemprop="title">R2013xxx</span></a></li>';
		for($i=0;$i<count($this->pathName);$i++)
		{
			$result .= '<li itemscope="" itemtype="http://www.data-vocabulary.org/Breadcrumb" itemprop="breadcrumb" class="breadcrumb_product"><a href="'.$this->pathUrl[$i].'" itemprop="url"><span itemprop="title">'.$this->pathName[$i].'</span></a></li>';
		}
		$result .= '</ul>';
		return $result;
	}
	
	function trasformFiles($current) {
		$temp = explode("#",$current["attrs"]["TARGET"]);
		if(file_exists($this->dirBase.$temp[0]))
		{
			$temp_file = $this->htmlTemplate;
			//echo "---------------------------------\n";
			//echo $temp[0]."\n";
			//echo "---------------------------------\n";
			$temp_file = str_replace("URL_CORRENTE",$temp[0],$temp_file);
			$content = implode("",file($this->dirBase.$temp[0]));
			if(preg_match('<!--FScategory:(.+)-->',$content,$matches))
			{
				//echo $matches[1]." ".$current["tagData"]."\n";
				$this->categories[$matches[1]][$current["tagData"]] = $temp[0];
				$this->alphabetical[$current["tagData"]] = $temp[0];
			}
			if(preg_match_all('/"purpose">(.+)<\/p>/Us',$content,$matches))
			{
				$this->pageDescriptions[$temp[0]] = $matches[1][0];
			}

		
			/*if (preg_match('/<title>(.+)<\/title>/',$content,$matches) && isset($matches[1] ))
			{
				$title = $matches[1];
				echo $title;*/
				if($temp[0] == "fsda_product_page.html")
				{
					$temp_file = str_replace("BGCOLOR","#436699",$temp_file);
				}
				else
				{
					$temp_file = str_replace("BGCOLOR","#FFFFFF",$temp_file);
				}
				$temp_file = str_replace("TITOLO",$current["tagData"],$temp_file);
				$temp_file = str_replace("DESCRIZIONE",$current["tagData"],$temp_file);
				//echo "---------------------------------\n";
			//}
			/*if (preg_match('/<body[^>]*>(.+)<\/body>/s',$content,$matches) && isset($matches[1] ))
			{	
				$body = $matches[1];
				//echo $body;
				$temp_file = str_replace("DOCUMENTO",$body,$temp_file);
			}*/
			//echo "\n";
			//$temp_file = str_replace("NAVIGATION",$this->getNavigation(),$temp_file);
			//$fd = fopen($this->dirTarget.$temp[0],"w+");
			//fwrite($fd,$temp_file);
			//fclose($fd);
			//copy($this->dirBase.$temp[0],$this->dirTarget.(str_replace(".html",".htm",$temp[0])));
			
			//copy($temp[0],$this->dirTarget.$temp[0]);
			//echo "---------------------------------\n";
		}
		else
		{
			//echo $temp[0]." KO\n";
		}
		//echo implode("->",$this->path)."\n";
		if(isset($current["children"]))
		{
			for($i = 0;$i<count($current["children"]);$i++)
			{
				$this->pathName[] = $current["children"][$i]["tagData"];
				$this->pathUrl[] = $current["children"][$i]["attrs"]["TARGET"];
				$this->trasformFiles($current["children"][$i]);
				array_pop($this->pathName);
				array_pop($this->pathUrl);
			}
		}
		
	}
	
	function createIndice()
	{
		$indice = implode("",file($this->fileIndex));
		$indice_one_line = str_replace("\t","",$indice);
		$indice_one_line = str_replace("\r\n","",$indice_one_line);
		$indice_one_line = str_replace("\n","",$indice_one_line);
		$indice_one_line = str_replace("<name>","",$indice_one_line);
		$indice_one_line = str_replace("</name>","",$indice_one_line);
		if(preg_match_all('/<tocreference target="(.+)"><\/tocreference>/',$indice,$matches))
		{
			for($i=0;$i<count($matches[1]);$i++)
			{
				$file_to_include = "";
				if(file_exists($this->dirBase.$matches[1][$i]))
				{
					$file_to_include = implode("",file($this->dirBase.$matches[1][$i]));
					$file_to_include = str_replace("\t","",$file_to_include);
					$file_to_include = str_replace("\r\n","",$file_to_include);
					$file_to_include = str_replace("\n","",$file_to_include);
					$file_to_include = str_replace("<name>","",$file_to_include);
					$file_to_include = str_replace("</name>","",$file_to_include);
					$file_to_include = preg_replace('/<purpose>.+?<\/purpose>/',"",$file_to_include);
					if( preg_match('/<toc version="2.0">(.+)<\/toc>/',$file_to_include,$matches_sec) && isset($matches_sec[1] ))
					{
						$indice_one_line = str_replace('<tocreference target="'.$matches[1][$i].'"></tocreference>',$matches_sec[1],$indice_one_line);
					}
					else
					{
						$indice_one_line = str_replace('<tocreference target="'.$matches[1][$i].'"></tocreference>',"",$indice_one_line);
					}
				}
			}
		}
		$this->indiceXML = $indice_one_line;
		$fd = fopen($this->dirTarget."indice.js","w+");
		fwrite($fd,"var indice = '".$indice_one_line."';");
		fclose($fd);
	}
	
	function generateCategorical()
	{
		//$temp_file = $this->htmlTemplate;
		$temp_file = $this->htmlTemplate;
		$temp_file = str_replace("TITOLO","Function Reference",$temp_file);
		$body = '<a name="top_of_page"></a>
<table border="0" cellpadding="0" cellspacing="0" class="nav" summary="Navigation aid" width="100%">
	<tr>
		<td valign="baseline"><b>FSDA Toolbox</b></td>
		<td align="right" valign="baseline"><a href="multivariate.html">
		<img align="bottom" alt="Capability Studies" border="0" src="images_help/b_prev.gif"></a>&nbsp;&nbsp;
		<a href="datasets.html">
		<img align="bottom" alt="Class Reference" border="0" src="images_help/b_next.gif"></a></td>
	</tr>
</table>
<table border="0" cellpadding="0" cellspacing="0" class="refpartnertable" width="100%">
	<tr>
		<td>
		<h1 class="refpartnerheading"><a name="bq_w_hm"></a>Function Reference</h1>
		</td>
		<td class="refpartnerlink"><a href="alphabetical.html">
		<img align="bottom" alt="Click for Alphabetical List" border="0" src="images_help/more_arrows.gif"> 
		Alphabetical List</a></td>
	</tr>
</table>
<hr>';
		
		
		foreach($this->categories as $key => $value)
		{
			$body = $body.'<h2 class="categorytitle">'.$key.'</h2>
			<a class="indexterm" name="zmw57dd0e24781"></a><a name="br0bjd1-1"></a>
			<p class="pagenavlink">
			<table border="1">
				<tr>
					<td>
						<table border="0" cellpadding="5" class="subcategorylist">';
			foreach($value as $key => $data)
			{
				$body = $body.'<tr>
					<td width="300"><a href="'.$data.'">'.$key.'</a></td>
					<td>'.$this->pageDescriptions[$data].'</td>
				</tr>';
			}
			//$body = $body."<tr><td style='height: 26px' valign='top'><a href='".$value."'>".$key."</a></td><td style=\"height: 26px\">".$this->pageDescriptions[$value]."</td></tr>";
			$body = $body.'</table>
					</td>
				</tr>
			</table><p class="pagenavlink">
			<a href="#top_of_page">Back to Top</a><p></p>';
		}

		$temp_file = str_replace("DOCUMENTO",$body,$temp_file);
		$fd = fopen($this->dirTarget."function_reference.html","w+");
		fwrite($fd,$temp_file);
		fclose($fd);
		
		
		
	}
	
	function generateAlphabetical()
	{
		$temp_file = $this->htmlTemplate;
		$temp_file = str_replace("TITOLO","Functions - Alphabetical List",$temp_file);
		asort($this->alphabetical);
		
		$body = '<a name="top_of_page"></a>
<table border="0" cellpadding="0" cellspacing="0" class="nav" summary="Navigation aid" width="100%">
	<tr>
		<td valign="baseline"><b>FSDA Toolbox</b></td>
		<td align="right" valign="baseline"><a href="function_reference.html">
		<img align="bottom" alt="Class Reference" border="0" src="images_help/b_prev.gif"></a>&nbsp;&nbsp;
		<a href="addt.html">
		<img align="bottom" alt="addt" border="0" src="images_help/b_next.gif"></a></td>
	</tr>
</table>
<table border="0" cellpadding="0" cellspacing="0" class="refpartnertable" width="100%">
	<tr>
		<td>
		<h1 class="refpartnerheading"><a name="bqzigzm-1"></a>Functions &minus; Alphabetical 
		List</h1>
		</td>
		<td class="refpartnerlink"><a href="function_reference.html">
		<img align="bottom" alt="Click for Categorical List" border="0" src="images_help/more_arrows.gif"> 
		Function Reference</a><br></td>
	</tr>
</table>
<hr>
<p></p><table border="0" cellpadding="2" cellspacing="0" width="100%">
		<colgroup>
		<col width="250"><col>
		</colgroup>';
	
		foreach($this->alphabetical as $key => $value)
		{
			$body = $body."<tr><td style='height: 26px' valign='top'><a href=\"".$value."\">".$key."</a></td><td style=\"height: 26px\">".$this->pageDescriptions[$value]."</td></tr>\n";
		}
		$body = $body."</table>";		
		$temp_file = str_replace("DOCUMENTO",$body,$temp_file);
		$fd = fopen($this->dirTarget."alphabetical.html","w+");
		fwrite($fd,$temp_file);
		fclose($fd);
	}


	function transcodeXML()
	{
		$files = glob($this->dirBase."/*.xml");
		foreach($files as $file)
		{
			$temp = explode("/",$file);
			$file_name = $temp[1];
			$file_raw = implode("",file($file));
			$file_raw = str_replace(".html",".htm",$file_raw);
			$fd = fopen($this->dirTarget.$file_name,"w+");
			fwrite($fd,$file_raw );
			fclose($fd);			
			
		}
	
	}
	
}

if(count($argv) == 1)
{
	echo "Questo programma serve aia trasformare i file html della vecchia documentazione\nmatlab nella nuova versione.\n";
	echo "Vengono richiesti 4 parametri:\n";
	echo "\n";
	echo "1) Il percorso alla directory della vecchia documentazione\n";
	echo "2) Il percorso alla directory in cui generare la nuova documentazione\n";
	echo "3) Il file indice.xml in versione toc\n";
	echo "4) Il percorso al file template della dcoumentazione nella nuova versione\n";
	echo "\n";
	echo "Esempio: trasform_help.exe ./old_helps ../new_helps indice.xml ../template.html\n";
	echo "\n";
	echo "Note: Questo file va eseguito nella directory in cui si trovano i vecchi file\ndella documentazione, insieme al file indice.xml in versione toc e\na tutti i file toc referenziati";
	die();
}

// Controllo se la directory in cui trovare i vecchi file esiste  
if(!file_exists($argv[1])) die("La Directory in cui trovare i vecchi files help non esiste");

// Controllo se la directory in cui copiare i file esiste  
if(!file_exists($argv[2])) die("La Directory in cui salvare i nuovi files help non esiste");

// Controllo che il file .xml dato come indice esista
if(!file_exists($argv[1].$argv[3])) die("L'indice XML non esiste");

// Controllo che il file .xml dato come indice esista
if(!file_exists($argv[4])) die("Il template html non esiste");


$objXML = new xml2Array();
$objXML->dirBase = $argv[1];
$objXML->dirTarget = $argv[2];
$objXML->fileIndex = $argv[1].$argv[3];
$objXML->htmlTemplate = implode("",file($argv[4]));

// Creazione indice.js

$objXML->createIndice();


$arrOutput = $objXML->parse($objXML->indiceXML);

//print_R($arrOutput);

/*''*/
$objXML->pathName[] = $arrOutput[0]["children"][0]["tagData"];
$objXML->pathUrl[] = $arrOutput[0]["children"][0]["attrs"]["TARGET"];
$objXML->trasformFiles($arrOutput[0]["children"][0]);

$objXML->generateAlphabetical();
$objXML->generateCategorical();
$objXML->transcodeXML();



//$alpha["attrs"]["TARGET"] = "alphabetical.html";
//$objXML->pathName[] = "alphabetical.html";
//$alpha["tagData"] = "Functions - Alphabetical List";
//$objXML->pathUrl[] = "Functions - Alphabetical List";
//$objXML->trasformFiles($alpha);
//print_R($objXML->pageDescriptions);
//print_R($objXML->categories);
?>
