#!/usr/local/bin/php
#$ -S /usr/local/bin/php
#$ -e log/
#$ -o log/
#$ -cwd
#$ -l s_vmem=4G,mem_req=4G,os7
<?php

require_once("./src/header.php");
ini_set("memory_limit","6144M");


function judge_concordance($ref_c, $alt_c, $ref_h, $alt_h){
	
	if( $ref_c + $alt_c ==0 |  $ref_h + $alt_h==0 | $ref_c + $ref_h==0 | $alt_c + $alt_h==0){return 0;}
	else{
		$chisqr_ue =( $ref_c + $alt_c + $ref_h + $alt_h) * ($ref_c * $alt_h -  $ref_h * $alt_c) * ($ref_c * $alt_h -  $ref_h * $alt_c);
		$chisqr_shita= ($ref_c + $alt_c) * ($ref_h + $alt_h) * ($ref_c + $ref_h) * ($alt_c + $alt_h);
		return $chisqr_ue/$chisqr_shita;
	}
}


function get_mutation($filename, $mode, $SampleID=null){
	$lines=file_get_lines($filename);

	// headerから、項目名と場所の対応付ける
	$header=array_shift($lines);
	$labels=explode("\t",rtrim($header));
	for($i=0;$i<count($labels);$i++)$pos[$labels[$i]]=$i;
	
	// ntog, gtonを保持する		
	static $ntog=null;if($ntog==null)$ntog=get_ntog();
	static $gton=null;if($gton==null)$gton=get_gton();
	
	// 各行を読み込み
	foreach($lines as $line){
		$d=explode("\t",rtrim($line));
		//
		$m=null;
		$m->line=rtrim($line);
		$m->mode=$mode;
		//
		$m->SampleID=$SampleID;
		if($d[$pos["SampleID"]]!=null)$m->SampleID=$d[$pos["SampleID"]];
		//
		$m->Func=$d[$pos["Func"]];
		//
		$m->GeneName=$d[$pos["Gene"]];
		if(preg_match("|^(.+)\(.+\)|",$m->GeneName,$match))$m->GeneName=$match[1];
		if(preg_match("|^(.+?);|is",$m->GeneName,$match))$m->GeneName=$match[1];
		if(preg_match("|^(.+?),(.+?)$|is",$m->GeneName,$match))$m->GeneName=strlen($match[1])<strlen($match[2])?$match[1]:$match[2];
		$m->GeneName=preg_replace("|\"|is","",$m->GeneName);	// "がついちゃってる名前がある
		//
		$m->ExonicFunc=$d[$pos["ExonicFunc"]];
		//
		$m->AAchange=$d[$pos["AAChange"]];
		$temp=explode(":",$m->AAchange);
		$m->nucleotideID=$temp[0];
		//
		$m->genomes1000=$d[$pos["1000G_ALL"]];
		$m->inhouseSNP=$d[$pos["inhouseSNP"]];
		$m->dbSNP131=$d[$pos["dbSNP131"]];
		$m->dbSNP132=$d[$pos["dbSNP132"]];
		//
		if($mode ==2){
			$m->chr=$d[$pos["Chr"]];
			$m->pos=$d[$pos["Start"]];
		}
		if($mode==1){
			$m->chr=$d[$pos["name"]];
			$m->pos=$d[$pos["position"]];
		}
		$m->ref=$d[$pos["Ref"]];
		$m->obs=$d[$pos["Obs"]];
		//
		$m->nmr=preg_replace("|\"|is","",$d[$pos["misRate_normal"]]);
		$m->tmr=preg_replace("|\"|is","",$d[$pos["misRate_tumor"]]);
		$m->p=preg_replace("|\"|is","",$d[$pos["p-value"]]);
		//
		$m->VariantGlue=$d[$pos["VariantGlue"]];
		//
		$m->indel=$d[$pos["indel"]];
		$m->misratio=preg_replace("|\"|","",$d[$pos["mismatch ratio"]]);
		//
		// GeneIDの取得、nucleotideIDで取れないときは名前で取る
		$m->GeneID						=$ntog[$m->nucleotideID];
		if($m->GeneID==null)$m->GeneID	=$ntog[$m->GeneName];
		// nucleotideIDが埋まってないときは、とりあえずGeneIDかGeneNameから埋める
		if($m->nucleotideID==null)$m->nucleotideID=$gton[$m->GeneID];
		if($m->nucleotideID==null)$m->nucleotideID=$gton[$m->GeneName];
		//
		/*
		if($m->VariantGlue!=null)$muts[$m->VariantGlue]=$m;
		else					 $muts[]=$m;
		*/
		$muts[]=$m;
	}
	return $muts;
}


$bamfile_list	=$argv[1];
$SNPs_list	=$argv[2];	if(!file_exists($SNPs_list))error("SNPs file list does not exist.");
$work_dir=$argv[3];
$mode =$argv[4];    // 1 for SNPs , 2 for Targeted seq

if( $mode ==1 ){
	$bam_path=$work_dir."/".$bamfile_list.".SNPbam";
}
if( $mode == 2){
	$bam_path=$work_dir."/".$bamfile_list.".Targetbam";
}
//echo $bam_path.PHP_EOL;
if(!file_exists($bam_path))error("bam file list does not exist.");


	

$bamfiles=file_get_lines($bam_path); 
array_shift($bamfiles);

//echo $bamfiles[0];



$samples[]=null;
$s=0;
foreach($bamfiles as $sample){
	$d = explode("\t", rtrim($sample));

	$m=null;
	$m->specimen = $d[0];
	$m->bam = $d[1];
	$samples[$s]=$m;

	++$s;
}

//echo "test".PHP_EOL;
//echo $s.PHP_EOL;


$calls[]=null;
$call_stat[]=null;
$n=0;
$first=0;
//$header_list="chr"."\t"."start"."\t"."end"."\t"."ref"."\t"."alt";
$header_list="chr"."\t"."start"."\t"."end"."\t"."ref"."\t"."alt";

//echo $header_list.PHP_EOL;


foreach(get_mutation($SNPs_list, $mode) as $m){
	$line1=$m->line;
	$line2=$m->line;


	for ($samp = 0; $samp < $s; $samp++) {
		if($first ==0 ){
			//$header_list.="\t".$samples[$samp]->specimen;
			$header_list.="\t".$samples[$samp]->specimen."_A"."\t".$samples[$samp]->specimen."_B";
		}
			
		$ref=$samples[$samp]->bam; 	$ref=rtrim($ref);
		if(!file_exists($ref))error("bam file does not exist.");


		//echo "$ref\t $m->chr \t $m->pos".PHP_EOL;

		if($m->mode==1){
			$p=get_pile_single2($ref, $m->chr, $m->pos);
		}
		if($m->mode==2){
			$p=get_pile_single($ref, $m->chr, $m->pos);
		}
		

		//echo "$p->pile".PHP_EOL;

		$ret=count_pile($p->pile,$m->ref,$m->obs);


		//$ref=$samples[$samp]->HLA_bam; ; 
		//if(!file_exists($ref))error("bam file does not exist.");

		//$p=get_pile_single($ref, $m->chr,$m->pos);
		//$ret_HLA=count_pile($p->pile,$m->ref,$m->obs);


		$line1.="\t".$ret["ref"]."\t".$ret["obs"];
		//$judge=judge_concordance($ret_conv["ref"], $ret_conv["obs"], $ret_HLA["ref"], $ret_HLA["obs"]);
		//$line2.="\t".$judge;
		print PHP_EOL;
	}
	print PHP_EOL;


	$calls[$first]=$line1.PHP_EOL;
	//$call_stat[$first]=$line2.PHP_EOL;
	++$first;
}


$header_list.=PHP_EOL;
$header_list2.=PHP_EOL;


//file_put_contents("output.txt", $header_list);
//file_put_contents("output.txt", implode("",$call_stat), $flags=FILE_APPEND);

$outpath=$work_dir."/".$bamfile_list."calls".$mode.".txt";
file_put_contents($outpath, $header_list);
file_put_contents($outpath, implode("",$calls), $flags=FILE_APPEND);




?>
