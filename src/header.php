#!/usr/bin/php
#$ -S /usr/bin/php
#$ -e log
#$ -o log
#$ -cwd

<?php

//------------------------------------------------------------------------------
// 初期化（・ｗ・）ノ
//------------------------------------------------------------------------------

set_time_limit(86400);	// 24h
ini_set("memory_limit","2048M");
ini_set("error_reporting","E_ALL & ~E_WARNING & ~E_NOTICE");
//ini_set("error_reporting","E_ALL & ~E_NOTICE");

//------------------------------------------------------------------------------
// 関数（・ｗ・）ノ
//------------------------------------------------------------------------------

// get_mutationlist 変異リスト読み込む、列の順番に依存しない、最強版の予定
function get_mutationlist($filename,$SampleID=null){
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
		$m->chr=$d[$pos["Chr"]];
		$m->pos=$d[$pos["Start"]];
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

/*
// get_mutationlist mutationlist.txtを開いて変数に格納（・ｗ・）GeneIDも読み込むよ。
function get_mutationlist($filename="mutationlist.txt"){
	$ntog=get_ntog();
	$gton=get_gton();
	foreach(file_get_lines($filename) as $line){
		$p=explode("\t",rtrim($line));
		// RCC102<>exonic<>FOXP1<>nonsynonymous SNV<>NM_001012505:c.G64C:p.G22R<>…
		$m=null;
		$m->SampleID		=$p[0];
		$m->Func			=$p[1];
		$m->GeneName		=$p[2];
		if(preg_match("|^(.+?);|is",$m->GeneName,$match))$m->GeneName=$match[1];
		if(preg_match("|^(.+?),(.+?)$|is",$m->GeneName,$match))$m->GeneName=strlen($match[1])<strlen($match[2])?$match[1]:$match[2];
		$m->GeneName=preg_replace("|\"|is","",$m->GeneName);	// "がついちゃってる名前がある
		$m->ExonicFunc		=$p[3];
		$m->AAchange		=$p[4];
		// $ref=$p[21];$obs=$p[22];$tumor=$p[23];$normal=$p[24];
		// nucleotideID,ref,obsの取得
		$d=explode(":",$m->AAchange);
		$m->nucleotideID	=$d[0];
		// ref/obs/posの取得
		$m->pos				=$p[20];
		$m->ref				=$p[22];
		$m->obs				=$p[23];
		$m->ref=$m->obs		='@';
		if(preg_match("/(synonymous|stopgain)/is",$m->ExonicFunc)){
			preg_match("|([ACGT])(\d+)([ACGT])|is",$d[1],$match);
			$m->ref			=$match[1];
			$m->pos			=$match[2];
			$m->obs			=$match[3];
			// deletionによるstopgainの場合はマッチしないので、再度当てる
			if($m->ref==null){$m->ref=$m->obs='@';}
		}
		// GeneIDの取得、nucleotideIDで取れないときは名前で取る
		$m->GeneID						=$ntog[$m->nucleotideID];
		if($m->GeneID==null)$m->GeneID	=$ntog[$m->GeneName];
		// nucleotideIDが埋まってないときは、とりあえずGeneIDかGeneNameから埋める
		if($m->nucleotideID==null)$m->nucleotideID=$gton[$m->GeneID];
		if($m->nucleotideID==null)$m->nucleotideID=$gton[$m->GeneName];
		//
		$muts[]=$m;
	}	
	return $muts;
}
*/
// get_mutationlist_s 佐藤先生の変異リストを開いて変数に格納（・ｗ・）GeneIDも読み込むよ。
function get_mutationlist_s($filename="RCCmutationlist_120106.txt"){
//function get_mutationlist_s($filename="RCC_mutation list.111219.txt"){
	$ntog=get_ntog();
	$gton=get_gton();
	foreach(file_get_lines($filename,1) as $line){
		$p=explode("\t",rtrim($line));
		// RCC102<>exonic<>FOXP1<>nonsynonymous SNV<>NM_001012505:c.G64C:p.G22R<>…
		$m=null;
		$m->SampleID		=$p[0];
		$m->Func			=$p[1];
		$m->GeneName		=$p[2];
		if(preg_match("|^(.+?);|is",$m->GeneName,$match))$m->GeneName=$match[1];
		if(preg_match("|^(.+?),(.+?)$|is",$m->GeneName,$match))$m->GeneName=strlen($match[1])<strlen($match[2])?$match[1]:$match[2];
		$m->GeneName=preg_replace("|\"|is","",$m->GeneName);	// "がついちゃってる名前がある
		$m->ExonicFunc		=$p[3];
		$m->AAchange		=$p[4];
		// $ref=$p[21];$obs=$p[22];$tumor=$p[23];$normal=$p[24];
		// nucleotideID,ref,obsの取得
		$d=explode(":",$m->AAchange);
		$m->nucleotideID	=$d[0];

		// ref/obs/posの取得
		$m->pos				=$p[20];
		$m->ref				=$p[22];
		$m->obs				=$p[23];
		
		//
		$m->Chr				=$p[19];
		// tumor mismatch ratioの取得
		$m->tmr				=$p[26];

		// 2012/02/21 （・ｗ・）祭
		$ACGT=array("A"=>0,"C"=>1,"G"=>2,"T"=>3);
		$temp=preg_replace("|\"|is","",$p[24]);
		//print $temp." ";
		$t=explode(",",$temp);
		$m->refcount_tumor	=$t[$ACGT[$m->ref]];
		$m->obscount_tumor	=$t[$ACGT[$m->obs]];

		/*
		$m->ref=$m->obs		='@';
		if(preg_match("/(synonymous|stopgain)/is",$m->ExonicFunc)){
			preg_match("|([ACGT])(\d+)([ACGT])|is",$d[1],$match);
			$m->ref			=$match[1];
			$m->pos			=$match[2];
			$m->obs			=$match[3];
			// deletionによるstopgainの場合はマッチしないので、再度当てる
			if($m->ref==null){$m->ref=$m->obs='@';}
		}
		*/
		// GeneIDの取得、nucleotideIDで取れないときは名前で取る
		$m->GeneID						=$ntog[$m->nucleotideID];
		if($m->GeneID==null)$m->GeneID	=$ntog[$m->GeneName];
		// nucleotideIDが埋まってないときは、とりあえずGeneIDかGeneNameから埋める
		if($m->nucleotideID==null)$m->nucleotideID=$gton[$m->GeneID];
		if($m->nucleotideID==null)$m->nucleotideID=$gton[$m->GeneName];
		//
		$muts[]=$m;
	}	
	return $muts;
}

/*
// 昔の手書き気味のリストを読んでいたころのやつ
// 佐藤先生のリスト読む
function get_mutationlist_s($filename="mutationlist_s.txt"){
	$ntog=get_ntog();
	$gton=get_gton();
	$lines=file_get_lines($filename);
	array_shift($lines);
	foreach($lines as $line){
		$p=explode("\t",rtrim($line));
		// RCC102<>exonic<>FOXP1<>nonsynonymous SNV<>NM_001012505:c.G64C:p.G22R<>…
		$m=null;
		$m->SampleID		=$p[0];
		$m->Func			="";
		$m->GeneName		=$p[1];
		if(preg_match("|^(.+?);|is",$m->GeneName,$match))$m->GeneName=$match[1];
		if(preg_match("|^(.+?),(.+?)$|is",$m->GeneName,$match))$m->GeneName=strlen($match[1])<strlen($match[2])?$match[1]:$match[2];
		$m->ExonicFunc		=$p[2];
		$m->AAchange		=$p[3];
		//
		// $ref=$p[21];$obs=$p[22];$tumor=$p[23];$normal=$p[24];
		// nucleotideID,ref,obsの取得
		$d=explode(":",$m->AAchange);
		$m->nucleotideID	=$d[0];
		$m->ref=$m->obs		='@';
		if(preg_match("/(synonymous|stopgain)/is",$m->ExonicFunc)){
			preg_match("|([ACGT])(\d+)([ACGT])|is",$d[1],$match);
			$m->ref			=$match[1];
			$m->pos			=$match[2];
			$m->obs			=$match[3];
			// deletionによるstopgainの場合はマッチしないので、再・x当てる
			if($m->ref==null){$m->ref=$m->obs='@';}
		}
		// GeneIDの取得、nucleotideIDで取れないときは名前で取る
		$m->GeneID						=$ntog[$m->nucleotideID];
		if($m->GeneID==null)$m->GeneID	=$ntog[$m->GeneName];
		// nucleotideIDが埋まってないときは、とりあえずGeneIDかGeneNameから埋める
		if($m->nucleotideID==null)$m->nucleotideID=$gton[$m->GeneID];
		if($m->nucleotideID==null)$m->nucleotideID=$gton[$m->GeneName];
		//
		// 手書きのやつ(Barcode由来)は適当に補完
		if(!preg_match("/(synonymous|stopgain|stoploss|insertion|deletion|---)/is",$m->ExonicFunc)){
			$m->AAchange=$m->ExonicFunc;	// 場所移動
			// nonsynonymousか、frameshiftか、stopgainか適当に判定
			$m->ExonicFunc="nonsynonymous";
			if(preg_match("|X$|is",$m->AAchange))$m->ExonicFunc="stopgain";
			if(preg_match("|fs$|is",$m->AAchange))$m->ExonicFunc="frameshift";
		}	
		//
		$muts[]=$m;
	}	
	return $muts;
}
*/

// get_folderList フォルダ内のファイル一覧を取得
function get_folderList($folder){$dir=opendir($folder);while($filename=readdir($dir))if($filename!=="."&&$filename!=="..")$filenames[]=$filename;closedir($dir);return $filenames;}

// convert_GOID 数字部分だけのGOID(NCBIが返してくる)をGO:XXXXXXXという形にする
function convert_GOID($GOID){$id=$GOID+10000000;$id=substr(strval($id),1,7);return "GO:".$id;}

// number_score関連
function calcPDF($lambda,$obs){$pdf=pow($lambda,$obs)*exp(-$lambda);while($obs>0){$pdf/=$obs;$obs--;}return $pdf;}
function calcp($lambda,$obs){$p=1;for($i=0;$i<$obs;$i++)$p-=calcPDF($lambda,$i);return $p;}
function sortbylength($a,$b){return strlen($b)-strlen($a);}
function number_score($mutrate,$length,$samples,$count){$prob=$length*$samples*$mutrate/1000000;return calcp($prob,$count);}

// file_get_lines ファイルを開き、行まで展開して返す（・ｗ・）手抜き用
function file_get_lines($filename,$headerlines=0){
	$file=file_get_contents($filename);
	$lines=explode(PHP_EOL,rtrim($file));
	if(count($lines)==1&&preg_match("|\n|",$lines[0]))$lines=explode("\n",rtrim($file));
	for($i=0;$i<$headerlines;$i++)array_shift($lines);
	return $lines;
}
//function file_get_lines_cache($filename,$folder){$file=$folder==null?file_get_contents_cache($filename):file_get_contents_cache($filename,$folder);$lines=explode(PHP_EOL,rtrim($file));return $lines;}

// get_geneinfo テスト中
function get_geneinfo($filename="cds.dat"){$geneinfo=null;foreach(file_get_lines($filename) as $line){
	list($nucleotideID,$GeneName,$GeneID,$length,$cds)=explode("\t",rtrim($line));
	if($nucleotideID!=null)$geneinfo[$nucleotideID]->length=$length;
	if($GeneName!=null)$geneinfo[$GeneName]->length=$length;
	if($GeneID!=null)$geneinfo[$GeneID]->length=$length;
}return $geneinfo;}

// get_ntog rna.dat開いてnucleotideID or GeneName → GeneIDの変換リスト作る（・ｗ・）大事、逆もしかり
function get_ntog($filename="rna.dat"){foreach(file_get_lines($filename) as $line){
	list($nucleotideID,$GeneName,$GeneID,$length)=explode("\t",rtrim($line));
	if($nucleotideID!=null)$ntog[$nucleotideID]=$GeneID;
	if($GeneName!=null)$ntog[$GeneName]=$GeneID;
}return $ntog;}
function get_gton($filename="rna.dat"){foreach(file_get_lines($filename) as $line){
	list($nucleotideID,$GeneName,$GeneID,$length)=explode("\t",rtrim($line));
	if($GeneID!=null)$gton[$GeneID]=$nucleotideID;
	if($GeneName!=null)$gton[$GeneName]=$nucleotideID;
}return $gton;}

// get_cdsna cds.datを開いてnucleotideID or GeneName → cds(大文字ATGC)の変換リスト作る
function get_cdsna($filename="cds.dat"){foreach(file_get_lines($filename) as $line){
	list($nucleotideID,$GeneName,$GeneID,$length,$cds)=explode("\t",rtrim($line));
	if($nucleotideID!=null)$cdsna[$nucleotideID]=$cds;
	if($GeneName!=null)$cdsna[$GeneName]=$cds;
}return $cdsna;}

/*
// get_pile_single 1塩基分のpileを返す（・ｗ・）構造体として。
// 2012/02/15 最強版（・ｗ・）９
// ↑単に、mpileupした方が強力なので、最強ではない（・ｗ||||
function get_pile_single($filename,$chr,$pos,$qual=0){
	exec("samtools view -b -h ".$filename." ".$chr.":".$pos."-".$pos." > temp.bam");
	//$lines=null;exec("samtools mpileup -q $qual temp.bam",$lines);
	$lines=null;exec("samtools mpileup temp.bam",$lines);
	//echo "get_pile_single1".PHP_EOL;
	foreach($lines as $line){
		if(preg_match("|^".$chr."\t".$pos."|is",$line)){
			$p=null;list($p->chr,$p->pos,$p->N,$p->total,$p->pile)=explode("\t",rtrim($line));
			return $p;
		}
	}
}
*/


function get_pile_single2($filename,$chr,$pos,$qual=0){
	exec("samtools view -b -h ".$filename." ".$chr.":".$pos."-".$pos." > temp.bam");
	//$lines=null;exec("samtools mpileup -q $qual temp.bam",$lines);
	$lines=null;exec("samtools mpileup temp.bam",$lines);
	//echo "get_pile_single1".PHP_EOL;
	foreach($lines as $line){
		if(preg_match("|^".$chr."\t".$pos."|is",$line)){
			$p=null;list($p->chr,$p->pos,$p->N,$p->total,$p->pile)=explode("\t",rtrim($line));
			return $p;
		}
	}
}

// get_pile_single 1塩基分のpileを返す（・ｗ・）構造体として。
function get_pile_single($filename,$chr,$pos,$param=null,$samtools="samtools"){
	$output=null;exec($samtools." mpileup ".$param." -r ".$chr.":".$pos."-".$pos." ".$filename,$output);
	//echo "get_pile_single2".PHP_EOL;
	
	$p=null;list($p->chr,$p->pos,$p->N,$p->total,$p->pile)=explode("\t",rtrim($output[0]));

	return $p;
}

// count_pile pile内の成分を数える（・ｗ・）ノてきとう。
function count_pile($pile,$ref="",$obs=""){
	$ret=array();
	// insertion,deletion数える
	$ret["-"]=substr_count($pile,"*");
	$ret["+"]=substr_count($pile,"+");
	// insertionのないpile作成、大文字に
	preg_match_all("|\+(\d+)|is",$pile,$matches);
	$array=array_unique($matches[1]);
	foreach($array as $match){$no=$match;$pile=preg_replace("/\+".$no."(A|T|G|C){".$no."}/is","",$pile);}
	$pile=strtoupper($pile);
	// 数える
	$ret["A"]=substr_count($pile,"A");
	$ret["C"]=substr_count($pile,"C");
	$ret["G"]=substr_count($pile,"G");
	$ret["T"]=substr_count($pile,"T");
	// ref,obsが与えられているときは、それを求める
	if($ref!=null&&$obs!=null){
		if($ref==="-"){			// insertion
			$ret["ref"]=max($ret["A"],$ret["C"],$ret["G"],$ret["T"]);
			$ret["obs"]=$ret["+"];
		}else if($obs==="-"){	// deletion
			$ret["ref"]=$ret[substr($ref,0,1)];
			$ret["obs"]=$ret["-"];
		}else{					// SNV
			$ret["ref"]=$ret[$ref];
			$ret["obs"]=$ret[$obs];
		}
	}
	//
	return $ret;
}

function get_BioGRIDDB($filename){
	$lines=file_get_lines($filename,1);
	$DB=null;
	foreach($lines as $line){
		list($temp,$GeneID1,$GeneID2)=explode("\t",$line);
		if($GeneID1!=7316&&$GeneID2!=7316&&$GeneID1!=$GeneID2){	// ユビキチン除去＋自分除去（・ｗ・）なぜかある
			$DB[$GeneID1]->ints[$GeneID2]=1;
			$DB[$GeneID2]->ints[$GeneID1]=1;
		}
	}
	return $DB;
}

function stats_standard_deviation($ary) {
	$avg = array_sum($ary)/count($ary);
	$diff_ary = array();
	foreach ($ary as $val) {$diff = $val-$avg;$diff_ary[] = pow($diff,2);}
	$diff_total = array_sum($diff_ary);
	$diff_avg   = $diff_total/count($diff_ary);
	$stdev = sqrt($diff_avg);
	return $stdev;
}

?>
