// this script manages the bottom bar behavior
// please note that: 	'bottomscript.js'   should be put in folder (...)\helpfiles\FSDA\includeFS
// 				   		'bottom.html'       should be put in folder (...)\helpfiles\FSDA\includeFS

var fname;

 
// get the name of the current PARENT file that include bottom.html removing the path with a regular expression
fname=document.location.pathname.match(/[^\/]+$/)[0];

// define new array that will be populated by publishfunctionAlpha.m MATLAB function
var fileArray = new Array( "index.html","ace.html","aceplot.html","add2spm.html","add2yX.html","addt.html","avas.html","barnardtest.html","basicPower.html","bc.html","boxplotb.html","boxtest.html","brushFAN.html","brushRES.html","brushROB.html","bwe.html","cabc.html","carbikeplot.html","cascade.html","cdsplot.html","clickableMultiLegend.html","ClusterRelabel.html","combsFS.html","CorAna.html","CorAnaplot.html","corrNominal.html","corrOrdinal.html","covplot.html","CressieRead.html","crosstab2datamatrix.html","ctsub.html","dempk.html","ellipse.html","exactcdf.html","existFS.html","fanplot.html","findDir.html","findFile.html","forecastTS.html","FowlkesMallowsIndex.html","FSM.html","FSMbonfbound.html","FSMbsb.html","FSMeda.html","FSMedaeasy.html","FSMenvmmd.html","FSMfan.html","FSMinvmmd.html","FSMmmd.html","FSMmmdeasy.html","FSMmmdrs.html","FSMmmdWITHRANKCHECK.html","FSMtra.html","FSR.html","FSRaddt.html","FSRB.html","FSRBbsb.html","FSRBeda.html","FSRBmdr.html","FSRbonfbound.html","FSRBr.html","FSRbsb.html","FSRcp.html","FSReda.html","FSRenvmdr.html","FSRfan.html","FSRH.html","FSRHbsb.html","FSRHeda.html","FSRHmdr.html","FSRinvmdr.html","FSRmdr.html","FSRmdrrs.html","FSRms.html","FSRr.html","FSRts.html","FSRtsbsb.html","FSRtsmdr.html","genr8.html","genSigmaGPCM.html","GowerIndex.html","GYfilt.html","HAbdp.html","HAc.html","HAeff.html","HApsi.html","HApsider.html","HApsix.html","HArho.html","HAwei.html","histFS.html","htmlwriteFS.html","HUeff.html","HUpsi.html","HUpsider.html","HUpsix.html","HUrho.html","HUwei.html","HYPbdp.html","HYPc.html","HYPck.html","HYPeff.html","HYPk.html","HYPpsi.html","HYPpsider.html","HYPpsix.html","HYPrho.html","HYPwei.html","inversegamcdf.html","inversegaminv.html","inversegampdf.html","isfunction.html","kdebiv.html","levfwdplot.html","lexunrank.html","lga.html","logmvnpdfFS.html","LTSts.html","LTStsVarSel.html","LXS.html","mahalFS.html","makecontentsfileFS.html","malfwdplot.html","malindexplot.html","mcd.html","mdrplot.html","MixSim.html","MixSimreg.html","mmdplot.html","mmdrsplot.html","MMmult.html","MMmultcore.html","MMmulteda.html","MMreg.html","MMregcore.html","MMregeda.html","mreadFS.html","Mscale.html","mtR.html","mve.html","mveeda.html","nchoosekFS.html","ncpci.html","ncx2mixtcdf.html","normBoxCox.html","normYJ.html","openMatlabFileFromHTML.html","OPTbdp.html","OPTc.html","OPTeff.html","OPTpsi.html","OPTpsider.html","OPTpsix.html","OPTrho.html","OPTwei.html","overlap.html","overlapmap.html","PoolClose.html","PoolPrepare.html","position.html","Powertra.html","publishBibliography.html","publishFS.html","publishFunctionAlpha.html","publishFunctionCate.html","Qn.html","qqplotFS.html","quickselectFS.html","RandIndexFS.html","randsampleFS.html","rcontFS.html","regressB.html","regressH.html","regressHart.html","regressHhar.html","regressts.html","removeExtraSpacesLF.html","repDupValWithMean.html","resfwdplot.html","resfwdplot.html","resindexplot.html","restrdeter.html","restrdeterGPCM.html","restrdeterGPCM.html","restreigen.html","restreigeneasy.html","restrshapeGPCM.html","restrshapeGPCM.html","RKbdp.html","RKeff.html","RKpsi.html","RKpsider.html","RKpsix.html","RKrho.html","RKwei.html","rlga.html","rlsmo.html","RobCov.html","RobRegrSize.html","rthin.html","Score.html","ScoreYJ.html","ScoreYJpn.html","SDest.html","shuffling.html","simdataset.html","simdatasetreg.html","simulateTS.html","smothr.html","Smult.html","Smulteda.html","Sn.html","SparseTableTest.html","spmplot.html","Sreg.html","Sregeda.html","subsets.html","suplabel.html","supsmu.html","tabulateFS.html","Taureg.html","TBbdp.html","TBc.html","TBeff.html","tBothSides.html","TBpsi.html","TBpsider.html","TBpsix.html","TBrho.html","TBwei.html","tclust.html","tclusteda.html","tclustIC.html","tclustICplot.html","tclustICsol.html","tclustreg.html","tclustregIC.html","tkmeans.html","triu2vec.html","unibiv.html","upperfracpos.html","verLessThanFS.html","vervaatrnd.html","vervaatsim.html","vervaatxdf.html","wedgeplot.html","winsor.html","WNChygepdf.html","wraptextFS.html","wthin.html","xmlcreateFS.html","yXplot.html","zscoreFS.html","function-cate.html" );

// create the array index of the current file
var ind=fileArray.indexOf(fname);



if (ind==1){
// substitute the following placeholder tags with the correct ones 
$("a[href='http://www.left.com/']").attr('href', fileArray [ind]);
$("a[href='http://www.right.com/']").attr('href', fileArray [ind+1]);
$("#left").text(fileArray[ind].substr(0,fileArray[ind].length-5));
$("#right").text(fileArray[ind+1].substr(0,fileArray[ind+1].length-5));
}
else if ((ind+1)==(fileArray.length-1)) 
{
// substitute the following placeholder tags with the correct ones 
$("a[href='http://www.left.com/']").attr('href', fileArray [ind-1]);
$("a[href='http://www.right.com/']").attr('href', fileArray [ind]);
$("#left").text(fileArray[ind-1].substr(0,fileArray[ind-1].length-5));
$("#right").text(fileArray[ind].substr(0,fileArray[ind].length-5));
}
// when the file is not included in the list 
else if (ind==-1) 
{
// substitute the following placeholder tags with the correct ones 
$("a[href='http://www.left.com/']").attr('href', 'function-alpha.html');
$("a[href='http://www.right.com/']").attr('href', 'function-cate.html');
$("#left").text('function-alpha.html');
$("#right").text('function-cate.html');
}

else{
// substitute the following placeholder tags with the correct ones 
$("a[href='http://www.left.com/']").attr('href', fileArray [ind-1]);
$("a[href='http://www.right.com/']").attr('href', fileArray [ind+1]);
$("#left").text(fileArray[ind-1].substr(0,fileArray[ind-1].length-5));
$("#right").text(fileArray[ind+1].substr(0,fileArray[ind+1].length-5));
}


