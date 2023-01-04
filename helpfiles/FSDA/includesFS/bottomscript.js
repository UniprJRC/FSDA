// this script manages the bottom bar behavior
// please note that: 	'bottomscript.js'   should be put in folder (...)\helpfiles\FSDA\includeFS
// 				   		'bottom.html'       should be put in folder (...)\helpfiles\FSDA\includeFS

var fname;

 
// get the name of the current PARENT file that include bottom.html removing the path with a regular expression
fname=document.location.pathname.match(/[^\/]+$/)[0];

// define new array that will be populated by publishfunctionAlpha.m MATLAB function
var fileArray = new Array( "index.html","ace.html","aceplot.html","add2boxplot.html","add2spm.html","add2yX.html","addt.html","avas.html","avasms.html","balloonplot.html","barnardtest.html","barVariableWidth.html","basicPower.html","bc.html","biplotFS.html","boxcoxR.html","boxplotb.html","boxtest.html","brushFAN.html","brushRES.html","brushROB.html","bwe.html","cabc.html","carbikeplot.html","carbikeplotGPCM.html","cascade.html","cdsplot.html","CEVmodel.html","clickableMultiLegend.html","ClusterRelabel.html","combsFS.html","Copy_of_ctlcurves.html","CorAna.html","CorAnaplot.html","corrcdf.html","corrNominal.html","corrOrdinal.html","corrpdf.html","covplot.html","CressieRead.html","crosstab2datamatrix.html","ctlcurves.html","ctlcurvesplot.html","ctsub.html","dempk.html","distribspec.html","ellipse.html","estregimeTAR.html","exactcdf.html","existFS.html","fanBIC.html","fanBICpn.html","fanplot.html","FE_int_vol.html","FE_int_vol_Fejer.html","FE_spot_vol.html","FE_spot_vol_FFT.html","findDir.html","findFile.html","forecastTS.html","FowlkesMallowsIndex.html","FSCorAnaeda.html","FSCorAnaenvmmd.html","FSM.html","FSMbonfbound.html","FSMbsb.html","FSMeda.html","FSMedaeasy.html","FSMenvmmd.html","FSMfan.html","FSMinvmmd.html","FSMmmd.html","FSMmmdeasy.html","FSMmmdrs.html","FSMtra.html","FSR.html","FSRaddt.html","FSRB.html","FSRBbsb.html","FSRBeda.html","FSRBmdr.html","FSRbonfbound.html","FSRBr.html","FSRbsb.html","FSRcp.html","FSReda.html","FSRenvmdr.html","FSRfan.html","FSRH.html","FSRHbsb.html","FSRHeda.html","FSRHmdr.html","FSRinvmdr.html","FSRmdr.html","FSRmdrrs.html","FSRms.html","FSRr.html","FSRts.html","FSRtsbsb.html","FSRtsmdr.html","funnelchart.html","genSigmaGPCM.html","GowerIndex.html","GUIconcentration.html","GUIcorr.html","GUIcov.html","GUImad.html","GUIpowermean.html","GUIquantile.html","GUIregress.html","GUIskewness.html","GUIstd.html","GUItrimmean.html","GUIvar.html","GYfilt.html","HAbdp.html","HAc.html","HAeff.html","HApsi.html","HApsider.html","HApsix.html","HArho.html","HAwei.html","histFS.html","htmlwriteFS.html","HUc.html","HUeff.html","HUpsi.html","HUpsider.html","HUpsix.html","HUrho.html","HUwei.html","HYPbdp.html","HYPc.html","HYPck.html","HYPeff.html","HYPk.html","HYPpsi.html","HYPpsider.html","HYPpsix.html","HYPrho.html","HYPwei.html","inversegamcdf.html","inversegaminv.html","inversegampdf.html","isfunction.html","kdebiv.html","levfwdplot.html","lexunrank.html","lga.html","logmvnpdfFS.html","LTSts.html","LTStsLSmult.html","LTStsVarSel.html","LXS.html","mahalCorAna.html","mahalFS.html","makecontentsfileFS.html","malfwdplot.html","malindexplot.html","mcd.html","mcdCorAna.html","mcdeda.html","mdpattern.html","mdpd.html","mdpdR.html","mdpdReda.html","mdrplot.html","mdrrsplot.html","medcouple.html","MixSim.html","MixSimreg.html","Mlocation.html","Mlocsca.html","mmdplot.html","mmdrsplot.html","MMmult.html","MMmultcore.html","MMmulteda.html","MMreg.html","MMregcore.html","MMregeda.html","moonplot.html","mreadFS.html","Mscale.html","mtR.html","mve.html","mveeda.html","nchoosekFS.html","ncpci.html","ncx2mixtcdf.html","normBoxCox.html","normYJ.html","normYJpn.html","openMatlabFileFromHTML.html","OPTbdp.html","OPTc.html","OPTeff.html","OptimalCuttingFrequency.html","OPTpsi.html","OPTpsider.html","OPTpsix.html","OPTrho.html","OPTwei.html","overlap.html","overlapmap.html","pcaFS.html","PDbdp.html","PDc.html","PDeff.html","PDpsi.html","PDpsider.html","PDpsix.html","PDrho.html","PDwei.html","position.html","Powertra.html","publishBibliography.html","publishFS.html","publishFunctionAlpha.html","publishFunctionCate.html","Qn.html","qqplotFS.html","quickselectFS.html","quickselectFS_demo.html","quickselectFSw.html","quickselectFSw_demo.html","RandIndexFS.html","randsampleFS.html","rcontFS.html","regressB.html","regressH.html","regressHart.html","regressHhar.html","regressts.html","removeExtraSpacesLF.html","repDupValWithMean.html","resfwdplot.html","resindexplot.html","restrdeter.html","restrdeterGPCM.html","restreigen.html","restreigeneasy.html","restreigenmemopt.html","restrshapeExact.html","restrshapeGPCM.html","restrSigmaGPCM.html","RhoPsiWei.html","RKbdp.html","RKeff.html","RKpsi.html","RKpsider.html","RKpsix.html","RKrho.html","RKwei.html","rlga.html","rlsmo.html","RobCov.html","RobRegrSize.html","rthin.html","scatterboxplot.html","Score.html","ScoreYJ.html","ScoreYJall.html","ScoreYJmle.html","ScoreYJpn.html","SDest.html","SETARX.html","setdiffFS.html","setToolboxStartEnd.html","shuffling.html","simdataset.html","simdatasetreg.html","simulateLM.html","simulateTS.html","smothr.html","Smult.html","Smulteda.html","Sn.html","SparseTableTest.html","spmplot.html","Sreg.html","Sregeda.html","subsets.html","suplabel.html","supsmu.html","tabulateFS.html","Taureg.html","TBbdp.html","TBc.html","TBeff.html","tBothSides.html","TBpsi.html","TBpsider.html","TBpsix.html","TBrho.html","TBwei.html","tclust.html","tclusteda.html","tclustIC.html","tclustICgpcm.html","tclustICplot.html","tclustICplotGPCM.html","tclustICsol.html","tclustICsolGPCM.html","tclustreg.html","tclustregeda.html","tclustregIC.html","tkmeans.html","triu2vec.html","twdcdf.html","twdpdf.html","twdrnd.html","unibiv.html","upperfracpos.html","verLessThanFS.html","vervaatrnd.html","vervaatsim.html","vervaatxdf.html","VIOM.html","waterfallchart.html","wedgeplot.html","winsor.html","WNChygecdf.html","WNChygepdf.html","wraptextFS.html","wthin.html","xmlcreateFS.html","yXplot.html","zscoreFS.html","function-cate.html" );

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


