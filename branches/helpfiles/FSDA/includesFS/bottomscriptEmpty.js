// this script manages the bottom bar behavior
// please note that: 	'bottomscript.js'   should be put in folder (...)\helpfiles\FSDA\includeFS
// 				   		'bottom.html'       should be put in folder (...)\helpfiles\FSDA\includeFS

var fname;

 
// get the name of the current PARENT file that include bottom.html removing the path with a regular expression
fname=document.location.pathname.match(/[^\/]+$/)[0];

// define new array that will be populated by publishfunctionAlpha.m MATLAB function
var fileArray = new Array( LIST_OF_FILES );

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


