// this script manages the bottom bar behavior
// plz. note that: 'function-alpha.txt' should be put in folder (...)\helpfiles\FSDA
// 				   'bottomscript.js'   should be put in folder (...)\helpfiles\FSDA\includeFS
// 				   'bottom.html'       should be put in folder (...)\helpfiles\FSDA\includeFS

var fname;

 
// get the name of the current file that include bottom.html removing the path with a regular expression
fname=document.location.pathname.match(/[^\/]+$/)[0];

// load the file list and store the big string in a variable by a callback f(x)='function(stuff)'
//$.get('function-alpha.txt', function(stuff){

// define new blank array and...
var fileArray = new Array( LIST_OF_FILES );

// ...populate it splitting the big string
//fileArray = stuff.split(",");
	//fileArray = ["index.html" ,"add2spm.html","FSRaddt.html","FSRB.html","FSRBbsb.html","FSRBeda.html"];
	

// array index of the current file
//var ind=$.inArray(fname, fileArray );

var ind=fileArray.indexOf(fname);


// alert(fileArray.length);
//alert(ind);


if (ind==1){
// substitute the following placeholder tags with the correct ones 
$("a[href='http://www.topleft.com/']").attr('href', fileArray [ind]);
$("a[href='http://www.topright.com/']").attr('href', fileArray [ind+1]);
$("#topleft").text(fileArray[ind].substr(0,fileArray[ind].length-5));
$("#topright").text(fileArray[ind+1].substr(0,fileArray[ind+1].length-5));
}
else if ((ind+1)==(fileArray.length-1)) 
{
// substitute the following placeholder tags with the correct ones 
$("a[href='http://www.topleft.com/']").attr('href', fileArray [ind-1]);
$("a[href='http://www.topright.com/']").attr('href', fileArray [ind]);
$("#topleft").text(fileArray[ind-1].substr(0,fileArray[ind-1].length-5));
$("#topright").text(fileArray[ind].substr(0,fileArray[ind].length-5));
}
// when the file is not included in the list 
else if (ind==-1) 
{
// substitute the following placeholder tags with the correct ones 
$("a[href='http://www.topleft.com/']").attr('href', 'function-alpha.html');
$("a[href='http://www.topright.com/']").attr('href', 'function-cate.html');
$("#topleft").text('function-alpha.html');
$("#topright").text('function-cate.html');
}

else{
// substitute the following placeholder tags with the correct ones 
$("a[href='http://www.topleft.com/']").attr('href', fileArray [ind-1]);
$("a[href='http://www.topright.com/']").attr('href', fileArray [ind+1]);
$("#topleft").text(fileArray[ind-1].substr(0,fileArray[ind-1].length-5));
$("#topright").text(fileArray[ind+1].substr(0,fileArray[ind+1].length-5));
}


