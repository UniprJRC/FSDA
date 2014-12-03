BEGIN {
    counter = 0;	
}

{    
	if (index($0, "Execution of") > 0) {		  
 		  s = substr($0, index($0, "EXAMPLES_test"));
 		  n = split(s, records, " ");
 		  testName = print substr(records[1], 15);		  
 		  print "Parsing: " testName;
	    counter++;
        if (index($0, "completed successfully")) {	    
	        print "ok " counter " - " testName " " completed successfully " >>"pippo.tap"
	    } else if (index($0, "FAILED") > 0) {
            errorMessage = substr($0, index($0, "FAILED") + 8);
	        print "not ok " counter " - " testName " " errorMessage >>"pippo.tap"
	    } else {
	        print "not ok " counter " - Test result not recognized" >>"pippo.tap";
	    }		
	}
}

END {
}

