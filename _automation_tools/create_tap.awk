BEGIN {
    counter = 0;	
}

{    
	if (index($0, "Execution of") > 0) {
		  n = split($0, records, "\\");
 		  s = substr($0, index($0, "EXAMPLES_test"));
 		  print substr(s, 14, index(s, " "));		  
	    counter++;
        if (index($0, "completed successfully")) {	    
	        print "ok " counter " - completed successfully " >>"pippo.tap"
	    } else if (index($0, "FAILED") > 0) {
            errorMessage = substr($0, index($0, "FAILED") + 8);
	        print "not ok " counter " - " errorMessage >>"pippo.tap"
	    } else {
	        print "not ok " counter " - Test result not recognized" >>"pippo.tap";
	    }		
	}
}

END {
}

