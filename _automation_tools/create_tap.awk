BEGIN {
    counter = 0;	
}

{    
    tapFile = "results.tap";
    
    if (index($0, "Execution of") > 0) {		  
        s = substr($0, index($0, "EXAMPLES_test"));
        n = split(s, records, " ");
        testName = substr(records[1], 15);		  
        print "Parsing: " testName;
        counter++;
	    
        if (index($0, "completed successfully")) {	    
            print "ok " counter " - " testName " completed successfully " >>tapFile;
        } else if (index($0, "FAILED") > 0) {
            errorMessage = substr($0, index($0, "FAILED") + 7);
            print "not ok " counter " - " testName " " errorMessage >>tapFile;
            print "  ---" >>tapFile;
            print "  message: " "\x22" errorMessage "\x22" >>tapFile;
            print "  severity: fail" >>tapFile;
            print "  ..." >>tapFile;
        } else {
            print "not ok " counter " - Test result not recognized" >>tapFile;
        }		
	  }
}

END {
}

