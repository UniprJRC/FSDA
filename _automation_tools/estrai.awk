BEGIN {
    trovato=0;
    num=0;
    p=split( ARGV[ARGC-1],tok,"/" );
    q=split( tok[p],nome,"." )
	
	/* BEGIN_TEST_CLASS */
	classFileName = sprintf("EXAMPLES_classes/Test_%s.m", nome[1]);	
	print "%% Extracted from " tok[p] > classFileName;
	print "classdef Test_" nome[1] " < matlab.unittest.TestCase" >> classFileName;
	print "methods(Test)" >> classFileName;
	/* END_TEST_CLASS */
}

END {
    /* BEGIN_TEST_CLASS */
    if (nomefile != "") {
        print "clear all;" >> nomefile;
        print "close all;" >> nomefile;    	
    }
    
    if (classFileName != null) {
        print "end % methods" >> classFileName;
        print "end % class" >> classFileName;    	
    }
    	
	/* END_TEST_CLASS */
}

/^%{/ {

	if ((getline tmp) > 0) {
	    if (index(tmp, "Interactive_example") > 0 || index(tmp, "example_producing_error") > 0 || index(tmp, "Skip_example")) {
		    print tmp " - Example with special tag found. Skipping...";
		} else {
			trovato=1;
			num++;
			nomefile=sprintf ("EXAMPLES_test/%s_ex%02d.m",nome[1],num) ;
			print nomefile ;
			print "%% Extracted from " tok[p] > nomefile ;	
			print tmp >> nomefile;

			/* BEGIN_TEST_CLASS */
			print "function test_" nome[1] "_ex" num "(testCase)" >> classFileName;	
			/* END_TEST_CLASS */		
		}
	}
	
	getline 
}

/^%}/ {
	trovato=0
	
	/* BEGIN_TEST_CLASS */
	print "end" >> classFileName;
	/* END_TEST_CLASS */
}

{
	if (trovato==1) {
		print $0 >> nomefile
	
    /* BEGIN_TEST_CLASS */	
		print $0 >> classFileName;
		/* END_TEST_CLASS */
	}   
}

