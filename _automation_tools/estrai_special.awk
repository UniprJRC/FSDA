BEGIN {
	skipping_example=0;
	
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
    if (nomefile != null) {
        print "clear all;" >> nomefile;
        print "close all;" >> nomefile;    	
    }
    
    if (classFileName != null) {
        print "end % methods" >> classFileName;
        print "end % class" >> classFileName;    	
    }    	
}

/^%%/ { 	
	if ((getline tmp) > 0) {
	  if (index(tmp, "Interactive_example") > 0 || index(tmp, "example_producing_error") > 0 || index(tmp, "examples_regression") > 0) {
		    print tmp " - Example with special tag found. Skipping...";
		    num++;
		    skipping_example = 1;
		} else {
      skipping_example = 0;
			num++;
			nomefile=sprintf ("EXAMPLES_test/%s_ex%02d.m",nome[1],num) ;
			print nomefile ;
			print "%% Extracted from " tok[p] > nomefile ;	
			print tmp >> nomefile;

			print "function test_" nome[1] "_ex" num "(testCase)" >> classFileName;	
		}
	}
	
	getline 
}

{
	if (nomefile != null) {
		if (skipping_example) {
			getline
		} else {
			print $0 >> nomefile;
		}
	}
	
  if (classFileName != null) {
		if (skipping_example) {
			getline
	  } else {
  	  print $0 >> classFileName;
	  }		
  }		
}

