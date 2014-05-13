BEGIN {
    primo=0;
    FS="href=" ;
}


/<a href="/ {
 	if (primo == 0) {
	    print FILENAME
	};
	primo++;
#	print $0;
	split ($2, tok,"\"");
	ref=tok[2];
	if ( ! match (ref, "mailto"))
		if ( ! match (ref, "^#"))
			if ( ! match (ref, "riani"))
				if ( ! match (ref, "matlab"))
					if ( ! match (ref, "matworks")){
						print tok[2]
					}
#	for (i=1;i<=NF;i++) {
#	print "il campo " i " : " $i 
#	} ;
}

