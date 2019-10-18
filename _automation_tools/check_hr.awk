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
	if ( ! match (ref,"www.eio.uva.es"))
	if ( ! match (ref,"doi.org"))
	if ( ! match (ref,"www.r-project.org"))
	 if ( ! match (ref, "mailto"))
		if ( ! match (ref, "^#"))
			if ( ! match (ref, "riani"))
				if ( ! match (ref, "matlab"))
				if ( ! match (ref, "javascript"))
				if ( ! match (ref, "#Example"))
				if ( ! match (ref, "support.sas"))
				if ( ! match (ref, "en.wikipedia"))
				if ( ! match (ref, "left.com"))
				if ( ! match (ref, "right.com"))
				if ( ! match (ref, "ec.europa.eu"))
				if ( ! match (ref, "www"))
				if ( ! match (ref, "rosa.unipr"))
				if ( ! match (ref, "users.ugent.be"))
				if ( ! match (ref, "wis.kuleuven.be"))
				if ( ! match (ref, "r-project.org"))
				
					if ( ! match (ref, "matworks")){
						print tok[2]
					}
#	for (i=1;i<=NF;i++) {
#	print "il campo " i " : " $i 
#	} ;
}

