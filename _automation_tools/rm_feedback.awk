BEGIN {
    trovato=0;
}

/<table border="0" cellpadding="0" cellspacing="0" class="feedbacklink/ {
	trovato=1;
	getline
}



/<\/table>/ {

	if (trovato==1) {
		trovato=0;
		getline
	}

}



{
	if (trovato==0) {
		print $0 
	}   
}

