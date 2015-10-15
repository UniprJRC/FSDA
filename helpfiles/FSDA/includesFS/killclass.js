// jq162 = jQuery.noConflict( true );
//jQuery(window).load(function($) {
//jQuery(document).ready(function($) {

// con questa sintassi, la funzione  viene eseguita una volta sola
jQuery(document).one('ready',function(){

// fa rientrare la TOC xml
window.parent.$('.row-offcanvas').toggleClass('active');

 // alert($('#3ptoc').text());
 
 // toglie il tasto sopra la TOC xml con close
 window.parent.$('#nav_toggle').detach();

// elimina completamente la class="toc_container toc_expanded" della TOC xml 
window.parent.$('#3ptoc').detach();


// window.parent.$('#doc_center_content').replaceWith('<iframe id="doc_center_content" onload="updateContentSize();" src="file:///E:/FSDA/helpfiles/FSDA/tclusttmp3.html?3pcontent=true" style="height: 76px;" allowTransparency="true" scrolling="yes" frameborder="1"></iframe>' );
// });

//<script src="includesFS/jquery-latest.js" type="text/javascript"></script>

//window.parent.$('#doc_center_content').replaceWith('<section id="doc_center_content" lang="en"> <div class="function_ref"> <object type="text/html" data="tclusttmp3.html"></object></div></section>');
 });

