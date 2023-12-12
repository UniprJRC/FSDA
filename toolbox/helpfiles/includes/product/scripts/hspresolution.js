//===============================================================================
////==================HSP link templates=========================================
window.JST = window.JST || {};
JST['installedHspTmpl'] = _.template(
       '<% for (var i = 0; i < installedhsps.length; i++) { %>' +
       '<% var hsp = installedhsps[i]; %>' +
       '<li style="display: list-item;"><a href="<%= hsp.helplocation %>"><%= hsp.displayname %></a></li>' +
       '<% } %>'
);

//===============================================================================