window.JST = window.JST || {};

function getReleaseNotesList() {
    var params = {};
    var services = {
        "messagechannel":"prodfilter",
        "requesthandler":"productfilter:handleSelectedProducts",
        "webservice":getReleaseNotesListWebServiceUrl()
    }
    requestHelpService(params, services, function(data) {
        var prodNavList = data.prodnavlist;
        if (typeof prodNavList === "string") {
            prodNavList = $.parseJSON(prodNavList);
        }
        handleReleaseNotes(prodNavList);
  	});
}

function getReleaseNotesListWebServiceUrl() {
    var release = getDocReleaseFromSearchBox();
    if (typeof getDocRelease === 'function') {
        release = getDocRelease();
    }
    return "/help/search/prodfilter/doccenter/en/" + release;
}

function handleSelectedProducts(prodList, prodNavList) {
    handleReleaseNotes(prodNavList);
}

function handleReleaseNotes(prodNavList) {
    var docRoot = findDocRoot();
    addLeftNav(prodNavList, docRoot);
    populateDropDown(prodNavList, docRoot);
}

function addLeftNav(prodNavList, docRoot) {
    var leftHtml = createLeftNavHtml(prodNavList, docRoot);
    $("#left_nav_categories").html(leftHtml);
    addSmoothScroll();
}

function createLeftNavHtml(prodNavList, docRoot) {
    return JST['leftnav_item']({"prodNavList" : prodNavList, "docRoot" : docRoot});
}

JST['leftnav_item'] = _.template(
  '<div class="search_refine">' +
    '<h3><%= getLocalizedString("topnav_header_category") %></h3>' +
    '<div id="nav_categories" style="max-height:350px;overflow-y:auto;">' +
    '<ul class="nav_toc" id="nav_siblings">' +
    '<% _.each(prodNavList, function(prod) { %>' +
      '<li style="display: list-item;"><a href="<%= docRoot %>/<%= getReleaseNotesPage(prod.helplocation) %>"><%= prod.displayname %></a></li>' +
    '<% }); %>' +
    '</ul>' +
    '</div>' +
  '</div>'
);

function populateDropDown(prodNavList, docRoot) {
    var dropDownHtml = createDropDownHtml(prodNavList, docRoot);
    $("#grn_to_prodrn_dropdown").html(dropDownHtml);    
}

function createDropDownHtml(prodNavList, docRoot) {
    return JST['dropdown_item']({"prodNavList" : prodNavList, "docRoot" : docRoot});
}

JST['dropdown_item'] = _.template(
    '<% _.each(prodNavList, function(prod) { %>' +
      '<li style="display: list-item;"><a href="<%= docRoot %>/<%= getReleaseNotesPage(prod.helplocation) %>"><%= prod.displayname %></a></li>' +
    '<% }); %>'
);

function getReleaseNotesPage(page) {
    return page.replace(/(index.html)$/, 'release-notes.html');
}