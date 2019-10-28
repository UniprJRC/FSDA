window.JST = window.JST || {};

function getProjectDoc() {
    var params = parseParams();
    var services = {
        "messagechannel":"projectdoc"
    }
    requestHelpService(params, services, addProjectDoc);        
}

function addProjectDoc(data) {
    switch (data.pagetype) {
        case "list" : 
            addListPage(data); 
            return;
        case "content" : 
            addContentPage(data); 
            return;
        default : ;
    }
}

function addListPage(data) {
    if (data.contenttype == "landing") {
        addLandingPage(data); 
    } else {
        addContentTypeListPage(data.contenttype, data); 
    }
}

function addLandingPage(data) {
    $("#pd_landing_container").show();
    addLandingPageHeading(data.projectName);
    var listData = data.listData;
    // listData is a structure containing data for each content type (function, example, concept)
    for (var i = 0; i < listData.length; i++) {
        var contentTypeData = listData[i]
        var sectionContentType = contentTypeData.contentType;
        var sectionListData = contentTypeData.listData;    
        if (sectionListData.length > 0) {            
            var h2Text = getContentTypeLabel(sectionContentType);
            addLandingPageSectionHeading(sectionContentType, h2Text);
            addLandingPageSectionContent(data.projectPath, data.projectName, data.projectContentTypes, sectionListData, sectionContentType);   
        }
    }
    addLeftNavHtml(data.projectName, data.projectContentTypes);    
}

function addLandingPageHeading(h1Text) {
    $("#pd_landing_container h1").html(h1Text);
    $("title").html(h1Text);
}

function addLandingPageSectionHeading(contentType, h2Text) {
    $("#pd_landing_container_" + contentType).show();
    $("#pd_landing_container_" + contentType + " h2").html(h2Text);
}

function addLandingPageSectionContent(projectPath, projectName, projectContentTypes, listData, contentType) {
    var contentHtml = createListPageContentHtml(projectPath, projectName, projectContentTypes, listData, contentType);
    $("#pd_landing_content_" + contentType).append(contentHtml);
}

function addContentTypeListPage(contentType, data) {
    $("#pd_list_container").show();
    var listData = data.listData;
    // listData is a structure containing data for each content type (function, example, concept)
    for (var i = 0; i < listData.length; i++) {
        var contentTypeData = listData[i]
        var sectionContentType = contentTypeData.contentType;
        var sectionListData = contentTypeData.listData;
        if ((sectionContentType == contentType) && (sectionListData.length > 0)) {            
            var h1Text = getListPageH1Text(contentType, data.projectName);
            addListPageHeading(h1Text);
            addListPageContent(data.projectPath, data.projectName, data.projectContentTypes, sectionListData, contentType);   
        }
    }
    addLeftNavHtml(data.projectName, data.projectContentTypes);    
}

function getListPageH1Text(contentType, projectName) {
    var headingLabel = getContentTypeLabel(contentType);
    var h1Text = (projectName && projectName.length > 0) ? projectName + " &mdash; " + headingLabel : headingLabel;
    return h1Text;
}

function getContentTypeLabel(contentType) {
    switch (contentType) {
        case "concept" : 
            return getLocalizedString("project_doc_concepts");
        case "example" : 
            return getLocalizedString("topnav_examples");
        case "function" : 
            return getLocalizedString("project_doc_functions");
        case "landing" : 
            return "";
        default :
            return "";
    }
}

function addListPageHeading(h1Text) {
    $("#pd_list_container h1").html(h1Text);
    $("title").html(h1Text);
}

function addListPageContent(projectPath, projectName, projectContentTypes, listData, contentType) {
    var contentHtml = createListPageContentHtml(projectPath, projectName, projectContentTypes, listData, contentType);
    $("#pd_list_content").append(contentHtml);
}

function createListPageContentHtml(projectPath, projectName, projectContentTypes, listData, contentType) {
    // Intentionally hard coding the pagetype to be content here. The action will be to display the content page.
    var projParams = {};
    projParams.projectPath = projectPath;
    projParams.projectName = projectName;

    var contentTypesString = '';
    for (var i = 0; i < projectContentTypes.length; i++) {
        var ctString = projectContentTypes[i];
        contentTypesString = contentTypesString + ctString + ',';
    }

    projParams.projectContentTypes = contentTypesString;
    projParams.pagetype = "content";
    projParams.contenttype = contentType;
    return JST['list_page_content_item']({"projParams" : projParams, "listData" : listData});
}

JST['list_page_content_item'] = _.template(
  '<div>' +
    '<ul>' +
    '<% _.each(listData, function(listItem) { %>' +
        '<% projParams.filePath = listItem.data.filePath; %>' +
        '<% var qs = $.param(projParams); %>' +
        '<li style="display: list-item;"> <a href="?<%= qs %>"> <%= listItem.name %> </a> </li>' +
    '<% }); %>' +
    '</ul>' +
  '</div>'
);

function addLeftNavHtml(projectName, projectContentTypes) {
    var leftNavHtml = createLeftNavHtml(projectName, projectContentTypes);
    $("#pd_left_nav").show();
    $("#pd_left_nav").append(leftNavHtml);
}

function createLeftNavHtml(projectName, projectContentTypes) {
    // Intentionally hard coding the pagetype to be list here. The action is to display the list page.
    return JST['left_nav']({"projectName" : projectName, "pagetype" : "list", "projectContentTypes" : projectContentTypes});
}

JST['left_nav'] = _.template(
  '<div>' +
    '<li style="display: list-item;"> <a href="?pagetype=list&contenttype=landing"> <%= projectName %> </a> </li>' +
    '<ul>' +
    '<% _.each(projectContentTypes, function(contentType) { %>' +
        '<% var projParams = {}; %>' +
        '<% projParams.pagetype = pagetype; %>' +
        '<% projParams.contenttype = contentType; %>' +
        '<% var qs = $.param(projParams); %>' +
        '<li style="display: list-item;"> <a href="?<%= qs %>"> <%= getContentTypeLabel(contentType) %> </a> </li>' +
    '<% }); %>' +
    '</ul>' +
  '</div>'
);

function addContentPage(data) {
    $(window).on('resize', resizeContentFrame);
    $("#pd_content_container").show();
    addContentPageContent(data.contenttype, data.html);
    addLeftNavHtml(data.projectName, data.projectContentTypes);    
}

function addContentPageContent(contentType, html) {
    $("#pd_frame").contents().find('html').html(html);
    resizeContentFrame();
}

function resizeContentFrame() {
    var iframe = document.getElementById("pd_frame");
    if (iframe) {
        var h = window.innerHeight;
        // The width is to wide. Not accounting for the left nav width.
        var w = window.innerWidth;
        iframe.style.height = h + "px";
        iframe.style.width = w + "px";    
    }
}
