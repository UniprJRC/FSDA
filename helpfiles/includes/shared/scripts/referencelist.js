window.JST = window.JST || {};

function getReferenceList(product) {
    var services = {
        "messagechannel":"reflist",
        "requesthandler":"reflist:getReferenceItems",
        "webservice":getRefListWebServiceUrl()
    }
    var params = parseParams();
    if (product) {
        params.product = product;
    }
    params.docroot = findDocRoot();
    if (! params.type) {
        params["type"] = "function";
    }
    
    selectTopNavLink(params.type);
    requestHelpService(params, services, function(data) {
        addRefItems(data);
    });


    var topNavParams = {
        "category":params.category
    }
    if (params.capability) topNavParams.capability = params.capability;
    if (params.listtype) topNavParams.listtype = params.listtype;
    
    if (product) {
        updateTopNavLinks(topNavParams, getCruxOverrides("gn", undefined, topNavParams["category"]));
    }
    
    if (params.type) {
        var typeLabel = getLabelForType(params.type);
        if (typeLabel) {
        var allParams = {};
        $.extend(allParams, params);
        var cruxParam = buildCruxParam(allParams, getCruxOverrides("lftnav", allParams["type"], undefined));
        var qs = $.param(cruxParam);
        var href = findDocRoot() + '/referencelist.html?type=' + params.type + "&" + qs;
            var allRefItemsLink = $('<a href="' + href + '">' + typeLabel + '</a>');
            $("ul.nav_breadcrumb li a").after(allRefItemsLink);
        }
    }
}

function getRefListWebServiceUrl() {
    var lang = getPageLanguage() || "en";

    var release = getDocReleaseFromSearchBox();
    if (typeof getDocRelease === 'function') {
        release = getDocRelease();
    }

    return "/help/search/reflist/doccenter/" + lang + "/" + release;
}

function createListTypeButtons(params) {
    var categoryLabel = getLocalizedString("topnav_by_category");
    var alphabeticalLabel = getLocalizedString("topnav_alphabetical_list");
    var html = '<div class="col-xs-12"><ul class="list-inline pull-right">';
    if (params.listtype && params.listtype === "alpha") {
        var link = buildTopNavUrl("referencelist.html", {"listtype":"cat"});
        html += '<li><a href="' + link + '">' + categoryLabel + '</a></li>';
        html += '<li class="add_list_separator_left">' + alphabeticalLabel + '</li>';
    } else {
        var link = buildTopNavUrl("referencelist.html", {"listtype":"alpha"});
        html += '<li>' + categoryLabel + '</span>';
        html += '<li class="add_list_separator_left">';
        html += '<a href="' + link + '">' + alphabeticalLabel + '</a>';
    }
    html += '</ul></div>';
    $("#listtypebutton").parent().replaceWith(html);
}

function getCategoryDisplayName(category) {
    if (category["label-key"]) {
        return getLocalizedString(category["label-key"]);
    } else {
        return category.name;
    }
}

JST['category_item'] = _.template(
    '<% if (level == 1) { %>' +
      '<% populateH1Tag(getCategoryDisplayName(category), typeLabel); %>' +
    '<% } else { %>' +
      '<h<%=level%> class="add_clear_both"><%= getCategoryDisplayName(category) %></h<%=level%>><span></span>' +
    '<% } %>' +
  '<% var childCats = category["child-categories"]; %>' +
    '<% var groups = category["grouped-leaf-items"]; %>' +
    '<% if (groups && groups.length > 0) { %>' +
      '<%= createAlphabetLinks(groups) %>' +
      '<% _.each(groups, function(group) { %>' +
        '<% if (group.id.length > 0) { %>' +
          '<h2 id="section_<%= group.id %>" class="anchor_link"><%= group.name %></h2><span></span>' +
        '<% } %>' +
        '<%= tableTemplate({"leafItems":group["leaf-items"],"listEltType":listEltType}) %>' +
      '<% }); %>' +
    '<% } %>' + 
    '<% var leafItems = category["leaf-items"]; %>' +
    '<% if (leafItems && leafItems.length > 0) { %>' +
      '<%= tableTemplate({"leafItems":leafItems,"listEltType":listEltType}) %>' +
    '<% } %>' +
  '<% if (childCats && childCats.length > 0) { %>' +
    '<% _.each(category["child-categories"], function(childCategory) { %>' +
      '<%= categoryTemplate({"category":childCategory,"level":level+1,"categoryTemplate":categoryTemplate,"tableTemplate":tableTemplate,"listEltType":listEltType}) %>' +
    '<% }); %>' +
  '<% } %>'
);


JST['reference_table'] = _.template(
  '<table class="table tablecondensed">' +
    '<% _.each(leafItems, function(leafItem) { %>' +
      '<tr><td class="term"><a href="<%= leafItem.path %>"><<%= listEltType %>><%= leafItem.name %></a></<%= listEltType %>></td>' +
      '<td class="description"><%= leafItem.purpose %></td></tr>' +
    '<% }); %>' +
  '</table>'
);

function createAlphabetLinks(groups) {
    html = '<div style="margin-bottom:10px">';
    var activeLinks = [];
    for (var i = 0; i < groups.length; i++) {
        var id = groups[i].id;
        if (id !== '') {
            activeLinks.push(id);
        }
    }
    
    
    for (var charCode = 'A'.charCodeAt(0); charCode <= 'Z'.charCodeAt(0); charCode++) {
        var letter = String.fromCharCode(charCode);
        if (activeLinks.length > 0 && letter === activeLinks[0]) {
            html += '<a href="#section_' + letter + '" style="margin-right:10px;" class="intrnllnk">' + letter + '</a>';
            activeLinks.shift();
        } else {
            html += '<span style="color:#999999;margin-right:10px;">' + letter + '</span>';
        }
    }
    
    html += "</div>";
    return html;
}

function addRefItems(data) {
    var params = parseParams();
    if (! params.type) {
        params["type"] = "function";
    }
    var typeLabel = getLabelForType(params.type);
    var templateData = {
        "category":data["category"],
        "level":1,
        "typeLabel":typeLabel,
        "categoryTemplate":JST["category_item"],
        "tableTemplate":JST["reference_table"],
        "listEltType":getListItemElement(params.type)
    }
    var cruxOverrides = getCruxOverrides("lftnav", params["type"], data["category"]);
    var categoryHtml;
    if (data.category.leafcount === 0) {
        populateH1Tag(data.category.name, typeLabel);
        categoryHtml = handleNoResults(data.category.name, data.relatedCategories);
    } else {
        if (params.type && params.type !== "app" && data.category.name !== "") {
            createListTypeButtons(params);
        }
        categoryHtml = JST["category_item"](templateData);
    }
    
    var messageHtml = createMessagesHtml(data.messages);
    
    $("#reflist_content").append(messageHtml);
    $("#reflist_content").append(categoryHtml);
    
    if (data.ancestors && data.ancestors.length > 0) {
        var ancestorHtml = createAncestorHtml(data, "referencelist.html", cruxOverrides);
        $("#left_nav_ancestors").replaceWith(ancestorHtml);
    } else {
        $("#left_nav_ancestors").remove();
    }

    var leftHtml = createLeftNavHtml(data, "referencelist.html", cruxOverrides);
    
    if (data.filtersDetail) {
        leftHtml += createFilterDetailHtml(data.filtersDetail);
    } else if (data.filters) {
        var filterNames = Object.getOwnPropertyNames(data.filters);
        for (var i = 0; i < filterNames.length; i++) {
            var filter = data.filters[filterNames[i]];
            leftHtml += createFilterHtml(filterNames[i], filter, {}, "referencelist.html");
        }
    }
    $("#left_nav_categories").html(leftHtml);

    var activeSibling = $("#nav_siblings li.active");
    if (activeSibling.length > 0) {
        var leftNavScrollTo = activeSibling.get(0).offsetTop;
        var categoriesTop = $("#nav_categories").get(0).offsetTop;
        $("#nav_categories").scrollTop(leftNavScrollTo-categoriesTop);
    }
    
    addSmoothScroll();
}

function getListItemElement(type) {
    var codeElementTypes = ['analysisopt', 'argument', 'class', 'com', 'constructor', 'cpp', 'function', 'matlabkeyword', 'method', 'mupadaxiom', 'mupadcat', 'mupadconstant', 'mupaddomain', 'mupadenvvar', 'mupadfunction', 'mupadgraphattr', 'mupadgraphprim', 'mupadkeyword', 'mupadlibrary', 'mupadoperator', 'net', 'object', 'package', 'property', 'reportgen', 'runtimecheck', 'simscapelanguage', 'sysobj'];
    return codeElementTypes.indexOf(type) >= 0 ? 'code' : 'span';
}
