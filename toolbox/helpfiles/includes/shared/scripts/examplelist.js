window.JST = window.JST || {};

function getExampleList(product) {
    var services = {
        "messagechannel":"examplelist",
        "requesthandler":"examplelist:getExamples",
        "webservice":getExampleListWebServiceUrl()
    }
    var params = parseParams();
    if (product) {
        params.product = product;
    }
    params.docroot = findDocRoot();
    if (! params.type) {
        params["type"] = "example";
    }
    
    requestHelpService(params, services, function(data) {
        addExampleItems(data);
    });
    
    if (product) {
        updateTopNavLinks({"category":params.category}, getCruxOverrides("gn", undefined, params.category));
    }

    var allParams = {};
    $.extend(allParams, params);
    var cruxParam = buildCruxParam(allParams, getCruxOverrides("lftnav", "example", undefined));
    var qs = $.param(cruxParam);
    var href = findDocRoot() + '/examples.html' + "?" + qs;
    var allExamplesLink = $('<a href="' + href + '"></a>');
    allExamplesLink.html(getLocalizedString("topnav_examples"));
    $("ul.nav_breadcrumb li a").after(allExamplesLink);
}

function getExampleListWebServiceUrl() {
    var lang = getPageLanguage() || "en";

    var release = getDocReleaseFromSearchBox();
    if (typeof getDocRelease === 'function') {
        release = getDocRelease();
    }


    return "/help/search/examplelist/doccenter/" + lang + "/" + release;
}

JST['category_item'] = _.template(
  '<div>' +
    '<% if (level == 1) { %>' +
      '<% populateH1Tag(category.name, typeLabel); %>' +
    '<% } else { %>' +
      '<h<%=level%> class="add_clear_both"><%= category.name %></h<%=level%>>' +
    '<% } %>' +
    '<% var leafItems = category["leaf-items"]; %>' +
    '<% if (leafItems && leafItems.length > 0) { %>' +
      '<section class="example_short_list">' +
      '<%_.each(category["leaf-items"], function(leafItem) { %>' +
        '<%= cardTemplate(leafItem) %>' +
      '<% }); %>' +
      '<% if (category.moreitems) { %>' +
        '<%= moreTemplate({"category":category}) %>' +
      '<% } %>' +
      '</section>' +
    '<% } %>'+
    '<% var childCats = category["child-categories"]; %>' +
    '<% if (childCats && childCats.length > 0) { %>' +
      '<% _.each(category["child-categories"], function(childCategory) { %>' +
        '<%= categoryTemplate({"category":childCategory,"level":level+1,"typeLabel":typeLabel,"categoryTemplate":categoryTemplate,"cardTemplate":cardTemplate,"moreTemplate":moreTemplate}) %>' +
      '<% }); %>' +
    '<% } else if (!leafItems || leafItems.length == 0) { %>' +
      'No Examples' +
    '<% } %>' +
  '</div>'
);

JST["more_examples_card"] = _.template(
    '<div class="card_container explorer_view add_long_title" data-ui-component="card">' +
        '<% var baseUrl = category.moreitemspath ? category.moreitemspath : "examples.html"; %>' +
        '<a href="<%= buildTopNavUrl(baseUrl, {"category":category.id}, {"area":"VM", "type":"example", "category":category.id}) %>">' +
          '<div class="card_body absolute_center">' +
            '<div class="panel-body"><h3><%=getViewMoreText(category.name)%></h3></div>' +
          '</div>' +
        '</a>' +
    '</div>'
);

function getViewMoreText(categoryName) {
    var viewMorePattern = getLocalizedString("topnav_view_more_examples");
    return viewMorePattern.replace(/\{0\}/, categoryName);
}

function addExampleItems(data) {
    var typeLabel = getLocalizedString("topnav_examples");
    var jsonData = {
        "category":data["category"],
        "level":1,
        "typeLabel":typeLabel,
        "categoryTemplate":JST["category_item"],
        "cardTemplate":JST["example_card"],
        "moreTemplate":JST["more_examples_card"]
    };

    var cruxOverrides = getCruxOverrides("lftnav", "example", data["category"]);
    var categoryHtml;
    if (data.category.leafcount === 0) {
        populateH1Tag(data.category.name, typeLabel);        
        categoryHtml = handleNoResults(data.category.name, data.relatedCategories);
    } else {
        var categoryHtml = JST['category_item'](jsonData);
    }
    $("#examplelist_content").html(categoryHtml);
    
    if (data.ancestors && data.ancestors.length > 0) {
        var ancestorHtml = createAncestorHtml(data, "examples.html", cruxOverrides);
        $("#left_nav_ancestors").replaceWith(ancestorHtml);
    } else {
        $("#left_nav_ancestors").remove();
    }

    var leftHtml = createLeftNavHtml(data, "examples.html", cruxOverrides);
    
    if (data.filtersDetail) {
        leftHtml += createFilterDetailHtml(data.filtersDetail);
    } else if (data.filters) {
        var filterNames = Object.getOwnPropertyNames(data.filters);
        for (var i = 0; i < filterNames.length; i++) {
            var filter = data.filters[filterNames[i]];
            leftHtml += createFilterHtml(filterNames[i], filter, {}, "examples.html");
        }
    }
    
    $("#left_nav_categories").html(leftHtml);
    
    var leftNavActive = $("ul.nav_toc li.active");
    var navCategories = $("#nav_categories");
    if (leftNavActive.length > 0 && navCategories.length > 0) {
        var leftNavScrollTo = leftNavActive.get(0).offsetTop;
        var categoriesTop = navCategories.get(0).offsetTop;
        navCategories.scrollTop(leftNavScrollTo-categoriesTop);
    }
    $(window).trigger('examples_cards_added');    
}
