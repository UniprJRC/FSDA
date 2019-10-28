$(function() {
    try {
        setUpTopNav();
    } catch (e) {
        // Swallow this error and move on.
    }
});

// Toggle a division as expanded or collapsed.
// Also toggle the arrow icon.
// Refer to the division and image by their IDs.
//
// "Collapsed" material is hidden using the
// display property in CSS.

// Used by adaptproduct function (see below)
// to support adaptive doc in the Windows
// version of the Help Browser.
var adaptiveIds = new Array();
function toggleexpander(blockid, arrowid) {
    arrow = document.getElementById(arrowid);
    block = document.getElementById(blockid);
    if (block.style.display == "none") {
        // Currently collapsed, so expand it.
        block.style.display = "block";
        arrow.src = arrow.src.replace("right", "down");
        arrow.title = getLocalizedString('click_to_collapse');
    }
    else {
        // Currently expanded, so collapse it.
        block.style.display = "none";
        arrow.src = arrow.src.replace("down", "right");
        arrow.title = getLocalizedString('click_to_expand');
    }
    return false; // Make browser ignore href.
}

// ===================================================
// Create and uniquely name two levels of upward navigation buttons
// for Functions -- By Category pages

var top_button_count = 0;
var current_section_id = 0;

function addTopOfPageButtons() {

    top_button_count = top_button_count + 1;

    var top_of_page_buttons =

        "<a class=\"pagenavimglink\" href=\"#top_of_page\" onMouseOver=\"document.images.uparrow" +
        top_button_count +
        ".src=\'doc_to_top_down.gif\'\;\" onMouseOut=\"document.images.uparrow" +
        top_button_count +
        ".src=\'doc_to_top_up.gif\'\;\"><img style=\"margin-top:0;margin-bottom:0px;padding-top:0;padding-bottom:0\" border=0 src=\"doc_to_top_up.gif\"  alt=\"" +
        getLocalizedString('back_to_top_of_page') +
        "\" title=\"" +
        getLocalizedString('back_to_top_of_page') +
        "\" name=\"uparrow" +
        top_button_count +
        "\">\&nbsp\;</a>";

    document.write(top_of_page_buttons);
}


function updateSectionId(id) {
    current_section_id = id;
}


function addTopOfSectionButtons() {

    top_button_count = top_button_count + 1;

    var top_of_page_buttons =

        "<a class=\"pagenavimglink\" href=" +
        "\"#" + current_section_id + "\"" +
        " onMouseOver=\"document.images.uparrow" +
        top_button_count +
        ".src=\'doc_to_section_down.gif\'\;\" onMouseOut=\"document.images.uparrow" +
        top_button_count +
        ".src=\'doc_to_section_up.gif\'\;\"><img style=\"margin-top:0;margin-bottom:0px;padding-top:0;padding-bottom:0\" border=0 src=\"doc_to_section_up.gif\"  alt=\"" +
        getLocalizedString('back_to_top_of_section') +
        "\" title=\"" +
        getLocalizedString('back_to_top_of_section') +
        "\" name=\"uparrow" +
        top_button_count +
        "\">\&nbsp\;</a>";

    document.write(top_of_page_buttons);
}

// ===================================================
// Create and write to the document stream HTML for
// the link to the Doc Feedback Survey site.
//
// Doing this through a JavaScript function is necessary
// to work around the an issue with pages that are found
// through the search facility of the help browser--
//
// When found as the result of a search,
// the document that is displayed in the Help browser
// is actually a temporary document with a trivial URL
// such as "text://5", not an actual page location.
//
// But the Help browser inserts a <BASE> element at the beginning
// of each such temporary page, and the <BASE> element stores the
// actual location.
//
// So this function tests the URL of the document for the expression "text://"
// and if that expression is found, attempts to use the URL stored in
// the <BASE> element.

function writeDocFeedbackSurveyLink() {
    var queryexpression = document.location.href;

    if (queryexpression.search(/text:\/\//) != -1) {
        var baseelement = document.getElementsByTagName("BASE")[0];
        queryexpression = baseelement.href;
    }
    survey_url_yes = "http://www.customersat3.com/TakeSurvey.asp?si=YU2FDmNEifg%3D&SF=" + queryexpression + "-YES";
    survey_url_no = "http://www.customersat3.com/TakeSurvey.asp?si=YU2FDmNEifg%3D&SF=" + queryexpression + "-NO";

    code = '<div style="padding-right:10px" class="feedbackblock">' + '<strong>' + getLocalizedString('was_this_topic_helpful') + '</strong> <input type="button" value="' + getLocalizedString('yes') + '" onClick="openWindow(\'' + survey_url_yes + '\',850,680, \'scrollbars=yes,resizable=yes\'); return false;"/>' + '&nbsp;&nbsp;' + '<input type="button" value="' + getLocalizedString('no') + '" onClick="openWindow(\'' + survey_url_no + '\',850,680, \'scrollbars=yes,resizable=yes\'); return false;"/>' + '</div>';
    document.write(code);
}


// Utility function replacing openWindow function used by the web-site survey link code.
// In the help browser, the original code would create a blank window before loading the URL into the system browser.
function openWindow(url, width, height, options, name) {
    // ignore the arguments, except url
    document.location = url;
} // end function openWindow


// Utility function for linking to feedback survey, as of R2012b.
function openFeedbackWindow(url) {
    window.open(url, "_blank");
} // end function openFeedbackWindow


// ===================================================
// Workaround for G801125.
// This global object check tests for IE8 or lower.
if (document.all && !document.getElementsByClassName) {
    document.createElement("section");
}


// ===================================================
// Function reference pages

$(window).load(function () {
    // Perform breadcrumb check in window load, since all the images in the breadcrumb
    // need to be loaded for correct width calculations.
    //getBreadcrumb().mwBreadcrumb();
    //delay the expanding to enable the page to load completely.
    setTimeout(function () {
        expandCollapsedContent();
    }, 0);
    $('body').scrollspy({ target: '.offcanvas_nav', offset: getScrollTopAdjustment() + 1});
});

$(document).ready(function () {
    // this function is derived from an earlier version of jquery
    // and we are using it here for backwards compatability.
    $.browserInfo = function () {
        var ua = navigator.userAgent.toLowerCase();
        var match = /(chrome)[ \/]([\w.]+)/.exec(ua) ||
            /(webkit)[ \/]([\w.]+)/.exec(ua) ||
            /(opera)(?:.*version|)[ \/]([\w.]+)/.exec(ua) ||
            /(msie) ([\w.]+)/.exec(ua) ||
            ua.indexOf("compatible") < 0 && /(mozilla)(?:.*? rv:([\w.]+)|)/.exec(ua) ||
            [];

        var browserName = match ? match [1] : navigator.appName;
        var browserVersion = match ? match[2] : navigator.appVersion;
        return {
            name: browserName,
            version: browserVersion
        }
    };
	
	$('select').each(function () {
		var select = $(this);
		var selectedValue = select.find('option[selected]').val();

		if (selectedValue) {
		  select.val(selectedValue);
		} else {
		  select.prop('selectedIndex', 0);
		}
	});


    $(".nav_scrollspy").addClass("nav");
    $(".nav_scrollspy ul").addClass("nav");

    if ($.getParameterByName("browser") !== "F1help") {
        if ($.fn.setupToc) {
            $('div.toc_container_wrapper').setupToc();
        }
    }

    //Perform JS code which has any user visible impact first.
    //Check image sizes. Do not scale animated images or any images with hotspots.

    var searchCrumbContainer = $('#search_crumb_container');
    //old template
    if (searchCrumbContainer.length > 0) {
        $('#doc_center_content img:not(".animated-image, [usemap]"), #content_container2 img:not(".animated-image, [usemap]")').scaleImage();
    }

    $('#doc_center_content img.animated-image, #content_container2 img.animated-image').animateImage();


    $('.collapse').on('hidden.bs.collapse shown.bs.collapse', function () {
        var id = $(this).attr("id");
        var expandParent = $("[data-target='#" + id +"']");
        checkExpandAllLinkState(expandParent);
        triggerContentResize();
    });

    addSmoothScroll();

    $('#content_container .expandAllLink').click(function (e) {
        e.stopPropagation();
        var link = $(this);
        var allExpanders = link.closest('.expandableContent').find("[data-toggle='collapse']");
        if (link.data('allexpanded')) {
            doCollapse(allExpanders);
            setExpandAllLinkState(link, "collapsed");
        } else {
            doExpand(allExpanders);
            setExpandAllLinkState(link, "expanded");
        }
    });

    $('#expandAllPage').click(function () {
        var allExpanders = $("#content_container").find("[data-toggle='collapse']");
        if ($(this).data('allexpanded')) {
            doCollapseAllInPage(allExpanders);
            $(this).data('allexpanded', false);
            $(this).html(getLocalizedString('expand_all_in_page'));
        } else {
            doExpandAllInPage(allExpanders);
            $(this).data('allexpanded', true);
            $(this).html(getLocalizedString('collapse_all_in_page'));
        }
    });

    applySearchHighlight();

    $(window).bind('toc_resize', function () {
        $('#search_crumb_container').width($('#content_container').width());
        getBreadcrumb().trigger('window_resize');
    });

    $(window).bind('resize', function (e) {
        $('#search_crumb_container').width($('#content_container').width());
        getBreadcrumb().trigger('window_resize');
        $('.toc_container_wrapper').trigger('window_resize');
        refreshScrollSpy();
    });

    $(window).bind('content_resize', function (e) {
        $('.toc_container_wrapper').trigger('window_resize');
        refreshScrollSpy();
    });

    $(window).bind('intrnllnk_clicked', function (e) {
        expandCollapsedContent();
    });

    $('#search_crumb_container').width($('#content_container').width());
    getBreadcrumb().trigger('window_resize');

    var $highlighted = $('[class^=highlight_]');
    $highlighted.each(function(){
        $(this).data("class-backup", $(this).attr("class"));
    });
    
    var oldOnpopstate = window.onpopstate;
    window.onpopstate = function(event) {
        if(event.state && event.state.highlight == "off")
            $highlighted.each(function(){
                $(this).attr("class","");
            });
        else
            $highlighted.each(function(){
                $(this).attr("class", $(this).data("class-backup"));
            });
        
        if(typeof oldOnpopstate === "function") oldOnpopstate(event);
    };

});

function setExpandAllLinkState(link, state) {
    if (state === 'expanded') {
        link.data('allexpanded', true);
        link.html(getLocalizedString('collapse_all'));
    } else if (state === 'collapsed') {
        link.data('allexpanded', false);
        link.html(getLocalizedString('expand_all'));
    }
}

var UrlParser = function (url){
    this.url = url;
    this.removeVar = function(name){
        name = name.replace(/[\-\[\]\/\{\}\(\)\*\+\?\.\\\^\$\|]/g, "\\$&");
        var expr = '([&\?])' + name + '=?[^&]*&?';
        var newurl = this.url.replace(new RegExp(expr),'$1').replace(/\?$/,'');
        return new UrlParser(newurl);
    }
}

function applySearchHighlight() {
    if ($.fn.highlight) {
        var highlighterCSSClass = ['highlight_01', 'highlight_02' , 'highlight_03', 'highlight_04', 'highlight_05'];
        var searchHighlightTerm = $.getParameterByName('searchHighlight');
        var additionalHighlight = '';
        if (localStorage) {
            additionalHighlight = localStorage[searchHighlightTerm];
        }
        var highlightExpand = {};
        if (additionalHighlight) {
            highlightExpand = JSON.parse(additionalHighlight);
        }
        
        if (searchHighlightTerm.length > 0) {
            var searchHighlightArray = searchHighlightTerm.match(/"[^"]+"|\S+/g);
            $.each(searchHighlightArray, function (index, value) {
                var searchTerm = value.replace(/^"|"$/g, '');
                var cssClass = highlighterCSSClass[index % highlighterCSSClass.length];
                $("#doc_center_content").highlight(searchTerm,
                    {className: cssClass, wordsOnly: true});

                if (!$.isEmptyObject(highlightExpand) && highlightExpand[searchTerm.toLowerCase()] !== undefined) {
                    $.each(highlightExpand[searchTerm.toLowerCase()], function (highlightIdx, highlightVal) {
                        var highlightTerm = highlightVal;
                        $("#doc_center_content").highlight(highlightTerm, {className: cssClass, wordsOnly: true});
                    });
                }

                var elements = $("#doc_center_content").find("." + cssClass);

                setTimeout(function () {
                    expandHighlightedElements(elements);
                }, 50);
            });

            $(document).keyup(function (e) {
                if (e.which === 27 && // escape key maps to keycode `27`
                    /([&\?])searchHighlight=?[^&]*&?/.test(document.location.search)) {
                    var classArray = $.map(highlighterCSSClass, function (value) {
                        return "." + value;
                    });
                    var highlightedEl = $(classArray.join(","));
                    $.each(highlightedEl, function () {
                        $(this).removeClass();
                    });
                    
                    var url = (new UrlParser(document.location.href))
                        .removeVar('searchHighlight').removeVar('s_tid').url;
            
                    // push this new url to the browser history
                    window.history.pushState({highlight: "off"}, "", url);
                }
            });
        }
    }
}

function expandHighlightedElements(elements) {
    $.each(elements, function () {
        if (!$(this).is(":visible")) {
            var collapsedParent = $(this).closest('.collapse');
            var id = collapsedParent.attr("id");
            var expandParent = $("[data-target='#" + id +"']");
            if (expandParent.length > 0) {
                prepareEltForExpansion(expandParent, true);
                doExpand(expandParent);
            }
        }
    });
}

//helper method to fetch the breadcrumb.
function getBreadcrumb() {
    var breadcrumb;
    if ($("#breadcrumbs").length != 0) {
        breadcrumb = $("#breadcrumbs");
    } else {
        breadcrumb = $(".breadcrumbs:first");
    }
    return breadcrumb;
}

function expandCollapsedContent() {
    if (location.hash.length > 0) {

        var target = getInternalLinkTarget(location.hash);
        if (target.length > 0) {
            var expandParent = getExpandParentForAnchorTarget(target);
            prepareEltForExpansion(expandParent);

            //scroll to the target first.
            var scrollParameter = getScrollParameter();
            var scrollTop = target.offset().top - getScrollTopAdjustment();
            $(scrollParameter).scrollTop(scrollTop);
            $('.anchor_hinting').removeClass('anchor_hinting');
            target.addClass('anchor_hinting');
            setTimeout(function () {
                target.removeClass('anchor_hinting');
            }, 5000);
            doExpand(expandParent);
        }
    }
}

function prepareEltForExpansion(elt, noAnimation) {
    doExpandNestedParent(elt, noAnimation);
}

function getExpandParentForAnchorTarget(target) {
    var collapseToggle = target.parents("[data-toggle='collapse']");
    if (collapseToggle.length > 0) {
        return collapseToggle.first();
    }
    return target;
}

function addSmoothScroll() {
    $(".intrnllnk").each(function () {
        var hash = this.hash;
        var target = getInternalLinkTarget(hash);
        if (target.length > 0) {
            $(this).click(function (evt) {
                evt.preventDefault();
                var expandParent = getExpandParentForAnchorTarget(target);
                prepareEltForExpansion(expandParent);

                var scrollParameter = getScrollParameter();
                var scrollTop = target.offset().top - getScrollTopAdjustment();
                $('.anchor_hinting').removeClass('anchor_hinting');
                $(scrollParameter).animate({scrollTop: scrollTop}, 700, function () {
                    target.addClass('anchor_hinting');
                    setTimeout(function () {
                        target.removeClass('anchor_hinting');
                    }, 5000);
                    doExpand(expandParent);
                });
                if(history.pushState) {
                    history.pushState(null, null, hash);
                }
                else {
                    location.hash = hash;
                }
                refreshScrollSpy();
            })
        }
    });
}

function getInternalLinkTarget(hash) {
    //search for anchor with given hash as "name" atrribute value;
    var target = [];

    // sanitizing for the url
    hash = hash.split('?')[0];

    //Remove the first '#' character from the name attribute. Escape any special character from the name/id.
    var escapedHash = hash.substring(1).replace(/([;&,.+*~':"!^#$%@\[\]\(\)=>\|])/g, '\\$1');

    target = $("#" + escapedHash);
    return target;
}

// Assumes that the elt passed in has data-toggle=collapse, and
// data-target specifying the element which needs to be expanded/collapsed
function findExpandableContent(elt) {
    var target = elt.data('target');
    return $(target);
}

function triggerContentResize() {
    $(window).trigger('content_resize');
}

function doExpand(elt_array, noAnimation) {
    $.each(elt_array, function(i, elt) {
        var expandable = findExpandableContent($(elt));
        expandable.collapse('show');
    });
}

function doCollapse(elt_array, noAnimation) {
    $.each(elt_array, function(i, elt) {
        var expandable = findExpandableContent($(elt));
        expandable.collapse('hide');
    });
}

function checkExpandAllLinkState(elt) {
    //Check if the expandable elt is nested within another expandable elt.
    var expandableParent = elt.closest('.expandableContent');

    // If element is not nested, or there is not expand all link, return
    if (expandableParent.length === 0 || expandableParent.find('.expandAllLink').length === 0) {
        return;
    }
    var expandAllLink = expandableParent.find('.expandAllLink:first');
    var expandableChildren = expandableParent.find("[data-toggle='collapse']");

    if (elt.hasClass('collapsed')) {
        var allChildrenCollapsed = true;
        expandableChildren.each(function () {
            if (!$(this).hasClass("collapsed")) {
                allChildrenCollapsed = false;
            }
        });
        if (allChildrenCollapsed) {
            setExpandAllLinkState(expandAllLink, "collapsed");
        }
    } else {
        var allChildrenExpanded = true;
        expandableChildren.each(function () {
            if ($(this).hasClass("collapsed")) {
                allChildrenExpanded = false;
            }
        });

        if (allChildrenExpanded) {
            setExpandAllLinkState(expandAllLink, "expanded");
        }
    }
}

function doExpandNestedParent(elt, noAnimation) {
    var expandableParent = elt.parents("[data-toggle='collapse']").first();
    if (expandableParent.length > 0) {
        var expandable = findExpandableContent(expandableParent);
        expandable.collapse('show');
    }
}

function doExpandAllInPage(elt_array) {
    $.each(elt_array, function(i, elt) {
        var expandable = findExpandableContent($(elt));
        expandable.collapse('show');
    });
    $.each($('#content_container .expandAllLink'), function (i, link) {
        setExpandAllLinkState($(link), "expanded");
    });
}

function doCollapseAllInPage(elt_array) {
    $.each(elt_array, function(i, elt) {
        var expandable = findExpandableContent($(elt));
        expandable.collapse('hide');
    });
    $.each($('#content_container .expandAllLink'), function (i, link) {
        setExpandAllLinkState($(link), "collapsed");
    });
}

function getScrollParameter() {
    return "html, body";
}

function getScrollTopAdjustment() {
    var scrollTop = 0;
    var searchCrumbContainer = $('#search_crumb_container');
    //old template
    if (searchCrumbContainer.length > 0) {
        if (searchCrumbContainer.css('position') === 'fixed') {
            scrollTop = 103;
        } else {
            scrollTop = 10;
        }
    } else {
        if ($('.sticky_header_container').is(":visible")) {
            scrollTop = $('.sticky_header_container').height() + 10;
        } else {
            scrollTop = 10;
        }

        
    }
    return scrollTop;
}

function refreshScrollSpy() {
    try {
        $('body#responsive_offcanvas').each(function () {
            $(this).data('bs.scrollspy').options.offset = getScrollTopAdjustment() + 1;
            var $spy = $(this).scrollspy('refresh');
        });
        $(window).scroll();
    } catch (e) {}
}

function setUpTopNav() {
    var topNavContainer = findTopNavContainer();
    doTopNavLayout(topNavContainer);
    $(window).resize(function() {
        adjustTopNavWidth();
    });
}

function adjustTopNavWidth() {
    var topNavContainer = findTopNavContainer();
    $("#topnav_more").hide();
    $("#topnav_more").removeClass("crux_nav_active");
    $("#topnav_more > ul > li").detach().appendTo(topNavContainer);
    var moreItem = $("#topnav_more").detach();
    topNavContainer.append(moreItem);

    doTopNavLayout(topNavContainer);
}

function doTopNavLayout(topNavContainer) {
    if (!window.topnav_link_positions || window.topnav_link_positions.length === 0) {
        if (!inspectTopNavLinks(topNavContainer)) {
            return;
        }
    }
    // The selected item should be selected if it gets added to the "More" menu.
    topNavContainer.children("li.crux_nav_active").addClass("active");

    var titleElt = findTitleElement();
    var headerArea = $("#section_header_title");
    var targetWidth = headerArea.width() - titleElt.width() - titleElt.offset().left - $("#topnav_more").width() - 5;
    
    var menuStart = -1;
    for (var i = 0; i < window.topnav_link_positions.length; i++) {
        if (window.topnav_link_positions[i] > targetWidth) {
            menuStart = i;
            break;
        }
    }
    
    if (menuStart >= 0) {
        var navItems = topNavContainer.children("li");
        navItems.slice(menuStart, navItems.length-1).detach().appendTo("#topnav_more > ul");
        if ($("#topnav_more > ul").find(".crux_nav_active").length > 0) {
            $("#topnav_more").addClass("crux_nav_active");
        }
        $("#topnav_more").show();
    }
}

function findTopNavContainer() {
    var topNav = $("#topnav_main");
    if (topNav.length === 0) {
        topNav = $("#crux_nav_desktop");
    }
    return topNav;
}

function findTitleElement() {
    var titleElement = $("div.section_header_title > h1");
    if (titleElement.length === 0) {
        titleElement = $("div.section_header_title > p.h1");
    }
    return titleElement;
}

function inspectTopNavLinks(topNavContainer) {
    var navItems = topNavContainer.children("li");
    if (!navItems || navItems.length === 0 || !navItems.is(":visible")) {
        // There are no nav items showing.
        window.topnav_link_positions = [];
        return false;
    }

    updateTopNavItemList(navItems);

    var itemsLeft = $(navItems[0]).offset().left;
    var linkPositions = [];
    
    navItems.each(function(idx,item) {
        item = $(item);
        if (item.attr("id") === "topnav_more") {
            return false;
        }
        
        var widthToItem = item.offset().left + item.width() - itemsLeft;
        linkPositions.push(widthToItem);
    });
    
    window.topnav_link_positions = linkPositions;
    return true;
}

function updateTopNavItemList(navItems) {
     navItems.each(function(idx,item) {
        
        item = $(item);

        var type = getTypeForId(item.attr("id"));
        var cruxParam = buildCruxParam(getCruxOverrides("gn", type));
        var qs = $.param(cruxParam);

        var aTags = item.children("a");
        var aTag = aTags[0];
        var href = aTag.getAttribute("href");

        if (href && qs) {
            // Remove the CRUX param if it already exists. A new one will be added.
            href = href.replace(/(&?s_tid=)(.[^&]*)/,"");

            href += href.indexOf('?') > 0 ? '&' : '?';
            href += qs;

            aTag.setAttribute("href", href);
        }
    });
}

function getTypeForId(id) {
    var type = id.replace("crux_nav_", "");
    return type;
}

function requestLeftNavItems(shortName, pathToDocRoot) {
    pathToDocRoot = pathToDocRoot.replace(/\/+$/,'');
    
    var services = {
      "messagechannel":"prodfilter",
      "requesthandler":"productfilter:handleSelectedProducts",
      "webservice": getProdFilterWebServiceUrl()
    }
    requestHelpService({}, services, function(data) {
        var prodlist = data.prodnavlist;
        if (typeof prodlist === "string") {
            prodlist = $.parseJSON(prodlist);
        }
        buildLeftNav(shortName, pathToDocRoot, prodlist);
    });
}

function getProdFilterWebServiceUrl() {
    var release = getDocReleaseFromSearchBox();
    if (typeof getDocRelease === 'function') {
        release = getDocRelease();
    }
    
    return "/help/search/prodfilter/doccenter/en/" + release;
}

function getDocReleaseFromSearchBox() {

    var localeEl = $("#docsearch_form");
    return localeEl.attr('data-release');
}

function buildLeftNav(shortName, pathToDocRoot, data) {
  var navElt = $("nav > ul.nav_toc");
  var catItems = navElt.find("li").detach();
  var selectedProd;
  
  
  $(data).each(function(idx,product) {
    var prodLineItem = $("<li></li>");
    var prodLink = $("<a></a>");
    prodLink.html(product.displayname);
    prodLink.attr("href", pathToDocRoot + "/" + product.helplocation);
    prodLineItem.append(prodLink);
    if (shortName === product.shortname) {
      selectedProd = prodLineItem;
      var catsContainer = $("<ul></ul>");
      catsContainer.append(catItems);
      prodLineItem.append(catsContainer);
      prodLineItem.addClass("active");
    }
    navElt.append(prodLineItem);
    
  });

  if (selectedProd.length > 0) {
      var prodListTop = navElt.get(0).offsetTop;
      var scrollToTop = selectedProd.get(0).offsetTop;
      
      // Use a window onload event here to make sure that the
      // CSS is loaded before attempting to scroll.
      $(window).on("load", function() {
        $("#sidebar").scrollTop(scrollToTop);
      });
  }
}
function buildCruxParam(cruxOverrides) {
    // Use allParams here so I don't modify cruxOverrides. 
    var allParams = {};
    $.extend(allParams, cruxOverrides);

    var cruxParam = {};
    var name = "s_tid"
    var value = "CRUX";
    if (allParams.area) {
        value = value + "_" + allParams.area;
    }
    if (allParams.type) {
        value = value + "_" + allParams.type;
    }
    cruxParam[name] = value;
    return cruxParam;
}

function getCruxOverrides(area, type) {
    var cruxOverrides = {};
    if (area) {
        cruxOverrides["area"] = area;
    }
    if (type) {
        cruxOverrides["type"] = type;
    }
    return cruxOverrides;
}

// Copyright 2002-2019 The MathWorks, Inc.
