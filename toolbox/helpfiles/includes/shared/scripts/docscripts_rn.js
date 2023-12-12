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
        expandCollapsedContent(true);
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

    addSmoothScroll();

    $("#content_container").delegate(".expandAllLink", "click", function (e) {
        e.stopPropagation();
        var link = $(this);
        if (link.data('allexpanded')) {
            doCollapse(link);
            setExpandAllLinkState(link, "collapsed");
        } else {
            doExpand(link);
            setExpandAllLinkState(link, "expanded");
        }
    });

    $("#content_container").delegate(".expand", "click", function (e) {
        e.stopPropagation();
        doToggle($(this));
        return false;
    });

    $('#expandAllPage').click(function () {
        if ($(this).data('allexpanded')) {
            doCollapseAllInPage($('#content_container .expand'));
            $(this).data('allexpanded', false);
            $(this).html(getLocalizedString('expand_all_in_page'));
        } else {
            doExpandAllInPage($('#content_container .expand'));
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
    });

    $(window).bind('content_resize', function (e) {
        $('.toc_container_wrapper').trigger('window_resize');
    });

    $(window).bind('intrnllnk_clicked', function (e) {
        expandCollapsedContent();
    });

    $('#search_crumb_container').width($('#content_container').width());
    getBreadcrumb().trigger('window_resize');

});

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
  var catItems = $("#rn-form-categories").detach();
  // Get rid of the header completely.
  $("#rn-category-group").remove();

  var navElt = $("#rn-nav-categories > ul.nav_toc");
  $(data).each(function(idx,product) {
    var prodLineItem = $("<li></li>");
    var prodLink = $("<a></a>");
    prodLink.html(product.displayname);
    prodLink.attr("href", pathToDocRoot + "/" + getReleaseNotesPage(product.helplocation));
    prodLineItem.append(prodLink);
    if (shortName === product.shortname) {
      var catsContainer = $("<div></div>");
      catsContainer.append(catItems);
      prodLineItem.append(catsContainer);
      prodLineItem.addClass("active");
    }
    navElt.append(prodLineItem);    
  });

  var activeSibling = $("#rn-nav-categories > ul.nav_toc > li.active");
  if (activeSibling.length > 0) {
    var leftNavScrollTo = activeSibling.get(0).offsetTop;
    var categoriesTop = $("#rn-nav-categories").get(0).offsetTop;
    $("#rn-nav-categories").scrollTop(leftNavScrollTo-categoriesTop);
  }
}

function getReleaseNotesPage(page) {
    return page.replace(/index\.html$/, 'release-notes.html');
}

function setExpandAllLinkState(link, state) {
    if (state === 'expanded') {
        link.data('allexpanded', true);
        link.html(getLocalizedString('collapse_all'));
    } else if (state === 'collapsed') {
        link.data('allexpanded', false);
        link.html(getLocalizedString('expand_all'));
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
                // apply additional highlight expansion
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
                if (e.which === 27) {
                    var classArray = $.map(highlighterCSSClass, function (value) {
                        return "." + value;
                    });
                    var highlightedEl = $(classArray.join(","));
                    $.each(highlightedEl, function () {
                        $(this).removeClass();
                    });
                }
            });
        }
    }
}

function expandHighlightedElements(elements) {
    $.each(elements, function () {
        if (!$(this).is(":visible")) {
            var collapsedParent = $(this).closest('.collapse').prev();
            if (collapsedParent.length > 0) {
                prepareEltForExpansion(collapsedParent, true);
                if (collapsedParent.hasClass('expand') && !collapsedParent.hasClass('expanded')) {
                    doToggle(collapsedParent, true);
                }
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

function expandCollapsedContent(noAnimation) {
    if (location.hash.length > 0) {

        var target = getInternalLinkTarget(location.hash);
        if (target.length > 0) {
            var nextSibling = getNextSiblingForAnchorTarget(target);
            prepareEltForExpansion(nextSibling, noAnimation);

            //scroll to the target first.
            var scrollParameter = getScrollParameter();
            var scrollTop = target.offset().top - getScrollTopAdjustment();
            $(scrollParameter).scrollTop(scrollTop);
            $('.anchor_hinting').removeClass('anchor_hinting');
            nextSibling.addClass('anchor_hinting');
            setTimeout(function () {
                nextSibling.removeClass('anchor_hinting');
            }, 5000);


            if (nextSibling.hasClass('expand') && !nextSibling.hasClass('expanded')) {
                doToggle(nextSibling);
            }
        }
    }
}

function prepareEltForExpansion(elt, noAnimation) {
    if (elt.hasClass('expand')) {
        doExpandNestedParent(elt, noAnimation);
    } else {
    	doExpandParent(elt, noAnimation);
    }    
}

function getDocReleaseFromSearchBox() {
    var localeEl = $("#docsearch_form");
    return localeEl.attr('data-release');
}

function getNextSiblingForAnchorTarget(target) {
    if (target.is('div') || target.is('span') || target.is('a')) {
        return target.next();
    } else if (target.parent().hasClass('expand')) {
        return target.parent();
    } else {
        return target;
    }
}

function addSmoothScroll() {
    $(".intrnllnk").each(function () {
        var hash = this.hash;
        var target = getInternalLinkTarget(hash);
        if (target.length > 0) {
            $(this).click(function (evt) {
                evt.preventDefault();
                var nextSibling = getNextSiblingForAnchorTarget(target);
                prepareEltForExpansion(nextSibling);

                var scrollParameter = getScrollParameter();
                var scrollTop = target.offset().top - getScrollTopAdjustment();
                $('.anchor_hinting').removeClass('anchor_hinting');
                $(scrollParameter).animate({scrollTop: scrollTop}, 700, function () {
                    nextSibling.addClass('anchor_hinting');
                    setTimeout(function () {
                        nextSibling.removeClass('anchor_hinting');
                    }, 5000);
                    if (nextSibling.hasClass('expand') && !nextSibling.hasClass('expanded')) {
                        doToggle(nextSibling);
                    }
                });
                if(history.pushState) {
                    history.pushState(null, null, hash);
                }
                else {
                    location.hash = hash;
                }
            })
        }
    });
}

function getInternalLinkTarget(hash) {
    //search for anchor with given hash as "name" atrribute value;
    var target = [];

    //Remove the first '#' character from the name attribute. Escape any special character from the name/id.
    var escapedHash = hash.substring(1).replace(/([;&,.+*~':"!^#$%@\[\]\(\)=>\|])/g, '\\$1');

    target = $("#" + escapedHash);
    return target;
}

function findExpandableContent(elt) {
    if (!elt.hasClass('expandableContent')) {
        elt = elt.closest('.expandableContent');
    }
    return elt;
}

function triggerContentResize() {
    $(window).trigger('content_resize');
}

function doExpand(elt, noAnimation) {
    var expandable = findExpandableContent(elt);
    var browserInfo = $.browserInfo();
    if ((browserInfo.name === 'msie' && browserInfo.version <= 8) || noAnimation) {
        expandable.find('.collapse').show();
        expandable.find('.expand').addClass('expanded');
        triggerContentResize();
        checkExpandAllLinkState(elt);
        return;
    }

    expandable.find('.collapse').slideDown(function () {
        expandable.find('.expand').addClass('expanded');
        triggerContentResize();
        checkExpandAllLinkState(elt);
    });
}

function doCollapse(elt, noAnimation) {
    var expandable = findExpandableContent(elt);
    //Before collapsing, check if the collapse child has any expandableContent children.
    //If it does, those divs need to be collapsed and not the parent.
    var collapsedChild = expandable.children('.collapse');
    if (collapsedChild.children('.expandableContent').length > 0) {
        expandable = expandable.children('.collapse').children('.expandableContent');
    }

    var browserInfo = $.browserInfo();
    if ((browserInfo.name === 'msie' && browserInfo.version <= 8) || noAnimation) {
        expandable.find('.collapse').hide();
        expandable.find('.expand').removeClass('expanded');
        triggerContentResize();
        checkExpandAllLinkState(elt);
        return;
    }

    expandable.find('.collapse').slideUp(function () {
        expandable.find('.expand').removeClass('expanded');
        triggerContentResize();
        checkExpandAllLinkState(elt);
    });
}

function doToggle(elt, noAnimation) {
    var expandable = findExpandableContent(elt);
    var browserInfo = $.browserInfo();
    if ((browserInfo.name === 'msie' && browserInfo.version <= 8) || noAnimation) {
        expandable.children('.collapse').toggle();
        expandable.children('.expand').toggleClass('expanded');
        triggerContentResize();
        checkExpandAllLinkState(elt);
        return;
    }

    expandable.children('.expand').toggleClass('expanded');
    expandable.children('.collapse').slideToggle(function () {
        triggerContentResize();
        checkExpandAllLinkState(elt);
    });
}

function checkExpandAllLinkState(elt) {
    //Check if the expandable elt is nested within another expandable elt.
    var expandableParent = elt.parents('.expandableContent:eq(1)');

    // If element is not nested, or there is not expand all link, return
    if (expandableParent.length === 0 || expandableParent.find('.expandAllLink').length === 0) {
        return;
    }
    var expandAllLink = expandableParent.find('.expandAllLink:first');
    var expandableChildren = expandableParent.find('.expandableContent');


    if (elt.hasClass('expanded')) {
        var allChildrenExpanded = true;
        expandableChildren.each(function () {
            var expand = $(this).find('.expand');
            if (!expand.hasClass("expanded")) {
                allChildrenExpanded = false;
            }
        });
        if (allChildrenExpanded) {
            setExpandAllLinkState(expandAllLink, "expanded");
        }
    } else {
        var allChildrenCollapsed = true;
        expandableChildren.each(function () {
            var expand = $(this).find('.expand');
            if (expand.hasClass("expanded")) {
                allChildrenCollapsed = false;
            }
        });

        if (allChildrenCollapsed) {
            setExpandAllLinkState(expandAllLink, "collapsed");
        }
    }
}

function doExpandNestedParent(elt, noAnimation) {
    var expandable = findExpandableContent(elt);
    var expandableParent = expandable.closest('.collapse').siblings('.expand');
    if (expandableParent.length > 0) {
        if (!expandableParent.hasClass('expanded')) {
            doToggle(expandableParent, noAnimation);
        }
    }
}

// This method is used to support the old HTML template for the ref pages. When all the ref pages have been
// converted to the new template, consolidate this method with the doExpandNestedParent method above.
function doExpandParent(elt, noAnimation) {
	var expandable = elt.closest('.collapse');
	while (!elt.is(":visible")) {
	    var expandableParent = expandable.siblings('.expand');
	    if (expandableParent.length > 0) {
	        if (!expandableParent.hasClass('expanded')) {
	            doToggle(expandableParent, noAnimation);
	        }
	    }
	    expandable = expandableParent.closest('.collapse');
	}
}

function doExpandAllInPage(elt) {
    var expandable = findExpandableContent(elt);
    expandable.find('.collapse').show();
    expandable.find('.expand').addClass('expanded');
    triggerContentResize();
    $.each($('#content_container .expandAllLink'), function (i, link) {
        setExpandAllLinkState($(link), "expanded");
    });
}

function doCollapseAllInPage(elt) {
    var expandable = findExpandableContent(elt);
    expandable.find('.collapse').hide();
    expandable.find('.expand').removeClass('expanded');
    triggerContentResize();
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
        scrollTop = $('.sticky_header_container').height() + 10;
    }
    return scrollTop;
}

// Copyright 2016-2019 The MathWorks, Inc.
