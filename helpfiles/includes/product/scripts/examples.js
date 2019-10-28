/*jslint nomen: true*/
/*global window, JST, _, $, document*/

window.JST = window.JST || {};

JST.example_card = _.template(
    '<div class="card_container explorer_view add_long_title" data-ui-component="card">' +
        '<a href="<%= url %>">' +
        '<div class="card_media" style="background-image:url(<%= thumbnail_url %>);">' +
        '<img alt="<%= title %>" src="<%= thumbnail_url %>" />' +
        '</div>' +
        '<div class="card_body">' +
        '<div class="panel panel_color_transparent panel_color_fill">' +
        '<div class="panel-heading">' +
        '<h3><%= title %></h3>' +
        '</div>' +
        '<div class="panel-body" style="overflow:hidden;">' +
        '<p class="card_description"><%= getDescription(description) %></p>' +
        '</div></div></div></a>' +
        '<div class="card_footer"><div class="row">' +
        '  <div class="col-xs-12"><div class="card_action_primary">' +
        '    <% if (open_command) { %><a href="<%= open_command %>"><%= getTranslateOpenExampleString(type) %></a><% } %>' +
        '  </div></div>' +
        '</div>' +
        '<div class="card_source">' +
        '<div class="card_type" data-toggle="tooltip" data-placement="top" title="MathWorks"> <span class="icon-membrane" aria-hidden="true"></span> </div>' +
        '</div>' +
        '</div>' +
        '</div>'
);

JST.example_cards  = _.template(
    '<div class="col-xs-12">' +
        '<section class="example_short_list">' +
        '<h2 class="add_clear_both"><a href="<%= product_example_location%>"><%= label %></a></h2>' +
            '<% _.each(examples, function(example) { %>' +
                '<%=    JST.example_card(example) %>' +
            '<% })%>' +
            '<a href="<%= product_example_location%>">' +
        '<div class="card_container explorer_view" data-ui-component="ui-example-card-view">' +
        '<div class="card_body absolute_center">' +
        '<div class="panel panel_color_transparent panel_color_fill">' +
        '<div class="panel-body" style="overflow:hidden;">' +
        '<h3><%= getTranslateViewMoreExampleString(label) %></h3>' +
        '</div></div></div>' +
        '</div>' +
        '</a>' +
        '</section>' +
        '</div>'
);


function convertExampleLinksToCards() {
    'use strict';
    var jsonSuffix = getJsonSuffix();
    var jsonUrl = "examples" + jsonSuffix + ".json";
    $.getJSON(jsonUrl)
        .done(function (data) {
            $.each(data, function () {
                var el = $('li a[href="' + this.url + '"]').parent();
                if (el.length === 0) {
                    el = $('li span[class="missinglink"]').parent();
                    el.replaceWith(JST.example_card(this));
                } else {
                    el.replaceWith(JST.example_card(this));
                }
            });
        }).fail(function () {
            $('.an-example').show();
        }).always(function () {
            $('h2').show();
            $('h3').show();
            $('h4').show();
            triggerExamplesCardsAdded();
        });
}

function loadExamplesLandingPageFromJSON(selectedProductList) {
    'use strict';
    $.getJSON('all_product_examples.json', function (allProductJson) {
        if (selectedProductList === null || selectedProductList.length === 0) {
            handleFeaturedExamples(allProductJson);
        } else {
            var filteredJson = allProductJson.filter(function (obj) {
                if ('shortname' in obj) {
                    for (var i = 0; i < selectedProductList.length; i++) {
                        if (selectedProductList[i] === obj.shortname) {
                            return true;
                        }
                    }
                    return false;
                } else {
                    return false;
                }
            });
            handleFeaturedExamples(filteredJson);
        }
    });
}

function triggerExamplesCardsAdded() {
    $(window).trigger('examples_cards_added');
}

function setExamplesVisibility() {
    'use strict';
    var prodFilterData = {};
    var services = {
        "messagechannel":"featuredexamples",
        "requesthandler":"featuredexamples:handleFeaturedExamples",
        "webservice":null
    };
    requestHelpService(prodFilterData, services, function(data) {
        var prodList = $.parseJSON(data.prodlist);
        handleFeaturedExamples(prodList);
    });
}

function handleFeaturedExamples(productList) {
    'use strict';
    loadExamples('matlab');
    if (productList.length > 0) {
        handleComingFromProductList(productList, 'all_product_list', 'all_products');
        checkAndLoadExamples(productList);
    }
}

function checkAndLoadExamples(productList) {
    'use strict';
    $.each(productList, function () {
        if (this.shortname === 'simulink') {
            loadExamples(this.shortname);
            return false;
        }
    });
}

function loadExamples(product) {
    'use strict';
    var landingPageExampleJsonUrl = product + "_featured_examples.json";
    $.getJSON(landingPageExampleJsonUrl, function (data) {
        var jsonData = data; 
        var productJsonUrl = product + "/examples" + getJsonSuffix() + ".json";
        $.getJSON(productJsonUrl, function (productData) {
                var augmentedExamplesList = [];
                $.each(productData, function(index, example) {
                    if ($.inArray(example.open_command, jsonData.examples) >= 0 ) {
                        example.url = product + "/" +  example.url;
                        example.thumbnail_url = example.thumbnail_url.replace(/^\.\.\//, '');
                        augmentedExamplesList.push(example);
                    }
                });
                jsonData.examples = augmentedExamplesList;
                $('#' + product + '-examples').append(JST.example_cards(jsonData));
                triggerExamplesCardsAdded();
        });
                
        
    });
}

function getDescription(description) {
    'use strict';
    var regex_array = [
        "This example shows how to ",
        "This example shows you how to ",
        "This example shows how you can ",
        "This example shows that ",
        "This example shows ",
        "This model shows ",
        "This is an example of ",
        "This example is ",
        "This example "
    ],
        regex = new RegExp("^(" + regex_array.join("|") + ")"),
        result = description.replace(regex, "");
    return result.charAt(0).toUpperCase() + result.slice(1);
}

function getJsonSuffix() {
   'use strict';
    var pageLanguage = getPageLanguage();
    var jsonSuffix = "";
    if (pageLanguage === undefined || pageLanguage === "en") {
        jsonSuffix = "";
    } else {
        jsonSuffix = "_" + pageLanguage;
    }
    return jsonSuffix;
}

function getTranslateOpenExampleString(string) {
    string = string || "Example";
    return getLocalizedString("Open " + string);
}

function getTranslateViewMoreExampleString(string) {
    return getLocalizedString("View more " + string + " Examples");
}