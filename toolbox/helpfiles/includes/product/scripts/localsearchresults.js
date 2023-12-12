var get_facet_label;
function populateResults(jsonObject) {
    var facetLabels = jsonObject.facetlabels ? jsonObject.facetlabels : {};
    populateLookupTable(facetLabels);
    
    if (jsonObject.responsetype === "noresults") {
        displayMessage(jsonObject);
        displaySpellCheck(jsonObject);
    } else if (jsonObject.responsetype === "error") {
        displayError(jsonObject);
    } else {
        populateResultsList(jsonObject.results, jsonObject.searchtext);
        populateResultData(jsonObject.pagedata);
        populateFacets(jsonObject.facets);
        displaySpellCheck(jsonObject);
        appendHighlightExpand(jsonObject.highlightexpand, jsonObject.searchtext);
    }
}

function populateResultsList(searchresults, searchTerm) {
    var highlightTerm = "";
    if (searchTerm && searchTerm.length > 0) {
        highlightTerm = searchTerm;
    } else if( $('#docsearch') && $('#docsearch').val() ){
        highlightTerm = $('#docsearch').val();
    }
    $('#wait').remove();

    var resultsHtml = '';
    resultsHtml = getSearchResultsHtml(searchresults, highlightTerm);

    var resultsDiv = $('#results_area');
    resultsDiv.html(resultsHtml);
}

function populateResultData(jsonObject) {
    var searchterm = jsonObject.searchterm;
    var productbreadcrumb = jsonObject.productbreadcrumb;
    var summarydata = jsonObject.summarydata;
    var footerdata = jsonObject.footerdata;
    var searchTextDesc = jsonObject.searchtext;

    $('#docsearch').val(searchterm);
    tokenizeSearchText();

    var summaryHtml = '';
    summaryHtml = getSearchSummaryHtml(footerdata);

    var searchingInfoDiv = $('#search_result_info_header');
    searchingInfoDiv.html(searchTextDesc);   

    var summaryDiv = $('#search_result_header');
    summaryDiv.html(summaryHtml);

    var footerHtml = '';
    footerHtml = getSearchFooterHtml(footerdata);

    var footerDiv = $('#search_result_footer');
    footerDiv.html(footerHtml);

    setPageTitle();
}

function populateFacets(facetJson) {
    var facetHtml = '';
    facetHtml = getFacetResultsHtml(facetJson);

    var facetDiv = $('#facets_area');
    facetDiv.html(facetHtml);
}

function displayError(error) {
    $('#docsearch').val(error.searchtext);
    var errorHtml = getErrorHtml(error.message);

    var errorDiv = $('#results_area');
    errorDiv.html(errorHtml);
    setPageTitle();
}

function displayMessage(message) {
    var messageDiv = $('#results_area');
    messageDiv.empty();
    $('#docsearch').val(message.searchtext);
    tokenizeSearchText();

    messageDiv.append(getSuggestionsListHtml(message));
    setPageTitle();
}

function setPageTitle() {
    document.title = "Search Results - " + $("#docsearch").val();
}

function tokenizeSearchText() {
    $('.tokenized_search_field').tokenize({
        fields: ["product", "category", "type"],
        remove_if_empty: true,
        label_function: get_facet_label
    });
}

function displaySpellCheck(jsonData) {
    if(jsonData === undefined) {
        return;
    }

    var spellcheckHtml = '';
    spellcheckHtml = getSpellCheckResultsHtml(jsonData);

    var messageDiv = $('#search_result_dym_header');
    messageDiv.html(spellcheckHtml);
}

function populateLookupTable(facetLabelJson) {

    var labels = {};
    for (var facetLabel in facetLabelJson) {
        if (facetLabelJson.hasOwnProperty(facetLabel)) {
            labels[facetLabelJson[facetLabel].field + ":" +
                    facetLabelJson[facetLabel].value] = facetLabelJson[facetLabel].label;
        }
    }
    get_facet_label = function (token) {
        var labelDataString = sessionStorage.getItem('facetlookuptable');
        if(labelDataString !== undefined) {
            var labelData = JSON.parse(labelDataString);
            return labelData[token.field + ":" + token.value];
        }
        return labels[token.field + ":" + token.value];
    }
    
    if(!$.isEmptyObject(labels)) {
        sessionStorage.setItem('facetlookuptable', JSON.stringify(labels));
    }
}

$(document).ready(function() {
  searchData = parseQueryString(); 
  var services = {"messagechannel":"docsearch"};
  requestHelpService(searchData, services, populateResults);
});


function parseQueryString() {
    var params = {};
    var qs = window.location.search;
    if (qs && qs.length > 0) {
        var paramsArray = qs.replace(/^\?/,"").split("&");
        for (var i = 0; i < paramsArray.length; i++) {
            var nameValPair = paramsArray[i].split("=");
            var name = nameValPair[0];
            var value = nameValPair.length > 1 ? nameValPair[1] : "";
            if (name && name.length > 0) {
                value = decodeURIComponent(value.replace(/\+/g," "));
                if (params[name]) {
                    params[name] += "," + value;
                } else {
                    params[name] = value;
                }
            }
        }
    }
    return params;
}