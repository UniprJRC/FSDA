window.JST = window.JST || {};

function init() {
    window.rnLocaleSuffix = getLocaleSuffix();
    handleReleaseNotes();
    $(window).on('hashchange', checkHash);
    $(document.forms.rnform).submit(function() {
        if (document.forms.rnform.searchHighlight) {
            document.forms.rnform.searchHighlight.value = document.forms.rnform.rntext.value;
        }
    });
}

function handleReleaseNotes() {
    id = parseHash();
    if (id && id.length > 0) {
        handleId(id);
    } else {
        var params = parseParams();
        populateDefaults(params);
        populateForm(params);
        populateSort(params);
        expandReleases(params);
        addReleaseRange(params);
        addNuggets(params);
        getAllReleaseNotes(params);
    }
}

function resetNotes() {
    document.location = document.location.pathname;
}

function checkHash() {
    var id = parseHash();
    if ($("#" + id).length == 0) {
        handleId(id);
    }
}

function parseHash() {
    var hash = window.location.hash;
    if (hash && hash.length > 0) {
        return hash.replace(/^#/,"");
    } else {
        return "";
    }
}

function handleId(id) {
    var params = {"id":id};
    populateDefaults(params);
    getAllReleaseNotes(params);
}

function parseParams() {
	var params = {};
    var qs = window.location.search;
    if (qs && qs.length > 0) {
        var paramsArray = qs.replace(/^\?/,"").split("&");
        for (var i = 0; i < paramsArray.length; i++) {
            var nameValPair = paramsArray[i].split("=");
            var name = nameValPair[0];
            var value = nameValPair.length > 1 ? nameValPair[1] : "";
            if (name === 'category') {
                value = '"' + value + '"';
            }
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
    
    if (params.rntext) {
        params.searchHighlight = params.rntext;
    } else if (params.searchHighlight) {
        params.rntext = params.searchHighlight;
    }

    return params;
}

function populateDefaults(params) {
    params.product = shortname;
    if (!params.groupby) {
        params.groupby = "release";
        params.sortby = "descending";
    }
}

function populateSort(params) {
    var sortvalue = params.groupby;
    if (params.sortby) {
        sortvalue += "-" + params.sortby;
    }
    var sortSelect = document.getElementById('selectsort');
    for (var i = 0; i < sortSelect.options.length; i++) {
    	var sortOption = sortSelect.options[i];
        if (sortOption.value === sortvalue) {
            sortOption.selected = true;
            break;
        }
    }
}

function handleSort(sortOption) {
    var optionParts = sortOption.split(/-/);
    var searchForm = document.forms.rnform;
    searchForm.groupby.value = optionParts[0];
    searchForm.sortby.value = optionParts[1];
    searchForm.submit();
}

function populateForm(params) {
    populateReleaseFields(params);
    populateTextFields(params, ["rntext","groupby","sortby","searchHighlight"]);
    document.forms.rnform.searchHighlight.value = document.forms.rnform.rntext.value;

    var searchForm = document.forms.rnform;
    searchForm.rntype.checked = params.rntype && params.rntype === "incompatibility";
    
    if (params.category) {
        var categories = getCategoriesFromParams(params);
        $("input[name='category']").each(function(i,checkboxElt) {
            checkboxElt.checked = categories.indexOf(checkboxElt.value) >= 0;
        });
    }
}

function getCategoriesFromParams(params) {
    var catString = params.category;
    catString = catString.replace(/^"(.*)"$/, "$1");
    return catString.split(/","/);
}

function populateTextFields(params, names) {
    var searchForm = document.forms.rnform;
    for (var i = 0; i < names.length; i++) {
        var fieldName = names[i];
        searchForm[fieldName].value = params[fieldName] ? params[fieldName] : "";
    }
}

function populateReleaseFields(params) {
    findSelectedReleases(params);
    var searchForm = document.forms.rnform;
    populateReleaseOptions(searchForm.startrelease, params.startrelease);
    populateReleaseOptions(searchForm.endrelease, params.endrelease);
}

function findSelectedReleases(params) {
    if (!params.endrelease) {
        params.endrelease = allRnReleases[allRnReleases.length-1];
    }
    
    if (!params.startrelease) {
        if (params.searchHighlight) {
            params.startrelease = allRnReleases[0];
        } else {
            var endindex = allRnReleases.indexOf(params.endrelease);
            var startindex = Math.max(0,endindex - 6);
            params.startrelease = allRnReleases[startindex];
        }
    }
}

function populateReleaseOptions(selectElt, selectedValue) {
    var found = false;
    for (var i = allRnReleases.length-1; i >= 0; i--) {
        var release = allRnReleases[i];
        var startOption = $("<option>").html(release);
        if (release === selectedValue) {
            startOption.attr("selected","selected");
            found = true;
        }
        $(selectElt).append(startOption);
    }
    
    if (!found) {
        $(selectElt).children().first().attr("selected","selected");
    }
}

function getAllReleaseNotes(params) {
    var d1 = getReleaseNotes(params);
    var d2 = getProductUpdates();

    $.when(d1,d2).done(function (v1,v2) {
        addReleaseNotes(v1,v2);
    }).fail(function (jqXHR, textStatus, error) {
        if (d1.state() === 'rejected') {
            showAllNotesForReleases(params.releases);
        } else if (d2.state() === 'rejected') {
            $.when(d1).done(function (v1) {
                addReleaseNotes(v1);
            }).fail(function (jqXHR, textStatus, error) {
                showAllNotesForReleases(params.releases);
            });
        }
    });
}

function getReleaseNotes(params) {
    $("#notes").empty();
    var deferred = $.Deferred(function() {});
    var services = {
        "messagechannel":"releasenotes",
        "requesthandler":"releasenotes://search",
        "webservice":getWebServiceUrl()
    }
    var errorhandler = function() {
        deferred.reject();
    }
    var successhandler = function(data) {
        deferred.resolve(data);
    }
    requestHelpService(params, services, successhandler, errorhandler);
    return deferred;
}

function getWebServiceUrl() {
    var lang = getPageLanguage() || "en";

    var release = getDocReleaseFromSearchBox();
    if (typeof getDocRelease === 'function') {
        release = getDocRelease();
    }
    
    return "/help/search/releasenotes/doccenter/" + lang + "/" + release;
}

function getProductUpdates() {
    var deferred = $.Deferred(function() {});
    var params = {};
    var url = "/support/bugreports/updates/" + basecode;
    var services = {
        "messagechannel":"productupdates",
        "webservice":url
    }
    var errorhandler = function() {
        deferred.reject();
    }
    var successhandler = function(data) {
        deferred.resolve(data);
    }
    requestHelpService(params, services, successhandler, errorhandler);
    return deferred;
}

function showAllNotesForReleases(releases) {
    var warning = '<div class="alert alert-warning" style="margin-top:10px">' +
        '<span class="alert_icon icon-alert-warning"></span>' +
        '<p>Could not retrieve release notes from server. Showing all release notes for the selected releases. ' +
        'Some functionality on this page might be unavailable.</p>' +
        '</div>';
    $("#notes").append(warning);
    getRnPage(releases.split(",").reverse());
}

function getRnPage(releases) {
    var onSuccess = function(data, release) {
        content = $(data).find("#doc_center_content");
        $("#notes").append(content.children());
    }
    var onComplete = function() {
        renderEquations();
    }
    
    getReleaseContent(releases, false, onSuccess, onComplete);
}

function expandReleases(params) {
    var searchForm = document.forms.rnform;
    var startElt = searchForm.startrelease;
    var endElt = searchForm.endrelease;
    
    var selectedReleases = [startElt.options[startElt.selectedIndex].value, endElt.options[endElt.selectedIndex].value];
    selectedReleases.sort(); // Releases sort alphabetically just fine.
    var startIndex = allRnReleases.indexOf(selectedReleases[0]);
    var endIndex = allRnReleases.indexOf(selectedReleases[1])+1;
	params.releases = allRnReleases.slice(startIndex, endIndex).join(",");
}

function addReleaseNotes(releaseNotes, productUpdates) {
    if (releaseNotes.filter) {
        populateReleaseFields(releaseNotes.filter);
        addReleaseRange(releaseNotes.filter);
    }
    storeHighlightExpand(document.forms.rnform.rntext.value, releaseNotes.highlightexpand);

    var totalNotes = 0;
	var notesDiv = $("#notes");
    var groups = releaseNotes.rngroups;
	for (var i = 0; i < groups.length; i++) {
		var groupObj = groups[i];
        var contentDiv = $('<div class="expandableContent"></div>');
        notesDiv.append(contentDiv);

		var hasContent = groupObj.notes || groupObj.rngroups;
        var expandedClass = hasContent ? "expanded" : "no_content expanded";
        var headerHtml = '<span id="' + groupObj.label + '" class="anchor_target"></span><div class="expand ' + expandedClass + '">' +
                         '<h2 id="' + groupObj.uniqueid + '">' + groupObj.label + '</h2>';
        var desc = groupObj.description ? groupObj.description : "&nbsp;";
        headerHtml += '<div class="doc_topic_desc">' + desc + '</div>';
        if (hasContent) {
            headerHtml += '<div class="switch"><a class="expandAllLink" href="javascript:void(0);"> expand all </a></div></div>';
            if (productUpdates) {
                var productUpdatesForRelease = productUpdates[groupObj.label];
                if (productUpdatesForRelease) {
                    var productUpdateMarkup = getProductUpdateMarkup(productUpdatesForRelease);
                    if (productUpdateMarkup && productUpdateMarkup.length > 0) {
                        headerHtml += productUpdateMarkup;  
                    }
                }                
            }
        }
		
        contentDiv.append(headerHtml);
        if (groupObj.notes) {
            var collapseDiv = $('<div class="collapse" style="display:block;"></div>');
            appendNotes(collapseDiv, groupObj.notes);
            contentDiv.append(collapseDiv);
            totalNotes += groupObj.notes.length;
		}
		
        var subGroups = groupObj.rngroups;
		if (subGroups) {
			for (var j = 0; j < subGroups.length; j++) {
				var subGroupObj = subGroups[j];
                collapseDiv = $('<div class="collapse" style="display:block;"></div>');
				collapseDiv.append('<h3 id="' + subGroupObj.uniqueid + '">' + subGroupObj.label + '</h3>');
                contentDiv.append(collapseDiv);
				appendNotes(collapseDiv, subGroupObj.notes);
                totalNotes += subGroupObj.notes.length;
			}
		}
	}
    
    $("#num-notes").html(totalNotes);
    $("#num-notes-container").show();
	
    var releases = Object.getOwnPropertyNames(releaseNotes.releases);
    
    var onSuccess = function(data, release) {
        var ids = releaseNotes.releases[release];
        handleRelease(data, ids, release);
    }
    
    var onComplete = function() {
        applySearchHighlight();
        renderEquations();
        addSmoothScroll();
        expandCollapsedContent(true);
    }
    
    getReleaseContent(releases, false, onSuccess, onComplete);
}

function getProductUpdateMarkup(productUpdatesForRelease) {
    if (productUpdatesForRelease) {
        var localizedText = getLocalizedText('bug_fixes', 'Bug Fixes');
        var updatesForRelease = [];
        var keys1 = Object.keys(productUpdatesForRelease);

        for (var ctr1 = 0; ctr1 < keys1.length; ctr1++) {
            var key1 = keys1[ctr1];
            var entry = productUpdatesForRelease[key1];

            var keys2 = Object.keys(entry);
            for (var ctr2 = 0; ctr2 < keys2.length; ctr2++) {
                updatesForRelease.push({"label":keys2[ctr2], "url":entry[keys2[ctr2]], "localizedtext":localizedText});
            }
        }
        var compiledTmpl = JST['updatesTmpl']({updates: updatesForRelease});
        return compiledTmpl;
    }
    return '';
}

JST['updatesTmpl'] = _.template(
       '<div class="collapse" style="display:block;">' +
       '<ul class="list-unstyled">' +
       '<% _.each(updates, function(update) { %>' +
         '<li><%= update.label %>: <a href="<%= update.url %>"><%= update.localizedtext %></a></li>' +
       '<% }); %>' +
       '</ul>' +
       '</div>'
);

function getLocalizedText(l10nKey, defaultString) {
    var localizedString = defaultString;
    if (typeof getLocalizedString !== "undefined") {
        localizedString = getLocalizedString(l10nKey);
    } 
    return localizedString;
}

function storeHighlightExpand(searchTerm, highlightExpand) {
    if (searchTerm && highlightExpand) {
        if (localStorage && localStorage.setItem) {
            localStorage.setItem(searchTerm, JSON.stringify(highlightExpand));    
        }
    }
}

function getLocaleSuffix() {
    var localeRegexp = /release-notes(_.*)\.html/;
    var localeMatch = localeRegexp.exec(window.location.pathname);
    if (localeMatch) {
        return localeMatch[1];
    } else {
        return "";
    }
}

// Retrieve the content for each release to display.
//   releases: an array containing all of the releases to retrieve.
//   fallback: a boolean indicating whether to fall back to English (used for installed
//             translated doc only.
//   onSuccess: a function called each time a release's content is retrieved.
//   onComplete: a function called when all release not content has been retrieved.
function getReleaseContent(releases, fallback, onSuccess, onComplete) {
    if (releases.length === 0) {
        onComplete();
        return;
    }
    
    var release = releases[0];
    var localeSuffix = fallback ? "" : window.rnLocaleSuffix;
    var url = "release-notes-" + release + localeSuffix + ".html";
    var jqxhr = $.get(url);
    jqxhr.done(function(data) {
        onSuccess(data, release);
        releases.shift();
        getReleaseContent(releases, false, onSuccess, onComplete);
    });
    jqxhr.fail(function() {
        if (localeSuffix.length > 0) {
            // Attempt to fall back to English.
            getReleaseContent(releases, true, onSuccess, onComplete);
        } else {
            releases.shift();
            getReleaseContent(releases, false, onSuccess, onComplete);
        }
    });
}

function renderEquations() {
    // These are the configuration parameters which are passed through to the
    // StandaloneEqnRenderer by default.
    var equationrendererDefaultConfig = {
            flavor: "MathType",
            equationFormat: "mathml",
            equationEncoding: "element",
            equationRootElement: "math",
            cacheFontMetrics: false
        };
        
    require([
        "MW/equations/renderer/StandaloneEqnRenderer", "dojo/domReady!"
    ], function (StandaloneEqnRenderer) {
        var ren = new StandaloneEqnRenderer(equationrendererDefaultConfig);
        ren.render();
    });
}

function appendNotes(contentDiv, notes) {
	for (var i = 0; i < notes.length; i++) {
		var note = notes[i];
        var containerDiv = $('<div class="expandableContent" id="container-' + note.id + '"></div>');
        var expandDiv = $('<div class="expand"></div>');
        var titleDiv = $('<h4 id="' + note.id + '">' + note.title + '</h4>');
        expandDiv.append(titleDiv);
        containerDiv.append(expandDiv);
        containerDiv.append('<div class="collapse" id="content-' + note.id + '"></div>');
        contentDiv.append(containerDiv);
	}
}

function handleRelease(data, ids, release) {
	var rn = $(data);
	for (var i = 0; i < ids.length; i++) {
		var id = ids[i];
		var noteContents = rn.find("#" + id).parent().next();
        if (noteContents && noteContents.length > 0) {
            $("#content-" + id).append(noteContents[0].innerHTML);
        }
	}
}

function addReleaseRange(params) {
    $('#start-release').html(params.startrelease);
    $('#end-release').html(params.endrelease);
}

function addNuggets(params) {
    var nuggetsDiv = $("#nugget-container");
    var numNuggets = 0;
    if (params.rntype && params.rntype === "incompatibility") {
        nuggetsDiv.show();
        // The incompatibility nugget is hard-coded on the page to avoid i18n issues.
        $('#nugget-incompat').show();
        $('#nugget-incompat').on("click", removeIncompatibilities);
        numNuggets++;
    }
    
    if (params.category) {
        var categories = getCategoriesFromParams(params);
        for (var i = 0; i < categories.length; i++) {
            var category = categories[i];
            addNugget(category, removeCategory);
            numNuggets++;
        }
    }
    
    if (params.rntext && params.rntext.length > 0) {
        addNugget(params.rntext, removeText);
        numNuggets++;
    }
    
    if (numNuggets > 1) {
        var nuggetSpan = nuggetsDiv.find("span[class='nuggets']");
        var removeAllSpan = $('<span class="nugget nugget_remove_all"><a href="javascript:void(0);"><span class="label">Remove All</span></a></span>');
        removeAllSpan.on("click", removeAllNuggets);
        nuggetSpan.append(removeAllSpan);
    }
}

function addNugget(label, removeFcn) {
    var nuggetsDiv = $("#nugget-container");
    nuggetsDiv.show();
    var nuggetSpan = nuggetsDiv.find("span[class='nuggets']");
    var newSpan = $('<span class="nugget"></span>');
    var newNugget = $('<a href="javascript:void(0);" class="icon-remove"><span class="label">' + label + "</span></a>");
    newNugget.on("click", {"label":label}, removeFcn);
    newSpan.append(newNugget);
    nuggetSpan.append(newSpan);
}

function removeCategory(event) {
    var checkbox = $("input[name='category'][value='" + event.data.label + "']");
    checkbox.attr('checked',false);
    document.forms.rnform.submit();
}

function removeIncompatibilities(event) {
    document.forms.rnform.rntype.checked = false;
    document.forms.rnform.submit();
}

function removeText(event) {
    document.forms.rnform.rntext.value = "";
    document.forms.rnform.searchHighlight.value = "";
    document.forms.rnform.submit();
}

function removeAllNuggets(event) {
    $("input[name='category']").attr("checked",false);
    var rnForm = document.forms.rnform;
    rnForm.rntype.checked = false;
    rnForm.rntext.value = "";
    rnForm.searchHighlight.value = "";
    rnForm.submit();
}

$(init);
