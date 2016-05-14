window.JST = window.JST || {};
var suggestionsObj;
$(function() {
    var searchField = $('#docsearch');
    suggestionsObj = new Suggestions(getSuggestions, getMore, searchField);
});

function getSuggestions() {
    var text = suggestionsObj.searchField.val();
    document.location='docsuggestion:' + text;
}

function getMore(type) {
    var text = suggestionsObj.searchField.val();
    document.location='docsuggestion:' + type + ':' + text;
}

function Suggestions(retrievalFunction, moreFunction, searchField) {

    JST['pagesuggestion'] = _.template(
		'<% var docroot = getDocRoot(); %>' +
		'<% _.each(suggestions, function (suggestion) { %>' +
		'<div tabindex="-1" id="<%= suggestion.type %>:<%= suggestion.shortName %>:<%= suggestion.pageId %>"' +
			 'class="suggestion <%= suggestion.type%>_suggestion" action="go(\'<%= docroot + suggestion.path %>\')">' +
			 '<span class="suggestiontext">' +
				 '<% _.each(suggestion.title, function(value, index) { %><%= displaySuggestionText(index, value) %><% }) %>' +
				 ' - ' +
				 '<span class="suggestionsummary">' +
					 '<% _.each(suggestion.summary, function(value, index) { %><%= displaySuggestionText(index, value) %><% }) %>' +
				 '</span>' +
			 '</span>' +
		 '<span class="suggestionprod">' +
			 '<%= suggestion.product %>' +
		 '</span>' +
		'</div>' +
		'<% }) %>'
    )


    JST['wordsuggestion'] = _.template(
        '<% _.each(wordlist, function(word) { %>' +
        '<div tabindex="-1" class="suggestion" data-searchword="<%= word.fullword %>">' +
        '<% _.each(word.splitwordlist, function(value, index) { %><%= displaySuggestionText(index, value) %><% }) %>' +
        '</div>' +
        '<% }) %>'
    )

    JST['morepages'] = _.template(
        '<div id="more<%= type %>" tabindex="-1" class="suggestion suggestion-more" data-pagetype="<%= type %>">' +
            '<span class="suggestmore">' +
            '<%= more %> more' +
            '</span>' +
            '</div>'
    )

    JST['morewords'] = _.template(
        '<div id="moreword" tabindex="-1" class="suggestion suggestion-more" data-pagetype="word"">' +
            '<span class="suggestmore">' +
            '<%= more %> more' +
            '</span>' +
            '</div>'
    );


    JST['suggestionsdropdown'] = _.template(
		'<div>' +
			'<% _.each(pages, function(pageGroup) { %>' +
			'<div class="suggestionheader"><%= pageGroup.header %></div>' +
			'<div id="<%= pageGroup.type %>area" class="suggestionarea" data-pagetype="<%= pageGroup.type %>">' +
				'<%= pageSuggestionTemplate({suggestions: pageGroup.suggestions, displaySuggestionText:displaySuggestionText})' +
				'%>' +
				'<% if (pageGroup["more"] && pageGroup["more"] > 0) { %>' +
				'<%= morePagesTemplate(pageGroup) %>' +
				'<% } %>' +
			'</div>' +
			'<% }) %>' +
			'<% var wordGroup = words; %>' +
			'<% if (wordGroup && wordGroup.wordlist.length > 0) { %>' +
			'<div class="suggestionheader">Search Suggestions</div>' +
			'<div id="wordarea" class="suggestionarea">' +
				'<%= wordSuggestionTemplate({wordlist: wordGroup.wordlist, displaySuggestionText: displaySuggestionText})' +
				'%>' +
				'<% if (wordGroup["more"] && wordGroup["more"] > 0) { %>' +
				'<%= moreWordsTemplate(wordGroup) %>' +
				'<% } %>' +
			'</div>' +
			'<% } %>' +
		'</div>'
    );


    if (!searchField) {
        // Try to find the search field if one wasn't provided.
        searchField = $('#docsearch');
        if (searchField.length == 0) {
            searchField = $('#searchfield');
        }
    }
    this.searchField = searchField;
    this.last = '';
    this.suggestions = {};
    var suggestionsObj = this;
    this.retrieveSuggestions = function() {
        var text = suggestionsObj.searchField.val();
        document.location = 'docsuggestion:' + encodeURIComponent(text);
    };

    this.retrieveMore = function(type) {
        var text = suggestionsObj.searchField.val();
        document.location = 'docsuggestion:' + type + '|' + encodeURIComponent(text);
    };
    this.selectionMode = 'keyboard';

    searchField.keyup(function(evt) {
        evt.stopPropagation();
        validateSearchForm($(this));
        // Show suggestions whenever the search text has changed or the user hits the down arrow.
        if (evt.keyCode == 40 || suggestionsObj.isSearchFieldChanged()) {
            suggestionsObj.last = searchField.val();
            suggestionsObj.findSuggestions();
        }
    });

    // Hide the suggestions for any click in the document, except clicks in the search field.  Later we will also
    // exempt clicks in the suggestion area itself.
    $(document).click(function(e) {
        if (e.button !== 0) {
            return true;
        }
        suggestionsObj.hidePopups(false);
    });
    searchField.click(function(evt) {
        validateSearchForm($(this));
        evt.stopPropagation();
    });

    $("#submitsearch").prop('disabled', true);

    $(window).resize(function() {
        var tokenizedDiv = searchField.closest('.tokenized');

        //Resize the div only if it has been tokenized, and the nuggets have been applied.
        if (tokenizedDiv.find('.tokens:first').children().length > 0) {
            tokenizedDiv.trigger('token.resize');
        }
    });
}


Suggestions.prototype.findSuggestions = function() {
    var searchText = this.searchField.val();
    if (searchText.match(/^(\w{2,}\s*)+$/) || searchText.match(/[^\w\s]+/)) {
        this.retrieveSuggestions();
    } else {
        this.hidePopups(false);
    }
};

Suggestions.prototype.isSearchFieldChanged = function() {
    return this.last != this.searchField.val();
};

Suggestions.prototype.displaySuggestions = function(suggestions) {
    var searchText = suggestions.searchtext;
    var currentSearchText = this.searchField.val();
    if (searchText === currentSearchText) {
        this.suggestions = suggestions;
        var html = this.formatSuggestions();
        if (html.length > 0) {
            this.displaySuggestionHtml(html);
        } else {
            this.hidePopups(false)
        }
    }
};

Suggestions.prototype.goToSelectedPage = function(group) {
    var index = $('.selected-suggestion').index();
    var selectedGroup = this.findGroup(group);
    var path = selectedGroup.suggestions[index].path;
    go(path);
};

Suggestions.prototype.findGroup = function(pageType) {
    for (var i = 0; i < this.suggestions.pages.length; i++) {
        var group = this.suggestions.pages[i];
        if (group.type === pageType) {
            return group;
        }
    }
    return null;
}

Suggestions.prototype.formatSuggestions = function() {
    var suggestionsObj = this;


    //if there are not suggestions, return empty.
    if (!((this.suggestions.pages && this.suggestions.pages.length > 0) || this.suggestions.words)) {
        return '';
    }

    var jsonData = _.extend({}, this.suggestions, {
        pageSuggestionTemplate: JST['pagesuggestion'],
        wordSuggestionTemplate: JST['wordsuggestion'],
        morePagesTemplate: JST['morepages'],
        moreWordsTemplate: JST['morewords']
    });




    var suggestionsElt = $(JST['suggestionsdropdown'](jsonData));
    //attach the event handlers.
    suggestionsElt.find('#wordarea').on('selection-action', function () {
            suggestionsObj.search($('.selected-suggestion').data('searchword'));
        });

    suggestionsElt.find('.suggestion-more').on('selection-action', function (evt) {
        evt.stopPropagation();
        suggestionsObj.more($(this).data('pagetype'));
    });

    suggestionsElt.find('.suggestionarea:not([wordarea])').on('selection-action', function () {
        suggestionsObj.goToSelectedPage($(this).data('pagetype'));
    });
    return suggestionsElt;

};

Suggestions.prototype.displaySuggestionHtml = function(html) {
    // TODO: is the logic here correct?  I feel like there's some asynchronous stuff happening that we're not taking
    // into account.
    if (!this.isSearchFieldChanged()) {
        var suggestionsElt = this.getCleanSuggestionsElement();
        if (html.length > 0) {
            suggestionsElt.append(html);
            if (suggestionsElt.is(':hidden')) {
                this.showSuggestions();
            }
            suggestionsElt.attr('suggestfor', this.searchField.val());
        } else {
            this.hidePopups(false);
        }
    }
};

Suggestions.prototype.getCleanSuggestionsElement = function() {
    var suggestionsElt = $('#suggestions');
    this.selectionMode = 'keyboard';
    if (suggestionsElt.length == 0) {
        suggestionsElt = $('<div id="suggestions"></div>');
        this.getSearchForm().after(suggestionsElt);
        this.addSuggestionHandlers(suggestionsElt);
    } else {
        suggestionsElt.empty();
    }

    var width = getPopupWidth(this.searchField);
    suggestionsElt.width(width);
    return suggestionsElt;
};

Suggestions.prototype.addSuggestionHandlers = function(suggestionsElt) {
    var suggestionsObj = this;
    suggestionsElt.keydown(function(evt) {
        suggestionsObj.selectionMode = 'keyboard';
        return suggestionsObj.handleKeyDown(evt);
    });
    suggestionsElt.keyup(function(evt) {
        suggestionsObj.selectionMode = 'keyboard';
        if (evt.keyCode == 13) {
            suggestionsObj.handleSelection();
        }
    });
    suggestionsElt.mousemove(function() {
        suggestionsObj.selectionMode = 'mouse';
    });
    suggestionsElt.mouseover(function(evt) {
        if (suggestionsObj.selectionMode === 'mouse') {
            selectFromMouseEvent(evt);
        }
    });
    suggestionsElt.click(function(evt) {
        evt.stopPropagation();
        selectFromMouseEvent(evt);
        suggestionsObj.handleSelection();
    });
};

Suggestions.prototype.showSuggestions = function() {
    var suggestionsElt = $('#suggestions');
    suggestionsElt.slideDown('fast');
    this.searchField.unbind('keydown.searchfield-suggestions');

    var suggestionsObj = this;
    this.searchField.bind('keydown.searchfield-suggestions', function(evt) {
        return suggestionsObj.handleKeyDown(evt);
    });
};

Suggestions.prototype.hidePopups = function(focusSearchField) {
    this.suggestions = {};
    $('#suggestions').remove();
    this.searchField.unbind('keydown.searchfield-suggestions');
    if (focusSearchField) {
        this.searchField.get(0).focus();
    }
};

Suggestions.prototype.handleKeyDown = function(evt) {
    var key = evt.keyCode? evt.keyCode : evt.charCode;
    if (key == 40) {
        // Down arrow
        this.selectNextSuggestion();
        return false;
    } else if (key == 38) {
        // Up arrow
        this.selectPrevSuggestion();
        return false;
    } else if (key == 27) {
        // Escape
        this.hidePopups(true);
        return false;
    } else {
        return true;
    }
};

Suggestions.prototype.handleSelection = function() {
    var triggerElt = $('.selected-suggestion');
    if (!triggerElt.hasClass('suggestion-more')) {
        triggerElt = triggerElt.closest('.suggestionarea');
    }
    triggerElt.trigger('selection-action');
};

Suggestions.prototype.selectNextSuggestion = function() {
    var newselection = getNextSelection();
    if (newselection != null && newselection.length > 0) {
        var oldselection = $('.selected-suggestion');
        this.changeSelection(oldselection, newselection);
    }
};

Suggestions.prototype.selectPrevSuggestion = function() {
    var newselection = getPreviousSelection();
    var oldselection = $('.selected-suggestion');
    if (newselection != null && newselection.length > 0) {
        this.changeSelection(oldselection, newselection);
    } else {
        oldselection.removeClass('selected-suggestion');
        this.searchField.get(0).focus();
    }
};

Suggestions.prototype.changeSelection = function(oldselection, newselection) {
    oldselection.removeClass('selected-suggestion');
    newselection.addClass('selected-suggestion');
    this.searchField.blur();
    newselection.get(0).focus();
};

Suggestions.prototype.getSearchForm = function() {
    return this.searchField.closest('form');
};

Suggestions.prototype.search = function(q) {
    this.searchField.val(q);
    this.getSearchForm().submit();
    this.hidePopups(false);
};

Suggestions.prototype.getSuggestionsUrl = function(action, params) {
    var formAction = this.getSearchForm().attr('action');
    params.width = this.searchField.width();
    return formAction.replace(/\/search(\/|$)/, '/search/' + action + '/') + '?' + $.param(params);
};

Suggestions.prototype.more = function(type) {
    var divId = type + 'area';
    var div = $('#' + divId);
    div.css('minHeight', div.outerHeight());
    div.css('maxHeight', div.outerHeight());

    var moreElt = $('#more' + type);
    moreElt.removeClass('suggestmore');
    moreElt.removeClass('suggestion');
    moreElt.removeClass('selected-suggestion');

    this.retrieveMore(type);
};

Suggestions.prototype.displayMoreSuggestions = function(type, template, suggestions) {
    $('#more' + type).remove();
	if (type !== 'word') {
        var toAdd = suggestions.suggestions;
        var group = this.findGroup(type);
        var suggestionsArray = group.suggestions;
        suggestionsArray = suggestionsArray.concat(toAdd);
        group.suggestions = suggestionsArray;
    }

    var templateToCall = JST[template];
    var suggestionElts = $(templateToCall(suggestions));
    suggestionElts.eq(0).addClass('selected-suggestion');
    suggestionElts.appendTo($('#' + type + 'area'));
    $('.selected-suggestion').focus();
};

/* ***** Utility functions start here ***** */
function selectFromMouseEvent(evt) {
    var newSelection = $(evt.target).closest('.suggestion');
    if (!newSelection.hasClass('selected-suggestion')) {
        var oldSelection = $('.selected-suggestion');
        oldSelection.removeClass('selected-suggestion');
        newSelection.addClass('selected-suggestion');
    }
}

function getNextSelection() {
    var selected = $('.selected-suggestion');
    var newselection = null;
    if (selected.length > 0) {
        newselection = selected.nextAll('.suggestion:first');
        if (newselection.length == 0) {
            // We're at the end of a section, jump to the next.
            var area = selected.closest('.suggestionarea').nextAll('.suggestionarea:first');
            return area.children('.suggestion:first');
        } else {
            return newselection;
        }
    } else {
        return $('.suggestion:first');
    }
}

function getPreviousSelection() {
    var selected = $('.selected-suggestion');
    var newselection = null;
    if (selected.length > 0) {
        newselection = selected.prevAll('.suggestion:first');
        if (newselection.length == 0) {
            // We're at the beginning of a selection, jump to the previous.
            var area = selected.closest('.suggestionarea').prevAll('.suggestionarea:first');
            return area.children('.suggestion:last');
        } else {
            return newselection;
        }
    } else {
        return null;
    }
}

function go(page) {
    document.location = page;
}

function getSectionHeader(title) {
    var header = $('<div class="suggestionheader"></div>');
    header.append(title);
    return header;
}

function getSectionDiv(type) {
    return $('<div id="' + type + 'area" class="suggestionarea"></div>');
}

function reportSuggestion(suiteName, suggestionsElt) {
    try {
        var suggId = suggestionsElt.attr('id');
        var parts = suggId.split(':');
        var type = parts[0];
        var searchTerm = parts[2];
        var elt = suggestionsElt.get(0);
        var s = s_gi(suiteName);
        s.linkTrackVars='prop3,eVar3,events';
        s.linkTrackEvents='event35';
        s.eVar3=searchTerm + '::lucene::' + type;
        s.prop3=searchTerm + '::lucene::' + type;
        s.events='event35';
        s.tl(elt,'o','Doc Suggestion');
    } catch (e) {
        // Don't report the suggestion.
    }
}

function validateSearchForm(searchField) {
    var val = searchField.val();
    if (val.length === 0) {
        $("#submitsearch").prop('disabled', true);
    } else {
        $("#submitsearch").prop('disabled', false);
    }
}

function getPopupWidth(searchField) {
    var width = searchField.outerWidth();
    var parentDiv = searchField.closest('.search_field');
    if (parentDiv.length > 0) {
        width = parentDiv.outerWidth();
    }
    return width;
}

function displaySuggestionText(index, value) {
    if (value.length > 0) {
        if (index % 2 === 0) {
            return value;
        } else {
            return "<span class='suggestionhighlight'>" + value + "</span>";
        }
    } else {
        return "";
    }
}

function getDocRoot() {
    // In the product we always use full paths for suggestions.
	return "";
}