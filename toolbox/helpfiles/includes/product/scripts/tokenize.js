/*
 * tokenize -- turn an text input field into a tokenized widget
 */
(function ($) {
    var parse_query = function (query, fields) {
        var result = { query: "", tokens: [] },
            regexp, matchdata;

        $.each(fields, function (i, field) {
            regexp = new RegExp(field + ":(?:\"([^\"]+)\"|'([^']+)'|([^ ]+))", "i");

            while (matchdata = query.match(regexp)) {
                result.tokens.push({
                    field: field,
                    value: (matchdata[1] || matchdata[2] || matchdata[3])
                });
                query = query.replace(matchdata[0], "");
            }
        });
        result.query = $.trim(query);

        return result;
    };

    var create_token = function (parent, searchInput, options) {
        var token = $('<span class="nugget ' + options.field + '"/>');

        token.append($('<a href="#remove-token"class="icon-remove" title="'+ getLocalizedNuggetsString('remove_facet_nuggets', 'Remove: ') + options.label + '"><span class="label">' + options.label + '</span></a>'));

        token.appendTo(parent).data({
            field: options.field,
            value: options.value
        });

        token.find("a").bind("click", function (event) {
            var form = searchInput.closest("form");
            token.trigger("nugget.remove");
            form.submit();
            event.preventDefault();
        });

        token.find("a").hover(function (event) {
                $(event.target).closest('.nugget').addClass('nugget_pending_delete');
            },
            function (event) {
                $(event.target).closest('.nugget').removeClass('nugget_pending_delete');
            });
    };


    var tokens_to_query = function (tokens) {
        return $.map(tokens, function (token) {
            return $(token).data("field") + ':' + $(token).data("value");
        }).join(" ");
    };

    var getLocalizedNuggetsString = function(l10nKey, defaultString) {
        var remove_facet_nuggets = defaultString;
        if (typeof getLocalizedString !== "undefined") {
            remove_facet_nuggets = getLocalizedString(l10nKey);
        }
        return remove_facet_nuggets == undefined ? defaultString : remove_facet_nuggets + ' ';
    };

    $.fn.tokenize = function (options) {
        var input = this.find("input[type=text]"),
            hidden = $('<input type="hidden" name="' + input.attr("name") + '">'),
            nugget_box = options.nugget_box || $('#nugget-box'),
            $nugget_container = $('#nugget-container'),
            tokens = $('<span class="nuggets"/>').prependTo(nugget_box),
            parsed_query;

//      nugget_box.addClass("tokenized");

        // use a hidden input to submit the actual field value
        this.closest("form").append(hidden).bind("submit", function (event) {
            var term = input.val() + " " + tokens_to_query(nugget_box.find(".nugget").not('.nugget_remove_all'));
            term = $.trim(term);

            if (term == "" && options.remove_if_empty) {
                hidden.remove();
            } else {
                hidden.val(term);
            }
        });
        input.removeAttr("name");

        // extract tokens from the query
        parsed_query = parse_query(input.val(), options.fields || []);

        // update the text input
        input.val(parsed_query.query);


        // turn the parsed tokens into elements
        $.each(parsed_query.tokens, function (i, token) {
            create_token(tokens, input, {
                field: token.field,
                value: token.value,
                label: options.label_function ? options.label_function(token) : token.value
            });
        });

        if (parsed_query.tokens.length > 1) {
            // Create remove link
            var remove_all_facet_nuggets = getLocalizedNuggetsString('remove_all_facet_nuggets', 'Remove All');

            var $remove_all = $('<span class="nugget nugget_remove_all"><a href="javascript:void(0)" title="' + remove_all_facet_nuggets + '"><span class="label">' + remove_all_facet_nuggets + '</span></a></span>');
            $remove_all.appendTo(tokens);
            // Attach event handler to remove first nugget
            $remove_all.on("click", function (event) {
                tokens.find('.nugget').trigger('nugget.remove_all');
                input.closest('form').submit();
            });
        }

        if (parsed_query.tokens.length > 0) {
            $nugget_container.show();
        }


        // If the user clicks the back button we need to reset the text field
//    var reinitialize = function() {
//      parsed_query = parse_query(input.val(), options.fields || []);
//      input.val(parsed_query.query);
//      clearInterval(reinitialize_watcher);
//    };
//
//    $(window).unload(function(event) {
//      var reinitialize_watcher = setInterval(reinitialize, 1000);
//    });


        nugget_box.find(".nugget").bind("nugget.remove_all", function (event) {
            var $target = $(event.target);
            $target.remove();
        });

        nugget_box.find(".nugget").bind("nugget.remove", function (event) {
            var $target = $(event.target);
            $target.remove();
        });


        return this;
    }
})(jQuery);


