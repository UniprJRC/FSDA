//Copyright 2015, The MathWorks Inc.
(function ($) {
    $.getParameterByName = function(name) {
        name = name.replace(/[\[]/, "\\\[").replace(/[\]]/, "\\\]");
        var regexS = "[\\?&]" + name + "=([^&#]*)";
        var regex = new RegExp(regexS, 'g');
        var results;
        var value = "";
        while (true) {
            results = regex.exec(window.location.href);
            if (results == null) {
                break;
            }
            value = decodeURIComponent(results[1].replace(/\+/g, " "));
        }
        return value;
    };


    $(document).ready(function () {
        if ($.getParameterByName("browser") === "F1help") {
            //body#responsive_offcanvas
            $('body#responsive_offcanvas').addClass('navigation_off').addClass('no_animate');
            $('.site_container:first').removeClass('site_toc_closed').addClass('navigation_off');
            // g1028547: Detect the index-not-found page through the occurrence or the alert container,
            // and add another class to the site container to show the banner inside the search_crumb_container
            if ($('.form_alert_container').length > 0) {
                $('.site_container:first').addClass('show_banner');
            }

            if ($.getParameterByName("client") !== "addons") {
                $("#doc_center_content").on('mouseenter', 'a', function () {
                    if ($(this).attr('href').indexOf("matlab:") === 0) {
                        return;
                    }
                    if ($(this).hasClass('corrected_url')) {
                        return;
                    }
                    $(this).addClass('corrected_url');
                    $(this).attr('href', function (i, h) {
                        if (h === undefined) {
                            return "";
                        }
                        var srcUrl, hash, hashIndex;
                        if (h.indexOf('#') > 0) {
                            hashIndex = h.indexOf('#');
                            hash = h.substring(hashIndex, h.length);
                        } else {
                            hash = '';
                            hashIndex = h.length;
                        }

                        srcUrl = h.substring(0, hashIndex);
                        return srcUrl +
                            (srcUrl.indexOf('?') != -1 ? "&browser=F1help" : "?browser=F1help") +
                            hash;
                    });
                });
            }
        }
    });
})(window.jQuery);


