/*global window, $, document */

function usePlatformScreenshots() {
    "use strict";
    var ua = window.navigator.userAgent,
        platformPattern = /\/examples\/(\w+)\/(glnxa64|maci64|win64)\//,
        platform = 'win64';
    if (ua.match(/Linux/)) {
        platform = 'glnxa64';
    } else if (ua.match(/Macintosh/)) {
        platform = 'maci64';
    }
    $('img').each(function () {
        var img = $(this),
            src = img.attr('src');
        src = src.replace(platformPattern, '/examples/$1/' + platform + '/');
        img.attr('src', src);
        if (img.parent('.card_media').length === 1) {
            var cardMediaDiv = img.parent('.card_media').first();
            var backgroundImgSrc = cardMediaDiv.css('background-image');
            backgroundImgSrc = backgroundImgSrc.replace(platformPattern, '/examples/$1/' + platform + '/');
            cardMediaDiv.css('background-image', backgroundImgSrc);
        }
    });

}

$(document).ready(usePlatformScreenshots);

$(window).bind('examples_cards_added', function(e) {
   usePlatformScreenshots();
});
