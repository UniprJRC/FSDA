function handleDocLinksClick3P(aTags) {
  var i;
  for (i = 0; i < aTags.length; i++) {
    aTags[i].onclick = function(evt) {
         var href = getHref(evt);
         if (href) {
             var hrefString = String(href);
             if (hrefString) {
                 var currentPageHost = window.location.host;
                 var currentPageProtocol = window.location.protocol;
                 var currentPage = window.location.href;
                 if (hrefString && hrefString.match(/^\s*matlab:.*/)) {
                     evt.stopImmediatePropagation();
                     var messageObj = {
                         "url" : hrefString,
                         "currenturl" : currentPage
                     };
                     var services = {
                        "messagechannel" : "matlab"
                     }
                     requestHelpService(messageObj, services, function() {});
                     return false;
                 } else if (hrefString
                     && (hrefString.indexOf(currentPageProtocol) < 0 || hrefString.indexOf(currentPageHost) < 0)
                     && hrefString.indexOf('http') === 0) {
                     evt.stopImmediatePropagation();
                     var messageObj = {
                         "url" : hrefString
                     };
                     var services = {
                        "messagechannel" : "externalclick"
                     }
                     requestHelpService(messageObj, services , function() {});
                     return false;
                 }
             }
         }
    }                        
  }
}

function getHref(evt) {
    if (evt.target && ($(evt.target).attr("href"))) {
        return evt.target;
    }

    if (evt.currentTarget && ($(evt.currentTarget).attr("href"))) {
        return evt.currentTarget;
    }

    return null;
}

function handleContextMenu3P(iframeElement) {
    iframeElement.contentDocument.oncontextmenu = function (e) {
        var id = "docinfo_" + new Date().getTime();

        var selectedText = getSelectedText3P(iframeElement.contentDocument);        
        var selectedUrl = getSelectedUrl(e.target);
        var pageUrl = window.location.href;

        var offsets = iframeElement.getBoundingClientRect();
        var clientX = e.clientX + offsets.left;
        var clientY = e.clientY + offsets.top;

        var messageObj = {
            "domchannel" : 'domeventinfo',
            "title" : window.document.title,
            "url": pageUrl,
            "domevent": 'contextmenu',
            "eventData" : {
                "pageX": e.pageX,
                "pageY": e.pageY,
                "clientX": clientX,
                "clientY": clientY,
                "selectedText": selectedText,
                "selectedUrl": selectedUrl
            },
            "id" : id
        };

        e.preventDefault();
        helpserviceIframePostMessage(messageObj);
    };
}

function getSelectedText3P(contentDocument) {
    var text = "";
    if (contentDocument.getSelection) {
        text = contentDocument.getSelection().toString();
    }    
    return text;
}