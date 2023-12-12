var pageLang = document.getElementById("doc_center_content").getAttribute("lang");

// g1884073
if ( (window.parent && window.parent.location !== window.location) || (window.location.href.indexOf('http') === 0) ) 
    {
        // g1930999: China site is under mathworks.cn
        if (location.hostname.indexOf('mathworks.cn') > 0)
        {
            $('#doc_survey').attr('src','https://ww2.mathworks.cn/programs/bounce_hub_help.html?s_cid=Help_Topic_Survey&surveyParams='+document.location.href+'-G11N-'+pageLang);
        } else if ((location.hostname.indexOf('mathworks.com') > 0) || (location.hostname.indexOf('localhost') > -1)) 
        {
            // g1917116: this line takes care of both web-format and in-product docs.
            $('#doc_survey').attr('src','https://www.mathworks.com/programs/bounce_hub_help.html?s_cid=Help_Topic_Survey&surveyParams='+document.location.href+'-G11N-'+pageLang);
        }
     
    }
var ifrm = document.getElementById("doc_survey");
ifrm.style.display = "none";
window.addEventListener("message", function(evt){
  if(evt.origin.toLowerCase().indexOf('www.getfeedback.com') >=0)
  {
  ifrm.style.display = "";
  }
});


// Copyright 2018-2019 The MathWorks, Inc.
