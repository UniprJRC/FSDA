// Display only the products not included in the "not coming from product" list.
function handleComingFromProductList(prodList, listID, divID) {
    'use strict';
    $.getJSON('not_coming_from_product.json', function (notComingFromProductJson) {
      var filteredList = prodList.filter(function (obj) {
          if (! ('shortname' in obj)) { 
              return true;
          }

          for (var i = 0; i < notComingFromProductJson.length; i++) {
              if (('shortname' in notComingFromProductJson[i]) && (notComingFromProductJson[i].shortname === obj.shortname)) {
                  return false;
              }
          }

          return true;
      });
      
      handleProductList(filteredList, listID, divID);      
    });
}

function handleProductList(docSetItemList, listID, divID) {
  var webOnlyProducts = getWebOnlyProducts();
  var pageLocation = window.location.href;
  var filtereddocSetItemList = docSetItemList.filter(function (obj) {
      if ('shortname' in obj) {
            if ((typeof isRequestFromArchiveArea === 'function') && isRequestFromArchiveArea(pageLocation) && webOnlyProducts.indexOf(obj.shortname) >= 0) {
                return false;
            }
          return true;
      } else {
          return false;
      }
  });

  if (filtereddocSetItemList !== undefined && filtereddocSetItemList.length > 0) {
    var compiledTmpl = JST.installedHspTmpl({installedhsps: filtereddocSetItemList});
    $('#' + listID).append(compiledTmpl);
    $('#' + divID).show();
  }
}

function getWebOnlyProducts() {
  if (typeof BaseCodeMap === 'function') {
    var baseCodeMap = new BaseCodeMap();
    var webOnlyProducts = baseCodeMap.global;
    webOnlyProducts.push("matlabdrive");
    webOnlyProducts.push("matlabgrader");
    var index = webOnlyProducts.indexOf("install");
    if (index != -1) {
      webOnlyProducts.splice(index, 1);
    }
    return webOnlyProducts; 
  } else {
    return [];
  }
  
}

function setVisibility() {
    getProductDataFromCookie();
}

function setProductsVisibility(cookieData) {
    if (cookieData == null || cookieData.length == 0) {
        handleProducts(null);
        handleAddOns(null);
        return;
    }
    var baseCodes = cookieData[0];
    handleProducts(baseCodes);
    if (cookieData.length > 1) {
        var addOns = cookieData[1];
        handleAddOns(addOns);
    } else {
        handleAddOns(null);
    }
}

function getProductDataFromCookie() {
   var cookieRegexp = /MW_Doc_Template="?([^;"]*)/;
   var cookies = document.cookie;
   var matches = cookieRegexp.exec(cookies);
   if (matches != null && matches.length > 0) {
       var docCookie = matches[1];
       var parts = docCookie.split(/\|\|/);
       if (parts[0].indexOf("ONLINE") !== -1) {
           handleMATLABOnlineProductsList(null);
           return;
       }
       if (parts.length > 3) {
           setProductsVisibility([parts[1].split(/\W+/), parts[3].split(/\W+/)]);
       } else if (parts.length > 1) {
           setProductsVisibility([parts[1].split(/\W+/)]);
       } else {
           setProductsVisibility(null);
       }
   } else {
       setProductsVisibility(null);
   }
}

function handleMATLABOnlineProductsList() {
    $.ajax({
            url: '/help/matlab_online_products.json',
            method: 'get',
            dataType: 'json',
            success: function(data) {
                setProductsVisibility(data.basecodes);
            },
            error: function (xhr, ajaxOptions, thrownError) {
                setProductsVisibility(null);
                console.log("Server respond status: " + xhr.status + " : " + thrownError);
            }     
        });
}                


function handleProducts(baseCodes) {
    if (baseCodes == null || baseCodes.length == 0 || (baseCodes.length == 1 && baseCodes[0].length == 0)) {
        handleSelectedProducts(null);
    } else {
        var shortNames = new BaseCodeMap().getProductShortNames(baseCodes);
        handleSelectedProducts(shortNames);
    }
}

function handleAddOns(addOns) {
    if (addOns == null || addOns.length == 0 || (addOns.length == 1 && addOns[0].length == 0)) {
        handleSelectedAddOns(null);
    } else {
        handleSelectedAddOns(addOns);
    }
}