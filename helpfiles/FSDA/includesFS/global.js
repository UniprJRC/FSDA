// JavaScript Document

// Variable to be used by various macros so functions can only be defined when used on the page
// Not using MW to avoid conflict with the MWA variable
var MWS = new Array();

$(document).ready( function() {
	// Sticks the banner to the top:
	setTimeout(function() { 
		$(".sticky_header_container").affix({offset: { top: function() { 
			return (this.bottom = $('.header').length > 0 ? Math.max($("#header_desktop").outerHeight(true), $("#header_mobile").outerHeight(true)) : 0)
		}}})
	}, 100);

	/* Go to top button */
	$("<a class=\"btn btn_primary\" id=\"go-top\"><span class=\"icon-arrow-open-up icon_24\"></span></a>").on("click", function(event) {
		event.preventDefault();
		$('html, body').animate({scrollTop: 0}, 300);
	}).appendTo("body");

	$(window).on("scroll", function() {
		if ($(this).scrollTop() > 200) {
			$('#go-top').fadeIn(200);
		} else {
			$('#go-top').fadeOut(200);
		}
	});

	/* Smooth scroll on all anchor links on the same page */
	$('a.add_smooth_scroll[href*=#]').on("click", function(e) {
		e.preventDefault();
		$('html, body').animate({ scrollTop: $( $.attr(this, 'href') ).offset().top }, 500);
		location.hash = $(this).attr("href").split("#")[1];
	});
	
	/* Swipe controls for the carousel */
	if ($(".carousel").length) {
		var carousels = $(".carousel");
		/* Swipe/Touch Support for Bootstrap Carousels */
		carousels.swiperight(function() { $(this).carousel('prev'); });  
		carousels.swipeleft(function()  { $(this).carousel('next'); });

		carousels.find(".carousel-control, .carousel-indicators li").on("click", function() {
			$(this).closest(".carousel").carousel("pause");
		});
	}

	/* Tabs Support */
	if ($(".tab-container").length) {
		fakewaffle.responsiveTabs(['xs', 'sm']);
		$(".tab-container .panel-group.responsive").each(function() {
			$(this).find(".accordion-toggle:first").click();
			$(this).find(".accordion-toggle:gt(0)").addClass("collapsed");
		});
	}

	/* Shadow(box) onload */
	if($("#shadowonload").length) {
		$('#shadowonload').modal('toggle').on("click", function(event) {
			if (event.target == this) {
				window.location.href = $(this).find(".modal-header a.close").attr("href");
			}
		});
	}


	// Expand Collapse
	if($(".expand_collapse").length) {
		MWS["expand"] = {
			"querystr": "",
			"curpage": "",
			"open_obj": "",
			"getExpandData": function() {
				return sessionStorage.getItem('expand')
			},
			"storeExpandData": function() {
				//session storage, will expire when broswer is closed
				sessionStorage.expand = JSON.stringify(MWS["expand"].open_obj);
			},
			"getUrlVars": function() {
				var vars = [], hash;
				var hashes = window.location.href.slice(window.location.href.indexOf('?') + 1).split('&');
				for (var i = 0; i < hashes.length; i++) {
						hash = hashes[i].split('=');
						vars.push(hash[0]);
						vars[hash[0]] = hash[1];
				}
				return vars;
			}						
		}
		MWS["expand"].querystr = window.location.href.slice(window.location.href.indexOf('?') + 1);
		if (window.location.href.indexOf('?') < 0) {
			MWS["expand"].querystr = "";
		}

		MWS["expand"].curpage = { "url": window.location.pathname, "query": MWS["expand"].querystr, "open": [] };
		MWS["expand"].open_obj = { "page": [ MWS["expand"].curpage ] };
		
		if (MWS["expand"].getExpandData()) {
			var found = 0;
			MWS["expand"].open_obj = JSON.parse(MWS["expand"].getExpandData());
			//find the page and assign to curpage var
			$.each(MWS["expand"].open_obj.page, function(index, value) {
				if(value.url == window.location.pathname) {
					MWS["expand"].curpage = value;
					found = 1;
				}
			});
			//if the page wasn't found then add it
			if (!found) {
				MWS["expand"].open_obj.page.push(MWS["expand"].curpage);
				MWS["expand"].storeExpandData();
			}
		} else {
			MWS["expand"].storeExpandData();
		}

		// Get the query string value and store in an array
		var expandparam = MWS["expand"].getUrlVars()["expand"];
		if (typeof(expandparam) != "undefined") {
			expandparam = expandparam.split(",");
		} else { 
			expandparam = [];
		}


		// Determine which way to use as the on page load open. Clear the open values if necessary.
		var toggletype = "default"; // Everything is closed unless the html wants a default open
			
		if (MWS["expand"].querystr == "" && MWS["expand"].curpage.query != "") {
			toggletype = "default";
			MWS["expand"].curpage.open = [];
		} else if (MWS["expand"].querystr != "" && MWS["expand"].curpage.query != MWS["expand"].querystr && expandparam.length > 0) { // Query string if different than the last used query string
			toggletype = "query";
			MWS["expand"].curpage.open = [];
		} else if ( MWS["expand"].curpage.open.length > 0 ) { // Cookie values if there is a cookie already set
			toggletype = "cookie";
		} else if (MWS["expand"].querystr != "" && expandparam.length > 0) { // Query string if there isn't a cookie
			toggletype = "query";
			MWS["expand"].curpage.open = [];
		} else {
			MWS["expand"].curpage.open = [];
		}
		// store the query string in the cookie for comparing on next page load
		MWS["expand"].curpage.query = MWS["expand"].querystr;


		$("body").on("click", ".expand_trigger", function() {
			var clickedexpander = $(this);
			index = clickedexpander.data("index");
			if(clickedexpander.hasClass("collapsed")) { // collapsing
				if($.inArray(index, MWS["expand"].curpage.open) < 0) {
					MWS["expand"].curpage.open.push(index);
					MWS["expand"].storeExpandData();
				}
			} else { //expanding
				MWS["expand"].curpage.open.splice($.inArray(index, MWS["expand"].curpage.open), 1);
				MWS["expand"].storeExpandData();
			}
			
			//find the collapse / expand that relates to this and toggle if necessary
			setTimeout(function() {
				var toggle = $("."+clickedexpander.closest(".expand_collapse").data("toggle"));
				var target = toggle.data("target");
				var totalcount = 1;
				var totalactive = 1;
				if (target == "all") {
					totalcount = $(".expand_trigger").length;
					totalactive = $(".expand_trigger:not(.collapsed)").length;
				} else {
					totalcount = $("."+target).find(">.expand_trigger + .collapse").length;
					totalactive = $("."+target).find(">.expand_trigger + .collapse.in").length;
				}
				console.log("1 " + totalcount + " : " + totalactive);
				if (totalcount == totalactive) {
					toggle.find(".expand_all").hide();
					toggle.find(".collapse_all").show();
				} else if (totalactive == 0) {
					toggle.find(".expand_all").show();
					toggle.find(".collapse_all").hide();
				}
			}, 500);
		});
		
		
						
		$("body").on("click", ".expand_all", function(event) {
			event.preventDefault();
			var target = $(this).hide().next().show().parent().data("target");
			var targetIds = new Array();
			if (target == "all") {
				$(".expand_trigger.collapsed").each(function() {
					$(this).click();
					targetIds.push($(this).parent().data("index"));
				});
			} else {
				$("."+target+" > .expand_trigger.collapsed").each(function() {
					$(this).click();
					targetIds.push($(this).parent().data("index"));
				});
			}
			
			$.each(targetIds, function(index, value) {
				if($.inArray(value, MWS["expand"].curpage.open) < 0) {
					MWS["expand"].curpage.open.push(value);
				}
			});
			MWS["expand"].storeExpandData();
		});

		$("body").on("click", ".collapse_all", function(event) {
			event.preventDefault();
			var target = $(this).hide().prev().show().parent().data("target");
			var targetIds = new Array();
			if (target == "all") {
				$(".expand_trigger:not(.collapsed)").each(function() {
					$(this).click();
					targetIds.push($(this).parent().data("index"));
				});
			} else {
				$("."+target+" > .expand_trigger:not(.collapsed)").each(function() {
					$(this).click();
					targetIds.push($(this).parent().data("index"));
				});
			}

			$.each(targetIds, function(index, value) {
				MWS["expand"].curpage.open.splice($.inArray(value, MWS["expand"].curpage.open), 1);
			});
			MWS["expand"].storeExpandData();

		});
	
		$(".expand_trigger").each(function(index) {
			var expander = $(this);
			expander.data("index", index);
			var expanderIds = new Array();
			// If it should be opened on page load then open it
			if ( (toggletype == "query" && $.inArray(expander.data("target").substr(1), expandparam) > -1) || 
					 (toggletype == "cookie" && $.inArray(index, MWS["expand"].curpage.open) > -1) || 
					 (toggletype == "default" && expander.hasClass("default_open")) ) {
				expander.click();
				expanderIds.push(index);
			}
			$.each(expanderIds, function(index, value) {
				if($.inArray(value, MWS["expand"].curpage.open) < 0) {
					MWS["expand"].curpage.open.push(value);
				}
			});
			MWS["expand"].storeExpandData();
		});
	
		switches = $(".switch");
		switchcount = switches.length;
		switches.each(function(index) {
			var switcher = $(this);
			var totalcount = 0;
			var totalactive = 0;
			if (switchcount == 1) {
				switcher.data("target", "all").addClass("mwToggle0");
				$(".expand_collapse").data("toggle", "mwToggle0");
				totalcount = $(".expand_trigger").length;
				totalactive = $(".expand_trigger:not(.collapsed)").length;
			} else {
				if (switcher.data("target") == "" || switcher.data("target") == null) {
					switcher.data("target", "mwExpandToggle"+index).addClass("mwToggle"+index);
					//Change from parent to closest
					switcher.parent().addClass("mwExpandToggle"+index).data("toggle", "mwToggle"+index);
				} else {
					var target = switcher.addClass("mwToggle"+index).data("target");
					$("."+target).data("toggle", "mwToggle"+index);
				}
				totalcount = $("."+switcher.data("target")).find(">.expand_trigger + .collapse").length;
				totalactive = $("."+switcher.data("target")).find(">.expand_trigger + .collapse.in").length;
			}
			console.log("2 " + totalcount + " : " + totalactive);
			if (totalcount == totalactive) {
				switcher.find(".expand_all").hide();
				switcher.find(".collapse_all").show();
			} else {
				switcher.find(".expand_all").show();
				switcher.find(".collapse_all").hide();
			}
		});
	}
	// End Expand All
	
	
	
});