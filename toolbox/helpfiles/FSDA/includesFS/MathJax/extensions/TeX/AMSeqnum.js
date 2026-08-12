/*
 *  /MathJax/extensions/TeX/AMSeqnum.js
 *
 *  Adds automatic equation numbering and cross-referencing (\label, \ref,
 *  \eqref) to the AMSmath extension of this (MathJax 1.1) build.
 *
 *  It is a *pure add-on*: no file of the stock distribution is modified.
 *  It patches AMSmath.js in memory once "TeX AMSmath Ready" has been posted,
 *  so loading it is enough to enable the feature.
 *
 *  ---------------------------------------------------------------------
 *  USAGE
 *  ---------------------------------------------------------------------
 *
 *      <script type="text/x-mathjax-config">
 *        MathJax.Hub.Config({
 *          jax: ["input/TeX","output/HTML-CSS"],
 *          extensions: ["tex2jax.js","TeX/AMSeqnum.js"],
 *          TeX: { equationNumbers: { autoNumber: "AMS" } }
 *        });
 *      </script>
 *      <script src="includesFS/MathJax/MathJax.js"></script>
 *
 *  Note that AMSeqnum.js goes in the *page* extension list, not in
 *  TeX.extensions: it has to be in place before tex2jax scans the page, so
 *  that a bare \ref{...} / \eqref{...} in running text is recognised.  It
 *  pulls in AMSmath.js by itself.  (Listing it in TeX.extensions also works,
 *  but then \ref must be written inside math delimiters, e.g. $\ref{eq:x}$.)
 *
 *  In the page:
 *
 *      \begin{equation}\label{eq:mcd} d^2_i=(x_i-\mu)'\Sigma^{-1}(x_i-\mu)\end{equation}
 *      ... see equation \eqref{eq:mcd} above ...
 *
 *  ---------------------------------------------------------------------
 *  WHAT YOU GET
 *  ---------------------------------------------------------------------
 *
 *   * automatic numbering of equation, align, gather, multline, alignat
 *     (starred forms stay unnumbered), or of *every* display when
 *     autoNumber is set to "all";
 *   * \label{key}   -- attaches a key (and a DOM id) to the current number;
 *   * \ref{key}     -- prints the bare number, hyperlinked to the equation;
 *   * \eqref{key}   -- same, but formatted, i.e. (3);
 *   * \notag / \nonumber  -- suppress the number of the current line;
 *   * \tag{...} / \tag*{...} -- manual tag, also usable with \label;
 *   * forward references: a \ref to an equation that appears *later* on the
 *     page is resolved by a second, automatic rendering pass;
 *   * re-rendering (font change from the MathJax menu, Hub.Rerender, ...)
 *     is idempotent: numbers do not drift.
 *
 *  ---------------------------------------------------------------------
 *  CONFIGURATION (TeX.equationNumbers)
 *  ---------------------------------------------------------------------
 *
 *   autoNumber    "AMS"  number only the AMS numbered environments (default)
 *                 "all"  number every displayed equation
 *                 "none" number nothing (only explicit \tag)
 *   tagPlacement  "adjacent" (default) put the number right after the
 *                 formula, the two being centred together;
 *                 "margin"   classic LaTeX placement, i.e. flush against
 *                 the right border of the page (this is what MathJax does
 *                 out of the box, and it can look very far away on a wide
 *                 browser window)
 *   tagSpacing    gap between formula and number, default "1.5em"
 *   formatNumber(n)  number -> string shown            default String(n)
 *   formatTag(s)     string -> tag as displayed        default "(s)"
 *   formatID(s)      string -> DOM id                  default "mjx-eqn-s"
 *   formatURL(id)    id     -> href of a \ref          default "#id"
 *   useLabelIds      true: DOM ids are built from the \label, not from the
 *                    number (default true)
 *
 *  Runtime helpers:
 *
 *   MathJax.InputJax.TeX.resetEquationNumbers(n, keepLabels)
 *
 *  This file is licensed under the same Apache License 2.0 as MathJax.
 */

MathJax.Extension["TeX/AMSeqnum"] = { version: "1.0" };

/*  =====================================================================
 *  Part 1 -- runs as soon as the file is loaded.
 *  ===================================================================== */

/*  Make sure AMSmath.js is loaded: everything below patches it.  */
MathJax.Hub.Register.StartupHook("TeX Jax Require", function () {
  var cfg = MathJax.InputJax.TeX.config;
  cfg.extensions = cfg.extensions || [];
  for (var i = 0, m = cfg.extensions.length; i < m; i++) {
    if (cfg.extensions[i] === "AMSmath.js") { return; }
  }
  cfg.extensions.push("AMSmath.js");
});

/*  \ref / \eqref in running text.
 *
 *  tex2jax of this build only recognises math delimiters and
 *  \begin{...}...\end{...}.  Teach it to pick up a bare \ref{...} or
 *  \eqref{...} as well (the "processRefs" option of MathJax 2.x).
 *  Set tex2jax: {processRefs: false} to switch it off.
 *
 *  This has to happen before the page is pre-processed, which is why this
 *  file belongs in the page-level "extensions" list.
 */
(function () {

  var patchTex2jax = function () {
    var T2J = MathJax.Extension.tex2jax;
    if (!T2J || T2J.AMSeqnumPatched) { return; }
    T2J.AMSeqnumPatched = true;
    if (typeof T2J.config.processRefs === "undefined") { T2J.config.processRefs = true; }

    var oldCreatePatterns = T2J.createPatterns;
    T2J.createPatterns = function () {
      oldCreatePatterns.call(this);
      if (this.config.processRefs) {
        this.start = new RegExp(this.start.source + "|\\\\(?:eq)?ref\\s*\\{[^}]*\\}", "g");
      }
    };

    var oldStartMatch = T2J.startMatch;
    T2J.startMatch = function (match, element) {
      if (/^\\(eq)?ref\s*\{/.test(match[0])) {
        this.search = {
          mode: "", end: "", open: element, olen: 0,
          opos: this.pattern.lastIndex - match[0].length,
          close: element, cpos: this.pattern.lastIndex, clen: 0,
          matched: true
        };
        element = this.encloseMath(element);
        this.switchPattern(this.start);
        return element;
      }
      return oldStartMatch.call(this, match, element);
    };
  };

  if (MathJax.Extension.tex2jax) {
    patchTex2jax();
  } else {
    MathJax.Ajax.LoadHook("[MathJax]/extensions/tex2jax.js", patchTex2jax);
  }

})();

MathJax.Ajax.Styles({
  ".MathJax_ref":   { "text-decoration": "none" },
  "a .MathJax_ref": { "text-decoration": "none" }
});

/*  =====================================================================
 *  Part 2 -- runs once AMSmath.js has installed itself.
 *  ===================================================================== */

MathJax.Hub.Register.StartupHook("TeX AMSmath Ready", function () {

  var MML       = MathJax.ElementJax.mml,
      TEX       = MathJax.InputJax.TeX,
      TEXDEF    = TEX.Definitions,
      STACKITEM = TEX.Stack.Item,
      HUB       = MathJax.Hub;

  /*  ---------------- configuration ---------------------------------- */

  var CONFIG = HUB.Insert({
    autoNumber:   "AMS",
    useLabelIds:  true,
    tagPlacement: "adjacent",   // "adjacent" | "margin"
    tagSpacing:   "1em",
    formatNumber: function (n)  { return String(n); },
    formatTag:    function (n)  { return "(" + n + ")"; },
    formatID:     function (n)  { return "mjx-eqn-" + String(n).replace(/\s/g, "_"); },
    formatURL:    function (id) { return "#" + escape(id); }
  }, (TEX.config.equationNumbers || {}));
  TEX.config.equationNumbers = CONFIG;

  /*  ---------------- state ------------------------------------------ */

  var AMS = MathJax.Extension["TeX/AMSeqnum"];

  AMS.number   = 0;     // current equation number
  AMS.labels   = {};    // label -> {tag:"3", id:"mjx-eqn-eq:foo"} (whole page)
  AMS.eqlabels = {};    // ... defined by the equation being parsed
  AMS.refs     = [];    // scripts holding an unresolved \ref
  AMS.badref   = false; // current equation contains an unresolved \ref
  AMS.reparse  = false; // we are re-parsing an already numbered script
  AMS.updating = false; // the \ref resolution pass is running
  AMS.display  = false; // current equation is in display mode

  /*  ---------------- tag construction -------------------------------- */

  AMS.tagMtd = function (text) {
    return MML.mtd(MML.mtext(MML.chars(text)).With({ displaystyle: false }));
  };

  //  Create the automatic tag for the equation/line being closed.
  AMS.autoTag = function (global) {
    if (global.tag || global.nonumber || global.notag) { return; }
    AMS.number++;
    global.tagID = CONFIG.formatNumber(AMS.number);
    global.tag   = AMS.tagMtd(CONFIG.formatTag(global.tagID));
  };

  //  Give the tag its DOM id and record the \label that points at it.
  AMS.useTag = function (global) {
    var tag = global.tag;
    if (!tag) { return; }
    global.tagged = true;
    var key = (CONFIG.useLabelIds && global.label) ? global.label : global.tagID;
    if (key != null) { tag.id = CONFIG.formatID(key); }
    if (global.label && !AMS.reparse) {
      AMS.eqlabels[global.label] = { tag: global.tagID, id: tag.id || "" };
    }
  };

  AMS.clearTag = function (global) {
    delete global.tag;
    delete global.tagID;
    delete global.label;
    delete global.nonumber;
  };

  /*  ---------------- "adjacent" tag placement -------------------------- */

  /*  With the classic (LaTeX) placement the number lives in a separate
   *  <mlabeledtr> label column which the HTML-CSS output pushes against the
   *  right border of the *container*, i.e. of the whole page.  With
   *  tagPlacement:"adjacent" the number is instead glued to the end of the
   *  formula, and formula + number are centred together.
   */
  var ADJACENT = (CONFIG.tagPlacement === "adjacent");

  //  the tag mtd, unwrapped into an <mrow> that carries the anchor id
  AMS.tagNode = function (tagMtd) {
    var tag = (tagMtd.data.length === 1 ? MML.mrow(tagMtd.data[0])
                                        : MML.mrow.apply(MML, tagMtd.data));
    if (tagMtd.id) { tag.id = tagMtd.id; }
    return tag;
  };

  //  content  ->  <mrow> content <mspace/> (n) </mrow>
  AMS.joinTag = function (content, tagMtd) {
    return MML.mrow(content,
                    MML.mspace().With({ width: CONFIG.tagSpacing }),
                    AMS.tagNode(tagMtd));
  };

  //  Add the pending tag as one extra cell at the end of a table row, so
  //  that the numbers of the different lines line up in a column of their
  //  own just after the formulas.
  AMS.pushTagCell = function (row, global) {
    row.push(MML.mtd(MML.mrow(MML.mspace().With({ width: CONFIG.tagSpacing }),
                              AMS.tagNode(global.tag)))
                .With({ columnalign: MML.ALIGN.LEFT }));
    delete global.tag;
  };

  //  ... but a multline table is a single 100%-wide column, so there the
  //  tag has to be merged into the content of the last line instead.
  AMS.joinTagToRow = function (row, global) {
    var cell = row[row.length - 1];
    if (!cell) { return; }
    row[row.length - 1] = MML.mtd(AMS.joinTag(cell.data[0] || MML.mrow(), global.tag))
                             .With({ columnalign: cell.columnalign });
    delete global.tag;
  };

  /*  ---------------- new / redefined macros --------------------------- */

  HUB.Insert(TEXDEF, {
    macros: {
      tag:       "HandleTag",
      notag:     "HandleNoTag",
      nonumber:  "HandleNoTag",
      label:     "HandleLabel",
      ref:       ["HandleRef", false],
      eqref:     ["HandleRef", true]
    },
    environment: {
      equation:    ["EquationBegin", "EquationEnd", true],
      "equation*": ["EquationBegin", "EquationEnd", false]
    }
  });

  TEX.Parse.Augment({

    //  \begin{equation} / \begin{equation*}
    EquationBegin: function (begin, numbered) {
      begin.isEquation = true;
      this.stack.global.numbered = numbered;
      return begin;
    },
    EquationEnd: function (begin, data) {
      return data;
    },

    //  \tag{...}  \tag*{...}
    HandleTag: function (name) {
      var global = this.stack.global;
      var arg = this.trimSpaces(this.GetArgument(name)), star = (arg === "*");
      if (star) { arg = this.trimSpaces(this.GetArgument(name)); }
      if (global.notag) { TEX.Error(name + " not allowed in " + global.notag + " environment"); }
      if (global.tag)   { TEX.Error("Multiple " + name); }
      global.tagID = arg;
      global.tag   = MML.mtd.apply(MML, this.InternalMath(star ? arg : CONFIG.formatTag(arg)));
      delete global.nonumber;
    },

    //  \notag  \nonumber
    HandleNoTag: function (name) {
      var global = this.stack.global;
      delete global.tag;
      delete global.tagID;
      global.nonumber = true;
    },

    //  \label{key}
    HandleLabel: function (name) {
      var global = this.stack.global;
      var label  = this.trimSpaces(this.GetArgument(name));
      if (label === "") { return; }
      if (global.label) { TEX.Error("Multiple " + name + "'s"); }
      global.label = label;
      if (!AMS.reparse && (AMS.labels[label] || AMS.eqlabels[label])) {
        TEX.Error("Label '" + label + "' multiply defined");
      }
    },

    //  \ref{key}  \eqref{key}
    HandleRef: function (name, eqref) {
      var label = this.trimSpaces(this.GetArgument(name));
      var ref   = AMS.eqlabels[label] || AMS.labels[label];
      if (!ref) { ref = { tag: "???", id: "" }; AMS.badref = true; }
      var tag = (eqref ? CONFIG.formatTag(ref.tag) : ref.tag);
      var mml = MML.mrow.apply(MML, this.InternalMath(tag)).With({ "class": "MathJax_ref" });
      if (ref.id) { mml.href = CONFIG.formatURL(ref.id); }
      this.Push(mml);
    }

  });

  /*  ---------------- numbered environments ---------------------------- */

  //  multline / multline*  (pass the "numbered" flag on to the stack item)
  var oldMultline = TEX.Parse.prototype.Multline;
  TEX.Parse.Augment({
    Multline: function (begin, numbered) {
      var item = oldMultline.call(this, begin, numbered);
      item.numbered = numbered;
      return item;
    }
  });

  var oldMultlineEndTable = STACKITEM.multline.prototype.EndTable;
  STACKITEM.multline.Augment({
    EndTable: function () {
      if (this.numbered) { AMS.autoTag(this.global); }
      AMS.useTag(this.global);
      if (ADJACENT && this.global.tag) {
        //  flush the pending line, then glue the number to the last one
        if (this.data.length) { this.EndEntry(); }
        if (this.row.length)  { this.EndRow(); }
        if (this.table.length) { AMS.joinTagToRow(this.table[this.table.length - 1], this.global); }
      }
      oldMultlineEndTable.call(this);
      AMS.clearTag(this.global);
    }
  });

  //  align, align*, gather, gather*, alignat, alignat*, split, ...
  var oldAMSarrayEndRow = STACKITEM.AMSarray.prototype.EndRow;
  STACKITEM.AMSarray.Augment({
    EndRow: function () {
      if (this.numbered && !this.global.notag) { AMS.autoTag(this.global); }
      AMS.useTag(this.global);
      if (ADJACENT && this.global.tag) { AMS.pushTagCell(this.row, this.global); }
      oldAMSarrayEndRow.call(this);
      AMS.clearTag(this.global);
    }
  });

  //  equation, equation* and (when autoNumber is "all") any display
  var oldStartCheckItem = STACKITEM.start.prototype.checkItem;
  STACKITEM.start.Augment({
    checkItem: function (item) {
      if (item.type === "stop") {
        var global = this.global;
        if (!global.tagged &&
            (global.numbered || (CONFIG.autoNumber === "all" && AMS.display))) {
          AMS.autoTag(global);
        }
        AMS.useTag(global);
        if (ADJACENT && global.tag) {
          this.data = [AMS.joinTag(this.mmlData(), global.tag)];
          delete global.tag;
        }
        var result = oldStartCheckItem.call(this, item);
        AMS.clearTag(global);
        return result;
      }
      return oldStartCheckItem.call(this, item);
    }
  });

  if (CONFIG.autoNumber === "none") {
    AMS.autoTag = function () {};
  }

  /*  ---------------- per-equation bookkeeping -------------------------- */

  var oldPrefilter  = TEX.prefilterMath,
      oldPostfilter = TEX.postfilterMath;

  TEX.Augment({

    prefilterMath: function (math, displaymode, script) {
      AMS.display  = displaymode;
      AMS.eqlabels = {};
      AMS.badref   = false;
      AMS.reparse  = false;
      if (script) {
        //  Remember where the counter stood the first time this script was
        //  seen, and rewind to it on every later pass, so that numbers are
        //  stable across re-renderings and across the \ref resolution pass.
        if (script.MathJaxEqnStart != null) {
          AMS.saved   = AMS.number;          // restored in postfilterMath
          AMS.number  = script.MathJaxEqnStart;
          AMS.reparse = true;
        } else {
          script.MathJaxEqnStart = AMS.number;
        }
      }
      return oldPrefilter.call(TEX, math, displaymode, script);
    },

    postfilterMath: function (data, displaymode, script) {
      var result = oldPostfilter.call(TEX, data, displaymode, script);
      for (var label in AMS.eqlabels) {
        if (AMS.eqlabels.hasOwnProperty(label)) { AMS.labels[label] = AMS.eqlabels[label]; }
      }
      if (AMS.badref && !AMS.reparse && script) { AMS.refs.push(script); }
      if (AMS.reparse) { AMS.number = AMS.saved; AMS.reparse = false; }
      return result;
    },

    //  Restart numbering (e.g. per chapter).  Call it before typesetting.
    resetEquationNumbers: function (n, keepLabels) {
      AMS.number = (n || 0);
      AMS.refs   = [];
      if (!keepLabels) { AMS.labels = {}; }
      var scripts = document.getElementsByTagName("script");
      for (var i = 0, m = scripts.length; i < m; i++) {
        if (scripts[i].MathJaxEqnStart != null) { delete scripts[i].MathJaxEqnStart; }
      }
    }

  });

  /*  ---------------- forward references: second pass -------------------- */

  HUB.Register.MessageHook("End Math", function () {
    if (AMS.updating || !AMS.refs.length) { return; }
    var scripts = AMS.refs; AMS.refs = [];
    AMS.updating = true;
    var queue = [];
    for (var i = 0, m = scripts.length; i < m; i++) {
      queue.push(["Reprocess", HUB, scripts[i]]);
    }
    queue.push(function () { AMS.updating = false; });
    HUB.Queue.apply(HUB, queue);
  });

  MathJax.Hub.Startup.signal.Post("TeX AMSeqnum Ready");

});

MathJax.Ajax.loadComplete("[MathJax]/extensions/TeX/AMSeqnum.js");
