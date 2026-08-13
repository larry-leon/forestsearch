/* =========================================================================
   forestsearch — Extreme Subgroup Effects Under a Uniform Treatment Benefit
   Fixed-Baseline (Conditional-on-X) Simulation Study — slide deck builder
   Build:  node build_deck.js   →  forestsearch_fixed_baseline_slides.pptx
   ========================================================================= */
const pptxgen = require("pptxgenjs");

/* ── Ocean Gradient palette (forestsearch brand) ─────────────────────── */
const C = {
  navy:  "0C2340",  // dominant dark
  deep:  "1B4F72",  // deep blue
  blue2: "2471A3",
  teal:  "148F77",
  cream: "F7F1E1",
  gold:  "C9A227",
  orange:"E67E22",  // random-benchmark orange (matches plots)
  red:   "B22222",  // true-HR reference red
  ink:   "1A2530",
  gray:  "5F6B7A",
  line:  "D7DEE6",
  tintB: "EEF3F8",  // light blue tint
  tintT: "E9F5F1",  // light teal tint
  tintG: "FBF4DC",  // light gold tint
  tintO: "FDEFE2",  // light orange tint
  white: "FFFFFF"
};
const HEAD = "Cambria";
const BODY = "Calibri";
const MONO = "Courier New";

/* ── Optional icons via react-icons + sharp (graceful fallback) ───────── */
const fs = require("fs");
const path = require("path");
let ICON = {};
async function makeIcons() {
  try {
    const React = require("react");
    const { renderToStaticMarkup } = require("react-dom/server");
    const fa = require("react-icons/fa");
    const sharp = require("sharp");
    const wanted = {
      lock: fa.FaLock, dice: fa.FaDice, db: fa.FaDatabase,
      clock: fa.FaStopwatch, search: fa.FaSearch, check: fa.FaCheck,
      scale: fa.FaBalanceScale, sliders: fa.FaSlidersH,
      chart: fa.FaChartBar, cogs: fa.FaCogs, flask: fa.FaFlask,
      calc: fa.FaCalculator
    };
    fs.mkdirSync(path.join(__dirname, "icons"), { recursive: true });
    for (const [k, Cmp] of Object.entries(wanted)) {
      const svg = renderToStaticMarkup(
        React.createElement(Cmp, { color: "#FFFFFF", size: 256 })
      );
      const out = path.join(__dirname, "icons", `${k}.png`);
      await sharp(Buffer.from(svg)).resize(256, 256).png().toFile(out);
      ICON[k] = out;
    }
  } catch (e) {
    console.log("icons skipped:", e.message);
    ICON = {};
  }
}

/* ── Deck setup ───────────────────────────────────────────────────────── */
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";                 // 13.33 x 7.5
pres.author = "Larry F. Leon";
pres.title  = "Extreme Subgroup Effects Under a Uniform Treatment Benefit — Fixed-Baseline Simulation";

const W = 13.33, H = 7.5, MX = 0.62;
let pageNo = 0;

/* ── Small helpers (fresh option objects every call) ──────────────────── */
function footer(slide, dark = false) {
  pageNo += 1;
  slide.addText("forestsearch  \u00B7  fixed-baseline simulation study", {
    x: MX, y: 7.14, w: 6.5, h: 0.3, fontFace: BODY, fontSize: 8.5,
    color: dark ? "8FA3B8" : C.gray, align: "left", margin: 0
  });
  slide.addText(String(pageNo), {
    x: W - 1.1, y: 7.14, w: 0.5, h: 0.3, fontFace: BODY, fontSize: 8.5,
    color: dark ? "8FA3B8" : C.gray, align: "right", margin: 0
  });
}
function contentSlide(kicker, title, titleW) {
  const s = pres.addSlide();
  s.background = { color: C.white };
  s.addText(kicker.toUpperCase(), {
    x: MX, y: 0.30, w: 9.5, h: 0.26, fontFace: BODY, fontSize: 10.5,
    bold: true, color: C.teal, charSpacing: 2, margin: 0
  });
  s.addText(title, {
    x: MX, y: 0.52, w: titleW || (W - 2 * MX), h: 0.62, fontFace: HEAD,
    fontSize: 28, bold: true, color: C.navy, margin: 0
  });
  footer(s, false);
  return s;
}
function card(s, x, y, w, h, fill, lineCol) {
  s.addShape(pres.ShapeType.roundRect, {
    x, y, w, h, fill: { color: fill },
    line: lineCol ? { color: lineCol, width: 0.75 } : { color: fill, width: 0 },
    rectRadius: 0.07
  });
}
function bullets(s, items, x, y, w, h, opts) {
  const o = opts || {};
  const baseColor = o.color || C.ink;
  const gap = o.gap !== undefined ? o.gap : 8;
  const bul = o.noBullet ? false : { code: "2022", indent: 12 };
  const runs = [];
  items.forEach((t) => {
    if (t && Array.isArray(t.text)) {
      t.text.forEach((r, j) => {
        const ro = Object.assign({ color: baseColor }, r.options || {});
        ro.bullet = bul; ro.paraSpaceAfter = gap;
        if (j === t.text.length - 1) ro.breakLine = true;
        runs.push({ text: r.text, options: ro });
      });
    } else {
      const txt = t && t.text !== undefined ? t.text : t;
      const ro = Object.assign({ color: baseColor }, (t && t.options) || {});
      ro.bullet = bul; ro.breakLine = true; ro.paraSpaceAfter = gap;
      runs.push({ text: txt, options: ro });
    }
  });
  s.addText(runs, {
    x, y, w, h, fontFace: BODY, fontSize: o.size || 13, align: "left",
    valign: o.valign || "top", margin: 0, lineSpacingMultiple: 1.04
  });
}
function iconCircle(s, x, y, d, fill, key) {
  s.addShape(pres.ShapeType.ellipse, {
    x, y, w: d, h: d, fill: { color: fill }, line: { color: fill, width: 0 }
  });
  if (ICON[key]) {
    const pad = d * 0.26;
    s.addImage({ path: ICON[key], x: x + pad, y: y + pad, w: d - 2 * pad, h: d - 2 * pad });
  }
}
function statBig(s, x, y, w, num, label, color) {
  s.addText(num, { x, y, w, h: 0.78, fontFace: HEAD, fontSize: 40, bold: true,
    color: color || C.deep, align: "left", margin: 0 });
  s.addText(label, { x, y: y + 0.74, w, h: 0.55, fontFace: BODY, fontSize: 11.5,
    color: C.gray, align: "left", margin: 0, lineSpacingMultiple: 1.0 });
}
function chip(s, x, y, w, text) {
  s.addShape(pres.ShapeType.roundRect, {
    x, y, w, h: 0.62, fill: { color: C.white, transparency: 88 },
    line: { color: "3A5573", width: 1 }, rectRadius: 0.10
  });
  s.addText(text, { x, y, w, h: 0.62, fontFace: BODY, fontSize: 12.5, bold: true,
    color: C.cream, align: "center", valign: "middle", margin: 0 });
}

/* =======================================================================
   SLIDE 1 — Title (dark)
   ======================================================================= */
function s01() {
  const s = pres.addSlide();
  s.background = { color: C.navy };
  s.addShape(pres.ShapeType.ellipse, { x: 10.4, y: -2.3, w: 6.4, h: 6.4,
    fill: { color: "12395F", transparency: 35 }, line: { color: "12395F", width: 0 } });
  s.addShape(pres.ShapeType.ellipse, { x: -1.8, y: 5.4, w: 4.6, h: 4.6,
    fill: { color: "0E6655", transparency: 55 }, line: { color: "0E6655", width: 0 } });

  s.addText("FORESTSEARCH SIMULATION STUDY", {
    x: MX, y: 0.85, w: 9, h: 0.32, fontFace: BODY, fontSize: 12, bold: true,
    color: C.gold, charSpacing: 3, margin: 0 });
  s.addText("Extreme Subgroup Effects Under a Uniform Treatment Benefit", {
    x: MX, y: 1.30, w: 11.9, h: 1.75, fontFace: HEAD, fontSize: 40, bold: true,
    color: C.white, margin: 0, lineSpacingMultiple: 1.02 });
  s.addText("A Fixed-Baseline (Conditional-on-X) Simulation Study Based on Observed Trial Data", {
    x: MX, y: 3.10, w: 11.6, h: 0.55, fontFace: BODY, fontSize: 18,
    color: "BFD3E6", margin: 0 });

  const cw = 2.86, gap = 0.22; let cx = MX;
  ["10,000 simulated trials", "56 pre-defined subgroups",
   "True HR = 0.70 for every patient", "Fixed panel:  N = 686"].forEach(t => {
    chip(s, cx, 4.05, cw, t); cx += cw + gap;
  });

  s.addText([
    { text: "Larry F. Le\u00F3n", options: { bold: true, color: C.white, breakLine: true } },
    { text: "forestsearch R package  \u00B7  larry-leon.github.io/forestsearch", options: { color: "9FB6CC", breakLine: true } },
    { text: "August 2026  \u00B7  companion to the random-X vignette (extreme_subgroups.qmd)", options: { color: "9FB6CC" } }
  ], { x: MX, y: 5.55, w: 10.5, h: 1.1, fontFace: BODY, fontSize: 14,
       margin: 0, lineSpacingMultiple: 1.25 });
  footer(s, true);
}

/* =======================================================================
   SLIDE 2 — The Core Question
   ======================================================================= */
function s02() {
  const s = contentSlide("Motivation", "The Core Question");

  bullets(s, [
    "Regulatory forest plots routinely display dozens of subgroup rows — and the smallest rows are the noisiest.",
    "A reviewer seeing HR = 1.4 with UB(HR) = 3.5 in 35 patients cannot tell a genuine signal from sampling noise.",
    "We answer the question by construction: build a DGM from real trial data where the truth is known exactly."
  ], MX, 1.45, 6.5, 2.2, { size: 13.5, gap: 10 });

  card(s, MX, 3.85, 6.5, 2.55, C.navy);
  s.addText("This study asks, directly:", { x: MX + 0.3, y: 4.05, w: 5.9, h: 0.3,
    fontFace: BODY, fontSize: 12, bold: true, color: C.gold, margin: 0 });
  s.addText("When every patient shares the same true HR = 0.70, how alarming can a stratified Cox subgroup analysis look by chance \u2014 given exactly these 686 patients?",
    { x: MX + 0.3, y: 4.40, w: 5.9, h: 1.85, fontFace: HEAD, fontSize: 17.5, bold: true,
      color: C.white, margin: 0, lineSpacingMultiple: 1.1 });

  const rx = 7.55, rw = 5.15;
  card(s, rx, 1.45, rw, 4.95, C.tintT, C.line);
  iconCircle(s, rx + 0.3, 1.75, 0.62, C.teal, "lock");
  s.addText("What is different in this deck", { x: rx + 1.05, y: 1.80, w: rw - 1.3, h: 0.5,
    fontFace: HEAD, fontSize: 16, bold: true, color: C.navy, margin: 0 });
  bullets(s, [
    { text: [{ text: "Every replicate re-runs the trial on the ", options: {} },
             { text: "same covariate panel", options: { bold: true } },
             { text: " \u2014 simulate_from_dgm(baseline = \"fixed\").", options: {} }] },
    "Subgroup definitions, memberships, and sizes are identical in all 10,000 trials \u2014 the N column is exact, not an average.",
    "Between-trial spread isolates outcome-generation noise, conditional on the enrolled patients (conditional-on-X).",
    "Only treatment assignment, entry times, latent event times, and censoring are re-drawn per trial."
  ], rx + 0.3, 2.55, rw - 0.6, 3.7, { size: 12.5, gap: 9 });
}

/* =======================================================================
   SLIDE 3 — Two Ways to Simulate a Trial
   ======================================================================= */
function s03() {
  const s = contentSlide("Design choice", "Two Ways to Simulate a Trial from One DGM");

  const rows = [
    [hdr(""), hdr("Random baseline  (unconditional)"), hdr("Fixed baseline  (this study)")],
    [lab("Covariate rows per trial"),
     td("Resample n rows with replacement from the 5,000-row super-population (df_super)"),
     tdB("dgm$df_source \u2014 every GBSG patient exactly once, in source order; no row sampling")],
    [lab("Per-trial N"),
     td("Free parameter (n = 686 chosen to match the source trial)"),
     tdB("nrow(df_source) = 686, fixed by construction \u2014 n is deliberately not passed")],
    [lab("Subgroup sizes"),
     td("Vary trial-to-trial (binomial sampling noise); reported N is a mean"),
     tdB("Exact and identical in every replicate")],
    [lab("Question answered"),
     td("\u201CWhat if this trial had enrolled a different sample from the same population?\u201D"),
     tdB("\u201CGiven exactly these patients, how variable are the subgroup analyses across re-runs?\u201D")],
    [lab("Sources of between-trial spread"),
     td("Covariate sampling + treatment + entry + event times + censoring"),
     tdB("Treatment + entry + event times + censoring only")]
  ];
  function hdr(t) { return { text: t, options: { fill: { color: C.navy }, color: C.white,
    bold: true, fontSize: 12.5, fontFace: BODY, valign: "middle" } }; }
  function lab(t) { return { text: t, options: { fill: { color: C.tintB }, color: C.navy,
    bold: true, fontSize: 11.5, fontFace: BODY, valign: "middle" } }; }
  function td(t)  { return { text: t, options: { fill: { color: C.white }, color: C.ink,
    fontSize: 11.5, fontFace: BODY, valign: "middle" } }; }
  function tdB(t) { return { text: t, options: { fill: { color: C.tintT }, color: C.ink,
    fontSize: 11.5, fontFace: BODY, valign: "middle" } }; }

  s.addTable(rows, {
    x: MX, y: 1.42, w: W - 2 * MX, colW: [2.35, 4.87, 4.87],
    border: { type: "solid", color: C.line, pt: 0.75 },
    margin: 0.07, rowH: [0.42, 0.78, 0.66, 0.62, 0.86, 0.66]
  });

  s.addText([
    { text: "Same DGM machinery either way \u2014 ", options: { color: C.ink } },
    { text: "only the row-sampling step changes.", options: { bold: true, color: C.navy } },
    { text: "  Both designs are exposed through the single baseline argument of simulate_from_dgm().", options: { color: C.gray } }
  ], { x: MX, y: 5.75, w: W - 2 * MX, h: 0.45, fontFace: BODY, fontSize: 13, margin: 0 });
}

/* =======================================================================
   SLIDE 4 — Pipeline at a glance
   ======================================================================= */
function s04() {
  const s = contentSlide("Workflow", "Simulation Pipeline at a Glance");

  const steps = [
    ["Time-scale prep", "rfstime \u2192 months (\u00F7 30.4375). Verify exp(\u03BC) against the KM median event time before anything else."],
    ["generate_aft_dgm_flex()", "Weibull AFT outcome model + AIC-selected censoring model. Stores df_super (5,000) and df_source (686 \u2014 the fixed panel)."],
    ["calibrate_k_treat()", "uniroot rescales the fitted treatment coefficient so the AHR is exactly 0.70 \u2014 the uniform truth."],
    ["Censoring calibration", "check_censoring_dgm() diagnoses; calibrate_cens_adjust() matches the simulated censoring rate to the observed 57.6%."],
    ["simulate_from_dgm(baseline = \"fixed\")", "One call = one trial on the fixed 686-patient panel: randomize, draw entry, event, and censoring times."],
    ["Cox loop + summaries", "Stratified Cox in 56 subgroups \u00D7 10,000 trials \u2192 ECI, Pr(HR<0.5), Pr(HR>1.0), Pr(UB\u22652), Pr(UB\u22653)."]
  ];
  const bw = 3.86, bh = 2.28, gx = 0.26, gy = 0.42;
  const x0 = MX, y0 = 1.55;
  steps.forEach((st, i) => {
    const r = Math.floor(i / 3), c = i % 3;
    const x = x0 + c * (bw + gx), y = y0 + r * (bh + gy);
    card(s, x, y, bw, bh, r === 0 ? C.tintB : C.tintT, C.line);
    s.addShape(pres.ShapeType.ellipse, { x: x + 0.22, y: y + 0.22, w: 0.5, h: 0.5,
      fill: { color: C.gold }, line: { color: C.gold, width: 0 } });
    s.addText(String(i + 1), { x: x + 0.22, y: y + 0.22, w: 0.5, h: 0.5, fontFace: HEAD,
      fontSize: 17, bold: true, color: C.navy, align: "center", valign: "middle", margin: 0 });
    s.addText(st[0], { x: x + 0.86, y: y + 0.20, w: bw - 1.05, h: 0.62, fontFace: BODY,
      fontSize: 13.5, bold: true, color: C.navy, margin: 0, valign: "top", lineSpacingMultiple: 0.98 });
    s.addText(st[1], { x: x + 0.24, y: y + 0.92, w: bw - 0.48, h: bh - 1.1, fontFace: BODY,
      fontSize: 11, color: C.ink, margin: 0, valign: "top", lineSpacingMultiple: 1.03 });
    if (c < 2) {
      s.addShape(pres.ShapeType.line, {
        x: x + bw + 0.03, y: y + bh / 2, w: gx - 0.06, h: 0,
        line: { color: C.gray, width: 1.5, endArrowType: "triangle" } });
    }
  });
}

/* =======================================================================
   SLIDE 5 — Source data
   ======================================================================= */
function s05() {
  const s = contentSlide("Ingredients", "Source Data \u2014 GBSG Breast Cancer Trial");

  iconCircle(s, MX, 1.5, 0.62, C.deep, "db");
  s.addText("survival::gbsg  \u00B7  hormone therapy vs control", {
    x: MX + 0.78, y: 1.56, w: 6.2, h: 0.5, fontFace: BODY, fontSize: 13.5,
    bold: true, color: C.deep, margin: 0 });

  statBig(s, MX,        2.35, 2.6, "686",  "patients\n(the fixed panel)", C.navy);
  statBig(s, MX + 2.65, 2.35, 2.6, "299",  "events (44%)", C.navy);
  statBig(s, MX,        3.85, 2.6, "59.4", "months \u2014 KM median\nevent time", C.navy);
  statBig(s, MX + 2.65, 3.85, 2.6, "57.6%", "censoring rate\n(calibration target)", C.navy);

  const rx = 6.35, rw = 6.35;
  card(s, rx, 1.45, rw, 2.5, C.tintB, C.line);
  s.addText("Covariates carried into the DGM", { x: rx + 0.28, y: 1.62, w: rw - 0.56, h: 0.35,
    fontFace: HEAD, fontSize: 15, bold: true, color: C.navy, margin: 0 });
  bullets(s, [
    { text: [{ text: "Continuous:  ", options: { bold: true } },
             { text: "age, size (mm), nodes, pgr, er", options: { fontFace: MONO, fontSize: 12 } }] },
    { text: [{ text: "Factors:  ", options: { bold: true } },
             { text: "meno, grade", options: { fontFace: MONO, fontSize: 12 } },
             { text: "   (grade also stratifies the analysis Cox model)", options: { fontSize: 11.5, color: C.gray } }] },
    { text: [{ text: "Treatment:  ", options: { bold: true } },
             { text: "hormon", options: { fontFace: MONO, fontSize: 12 } },
             { text: "   \u00B7   Outcome: recurrence-free survival, in months", options: { fontSize: 11.5 } }] }
  ], rx + 0.28, 2.05, rw - 0.56, 1.8, { size: 12.5, gap: 7 });

  card(s, rx, 4.15, rw, 2.25, C.tintG, C.line);
  iconCircle(s, rx + 0.26, 4.4, 0.56, C.gold, "search");
  s.addText("Time-scale consistency \u2014 check it first", { x: rx + 0.96, y: 4.44, w: rw - 1.2,
    h: 0.35, fontFace: HEAD, fontSize: 14.5, bold: true, color: C.navy, margin: 0 });
  bullets(s, [
    "The most common simulation failure: DGM built in days, analysis_time given in months \u21D2 universal administrative censoring.",
    { text: [{ text: "Diagnostic: exp(\u03BC) \u2248 observed KM median.  Here: ", options: {} },
             { text: "55.5 vs 59.4 months \u2713", options: { bold: true, color: C.teal } }] }
  ], rx + 0.28, 4.95, rw - 0.56, 1.4, { size: 12, gap: 7 });
}

/* =======================================================================
   SLIDE 6 — Outcome model
   ======================================================================= */
function s06() {
  const s = contentSlide("DGM \u00B7 outcome model", "Weibull AFT with model = \"null\" \u2014 a Uniform Truth");

  card(s, MX, 1.45, 7.0, 1.35, C.navy);
  s.addText([
    { text: "log T\u1D62  =  \u03BC  +  x\u1D62\u2032\u03B2  +  \u03C4\u00B7\u03B5\u1D62 ,      \u03B5\u1D62 ~ minimum extreme value", options: { breakLine: true } },
    { text: "\u03C8\u2070(L)  =  log(0.70)   for every covariate profile L", options: { color: C.gold } }
  ], { x: MX + 0.3, y: 1.62, w: 6.4, h: 1.05, fontFace: HEAD, fontSize: 16.5, bold: true,
       color: C.white, margin: 0, lineSpacingMultiple: 1.25 });

  bullets(s, [
    { text: [{ text: "Linear predictor: ", options: { bold: true } },
             { text: "treatment + all baseline covariates; a treatment \u00D7 subgroup interaction exists only under model = \"alt\".", options: {} }] },
    { text: [{ text: "model = \"null\"", options: { fontFace: MONO, fontSize: 12.5, bold: true } },
             { text: " removes the interaction entirely \u21D2 one common treatment effect. Any subgroup finding downstream is a false positive by construction.", options: {} }] },
    { text: [{ text: "k_treat", options: { fontFace: MONO, fontSize: 12.5, bold: true } },
             { text: " rescales the fitted treatment coefficient before sampling \u2014 the calibration hook used on the next slide (k_inter plays the same role for interactions under \"alt\").", options: {} }] },
    { text: [{ text: "Latent times are quantile-transform draws", options: { bold: true } },
             { text: " from the fitted distribution (via log(rexp(n))) \u2014 nothing is re-fitted per replicate.", options: {} }] }
  ], MX, 3.05, 7.0, 3.3, { size: 13, gap: 10 });

  const rx = 8.0, rw = 4.7;
  card(s, rx, 1.45, rw, 4.95, C.tintB, C.line);
  s.addText("The constructor returns two frames", { x: rx + 0.28, y: 1.65, w: rw - 0.56, h: 0.6,
    fontFace: HEAD, fontSize: 15, bold: true, color: C.navy, margin: 0 });
  bullets(s, [
    { text: [{ text: "df_super", options: { fontFace: MONO, fontSize: 12.5, bold: true, color: C.deep } },
             { text: " \u2014 5,000 rows resampled from the source; the sampling frame for the random-X design.", options: {} }] },
    { text: [{ text: "df_source", options: { fontFace: MONO, fontSize: 12.5, bold: true, color: C.teal } },
             { text: " \u2014 the prepared source data, every patient exactly once, carrying per-row linear predictors for outcome and censoring under both arms. This is the fixed-baseline panel.", options: {} }] },
    "Build settings here:  n_super = 5000,  seed = 99,  standardize = FALSE.",
    { text: [{ text: "Fail fast: ", options: { bold: true } },
             { text: "the vignette hard-stops at build time if df_source is missing (older package version) \u2014 not 1,000 lines later inside the parallel loop.", options: {} }] }
  ], rx + 0.28, 2.30, rw - 0.56, 4.0, { size: 12, gap: 9 });
}

/* =======================================================================
   SLIDE 7 — k_treat calibration
   ======================================================================= */
function s07() {
  const s = contentSlide("DGM \u00B7 treatment effect", "Calibrating k_treat \u2014 the AHR Is the Right Target");

  bullets(s, [
    { text: [{ text: "With k_treat = 1 the truth is whatever the GBSG fit produced. ", options: {} },
             { text: "calibrate_k_treat()", options: { fontFace: MONO, fontSize: 12.5, bold: true } },
             { text: " finds the multiplier by uniroot (range [\u22125, 5], tol 1e\u22126) so a chosen HR summary hits 0.70 exactly.", options: {} }] },
    { text: [{ text: "Two possible targets: ", options: {} },
             { text: "AHR", options: { bold: true } },
             { text: " = exp(mean log HR of potential outcomes) \u2014 equals the patient-level HR under \"null\" \u2014 or the ", options: {} },
             { text: "marginal Cox HR", options: { bold: true } },
             { text: " from a no-covariate fit on stacked potential outcomes.", options: {} }] },
    { text: [{ text: "The Part IV analysis is a Cox model stratified by grade. ", options: {} },
             { text: "By Cox non-collapsibility", options: { bold: true } },
             { text: ", the marginal HR is attenuated toward 1; a stratified fit recovers the patient-level HR.", options: {} }] },
    { text: [{ text: "So we set use_ahr = TRUE. ", options: { bold: true, color: C.navy } },
             { text: "Calibrating the marginal HR to 0.70 instead would drag the simulation median below 0.70 by the same non-collapsibility shift.", options: {} }] }
  ], MX, 1.5, 6.9, 4.6, { size: 13, gap: 12 });

  const rx = 7.95, rw = 4.75;
  card(s, rx, 1.5, rw, 3.5, C.navy);
  s.addText("Calibration result", { x: rx + 0.3, y: 1.68, w: rw - 0.6, h: 0.35,
    fontFace: BODY, fontSize: 12, bold: true, color: C.gold, margin: 0 });
  const kv = [
    ["k_treat", "0.956"],
    ["AHR  (target 0.70)", "0.700  \u2713"],
    ["Marginal Cox HR (attenuated)", "0.766"],
    ["Implied median event time  exp(\u03BC)", "55.5 mo"]
  ];
  kv.forEach((r, i) => {
    s.addText(r[0], { x: rx + 0.3, y: 2.12 + i * 0.68, w: 3.15, h: 0.6, fontFace: BODY,
      fontSize: 12.5, color: "CFE0EE", margin: 0, valign: "middle" });
    s.addText(r[1], { x: rx + 3.35, y: 2.12 + i * 0.68, w: 1.15, h: 0.6, fontFace: HEAD,
      fontSize: 15, bold: true, color: C.white, margin: 0, valign: "middle", align: "right" });
  });
  card(s, rx, 5.2, rw, 1.2, C.tintT, C.line);
  s.addText([
    { text: "One uniroot pattern, three calibrators:  ", options: { bold: true, color: C.navy } },
    { text: "calibrate_k_treat() \u00B7 calibrate_k_inter() \u00B7 calibrate_cens_adjust()", options: { fontFace: MONO, fontSize: 11.5, color: C.ink } }
  ], { x: rx + 0.25, y: 5.32, w: rw - 0.5, h: 0.95, fontFace: BODY, fontSize: 12,
       margin: 0, valign: "middle", lineSpacingMultiple: 1.1 });
}

/* =======================================================================
   SLIDE 8 — Censoring model selection
   ======================================================================= */
function s08() {
  const s = contentSlide("DGM \u00B7 censoring model", "Censoring \u2014 AIC-Selected AFT, Plus Two Manual Modes");

  bullets(s, [
    { text: [{ text: "Fitted symmetrically: ", options: { bold: true } },
             { text: "reverse the event indicator so censored observations become the \u201Cevent\u201D, then fit an AFT model to censoring times.", options: {} }] },
    { text: [{ text: "select_censoring = TRUE", options: { fontFace: MONO, fontSize: 12.5, bold: true } },
             { text: " (default): four candidates compete on AIC. Intercept-only candidates guard against overfitting when censoring is unrelated to covariates.", options: {} }] },
    "Censoring covariates default to the outcome model\u2019s; override via continuous_vars_cens / factor_vars_cens."
  ], MX, 1.48, 6.75, 2.5, { size: 12.5, gap: 9 });

  const rows = [
    [h("Candidate"), h("Distribution"), h("Covariates")],
    [m("Weibull"),    t("Weibull"),    t("treat + censoring covariates")],
    [m("LogNormal"),  t("Log-normal"), t("treat + censoring covariates")],
    [m("Weibull0"),   t("Weibull"),    t("intercept only")],
    [m("LogNormal0"), t("Log-normal"), t("intercept only")]
  ];
  function h(x) { return { text: x, options: { fill: { color: C.navy }, color: C.white, bold: true, fontSize: 11.5, fontFace: BODY, valign: "middle" } }; }
  function m(x) { return { text: x, options: { fill: { color: x === "Weibull" ? C.tintT : C.white }, color: C.ink, fontFace: MONO, fontSize: 11, valign: "middle", bold: x === "Weibull" } }; }
  function t(x) { return { text: x, options: { fill: { color: C.white }, color: C.ink, fontSize: 11.5, fontFace: BODY, valign: "middle" } }; }
  s.addTable(rows, { x: MX, y: 4.05, w: 6.75, colW: [1.75, 1.75, 3.25],
    border: { type: "solid", color: C.line, pt: 0.75 }, margin: 0.06,
    rowH: [0.36, 0.4, 0.4, 0.4, 0.4] });
  s.addText([
    { text: "Selected here:  Weibull, covariate-adjusted   \u00B7   ", options: { bold: true, color: C.teal } },
    { text: "\u03BC_c = 4.031 (exp = 56.3 mo),  \u03C4_c = 0.425", options: { color: C.ink } }
  ], { x: MX, y: 6.25, w: 6.75, h: 0.35, fontFace: BODY, fontSize: 12.5, margin: 0 });

  const rx = 7.75, rw = 4.95;
  card(s, rx, 1.48, rw, 2.4, C.tintB, C.line);
  iconCircle(s, rx + 0.26, 1.72, 0.56, C.deep, "sliders");
  s.addText("Manual modes  (select_censoring = FALSE)", { x: rx + 0.96, y: 1.76, w: rw - 1.15,
    h: 0.55, fontFace: HEAD, fontSize: 14, bold: true, color: C.navy, margin: 0 });
  bullets(s, [
    { text: [{ text: "Analytical: ", options: { bold: true } },
             { text: "supply cens_params$mu and $tau directly \u2014 no model is fitted.", options: {} }] },
    { text: [{ text: "Force-fit: ", options: { bold: true } },
             { text: "leave cens_params empty; one AFT with your family (cens_type) and formula (cens_intercept_only).", options: {} }] }
  ], rx + 0.28, 2.42, rw - 0.56, 1.4, { size: 12, gap: 7 });

  card(s, rx, 4.05, rw, 2.55, C.tintG, C.line);
  s.addText("Why a good fit still is not enough", { x: rx + 0.28, y: 4.22, w: rw - 0.56, h: 0.35,
    fontFace: HEAD, fontSize: 14, bold: true, color: C.navy, margin: 0 });
  s.addText("The fitted censoring distribution interacts with the administrative cutoff (analysis_time), so the realised censoring rate in simulated trials drifts away from the observed 57.6%. Part II diagnoses the drift and calibrates it out with a single log-scale shift, cens_adjust.",
    { x: rx + 0.28, y: 4.62, w: rw - 0.56, h: 1.85, fontFace: BODY, fontSize: 12,
      color: C.ink, margin: 0, lineSpacingMultiple: 1.12 });
}

/* =======================================================================
   SLIDE 9 — Entry times & observation window
   ======================================================================= */
function s09() {
  const s = contentSlide("DGM \u00B7 follow-up structure", "Entry Times and the Observation Window");

  // Timeline geometry: months 0..92 mapped to x
  const x0 = 1.0, x1 = 12.3, m0 = 0, m1 = 92;
  const mx = (m) => x0 + (m - m0) / (m1 - m0) * (x1 - x0);
  const yAx = 4.35;

  // recruitment band
  s.addShape(pres.ShapeType.rect, { x: mx(0), y: 1.7, w: mx(24) - mx(0), h: yAx - 1.7,
    fill: { color: C.tintT }, line: { color: C.tintT, width: 0 } });
  s.addText("recruitment window \u2014 entry ~ Uniform(0, max_entry = 24)", {
    x: mx(0) + 0.05, y: 1.76, w: mx(24) - mx(0) + 2.6, h: 0.28, fontFace: BODY,
    fontSize: 10.5, color: C.teal, bold: true, align: "left", margin: 0 });

  // axis
  s.addShape(pres.ShapeType.line, { x: x0, y: yAx, w: x1 - x0, h: 0,
    line: { color: C.ink, width: 1.75, endArrowType: "triangle" } });
  [[0, "0"], [24, "24"], [84, "84 = analysis_time"]].forEach(([m, lab]) => {
    s.addShape(pres.ShapeType.line, { x: mx(m), y: yAx - 0.06, w: 0, h: 0.12,
      line: { color: C.ink, width: 1.5 } });
    s.addText(lab, { x: mx(m) - 0.9, y: yAx + 0.10, w: 1.8, h: 0.28, fontFace: BODY,
      fontSize: 10.5, color: C.ink, align: "center", margin: 0, bold: m === 84 });
  });
  s.addText("calendar time (months)", { x: x1 - 2.2, y: yAx + 0.38, w: 2.2, h: 0.26,
    fontFace: BODY, fontSize: 10, italic: true, color: C.gray, align: "right", margin: 0 });

  // administrative cutoff dashed line
  s.addShape(pres.ShapeType.line, { x: mx(84), y: 1.62, w: 0, h: yAx - 1.55,
    line: { color: C.red, width: 1.5, dashType: "dash" } });

  // patients
  function patient(y, entry, end, kind, label) {
    s.addShape(pres.ShapeType.line, { x: mx(entry), y, w: mx(end) - mx(entry), h: 0,
      line: { color: C.deep, width: 2.5 } });
    s.addShape(pres.ShapeType.line, { x: mx(entry), y: y - 0.07, w: 0, h: 0.14,
      line: { color: C.deep, width: 2.5 } });
    if (kind === "event") {
      s.addShape(pres.ShapeType.ellipse, { x: mx(end) - 0.055, y: y - 0.055, w: 0.11, h: 0.11,
        fill: { color: C.red }, line: { color: C.red, width: 0 } });
    } else if (kind === "cens") {
      s.addShape(pres.ShapeType.ellipse, { x: mx(end) - 0.055, y: y - 0.055, w: 0.11, h: 0.11,
        fill: { color: C.white }, line: { color: C.deep, width: 1.75 } });
    } else {
      s.addShape(pres.ShapeType.rect, { x: mx(end) - 0.05, y: y - 0.07, w: 0.1, h: 0.14,
        fill: { color: C.orange }, line: { color: C.orange, width: 0 } });
    }
    s.addText(label, { x: mx(entry) + 0.06, y: y - 0.34, w: 6.4, h: 0.26, fontFace: BODY,
      fontSize: 10, color: C.gray, margin: 0 });
  }
  patient(2.72, 3, 54,  "event", "enters month 3 \u00B7 event observed  (\u25CF)");
  patient(3.30, 20, 71, "cens",  "enters month 20 \u00B7 latent censoring time reached  (\u25CB)");
  patient(3.88, 10, 84, "admin", "enters month 10 \u00B7 follow-up truncated at analysis_time  (administrative censoring)");

  bullets(s, [
    { text: [{ text: "entry_var = NULL", options: { fontFace: MONO, fontSize: 12, bold: true } },
             { text: " (default) draws synthetic uniform entry times; alternatively pass a column of observed entry times.", options: {} }] },
    { text: [{ text: "Observed time = min(latent event, latent censoring, analysis_time \u2212 entry)", options: { bold: true } },
             { text: " \u2014 so each patient\u2019s follow-up window is analysis_time \u2212 entry.", options: {} }] },
    "The admin cutoff interacts with the fitted censoring model \u21D2 the realised censoring rate drifts \u2014 the reason cens_adjust calibration exists (next slide).",
    "time_eos offers a further end-of-study truncation; rand_ratio / strata_rand control the randomisation itself."
  ], MX, 5.0, W - 2 * MX, 1.95, { size: 12.5, gap: 7 });
}

/* =======================================================================
   SLIDE 10 — Censoring diagnostics & calibration
   ======================================================================= */
function s10() {
  const s = contentSlide("DGM \u00B7 calibration", "Censoring Diagnostics and the cens_adjust Shift");

  const lw = 5.9;
  card(s, MX, 1.45, lw, 2.45, C.tintB, C.line);
  s.addText("1 \u00B7 Diagnose \u2014 check_censoring_dgm()", { x: MX + 0.28, y: 1.6, w: lw - 0.56,
    h: 0.35, fontFace: HEAD, fontSize: 14.5, bold: true, color: C.navy, margin: 0 });
  bullets(s, [
    "Three checks vs the DGM reference: overall censoring rate (tol 10 pp), censoring-time quantiles, KM median censoring time (tol 25%).",
    { text: [{ text: "Pilot trial at cens_adjust = 0:  59.2% vs 57.6%  (+1.6 pp) \u2014 passes, ", options: {} },
             { text: "but we calibrate anyway.", options: { bold: true } }] }
  ], MX + 0.28, 2.02, lw - 0.56, 1.8, { size: 12, gap: 7 });

  card(s, MX, 4.1, lw, 2.5, C.tintT, C.line);
  s.addText("2 \u00B7 Calibrate \u2014 calibrate_cens_adjust()", { x: MX + 0.28, y: 4.25, w: lw - 0.56,
    h: 0.35, fontFace: HEAD, fontSize: 14.5, bold: true, color: C.navy, margin: 0 });
  bullets(s, [
    "uniroot on a log-scale shift of the censoring distribution; target = \"rate\" (or \"km_median\").",
    "Objective: simulate \u2192 extract metric \u2192 return sim \u2212 ref. A fixed inner seed keeps the root-finder off Monte Carlo noise.",
    "Monotone: larger shift \u21D2 longer censoring times \u21D2 lower rate. If uniroot fails, widen interval = c(\u22123, 3)."
  ], MX + 0.28, 4.67, lw - 0.56, 1.85, { size: 12, gap: 6 });

  const rx = 6.85, rw = 5.85;
  card(s, rx, 1.45, rw, 2.9, C.navy);
  s.addText("Calibration result (target = \"rate\")", { x: rx + 0.3, y: 1.62, w: rw - 0.6, h: 0.32,
    fontFace: BODY, fontSize: 12, bold: true, color: C.gold, margin: 0 });
  const kv = [
    ["cens_adjust", "0.0412"],
    ["uniroot iterations  \u00B7  residual", "18  \u00B7  0.0002"],
    ["Reference \u2192 achieved censoring rate", "57.6% \u2192 57.6%"],
    ["Re-check at tighter tolerances (5 pp / 15%)", "+0.3 pp  \u2014 OK"]
  ];
  kv.forEach((r, i) => {
    s.addText(r[0], { x: rx + 0.3, y: 2.0 + i * 0.56, w: 3.65, h: 0.5, fontFace: BODY,
      fontSize: 12, color: "CFE0EE", margin: 0, valign: "middle" });
    s.addText(r[1], { x: rx + 3.9, y: 2.0 + i * 0.56, w: 1.75, h: 0.5, fontFace: HEAD,
      fontSize: 13.5, bold: true, color: C.white, margin: 0, valign: "middle", align: "right" });
  });

  card(s, rx, 4.55, rw, 2.05, C.tintG, C.line);
  iconCircle(s, rx + 0.26, 4.8, 0.56, C.gold, "lock");
  s.addText("Fixed-baseline detail", { x: rx + 0.96, y: 4.84, w: rw - 1.2, h: 0.35,
    fontFace: HEAD, fontSize: 14, bold: true, color: C.navy, margin: 0 });
  s.addText("The calibration objective forwards baseline = \"fixed\" through \u201C\u2026\u201D into simulate_from_dgm() \u2014 so n and n_eval must both equal the panel size (686). They are pinned, not free Monte Carlo sizes as in the random-X vignette.",
    { x: rx + 0.28, y: 5.3, w: rw - 0.56, h: 1.2, fontFace: BODY, fontSize: 12,
      color: C.ink, margin: 0, lineSpacingMultiple: 1.1 });
}

/* =======================================================================
   SLIDE 11 — baseline="fixed" mechanics  (FOCUS)
   ======================================================================= */
function s11() {
  const s = contentSlide("Focus \u00B7 fixed covariates", "baseline = \"fixed\" \u2014 What One Trial Draw Actually Does");

  const steps = [
    ["Take the panel as-is", "dgm$df_source: every source patient exactly once, in source order. No row sampling of any kind."],
    ["Randomise treatment", "treat_sim ~ 1:1 (rand_ratio; strata_rand available for stratified randomisation)."],
    ["Draw entry times", "entry ~ Uniform(0, max_entry = 24), or your observed entry column via entry_var."],
    ["Generate latent event times", "Quantile-transform draw from the fitted Weibull AFT under the assigned arm \u2014 per-row linear predictors were pre-computed at build time."],
    ["Censor and truncate", "Calibrated censoring model (+ cens_adjust = 0.0412), then the administrative cutoff at analysis_time = 84."]
  ];
  let y = 1.5;
  steps.forEach((st, i) => {
    s.addShape(pres.ShapeType.ellipse, { x: MX, y: y + 0.02, w: 0.44, h: 0.44,
      fill: { color: C.teal }, line: { color: C.teal, width: 0 } });
    s.addText(String(i + 1), { x: MX, y: y + 0.02, w: 0.44, h: 0.44, fontFace: HEAD,
      fontSize: 14, bold: true, color: C.white, align: "center", valign: "middle", margin: 0 });
    s.addText([
      { text: st[0] + "   ", options: { bold: true, color: C.navy } },
      { text: st[1], options: { color: C.ink } }
    ], { x: MX + 0.62, y, w: 6.35, h: 0.92, fontFace: BODY, fontSize: 12,
         margin: 0, valign: "top", lineSpacingMultiple: 1.05 });
    y += 0.97;
  });
  s.addText([
    { text: "Nothing is re-fitted per replicate", options: { bold: true, color: C.teal } },
    { text: " \u2014 each call is a draw from distributions stored in the DGM object.", options: { color: C.ink } }
  ], { x: MX, y: 6.42, w: 7.0, h: 0.4, fontFace: BODY, fontSize: 12.5, margin: 0 });

  const rx = 7.85, rw = 4.85;
  card(s, rx, 1.5, rw, 3.1, C.navy);
  s.addText("One replicate  (inside run_one_sim)", { x: rx + 0.28, y: 1.64, w: rw - 0.56,
    h: 0.3, fontFace: BODY, fontSize: 11.5, bold: true, color: C.gold, margin: 0 });
  s.addText(
    "df_s <- simulate_from_dgm(\n" +
    "  dgm           = dgm_uniform,\n" +
    "  baseline      = \"fixed\",  # n implied\n" +
    "  analysis_time = 84,\n" +
    "  max_entry     = 24,\n" +
    "  cens_adjust   = 0.0412,\n" +
    "  seed          = seed_base + ss\n" +
    ")",
    { x: rx + 0.28, y: 2.0, w: rw - 0.56, h: 2.5, fontFace: MONO, fontSize: 12,
      color: "E8F0F7", margin: 0, lineSpacingMultiple: 1.08 });

  card(s, rx, 4.85, rw, 1.75, C.tintB, C.line);
  s.addText("Other knobs (defaults shown in the vignette)", { x: rx + 0.28, y: 5.0, w: rw - 0.56,
    h: 0.3, fontFace: HEAD, fontSize: 13, bold: true, color: C.navy, margin: 0 });
  s.addText("rand_ratio \u00B7 draw_treatment \u00B7 entry_var \u00B7 strata_rand \u00B7 time_eos \u00B7 keep_rand \u00B7 hrz_crit",
    { x: rx + 0.28, y: 5.38, w: rw - 0.56, h: 1.1, fontFace: MONO, fontSize: 11.5,
      color: C.ink, margin: 0, lineSpacingMultiple: 1.25 });
}

/* =======================================================================
   SLIDE 12 — Fixed vs re-drawn  (FOCUS 2)
   ======================================================================= */
function s12() {
  const s = contentSlide("Focus \u00B7 fixed covariates", "What Is Frozen \u2014 and What Is Re-Drawn \u2014 Across 10,000 Trials");

  const cw = 5.9, ch = 3.9, y0 = 1.5;
  // FIXED card
  card(s, MX, y0, cw, ch, C.tintB, C.line);
  iconCircle(s, MX + 0.28, y0 + 0.24, 0.62, C.navy, "lock");
  s.addText("FIXED in every replicate", { x: MX + 1.04, y: y0 + 0.32, w: cw - 1.3, h: 0.45,
    fontFace: HEAD, fontSize: 16.5, bold: true, color: C.navy, margin: 0 });
  bullets(s, [
    "The covariate matrix: the 686 GBSG rows, in source order.",
    { text: [{ text: "Every covariate-defined subgroup \u2014 definition, membership, and size. ", options: {} },
             { text: "The N column is exact, not a mean.", options: { bold: true } }] },
    "Per-trial N itself: not a parameter \u2014 n is deliberately absent from sim_config.",
    { text: [{ text: "Structural zeros stay zero: ", options: { bold: true } },
             { text: "the panel has no node-negative patients, so 5 node-defined subgroups are empty in every trial (51 of 56 estimable).", options: {} }] }
  ], MX + 0.3, y0 + 1.0, cw - 0.6, ch - 1.15, { size: 12, gap: 8 });

  // RE-DRAWN card
  const rx = MX + cw + 0.3;
  card(s, rx, y0, cw, ch, C.tintO, C.line);
  iconCircle(s, rx + 0.28, y0 + 0.24, 0.62, C.orange, "dice");
  s.addText("RE-DRAWN each replicate", { x: rx + 1.04, y: y0 + 0.32, w: cw - 1.3, h: 0.45,
    fontFace: HEAD, fontSize: 16.5, bold: true, color: "B35A0F", margin: 0 });
  bullets(s, [
    "Treatment assignment (1:1 randomisation).",
    "Entry times ~ Uniform(0, 24 months).",
    "Latent event times from the Weibull AFT; latent censoring times from the calibrated censoring model; administrative cutoff at 84 months.",
    { text: [{ text: "random* membership", options: { bold: true } },
             { text: " \u2014 re-drawn over the fixed panel, so sizes are exactly 60 / 40 / 20 / 15: pure membership noise at fixed N.", options: {} }] }
  ], rx + 0.3, y0 + 1.0, cw - 0.6, ch - 1.15, { size: 12, gap: 8 });

  card(s, MX, 5.65, W - 2 * MX, 1.0, C.navy);
  s.addText([
    { text: "Interpretation shift:  ", options: { bold: true, color: C.gold } },
    { text: "\u201Cwhat if the trial had enrolled a different sample?\u201D  \u2192  \u201Cgiven exactly these patients, re-run the trial.\u201D  ", options: { color: C.white } },
    { text: "All between-trial spread is outcome-generation noise, conditional on X.", options: { bold: true, color: C.white } }
  ], { x: MX + 0.3, y: 5.65, w: W - 2 * MX - 0.6, h: 1.0, fontFace: BODY, fontSize: 13.5,
       margin: 0, valign: "middle", lineSpacingMultiple: 1.1 });
}

/* =======================================================================
   SLIDE 13 — Implementation notes
   ======================================================================= */
function s13() {
  const s = contentSlide("Focus \u00B7 fixed covariates", "Implementation Notes \u2014 for Analysts Reusing the Workflow");

  const items = [
    ["cogs", C.deep, "One sim_config, everywhere",
     "baseline, analysis_time = 84, max_entry = 24, and the seed scheme live in a single list reused by calibration, the single-trial demo, and the main loop \u2014 no duplication, no drift."],
    ["check", C.teal, "Fail fast on df_source",
     "The vignette hard-stops at DGM build time if df_source is missing (older forestsearch), rather than deep inside the parallel loop."],
    ["dice", C.orange, "Reproducible under any schedule",
     "Trial seed = seed_base + ss; random-benchmark seed adds a 1e6 offset. Iteration-local seeds \u21D2 bit-for-bit results regardless of worker execution order."],
    ["lock", C.navy, "Calibrate under the same design",
     "Pass baseline = \"fixed\" through calibrate_cens_adjust()\u2019s \u201C\u2026\u201D and pin n = n_eval = nrow(df_source) = 686."],
    ["sliders", C.gold, "Swap the analysis, not the loop",
     "The per-subgroup fit is a cox_fn argument (here: treat_sim + strata(grade)). Substitute an unstratified or adjusted model for sensitivity analyses."],
    ["db", C.blue2, "Bring your own trial",
     "Replace gbsg and the subgroup list; the pipeline \u2014 DGM, calibration, fixed-panel loop, summaries \u2014 is unchanged."]
  ];
  const cw2 = 5.9, ch2 = 1.62, gy = 0.16;
  items.forEach((it, i) => {
    const col = i % 2, row = Math.floor(i / 2);
    const x = MX + col * (cw2 + 0.3), y = 1.5 + row * (ch2 + gy);
    card(s, x, y, cw2, ch2, col === 0 ? C.tintB : C.tintT, C.line);
    iconCircle(s, x + 0.2, y + 0.2, 0.5, it[1], it[0]);
    s.addText(it[2], { x: x + 0.84, y: y + 0.16, w: cw2 - 1.0, h: 0.34, fontFace: BODY,
      fontSize: 13, bold: true, color: C.navy, margin: 0 });
    s.addText(it[3], { x: x + 0.84, y: y + 0.5, w: cw2 - 1.05, h: ch2 - 0.62, fontFace: BODY,
      fontSize: 10.6, color: C.ink, margin: 0, lineSpacingMultiple: 1.04, valign: "top" });
  });
}

/* =======================================================================
   SLIDE 14 — Execution & timing
   ======================================================================= */
function s14() {
  const s = contentSlide("Execution", "10,000 Trials \u00D7 56 Subgroups in 46 Seconds");

  const rows = [
    [h2("Stage"), h2("Elapsed"), h2("Per trial")],
    [t2("DGM construction  (generate_aft_dgm_flex)"), t2r("1.17 s"), t2r("\u2014")],
    [t2("Censoring calibration  (calibrate_cens_adjust)"), t2r("0.43 s"), t2r("\u2014")],
    [tb("Simulation loop  (10,000 \u00D7 56 subgroups, 116 workers)"), tbr("44.3 s"), tbr("\u22484 ms")],
    [t2b("Total"), t2br("45.9 s"), t2br("\u2014")]
  ];
  function h2(x)  { return { text: x, options: { fill: { color: C.navy }, color: C.white, bold: true, fontSize: 12, fontFace: BODY, valign: "middle" } }; }
  function t2(x)  { return { text: x, options: { fill: { color: C.white }, color: C.ink, fontSize: 12, fontFace: BODY, valign: "middle" } }; }
  function t2r(x) { return { text: x, options: { fill: { color: C.white }, color: C.ink, fontSize: 12, fontFace: BODY, valign: "middle", align: "right" } }; }
  function tb(x)  { return { text: x, options: { fill: { color: C.tintT }, color: C.ink, fontSize: 12, fontFace: BODY, valign: "middle", bold: true } }; }
  function tbr(x) { return { text: x, options: { fill: { color: C.tintT }, color: C.ink, fontSize: 12, fontFace: BODY, valign: "middle", align: "right", bold: true } }; }
  function t2b(x) { return { text: x, options: { fill: { color: C.tintB }, color: C.navy, fontSize: 12, fontFace: BODY, valign: "middle", bold: true } }; }
  function t2br(x){ return { text: x, options: { fill: { color: C.tintB }, color: C.navy, fontSize: 12, fontFace: BODY, valign: "middle", align: "right", bold: true } }; }

  s.addTable(rows, { x: MX, y: 1.55, w: 6.7, colW: [4.4, 1.2, 1.1],
    border: { type: "solid", color: C.line, pt: 0.75 }, margin: 0.07,
    rowH: [0.4, 0.5, 0.5, 0.62, 0.46] });

  iconCircle(s, MX, 4.6, 0.62, C.teal, "clock");
  statBig(s, MX + 0.85, 4.35, 3.0, "225.8", "trials / second  (wall-clock)", C.teal);
  statBig(s, MX + 3.9,  4.35, 3.0, "\u22484 ms", "per trial \u00B7 56 Cox fits each", C.deep);

  const rx = 7.6, rw = 5.1;
  card(s, rx, 1.55, rw, 4.85, C.tintB, C.line);
  s.addText("Parallel design", { x: rx + 0.28, y: 1.72, w: rw - 0.56, h: 0.35,
    fontFace: HEAD, fontSize: 15, bold: true, color: C.navy, margin: 0 });
  bullets(s, [
    "doFuture + multisession: 116 of 128 cores (91%; \u226590% target), portable across Windows / macOS / Linux.",
    "Embarrassingly parallel \u2014 iteration-local seeds, so results are identical under any execution order.",
    { text: [{ text: "Convergence: ", options: { bold: true } },
             { text: "\u2265 99.7% everywhere. random15: 9,971 / 10,000; all other subgroups \u2248 100%. Subgroups with < 5 rows are skipped by guard.", options: {} }] },
    { text: [{ text: "Structural note: ", options: { bold: true } },
             { text: "5 node-negative subgroups are empty in every replicate (all-node-positive panel) \u2192 51 estimable subgroups.", options: {} }] },
    "Projected budgets at this throughput: 1,000 sims \u2248 6 s; 5,000 sims \u2248 24 s."
  ], rx + 0.28, 2.15, rw - 0.56, 4.15, { size: 12, gap: 9 });
}

/* =======================================================================
   SLIDE 15 — Chart: Pr(UB>=2) vs subgroup size
   ======================================================================= */
function s15() {
  const s = contentSlide("Results", "Sample Size \u2014 Not Clinical Meaning \u2014 Drives Chance Extremes");

  const cats = [
    "random15  (N=15)", "random20  (N=20)", "Grade 3 / PGR high  (N=34)",
    "Pre-meno / Age>50  (N=35)", "random40  (N=40)", "Grade 3 / High nodes  (N=46)",
    "Grade 3 / ER-high  (N=55)", "Post-meno / G3 / ER-low  (N=58)",
    "random60  (N=60)", "Grade 3  (N=161)", "All Patients  (N=686)"
  ];
  const clin = [null, null, 55.2, 62.3, null, 22.2, 27.9, 17.3, null, 0.3, 0.0];
  const rand = [90.1, 80.4, null, null, 50.6, null, null, null, 31.1, null, null];

  s.addChart(pres.ChartType.bar, [
    { name: "Clinical / interaction subgroups", labels: cats, values: clin },
    { name: "Random benchmarks",                labels: cats, values: rand }
  ], {
    x: MX, y: 1.4, w: 8.35, h: 5.25,
    barDir: "bar", barGrouping: "clustered", barOverlapPct: 100, barGapWidthPct: 55,
    chartColors: [C.deep, C.orange],
    catAxisOrientation: "maxMin",
    showTitle: true, title: "Pr( UB(HR) \u2265 2.0 )  across 10,000 trials \u2014 true HR = 0.70 everywhere",
    titleFontSize: 13, titleColor: C.navy, titleFontFace: BODY,
    showValue: true, dataLabelPosition: "outEnd", dataLabelFormatCode: "0.#\"%\"",
    dataLabelFontSize: 10, dataLabelColor: C.ink, dataLabelFontFace: BODY,
    valAxisMinVal: 0, valAxisMaxVal: 100, valAxisMajorUnit: 20,
    valAxisLabelFormatCode: "0\"%\"", valAxisLabelFontSize: 10, valAxisLabelColor: C.gray,
    catAxisLabelFontSize: 10.5, catAxisLabelColor: C.ink, catAxisLabelFontFace: BODY,
    valGridLine: { color: "E4EAF0", size: 0.75 }, catGridLine: { style: "none" },
    showLegend: true, legendPos: "b", legendFontSize: 11, legendColor: C.ink
  });

  const rx = 9.25, rw = 3.45;
  card(s, rx, 1.55, rw, 2.3, C.tintT, C.line);
  s.addText("The Cox model is fine", { x: rx + 0.24, y: 1.7, w: rw - 0.48, h: 0.32,
    fontFace: HEAD, fontSize: 13.5, bold: true, color: C.navy, margin: 0 });
  s.addText("Median HR = 0.69\u20130.73 in every estimable subgroup \u2014 unbiased throughout. What changes with N is the spread, not the centre.",
    { x: rx + 0.24, y: 2.06, w: rw - 0.48, h: 1.7, fontFace: BODY, fontSize: 11.5,
      color: C.ink, margin: 0, lineSpacingMultiple: 1.12 });

  card(s, rx, 4.05, rw, 2.45, C.tintO, C.line);
  s.addText("Orange = noise floor", { x: rx + 0.24, y: 4.2, w: rw - 0.48, h: 0.32,
    fontFace: HEAD, fontSize: 13.5, bold: true, color: "B35A0F", margin: 0 });
  s.addText("The random benchmarks have no clinical meaning, yet at each size they bracket the clinical subgroups\u2019 tail risk. Reading a forest plot without this calibration invites over-interpretation.",
    { x: rx + 0.24, y: 4.56, w: rw - 0.48, h: 1.85, fontFace: BODY, fontSize: 11.5,
      color: C.ink, margin: 0, lineSpacingMultiple: 1.12 });
}

/* =======================================================================
   SLIDE 16 — Results table
   ======================================================================= */
function s16() {
  const s = contentSlide("Results", "Six Summary Statistics \u2014 Selected Subgroups");

  const data = [
    ["All Patients (ITT)",   "686", "0.1%",  "0.3%",  "0.72", "0%",    "0%",    "0.90", false],
    ["Grade 3",              "161", "5.1%",  "6.4%",  "0.72", "0.3%",  "0%",    "1.11", false],
    ["random60",             "60",  "20.7%", "21.4%", "0.71", "31.1%", "7.8%",  "1.63", true],
    ["Grade 3 / ER-high",    "55",  "18.3%", "20.5%", "0.73", "27.9%", "5.9%",  "1.59", false],
    ["random40",             "40",  "25.7%", "25.9%", "0.72", "50.6%", "22.8%", "2.01", true],
    ["Pre-meno / Age>50",    "35",  "29.9%", "28.7%", "0.71", "62.3%", "36.0%", "2.41", false],
    ["Grade 3 / PGR high",   "34",  "27.3%", "27.4%", "0.72", "55.2%", "27.1%", "2.14", false],
    ["random20",             "20",  "34.1%", "33.7%", "0.71", "80.4%", "60.5%", "3.67", true],
    ["random15",             "15",  "37.6%", "35.7%", "0.69", "90.1%", "76.7%", "5.83", true]
  ];
  const header = ["Subgroup", "N (exact)", "Pr(HR<0.5)", "Pr(HR>1.0)", "mHR",
                  "Pr(UB\u22652)", "Pr(UB\u22653)", "mUB"].map(t => ({
    text: t, options: { fill: { color: C.navy }, color: C.white, bold: true,
      fontSize: 11.5, fontFace: BODY, valign: "middle", align: t === "Subgroup" ? "left" : "center" } }));
  const rows = [header].concat(data.map(r => {
    const fill = r[8] ? C.tintO : C.white;
    return r.slice(0, 8).map((v, j) => ({
      text: v, options: { fill: { color: fill }, color: r[8] && j === 0 ? "B35A0F" : C.ink,
        bold: j === 0, fontSize: 11.5, fontFace: j === 0 ? BODY : BODY, valign: "middle",
        align: j === 0 ? "left" : "center" } }));
  }));
  s.addTable(rows, { x: MX, y: 1.45, w: W - 2 * MX, colW: [2.9, 1.15, 1.35, 1.35, 1.0, 1.35, 1.35, 1.64],
    border: { type: "solid", color: C.line, pt: 0.75 }, margin: 0.06,
    rowH: 0.40 });

  card(s, MX, 5.92, W - 2 * MX, 0.9, C.navy);
  s.addText([
    { text: "Read one row:  ", options: { bold: true, color: C.gold } },
    { text: "in a 15-patient subgroup, 9 of 10 trials produce UB \u2265 2 and 3 of 4 produce UB \u2265 3 \u2014 under a uniform benefit of HR = 0.70.  ", options: { color: C.white } },
    { text: "mHR \u2248 0.70 throughout: the estimator is unbiased; the tails are pure sampling noise.", options: { color: "BFD3E6" } }
  ], { x: MX + 0.3, y: 5.92, w: W - 2 * MX - 0.6, h: 0.9, fontFace: BODY, fontSize: 12.5,
       margin: 0, valign: "middle", lineSpacingMultiple: 1.08 });
  s.addText("Orange rows: random benchmarks (membership re-drawn per trial over the fixed panel; sizes exact by construction).",
    { x: MX, y: 5.52, w: W - 2 * MX, h: 0.3, fontFace: BODY, fontSize: 10.5, italic: true,
      color: C.gray, margin: 0 });
}

/* =======================================================================
   SLIDE 17 — Clinical rationale is no protection
   ======================================================================= */
function s17() {
  const s = contentSlide("Results", "Clinical Rationale Provides No Statistical Protection");

  const pairs = [
    ["Pre-meno / Age > 50", "N = 35", "62.3%", "random40  (N = 40)", "50.6%", "random20  (N = 20)", "80.4%"],
    ["Grade 3 / PGR high",  "N = 34", "55.2%", "random40  (N = 40)", "50.6%", "random20  (N = 20)", "80.4%"],
    ["Grade 3 / ER-high",   "N = 55", "27.9%", "random60  (N = 60)", "31.1%", "random40  (N = 40)", "50.6%"]
  ];
  const cw = 3.93, gy = 0.24; let x = MX;
  pairs.forEach(p => {
    card(s, x, 1.5, cw, 3.35, C.tintB, C.line);
    s.addText(p[0], { x: x + 0.24, y: 1.66, w: cw - 0.48, h: 0.36, fontFace: HEAD,
      fontSize: 14.5, bold: true, color: C.navy, margin: 0 });
    s.addText(p[1] + "  \u00B7  clinically defined", { x: x + 0.24, y: 2.04, w: cw - 0.48, h: 0.28,
      fontFace: BODY, fontSize: 10.5, color: C.gray, margin: 0 });
    s.addText(p[2], { x: x + 0.24, y: 2.34, w: cw - 0.48, h: 0.66, fontFace: HEAD,
      fontSize: 30, bold: true, color: C.deep, margin: 0 });
    s.addText("of trials reach UB(HR) \u2265 2", { x: x + 0.24, y: 2.98, w: cw - 0.48, h: 0.26,
      fontFace: BODY, fontSize: 10.5, color: C.gray, margin: 0 });
    s.addShape(pres.ShapeType.line, { x: x + 0.24, y: 3.38, w: cw - 0.48, h: 0,
      line: { color: C.line, width: 1 } });
    s.addText([
      { text: "size-matched noise floor\n", options: { fontSize: 10, color: C.gray } },
      { text: p[3] + ":  ", options: { fontSize: 11.5, color: C.ink } },
      { text: p[4] + "\n", options: { fontSize: 12.5, bold: true, color: C.orange } },
      { text: p[5] + ":  ", options: { fontSize: 11.5, color: C.ink } },
      { text: p[6], options: { fontSize: 12.5, bold: true, color: C.orange } }
    ], { x: x + 0.24, y: 3.5, w: cw - 0.48, h: 1.25, fontFace: BODY, margin: 0,
         lineSpacingMultiple: 1.15 });
    x += cw + gy;
  });

  bullets(s, [
    { text: [{ text: "At matched N, clinical and random subgroups share the same chance distribution", options: { bold: true, color: C.navy } },
             { text: " \u2014 a clinical label does not change the sampling behaviour of the estimator.", options: {} }] },
    { text: [{ text: "High-risk screen:", options: { bold: true } },
             { text: " a Pr(UB\u22652) \u2265 10% filter catches 11 of the 51 estimable subgroups \u2014 7 clinical / interaction subgroups plus all 4 random benchmarks \u2014 every one with N \u2264 87.", options: {} }] }
  ], MX, 5.2, W - 2 * MX, 1.6, { size: 13, gap: 9 });
}

/* =======================================================================
   SLIDE 18 — Does the baseline design matter?
   ======================================================================= */
function s18() {
  const s = contentSlide("Design robustness", "Fixed vs Random Baseline \u2014 Does the Choice Matter?");

  iconCircle(s, MX, 1.55, 0.62, C.deep, "scale");
  s.addText("Companion comparison: same DGM, same 56 subgroups, random-X vignette vs this fixed-X study",
    { x: MX + 0.8, y: 1.6, w: 10.8, h: 0.5, fontFace: BODY, fontSize: 13.5, bold: true,
      color: C.deep, margin: 0 });

  bullets(s, [
    { text: [{ text: "Inferential conclusions are essentially identical: ", options: { bold: true } },
             { text: "across the 56 subgroups, mean \u0394Pr(UB\u22652) \u2248 \u22120.3 percentage points between designs.", options: {} }] },
    { text: [{ text: "Residual differences trace to exact-N vs mean-N composition: ", options: { bold: true } },
             { text: "under random X, an \u201CN \u2248 35\u201D subgroup is sometimes 28 and sometimes 42 patients; under fixed X it is 35 in every replicate.", options: {} }] },
    "The extreme-subgroup phenomenon is a property of subgroup size and outcome noise \u2014 not of how the covariates were sampled."
  ], MX, 2.35, 11.9, 2.3, { size: 13, gap: 10 });

  const cw = 5.9;
  card(s, MX, 4.75, cw, 1.85, C.tintT, C.line);
  s.addText("Choose fixed when\u2026", { x: MX + 0.26, y: 4.9, w: cw - 0.52, h: 0.32,
    fontFace: HEAD, fontSize: 13.5, bold: true, color: C.teal, margin: 0 });
  s.addText("the audience is reviewing this trial: same patients every replicate, exact subgroup Ns in every output, and the clean conditional message \u2014 \u201Cthese patients, re-run.\u201D",
    { x: MX + 0.26, y: 5.26, w: cw - 0.52, h: 1.25, fontFace: BODY, fontSize: 12,
      color: C.ink, margin: 0, lineSpacingMultiple: 1.12 });
  card(s, MX + cw + 0.3, 4.75, cw, 1.85, C.tintB, C.line);
  s.addText("Choose random when\u2026", { x: MX + cw + 0.56, y: 4.9, w: cw - 0.52, h: 0.32,
    fontFace: HEAD, fontSize: 13.5, bold: true, color: C.deep, margin: 0 });
  s.addText("the question is about the enrolment process itself: population-sampling uncertainty is added on top, answering the unconditional \u201Ca different sample from the same population\u201D question.",
    { x: MX + cw + 0.56, y: 5.26, w: cw - 0.52, h: 1.25, fontFace: BODY, fontSize: 12,
      color: C.ink, margin: 0, lineSpacingMultiple: 1.12 });
}

/* =======================================================================
   SLIDE 19 — Conclusions (dark)
   ======================================================================= */
function s19() {
  const s = pres.addSlide();
  s.background = { color: C.navy };
  s.addShape(pres.ShapeType.ellipse, { x: 11.0, y: 4.9, w: 5.0, h: 5.0,
    fill: { color: "12395F", transparency: 40 }, line: { color: "12395F", width: 0 } });
  s.addText("CONCLUSIONS", { x: MX, y: 0.42, w: 6, h: 0.3, fontFace: BODY, fontSize: 11,
    bold: true, color: C.gold, charSpacing: 3, margin: 0 });
  s.addText("Five Takeaways", { x: MX, y: 0.72, w: 9, h: 0.7, fontFace: HEAD, fontSize: 30,
    bold: true, color: C.white, margin: 0 });

  const items = [
    ["Small subgroups produce extreme results by chance.",
     "At N \u2248 15, Pr(UB \u2265 2) = 90% and the median upper bound is 5.8 \u2014 under a uniform true HR of 0.70."],
    ["Clinical rationale is not statistical protection.",
     "Size-matched random benchmarks reproduce the clinical subgroups\u2019 HR and UB(HR) behaviour."],
    ["The Cox model is working correctly.",
     "Median HR \u2248 0.70 in every estimable subgroup; wide intervals in small subgroups are honest uncertainty."],
    ["The fixed-baseline design isolates outcome noise.",
     "Same 686 patients every replicate; subgroup Ns exact; conclusions match the random-X design (mean \u0394Pr(UB\u22652) \u2248 \u22120.3 pp)."],
    ["Simulation-based calibration is essential \u2014 and cheap.",
     "A full 10,000-trial, 56-subgroup benchmark runs in under a minute with forestsearch."]
  ];
  let y = 1.7;
  items.forEach((it, i) => {
    s.addShape(pres.ShapeType.ellipse, { x: MX, y: y + 0.03, w: 0.46, h: 0.46,
      fill: { color: C.gold }, line: { color: C.gold, width: 0 } });
    s.addText(String(i + 1), { x: MX, y: y + 0.03, w: 0.46, h: 0.46, fontFace: HEAD,
      fontSize: 15, bold: true, color: C.navy, align: "center", valign: "middle", margin: 0 });
    s.addText([
      { text: it[0] + "  ", options: { bold: true, color: C.white } },
      { text: it[1], options: { color: "BFD3E6" } }
    ], { x: MX + 0.68, y, w: 11.6, h: 0.9, fontFace: BODY, fontSize: 13.5, margin: 0,
         valign: "top", lineSpacingMultiple: 1.05 });
    y += 0.95;
  });

  s.addText([
    { text: "pak::pak(\"larry-leon/forestsearch\")", options: { fontFace: MONO, fontSize: 13, color: C.gold } },
    { text: "      \u00B7      larry-leon.github.io/forestsearch", options: { fontFace: BODY, fontSize: 13, color: "9FB6CC" } }
  ], { x: MX, y: 6.6, w: 11.5, h: 0.4, margin: 0 });
  footer(s, true);
}

/* =======================================================================
   SLIDE 20 — Appendix: HTE variant
   ======================================================================= */
function s20() {
  const s = contentSlide("Appendix", "The Same DGM with Genuine Heterogeneity \u2014 model = \"alt\"");

  bullets(s, [
    { text: [{ text: "Three extra knobs turn the null DGM into an HTE DGM: ", options: { bold: true } },
             { text: "subgroup_vars = c(\"er\", \"meno\"),  subgroup_cuts = list(er = 20, meno = 0),  k_inter = 1.5.", options: { fontFace: MONO, fontSize: 12 } }] },
    "This embeds a harm subgroup \u2014 ER-low and pre-menopausal \u2014 covering 18.8% of the panel; flag_harm marks membership in every simulated dataset.",
    { text: [{ text: "Truth by construction: ", options: { bold: true } },
             { text: "marginal HR \u2248 1.21 inside the harm subgroup vs a clear benefit in the complement \u2014 spline_spec and set_beta_spec allow still finer control of the effect surface.", options: {} }] },
    "The entire fixed-baseline pipeline is reused unchanged: shared sim_config, its own calibrated cens_adjust = 0.041, same panel of 686 patients.",
    { text: [{ text: "Why it matters: ", options: { bold: true, color: C.navy } },
             { text: "a known-truth HTE target for evaluating subgroup-discovery methods such as forestsearch() \u2014 the flip side of the false-positive benchmark in this deck.", options: {} }] }
  ], MX, 1.55, 7.3, 4.9, { size: 13, gap: 12 });

  const rx = 8.35, rw = 4.35;
  card(s, rx, 1.55, rw, 3.2, C.navy);
  s.addText("HTE build (delta vs null)", { x: rx + 0.26, y: 1.7, w: rw - 0.52, h: 0.3,
    fontFace: BODY, fontSize: 11.5, bold: true, color: C.gold, margin: 0 });
  s.addText(
    "generate_aft_dgm_flex(\n" +
    "  \u2026,\n" +
    "  model         = \"alt\",\n" +
    "  subgroup_vars = c(\"er\",\"meno\"),\n" +
    "  subgroup_cuts = list(er = 20,\n" +
    "                       meno = 0),\n" +
    "  k_inter       = 1.5, \u2026)",
    { x: rx + 0.26, y: 2.05, w: rw - 0.52, h: 2.6, fontFace: MONO, fontSize: 11.5,
      color: "E8F0F7", margin: 0, lineSpacingMultiple: 1.1 });

  card(s, rx, 4.95, rw, 1.5, C.tintT, C.line);
  s.addText([
    { text: "Harm-subgroup prevalence: ", options: { bold: true, color: C.navy } },
    { text: "18.8% of the panel (134 patients in the demo trial).", options: { color: C.ink } }
  ], { x: rx + 0.26, y: 5.05, w: rw - 0.52, h: 1.3, fontFace: BODY, fontSize: 12.5,
       margin: 0, valign: "middle", lineSpacingMultiple: 1.15 });
}

/* ── Build ────────────────────────────────────────────────────────────── */
(async () => {
  await makeIcons();
  [s01, s02, s03, s04, s05, s06, s07, s08, s09, s10,
   s11, s12, s13, s14, s15, s16, s17, s18, s19, s20].forEach(f => f());
  await pres.writeFile({ fileName: "forestsearch_fixed_baseline_slides.pptx" });
  console.log("Wrote forestsearch_fixed_baseline_slides.pptx with", pageNo, "slides");
})();

