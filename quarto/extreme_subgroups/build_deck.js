/**
 * build_deck.js
 *
 * Section 6.6.4 summary deck for the extreme_subgroups vignette.
 *
 * Renders 9 slides to extreme_subgroups_section_6_6_4.pptx.
 *
 * Topics:
 *   1. Title
 *   2. Statistical / clinical objective
 *   3. Simulation setup (data, DGM, calibration)
 *   4. Subgroup taxonomy (56 total, bar chart)
 *   5. Section 6.6.4 intro: high-risk filter
 *   6. UB(HR) distribution figure
 *   7. HR distribution figure
 *   8. Key takeaways (four-card grid)
 *   9. Resources
 */

const pptxgen = require("pptxgenjs");
const path = require("path");

const SW = 13.333;
const SH = 7.5;
const MARGIN = 0.5;

// Ocean Gradient palette (consistent with prior project decks)
const C = {
  navy:    "21295C",
  deep:    "065A82",
  teal:    "1C7293",
  cream:   "F8F5F0",
  white:   "FFFFFF",
  gold:    "C9A227",
  text:    "2C2C2C",
  muted:   "6B7280",
  random:  "E67E22"
};

const FONT = {
  header: "Georgia",
  body:   "Calibri"
};

const ASSETS = path.join(__dirname, "assets");

let pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.author = "Larry F. Leon";
pres.title  = "Extreme Subgroups — Section 6.6.4 Summary";

// ── Helpers ────────────────────────────────────────────────────────────────
const makeShadow = () => ({
  type: "outer", color: "000000", blur: 8, offset: 2, angle: 90, opacity: 0.10
});

function addStandardHeader(s, title, subtitle) {
  s.addText(title, {
    x: MARGIN, y: 0.35, w: SW - 2 * MARGIN, h: 0.65,
    fontSize: 28, fontFace: FONT.header, color: C.navy, bold: true,
    align: "left", valign: "top", margin: 0
  });
  if (subtitle) {
    s.addText(subtitle, {
      x: MARGIN, y: 1.0, w: SW - 2 * MARGIN, h: 0.4,
      fontSize: 14, fontFace: FONT.body, color: C.muted, italic: true,
      align: "left", valign: "top", margin: 0
    });
  }
}

// =============================================================================
// SLIDE 1 — Title
// =============================================================================
function buildSlide1() {
  let s = pres.addSlide();
  s.background = { color: C.navy };

  s.addText("Extreme Subgroup Effects", {
    x: 1.0, y: 2.0, w: SW - 2.0, h: 0.9,
    fontSize: 44, fontFace: FONT.header, color: C.white, bold: true,
    align: "left", valign: "top", margin: 0
  });
  s.addText("Under a Uniform Treatment Benefit", {
    x: 1.0, y: 2.9, w: SW - 2.0, h: 0.7,
    fontSize: 30, fontFace: FONT.header, color: C.cream, italic: true,
    align: "left", valign: "top", margin: 0
  });
  s.addText("Section 6.6.4 — High-Risk Subgroup Analysis", {
    x: 1.0, y: 4.3, w: SW - 2.0, h: 0.5,
    fontSize: 18, fontFace: FONT.body, color: C.gold, bold: true,
    align: "left", valign: "top", margin: 0
  });

  s.addText("Larry F. Leon", {
    x: 1.0, y: 6.3, w: SW - 2.0, h: 0.4,
    fontSize: 14, fontFace: FONT.body, color: C.cream,
    align: "left", valign: "top", margin: 0
  });
  s.addText("forestsearch", {
    x: 1.0, y: 6.7, w: SW - 2.0, h: 0.4,
    fontSize: 12, fontFace: FONT.body, color: C.cream, italic: true,
    align: "left", valign: "top", margin: 0
  });
}

// =============================================================================
// SLIDE 2 — Objective
// =============================================================================
function buildSlide2() {
  let s = pres.addSlide();
  addStandardHeader(s, "The Question",
    "Statistical and clinical motivation");

  const cardY = 1.85, cardH = 4.5;
  const colGap = 0.4;
  const colW = (SW - 2 * MARGIN - colGap) / 2;

  // Left: statistical objective
  s.addShape(pres.shapes.RECTANGLE, {
    x: MARGIN, y: cardY, w: colW, h: cardH,
    fill: { color: C.deep }, line: { type: "none" },
    shadow: makeShadow()
  });
  s.addText("STATISTICAL OBJECTIVE", {
    x: MARGIN + 0.4, y: cardY + 0.35, w: colW - 0.8, h: 0.4,
    fontSize: 12, fontFace: FONT.body, color: C.gold, bold: true,
    charSpacing: 4, margin: 0
  });
  s.addText(
    "When a treatment has a perfectly uniform benefit across all patients, how often do standard subgroup analyses produce extreme-looking results by chance alone?",
    { x: MARGIN + 0.4, y: cardY + 0.95, w: colW - 0.8, h: cardH - 1.3,
      fontSize: 19, fontFace: FONT.header, color: C.white,
      valign: "top", margin: 0, paraSpaceAfter: 6 }
  );

  // Right: clinical motivation
  const rightX = MARGIN + colW + colGap;
  s.addShape(pres.shapes.RECTANGLE, {
    x: rightX, y: cardY, w: colW, h: cardH,
    fill: { color: C.teal }, line: { type: "none" },
    shadow: makeShadow()
  });
  s.addText("WHY IT MATTERS", {
    x: rightX + 0.4, y: cardY + 0.35, w: colW - 0.8, h: 0.4,
    fontSize: 12, fontFace: FONT.body, color: C.gold, bold: true,
    charSpacing: 4, margin: 0
  });
  s.addText(
    "Regulatory reviewers see forest plots with apparently alarming subgroup hazard ratios or upper confidence bounds. Are these real signals — or sampling noise from small denominators?",
    { x: rightX + 0.4, y: cardY + 0.95, w: colW - 0.8, h: cardH - 1.3,
      fontSize: 16, fontFace: FONT.body, color: C.white,
      valign: "top", margin: 0, paraSpaceAfter: 6 }
  );

  // Headline footer
  s.addText(
    "Answer: extreme results are routine in small subgroups — clinical rationale provides no statistical protection.",
    { x: MARGIN, y: 6.65, w: SW - 2 * MARGIN, h: 0.45,
      fontSize: 14, fontFace: FONT.body, color: C.navy,
      bold: true, italic: true, align: "center", margin: 0 }
  );
}

// =============================================================================
// SLIDE 3 — Simulation setup
// =============================================================================
function buildSlide3() {
  let s = pres.addSlide();
  addStandardHeader(s, "The Simulation Setup",
    "Real-data DGM with a perfectly uniform treatment effect");

  const steps = [
    { num: "1", title: "Source data",
      body: "GBSG breast cancer trial: N = 686 patients, 299 events" },
    { num: "2", title: "Weibull AFT model",
      body: "Fit to the source data; build a super-population (N = 5,000)" },
    { num: "3", title: "Uniform treatment effect",
      body: "True hazard ratio = 0.70 for every patient — no heterogeneity, by construction" },
    { num: "4", title: "Censoring calibration",
      body: "uniroot-based search matches simulated censoring to the source data" },
    { num: "5", title: "10,000 synthetic trials",
      body: "Each trial: 686 patients × 56 pre-defined subgroups × Cox stratified by grade" }
  ];

  const startY = 1.85, rowH = 0.95, badge = 0.6;
  steps.forEach((st, i) => {
    const y = startY + i * rowH;

    s.addShape(pres.shapes.OVAL, {
      x: MARGIN + 0.1, y: y, w: badge, h: badge,
      fill: { color: i < 3 ? C.deep : C.teal }, line: { type: "none" }
    });
    s.addText(st.num, {
      x: MARGIN + 0.1, y: y, w: badge, h: badge,
      fontSize: 22, fontFace: FONT.header, color: C.white, bold: true,
      align: "center", valign: "middle", margin: 0
    });

    const textX = MARGIN + 0.1 + badge + 0.3;
    s.addText(st.title, {
      x: textX, y: y, w: SW - textX - MARGIN, h: 0.35,
      fontSize: 16, fontFace: FONT.body, color: C.navy, bold: true,
      align: "left", valign: "top", margin: 0
    });
    s.addText(st.body, {
      x: textX, y: y + 0.35, w: SW - textX - MARGIN, h: 0.55,
      fontSize: 13, fontFace: FONT.body, color: C.text,
      align: "left", valign: "top", margin: 0
    });
  });
}

// =============================================================================
// SLIDE 4 — Subgroup taxonomy
// =============================================================================
function buildSlide4() {
  let s = pres.addSlide();
  addStandardHeader(s, "Subgroups Under Study",
    "56 subgroups examined in every simulated trial");

  // Bar chart (left)
  const chartX = MARGIN, chartY = 1.85, chartW = 7.2, chartH = 4.8;
  const chartData = [{
    name: "Count",
    labels: ["Two-way interactions", "Continuous", "3-way interactions",
             "Clinical", "Random benchmarks", "ITT"],
    values: [29, 13, 5, 4, 4, 1]
  }];
  s.addChart(pres.charts.BAR, chartData, {
    x: chartX, y: chartY, w: chartW, h: chartH, barDir: "bar",
    chartColors: [C.teal],
    chartArea: { fill: { color: C.white } },
    catAxisLabelColor: C.text, catAxisLabelFontSize: 11, catAxisLabelFontFace: FONT.body,
    valAxisLabelColor: C.muted, valAxisLabelFontSize: 10,
    valGridLine: { color: "E5E7EB", size: 0.5 },
    catGridLine: { style: "none" },
    showValue: true, dataLabelPosition: "outEnd",
    dataLabelColor: C.navy, dataLabelFontSize: 12, dataLabelFontBold: true,
    showLegend: false,
    showTitle: true, title: "Subgroups by category (total: 56)",
    titleColor: C.navy, titleFontSize: 14, titleFontFace: FONT.header
  });

  // Annotation column (right)
  const annX = MARGIN + chartW + 0.5;
  const annW = SW - annX - MARGIN;

  s.addText("Why so many — and what each adds", {
    x: annX, y: 1.85, w: annW, h: 0.4,
    fontSize: 14, fontFace: FONT.header, color: C.navy, bold: true,
    margin: 0
  });

  const items = [
    { label: "ITT (1)",            body: "Full trial — anchors the comparison" },
    { label: "Clinical (4)",       body: "Standard categorical splits (meno, grade)" },
    { label: "Continuous (13)",    body: "Median / quartile cut-points on age, size, nodes, PGR, ER" },
    { label: "Two-way (29)",       body: "Pairwise clinical x biomarker combinations" },
    { label: "3-way (5)",          body: "Three-way intersections — progressively smaller N" },
    { label: "Random (4)",         body: "Pure-noise benchmarks: random60 / 40 / 20 / 15" }
  ];
  const itemY0 = 2.4, itemH = 0.7;
  items.forEach((it, i) => {
    const y = itemY0 + i * itemH;
    s.addText([
      { text: it.label, options: { bold: true, color: C.deep } },
      { text: "  " + it.body, options: { color: C.text } }
    ], {
      x: annX, y: y, w: annW, h: 0.6,
      fontSize: 12, fontFace: FONT.body, valign: "top", margin: 0
    });
  });
}

// =============================================================================
// SLIDE 5 — Section 6.6.4 intro
// =============================================================================
function buildSlide5() {
  let s = pres.addSlide();
  addStandardHeader(s, "Section 6.6.4 — High-Risk Subgroups",
    "Filtering to subgroups most prone to extreme upper bounds");

  // Filter rule callout
  s.addShape(pres.shapes.RECTANGLE, {
    x: MARGIN, y: 1.85, w: SW - 2 * MARGIN, h: 1.55,
    fill: { color: C.deep }, line: { type: "none" }, shadow: makeShadow()
  });
  s.addText("FILTER RULE", {
    x: MARGIN + 0.4, y: 1.95, w: SW - 2 * MARGIN - 0.8, h: 0.35,
    fontSize: 12, fontFace: FONT.body, color: C.gold, bold: true,
    charSpacing: 4, margin: 0
  });
  s.addText("Pr( UB(HR) >= 2.0 )  >=  10%   across 10,000 simulated trials", {
    x: MARGIN + 0.4, y: 2.35, w: SW - 2 * MARGIN - 0.8, h: 0.6,
    fontSize: 20, fontFace: FONT.header, color: C.white, bold: true,
    margin: 0
  });
  s.addText("All Patients (ITT) included for comparison regardless of whether it meets the filter.", {
    x: MARGIN + 0.4, y: 3.0, w: SW - 2 * MARGIN - 0.8, h: 0.35,
    fontSize: 13, fontFace: FONT.body, color: C.cream, italic: true,
    margin: 0
  });

  // Two info cards
  const colY = 3.75, colH = 3.2, colGap = 0.4;
  const colW = (SW - 2 * MARGIN - colGap) / 2;

  s.addShape(pres.shapes.RECTANGLE, {
    x: MARGIN, y: colY, w: colW, h: colH,
    fill: { color: C.cream }, line: { color: C.teal, width: 0.75 }
  });
  s.addText("Two views, one set of subgroups", {
    x: MARGIN + 0.3, y: colY + 0.2, w: colW - 0.6, h: 0.4,
    fontSize: 14, fontFace: FONT.header, color: C.navy, bold: true,
    margin: 0
  });
  s.addText([
    { text: "UB(HR) panel  ", options: { bold: true, color: C.deep } },
    { text: "asks: how alarming can the upper confidence bound look by chance?", options: { color: C.text, breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "HR panel  ", options: { bold: true, color: C.deep } },
    { text: "asks: how misleading can the point estimate itself be?", options: { color: C.text } }
  ], {
    x: MARGIN + 0.3, y: colY + 0.7, w: colW - 0.6, h: colH - 0.9,
    fontSize: 13, fontFace: FONT.body, valign: "top", margin: 0,
    paraSpaceAfter: 4
  });

  const rightX = MARGIN + colW + colGap;
  s.addShape(pres.shapes.RECTANGLE, {
    x: rightX, y: colY, w: colW, h: colH,
    fill: { color: C.cream }, line: { color: C.teal, width: 0.75 }
  });
  s.addText("11 subgroups + ITT comparator", {
    x: rightX + 0.3, y: colY + 0.2, w: colW - 0.6, h: 0.4,
    fontSize: 14, fontFace: FONT.header, color: C.navy, bold: true,
    margin: 0
  });
  s.addText([
    { text: "Random benchmarks (random15 / 20 / 40 / 60) lead the list — pure sampling noise.", options: { breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "Clinical and biomarker-defined subgroups of similar size are statistically indistinguishable from them.", options: {} }
  ], {
    x: rightX + 0.3, y: colY + 0.7, w: colW - 0.6, h: colH - 0.9,
    fontSize: 13, fontFace: FONT.body, color: C.text, valign: "top",
    margin: 0, paraSpaceAfter: 4
  });
}

// =============================================================================
// SLIDE 6 — UB(HR) distribution figure
// =============================================================================
function buildSlide6() {
  let s = pres.addSlide();
  s.addText("UB(HR) Distribution — Section 6.6.4 Panel A", {
    x: MARGIN, y: 0.3, w: SW - 2 * MARGIN, h: 0.55,
    fontSize: 24, fontFace: FONT.header, color: C.navy, bold: true,
    align: "left", valign: "top", margin: 0
  });
  s.addText("What a safety reviewer sees: the upper 95% confidence bound, across simulations", {
    x: MARGIN, y: 0.85, w: SW - 2 * MARGIN, h: 0.35,
    fontSize: 12, fontFace: FONT.body, color: C.muted, italic: true,
    align: "left", valign: "top", margin: 0
  });

  // Figure (native 1500x750 → 2:1; render at 10" x 5")
  const imgW = 10, imgH = 5;
  const imgX = (SW - imgW) / 2;
  s.addImage({
    path: path.join(ASSETS, "section_6_6_4_fig1.png"),
    x: imgX, y: 1.3, w: imgW, h: imgH
  });

  // Bullets at bottom
  s.addText([
    { text: "random15 / random20 produce UB >= 2 in 75-89% of trials and UB >= 3 in 59-75%.", options: { bullet: true, breakLine: true } },
    { text: "Clinical / biomarker subgroups of similar size (Pre-meno/Age>50, Grade 3/PGR high, ...) show indistinguishable distributions.", options: { bullet: true, breakLine: true } },
    { text: "Arrowheads mark whiskers whose 99th percentile extends beyond the visible axis.", options: { bullet: true } }
  ], {
    x: MARGIN + 0.3, y: 6.4, w: SW - 2 * MARGIN - 0.3, h: SH - 6.4 - 0.3,
    fontSize: 11, fontFace: FONT.body, color: C.text, valign: "top",
    paraSpaceAfter: 3
  });
}

// =============================================================================
// SLIDE 7 — HR distribution figure
// =============================================================================
function buildSlide7() {
  let s = pres.addSlide();
  s.addText("HR Distribution — Section 6.6.4 Panel B", {
    x: MARGIN, y: 0.3, w: SW - 2 * MARGIN, h: 0.55,
    fontSize: 24, fontFace: FONT.header, color: C.navy, bold: true,
    align: "left", valign: "top", margin: 0
  });
  s.addText("Same subgroups, point-estimate view: how misleading can the headline HR be?", {
    x: MARGIN, y: 0.85, w: SW - 2 * MARGIN, h: 0.35,
    fontSize: 12, fontFace: FONT.body, color: C.muted, italic: true,
    align: "left", valign: "top", margin: 0
  });

  const imgW = 10, imgH = 5;
  const imgX = (SW - imgW) / 2;
  s.addImage({
    path: path.join(ASSETS, "section_6_6_4_fig2.png"),
    x: imgX, y: 1.3, w: imgW, h: imgH
  });

  s.addText([
    { text: "Every row's median HR is approximately 0.70 — the Cox model is unbiased even at N = 15.", options: { bullet: true, breakLine: true } },
    { text: "But 27-37% of trials in small subgroups produce HR < 0.5 (looks like a huge benefit); 25-36% produce HR > 1.0 (looks like harm).", options: { bullet: true, breakLine: true } },
    { text: "All this variability is sampling noise — the true HR is 0.70 by construction.", options: { bullet: true } }
  ], {
    x: MARGIN + 0.3, y: 6.4, w: SW - 2 * MARGIN - 0.3, h: SH - 6.4 - 0.3,
    fontSize: 11, fontFace: FONT.body, color: C.text, valign: "top",
    paraSpaceAfter: 3
  });
}

// =============================================================================
// SLIDE 8 — Key takeaways
// =============================================================================
function buildSlide8() {
  let s = pres.addSlide();
  addStandardHeader(s, "Key Takeaways",
    "Reading the UB and HR panels together");

  const startY = 1.85, gridGap = 0.3;
  const colW = (SW - 2 * MARGIN - gridGap) / 2;
  const rowH = (SH - startY - MARGIN - gridGap) / 2;

  const cards = [
    { x: MARGIN, y: startY, color: C.deep, num: "1",
      title: "Extremes are routine in small subgroups",
      body: "Under a uniform HR = 0.70, random15 trials produce UB(HR) > 5 about three-quarters of the time. Same true effect for everyone — the apparent heterogeneity is sampling noise.",
      bodyColor: C.cream },
    { x: MARGIN + colW + gridGap, y: startY, color: C.teal, num: "2",
      title: "Clinical rationale != protection",
      body: "Pre-meno/Age>50 (N approx 38) behaves like random40 (N approx 40). A defined subgroup of size N has the same chance distribution as a random pick of N.",
      bodyColor: C.cream },
    { x: MARGIN, y: startY + rowH + gridGap, color: C.navy, num: "3",
      title: "The Cox model is working correctly",
      body: "Median HR is approximately 0.70 across every row. The wide CIs accurately reflect uncertainty — the problem is interpretation, not the statistics.",
      bodyColor: C.cream },
    { x: MARGIN + colW + gridGap, y: startY + rowH + gridGap, color: C.gold, num: "4",
      title: "Simulation-based calibration is the fix",
      body: "Without a null benchmark, no reviewer can quantify how 'extreme' UB = 3.5 in a 40-patient subgroup actually is. The random* anchors provide that calibration.",
      bodyColor: C.navy }
  ];

  cards.forEach(c => {
    const numColor = c.color === C.gold ? C.navy : C.gold;
    s.addShape(pres.shapes.RECTANGLE, {
      x: c.x, y: c.y, w: colW, h: rowH,
      fill: { color: c.color }, line: { type: "none" },
      shadow: makeShadow()
    });
    s.addText(c.num, {
      x: c.x + 0.3, y: c.y + 0.2, w: 0.7, h: 0.7,
      fontSize: 36, fontFace: FONT.header, color: numColor, bold: true,
      align: "left", valign: "top", margin: 0
    });
    s.addText(c.title, {
      x: c.x + 1.0, y: c.y + 0.3, w: colW - 1.2, h: 0.55,
      fontSize: 16, fontFace: FONT.header, color: c.color === C.gold ? C.navy : C.white, bold: true,
      align: "left", valign: "top", margin: 0
    });
    s.addText(c.body, {
      x: c.x + 0.4, y: c.y + 1.05, w: colW - 0.8, h: rowH - 1.2,
      fontSize: 13, fontFace: FONT.body, color: c.bodyColor,
      align: "left", valign: "top", margin: 0, paraSpaceAfter: 4
    });
  });
}

// =============================================================================
// SLIDE 9 — Resources
// =============================================================================
function buildSlide9() {
  let s = pres.addSlide();
  s.background = { color: C.navy };

  s.addText("Resources", {
    x: 1.0, y: 1.6, w: SW - 2.0, h: 0.8,
    fontSize: 36, fontFace: FONT.header, color: C.white, bold: true,
    align: "left", valign: "top", margin: 0
  });

  const items = [
    { label: "Package",       url: "github.com/larry-leon/forestsearch" },
    { label: "Documentation", url: "larry-leon.github.io/forestsearch/" },
    { label: "Vignette",      url: "larry-leon.github.io/forestsearch/articles/extreme_subgroups.html" },
    { label: "Install",       url: 'pak::pak("larry-leon/forestsearch")' }
  ];
  items.forEach((it, i) => {
    const y = 3.3 + i * 0.65;
    s.addText(it.label, {
      x: 1.0, y: y, w: 3.5, h: 0.45,
      fontSize: 16, fontFace: FONT.body, color: C.gold, bold: true,
      align: "left", valign: "top", margin: 0
    });
    s.addText(it.url, {
      x: 4.5, y: y, w: 8.5, h: 0.45,
      fontSize: 15, fontFace: FONT.body, color: C.cream,
      align: "left", valign: "top", margin: 0
    });
  });

  s.addText("Larry F. Leon  -  forestsearch", {
    x: 1.0, y: 6.6, w: 8, h: 0.4,
    fontSize: 12, fontFace: FONT.body, color: C.cream, italic: true,
    margin: 0
  });
}

// ── Build all & write ──────────────────────────────────────────────────────
buildSlide1();
buildSlide2();
buildSlide3();
buildSlide4();
buildSlide5();
buildSlide6();
buildSlide7();
buildSlide8();
buildSlide9();

const outPath = path.join(__dirname, "extreme_subgroups_section_6_6_4.pptx");
pres.writeFile({ fileName: outPath })
  .then(name => console.log("Wrote:", name));
