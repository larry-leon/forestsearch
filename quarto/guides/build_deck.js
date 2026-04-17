// ============================================================================
// Extreme Subgroup Effects — PowerPoint Deck
// Based on Larry F. León's forestsearch vignette
// ============================================================================

const pptxgen = require("pptxgenjs");
const pres = new pptxgen();

pres.layout = "LAYOUT_WIDE";           // 13.3" × 7.5"
pres.author = "Larry F. León";
pres.title  = "Extreme Subgroup Effects Under a Uniform Treatment Benefit";

// ----------------------------------------------------------------------------
// Palette  —  "Ocean Gradient" with a clinical accent
// ----------------------------------------------------------------------------
const C = {
  navy:     "0C2344",   // primary dark
  deep:     "065A82",   // deep blue
  teal:     "1C7293",   // mid teal
  midnight: "21295C",
  ice:      "E8F0F7",   // very light blue
  pale:     "F5F8FC",   // near white
  white:    "FFFFFF",
  ink:      "1A2A3A",   // body text
  muted:    "5C6F83",   // secondary text
  line:     "C9D3DD",   // hairlines
  accent:   "C0392B",   // accent red (true-HR line, alerts)
  gold:     "D4A017",   // highlight
  green:    "1E8449",   // positive / go
  orange:   "E67E22",   // random-benchmark / caution
  rand:     "E67E22",
  clinical: "1C7293",
  itt:      "0C2344"
};

// Geometric / spacing
const SW = 13.3, SH = 7.5;
const MARGIN = 0.55;

// ----------------------------------------------------------------------------
// Helper: title band
// ----------------------------------------------------------------------------
function titleBand(slide, title, eyebrow) {
  // thin left accent strip
  slide.addShape(pres.shapes.RECTANGLE, {
    x: 0, y: 0, w: 0.12, h: SH,
    fill: { color: C.deep }, line: { color: C.deep, width: 0 }
  });
  if (eyebrow) {
    slide.addText(eyebrow.toUpperCase(), {
      x: MARGIN, y: 0.35, w: SW - 2*MARGIN, h: 0.32,
      fontFace: "Calibri", fontSize: 11, bold: true,
      color: C.teal, charSpacing: 4, margin: 0
    });
  }
  slide.addText(title, {
    x: MARGIN, y: eyebrow ? 0.67 : 0.45, w: SW - 2*MARGIN, h: 0.75,
    fontFace: "Calibri", fontSize: 30, bold: true,
    color: C.navy, margin: 0
  });
}

// Helper: page footer (subtle)
function footer(slide, pageNum, totalPages) {
  slide.addText(
    `forestsearch  ·  Extreme Subgroup Effects  ·  L. F. León`,
    { x: MARGIN, y: SH - 0.35, w: 8, h: 0.25,
      fontFace: "Calibri", fontSize: 9, color: C.muted, margin: 0 }
  );
  slide.addText(
    `${pageNum} / ${totalPages}`,
    { x: SW - 1.2, y: SH - 0.35, w: 0.6, h: 0.25,
      fontFace: "Calibri", fontSize: 9, color: C.muted,
      align: "right", margin: 0 }
  );
}

// Helper: a "callout card" with coloured left stripe
function card(slide, opts) {
  const { x, y, w, h, stripe = C.deep, fill = C.pale, pad = 0.2 } = opts;
  slide.addShape(pres.shapes.RECTANGLE, {
    x, y, w, h, fill: { color: fill }, line: { color: C.line, width: 0.5 }
  });
  slide.addShape(pres.shapes.RECTANGLE, {
    x, y, w: 0.08, h,
    fill: { color: stripe }, line: { color: stripe, width: 0 }
  });
  return { innerX: x + pad + 0.1, innerY: y + pad, innerW: w - pad*2 - 0.1, innerH: h - pad*2 };
}

// Helper: numbered step badge (circle w/ number)
function stepBadge(slide, x, y, n, color = C.deep) {
  slide.addShape(pres.shapes.OVAL, {
    x, y, w: 0.55, h: 0.55,
    fill: { color: color }, line: { color: color, width: 0 }
  });
  slide.addText(String(n), {
    x, y, w: 0.55, h: 0.55,
    fontFace: "Calibri", fontSize: 20, bold: true,
    color: C.white, align: "center", valign: "middle", margin: 0
  });
}

// ============================================================================
// SLIDE 1 — Title slide
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.navy };

  // big accent block on right
  s.addShape(pres.shapes.RECTANGLE, {
    x: SW - 4.2, y: 0, w: 4.2, h: SH,
    fill: { color: C.deep }, line: { color: C.deep, width: 0 }
  });
  s.addShape(pres.shapes.RECTANGLE, {
    x: SW - 4.2, y: 0, w: 0.12, h: SH,
    fill: { color: C.gold }, line: { color: C.gold, width: 0 }
  });

  // tiny eyebrow
  s.addText("FORESTSEARCH  ·  VIGNETTE", {
    x: 0.8, y: 0.9, w: 8, h: 0.35,
    fontFace: "Calibri", fontSize: 12, bold: true, color: C.gold,
    charSpacing: 6, margin: 0
  });

  // main title
  s.addText("Extreme Subgroup Effects", {
    x: 0.8, y: 1.35, w: 8, h: 1.1,
    fontFace: "Calibri", fontSize: 46, bold: true,
    color: C.white, margin: 0
  });
  s.addText("Under a Uniform Treatment Benefit", {
    x: 0.8, y: 2.45, w: 8, h: 0.8,
    fontFace: "Calibri", fontSize: 30,
    color: C.white, italic: true, margin: 0
  });

  // divider line
  s.addShape(pres.shapes.LINE, {
    x: 0.8, y: 3.5, w: 1.5, h: 0,
    line: { color: C.gold, width: 3 }
  });

  // subtitle
  s.addText("A Simulation Study Based on Observed Trial Data", {
    x: 0.8, y: 3.75, w: 8, h: 0.5,
    fontFace: "Calibri", fontSize: 18,
    color: C.ice, margin: 0
  });

  // "pull-quote" style tagline on right
  s.addText([
    { text: "Before concluding a subgroup finding is real, compare it against what chance alone produces for a subgroup of ",
      options: { color: C.white, fontSize: 16 } },
    { text: "the same size.",
      options: { color: C.gold, fontSize: 16, bold: true, italic: true } }
  ], {
    x: SW - 3.9, y: 2.7, w: 3.4, h: 3.3,
    fontFace: "Georgia", italic: true,
    valign: "top", margin: 0
  });

  // author + affiliation
  s.addText("Larry F. León", {
    x: 0.8, y: SH - 1.6, w: 6, h: 0.4,
    fontFace: "Calibri", fontSize: 18, bold: true, color: C.white, margin: 0
  });
  s.addText("forestsearch R package  ·  github.com/larry-leon/forestsearch", {
    x: 0.8, y: SH - 1.2, w: 8, h: 0.3,
    fontFace: "Calibri", fontSize: 12, color: C.ice, margin: 0
  });

  // date
  s.addText("April 2026", {
    x: 0.8, y: SH - 0.8, w: 6, h: 0.3,
    fontFace: "Calibri", fontSize: 11, color: C.gold,
    charSpacing: 4, margin: 0
  });
}

// ============================================================================
// SLIDE 2 — The Central Question
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleBand(s, "The Central Question", "Primary Question");

  // Framing question — big pull-quote
  const q = card(s, { x: MARGIN, y: 1.55, w: SW - 2*MARGIN, h: 1.35,
                      stripe: C.accent, fill: C.pale });
  s.addText("Does an apparent subgroup finding reflect", {
    x: q.innerX, y: q.innerY + 0.05, w: q.innerW, h: 0.42,
    fontFace: "Georgia", fontSize: 20, italic: true, color: C.ink, margin: 0
  });
  s.addText([
    { text: "genuine heterogeneity", options: { bold: true, color: C.accent } },
    { text: "  or is it  ", options: { color: C.ink } },
    { text: "noise amplified by small sample size", options: { bold: true, color: C.accent } },
    { text: "?", options: { color: C.ink } }
  ], {
    x: q.innerX, y: q.innerY + 0.5, w: q.innerW, h: 0.7,
    fontFace: "Georgia", fontSize: 22, italic: true, margin: 0
  });

  // Three framing cards: What reviewers see / Why it matters / What chance does
  const cardY = 3.15;
  const cardH = 3.55;
  const cardW = (SW - 2*MARGIN - 0.3*2) / 3;

  const cards = [
    {
      stripe: C.deep,
      title: "What reviewers see",
      body: [
        "Forest plots with dozens of subgroup HRs and 95 % CIs.",
        "Small subgroups routinely show dramatic point estimates and wide intervals.",
        "Country-level analyses, biomarker strata, and post-hoc cuts amplify the problem."
      ]
    },
    {
      stripe: C.teal,
      title: "Why it matters",
      body: [
        "Regulatory and HTA reviews increasingly scrutinize subgroups.",
        "Apparent \u201cno-benefit\u201d or \u201charm\u201d signals in small groups can alter approval decisions.",
        "Clinical rationale does not automatically translate into statistical protection."
      ]
    },
    {
      stripe: C.gold,
      title: "What chance alone does",
      body: [
        "Even when the treatment helps every patient equally, small subgroups produce extreme estimates.",
        "Wide CIs in small N are the Cox model working correctly.",
        "Without calibration, \u201cextreme\u201d cannot be interpreted."
      ]
    }
  ];
  cards.forEach((c, i) => {
    const x = MARGIN + i * (cardW + 0.3);
    const inner = card(s, { x, y: cardY, w: cardW, h: cardH, stripe: c.stripe, fill: C.pale });
    s.addText(c.title, {
      x: inner.innerX, y: inner.innerY, w: inner.innerW, h: 0.45,
      fontFace: "Calibri", fontSize: 16, bold: true, color: c.stripe, margin: 0
    });
    s.addText(
      c.body.map((line, idx) => ({
        text: line,
        options: { bullet: { code: "25A0" }, breakLine: idx < c.body.length - 1, color: C.ink }
      })),
      { x: inner.innerX, y: inner.innerY + 0.55, w: inner.innerW, h: inner.innerH - 0.55,
        fontFace: "Calibri", fontSize: 13, paraSpaceAfter: 8, margin: 0 }
    );
  });
}

// ============================================================================
// SLIDE 3 — Strategy: use the trial data itself
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleBand(s, "Strategy: Let the Trial Data Define the Null", "Strategy");

  // Left: conceptual diagram (three stacked nodes with arrows)
  const colX = MARGIN + 0.1;
  const boxW = 4.4, boxH = 1.3, gap = 0.3;

  const blocks = [
    { label: "Observed trial data",
      sub: "GBSG  ·  N = 686  ·  all baseline covariates, censoring, outcomes",
      color: C.deep },
    { label: "Weibull AFT model + super-population",
      sub: "generate_aft_dgm_flex()  ·  fits the source structure",
      color: C.teal },
    { label: "Synthetic trials with a uniform benefit",
      sub: "simulate_from_dgm()  ·  true HR = 0.70 for every patient",
      color: C.gold }
  ];

  blocks.forEach((b, i) => {
    const y = 1.65 + i * (boxH + gap);
    s.addShape(pres.shapes.RECTANGLE, {
      x: colX, y, w: boxW, h: boxH,
      fill: { color: C.white }, line: { color: b.color, width: 2 }
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x: colX, y, w: 0.15, h: boxH,
      fill: { color: b.color }, line: { color: b.color, width: 0 }
    });
    s.addText(b.label, {
      x: colX + 0.3, y: y + 0.18, w: boxW - 0.4, h: 0.45,
      fontFace: "Calibri", fontSize: 16, bold: true, color: C.navy, margin: 0
    });
    s.addText(b.sub, {
      x: colX + 0.3, y: y + 0.65, w: boxW - 0.4, h: 0.55,
      fontFace: "Calibri", fontSize: 12, color: C.muted, italic: true, margin: 0
    });

    // arrow to next box
    if (i < blocks.length - 1) {
      const ax = colX + boxW / 2 - 0.12;
      const ay = y + boxH + 0.03;
      s.addShape(pres.shapes.RECTANGLE, {
        x: ax, y: ay, w: 0.24, h: gap - 0.06,
        fill: { color: C.teal }, line: { color: C.teal, width: 0 }
      });
    }
  });

  // Right: explanatory text
  const rx = colX + boxW + 0.6;
  const rw = SW - rx - MARGIN;

  s.addText("The data-generating mechanism (DGM) is built from the trial itself", {
    x: rx, y: 1.65, w: rw, h: 0.45,
    fontFace: "Calibri", fontSize: 18, bold: true, color: C.navy, margin: 0
  });

  s.addText([
    { text: "Observed covariate distributions, censoring process, and individual patient characteristics all carry over into the simulation.",
      options: { breakLine: true, color: C.ink } },
    { text: "", options: { breakLine: true } },
    { text: "Treatment effect is overridden to be ",
      options: { color: C.ink } },
    { text: "perfectly uniform",
      options: { bold: true, color: C.accent, breakLine: true } },
    { text: "   \u03C8\u2070(L) = log(HR\u2080)  for all L.",
      options: { fontFace: "Cambria", italic: true, color: C.midnight, breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "Any apparent subgroup finding is, by construction, ", options: { color: C.ink } },
    { text: "a false positive", options: { bold: true, color: C.accent } },
    { text: ".", options: { color: C.ink } }
  ], {
    x: rx, y: 2.2, w: rw, h: 2.5,
    fontFace: "Calibri", fontSize: 14, margin: 0
  });

  // Bottom: key idea callout
  const cal = card(s, { x: rx, y: 5.0, w: rw, h: 1.85, stripe: C.green, fill: C.pale });
  s.addText("Key design feature", {
    x: cal.innerX, y: cal.innerY, w: cal.innerW, h: 0.4,
    fontFace: "Calibri", fontSize: 14, bold: true, color: C.green, margin: 0
  });
  s.addText(
    "The DGM is study-specific, not generic. Covariate correlations, event rates, " +
    "censoring patterns — everything that makes real subgroups correlated with " +
    "each other — is preserved. This is what separates the framework from " +
    "abstract power calculations.",
    { x: cal.innerX, y: cal.innerY + 0.4, w: cal.innerW, h: cal.innerH - 0.4,
      fontFace: "Calibri", fontSize: 13, color: C.ink, margin: 0 }
  );
}

// ============================================================================
// SLIDE 4 — The gbsg case study (dataset summary)
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleBand(s, "Case Study: the GBSG Breast Cancer Trial", "Source Data");

  // Left column: key stats as big callouts
  const lx = MARGIN, lw = 5.3;

  const stats = [
    { big: "686",     label: "patients with complete prognostic data", size: 40, w: 2.0 },
    { big: "299",     label: "recurrence-free survival events",        size: 40, w: 2.0 },
    { big: "~33",     label: "months median event time",               size: 40, w: 2.0 },
    { big: "1984–89", label: "enrollment period  ·  Schumacher et al., JCO 1994", size: 26, w: 2.4 }
  ];

  stats.forEach((st, i) => {
    const y = 1.65 + i * 1.15;
    s.addShape(pres.shapes.RECTANGLE, {
      x: lx, y, w: lw, h: 1.02,
      fill: { color: C.pale }, line: { color: C.line, width: 0.5 }
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x: lx, y, w: 0.08, h: 1.02,
      fill: { color: C.deep }, line: { color: C.deep, width: 0 }
    });
    s.addText(st.big, {
      x: lx + 0.25, y: y + 0.05, w: st.w, h: 0.9,
      fontFace: "Calibri", fontSize: st.size, bold: true, color: C.navy,
      valign: "middle", margin: 0
    });
    s.addText(st.label, {
      x: lx + 0.25 + st.w + 0.1, y: y + 0.05, w: lw - 0.4 - st.w, h: 0.9,
      fontFace: "Calibri", fontSize: 13, color: C.ink,
      valign: "middle", margin: 0
    });
  });

  // Right column: variables table-ish
  const rx = lx + lw + 0.4, rw = SW - rx - MARGIN;
  s.addText("Key variables", {
    x: rx, y: 1.55, w: rw, h: 0.4,
    fontFace: "Calibri", fontSize: 16, bold: true, color: C.teal, margin: 0
  });

  const vars = [
    ["hormon",    "treatment: hormonal therapy (0 / 1)"],
    ["rfstime",   "recurrence-free survival time (days)"],
    ["status",    "event indicator (0 = censored, 1 = event)"],
    ["age",       "age in years"],
    ["meno",      "menopausal status (0 = pre, 1 = post)"],
    ["size",      "tumour size (mm)"],
    ["grade",     "tumour grade (1, 2, 3)"],
    ["nodes",     "number of positive nodes"],
    ["pgr",       "progesterone receptor (fmol)"],
    ["er",        "oestrogen receptor (fmol)"]
  ];

  const tableRows = vars.map((v) => [
    { text: v[0], options: { fontFace: "Consolas", fontSize: 11, color: C.deep, bold: true } },
    { text: v[1], options: { fontFace: "Calibri", fontSize: 11, color: C.ink } }
  ]);
  s.addTable(tableRows, {
    x: rx, y: 2.0, w: rw,
    colW: [1.3, rw - 1.3],
    rowH: 0.32,
    border: { type: "solid", pt: 0.5, color: C.line },
    fill: { color: C.white },
    align: "left", valign: "middle",
    margin: 0.08
  });

  // Source note
  s.addText(
    "Randomised 2×2 trial of hormonal therapy and chemotherapy duration in node-positive breast cancer.  " +
    "In this vignette we use hormonal therapy as the treatment factor.  " +
    "Available as  survival::gbsg.",
    { x: MARGIN, y: SH - 0.75, w: SW - 2*MARGIN, h: 0.4,
      fontFace: "Calibri", fontSize: 10, italic: true, color: C.muted, margin: 0 }
  );
}

// ============================================================================
// SLIDE 5 — Setting: a perfectly uniform treatment effect
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleBand(s, "Setting: a Perfectly Uniform Treatment Effect", "Null DGM");

  // Equation band
  s.addShape(pres.shapes.RECTANGLE, {
    x: MARGIN, y: 1.55, w: SW - 2*MARGIN, h: 1.2,
    fill: { color: C.navy }, line: { color: C.navy, width: 0 }
  });
  s.addText("True log hazard ratio is constant across every patient:", {
    x: MARGIN + 0.3, y: 1.65, w: SW - 2*MARGIN - 0.6, h: 0.4,
    fontFace: "Calibri", fontSize: 14, color: C.ice, italic: true, margin: 0
  });
  s.addText([
    { text: "\u03C8\u2070(L)", options: { fontFace: "Cambria", italic: true, color: C.white } },
    { text: "  =  ", options: { color: C.white } },
    { text: "log(HR", options: { color: C.white } },
    { text: "true", options: { color: C.white, subscript: true } },
    { text: ")", options: { color: C.white } },
    { text: "     for every covariate profile  ", options: { color: C.ice, fontSize: 18 } },
    { text: "L", options: { fontFace: "Cambria", italic: true, color: C.ice, fontSize: 18 } },
    { text: ".", options: { color: C.ice, fontSize: 18 } }
  ], {
    x: MARGIN + 0.3, y: 2.05, w: SW - 2*MARGIN - 0.6, h: 0.55,
    fontFace: "Cambria", fontSize: 24, bold: true, margin: 0
  });

  // Three column explanation
  const y0 = 3.1;
  const h0 = 3.45;
  const w0 = (SW - 2*MARGIN - 0.3 * 2) / 3;

  const cols = [
    {
      title: "HR  =  0.70",
      stripe: C.accent,
      body: "A moderate treatment benefit — identical for every subject. " +
            "Specified via the spline anchors:\n\nlog.hrs = log(c(0.70, 0.70, 0.70))"
    },
    {
      title: "No embedded subgroup",
      stripe: C.deep,
      body: "model = \"null\" in generate_aft_dgm_flex().\n\n" +
            "There is no harm subgroup, no effect modifier, no heterogeneity. " +
            "Every subgroup shares exactly the same HR."
    },
    {
      title: "Any finding is a false positive",
      stripe: C.gold,
      body: "By construction. What we observe across simulated trials is the " +
            "pure sampling-variability distribution of subgroup estimates."
    }
  ];
  cols.forEach((c, i) => {
    const x = MARGIN + i * (w0 + 0.3);
    const inner = card(s, { x, y: y0, w: w0, h: h0, stripe: c.stripe, fill: C.pale });
    s.addText(c.title, {
      x: inner.innerX, y: inner.innerY, w: inner.innerW, h: 0.5,
      fontFace: "Calibri", fontSize: 17, bold: true, color: c.stripe, margin: 0
    });
    s.addText(c.body, {
      x: inner.innerX, y: inner.innerY + 0.55, w: inner.innerW, h: inner.innerH - 0.6,
      fontFace: (c.title.indexOf("=") >= 0 || c.body.indexOf("log.hrs") >= 0) ? "Calibri" : "Calibri",
      fontSize: 13, color: C.ink, margin: 0
    });
  });
}

// ============================================================================
// SLIDE 6 — Installation & accessing the vignette
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleBand(s, "Installing forestsearch & Accessing the Vignette", "Getting Started");

  // Install block (dark)
  const ix = MARGIN, iy = 1.55, iw = 6.25, ih = 2.6;
  s.addShape(pres.shapes.RECTANGLE, {
    x: ix, y: iy, w: iw, h: ih,
    fill: { color: C.navy }, line: { color: C.navy, width: 0 }
  });
  s.addText("Installation", {
    x: ix + 0.3, y: iy + 0.15, w: iw - 0.6, h: 0.4,
    fontFace: "Calibri", fontSize: 16, bold: true, color: C.gold, margin: 0
  });
  s.addText("1.  Install pak (one time)", {
    x: ix + 0.3, y: iy + 0.55, w: iw - 0.6, h: 0.3,
    fontFace: "Calibri", fontSize: 11, color: C.ice, italic: true, margin: 0
  });
  s.addText('install.packages("pak")', {
    x: ix + 0.3, y: iy + 0.82, w: iw - 0.6, h: 0.32,
    fontFace: "Consolas", fontSize: 13, color: C.white, margin: 0
  });
  s.addText("2.  Install forestsearch from GitHub", {
    x: ix + 0.3, y: iy + 1.2, w: iw - 0.6, h: 0.3,
    fontFace: "Calibri", fontSize: 11, color: C.ice, italic: true, margin: 0
  });
  s.addText('pak::pak("larry-leon/forestsearch")', {
    x: ix + 0.3, y: iy + 1.47, w: iw - 0.6, h: 0.32,
    fontFace: "Consolas", fontSize: 13, color: C.white, margin: 0
  });
  s.addText("3.  Load the package", {
    x: ix + 0.3, y: iy + 1.85, w: iw - 0.6, h: 0.3,
    fontFace: "Calibri", fontSize: 11, color: C.ice, italic: true, margin: 0
  });
  s.addText('library(forestsearch)', {
    x: ix + 0.3, y: iy + 2.12, w: iw - 0.6, h: 0.32,
    fontFace: "Consolas", fontSize: 13, color: C.white, margin: 0
  });

  // Docs + vignette (right)
  const rx = ix + iw + 0.4, rw = SW - rx - MARGIN;
  s.addText("Documentation site", {
    x: rx, y: iy, w: rw, h: 0.4,
    fontFace: "Calibri", fontSize: 16, bold: true, color: C.teal, margin: 0
  });
  s.addText([
    { text: "larry-leon.github.io/forestsearch",
      options: { hyperlink: { url: "https://larry-leon.github.io/forestsearch/" },
                 color: C.deep, bold: true } }
  ], {
    x: rx, y: iy + 0.4, w: rw, h: 0.4,
    fontFace: "Consolas", fontSize: 14, margin: 0
  });
  s.addText("includes function reference and vignettes.", {
    x: rx, y: iy + 0.75, w: rw, h: 0.3,
    fontFace: "Calibri", fontSize: 11, color: C.muted, italic: true, margin: 0
  });

  s.addText("Extreme Subgroups vignette", {
    x: rx, y: iy + 1.15, w: rw, h: 0.4,
    fontFace: "Calibri", fontSize: 16, bold: true, color: C.teal, margin: 0
  });
  s.addText([
    { text: "larry-leon.github.io/forestsearch/articles/",
      options: { color: C.deep, bold: true, breakLine: true } },
    { text: "extreme_subgroups.html",
      options: { color: C.deep, bold: true,
                 hyperlink: { url: "https://larry-leon.github.io/forestsearch/articles/extreme_subgroups.html" } } }
  ], {
    x: rx, y: iy + 1.55, w: rw, h: 0.7,
    fontFace: "Consolas", fontSize: 12, margin: 0
  });

  // Bottom: companion vignettes
  const cy = 4.4, ch = 2.55;
  s.addText("Related vignettes on the pkgdown site", {
    x: MARGIN, y: cy, w: SW - 2*MARGIN, h: 0.4,
    fontFace: "Calibri", fontSize: 16, bold: true, color: C.navy, margin: 0
  });
  const items = [
    ["Getting Started",        "End-to-end walkthrough with the GBSG trial"],
    ["Methodology",            "Statistical framework, algorithm details, simulation results"],
    ["Treatment Effect Definitions", "Marginal HR, AHR, and CDE estimands"],
    ["Simulation Studies",     "Operating characteristics and benchmarks"]
  ];
  items.forEach((it, i) => {
    const col = i % 2, row = Math.floor(i / 2);
    const x = MARGIN + col * ((SW - 2*MARGIN) / 2 + 0.1);
    const y = cy + 0.55 + row * 0.85;
    const w = (SW - 2*MARGIN) / 2 - 0.15;
    s.addShape(pres.shapes.RECTANGLE, {
      x, y, w, h: 0.75,
      fill: { color: C.pale }, line: { color: C.line, width: 0.5 }
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x, y, w: 0.06, h: 0.75,
      fill: { color: C.teal }, line: { color: C.teal, width: 0 }
    });
    s.addText(it[0], {
      x: x + 0.18, y: y + 0.08, w: w - 0.25, h: 0.32,
      fontFace: "Calibri", fontSize: 13, bold: true, color: C.navy, margin: 0
    });
    s.addText(it[1], {
      x: x + 0.18, y: y + 0.38, w: w - 0.25, h: 0.32,
      fontFace: "Calibri", fontSize: 11, color: C.muted, margin: 0
    });
  });
}

// ============================================================================
// SLIDE 7 — The simulation workflow (8-step pipeline)
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleBand(s, "Simulation Workflow", "Pipeline");

  const steps = [
    { n: 1, title: "Prepare data",      fn: "—",                            purpose: "Convert source data to a consistent time scale (e.g. days → months)" },
    { n: 2, title: "Build the DGM",     fn: "generate_aft_dgm_flex()",      purpose: "Fit a Weibull AFT model; create a super-population (N = 5 000)" },
    { n: 3, title: "Diagnose censoring",fn: "check_censoring_dgm()",        purpose: "Compare simulated vs. observed censoring rate and KM median" },
    { n: 4, title: "Calibrate censoring",fn: "calibrate_cens_adjust()",     purpose: "Find a log-scale shift so simulated censoring matches the trial" },
    { n: 5, title: "Simulate trials",   fn: "simulate_from_dgm()",          purpose: "Draw synthetic RCTs from the calibrated DGM (one call per trial)" },
    { n: 6, title: "Fit Cox in subgroups", fn: "coxph()  (stratified)",     purpose: "For each pre-defined subgroup, store HR point estimate and UB(HR)" },
    { n: 7, title: "Summarise",         fn: "summarize_extreme_sims()",     purpose: "Empirical 99 % CI (ECI), Pr(UB ≥ 1), conditional UB quantiles" },
    { n: 8, title: "Visualise",         fn: "plot_hr_forest() / plot_ub_forest()", purpose: "Forest plots with random-benchmark anchors by subgroup size" }
  ];

  const rowH = 0.56;
  const y0 = 1.75;
  steps.forEach((st, i) => {
    const y = y0 + i * rowH;

    // Step badge
    stepBadge(s, MARGIN + 0.05, y + 0.01, st.n, i < 4 ? C.deep : C.teal);

    // title
    s.addText(st.title, {
      x: MARGIN + 0.75, y: y + 0.02, w: 2.4, h: 0.5,
      fontFace: "Calibri", fontSize: 14, bold: true, color: C.navy,
      valign: "middle", margin: 0
    });
    // function
    s.addText(st.fn, {
      x: MARGIN + 3.25, y: y + 0.02, w: 3.3, h: 0.5,
      fontFace: "Consolas", fontSize: 12, color: C.accent,
      valign: "middle", margin: 0
    });
    // purpose
    s.addText(st.purpose, {
      x: MARGIN + 6.65, y: y + 0.02, w: SW - MARGIN - 6.7 - 0.1, h: 0.5,
      fontFace: "Calibri", fontSize: 12, color: C.ink,
      valign: "middle", margin: 0
    });

    // hairline
    if (i < steps.length - 1) {
      s.addShape(pres.shapes.LINE, {
        x: MARGIN + 0.7, y: y + rowH - 0.03, w: SW - 2 * MARGIN - 0.7, h: 0,
        line: { color: C.line, width: 0.5 }
      });
    }
  });

  // column headers
  s.addText("Step", {
    x: MARGIN + 0.75, y: y0 - 0.32, w: 2.4, h: 0.28,
    fontFace: "Calibri", fontSize: 10, bold: true, color: C.muted,
    charSpacing: 4, margin: 0
  });
  s.addText("Function", {
    x: MARGIN + 3.25, y: y0 - 0.32, w: 3.3, h: 0.28,
    fontFace: "Calibri", fontSize: 10, bold: true, color: C.muted,
    charSpacing: 4, margin: 0
  });
  s.addText("Purpose", {
    x: MARGIN + 6.65, y: y0 - 0.32, w: SW - MARGIN - 6.7 - 0.1, h: 0.28,
    fontFace: "Calibri", fontSize: 10, bold: true, color: C.muted,
    charSpacing: 4, margin: 0
  });
}

// ============================================================================
// SLIDE 8 — Subgroups under study (56 total)
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleBand(s, "Subgroups Under Study — 56 in Total", "Study Design");

  // 6 category cards in 2 rows × 3 cols
  const catX = MARGIN;
  const totalW = SW - 2*MARGIN;
  const gap = 0.25;
  const cw = (totalW - gap * 2) / 3;
  const ch = 1.9;
  const startY = 1.55;

  const cats = [
    { color: C.navy,     label: "ITT",          n: "1",  ex: "All Patients  (N \u2248 686)" },
    { color: C.deep,     label: "Clinical",     n: "4",  ex: "Pre-/Post-meno  ·  Grade 3  ·  Grade 1/2" },
    { color: C.teal,     label: "Continuous",   n: "13", ex: "Age  ·  Size  ·  Nodes  ·  PGR  ·  ER (various cuts)" },
    { color: "9B5094",   label: "2-way interactions", n: "29", ex: "Meno \u00D7 Grade / Age / ER  ·  Grade \u00D7 Nodes / ER / PGR  ·  Size \u00D7 Nodes" },
    { color: "B03A7A",   label: "3-way interactions", n: "5",  ex: "Pre-meno/Young/ER-low  ·  Post-meno/G3/Node+  ·  \u2026" },
    { color: C.orange,   label: "Random benchmarks", n: "4",   ex: "random60  ·  random40  ·  random20  ·  random15   (re-sampled each sim)" }
  ];

  cats.forEach((c, i) => {
    const col = i % 3, row = Math.floor(i / 3);
    const x = catX + col * (cw + gap);
    const y = startY + row * (ch + gap);

    s.addShape(pres.shapes.RECTANGLE, {
      x, y, w: cw, h: ch,
      fill: { color: C.white }, line: { color: c.color, width: 1.5 }
    });
    // coloured header strip
    s.addShape(pres.shapes.RECTANGLE, {
      x, y, w: cw, h: 0.55,
      fill: { color: c.color }, line: { color: c.color, width: 0 }
    });

    s.addText(c.label, {
      x: x + 0.2, y: y + 0.05, w: cw - 0.9, h: 0.45,
      fontFace: "Calibri", fontSize: 14, bold: true, color: C.white,
      valign: "middle", margin: 0
    });
    s.addText(c.n, {
      x: x + cw - 0.85, y: y + 0.05, w: 0.7, h: 0.45,
      fontFace: "Calibri", fontSize: 18, bold: true, color: C.white,
      align: "right", valign: "middle", margin: 0
    });

    s.addText(c.ex, {
      x: x + 0.2, y: y + 0.7, w: cw - 0.4, h: ch - 0.8,
      fontFace: "Calibri", fontSize: 12, color: C.ink, margin: 0
    });
  });

  // Random benchmarks explainer card at the bottom
  const by = startY + 2 * (ch + gap) + 0.1;
  const bh = SH - by - 0.7;
  const inner = card(s, { x: MARGIN, y: by, w: SW - 2*MARGIN, h: bh,
                          stripe: C.orange, fill: C.pale });
  s.addText("The random* subgroups are the calibration anchor", {
    x: inner.innerX, y: inner.innerY, w: inner.innerW, h: 0.4,
    fontFace: "Calibri", fontSize: 14, bold: true, color: C.orange,
    valign: "top", margin: 0
  });
  s.addText(
    "They have no clinical meaning whatsoever — their true HR equals the ITT HR of 0.70. " +
    "Any extreme finding is pure sampling variability. Placing them alongside clinical " +
    "subgroups of similar size makes the noise floor explicit for each size class.",
    { x: inner.innerX, y: inner.innerY + 0.5, w: inner.innerW, h: inner.innerH - 0.5,
      fontFace: "Calibri", fontSize: 12, color: C.ink, valign: "top", margin: 0 }
  );
}

// ============================================================================
// SLIDE 9 — What is measured per subgroup per simulation
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleBand(s, "What Is Measured Per Subgroup, Per Simulation", "Metrics");

  // Cox model specification on top
  const tx = MARGIN, tw = SW - 2*MARGIN;
  s.addShape(pres.shapes.RECTANGLE, {
    x: tx, y: 1.55, w: tw, h: 0.9,
    fill: { color: C.navy }, line: { color: C.navy, width: 0 }
  });
  s.addText("Primary analysis", {
    x: tx + 0.3, y: 1.6, w: tw - 0.6, h: 0.3,
    fontFace: "Calibri", fontSize: 11, bold: true, color: C.gold,
    charSpacing: 4, margin: 0
  });
  s.addText("coxph(Surv(y_sim, event_sim) ~ treat_sim + strata(grade), data = subgroup_data)", {
    x: tx + 0.3, y: 1.9, w: tw - 0.6, h: 0.45,
    fontFace: "Consolas", fontSize: 14, color: C.white, margin: 0
  });

  // Two metric cards
  const my = 2.7;
  const mh = SH - my - 0.8;
  const mw = (SW - 2*MARGIN - 0.4) / 2;

  // HR card
  const hrInner = card(s, { x: MARGIN, y: my, w: mw, h: mh, stripe: C.deep, fill: C.pale });
  s.addText("HR point estimate distribution", {
    x: hrInner.innerX, y: hrInner.innerY, w: hrInner.innerW, h: 0.45,
    fontFace: "Calibri", fontSize: 16, bold: true, color: C.deep, margin: 0
  });
  s.addText([
    { text: "Across all simulated trials, record the HR estimate for each subgroup.",
      options: { breakLine: true, color: C.ink } },
    { text: "", options: { breakLine: true } },
    { text: "Empirical 99 % CI (ECI)",
      options: { bold: true, color: C.navy, breakLine: true } },
    { text: "    median + 1st and 99th percentiles",
      options: { color: C.ink, breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "Pr(HR < 0.80)",
      options: { bold: true, color: C.navy, fontFace: "Consolas", breakLine: true } },
    { text: "    frequency of \u201calarming\u201d point estimate",
      options: { color: C.ink, breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "Under truth the median should sit at 0.70 for every subgroup. The spread is what changes.",
      options: { italic: true, color: C.muted } }
  ], {
    x: hrInner.innerX, y: hrInner.innerY + 0.55, w: hrInner.innerW, h: hrInner.innerH - 0.6,
    fontFace: "Calibri", fontSize: 13, margin: 0
  });

  // UB card
  const ubInner = card(s, { x: MARGIN + mw + 0.4, y: my, w: mw, h: mh, stripe: C.accent, fill: C.pale });
  s.addText("Upper 95 % CI bound  UB(HR)", {
    x: ubInner.innerX, y: ubInner.innerY, w: ubInner.innerW, h: 0.45,
    fontFace: "Calibri", fontSize: 16, bold: true, color: C.accent, margin: 0
  });
  s.addText([
    { text: "What a safety reviewer looks at when asking whether harm can be ruled out.",
      options: { breakLine: true, color: C.ink } },
    { text: "", options: { breakLine: true } },
    { text: "Pr(UB < 1.0)",
      options: { bold: true, color: C.navy, fontFace: "Consolas", breakLine: true } },
    { text: "    power to demonstrate benefit; drops sharply with N",
      options: { color: C.ink, breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "Median UB | UB \u2265 1.0",
      options: { bold: true, color: C.navy, breakLine: true } },
    { text: "    worst-case bound when the CI fails to exclude 1.0",
      options: { color: C.ink, breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "The key safety quantity — how large the upper bound can grow by chance alone.",
      options: { italic: true, color: C.muted } }
  ], {
    x: ubInner.innerX, y: ubInner.innerY + 0.55, w: ubInner.innerW, h: ubInner.innerH - 0.6,
    fontFace: "Calibri", fontSize: 13, margin: 0
  });
}

// ============================================================================
// SLIDE 10 — Headline result: UB(HR) by size (chart)
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleBand(s, "What Chance Alone Produces: 99th Percentile of UB(HR)", "Headline Result");

  // Chart: 99th percentile UB by approximate subgroup size
  // Illustrative values from the study's summary table
  const labels = ["ITT (~686)", "~200", "~60", "random40 (~40)", "random20 (~20)", "random15 (~15)"];
  const values = [1.1, 1.5, 2.5, 3.5, 5.0, 5.8];

  s.addChart(pres.charts.BAR, [{
    name: "99th percentile UB(HR)",
    labels: labels,
    values: values
  }], {
    x: MARGIN, y: 1.55, w: 8.3, h: 5.3,
    barDir: "col",
    chartColors: [C.deep],
    chartArea: { fill: { color: C.white }, roundedCorners: false },
    plotArea:  { fill: { color: C.white } },
    catAxisLabelColor: C.ink,
    catAxisLabelFontSize: 11,
    catAxisLabelFontFace: "Calibri",
    valAxisLabelColor: C.muted,
    valAxisLabelFontSize: 10,
    valGridLine: { color: C.line, size: 0.5 },
    catGridLine: { style: "none" },
    showValue: true,
    dataLabelPosition: "outEnd",
    dataLabelColor: C.navy,
    dataLabelFontFace: "Calibri",
    dataLabelFontSize: 11,
    dataLabelFontBold: true,
    dataLabelFormatCode: "0.0",
    showLegend: false,
    showTitle: false,
    valAxisMinVal: 0,
    valAxisMaxVal: 7,
    valAxisMajorUnit: 1
  });

  // Right side explanation
  const rx = 9.2, rw = SW - rx - MARGIN;

  s.addText("Reading the chart", {
    x: rx, y: 1.6, w: rw, h: 0.4,
    fontFace: "Calibri", fontSize: 16, bold: true, color: C.teal, margin: 0
  });
  s.addText(
    "Each bar is the 99th percentile of the Cox 95 % upper confidence bound " +
    "across thousands of simulated trials where the true HR is 0.70 everywhere.",
    { x: rx, y: 2.05, w: rw, h: 1.3,
      fontFace: "Calibri", fontSize: 12, color: C.ink, margin: 0 }
  );

  // Key insight card
  const ki = card(s, { x: rx, y: 3.55, w: rw, h: 2.2, stripe: C.accent, fill: C.pale });
  s.addText("In 1 % of trials", {
    x: ki.innerX, y: ki.innerY, w: ki.innerW, h: 0.35,
    fontFace: "Calibri", fontSize: 12, bold: true, color: C.accent, margin: 0
  });
  s.addText([
    { text: "UB(HR) > 3.5", options: { bold: true, color: C.navy } },
    { text: "  in a 40-patient subgroup.", options: { color: C.ink, breakLine: true } },
    { text: "UB(HR) > 5.0", options: { bold: true, color: C.navy } },
    { text: "  in a 15-patient subgroup.", options: { color: C.ink, breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "…even though the treatment benefits every patient equally.",
      options: { italic: true, color: C.muted } }
  ], {
    x: ki.innerX, y: ki.innerY + 0.35, w: ki.innerW, h: ki.innerH - 0.35,
    fontFace: "Calibri", fontSize: 13, margin: 0
  });

  s.addText(
    "Illustrative quantitative anchors drawn from the null-DGM simulation summary.",
    { x: MARGIN, y: SH - 0.75, w: SW - 2*MARGIN, h: 0.3,
      fontFace: "Calibri", fontSize: 9, italic: true, color: C.muted, margin: 0 }
  );
}

// ============================================================================
// SLIDE 11 — Interpretation: the central message
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleBand(s, "Interpreting Subgroup Findings with the Benchmarks", "Interpretation");

  // Framed statement
  const f = card(s, { x: MARGIN, y: 1.55, w: SW - 2*MARGIN, h: 1.5,
                      stripe: C.gold, fill: C.pale });
  s.addText(
    "Before concluding that a subgroup finding reflects real heterogeneity, " +
    "compare it against what chance alone produces for a subgroup of the same size.",
    { x: f.innerX, y: f.innerY + 0.05, w: f.innerW, h: f.innerH - 0.1,
      fontFace: "Georgia", fontSize: 18, italic: true, color: C.navy,
      valign: "middle", margin: 0 }
  );

  // Three takeaways
  const y0 = 3.3;
  const h0 = 3.4;
  const w0 = (SW - 2*MARGIN - 0.3 * 2) / 3;
  const takes = [
    {
      n: "1",
      color: C.deep,
      title: "Wide CIs are correct",
      body: "The Cox model is working exactly as intended. A 95 % CI from 0.3 to 4.5 " +
            "in a 15-patient subgroup is a faithful summary of the evidence — not a bug."
    },
    {
      n: "2",
      color: C.teal,
      title: "Clinical rationale \u2260 statistical protection",
      body: "Clinical subgroups of N \u2248 40 have the same HR / UB distribution as " +
            "random benchmarks of N \u2248 40. The naming and motivation do not reduce variability."
    },
    {
      n: "3",
      color: C.accent,
      title: "Calibrate, then interpret",
      body: "A UB(HR) of 3.5 in a 40-patient subgroup is unremarkable against random40. " +
            "The same value in a 300-patient subgroup is genuinely extreme. Context is everything."
    }
  ];
  takes.forEach((t, i) => {
    const x = MARGIN + i * (w0 + 0.3);
    const inner = card(s, { x, y: y0, w: w0, h: h0, stripe: t.color, fill: C.pale });
    // big number
    s.addText(t.n, {
      x: inner.innerX, y: inner.innerY, w: 1.0, h: 1.0,
      fontFace: "Georgia", fontSize: 42, bold: true, color: t.color, margin: 0
    });
    s.addText(t.title, {
      x: inner.innerX + 0.9, y: inner.innerY + 0.1, w: inner.innerW - 0.9, h: 0.9,
      fontFace: "Calibri", fontSize: 15, bold: true, color: C.navy,
      valign: "middle", margin: 0
    });
    s.addText(t.body, {
      x: inner.innerX, y: inner.innerY + 1.15, w: inner.innerW, h: inner.innerH - 1.2,
      fontFace: "Calibri", fontSize: 13, color: C.ink, margin: 0
    });
  });
}

// ============================================================================
// SLIDE 12 — Swap in your own trial data (part 1 of 2: Conceptual)
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleBand(s, "Swap in Your Own Trial Data", "Adapting the Framework");

  // Intro band
  s.addText(
    "The framework is designed so that GBSG is an example, not a requirement. " +
    "Any trial with a time-to-event outcome, a treatment indicator, and baseline " +
    "covariates can be plugged in.",
    { x: MARGIN, y: 1.55, w: SW - 2*MARGIN, h: 0.75,
      fontFace: "Calibri", fontSize: 14, color: C.ink, italic: true, margin: 0 }
  );

  // 4 prerequisite tiles
  const pY = 2.45;
  const pH = 1.9;
  const pGap = 0.25;
  const pW = (SW - 2*MARGIN - pGap * 3) / 4;

  const prereqs = [
    { title: "Outcome",   detail: "A time variable and an event indicator (0/1)." },
    { title: "Treatment", detail: "A binary treatment arm variable." },
    { title: "Covariates",detail: "Any continuous or categorical baseline variables — age, biomarkers, stratification factors." },
    { title: "Subgroups", detail: "Your pre-defined subgroups of interest, expressed as logical / subset conditions." }
  ];

  prereqs.forEach((p, i) => {
    const x = MARGIN + i * (pW + pGap);
    s.addShape(pres.shapes.RECTANGLE, {
      x, y: pY, w: pW, h: pH,
      fill: { color: C.pale }, line: { color: C.line, width: 0.5 }
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x, y: pY, w: pW, h: 0.45,
      fill: { color: C.teal }, line: { color: C.teal, width: 0 }
    });
    s.addText(p.title, {
      x: x + 0.15, y: pY + 0.05, w: pW - 0.3, h: 0.35,
      fontFace: "Calibri", fontSize: 13, bold: true, color: C.white,
      valign: "middle", margin: 0
    });
    s.addText(p.detail, {
      x: x + 0.15, y: pY + 0.55, w: pW - 0.3, h: pH - 0.6,
      fontFace: "Calibri", fontSize: 12, color: C.ink, margin: 0
    });
  });

  // Bottom: key substitutions in code
  const sy = pY + pH + 0.35;
  const sh = SH - sy - 0.75;

  s.addShape(pres.shapes.RECTANGLE, {
    x: MARGIN, y: sy, w: SW - 2*MARGIN, h: sh,
    fill: { color: C.navy }, line: { color: C.navy, width: 0 }
  });
  s.addText("What actually changes in the code", {
    x: MARGIN + 0.3, y: sy + 0.1, w: SW - 2*MARGIN - 0.6, h: 0.35,
    fontFace: "Calibri", fontSize: 13, bold: true, color: C.gold,
    charSpacing: 3, margin: 0
  });
  s.addText([
    { text: "data            = ", options: { color: C.ice, fontFace: "Consolas" } },
    { text: "my_trial_data", options: { color: C.gold, fontFace: "Consolas", bold: true } },
    { text: "              # replaces gbsg", options: { color: C.muted, fontFace: "Consolas", breakLine: true } },

    { text: "outcome_var     = ", options: { color: C.ice, fontFace: "Consolas" } },
    { text: "\"time_months\"", options: { color: C.gold, fontFace: "Consolas", bold: true } },
    { text: "                # your time column", options: { color: C.muted, fontFace: "Consolas", breakLine: true } },

    { text: "event_var       = ", options: { color: C.ice, fontFace: "Consolas" } },
    { text: "\"event\"", options: { color: C.gold, fontFace: "Consolas", bold: true } },
    { text: "                      # your event indicator", options: { color: C.muted, fontFace: "Consolas", breakLine: true } },

    { text: "treatment_var   = ", options: { color: C.ice, fontFace: "Consolas" } },
    { text: "\"arm\"", options: { color: C.gold, fontFace: "Consolas", bold: true } },
    { text: "                        # your treatment indicator", options: { color: C.muted, fontFace: "Consolas", breakLine: true } },

    { text: "continuous_vars = ", options: { color: C.ice, fontFace: "Consolas" } },
    { text: "c(\"age\", \"biomarker\", \u2026)", options: { color: C.gold, fontFace: "Consolas", bold: true } },
    { text: "   # your continuous covariates", options: { color: C.muted, fontFace: "Consolas", breakLine: true } },

    { text: "factor_vars     = ", options: { color: C.ice, fontFace: "Consolas" } },
    { text: "c(\"region\", \"ECOG\", \u2026)", options: { color: C.gold, fontFace: "Consolas", bold: true } },
    { text: "     # your factor covariates", options: { color: C.muted, fontFace: "Consolas" } }
  ], {
    x: MARGIN + 0.3, y: sy + 0.5, w: SW - 2*MARGIN - 0.6, h: sh - 0.6,
    fontFace: "Consolas", fontSize: 12, margin: 0
  });
}

// ============================================================================
// SLIDE 13 — Swap-in steps (part 2 of 2: Concrete workflow)
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleBand(s, "Six Concrete Steps to Adapt the Vignette", "Adapting the Framework");

  const steps = [
    { n: 1, title: "Align the time scale",
      body: "Convert event times to the units the DGM will use (days \u2192 months is typical). Time-scale mismatch is the single most common cause of universal administrative censoring." },
    { n: 2, title: "Swap the data call",
      body: "Replace data = gbsg with your dataset. Map your column names into outcome_var, event_var, treatment_var, continuous_vars, factor_vars." },
    { n: 3, title: "Set model = \"null\"",
      body: "For the extreme-subgroup evaluation, use subgroup_vars = NULL and model = \"null\" so every patient shares the same true HR." },
    { n: 4, title: "Redefine the subgroup list",
      body: "Edit the subgroups list object. Use the same structure (id / name / grp) with expressions tailored to your covariates. Keep the random15/20/40/60 benchmarks — they are the anchor." },
    { n: 5, title: "Calibrate censoring",
      body: "Run check_censoring_dgm() then calibrate_cens_adjust() so synthetic trials reproduce your study\u2019s censoring rate. Set analysis_time and max_entry to match your trial." },
    { n: 6, title: "Run, summarise, plot",
      body: "Use the same simulation loop: simulate_from_dgm() \u2192 coxph() \u2192 summarize_extreme_sims() \u2192 plot_ub_forest(). All tables and plots update automatically." }
  ];

  // 2 columns × 3 rows
  const y0 = 1.55;
  const rowH = 1.75;
  const colW = (SW - 2*MARGIN - 0.3) / 2;

  steps.forEach((st, i) => {
    const col = i % 2, row = Math.floor(i / 2);
    const x = MARGIN + col * (colW + 0.3);
    const y = y0 + row * rowH;

    s.addShape(pres.shapes.RECTANGLE, {
      x, y, w: colW, h: rowH - 0.15,
      fill: { color: C.pale }, line: { color: C.line, width: 0.5 }
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x, y, w: 0.08, h: rowH - 0.15,
      fill: { color: C.deep }, line: { color: C.deep, width: 0 }
    });

    // badge
    stepBadge(s, x + 0.25, y + 0.15, st.n, C.deep);

    s.addText(st.title, {
      x: x + 0.95, y: y + 0.15, w: colW - 1.1, h: 0.5,
      fontFace: "Calibri", fontSize: 15, bold: true, color: C.navy,
      valign: "middle", margin: 0
    });
    s.addText(st.body, {
      x: x + 0.25, y: y + 0.72, w: colW - 0.4, h: rowH - 0.9,
      fontFace: "Calibri", fontSize: 12, color: C.ink, margin: 0
    });
  });

  // Footer tip
  s.addText(
    "Tip · start by replacing only data and the column mappings; keep all tuning parameters at their defaults for the first run.",
    { x: MARGIN, y: SH - 0.7, w: SW - 2*MARGIN, h: 0.3,
      fontFace: "Calibri", fontSize: 11, italic: true, color: C.muted, margin: 0 }
  );
}

// ============================================================================
// SLIDE 14 — Conclusions
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.white };
  titleBand(s, "Conclusions", "Wrap-up");

  const concl = [
    { n: "01", title: "Small subgroups produce extreme results by chance",
      body: "Under a uniform HR = 0.70, subgroups with N \u2272 40 regularly yield UB(HR) > 3.0 — purely from sampling variability." },
    { n: "02", title: "Clinical rationale offers no statistical protection",
      body: "Clinical and random subgroups of the same size have statistically indistinguishable HR and UB distributions." },
    { n: "03", title: "The Cox model is correct",
      body: "Wide intervals in small subgroups are faithful reflections of genuine uncertainty, not a modelling failure." },
    { n: "04", title: "Simulation-based calibration is essential",
      body: "Without a null-DGM benchmark, \u201cextreme\u201d cannot be meaningfully quantified for a given subgroup size." },
    { n: "05", title: "The framework is portable",
      body: "Any trial can be swapped in for GBSG. The forestsearch package provides the full, reproducible pipeline." }
  ];

  // Stacked list of wide cards
  const y0 = 1.55;
  const rowH = 0.95;
  concl.forEach((c, i) => {
    const y = y0 + i * rowH;
    s.addShape(pres.shapes.RECTANGLE, {
      x: MARGIN, y, w: SW - 2*MARGIN, h: rowH - 0.12,
      fill: { color: C.pale }, line: { color: C.line, width: 0.5 }
    });
    s.addShape(pres.shapes.RECTANGLE, {
      x: MARGIN, y, w: 0.08, h: rowH - 0.12,
      fill: { color: i < 3 ? C.deep : C.gold }, line: { color: C.deep, width: 0 }
    });

    // number
    s.addText(c.n, {
      x: MARGIN + 0.25, y: y + 0.08, w: 1.0, h: rowH - 0.25,
      fontFace: "Georgia", fontSize: 28, bold: true, color: i < 3 ? C.deep : C.gold,
      valign: "middle", margin: 0
    });
    // title
    s.addText(c.title, {
      x: MARGIN + 1.35, y: y + 0.05, w: 4.4, h: rowH - 0.2,
      fontFace: "Calibri", fontSize: 15, bold: true, color: C.navy,
      valign: "middle", margin: 0
    });
    // body
    s.addText(c.body, {
      x: MARGIN + 5.85, y: y + 0.05,
      w: SW - 2*MARGIN - 5.85 - 0.2, h: rowH - 0.2,
      fontFace: "Calibri", fontSize: 12, color: C.ink,
      valign: "middle", margin: 0
    });
  });
}

// ============================================================================
// SLIDE 15 — Closing / resources
// ============================================================================
{
  const s = pres.addSlide();
  s.background = { color: C.navy };

  // accent vertical bar
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0, y: 0, w: 0.25, h: SH,
    fill: { color: C.gold }, line: { color: C.gold, width: 0 }
  });

  s.addText("Thank you", {
    x: 0.9, y: 1.1, w: SW - 1.5, h: 1.2,
    fontFace: "Georgia", fontSize: 54, bold: true, color: C.white, margin: 0
  });
  s.addShape(pres.shapes.LINE, {
    x: 0.9, y: 2.45, w: 1.2, h: 0, line: { color: C.gold, width: 3 }
  });

  s.addText([
    { text: "forestsearch", options: { bold: true, color: C.gold } },
    { text: " is an open-source R package providing a reproducible pipeline for simulation-based subgroup calibration.",
      options: { color: C.ice } }
  ], {
    x: 0.9, y: 2.7, w: SW - 1.8, h: 0.9,
    fontFace: "Calibri", fontSize: 17, margin: 0
  });

  // Resource cards row
  const ry = 4.1;
  const rh = 2.3;
  const rw = (SW - 1.8 - 0.5) / 2;

  // Package card
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.9, y: ry, w: rw, h: rh,
    fill: { color: C.deep }, line: { color: C.deep, width: 0 }
  });
  s.addText("Package", {
    x: 1.1, y: ry + 0.2, w: rw - 0.4, h: 0.35,
    fontFace: "Calibri", fontSize: 12, bold: true, color: C.gold,
    charSpacing: 4, margin: 0
  });
  s.addText([
    { text: 'pak::pak("larry-leon/forestsearch")',
      options: { color: C.white, fontFace: "Consolas", bold: true, breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "library(forestsearch)",
      options: { color: C.white, fontFace: "Consolas", bold: true } }
  ], {
    x: 1.1, y: ry + 0.7, w: rw - 0.4, h: rh - 0.9,
    fontFace: "Consolas", fontSize: 14, margin: 0
  });

  // Site / vignette card
  s.addShape(pres.shapes.RECTANGLE, {
    x: 0.9 + rw + 0.5, y: ry, w: rw, h: rh,
    fill: { color: C.teal }, line: { color: C.teal, width: 0 }
  });
  s.addText("Docs & vignette", {
    x: 0.9 + rw + 0.7, y: ry + 0.2, w: rw - 0.4, h: 0.35,
    fontFace: "Calibri", fontSize: 12, bold: true, color: C.gold,
    charSpacing: 4, margin: 0
  });
  s.addText([
    { text: "larry-leon.github.io/forestsearch",
      options: { color: C.white, fontFace: "Consolas", bold: true,
                 hyperlink: { url: "https://larry-leon.github.io/forestsearch/" },
                 breakLine: true } },
    { text: "", options: { breakLine: true } },
    { text: "…/articles/extreme_subgroups.html",
      options: { color: C.white, fontFace: "Consolas", bold: true,
                 hyperlink: { url: "https://larry-leon.github.io/forestsearch/articles/extreme_subgroups.html" } } }
  ], {
    x: 0.9 + rw + 0.7, y: ry + 0.7, w: rw - 0.4, h: rh - 0.9,
    fontFace: "Consolas", fontSize: 13, margin: 0
  });

  // Tiny author footer
  s.addText("Larry F. León  ·  April 2026", {
    x: 0.9, y: SH - 0.6, w: SW - 1.8, h: 0.3,
    fontFace: "Calibri", fontSize: 11, italic: true, color: C.ice,
    charSpacing: 3, margin: 0
  });
}

// ============================================================================
// Write output
// ============================================================================
pres.writeFile({ fileName: "/home/claude/deck/extreme_subgroups_deck.pptx" })
  .then(fileName => { console.log("WROTE: " + fileName); })
  .catch(err => { console.error("ERROR: " + err); process.exit(1); });
