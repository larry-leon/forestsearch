const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  Header, Footer, AlignmentType, LevelFormat,
  HeadingLevel, BorderStyle, WidthType, ShadingType,
  PageNumber, PageBreak
} = require("docx");

const FONT = "Arial";
const PAGE_W = 12240;
const PAGE_H = 15840;
const MARGIN = 1440;
const CONTENT_W = PAGE_W - 2 * MARGIN;

const HEADING_COLOR = "1F3864";
const ACCENT_COLOR = "2E75B6";
const TABLE_HEADER = "D5E8F0";
const TABLE_ALT = "F2F7FA";

const border = { style: BorderStyle.SINGLE, size: 1, color: "CCCCCC" };
const borders = { top: border, bottom: border, left: border, right: border };
const cellMargins = { top: 60, bottom: 60, left: 100, right: 100 };

function headerCell(text, width) {
  return new TableCell({
    borders,
    width: { size: width, type: WidthType.DXA },
    shading: { fill: TABLE_HEADER, type: ShadingType.CLEAR },
    margins: cellMargins,
    children: [new Paragraph({
      spacing: { before: 40, after: 40 },
      children: [new TextRun({ text, bold: true, font: FONT, size: 18 })]
    })]
  });
}

function dataCell(text, width, opts = {}) {
  return new TableCell({
    borders,
    width: { size: width, type: WidthType.DXA },
    shading: opts.shading ? { fill: opts.shading, type: ShadingType.CLEAR } : undefined,
    margins: cellMargins,
    verticalAlign: "top",
    children: [new Paragraph({
      spacing: { before: 40, after: 40 },
      children: [new TextRun({ text, font: FONT, size: 18, ...(opts.italic ? { italics: true } : {}), ...(opts.bold ? { bold: true } : {}) })]
    })]
  });
}

function p(children, opts = {}) {
  const runs = typeof children === "string"
    ? [new TextRun({ text: children, font: FONT, size: opts.size || 22, ...(opts.bold ? { bold: true } : {}), ...(opts.italic ? { italics: true } : {}) })]
    : children;
  return new Paragraph({
    spacing: { before: opts.before || 120, after: opts.after || 120, line: opts.line || 276 },
    alignment: opts.align || AlignmentType.LEFT,
    children: runs,
    ...(opts.heading ? { heading: opts.heading } : {}),
  });
}

function heading(text, level, opts = {}) {
  const hl = level === 1 ? HeadingLevel.HEADING_1 : level === 2 ? HeadingLevel.HEADING_2 : HeadingLevel.HEADING_3;
  const sz = level === 1 ? 32 : level === 2 ? 26 : 22;
  return new Paragraph({
    heading: hl,
    spacing: { before: level === 1 ? 360 : 240, after: level === 1 ? 200 : 140 },
    ...(opts.pageBreakBefore ? { pageBreakBefore: true } : {}),
    children: [new TextRun({ text, bold: true, font: FONT, size: sz, color: HEADING_COLOR })]
  });
}

function bullet(children, level = 0) {
  const runs = typeof children === "string"
    ? [new TextRun({ text: children, font: FONT, size: 22 })]
    : children;
  return new Paragraph({
    numbering: { reference: "bullets", level },
    spacing: { before: 60, after: 60, line: 276 },
    children: runs
  });
}

function bulletBold(boldText, normalText, level = 0) {
  return new Paragraph({
    numbering: { reference: "bullets", level },
    spacing: { before: 60, after: 60, line: 276 },
    children: [
      new TextRun({ text: boldText, bold: true, font: FONT, size: 22 }),
      new TextRun({ text: normalText, font: FONT, size: 22 }),
    ]
  });
}

function subBullet(children) {
  return bullet(children, 1);
}

function subBulletBold(boldText, normalText) {
  return bulletBold(boldText, normalText, 1);
}

const doc = new Document({
  styles: {
    default: { document: { run: { font: FONT, size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 32, bold: true, font: FONT, color: HEADING_COLOR },
        paragraph: { spacing: { before: 360, after: 200 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 26, bold: true, font: FONT, color: HEADING_COLOR },
        paragraph: { spacing: { before: 240, after: 140 }, outlineLevel: 1 } },
      { id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 22, bold: true, font: FONT, color: HEADING_COLOR },
        paragraph: { spacing: { before: 200, after: 120 }, outlineLevel: 2 } },
    ]
  },
  numbering: {
    config: [
      { reference: "bullets",
        levels: [
          { level: 0, format: LevelFormat.BULLET, text: "\u2022", alignment: AlignmentType.LEFT,
            style: { paragraph: { indent: { left: 720, hanging: 360 } } } },
          { level: 1, format: LevelFormat.BULLET, text: "\u2013", alignment: AlignmentType.LEFT,
            style: { paragraph: { indent: { left: 1080, hanging: 360 } } } },
          { level: 2, format: LevelFormat.BULLET, text: "\u25E6", alignment: AlignmentType.LEFT,
            style: { paragraph: { indent: { left: 1440, hanging: 360 } } } },
        ] },
    ]
  },
  sections: [{
    properties: {
      page: {
        size: { width: PAGE_W, height: PAGE_H },
        margin: { top: MARGIN, right: MARGIN, bottom: MARGIN, left: MARGIN }
      }
    },
    headers: {
      default: new Header({
        children: [new Paragraph({
          alignment: AlignmentType.RIGHT,
          children: [new TextRun({ text: "Guidance on Subgroup Analyses in Clinical Trials \u2014 Appendix: Scenario D", font: FONT, size: 16, italics: true, color: "888888" })],
          border: { bottom: { style: BorderStyle.SINGLE, size: 4, color: ACCENT_COLOR, space: 4 } }
        })]
      })
    },
    footers: {
      default: new Footer({
        children: [new Paragraph({
          alignment: AlignmentType.CENTER,
          border: { top: { style: BorderStyle.SINGLE, size: 2, color: "CCCCCC", space: 4 } },
          children: [
            new TextRun({ text: "Page ", font: FONT, size: 16, color: "888888" }),
            new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 16, color: "888888" }),
          ]
        })]
      })
    },
    children: [

      // ===== TITLE =====
      heading("Scenario D: MRCT Regional Consistency via Subgroup Identification", 1),

      p([
        new TextRun({ text: "Technical Report Summary: ", font: FONT, size: 22, bold: true }),
        new TextRun({ text: "MRCT Consistency Evaluation via Subgroup Identification \u2014 A ForestSearch Framework with a View Towards Biomarker Causal Effects", font: FONT, size: 22, italics: true }),
      ]),
      p([
        new TextRun({ text: "Authors: ", font: FONT, size: 20, bold: true }),
        new TextRun({ text: "Larry Leon, Masashi Shimura, Xiaojing Cui, Kenichi Takahashi, Shuping Jiang, and William Wang", font: FONT, size: 20 }),
      ], { before: 40, after: 40 }),
      p([
        new TextRun({ text: "R package: ", font: FONT, size: 20, bold: true }),
        new TextRun({ text: "forestsearch (github.com/larry-leon/forestsearch)", font: FONT, size: 20 }),
      ], { before: 40, after: 160 }),


      // ===== PROBLEM STATEMENT =====
      heading("Motivation and Problem Statement", 2),

      bulletBold("The MRCT regional consistency problem: ",
        "A global oncology trial achieves statistical significance for the overall (ITT) population, but a regional subpopulation \u2014 such as Asia-Pacific (AP), comprising roughly 15\u201320% of the trial \u2014 shows a hazard ratio near 1.0, failing to demonstrate benefit."),
      bulletBold("Regulatory requirement: ",
        "ICH E17 guidelines require that regional subpopulations show evidence of treatment consistency with the overall result."),
      bulletBold("Central question: ",
        "Can a biomarker-defined subgroup be identified in the non-AP (rest-of-world) training data such that, when applied to the AP testing data, it recovers regional consistency?"),
      bulletBold("Approach: ",
        "The ForestSearch algorithm (L\u00E9on et al., 2024, Statistics in Medicine) is applied within a Weibull AFT/Cox dual modeling framework, with rigorous simulation-based evaluation of operating characteristics."),


      // ===== STATISTICAL FRAMEWORK =====
      heading("Statistical Framework", 2),

      p([
        new TextRun({ text: "Weibull AFT / Cox Dual Representation", font: FONT, size: 22, bold: true, color: ACCENT_COLOR }),
      ], { before: 160, after: 80 }),

      bulletBold("Core modeling strategy: ",
        "Exploits the mathematical equivalence between the Weibull Accelerated Failure Time (AFT) model and the Cox proportional hazards model under the Weibull assumption."),
      subBulletBold("AFT model: ",
        "log(T) = log(\u03B8) + L\u2032\u03B3 + \u03C4\u03B5, where \u03C4 = 1/\u03BD is the scale parameter and \u03B5 follows an extreme-value distribution."),
      subBulletBold("Cox model: ",
        "\u03BB(t; L) = \u03BB\u2080(t) exp(L\u2032\u03B2), with Weibull baseline hazard."),
      subBulletBold("Bridging identity: ",
        "\u03B2 = \u2212\u03B3/\u03C4 converts AFT coefficients to hazard-ratio coefficients."),
      bulletBold("Dual interpretation: ",
        "Potential outcomes are simulated on the AFT scale (clean causal interpretation as geometric mean ratios of survival times), while treatment effects are interpreted on the hazard-ratio scale (standard clinical reporting currency)."),

      p([
        new TextRun({ text: "Two-Phase Spline Model for Biomarker Heterogeneity", font: FONT, size: 22, bold: true, color: ACCENT_COLOR }),
      ], { before: 160, after: 80 }),

      bulletBold("Piecewise-linear Cox linear predictor: ",
        "Single knot at biomarker value k, with five parameters capturing the main treatment effect, prognostic biomarker effect, and biomarker \u00D7 treatment interactions above and below the knot."),
      bulletBold("Working example: ",
        "With knot k = 5 and upper anchor \u03B6 = 10, parameterized so that HR(z=0) = 3.00 (strong harm), HR(z=5) = 1.25 (marginal harm), and HR(z=10) = 0.50 (strong benefit)."),
      bulletBold("Clinical interpretation: ",
        "Creates a realistic scenario where treatment harms patients with low biomarker values and benefits those with high values, with the treatment effect continuously varying across the biomarker spectrum."),

      p([
        new TextRun({ text: "Three Complementary Estimands", font: FONT, size: 22, bold: true, color: ACCENT_COLOR }),
      ], { before: 160, after: 80 }),

      // Estimands table
      new Table({
        width: { size: CONTENT_W, type: WidthType.DXA },
        columnWidths: [2200, 3080, 4080],
        rows: [
          new TableRow({ children: [
            headerCell("Estimand", 2200),
            headerCell("Definition", 3080),
            headerCell("Interpretation", 4080),
          ]}),
          new TableRow({ children: [
            dataCell("Average Hazard Ratio (AHR)", 2200, { bold: true }),
            dataCell("Exponentiated average of individual causal log-HRs over subgroup S", 3080),
            dataCell("Clean causal average. Inherits dual causal status under Weibull AFT: simultaneously individual-level and population-level.", 4080),
          ]}),
          new TableRow({ children: [
            dataCell("Controlled Direct Effect (CDE)", 2200, { bold: true, shading: TABLE_ALT }),
            dataCell("Ratio of average exponentiated hazards (treated vs. control) within subgroup", 3080, { shading: TABLE_ALT }),
            dataCell("Accounts for prognostic covariate distribution within subgroup. Robust to within-subgroup confounding.", 4080, { shading: TABLE_ALT }),
          ]}),
          new TableRow({ children: [
            dataCell("Marginal HR", 2200, { bold: true }),
            dataCell("Cox model coefficient on stacked potential outcomes", 3080),
            dataCell("Clinically familiar Cox coefficient. Standard reporting currency for regulatory communication.", 4080),
          ]}),
        ]
      }),

      p("These three estimands complement one another: AHR gives a clean causal average, CDE accounts for covariate confounding within the subgroup, and the Marginal HR provides the clinically familiar Cox coefficient.", { before: 80, after: 80 }),

      new Paragraph({ children: [new PageBreak()] }),

      // ===== FORESTSEARCH ALGORITHM =====
      heading("The ForestSearch Algorithm", 2),

      p("ForestSearch (L\u00E9on et al., 2024) is a three-stage exploratory subgroup identification algorithm:"),

      p([
        new TextRun({ text: "Stage 1 \u2014 Candidate Feature Selection", font: FONT, size: 22, bold: true, color: ACCENT_COLOR }),
      ], { before: 120, after: 60 }),

      bulletBold("Generalized Random Forests (GRF): ",
        "Causal survival forests targeting RMST differences to detect treatment effect heterogeneity. GRF suggests data-driven cut-points for continuous covariates."),
      bulletBold("LASSO regularization: ",
        "Cox penalized regression performs variable screening, mitigating false discovery from noise variables."),
      bulletBold("Forced cut-points: ",
        "Domain-guided biomarker cuts ensure clinically meaningful thresholds are always evaluated (e.g., biomarker \u2264 0, 1, 2, 5)."),
      bullet("Continuous covariates are converted into binary indicator splits (medians, quartiles, GRF-derived cuts)."),

      p([
        new TextRun({ text: "Stage 2 \u2014 Exhaustive Combinatorial Search", font: FONT, size: 22, bold: true, color: ACCENT_COLOR }),
      ], { before: 120, after: 60 }),

      bullet("All single and pairwise combinations of L candidate indicators are enumerated (up to L(L\u22121)/2 + L candidates)."),
      bullet("Candidates filtered by minimum sample size (n \u2265 60) and minimum events per arm (d \u2265 12)."),
      bullet("Candidates screened for signal of harm: Cox HR \u2265 threshold (e.g., 0.90)."),

      p([
        new TextRun({ text: "Stage 3 \u2014 Splitting Consistency Validation", font: FONT, size: 22, bold: true, color: ACCENT_COLOR }),
      ], { before: 120, after: 60 }),

      bullet("Each surviving candidate undergoes R = 1,000 random 50/50 splits."),
      bulletBold("Consistency criterion: ",
        "A candidate is \u201Cconsistent with harm\u201D if both halves yield HR \u2265 1.0 in at least 80\u201390% of splits."),
      bulletBold("Selection: ",
        "The algorithm selects the subgroup \u201Cmaximally consistent with harm\u201D per the chosen focus (sg_focus)."),

      p([
        new TextRun({ text: "Three Subgroup Identification Objectives", font: FONT, size: 22, bold: true, color: ACCENT_COLOR }),
      ], { before: 160, after: 80 }),

      new Table({
        width: { size: CONTENT_W, type: WidthType.DXA },
        columnWidths: [1200, 2400, 2400, 3360],
        rows: [
          new TableRow({ children: [
            headerCell("Objective", 1200),
            headerCell("Goal", 2400),
            headerCell("HR Criterion", 2400),
            headerCell("Strategy", 3360),
          ]}),
          new TableRow({ children: [
            dataCell("A", 1200, { bold: true }),
            dataCell("Largest detrimental subgroup", 2400),
            dataCell("HR \u2265 0.90", 2400),
            dataCell("Identify harm subgroup; complement = benefit population", 3360),
          ]}),
          new TableRow({ children: [
            dataCell("B (Primary)", 1200, { bold: true, shading: TABLE_ALT }),
            dataCell("Smallest harmful subgroup", 2400, { shading: TABLE_ALT }),
            dataCell("HR \u2265 0.90", 2400, { shading: TABLE_ALT }),
            dataCell("Minimize harm subgroup so that complement (\u201Cgood\u201D subgroup) is maximally large with enhanced benefit", 3360, { shading: TABLE_ALT }),
          ]}),
          new TableRow({ children: [
            dataCell("C", 1200, { bold: true }),
            dataCell("Largest beneficial subgroup", 2400),
            dataCell("HR \u2264 0.60", 2400),
            dataCell("Directly target strong benefit", 3360),
          ]}),
        ]
      }),

      p([
        new TextRun({ text: "In this application, Objective B with sg_focus = \u201CminSG\u201D is the primary strategy: ", font: FONT, size: 22 }),
        new TextRun({ text: "find the smallest group with attenuated/harmful treatment effect so that the complementary \u201Cgood\u201D subgroup is as large as possible while still showing enhanced benefit.", font: FONT, size: 22, italics: true }),
      ], { before: 80 }),


      // ===== TRAIN/TEST PARADIGM =====
      heading("Train/Test Paradigm: Regional Separation", 2),

      bulletBold("Training set (Non-AP): ",
        "~80\u201385% of the trial. ForestSearch runs exclusively on this population to identify subgroups."),
      bulletBold("Testing set (AP): ",
        "~15\u201320% of the trial. The identified subgroup definition is applied to AP to evaluate whether regional consistency is recovered."),
      bulletBold("Overfitting prevention: ",
        "The subgroup is never \u201Ctuned\u201D to the AP data. The AP region is given a substantially different baseline prognosis from the rest of the world, making the consistency challenge realistic."),
      bulletBold("Rationale: ",
        "This mirrors the regulatory use case where a globally identified biomarker signature must independently demonstrate consistency in a smaller regional population."),

      new Paragraph({ children: [new PageBreak()] }),

      // ===== SIMULATION DESIGN =====
      heading("Simulation Design and Data-Generating Mechanism", 2),

      p([
        new TextRun({ text: "Case-Study Dataset", font: FONT, size: 22, bold: true, color: ACCENT_COLOR }),
      ], { before: 120, after: 60 }),

      bullet("Synthetic but realistic dataset (N \u2248 2,000) mimicking a global oncology trial."),
      bullet("Preserves realistic covariate distributions, censoring patterns, and regional composition."),
      bullet("The data-generating mechanism (DGM) is built by fitting a Weibull AFT model to this dataset via generate_aft_dgm_flex(), incorporating the two-phase spline for biomarker-driven heterogeneous treatment effects and a separate censoring model selected by AIC/BIC."),

      p([
        new TextRun({ text: "Two Simulation Scenarios", font: FONT, size: 22, bold: true, color: ACCENT_COLOR }),
      ], { before: 160, after: 80 }),

      new Table({
        width: { size: CONTENT_W, type: WidthType.DXA },
        columnWidths: [2000, 3380, 3980],
        rows: [
          new TableRow({ children: [
            headerCell("Scenario", 2000),
            headerCell("Spline Log-HRs", 3380),
            headerCell("Description", 3980),
          ]}),
          new TableRow({ children: [
            dataCell("HTE (Alternative)", 2000, { bold: true }),
            dataCell("log(3.00, 1.25, 0.50)", 3380),
            dataCell("Strong biomarker-driven heterogeneity: harm at low biomarker, benefit at high biomarker.", 3980),
          ]}),
          new TableRow({ children: [
            dataCell("Null (Uniform Benefit)", 2000, { bold: true, shading: TABLE_ALT }),
            dataCell("log(0.70, 0.70, 0.70)", 3380, { shading: TABLE_ALT }),
            dataCell("Flat HR = 0.70 across all biomarker levels. No heterogeneity; treatment benefits everyone equally.", 3980, { shading: TABLE_ALT }),
          ]}),
        ]
      }),

      p([
        new TextRun({ text: "Per-Replicate Simulation Workflow (N = 1,000\u20135,000 replicates)", font: FONT, size: 22, bold: true, color: ACCENT_COLOR }),
      ], { before: 160, after: 80 }),

      bulletBold("Step 1 \u2014 Simulate: ",
        "Draw n = 500 patients from the DGM; generate survival times on the AFT scale; apply administrative censoring at analysis_time = 60 months."),
      bulletBold("Step 2 \u2014 Split: ",
        "Divide into Non-AP (training, ~80\u201385%) and AP (testing, ~15\u201320%)."),
      bulletBold("Step 3 \u2014 ForestSearch: ",
        "Run the full GRF \u2192 LASSO \u2192 enumeration \u2192 screening \u2192 consistency pipeline on Non-AP training data."),
      bulletBold("Step 4 \u2014 Decision: ",
        "If a subgroup is found, apply its definition to AP and compute HR(AP, subgroup). If not found, record any_found = 0."),
      bulletBold("Step 5 \u2014 Record: ",
        "Per-replicate results including subgroup definition, HR(non-AP), HR(AP), HR(ITT), and HR(AP, subgroup)."),

      p([
        new TextRun({ text: "Key simulation parameters: ", font: FONT, size: 22, bold: true }),
        new TextRun({ text: "hr.threshold = 0.90, hr.consistency = 0.80, pconsistency.threshold = 0.80, use_grf = TRUE, use_lasso = TRUE, minimum sizes n \u2265 60 and events d \u2265 12, forced biomarker cuts at z_bm \u2264 0, 1, 2, 5.", font: FONT, size: 22 }),
      ]),


      // ===== OPERATING CHARACTERISTICS =====
      heading("Operating Characteristics Evaluated", 2),

      bulletBold("Subgroup discovery rate: ",
        "How often does ForestSearch identify a qualifying subgroup? (any_found across replicates)"),
      bulletBold("HR distributions: ",
        "Density/histogram plots of HR(Non-AP, ITT), HR(Overall, ITT), HR(AP, ITT), and HR(AP, subgroup) across replicates."),
      bulletBold("Subgroup composition: ",
        "Which subgroup definitions are most frequently selected (visualized by plot_sg_distribution)."),
      bulletBold("Summary metrics: ",
        "Combined HTE and Null scenario metrics including false-discovery rates, estimation bias, and CI coverage."),

      p([
        new TextRun({ text: "Expected results:", font: FONT, size: 22, bold: true }),
      ], { before: 120, after: 40 }),
      bulletBold("Under HTE scenario: ",
        "ForestSearch should successfully identify biomarker-high subgroups where AP consistency is recovered."),
      bulletBold("Under Null scenario: ",
        "ForestSearch should have a low false-discovery rate (i.e., rarely claim a subgroup exists when treatment benefit is truly uniform)."),


      // ===== CROSS-VALIDATION & INFERENCE =====
      heading("Cross-Validation and Bootstrap Inference", 2),

      bulletBold("Cross-validation: ",
        "Leave-one-out and 10-fold CV with agreement metrics (sensitivity, PPV, subgroup match) between CV-test and full-data results."),
      bulletBold("Parallelized bootstrap: ",
        "B = 300\u20132,000 replicates via doFuture, where the full ForestSearch algorithm is re-run within each replicate."),
      bulletBold("De-biased estimation: ",
        "Bootstrap bias correction with infinitesimal jackknife (IJ) variance estimation for valid confidence intervals."),


      new Paragraph({ children: [new PageBreak()] }),

      // ===== KEY TAKEAWAYS =====
      heading("Key Findings and Takeaways", 2),

      bulletBold("ForestSearch provides a principled approach ",
        "to subgroup identification that directly addresses MRCT regional consistency challenges."),
      bulletBold("The Weibull AFT/Cox dual framework enables causal interpretation: ",
        "simulated on the AFT scale, interpreted via hazard ratios. The scale-change parameter simultaneously captures individual-level and population-level causal effects."),
      bulletBold("Three complementary estimands (AHR, CDE, Marginal HR) ",
        "robustly characterize treatment effect heterogeneity from different perspectives, strengthening the evidence basis for subgroup claims."),
      bulletBold("The train/test paradigm (non-AP \u2192 AP) avoids overfitting: ",
        "subgroup definitions identified in the rest-of-world are independently validated in the regional population of interest."),
      bulletBold("Data-driven simulation evaluates full operating characteristics: ",
        "false-discovery rates, estimation bias, and CI coverage, all grounded in realistic trial data via calibrated data-generating mechanisms."),
      bulletBold("Simulation results demonstrate ",
        "the approach can recover AP consistency under heterogeneous effects while controlling false-discovery under the null."),


      // ===== REGULATORY CONTEXT =====
      heading("Relationship to Other Technical Reports", 2),

      p("This work directly complements the Extreme Subgroups simulation study (Scenarios A/B/E), which demonstrates that standard subgroup analyses can produce extreme-looking hazard ratios purely from sampling variability. Together, the two bodies of work address different sides of the MRCT problem:"),

      bulletBold("Extreme Subgroups study (Scenarios A/B/E): ",
        "Establishes calibration benchmarks showing how extreme small-subgroup results can appear under a uniform treatment effect, providing null-hypothesis context for interpreting subgroup findings."),
      bulletBold("MRCT Consistency study (Scenario D, this report): ",
        "Proposes a constructive solution \u2014 when regional consistency appears to fail, ForestSearch can identify biomarker-defined subgroups that recover consistency, with rigorous simulation-based validation."),

      p([
        new TextRun({ text: "Key regulatory applications:", font: FONT, size: 22, bold: true }),
      ], { before: 120, after: 40 }),
      bullet("Proactive risk quantification for sponsors preparing for regulatory review of MRCT regional consistency."),
      bullet("Calibration tools for regulatory reviewers assessing whether apparent regional inconsistency reflects genuine heterogeneity or sampling variability."),
      bullet("HTA subgroup evaluations where small country-level populations face elevated risk of extreme point estimates."),

      p([
        new TextRun({ text: "The overarching message is that principled, simulation-calibrated subgroup identification can transform the MRCT consistency problem from an apparent failure into a scientifically grounded refinement of the treatment population, while maintaining rigorous control of false discovery.", font: FONT, size: 22, italics: true }),
      ], { before: 160, after: 120 }),


      // ===== KEY R FUNCTIONS =====
      heading("Key R Packages and Functions", 2),

      new Table({
        width: { size: CONTENT_W, type: WidthType.DXA },
        columnWidths: [3200, 6160],
        rows: [
          new TableRow({ children: [
            headerCell("Package / Function", 3200),
            headerCell("Role", 6160),
          ]}),
          new TableRow({ children: [
            dataCell("forestsearch", 3200, { bold: true }),
            dataCell("Core subgroup identification algorithm (L\u00E9on et al., 2024)", 6160),
          ]}),
          new TableRow({ children: [
            dataCell("weightedsurv", 3200, { bold: true, shading: TABLE_ALT }),
            dataCell("Kaplan-Meier plotting with weighted survival", 6160, { shading: TABLE_ALT }),
          ]}),
          new TableRow({ children: [
            dataCell("grf", 3200, { bold: true }),
            dataCell("Generalized Random Forests \u2014 causal survival forests for HTE detection", 6160),
          ]}),
          new TableRow({ children: [
            dataCell("survival", 3200, { bold: true, shading: TABLE_ALT }),
            dataCell("survreg() for Weibull AFT, coxph() for Cox PH", 6160, { shading: TABLE_ALT }),
          ]}),
          new TableRow({ children: [
            dataCell("generate_aft_dgm_flex()", 3200, { bold: true }),
            dataCell("Builds the Weibull AFT DGM with spline biomarker effects and censoring model", 6160),
          ]}),
          new TableRow({ children: [
            dataCell("simulate_from_dgm()", 3200, { bold: true, shading: TABLE_ALT }),
            dataCell("Draws simulated trial data from the DGM", 6160, { shading: TABLE_ALT }),
          ]}),
          new TableRow({ children: [
            dataCell("mrct_region_sims()", 3200, { bold: true }),
            dataCell("Orchestrates the full simulation study across replicates (parallelized via doFuture)", 6160),
          ]}),
          new TableRow({ children: [
            dataCell("cox_cs_fit()", 3200, { bold: true, shading: TABLE_ALT }),
            dataCell("Fits the two-phase spline Cox model and plots biomarker heterogeneity", 6160, { shading: TABLE_ALT }),
          ]}),
          new TableRow({ children: [
            dataCell("cox_ahr_cde_analysis()", 3200, { bold: true }),
            dataCell("Computes AHR and CDE estimands across biomarker thresholds", 6160),
          ]}),
          new TableRow({ children: [
            dataCell("SGplot_estimates()", 3200, { bold: true, shading: TABLE_ALT }),
            dataCell("Plots HR distributions across replicates and populations", 6160, { shading: TABLE_ALT }),
          ]}),
          new TableRow({ children: [
            dataCell("summaryout_mrct()", 3200, { bold: true }),
            dataCell("Produces combined summary tables for HTE and Null scenario operating characteristics", 6160),
          ]}),
        ]
      }),

    ]
  }]
});

Packer.toBuffer(doc).then(buffer => {
  fs.writeFileSync("/home/claude/mrct_appendix_summary.docx", buffer);
  console.log("MRCT appendix summary created successfully");
});
