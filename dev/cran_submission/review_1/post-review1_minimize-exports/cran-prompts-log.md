# CRAN Submission Prompts Log

> **What is this file?**
> A simple running log of prompts I've used with Claude to resolve CRAN submission
> issues for the ForestSearch package. Each entry records what the problem was,
> what I asked, and what worked. This way I don't have to re-discover solutions.

---

## How to Add a New Entry

Copy the blank template below, paste it at the top of the "Log Entries" section,
and fill in the fields. Delete any optional fields you don't need.

### Blank Template

```
---

### Entry: [short title]

- **Date:** YYYY-MM-DD
- **CRAN issue:** [What CRAN flagged or what R CMD check reported]
- **Package version:** [e.g., 0.2.1]
- **Prompt to Claude:**

> [Paste your prompt here. Keep the ">" at the start of each line
> so it shows up as a quote block.]

- **Key takeaway:** [One or two sentences on what fixed it]
- **Files changed:** [List any R scripts, DESCRIPTION, etc. that were modified]
- **Status:** ✅ Resolved / ⏳ In progress / ❌ Did not work
- **Notes:** [Optional — anything else worth remembering]
```

---

## Log Entries

*(Newest first — paste new entries right below this line)*

---

### Entry: Example — Global variable NOTE in R CMD check

- **Date:** 2026-03-18
- **CRAN issue:** `R CMD check` NOTE: "no visible binding for global variable 'x'"
  caused by `dplyr` column references inside `mutate()`, `filter()`, etc.
- **Package version:** 0.2.0
- **Prompt to Claude:**

> My R package ForestSearch gets this NOTE from R CMD check:
> "no visible binding for global variable 'trt'". It comes from
> dplyr code like `filter(data, trt == 1)`. What is the
> standard way to fix this for CRAN?

- **Key takeaway:** Add a line like `trt <- NULL` (or use
  `utils::globalVariables("trt")`) before or inside the function. The
  `.data$trt` pronoun from `rlang` is another option. CRAN prefers any
  of these over ignoring the NOTE.
- **Files changed:** `R/forestsearch.R`, `R/globals.R` (new file)
- **Status:** ✅ Resolved
- **Notes:** Created a dedicated `R/globals.R` with one
  `utils::globalVariables(c(...))` call listing all column names.

---

### Entry: Example — DESCRIPTION formatting

- **Date:** 2026-03-15
- **CRAN issue:** CRAN rejected submission because the Description field in
  DESCRIPTION used the package name without single quotes.
- **Package version:** 0.2.0
- **Prompt to Claude:**

> CRAN says my DESCRIPTION file needs changes. The email says:
> "Please always write package names, software names and API names
> in single quotes in Title and Description." How should I format
> the Description field?

- **Key takeaway:** Wrap package names in single quotes (e.g., 'ForestSearch',
  'survival') and spell out acronyms on first use. Don't start the Description
  with "This package" — CRAN prefers starting with a direct statement.
- **Files changed:** `DESCRIPTION`
- **Status:** ✅ Resolved
- **Notes:** Also had to remove a period from the Title field.

---

### Entry: Example — Rd examples running too long

- **Date:** 2026-03-10
- **CRAN issue:** `R CMD check` flagged examples taking > 5 seconds, which
  CRAN may reject.
- **Package version:** 0.1.9
- **Prompt to Claude:**

> Some of my roxygen examples in ForestSearch run bootstrap
> simulations and take over 10 seconds. CRAN check warns about
> this. What are my options to keep the examples but pass the
> check?

- **Key takeaway:** Wrap slow examples in `\donttest{}` (preferred over
  `\dontrun{}`). CRAN accepts `\donttest` for examples that are valid but
  slow. Use `\dontrun` only for examples that truly cannot be executed
  (e.g., need credentials). Alternatively, reduce `n_bootstrap` to a tiny
  value in the example.
- **Files changed:** `R/forestsearch_bootstrap_dofuture.R`
- **Status:** ✅ Resolved
- **Notes:** Kept a fast 2-iteration example as a runnable default,
  wrapped the full demo in `\donttest{}`.
  
  
  
  ```
---

### Entry: [short title]

- **Date:** YYYY-MM-DD
- **CRAN issue:** [What CRAN flagged or what R CMD check reported]
- **Package version:** [e.g., 0.2.1]
- **Prompt to Claude:**

> I am preparing an R package called forestsearch for CRAN submission.
> Please fetch the most recent forestsearch codebase https://github.com/larry-leon/forestsearch. 
>
> An extensive review has been conducted by Claude, however an independent reviewer has raised the following concerns:
> (A) I am concerned about the descriptive stubs. These are `@examples` sections which are comments-only. I suspect these will
> be flagged by the CRAN reviewer as insufficient. Do these need to be exported functions? For example,
> `extract_selected_tree_cuts()` seems like an internal-only function;
> (B)  If the example is short and runs quickly, there is no harm in keeping it.
> But it doesn't make much sense to include complex, long-running examples that
> have to be wrapped in \dontrun{} if an end user is unlikely to ever run the example code anways;
> (C)  In the branch cran-review-1, {forestsearch} currently has 118 exported functions in NAMESPACE and
> 255 manual pages in man/. However, only 42 out of the 118 are actually used in the vignettes/articles
> under vignettes/. This suggests that only those 42 functions are required to be used directly
> by an end user performing an analysis described in the vignettes. I recommend focusing your attention on writing thorough
> documentation for only important functions. All the other functions can be internal to the package and only receive a bare minimum of documentation;
> (D) Four examples in R/oc_analyses_gbsg.R still use \dontrun{}. Please conduct a comprehensive scan to ensure \dontrun{} is NOT used.
>
> Before proceeding please create a step-by-step strategy to address (A) - (D).  In particular, the main concern is to only export necessary functions 
> and for internal functions ensure that example are completely omitted.
> Please consult the vignettes and quarto documents in the KB to determine which functions should be internal, and for such internal functions ensure 
> examples can be completely eliminated.




- **Key takeaway:** [One or two sentences on what fixed it]
- **Files changed:** [List any R scripts, DESCRIPTION, etc. that were modified]
- **Status:** ✅ Resolved / ⏳ In progress / ❌ Did not work
- **Notes:** [Optional — anything else worth remembering]
```

  
  
