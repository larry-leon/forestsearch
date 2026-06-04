# Session Directive — Read Before Anything Else

**Acknowledge these rules in your first response by listing them back in
one short sentence each. Do not skip this.**

---

## 1. Source of truth for files

For the entire duration of this chat:

- The current version of any file is whichever is **most recent in this
  conversation** — an artifact you produced, a file I pasted, or a file I
  uploaded mid-chat. Project knowledge (`/mnt/project/`) is a stale
  snapshot and is **not** the source of truth for any file we have
  already touched.
- Do **NOT** call `project_knowledge_search` to re-fetch a file that has
  already appeared as an artifact in this conversation, even partially.
- You may use `project_knowledge_search` only for files that have **not
  yet appeared** in this chat in any form.

## 2. Before producing any file revision

1. Scan the conversation history for the most recent version of the
   target file.
2. State explicitly which version you are building on, e.g.: *"Building
   on the artifact `forestsearch_main.R` you and I last revised in turn
   N, where the `merge(by.y = "id")` call was confirmed intentional."*
3. If you cannot determine which version is current — artifact vs.
   project knowledge vs. pasted snippet — **STOP and ask** before
   touching the file. Do not guess. Do not search.

## 3. Scope discipline

- Modify **only** what I explicitly ask you to modify.
- Do NOT opportunistically "also fix" unrelated code in the same file.
- Do NOT re-introduce code, comments, or formatting that was explicitly
  removed in a prior turn.
- Do NOT overwrite my edits silently. If your change conflicts with
  something I wrote, flag it and ask.

## 4. Ask for files instead of guessing

If answering, verifying, or revising something requires the **actual
contents** of a file you have not been shown in this chat — a sourced
params bundle (e.g. `_fsparams/*.R`), a helper, a config, a data
definition, a function body — **STOP and ask me to provide that file**
before proceeding.

- Do **NOT** infer, guess, or reason from the *assumed* contents of a
  file you have not seen. This applies even when you believe you can
  predict what the file contains.
- A correct-looking answer built on a guessed file is worse than a
  one-line request for the file: guessing file contents has repeatedly
  cost hours of rework and failed renders.
- If the file *has* already appeared in this chat, Rule 1 governs — use
  the most recent in-chat version; do not re-fetch or guess.
- "I need to see `X` before I can verify this" is always an acceptable
  and preferred response.

## 5. When in doubt

Ask. Always ask before searching project knowledge for a file we've
already worked on. The cost of one clarifying question is trivial; the
cost of reverting an hour of edits is not.

## 6. Correction phrase

If I say:

> **"Stale revert — re-read the most recent artifact for `<filename>`
> before continuing."**

Then immediately:
1. Locate the most recent artifact for that file in this conversation.
2. Re-read it in full.
3. Rebuild from that version.
4. Do **NOT** call `project_knowledge_search` for that file.

---

## Confirmation required

In your first response, list these six rules in one sentence each, then
proceed with my actual request. This confirmation is non-negotiable for
this session.
