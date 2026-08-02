
## Documentation Standardss

1. **Root-level main markdown**: Any functional, architectural, or code updates MUST update relevant subdirectory documentation after work completion.

2. **Every folder MUST have a README.md with**:
   - Concise architecture description (≤3 lines)
   - List of each file with: name, status, and function
   - Header declaration: "Once my folder has changes, please update me."

3. **Every file MUST have header comments**:
   - **Input**: External dependencies (what this file needs)
   - **Output**: What it provides externally (what others use from this)
   - **Pos**: Role/position in the local system architecture
   - Reminder: "Once I am updated, update my header comments and folder's md."


## Lab Notebook Protocol
1. **`docs/metalabbook.md`** — single index table for ALL studies. Columns:
   `Slug | Status (active/done/abandoned) | Started | Last Updated | One-line Summary`.
   Rows sorted by `Last Updated` descending.

2. **`docs/labbooks/<NN>_<slug>.md`** — one file per study. `<NN>` = two-digit study number
   (next unused number in `docs/labbooks/`); the study's report reuses the same `<NN>_<slug>` stem.

   Top of file:
   ```
   Hypothesis: <one falsifiable sentence>
   Status: active          (active | done | abandoned — keep in sync with the metalabbook row)
   Started: YYYY-MM-DD
   ```

   Below the header: one `## YYYY-MM-DD` entry per day, **newest first** — prepend new entries
   directly under the header, never append at the bottom. If today already has an entry, extend
   it rather than adding a second heading for the same date.

   Entry template (use these exact bold labels):
   ```
   ## YYYY-MM-DD
   **What I did** — actions, commands/configs run, code touched, with paths.
   **What I observed** — concrete numbers, output/plot paths, error messages. Facts only; no interpretation here.
   **What I think it means** — interpretation. Confidence: low / medium / high. One line on what evidence would change this reading.
   **Next** — (optional) the immediate follow-up this implies.
   ```

3. **`docs/features/<NN>_<slug>.md`** — one file per substantial new feature (new behaviour or
   algorithm — not bugfixes or refactors), written when the feature lands. `<NN>` = two-digit
   feature number, the next unused one in `docs/features/` (numbered independently of the labbooks).

   Template (use these exact headings):
   ```
   # <Feature name>
   Added: YYYY-MM-DD
   Code: <main files/modules involved>

   ## What it does
   One paragraph: the problem it solves and the user-visible behaviour.

   ## Design choices
   One bullet per significant choice: what was chosen, why, and what was rejected.

   ## Algorithm details
   How it works — the steps, key equations/parameters/thresholds, and where each lives in the code.

   ## Caveats
   Known limitations, assumptions that can break, inputs it handles badly.
   ```

### When I report progress on a study
- Locate the matching `docs/labbooks/<NN>_<slug>.md`.
- Prepend a new `## YYYY-MM-DD` entry (do NOT append at the bottom — newest entries go on top, under the hypothesis/status header).
- Update that study's `Last Updated` in `docs/metalabbook.md` and re-sort the table so it stays in `Last Updated` descending order.
- If `Status` should change (e.g. `active` → `done`/`abandoned`), update both the labbook file header and the metalabbook row.

### When the work looks like a NEW study
- Do NOT silently create a new slug or new labbook file. Ask me first: propose a slug, the hypothesis, and confirm before creating `docs/labbooks/<slug>.md` and adding the row to `docs/metalabbook.md`.

### When finishing a study
- Write a detailed report including the intro, method, results and discussion in `docs/report/<NN>_<slug>.md`

### When finishing adding a new feature
- Write (or update) `docs/features/<NN>_<slug>.md` following the template in item 3 above.
- Add a one-line mention of the feature to `README.md`.


# General 
- use agents, and different agents for different levels of tasks.
- remove the short investigation files after using them.
- never write implementation/diagnosis reports .md until the finalisation or with my clear instructions, just reply in the chat briefly.
- when I say discuss interactively on something (like research plan), use more AskUserQuestion.
- Try to use straightforward and clear language when I'm asking you to explain anything or when you are reporting things to me. Make some efforts to make sure Your output is reasonable in terms of the language use. 
- Always update labbook when you have a long task and you are waiting for the queue or the task to finish.
- Refer to labbook for information first.
- refer to 'docs/metalabbook.md' for the project plan and context.
- ALWAYS update the read me, docs/metalabbook.md and docs/labbook (or the dedicated version for that task) after any big progress made. in plan just write a small sentence to say what is done and what is found. in labbook, include more details such as numerical results and folders/directories involved.
- use uv to manage the environment



