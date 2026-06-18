# Implementation plan

Short-medium term, the goal is a small CLI + library that the SCBL (and
eventually a user) runs on a panel submission (gene list + reference dataset).
It reads a target list and a reference dataset, checks them against a pinned
genome, and hands back a list of things to fix. The point is to move routine
sanity checks to submission time so we stop spending several weeks back and
forth via email about symbol/ID mismatch and incorrectly formatted references.

Medium-medium term, we add separate validation and internal representation for 
custom targets (distinguishing isoforms, capturing common sub-gene sequences,
exogeneous targets, etc.).  More thought needed here.

Long-term, transition to a web UI and reuse this core rust lib as the backend
for gene list and reference dataset validation.

## What's already here

- `target_list.rs` reads a CSV target list and checks each target's Ensembl ID against
  the genome, plus that the name in the list matches the name the genome has for that ID.
  The `Target` struct already has `name`, `ensembl_id`, `group`, `is_backup`, and uses
  serde `alias` for a couple of column spellings.
- `gene_annotation.rs` reads a GTF with `noodles` into an id -> name map.
- `main.rs` exposes one `validate-targets` subcommand that writes errors out as JSON.

## Design decisions

Foundational assumptions:

- Pin to 2020-A. 
    - All IDs validate against the Panel Designer's reference
      (`refdata-gex-GRCh38-2020-A` for human, `refdata-gex-mm10-2020-A` for mouse,
      both Ensembl 98), not whatever release the user's analysis was on. 
    - An ID can be valid in a user's 2024-A analysis and not exist in the XPD universe
    - Gene symbols drift between releases, so a symbol can attached to the wrong ID 
    - IDs are versionless (`ENSG00000102755`, not `ENSG00000102755.6`).

- Permissive input, strict internal representation.
    - We cannot expect users to always provide us exactly what we ask for.
    - We will accept stale CSV/XLSX templates and normalize them to one
      canonical internal shape. 
    - Column aliases map to canonical names; we return the canonical name back
      (so hopefully the user gets the message). 
    - We never hand-fix a submission for the user; the tool tells them what's
      wrong and they fix it.

- Output is a report of "findings", and don't stop at the first issue.
    - The output is a list of findings, each with a stable name, code, severity,
      and location in the input. 
    - Some findings are deterministically correctable and we suggest the fix;
      some aren't. The report encompasses all findings.

- Use built-in name<->id mapping. 
  -  The gene_id -> gene_name mapping ships inside the binary.

- Read the `.h5ad` reference datasets in Rust.
  -  We can compile libhd5f within a static rust binary and can read an AnnData
     file directly.  No Python needed like we originally thought.

## Plan

Roughly in order. Each step should be usable on its own.

### 1. Finish the gene-list checks

Update `target_list.rs`

- Add `must_have` boolean to `Target`.
- Required columns: 
  ```{rust}
  gene_name 
  gene_id
  is_backup
  group 
  must_have 
  ```
- Keep a known set of aliases per column and accept those, but normalize to the
  canonical name. The alias map should be data (a table we edit), not a wall of serde attributes.
- Screen for two identical gene_name/gene_ids are an error. (Targeting the same
  gene in two genuinely distinct ways is a custom-target thing, which we're not
  doing yet, so for now same-gene-twice is just a duplicate.)
- Screen for disallowed-targets list, see below.

### 1a. The disallowed-targets list

10x publishes lists of genes that aren't covered by their probes and shouldn't go in a
custom design. We reject anything on the matching list. The lists live under "Genes not
available for design" in the Panel Designer docs:

- Human, Xenium v1:
  https://cdn.10xgenomics.com/raw/upload/v1702081344/software-support/Xenium-panels/human_2020-A-ref_noprobe_genesv2.csv
- Mouse, Xenium v1:
  https://cdn.10xgenomics.com/raw/upload/v1702081344/software-support/Xenium-panels/mouse_2020-A-ref_noprobe_genesv2.csv
- Human, Xenium Prime:
  https://cdn.10xgenomics.com/raw/upload/v1715034995/software-support/Xenium-panels/human-2020-A-ref-noprobe-genes-xenium-prime.csv
- Mouse, Xenium Prime:
  https://cdn.10xgenomics.com/raw/upload/v1715034995/software-support/Xenium-panels/mouse-2020-A-ref-noprobe-genes-xenium-prime.csv

We fetch these once and embed them the same way as the genome mapping. The right
list depends on species (human/mouse) and panel type (v1/Prime), so the submission
has to tell us which panel type it's for, or we check against the relevant pair.

**See open questions re: how to handle accepting species & chemistry type**

### 2. Built-in genome mapping

Right now the genome is a GTF passed in at runtime and parsed on every run.
Instead, we'll embed the gene_id -> gene_name mapping inside the binary so
there's nothing to download or point at.

- An offline step (not part of the shipped binary) parses the pinned 2020-A GTF with
  `noodles` into the id <-> name maps we need, serializes them, and writes them next to
  the crate to be embedded.  We could think about adding this later to the
  shipped binary for users who want to add their own reference, but hold for
  now.
- The binary loads the embedded mapping once at startup. At ~60k genes across 2
  species, this is small.
- One mapping per pinned reference (human, mouse), committed alongside the code.

This replaces the runtime GTF read in `gene_annotation.rs`. Keep the GTF-parsing code
around only as the offline build step that produces the embedded mapping.

### 3. Reference dataset checks (`.h5ad`)

Add a second mode that takes the reference dataset and checks it. From GOALS, it checks:

- Raw counts exist and are integer-valued.
- `var` has Ensembl IDs (and symbols).
- Warn about targets that aren't present in the reference dataset. Report the
  missing set.
- Warn that the dataset may be on a different genome release: find IDs whose
  gene name in the dataset disagrees with the name 2020-A has for that ID. That
  mismatch is the tell.
- Ensure there are valid observation annotations.

Read the file in Rust with included libhdf5. The user tells us which `.obs`
column is the annotation and where raw counts live (flags or a small config).
In the future, we can also add functionality to save the reference dataset in
10x's prefered format.

### 4. The report

Pull the findings into one structured object instead of the current flat error list.

- Each finding: code, severity, locus (file/row/field), message, optional suggested fix.
- The report also carries the summary GOALS asks for: cell types and cell counts
  per type from the reference dataset, plus a short reference-dataset summary.
- Render it two ways at least: human-readable to the terminal, and JSON for anything
  downstream.

### 5. XLSX input

We need to accept both CSV and XLSX input. CSV works today; add XLSX reading so
people can submit the spreadsheet template directly. A web UI comes later and is
out of scope here, but keeping the validation logic in the library (not the CLI)
is what makes that later UI cheap, so keep the split clean.

## Not doing yet

- Nonstandard / custom targets (isoforms, transgenes, microbial, cross-species). 
- Probeset-count suggestions from the reference dataset.
- Saving the reference dataset reformmated for 10x.
- The web UI.
- Any distribution/installer work.

## Current Open questions

- Panel type: do we make the submission declare v1 vs Prime, or check against both
  no-probe lists? Same question decides which disallowed list applies.
- To what level are we going to use LLMs here? 