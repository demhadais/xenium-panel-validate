# Goals

## Standard Genes List Validation
- Every gene name + gene ID in target list is in the reference genome with the same name/ID combo
- No duplicates
- For simple genes, we need the following columns (snake_case, whitespace-stripped):
  - Gene name
  - Gene ID
  - Is backup
  - Group
  - Must have
- *For later*: A UI that allows a user to specify the minimum number of genes per group to try to retain
- Not in the disallowed targets-list
- Define a set of known aliases for each column that are also accepted, but when outputting, output the "canonical" column name

## Nonstandard Targets List Validation
TBD - we need to think about this more, but here are some ramblings
- For each target, we need the following three pieces of information (column names pending)
  - What do you want to detect?
  - Granularity of what you want to detect
  - What do you not want to detect
- maybe we just require different CSVs for different types of nonstandard targets

## Reference Dataset Validation
- Raw counts exist
- Gene IDs are present
- Warn about targets not in reference dataset
- Warn that the dataset may have been aligned against a different reference genome. Can detect this by finding gene IDs whose gene name don't match the mapping in the reference genome

## Report-generation
- Generate some kind of structured object that details what the user gave us, what we did, and what the output was. The overall paradigm of this report should be a list of errors or "findings" that may or may not be correctable. For the correctable case, we present the user with our best guess(es) of how the error could be fixed and ask them to choose. For example:
  ```json
  {
    "gene_list_errors": [
      {
        "submitted_gene_name": "gene",
        "submitted_gene_id": "ID",
        "error_type": "gene_name_id_mismatch"
      }
    ]
  }
  ```

## For the future: Probeset calculation
Eventually, do some calculations on the reference dataset so as to pre-calculate a "suggested" number of probesets for the panel designer

# Testing the tool
TBD, but we can definitely engineer some gene lists with failure modes and use real reference datasets with known failure modes?

# Outputs
- A report of everything to fix
  - In this report, a list of cell types and the number of cells per cell type
- A summary of the reference dataset
