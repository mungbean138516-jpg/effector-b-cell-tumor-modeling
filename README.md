# GitHub-Ready MacLean Lab Package

This folder is a curated subset of the larger `MACLEAN_LAB` workspace, assembled specifically for GitHub upload and portfolio-style presentation.

It keeps the strongest representative materials from the project:

- core Julia modeling scripts
- post-meeting refinement scripts, including the updated `step1_d2`, `step2_initialization_d2`, and `step3_heatmap` work
- representative output figures and CSV summaries
- presentation builder scripts
- a small set of summary documents and PDFs

## Why this folder exists

The full project workspace contains many intermediate outputs, repeated drafts, zip archives, and presentation artifacts. This package trims that down to a cleaner set of materials that still shows the scope of the work:

- model development
- parameter sweeps
- stochastic and deterministic analyses
- post-meeting refinement across dense sweeps, initialization robustness, and the refined `b5 x b6` heatmap analysis
- presentation-ready summaries

## Folder structure

- `code/core/`: main Julia analysis scripts covering CTL comparison, beta sweeps, parameter sweeps, and SDE refinement.
- `code/post_meeting/step1_d2/`: updated dense-sweep refinement script and associated Step 1 follow-up materials.
- `code/post_meeting/step2_initialization_d2/`: updated initialization-robustness script for the `d2` round.
- `code/post_meeting/step3_heatmap/`: refined `b5 x b6` phase-diagram script plus the follow-up plotting-fix script.
- `code/presentation_builders/`: Python scripts used to generate or assemble presentation materials.
- `results/ctl_compare/`: representative outputs from the direct-vs-coupled CTL comparison.
- `results/beta5/`: representative beta-5 deterministic and stochastic outputs.
- `results/beta6/`: representative refined beta-6 outputs and heatmap summaries.
- `results/parameter_sweeps/`: representative wide-sweep figures and CSV outputs.
- `results/post_meeting/step1_d2/`: updated Step 1 dense-sweep figures and CSV outputs.
- `results/post_meeting/step2_initialization_d2/`: updated Step 2 robustness tables, heatmaps, and supplementary state plots.
- `results/post_meeting/step3_heatmap/raw/`: original Step 3 heatmap outputs and boundary CSVs.
- `results/post_meeting/step3_heatmap/fixed/`: corrected Step 3 presentation-ready plots and cleaned boundary CSV.
- `docs/`: concise written summaries and the equation reference PDF.
- `presentations/`: slide decks and PDFs that capture the project narrative in a presentation-ready format.

## Recommended upload scope

If you want a clean GitHub repository, this folder is already organized for that purpose. A good first upload would include everything here.

If you want an even lighter portfolio version, prioritize:

- `README.md`
- `Project.toml`
- `code/`
- `results/beta6/`
- `results/post_meeting/step2_initialization_d2/`
- `results/post_meeting/step3_heatmap/fixed/`
- `docs/`

## Notes

- This package intentionally omits zip archives, virtual environments, and a large number of duplicate outputs.
- The post-meeting materials are organized by update round so newer files replace older flat duplicates.
- The included figures and CSV files are selected for representativeness, not completeness.
- The original full workspace remains unchanged outside this curated folder.
