# GitHub-Ready MacLean Lab Package

This folder is a curated subset of the larger `MACLEAN_LAB` workspace, assembled specifically for GitHub upload and portfolio-style presentation.

It keeps the strongest representative materials from the project:

- core Julia modeling scripts
- post-meeting refinement scripts
- representative output figures and CSV summaries
- presentation builder scripts
- a small set of summary documents and PDFs

## Why this folder exists

The full project workspace contains many intermediate outputs, repeated drafts, zip archives, and presentation artifacts. This package trims that down to a cleaner set of materials that still shows the scope of the work:

- model development
- parameter sweeps
- stochastic and deterministic analyses
- post-meeting refinement
- presentation-ready summaries

## Folder structure

- `code/core/`: main Julia analysis scripts covering CTL comparison, beta sweeps, parameter sweeps, and SDE refinement.
- `code/post_meeting/`: scripts from the post-meeting dense sweep and initialization-robustness follow-up.
- `code/presentation_builders/`: Python scripts used to generate or assemble presentation materials.
- `results/ctl_compare/`: representative outputs from the direct-vs-coupled CTL comparison.
- `results/beta5/`: representative beta-5 deterministic and stochastic outputs.
- `results/beta6/`: representative refined beta-6 outputs and heatmap summaries.
- `results/parameter_sweeps/`: representative wide-sweep figures and CSV outputs.
- `results/post_meeting/`: representative figures and tables from the post-meeting refinement phase.
- `docs/`: concise written summaries and the equation reference PDF.
- `presentations/`: presentation PDFs that capture the project narrative in slide format.

## Recommended upload scope

If you want a clean GitHub repository, this folder is already organized for that purpose. A good first upload would include everything here.

If you want an even lighter portfolio version, prioritize:

- `README.md`
- `Project.toml`
- `code/`
- `results/beta6/`
- `results/post_meeting/`
- `docs/`

## Notes

- This package intentionally omits zip archives, virtual environments, and a large number of duplicate outputs.
- The included figures and CSV files are selected for representativeness, not completeness.
- The original full workspace remains unchanged outside this curated folder.
