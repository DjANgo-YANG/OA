# Contributing

## Script changes

- Keep input, configuration, analysis and output sections clearly separated.
- Read paths from `config/paths.R` or environment variables.
- Write generated files only under `results/`.
- Set and document random seeds where applicable.
- Do not commit raw data, RDS objects, patient identifiers or credentials.
- Update `script_manifest.csv` when replacing a source snapshot.

## Pull requests

Describe the affected workflow stage, required inputs, expected outputs and validation performed. For analysis changes, include a concise comparison of the old and new result rather than committing all generated figures.
