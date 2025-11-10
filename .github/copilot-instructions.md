# Copilot Instructions for unstudied_genes

## Project Overview
Early-stage bioinformatics analysis project focused on protein structure predictions for unstudied genes (PUGs) using PLDDT (predicted local distance difference test) scoring from AlphaFold API.

## Key Components & Data Flow

### Main Workflow: `plddt_for_PUGs.ipynb`
Single Jupyter notebook coordinating:
1. **PUG Dataset** - Load S. pombe unstudied genes (lacks functional annotation)
2. **PLDDT Scoring** - Query AlphaFold API to retrieve predicted protein confidence scores
3. **Control Genes** - Load reference gene set (all S. pombe genes or curated controls)
4. **Comparative Analysis** - Compare PLDDT distributions between PUGs and controls

**Critical**: PLDDT scores indicate protein structure prediction confidence (0-100); higher = more reliable predictions. This is the core metric for downstream prioritization.

## Environment & Dependencies

### Python Environment
- **Manager**: Conda (configured in `.vscode/settings.json`)
- **Python Version**: (Configure as needed; recommend 3.9+)
- **Key Dependencies** (expected):
  - `requests` - for AlphaFold API calls
  - `pandas` - data manipulation & analysis
  - `matplotlib` / `seaborn` - visualization of PLDDT distributions
  - `numpy` - numerical computations

### Setup
When adding new cells to the notebook, ensure API credentials/endpoints are clearly documented and never hardcoded.

## Conventions & Patterns

### Notebook Structure
- Comment sections clearly mark major workflow stages (see existing `#Load` comments)
- Each section handles: data loading → processing → analysis/visualization
- Cell outputs should produce: summary statistics, plots, or filtered datasets

### S. pombe Genomics Context
- Organism: *Saccharomyces pombe* (fission yeast)
- Use standard gene identifiers (e.g., SPAC####, SPBC####)
- "Unknown genes" = functionally unannotated genes; PLDDT may correlate with protein conservation/function

## Integration Points

- **AlphaFold API**: Primary external dependency for PLDDT retrieval
  - Handle rate limits & timeouts gracefully
  - Cache results to avoid redundant calls
- **Gene Databases**: S. pombe genome reference (PomBase or local annotation files)

## Development Tips
- Incrementally test API calls on small gene subsets before full runs
- Document any gene ID mappings or data transformations clearly
- Use descriptive variable names reflecting biological entities (e.g., `pug_plddt_scores` not `data1`)
