# agem — AI-driven Genome-scale Metabolic Model Reconstruction

A comprehensive survey and technical design for next-generation AI-driven GEM reconstruction.

## Repository Contents

### Research Report
- **[GEM-重建方法调研报告.md](GEM-重建方法调研报告.md)** — Detailed survey of 9 GEM reconstruction tools:
  CarveMe, ModelSEED/KBase/RAST, gapseq, AGORA/DEMETER, GEMsembler,
  gempipe, pyFBA, Bactabolize
  - Methodological comparison
  - Output quality benchmarks
  - LLM-improvement opportunity matrix

### Technical Designs
- **[agem_llm-技术路线-v4.md](agem_llm-技术路线-v4.md)** — Final technical design for the agem_llm pipeline:
  - 7-stage pipeline: protein pre-filtering → multi-tool annotation → multi-database GPR assembly
    → SBML construction → topology + template comparison → GPR metadata store → LLM refinement
  - 8-database collaborative decision system
  - Synthesizes strengths of CarveMe (template completeness) + gapseq (GPR coverage)

### Implementation Plans
- **[docs/superpowers/plans/](docs/superpowers/plans/)** — 3-phase implementation plan:
  - Phase 1: Tool layer + database clients (9 tasks)
  - Phase 2: Pipeline core — GPR builder, SBML builder, topology, LLM refinement (7 tasks)
  - Phase 3: Benchmark framework vs CarveMe/gapseq (2 tasks)

## Key Design Principles

1. **LLM as post-hoc refiner, not upfront filter** — Deterministic multi-database GPR assembly first; LLM polishes at the end
2. **No KEGG dependency** — Primary reaction database is Rhea (free REST API, 17,783 reactions)
3. **Multi-evidence fusion** — Every GPR entry backed by up to 7 evidence dimensions
4. **Fills known tool gaps** — CarveMe universal template comparison finds missing reactions
5. **Experimental validation** — BacDive phenotype data used for both input constraints AND quality control
