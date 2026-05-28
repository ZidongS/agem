# agem — AI-driven Genome-scale Metabolic Model Reconstruction

A comprehensive survey and technical design for next-generation AI-driven GEM reconstruction.

## Repository Contents

### Research Report
- **GEM-重建方法调研报告.md** — Detailed survey of 9 GEM reconstruction tools/platforms:
  CarveMe, ModelSEED/KBase/RAST, gapseq, AGORA/DEMETER, GEMsembler,
  gempipe, pyFBA, Bactabolize
  - For each tool: databases used (exact entry counts), annotation tools (input/output/method), pipeline steps
  - Cross-tool database comparison (15 databases, entry counts, which tool uses which)
  - Annotation tool comparison (12 tools, input/output/method/workflow association)
  - Longitudinal analysis: 16-step genome-to-GEM workflow with SOTA tools per step
  - LLM improvement opportunity matrix

### Technical Designs
- **agem_llm-技术路线-v4.md** — Final technical design for the agem_llm pipeline (7-stage, 8-database collaborative system)

### Implementation Plans
- **docs/superpowers/plans/** — 3-phase implementation plan (18 tasks total)

## Key Clarification
- **KBase** is a workflow **platform** (web GUI for genome→GEM)
- **ModelSEED** is both a biochemistry **database** (33,978 compounds, 36,645 reactions) AND a reconstruction **workflow** (runs as an App within KBase)
