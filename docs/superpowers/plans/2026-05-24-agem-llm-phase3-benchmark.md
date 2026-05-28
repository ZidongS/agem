# agem_llm Phase 3: Benchmark 评估框架

**Goal:** 构建 agem_llm vs CarveMe vs gapseq 三工具对比 benchmark 框架。

**Architecture:** 对每个测试基因组跑三个工具，统一评估 MEMOTE、基因必要性 AUCPR、生长表型、BacDive 一致性、多工具 Jaccard。

**Tech Stack:** cobrapy, memote, numpy

---

## Task 1: 评估框架
- Create: agem_llm/benchmark/eval_framework.py
- GEMEvaluator: evaluate_gene_essentiality(), evaluate_growth_phenotypes(), evaluate_annotation_quality(), evaluate_bacdive_consistency() [新]

## Task 2: 批量 Benchmark 运行脚本
- Create: agem_llm/benchmark/run_benchmarks.py
- 对 Tier 1 基因组（E. coli, B. subtilis, P. aeruginosa, S. aureus）：
  1. 运行 CarveMe（CLI）
  2. 运行 gapseq（CLI）
  3. 运行 agem_llm pipeline
  4. 汇总对比报告（CSV + Markdown）

## 评估指标

| 层级 | 方法 | 目标 |
|------|------|------|
| L0 | SBML 语法 | libsbml 解析验证 |
| L1 | MEMOTE 结构 | 质量平衡、死胡同、GPR 覆盖率 |
| L2 | FBA 生长 | 生物量生产 |
| L3 | BacDive 碳源 | 预测 vs 实验碳源利用一致性 |
| L4 | 基因必要性 | AUCPR（需要 Keio 数据） |
| L5 | 多工具 Jaccard | CarveMe/gapseq/agem_llm 结构对比 |
