# GEM 自动重建方法详细调研报告

> 调研日期：2026-05-24  
> 目标：为 agem 项目提供方法论基础，识别 LLM 可以改进的模块

---

## 目录

1. [CarveMe — 自上而下雕刻法](#1-carveme)
2. [ModelSEED — 自下而上组装法](#2-modelseed)
3. [gapseq — 权重证据预测法](#3-gapseq)
4. [AGORA/DEMETER — 半自动化精炼法](#4-agorademeter)
5. [GEMsembler — 多工具共识组装法](#5-gemsembler)
6. [gempipe — 混合重建 + pan/multi-strain 分析](#6-gempipe)
7. [pyFBA — 功能角色→反应映射 + 迭代 gap-filling](#7-pyfba)
8. [ModelSEED/KBase/RAST — 完整平台生态系统](#8-modelseedkbacerast)
9. [Bactabolize — 参考模型驱动的快速高性能重建](#9-bactabolize)
10. [九工具系统性对比](#10-九工具系统性对比)
11. [LLM 可改进模块矩阵](#11-llm-可改进模块矩阵)

---

## 1. CarveMe — 自上而下雕刻法

**论文**: Machado et al. (2018), *Nucleic Acids Research*  
**代码**: https://github.com/cdanielmachado/carveme

### 1.1 核心方法论

CarveMe 采用独特的**自上而下（top-down）雕刻法**：从 BiGG 数据库整合的通用代谢模板出发，用 DIAMOND 将基因组蛋白序列与 BiGG 蛋白数据库比对，将比特分数经 GPR 映射转化为反应分数，最后通过 MILP 优化移除无遗传证据的反应，保留能产生生物量的最小网络。

### 1.2 关键特点

- 命名空间：BiGG
- 速度：极快（~30秒/10个模型）
- 假阴性率：高（32%）
- 死胡同代谢物：低（模板已预清理）

### 1.3 主要局限

1. 模板盲区：不能发现模板中不存在的新反应
2. GPR 映射粗糙：仅依赖 DIAMOND 同源性
3. 无 EC/KO 验证：不检查酶功能的生物学一致性

---

## 2. ModelSEED — 自下而上组装法

**论文**: Henry CS et al. (2010), *Nature Biotechnology*  
**平台**: https://modelseed.org / KBase

### 2.1 核心方法论

RAST 注释 → GPR 关联（基于 33,978 化合物、36,645 反应的生物化学数据库）→ 草稿网络组装 → 生物量方程生成 → MILP gap-filling → GrowMatch 优化（准确率从 66%→87%）

### 2.2 关键特点

- 数据库最丰富：33,978 化合物 + 36,645 反应
- 热力学参数：内置基团贡献法 ΔG
- 命名空间：ModelSEED

### 2.3 主要局限

- 数据库更新频率低于 KEGG/BiGG
- MILP 可能引入生物学上不合理的反应

---

## 3. gapseq — 权重证据预测法

**论文**: Zimmermann J et al. (2021), *Genome Biology*  
**代码**: https://github.com/jotech/gapseq

### 3.1 核心方法论

四步流程：通路与反应预测（bitscore≥200, coverage≥75%）→ 转运体预测 → 加权草稿重建（bitscore→权重线性转换）→ 生长培养基预测（74种化合物）→ 多步 gap-filling

### 3.2 核心创新

不使用 hard cutoff，而是根据 bitscore 持续分配权重：bitscore≤l(50)→最大权重，bitscore≥u(200)→最小权重

### 3.3 关键特点

- 假阴性率最低（6%）
- GPR 覆盖率最高（79-86%）
- 速度慢（~5.5h/10模型）

---

## 4. AGORA/DEMETER — 半自动化精炼法

**论文**: Heinken A et al. (2023), *Nature Biotechnology*  
**平台**: https://www.vmh.life

### 4.1 核心方法论

KBase 草稿重建 → DEMETER 精炼：命名翻译 → BOF 优化（gram 状态）→ 实验数据整合（>1500物种）→ PubSEED 比较基因组学 → 数据驱动 gap-filling → 迭代 QC

### 4.2 质量指标

- 7,302 菌株、1,738 物种、25 门
- 厌氧生长率 100%（精炼后）vs 27%（草稿）
- 实验数据一致性 >99%

---

## 5. GEMsembler — 多工具共识组装法

**论文**: Matveishina EK et al. (2025), *mSystems*  
**代码**: https://github.com/zimmmermann-kogadeeva-group/GEMsembler

### 5.1 核心方法论

4工具输入 → 特征转换到 BiGG 统一命名空间 → supermodel 组装 → coreX 共识模型（core2/core3/core4）→ GPR 组合优化（SA + GA）

### 5.2 关键发现

- ~25% 特征被 4 工具共同支持
- core3 超越黄金标准（L. plantarum）
- GA 优化 AUCPR: 0.556→0.711

---

## 6. gempipe — 混合重建 + pan/multi-strain 分析

**论文**: Lazzari G et al. (2025), *mSystems*  
**代码**: https://github.com/gempipe/gempipe

### 6.1 核心方法论

双路重建：参考模型扩展（BLAST 直系同源）+ 非参考宇宙重建（BiGG+ModelSEED+...全宇宙模板），自动合并两路结果。

### 6.2 性能基准

- 参考覆盖率 91.2%（vs Bactabolize 94.5%）
- 代谢预测准确率 92.3%（最高）
- 基因恢复完整度 97.8%（最高）
- 孤儿反应率 4.1%（最低）

---

## 7. pyFBA — 功能角色→反应映射 + 迭代 gap-filling

**论文**: Cuevas DA et al. (2016), *Frontiers in Microbiology*  
**代码**: https://github.com/linsalrob/PyFBA

### 7.1 核心方法论

功能角色→酶复合体→反应的三层映射体系（many-to-many），6 个独立 gap-filling 模块有序组合：universal reactions, transport, orphan compounds, subsystem completion, comparative genomics, reaction reduction

### 7.2 关键特点

- 过程透明，每步可审查
- gap-filling 模块可定制
- 依赖 Model SEED 数据库

---

## 8. ModelSEED/KBase/RAST — 完整平台生态系统

### 8.1 KBase 代谢模型构建工作流

上传基因组 → RASTtk 注释（子系统上下文）→ Build Metabolic Model → Gapfill → FBA/模型分析 → 导出 SBML

### 8.2 RAST 注释引擎

- 基于 >1,500 curated subsystems 的子系统注释
- 输出：功能角色、EC、GO、FIGfam
- 独特优势：通路上下文——不是单基因注释

---

## 9. Bactabolize — 参考驱动快速高性能重建

**论文**: Vezina B, Watts SC et al. (2023), *eLife*  
**代码**: https://github.com/kelwyres/Bactabolize

### 9.1 核心方法论

纯参考模型驱动的还原方法：高质量 pan-genome 参考模型 → BBH 找直系同源 → 参考基因缺失则移除对应反应 → 自动 gap-filling → 底物预测/基因必要性分析

### 9.2 关键特点

- 底物预测准确率 0.97（最高于 KpSC 内）
- 速度 <3 min/基因组
- 假阳性最少
- 局限：无法发现参考模型外的新反应

---

## 10. 九工具系统性对比

### 10.1 范式对比

| 维度 | CarveMe | ModelSEED/KBase | gapseq | AGORA/DEMETER | GEMsembler | gempipe | pyFBA | Bactabolize |
|------|---------|-----------------|--------|--------------|-----------|---------|-------|------------|
| 范式 | 自上而下 | 自下而上 | 加权自下而上 | 草稿+精炼 | 多工具共识 | 混合（参考+宇宙） | 功能角色→反应 | 参考驱动还原 |
| gap-filling | MILP | MILP+GrowMatch | 加权多步 | 数据驱动 | 无需 | 双路补全 | 多模块迭代 | 参考 gap-fill |
| 速度 | ~30s | ~3min | ~5.5h | 数小时 | 取决于输入 | 中 | 快 | <3min |
| pan/多菌株 | 无 | 无 | 无 | 无 | 无 | 内置 | 无 | 批量设计 |

### 10.2 输出质量对比

| 指标 | CarveMe | gapseq | Bactabolize | gempipe |
|------|---------|--------|-----------|----------|
| GPR 覆盖率 | 58-78% | 79-86% | 92.4% | **97.8%** |
| 底物预测率 | 0.78 | 0.82 | **0.97** | 0.92 |
| 死胡同代谢物 | **最低** | 最高 | 低 | 最低 |
| 孤儿反应率 | 7.3% | 高 | 5.9% | **4.1%** |

---

## 11. LLM 可改进模块矩阵

### 11.1 改进潜力排序

| 排名 | 模块 | LLM 优势 | 来源 |
|------|------|---------|------|
| 1 | 注释冲突裁决 | 语义理解超越序列分数 | 全工具 |
| 2 | gap-filling 生物学判断 | 通路/分类群上下文 | pyFBA 子系统补全 |
| 3 | MEMOTE 迭代修复 | 诊断+修复闭环 | 全工具 |
| 4 | 代谢物映射 | 化学结构语义匹配 | 全工具 |
| 5 | 多工具语义加权 | 上下文感知权重 | GEMsembler+gempipe |
| 6 | 通路完整性补全 | 通路上下文+补全推理 | pyFBA 子系统 |
| 7 | 参考模型质量验证 | 自动化 DEMETER 式测试 | Bactabolize+AGORA |
| 8 | pan-strain 差异解释 | 代谢分歧→表型推理 | gempipe |

### 11.2 新工具带来的 LLM 机会

| 来源 | 可借鉴模块 | LLM 可实现方式 |
|------|-----------|--------------|
| gempipe | 双路合并（参考 + 宇宙） | LLM 决定参考模型反应保留 + 宇宙额外反应 |
| pyFBA | 功能角色→复合体→反应映射 | LLM 裁决多对多映射争议 |
| pyFBA | 子系统补全 gap-filling | LLM 查询通路完整性 |
| Bactabolize | 参考模型质量控制 | LLM 验证成分一致性 |

---

## 参考文献

1. Machado D et al. (2018). *NAR*. [CarveMe]
2. Henry CS et al. (2010). *Nat Biotechnol*. [ModelSEED]
3. Seaver SMD et al. (2021). *NAR*. [ModelSEED DB]
4. Zimmermann J et al. (2021). *Genome Biol*. [gapseq]
5. Heinken A et al. (2023). *Nat Biotechnol*. [AGORA2/DEMETER]
6. Matveishina EK et al. (2025). *mSystems*. [GEMsembler]
7. Lazzari G et al. (2025). *mSystems*. [gempipe]
8. Cuevas DA et al. (2016). *Front Microbiol*. [pyFBA]
9. Vezina B, Watts SC et al. (2023). *eLife*. [Bactabolize]
10. Hsieh YE et al. (2024). *npj Syst Biol Appl*. [comparison]
11. Bansal P et al. (2024). *NAR*. [Rhea]
12. Cui H et al. (2025). *Nature*. [EZSpecificity]
13. Yu T et al. (2023). *Science*. [CLEAN]
14. RAST/RASTtk, KBase Documentation. [KBase platform]
