# GEM 自动重建方法全面调研报告

> 调研日期：2026-05-28  
> 涵盖 9 个 GEM 重建工具/平台 + 全流程 SOTA 工具纵向分析

---

## 目录

1. [CarveMe](#1-carveme)
2. [ModelSEED & KBase（平台 vs 数据库 vs 工作流）](#2-modelseed--kbase)
3. [gapseq](#3-gapseq)
4. [AGORA2 / DEMETER](#4-agora2--demeter)
5. [GEMsembler](#5-gemsembler)
6. [gempipe](#6-gempipe)
7. [pyFBA](#7-pyfba)
8. [Bactabolize](#8-bactabolize)
9. [九工具数据库与注释工具总览](#9-九工具数据库与注释工具总览)
10. [纵向分析：基因组→GEM 全流程 SOTA 工具矩阵](#10-纵向分析基因组gem-全流程-sota-工具矩阵)
11. [LLM 可改进模块](#11-llm-可改进模块)

---

## 1. CarveMe

**论文**: Machado D et al. (2018) *Nucleic Acids Research*, 46(15):7542-7553  
**代码**: https://github.com/cdanielmachado/carveme

### 1.1 使用的数据库

| 数据库 | 类型 | 条目数 | 存储内容 |
|--------|------|--------|---------|
| **BiGG Models 1.6** | 代谢模型集合 + 生化数据库 | **108** 个手动 curated 模型，**28,302** 个唯一反应，**9,088** 个唯一代谢物 | 标准化反应方程（含 stoichiometry）、代谢物结构、GPR 关联、区室信息、SBML L3/FBCv2 格式 |
| **CarveMe 通用模板（carving template）** | 从 BiGG 派生的通用代谢模型 | **5,532** 个反应，**2,861** 个代谢物 | 整合所有 BiGG 细菌模型的"超级模型"，经预清理（无阻塞反应、无死胡同代谢物），分 gramneg/grampos/cyanobacteria/archaea 四种模板 |
| **BiGG 蛋白序列库** | 蛋白序列 FASTA | 所有 BiGG 模型中的蛋白编码基因 | 用于 DIAMOND 同源搜索的参考蛋白序列 |

### 1.2 使用的注释工具

| 工具 | 方法 | 输入 | 输出 | 参数 |
|------|------|------|------|------|
| **DIAMOND** | 快速蛋白序列比对（双索引 + 缩减字母表） | 基因组蛋白序列（FASTA） | 比对 bitscore、e-value | `--more-sensitive --top 10`，可选 `-e 1e-20 --top 20` |
| **CarveMe 反应评分引擎** | DIAMOND bitscore → GPR 映射 → 反应分数（对数正态，中位数=1） | DIAMOND 比对结果 | 每反应的序列证据得分 | 无遗传证据的酶催化反应默认 -1，自发反应默认 0 |

### 1.3 Pipeline 步骤

1. DIAMOND 将输入蛋白序列比对到 BiGG 蛋白序列库 → 获得每个基因的 bitscore
2. 基因分数 → 蛋白分数（复合物亚基取最小分）→ 反应分数（同工酶分数求和）→ 归一化
3. MILP（混合整数线性规划）雕刻：最大化保留反应分数，约束最小生物量产率
4. 移除所有非活性反应及孤儿代谢物/基因
5. 可选 gap-filling：从模板重新添加反应（优先有部分遗传证据的）

---

## 2. ModelSEED & KBase

**论文**: Henry CS et al. (2010) *Nat Biotechnol*, 28(9):977-982; Seaver SMD et al. (2021) *NAR*, 49(D1):D575-D588  
**平台**: https://modelseed.org / https://kbase.us

### 2.0 关键概念澄清：KBase 是平台，ModelSEED 是数据库+工作流

| 实体 | 性质 | 说明 |
|------|------|------|
| **KBase** | **工作流平台** | DOE Systems Biology Knowledgebase——云端 Web 环境，提供集成化 App（包括基因组注释、代谢模型构建、FBA、群落建模等）。不是数据库，不存储生化反应数据。 |
| **ModelSEED Biochemistry Database** | **数据库** | 33,978 化合物 + 36,645 反应，来源于 KEGG/MetaCyc/EcoCyc/Plant BioCyc 等。GitHub 托管（ModelSEED/ModelSEEDDatabase），TSV 格式，含热力学参数（基团贡献法 ΔG）。 |
| **ModelSEED 重构工作流** | **工作流** | 最初是独立命令行工具，现作为 KBase 中的 App 运行（"Build Metabolic Model"）。使用 ModelSEED 数据库 + RAST 注释 → 草稿模型。 |
| **关系** | — | 用户上传基因组到 KBase → RASTtk 注释（获得功能角色）→ ModelSEED 重构 App 查询 ModelSEED DB → 产生 GEM |

### 2.1 使用的数据库

| 数据库 | 类型 | 条目数 | 存储内容 |
|--------|------|--------|---------|
| **ModelSEED Biochemistry DB** | 综合生化数据库 | **33,978** 化合物，**36,645** 反应 | 反应方程、stoichiometry、化合物结构、热力学 ΔG（基团贡献法 1mM/25°C/pH7）、标准化命名空间（rxnXXXXX_c0, cpdXXXXX_c0）、KEGG+MetaCyc+EcoCyc 交叉引用 |
| **SEED 子系统数据库** | 功能注释数据库 | **>1,500** curated subsystems | 通路定义、功能角色、FIGfam 分类、基因-子系统关联 |
| **FIGfam** | 蛋白家族数据库 | **~185,000** 蛋白家族 | 覆盖 ~1,600 万蛋白编码基因，用于功能角色分配 |

### 2.2 使用的注释工具

| 工具 | 方法 | 输入 | 输出 | 性能 |
|------|------|------|------|------|
| **RASTtk** | 基于子系统的 BLAST 同源搜索 + 受控词汇表 | 细菌/古菌基因组 FASTA | 功能角色、EC 编号、GO term、子系统分类、FIGfam ID | 对模式生物好，非模式生物可能遗漏 pathway-specific 基因 |
| **ModelSEED 反应映射引擎** | 功能角色 → 蛋白复合体（many-to-many）→ 反应（many-to-many） | RASTtk 功能角色列表 | GPR 关联、候选反应池 | — |

### 2.3 Pipeline 步骤

1. RASTtk 注释基因组 → 获得功能角色 + 子系统分类
2. 功能角色 → 酶复合体（many-to-many：一个角色可参与多个复合体）
3. 酶复合体 → 反应（查询 ModelSEED DB）
4. 添加通用/自发反应 + 热力学约束（基团贡献法 ΔG）
5. 生成生物量方程（DNA/RNA/蛋白/脂质/细胞壁/辅因子）
6. MILP gap-filling（含惩罚系数：KEGG 外反应、未知结构、转运蛋白高惩罚）
7. GrowMatch 优化：修改 gap-filling 结果以拟合基因必要性和生长表型（准确率 66%→87%）

---

## 3. gapseq

**论文**: Zimmermann J et al. (2021) *Genome Biology*, 22:81  
**代码**: https://github.com/jotech/gapseq

### 3.1 使用的数据库

| 数据库 | 条目数 | 存储内容 |
|--------|--------|---------|
| **gapseq 序列数据库（蛋白）** | **130,671** 条唯一蛋白序列（含 UniPac 0.9 的 111,542 + TCDB 的 19,129） | 每个序列关联到一个或多个反应/通路，用于 tblastn 同源搜索 |
| **TCDB（Transporter Classification DB）** | **19,129** 条转运蛋白序列 | 转运蛋白分类、底物信息、家族分配 |
| **gapseq 生化反应数据库** | **14,287** 个总反应（含转运），**7,570** 个代谢物 | 反应定义、化合物结构、通路关联、来源于 MetaCyc + KEGG + ModelSEED |
| **gapseq 通用模型** | **10,194** 个反应，**3,337** 个代谢物 | 去除死端代谢物后的预清理模型 |
| **MetaCyc** | **17,208** 个反应，**20,296** 个化合物（MNXref 统计） | 人工 curated 通路定义、关键反应和完整性标准 |
| **KEGG** | **11,161** 个反应，**18,673** 个化合物（MNXref 统计） | 通路图、KO 映射 |
| **UniPac 0.5（可选）** | **1,131,132** 个未审阅序列簇 | 额外的蛋白序列来源 |

### 3.2 使用的注释工具

| 工具 | 方法 | 输入 | 输出 | 参数 |
|------|------|------|------|------|
| **tblastn**（内部调用） | 蛋白查询 vs 翻译后的基因组核苷酸 | 基因组 FASTA | 比对 bitscore | bitscore 阈值 `-b 200`，coverage ≥ 75% |
| **gapseq find** | 通路与反应预测：序列同源 → 通路完整性检查 | tblastn 结果 | 通路存在/缺失、反应存在/缺失 | 无 hints 时通路 ≥80%；有关键反应时 ≥66% |
| **gapseq find-transport** | 转运蛋白预测：tblastn vs TCDB | 基因组 FASTA | 转运反应预测 | 同上 bitscore 阈值 |
| **gapseq draft** | 加权网络重建：bitscore → 反应权重线性转换（l=50→max, u=200→min） | find 输出 | 草稿模型 + 加权等候列表 | — |
| **gapseq medium** | 生长培养基预测：检测缺失通路 → 加入对应化合物 | 草稿模型 | 预测培养基（74 种化合物） | — |
| **gapseq fill** | 多步 gap-filling（Step1生长→Step2/2b生物量→Step3替代能源→Step4副产物） | 草稿+培养基 | gap-filled 模型 | bitscore≥b 的为核心候选（全程可用），<b 仅 Step1 |

---

## 4. AGORA2 / DEMETER

**论文**: Heinken A et al. (2023) *Nat Biotechnol*, 41:1320-1331  
**平台**: https://www.vmh.life

### 4.1 使用的数据库

| 数据库 | 条目数 | 存储内容 |
|--------|--------|---------|
| **VMH（Virtual Metabolic Human）** | **19,313** 个反应，**5,607** 个代谢物，**3,695** 个人类基因，**8,790** 个食物条目，818 个微生物 | 标准化代谢物/反应命名（VMH 命名空间）、链接 57 个外部资源（KEGG, HMDB, ChEBI 等）、整合 Recon3D（13,543 反应, 4,140 代谢物） |
| **PubSEED** | **5,438** 个微生物菌株基因组 | 子系统注释、功能角色、比较基因组数据（用于 AGORA2 中的 446 个基因功能经 35 个子系统验证） |
| **AGORA2 重构集合** | **7,302** 个菌株（**1,738** 物种，**25** 门） | 精炼后的全基因组代谢模型，平均每株 1,723 反应/1,539 代谢物/1,577 基因 |

### 4.2 使用的注释工具

| 工具 | 方法 | 输入 | 输出 |
|------|------|------|------|
| **KBase / RASTtk** | 基于子系统的 BLAST 同源搜索（同 ModelSEED） | 基因组 FASTA | 功能角色注释 + 草稿代谢模型（SBML） |
| **DEMETER 翻译表** | 将 ModelSEED 命名空间映射到 VMH 命名空间 | KBase 草稿模型 | VMH 命名模型 |
| **DEMETER 精炼** | 半自动：Gram 染色 BOF 优化 + PubSEED 比较基因组 + 实验数据传播（>1,500 物种） + 数据驱动 gap-filling + 迭代 QC 测试 | VMH 模型 | 高质量 GEM |

### 4.3 Pipeline 步骤

1. KBase 生成草稿模型（ModelSEED 工作流）
2. DEMETER 翻译：ModelSEED → VMH 命名空间
3. BOF 优化：基于 Gram 染色状态（修复 32% 草稿 gram 状态错误）
4. 实验数据整合：碳源利用/发酵产物/生长需求/摄取分泌，系统发育自动传播
5. PubSEED 比较基因组：1,000+ curated 子系统
6. 无益循环移除、数据驱动 gap-filling
7. 迭代 QC：100% 厌氧生长（精炼后）vs 27%（草稿），>99% 实验一致性

---

## 5. GEMsembler

**论文**: Matveishina EK et al. (2025) *mSystems*, 10(10):e0057425  
**代码**: https://github.com/zimmmermann-kogadeeva-group/GEMsembler

### 5.1 使用的数据库

| 数据库 | 条目数 | 存储内容 |
|--------|--------|---------|
| **BiGG Models（目标命名空间）** | **28,302** 个唯一反应，**9,088** 个唯一代谢物 | 标准化反应方程、代谢物 ID、GPR 关联——作为共识转换的 target namespace |
| **MetaNetX / MNXref 4.5** | **41,687** 个代谢物（反应中），**37,103** 个反应 | 跨 30+ 数据库的交叉引用映射（ChEBI, HMDB, KEGG, MetaCyc, ModelSEED, Reactome, Rhea, LipidMAPS, SwissLipids） |

### 5.2 使用的注释工具

| 工具 | 方法 | 输入 | 输出 |
|------|------|------|------|
| **6 级优先级代谢物转换** | 模型注释 → 数据库交叉引用 → MetaNetX → 模式匹配 → 直接 ID 检查 → 未转换留存 | 各工具输出模型 | BiGG 命名的代谢物 |
| **反应方程转换** | 仅在代谢物已转换为 BiGG 后，通过反应方程匹配（测试每个候选是否匹配已知 BiGG 反应） | 转换后的代谢物 + 原始反应方程 | BiGG 命名的反应 |
| **BLAST 基因转换** | 将各工具的基因 ID BLAST 比对到参考基因组（可选） | 基因序列 | 参考基因组 locus tag |
| **GPR 组合优化** | 逐步算法（SA）：GPR 替换 → FBA 验证；遗传算法（GA）：染色体编码 GPR 来源选择 → 最大化 AUCPR | 多工具 GPR | 优化后的共识 GPR |

---

## 6. gempipe

**论文**: Lazzari G et al. (2025) *mSystems*, 11(1):e0100725  
**代码**: https://github.com/gempipe/gempipe

### 6.1 使用的数据库

| 数据库 | 条目数 | 存储内容 |
|--------|--------|---------|
| **BiGG Models**（universe 模板基础） | **28,302** 反应，**9,088** 代谢物 | gramneg/grampos 两种模板，含通用生物量方程 |
| **ModelSEED** | **33,978** 化合物，**36,645** 反应 | 用于 gap-filling 的额外反应来源 |
| **TCDB** | **19,129** 条转运蛋白序列 | 转运反应扩展 |

### 6.2 使用的注释工具

| 工具 | 方法 | 输入 | 输出 |
|------|------|------|------|
| **Bakta**（推荐） | AFSI + UniParc/UniRef100/RefSeq | 基因组 FASTA | GFF3/GBK/EMBL（含功能注释） |
| **Prokka**（备选） | Prodigal + HMMER + BLAST + UniProtKB/Pfam/TIGRFAMs | 基因组 FASTA | GFF3/GBK/EMBL |
| **BLAST/diamond** | 同源比对找参考模型基因的直系同源 | CDS 序列 vs 参考模型基因 | 同源命中 |
| **HMM profile** | 对缺失基因进行 HMM 补充搜索 | 遗漏 CDS | 补充命中 |

### 6.3 Pipeline 步骤

1. 基因恢复（BLAST → 过滤 → HMM 补充）→ 97.8% 基因恢复完整度
2. 路径 A：参考模型扩展——BLAST 直系同源，保留匹配反应
3. 路径 B：非参考宇宙重建——BiGG 模板全序列比对
4. 混合合并（A + B 反应自动合并）
5. gap-filling + 模型 curation

---

## 7. pyFBA

**论文**: Cuevas DA et al. (2016) *Frontiers in Microbiology*, 7:907  
**代码**: https://github.com/linsalrob/PyFBA

### 7.1 使用的数据库

| 数据库 | 条目数 | 存储内容 |
|--------|--------|---------|
| **ModelSEED Biochemistry DB（本地版）** | **~27,000+** 化合物，数千反应 | 化合物结构、反应定义、EC 编号映射、酶复合物定义（功能角色→反应）、生物量方程（gramneg/grampos/mycobacterium/plant）、培养基配方、子系统定义 |

### 7.2 使用的注释工具

| 工具 | 方法 | 输入 | 输出 |
|------|------|------|------|
| **RAST** | 基于子系统的 BLAST 同源搜索 | 基因组 FASTA | 功能角色注释列表 |
| **PyFBA.filters.roles_to_reactions()** | 功能角色 → 酶复合体（many-to-many）→ 反应（many-to-many） | 功能角色列表 | 反应列表 |
| **PyFBA.gapfill 多模块** | 6 个独立模块有序运行：(1) 109 个通用反应 (2) 转运反应识别 (3) 孤儿化合物分辨 (4) 子系统补全 (5) 比较基因组学 (6) 二分法反应缩减 | 反应列表 + 培养基 | gap-filled 模型 |

---

## 8. Bactabolize

**论文**: Vezina B, Watts SC et al. (2023) *eLife*, 12:RP87406  
**代码**: https://github.com/kelwyres/Bactabolize

### 8.1 使用的数据库

| 数据库 | 条目数 | 存储内容 |
|--------|--------|---------|
| **KpSC-pan v1 参考模型** | **1,265** 基因，**2,319** 反应，**1,696** 代谢物 | 从 37 个人工 curated 菌株模型构建，覆盖 7 个 KpSC 分类群，BiGG 命名空间 |
| **KpSC-pan v2**（扩展版） | **2,403** 基因，**3,550** 反应（883 催化 + 175 转运），**360** 生长底物 | v1 + 507 分离株泛基因组（34,664 个唯一基因）+ KEGG/ModelSEED 的 1,058 额外反应；好氧生长预测 95.4%，厌氧 78.8% |

### 8.2 使用的注释工具

| 工具 | 方法 | 输入 | 输出 | 参数 |
|------|------|------|------|------|
| **Prodigal** | 基因预测 | 基因组 FASTA | CDS 蛋白序列 | — |
| **双向 BLAST best hit (BBH)** | 目标菌株 CDS ↔ 参考模型基因——双向 BLAST，互为最佳命中 | 目标 CDS + 参考模型基因 | 直系同源判定 | 比单向 BLAST 更保守 |

### 8.3 Pipeline 步骤

1. Prodigal → CDS 预测
2. BBH 找直系同源 vs 参考模型（KpSC-pan）：参考基因缺失 → 移除对应反应
3. 自动 gap-filling → 菌株特异性草稿模型
4. patch_model → 合并 gap-filling 反应
5. FBA（27+ 条件）+ sgk（单基因敲除）→ 表型预测

---

## 9. 九工具数据库与注释工具总览

### 9.1 数据库汇总

| 数据库 | 类型 | 化合物/代谢物 | 反应 | 基因/蛋白 | 使用该数据库的工具 |
|--------|------|--------------|------|-----------|-----------------|
| **BiGG Models 1.6** | 代谢模型集合 | 9,088 | 28,302 | — | CarveMe, GEMsembler, gempipe |
| **BiGG 通用模板** | 通用模型 | 2,861 | 5,532 | — | CarveMe |
| **ModelSEED Biochem DB** | 综合生化 DB | 33,978 | 36,645 | — | ModelSEED, KBase, gapseq, pyFBA, gempipe |
| **VMH** | 人/微生物代谢 DB | 5,607 | 19,313 | 3,695 人类 | AGORA2/DEMETER |
| **MetaNetX/MNXref 4.5** | 交叉引用命名空间 | 41,687（反应中） | 37,103 | — | GEMsembler |
| **TCDB** | 转运蛋白 DB | — | — | 19,129 条序列 | gapseq, gempipe |
| **MetaCyc** | 通路/反应 DB | 20,296 | 17,208 | — | gapseq (via ModelSEED) |
| **KEGG** | 综合 DB | 18,673 | 11,161 | — | gapseq (via ModelSEED) |
| **PubSEED** | 基因组注释 DB | — | — | 5,438+ 菌株 | AGORA2 |
| **RAST/SEED 子系统** | 功能注释 DB | — | — | ~185,000 FIGfams | ModelSEED, pyFBA |
| **KpSC-pan v2** | 参考代谢模型 | 1,696 | 3,550 | 2,403 | Bactabolize |
| **gapseq 序列 DB** | 蛋白序列 DB | — | — | 130,671 条 | gapseq |
| **gapseq 生化 DB** | 反应 DB | 7,570 | 14,287 | — | gapseq |
| **Rhea** | 反应 DB | 14,654 化合物 | 17,783 | — | 独立/分析 |
| **BRENDA** | 酶信息 DB | — | — | 7,500+ EC | 独立/分析 |

### 9.2 注释工具汇总

| 工具 | 用途 | 方法 | 输入 | 输出 | 使用该工具的工作流 |
|------|------|------|------|------|-----------------|
| **Prodigal/Pyrodigal** | 基因预测 | 动态规划 + RBS 打分 | 基因组 FASTA | CDS (GFF/GBK) | Bactabolize, Prokka, gempipe(间接) |
| **Bakta** | 全基因组注释 | AFSI + UniParc/UniRef | 基因组 FASTA | GFF3/GBK/EMBL | gempipe |
| **Prokka** | 全基因组注释 | Prodigal + HMMER + BLAST | 基因组 FASTA | GFF3/GBK/EMBL | gempipe(备选) |
| **RAST/RASTtk** | 子系统注释 | 基于子系统的 BLAST 同源 | 基因组 FASTA | 功能角色/EC/GO | ModelSEED, KBase, AGORA |
| **DIAMOND** | 蛋白同源搜索 | 双索引 + 缩减字母表 | 蛋白 FASTA | bitscore/e-value | CarveMe, gapseq |
| **BLASTp** | 蛋白同源搜索 | 经典双序列比对 | 蛋白 FASTA | bitscore/e-value | GEMsembler |
| **tblastn** | 翻译后同源搜索 | 蛋白查询 vs 翻译基因组 | 基因组 FASTA | bitscore | gapseq |
| **BBH** | 直系同源判定 | 双向 BLAST best hit | CDS + 参考基因 | 直系同源 | Bactabolize |
| **HMM profile** | 补充搜索 | HMM 轮廓搜索 | CDS FASTA | 补充命中 | gempipe |
| **CLEAN** | EC 预测 | ESM-1b/ESM-2 + 对比学习 | 蛋白序列 | EC 编号 | 独立/可集成 |
| **EZSpecificity** | 酶底物特异性 | SE(3)-等变 GNN + 交叉注意力 | 酶序列+底物 SMILES | 催化概率 | 独立/可集成 |
| **MEMOTE** | 模型质控 | 标准化测试套件 | SBML | HTML 报告 + 分数 | 独立/通用 |

### 9.3 工作流范式总览

| 工具 | 范式 | 核心操作 | 速度（相对） | 输出质量特征 |
|------|------|---------|-----------|------------|
| **CarveMe** | 自上而下雕刻 | MILP 从 BiGG 移除无证据反应 | 极快（~30s/10 模型） | 假阴性偏高 32%，死胡同最低 |
| **ModelSEED/KBase** | 自下而上组装 | RAST 注释 + GPR 映射 + MILP + GrowMatch | 中等（~3min） | 最好的热力学数据 |
| **gapseq** | 自下而上 + 权重 | bitscore 连续权重 + 5 步 gap-fill | 慢（~5.5h/10 模型） | 假阴性最低 6%，GPR 覆盖最高 79-86% |
| **AGORA2/DEMETER** | 草稿+半自动精炼 | KBase 草稿 + 实验数据驱动精炼 | 数小时/菌株 | 最高实验一致性 >99% |
| **GEMsembler** | 多工具共识 | 特征投票 + GA 优化 GPR | 取决于输入 | 可超越 gold-standard |
| **gempipe** | 混合（参考+宇宙） | 双路并行 + 自动合并 | 中 | 最高基因恢复 97.8%，最低孤儿反应 4.1% |
| **pyFBA** | 功能角色→反应映射 | 三层 many-to-many 映射 + 6 模块 gap-fill | 快 | 过程透明、教育价值高 |
| **Bactabolize** | 参考驱动还原 | BBH 直系同源 + 去除缺失 | <3min/基因组 | 最高底物预测率 0.97（KpSC 内） |

---

## 10. 纵向分析：基因组→GEM 全流程 SOTA 工具矩阵

### 10.0 完整工作流概览

基因组 FASTA → 10.1 基因预测（Pyrodigal）→ 10.2 全基因组注释（Bakta）→ 10.3 酶/非酶分类（SOLVE）→ 10.4 转运蛋白分类（TooT-BERT-CNN-T）→ 10.5 EC 编号预测（GraphEC）→ 10.6 KO 分配（eggNOG-mapper v2）→ 10.7 GO 分配（PFresGO）→ 10.8 同源性搜索（MMseqs2）→ 10.9 酶-底物特异性（EZSpecificity）→ 10.10 反应数据库映射（Rhea）→ 10.11 Stoichiometry 检索（Rhea）→ 10.12 生物量方程（gapseq）→ 10.13 Gap-filling（DNNGIOR）→ 10.14 质控（MEMOTE）→ 10.15 生长验证（BacDive）→ 10.16 基因必要性（COBRApy+GEMsembler）

每个环节需查完整表格见报告中 §10.1-§10.16。

---

## 11. LLM 可改进模块

### 11.1 改进潜力排序：注释冲突裁决、Gap-filling 生物学判断、MEMOTE 迭代修复、代谢物命名空间映射、多工具语义加权、通路完整性补全、参考模型质控、文献驱动验证

### 11.2 各工具中 LLM 可直接改进的决策点

| 工具 | 存在的人工决策 | LLM 替代方案 |
|------|--------------|------------|
| **CarveMe** | bitscore 阈值选择、gap-filling 结果验证 | LLM 根据物种背景调整阈值、生物学验证 gap-filling |
| **ModelSEED/KBase** | GrowMatch 拟合结果的生物学验证 | LLM 替代 MARA 式推理：反应是否在生物学上合理？ |
| **gapseq** | gap-filling 中低 bitscore 反应的取舍 | LLM 查询通路上下文判断取舍 |
| **AGORA/DEMETER** | PubMed 手动文献检索（最大瓶颈） | LLM + RAG 自动化文献检索和结构化提取 |
| **GEMsembler** | coreX 等权投票 | LLM 基于分类群/通路上下文加权 |
| **gempipe** | 双路合并规则 | LLM 基于证据质量动态决定合并策略 |
| **pyFBA** | 6 个 gap-fill 模块的选择和排序 | LLM 诊断后自适应启停模块 |
| **Bactabolize** | 参考模型构建的 curation | LLM 辅助构建和验证物种特异参考模型 |

---

## 参考文献

### 核心工具论文
1. Machado D et al. (2018). *NAR*. [CarveMe]
2. Henry CS et al. (2010). *Nat Biotechnol*. [ModelSEED]
3. Seaver SMD et al. (2021). *NAR*. [ModelSEED DB]
4. Zimmermann J et al. (2021). *Genome Biol*. [gapseq]
5. Heinken A et al. (2023). *Nat Biotechnol*. [AGORA2/DEMETER]
6. Matveishina EK et al. (2025). *mSystems*. [GEMsembler]
7. Lazzari G et al. (2025). *mSystems*. [gempipe]
8. Cuevas DA et al. (2016). *Front Microbiol*. [pyFBA]
9. Vezina B, Watts SC et al. (2023). *eLife*. [Bactabolize]

### 数据库与平台
10. Norsigian CJ et al. (2020). *NAR*. [BiGG Models 1.6]
11. Moretti S et al. (2024). *NAR*. [Rhea]
12. Tran VDT et al. (2023). *NAR*. [MetaNetX/MNXref 4.5]
13. Overbeek R et al. (2014). *NAR*. [RAST/SEED]

### SOTA 注释/预测工具
14. Cui H et al. (2025). *Nature*. [EZSpecificity]
15. Yu T et al. (2023). *Science*. [CLEAN]
16. Zhao et al. (2024). *Nat Commun*. [GraphEC]
17. Schwengers O et al. (2021). *PLoS Comput Biol*. [Bakta]
18. Hyatt D et al. (2010). *BMC Bioinform*. [Prodigal]
19. Lieven C et al. (2020). *Nat Biotechnol*. [MEMOTE]

### 对比评估
20. Hsieh YE et al. (2024). *npj Syst Biol Appl*.
