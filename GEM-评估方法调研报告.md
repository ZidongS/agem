# GEM 评估方法论全面调研报告

> 调研日期：2026-05-28
> 涵盖基因必要性、生长表型、MEMOTE 测试套件、生物量一致性、多工具共识、实验数据交叉验证，以及 9 个工具的自我评估方法

---

## 目录

1. [基因必要性预测评估](#1-基因必要性预测评估)
2. [生长表型预测评估](#2-生长表型预测评估)
3. [MEMOTE 测试套件](#3-memote-测试套件)
4. [生物量一致性检查](#4-生物量一致性检查)
5. [多工具共识指标](#5-多工具共识指标)
6. [实验数据交叉验证](#6-实验数据交叉验证)
7. [各工具自我评估方法](#7-各工具自我评估方法)

---

## 1. 基因必要性预测评估

### 1.1 实验数据集

| 数据集 | 物种 | 必要基因数 | 非必要基因数 | 总基因数 | 方法 | 来源 |
|--------|------|-----------|------------|---------|------|------|
| **Keio Collection** | *E. coli* K-12 BW25113 | 303（候选） | 3,985（成功敲除） | 4,288（目标） | 单基因框内敲除（Kan 盒+FLP 位点） | Baba et al. 2006, *Mol Syst Biol* |
| **DEG 15.0** | 78 细菌 + 35 真核 + 2 古菌 | 细菌平均 ~300-600/物种 | 各物种不等 | 全基因组 | 单基因敲除/转座子诱变/RNAi/CRISPR-Cas9 | Zhang R, *NAR* 持续更新 |
| **OGEE v3** | 91 物种（含 581 人类细胞系） | 细菌 21,914（v2） | 细菌 78,075（v2） | v3 接近翻倍 | 实验+文本挖掘；条件性必要标记 | Gurumayum et al. 2021, *NAR* |
| **BacDive 生理数据** | 20,433+ 物种/100,866 菌株 | — | — | 1,721,500+ 生理条目 | 实验测定+文献整合 | DSMZ, 2025 版 |

### 1.2 预测准确性指标

| 指标 | 定义 | 公式 | 适用场景 | 特点 |
|------|------|------|---------|------|
| **AUROC** | TPR vs FPR 曲线下面积 | TPR vs FPR 积分 | 类别平衡时 | 对不平衡数据偏乐观 |
| **AUCPR** | 精确率 vs 召回率曲线下面积 | Precision vs Recall 积分 | **必要性数据强烈推荐** | 对不平衡数据更准确（必要基因 ~5-15%） |
| **Accuracy** | 总体正确率 | (TP+TN)/(TP+TN+FP+FN) | 基础比较 | 可能因高 TN 而虚高 |
| **F1** | 精确率与召回率的调和平均 | 2×P×R/(P+R) | 平衡 P 和 R | 对假阳性/假阴性等权 |
| **MCC** | Matthews 相关系数 | (TP×TN-FP×FN)/sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN)) | **金标准** | 唯一对 4 个象限等权，-1 到 +1 |

### 1.3 各工具基因必要性性能

| 工具 | 数据集 | AUCPR | Accuracy | F1 | 备注 |
|------|--------|-------|----------|-----|------|
| **GEMsembler Core3+GA** | E. coli Keio | 0.731 | 0.80 | — | GA 优化 50 代，15 碳源 |
| **GEMsembler Core3** | E. coli Keio | 0.556 | — | — | 无优化 |
| **gapseq** | E. coli DEG | — | 0.80 | — | 自动重建模型 |
| **ModelSEED+GrowMatch** | 多种 | — | ~0.87 | — | GrowMatch 优化后（66%→87%） |
| **iML1515a+GA** | E. coli Keio | 0.771 | — | — | 从 0.754 改进 |
| **CarveMe** | E. coli | — | — | — | 基因 ID 命名空间不兼容（cds_XXXXX vs b-number），实验数据不适用 |

### 1.4 预测阈值

- 必要基因：敲除后生长率 < 1% 野生型（growth < 0.01 × WT）
- 非必要基因：敲除后生长率 ≥ 1% 野生型

---

## 2. 生长表型预测评估

### 2.1 实验数据集

| 数据集 | 规模 | 内容 | 评估指标 |
|--------|------|------|---------|
| **Biolog PM1-PM10** | ~950 条件/板 | 190 碳源（PM1-2）+ 95 氮源（PM3）+ 59 磷源（PM4）+ 95 硫源 | Accuracy, TPR, FPR, No-growth accuracy |
| **BacDive 碳源** | 4,349 菌株 × 最多 49 碳源 | 碳源利用表型（+/−/弱） | 与 FBA 模拟对比 |
| **gapseq 表型基准** | 14,931 调查项（48 碳源 × ~310 物种） | 从文献/BacDive/ProTraits 整合 | gapseq vs CarveMe vs ModelSEED 三工具对比 |
| **AGORA2 验证集** | 3 个独立实验数据集 | 摄取/分泌代谢物（>99% 一致性） | Sensitivity, Specificity |

### 2.2 各工具生长表型性能

| 工具 | 碳源 TP 率 | 碳源假阴性率 | 碳源 Accuracy | 数据来源 |
|------|----------|------------|-------------|---------|
| **gapseq** | 47% | 6%（酶活性） | 80% | Hsieh 2024 + Zimmermann 2021 |
| **Bactabolize** | 85.79%（C 源） | — | 97% 临床株 | KpSC 内 10 菌株验证 |
| **CarveMe** | 24% | 32%（酶活性） | 66% | Hsieh 2024 |
| **ModelSEED/KBase** | 31% | 28%（酶活性） | — | Hsieh 2024 |
| **AGORA2** | — | — | >99%（摄取/分泌） | 3 独立验证集 |
| **gempipe** | — | — | 92.3%（代谢预测） | mSystems 2025 基准 |

### 2.3 FBA 模拟条件

- 碳源利用：设定唯一碳源交换反应 lower bound = -10 mmol/gDW/h，其余碳源关闭
- 生长阈值：growth > 1e-6 mmol/gDW/h → 预测生长
- 常用碳源测试集（E. coli）：20 种正例（glucose, fructose, galactose, glycerol, succinate...）+ 7 种负例（methanol, sucrose, lactose...）

---

## 3. MEMOTE 测试套件

### 3.1 测试类别

| 测试类别 | 测试数 | 权重 | 测试内容 | SBO 术语检查 |
|---------|--------|------|---------|------------|
| **Consistency** | 39 | 3× | 质量/电荷平衡、化学式一致性、反应循环检测 | — |
| **SBO Terms** | 15 | 2× | 反应/代谢物/基因的 SBO 术语标注 | SBO:0000627(交换)、SBO:0000176(生化)、SBO:0000185(转运)、SBO:0000247(代谢物)、SBO:0000243(基因) |
| **Metabolite Annotation** | 18 | 2× | InChI/SMILES/CHEBI/公式存在性 | — |
| **Reaction Annotation** | 7 | 2× | EC 编号/GPR 规则存在性 | — |
| **Gene Annotation** | 2 | 2× | 基因产物名称、标准注释 | — |
| **Biomass** | 7 | 1× | 生物量反应元素检查 | — |
| **Basic Tests** | 5 | 1× | 模型可解性、代谢物/反应数量 | — |
| **总计** | **93** | — | — | — |

### 3.2 评分系统

- 基础分数 = (passed × weight) / ((passed + failed) × weight)
- 新模型评分均值 ~50-60%
- 高质量模型 ≥ 80%
- CarveMe 在 E. coli 上达 85.2%（在所有基因组中最高）

### 3.3 各工具 MEMOTE 分数

| 基因组 | CarveMe | gapseq | GEM-Agent | Bactabolize |
|--------|---------|--------|-----------|------------|
| E. coli K-12 | 85.2% | — | 85.2% | — |
| P. aeruginosa PAO1 | 68.7% | — | 68.7% | — |
| B. subtilis 168 | 70.9% | — | 70.9% | — |
| S. aureus 8325 | 69.3% | — | 69.3% | — |
| K. pneumoniae | — | — | — | 高（含 MEMOTE 报告输出） |

---

## 4. 生物量一致性检查

### 4.1 生物量组分标准比例

| 大分子 | E. coli（%） | B. subtilis（%） | K. pneumoniae（%） |
|--------|------------|----------------|-------------------|
| Protein | 55.0 | 45.0 | 53.0 |
| RNA | 20.5 | 18.0 | 18.9 |
| DNA | 3.1 | 2.0 | 3.4 |
| Lipid | 9.1 | 9.0 | 9.8 |
| Cell Wall | 2.5（LPS/肽聚糖） | 5.0（teichoic acid） | 3.2 |
| Cofactors | 2.8 | 2.5 | 3.0 |
| Ion | 1.0 | 1.0 | 1.0 |

### 4.2 检查项

- 质量平衡：每个反应 + 生物量反应的 H/C/N/O/P/S 元素平衡
- 电荷平衡：每个反应净电荷为 0
- ATP 产量检查（GAM/NGAM）：GAM = Growth-Associated Maintenance, NGAM = Non-Growth-Associated Maintenance
- 生物量组分生产检查：每个生物量前体都必须有合成路径

### 4.3 各工具的生物量方程策略

| 工具 | 策略 |
|------|------|
| **CarveMe** | BiGG 通用模板的预定义生物量反应，含菌种特异裁剪 |
| **gapseq** | 基于 16S 预测 Gram 染色状态 → Gram+/Gram- 预设生物量方程 |
| **ModelSEED/KBase** | 含 DNA/RNA/蛋白/脂质/细胞壁/辅因子的通用方程 |
| **AGORA2/DEMETER** | Gram 状态修复（修复 32% 草稿错误）+ 菌种特异调整 |
| **pyFBA** | 多种生物类型的预定义方程（gramneg/grampos/mycobacterium/plant） |

---

## 5. 多工具共识指标

### 5.1 Hsieh et al. 2024 基准

| 基因组 | CM rxns | GS rxns | Shared BiGG | Jaccard | GPR 覆盖（CM/GS） |
|--------|---------|---------|------------|---------|-----------------|
| E. coli K-12 | 2,364 | 3,012 | 682 | 0.183 | 58%/79% |
| P. aeruginosa PAO1 | 2,471 | 2,775 | 586 | 0.152 | — |
| B. subtilis 168 | 1,141 | 2,758 | 444 | 0.166 | — |
| S. aureus 8325 | 1,267 | 2,265 | 356 | 0.135 | — |

关键发现：仅 ~25% 特征被 4 工具共同支持；~50% 仅被单一工具支持；Biomass 成分 Jaccard 仅 0.23-0.24

### 5.2 GEMsembler CoreX 系统

| 共识水平 | 定义 | 特征可信度 | 典型应用 |
|---------|------|----------|---------|
| **Core4** | 4 工具一致 | 最高（4/4） | 核心代谢功能 |
| **Core3** | ≥3 工具一致 | 高（3/4） | 通路分析与 curation |
| **Core2** | ≥2 工具一致 | 中（2/4） | 扩展现有知识 |
| **Core1** | ≥1 工具 | 低（1/4） | Assembly/Union |

### 5.3 GPR 优化算法对比

| 算法 | 修改 GPR 数 | 添加基因数 | AUCPR 改进 | 计算复杂度 |
|------|----------|----------|----------|----------|
| **逐步算法 (SA)** | 3-34 | 1-10 | +13.5%（core3: 0.556→0.691） | 低 |
| **遗传算法 (GA)** | 数百 | 6-248 | +15.5%（core3: 0.556→0.711） | 高（200×50 代 ×15 碳源） |

---

## 6. 实验数据交叉验证

### 6.1 验证数据集推荐

| 验证类型 | 推荐数据集 | 规模 | 适用物种 | 指标 |
|---------|----------|------|---------|------|
| **基因必要性** | Keio Collection | 4,288 基因 | E. coli K-12 | AUCPR, F1, MCC |
| **碳源利用** | BacDive（碳源表型） | 4,349 菌株 × 49 碳源 | 多种 | 与 FBA 预测比较 |
| **酶活性** | BacDive（酶活性表型） | 1,721,500+ 生理条目 | 多种 | TP/FP 率 |
| **生长表型** | gapseq 综合基准 | 14,931 调查项 | ~310 物种 | Accuracy, TP/FP |
| **摄取/分泌** | AGORA2 验证集 | 3 独立实验数据集 | 人类肠道微生物 | Sensitivity, Specificity |

### 6.2 BacDive 数据使用说明

- 碳源利用：+/-/weak → 用于验证 FBA 预测的碳源生长能力
- 酶活性：含 EC 编号的酶活性检测 → 验证对应蛋白功能预测
- 发酵产物：化合物列表 → 验证模型分泌反应
- 氧气需求：aerobic/anaerobic/facultative → 验证电子传递链设置
- 温度/pH 范围：min/opt/max → 用于 gap-filling 条件

---

## 7. 各工具自我评估方法

### 7.1 CarveMe
- 14,931 碳源利用查证
- 酶活性验证：假阴性 32%（偏高）
- Accuracy 0.66（低于 gapseq 0.80）
- 死胡同代谢物最少
- 无专门基因必要性验证（基因 ID 命名空间不兼容）

### 7.2 ModelSEED / KBase
- **GrowMatch**：基因必要性准确率 66% → 87%（MILP 搜索最优修改方案）
- FBA 一致性检验
- ModelSEED v2 基准：增强的 ATP 生物合成通路预测

### 7.3 gapseq
- 酶活性检验：10,538 测试项，假阴性 6%（三工具最优）
- 碳源利用：48 碳源 × ~310 物种
- 碳源真阳性率 47%（三工具最优，对比 CarveMe 24%, ModelSEED 31%）
- 基因必要性 Accuracy 0.80

### 7.4 AGORA2 / DEMETER
- **全面测试套件**：合格标准 100% 厌氧生长、真实 ATP 产量
- 74/74 定义培养基正确生长预测
- 3 独立验证集：摄取/分泌一致性 >99%
- 预测药物代谢准确率 0.81
- 精炼前：27% 草稿模型可厌氧生长；精炼后：100%

### 7.5 GEMsembler
- **唯一使用 AUCPR 作为优化目标**
- Core3 在 L. plantarum auxotrophy 预测上超越人工 gold-standard iLP728
- 基因必要性 AUCPR：core3+GA 0.711（E. coli）
- 发现并修正 gold-standard 错误（DBTS: b0778 → b0778 OR b1593）

### 7.6 gempipe
- 基因恢复完整度：97.8%
- 参考覆盖率：91.2%
- 代谢预测准确率：92.3%
- 孤儿反应率：4.1%（最低）
- 额外验证：Biolog 碳源利用在 K. pneumoniae 37 菌株上

### 7.7 pyFBA
- 定性验证（生长/不生长）——适合教学项目
- 二分法反应缩减（bisection）：验证最小反应集是否正确

### 7.8 Bactabolize
- 碳/氮/磷/硫源利用预测（C 源 85.79%, N 源 84.09%, P 源 89.47%）
- 基因必要性：F1 0.48（KpSC 10 株临床株）
- 参考模型基因/反应捕获率 >99%
- MEMOTE 报告自动生成

---

## 总结：推荐的评估矩阵

| 评估层级 | 方法 | 数据集 | 指标 | 适用场景 |
|---------|------|------|------|---------|
| **L0: 语法** | libsbml 解析 | — | 解析成功/失败 | 全部 |
| **L1: 结构** | MEMOTE | — | Score 0-100% | 全部 |
| **L2: 功能** | FBA 生长 | — | 生长率 > 0 | 全部 |
| **L3: 碳源** | FBA + BacDive | 4,349 菌株 × 49 碳源 | Accuracy, TP/FP/TN/FN | 有 BacDive 数据 |
| **L4: 必要性** | FBA 基因敲除 | Keio/DEG/OGEE | AUCPR, F1, MCC | 有同源必要数据 |
| **L5: 代谢物** | 摄取/分泌 | AGORA2 验证集 | Sensitivity, Specificity | 肠道微生物 |
| **L6: 共识** | 多工具 Jaccard | 同基因组跑多工具 | Jaccard 相似度 | 跨工具对比 |

## 参考文献

1. Baba T et al. (2006). *Mol Syst Biol* 2:2006.0008. [Keio Collection]
2. Zhang R et al. (2021). *Nucleic Acids Research*. [DEG 15.0]
3. Gurumayum S et al. (2021). *Nucleic Acids Research* 49(D1):D998-D1006. [OGEE v3]
4. Lieven C et al. (2020). *Nature Biotechnology* 38:272-276. [MEMOTE]
5. Hsieh YE et al. (2024). *npj Systems Biology and Applications* 10:54. [三工具对比]
6. Matveishina EK et al. (2025). *mSystems* 10(10):e0057425. [GEMsembler]
7. Zimmermann J et al. (2021). *Genome Biology* 22:81. [gapseq]
8. Heinken A et al. (2023). *Nature Biotechnology* 41:1320-1331. [AGORA2/DEMETER]
9. Vezina B, Watts SC et al. (2023). *eLife* 12:RP87406. [Bactabolize]
10. Lazzari G et al. (2025). *mSystems* 11(1):e0100725. [gempipe]
11. Machado D et al. (2018). *Nucleic Acids Research* 46(15):7542-7553. [CarveMe]
12. Henry CS et al. (2010). *Nature Biotechnology* 28(9):977-982. [ModelSEED]
13. Schober I et al. (2025). *Nucleic Acids Research*. [BacDive 2025]
