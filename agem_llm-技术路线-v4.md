# agem_llm 技术路线 v4

> 核心变更：向量库简化为元数据索引；新增拓扑连通性评估；新增 CarveMe 通用模板对比找缺；综合 CarveMe + gapseq 优势

---

## 系统架构（7阶段）

1. **蛋白质预过滤** — CLEAN-ESM2（酶分类）+ GenEfflux（转运体分类）
2. **多工具深度注释** — KofamScan(KO) + eggNOG(EC/GO) + MMseqs2(UniProt) + EZSpecificity(底物) + BioRegistry 统一 ID
3. **跨数据库 GPR 全量构建** — EC/KO → KEGG + BRENDA + Rhea + ModelSEED2 并行查询；UniProt 强匹配→Rhea 精确反应；EZSpecificity 打分
4. **全 GPR 直接构建 GEM** — Rhea stoichiometry 优先 → SBML L3V1 FBCv3（含生物量方程+转运+交换反应）
5. **拓扑连通性 + CarveMe 模板对比找缺**
   - 5a: 生物量前体可达性、死胡同分析、断连子网检测、hub 代谢物检查
   - 5b: 加载 gramneg/grampos 模板，求差集，回溯 GPR 证据
   - 5c: 综合 gapseq 级 GPR 覆盖率 + CarveMe 级拓扑完整性
6. **GPR 元数据库** — 纯元数据索引（CDS, EC, reaction, pathway, evidence, topology_status, template_missing），多维条件查询
7. **LLM 精修** — 模板找缺补全 + 拓扑修复 + 碳源表型修复，tool calling 检索 GPR 库，libsbml 执行修复

## 数据库角色分工

| 数据库 | 角色 | 提供信息 |
|--------|------|---------|
| Rhea | 主反应库 | reaction stoichiometry, ChEBI 化合物, EC×反应映射（REST API, 17,783 反应） |
| BRENDA | 有机体特异性验证 | 某 EC 在某物种/近缘物种中的实验证据 |
| EZSpecificity | CDS→底物特异性 | 给定 CDS 序列 + 候选底物 → 催化概率（Nature 2025） |
| BacDive | 实验表型（输入+QC） | 碳源利用, 发酵产物, 需氧性（Python REST 客户端） |
| ProTraits/细菌-古菌性状 | 预测表型（输入+QC） | 420+ 性状, 孢子形成, 运动性, Gram 染色 |
| KEGG | KO/EC→反应 | 本地预构建 SQLite |
| ModelSEED2 | 补充反应来源 | TSV 解析 EC→reaction |
| BioRegistry | 统一 ID 转换 | ChEBI⇄KEGG⇄BiGG⇄Rhea |

## 综合 CarveMe + gapseq 优势

| 质量维度 | CarveMe | gapseq | agem_llm 策略 |
|---------|---------|--------|-------------|
| GPR 覆盖率 | 58-78% | **79-86%** | 全量多数据库 GPR 汇编 |
| 死胡同代谢物 | **最低** | 高 | 连通性分析 + LLM 修复 |
| 拓扑完整性 | **高** | 中 | 模板对比 + LLM 补全 |
| 转运蛋白 | 无 | **专门模块** | GenEfflux + Rhea transport |
| 酶底物特异性 | DIAMOND | blast | **EZSpecificity (2025)** |
| 实验验证 | 无 | 无 | **BacDive + BRENDA** |

## 与其他工具对比

| 维度 | CarveMe | gapseq | GEMsembler | agem_llm |
|------|---------|--------|-----------|----------|
| 反应来源 | BiGG 模板 | 同源+通路 | 组合 4 工具 | Rhea+KEGG+BRENDA+ModelSEED2 |
| GPR 构建 | DIAMOND→评分 | bitscore 权重 | 继承各工具 | 多数据库 + EZSpec 打分 |
| 模板对比 | — | — | — | **CarveMe 通用模板对比** |
| 拓扑评估 | MILP 隐含 | gapfill 隐含 | 通路分析 | **显式连通性评估** |
| 决策引擎 | MILP | 加权 | 多数投票 | **LLM（后置精修）** |
| 可解释性 | 低 | 中 | 中 | **高（每 GPR 有完整证据链）** |

## 实施顺序

| Phase | 内容 |
|------|------|
| Phase 1 | 搬运 gem_agent 工具层 + BaseAgent |
| Phase 2 | CLEAN + GenEfflux 预过滤 |
| Phase 3a-f | 数据库客户端（Rhea, BRENDA, KEGG, ModelSEED2, BioRegistry, EZSpecificity） |
| Phase 4 | 多源 GPR 汇编 + SBML 构建 |
| Phase 5a-c | 拓扑评估 + 模板对比 + 综合评估 |
| Phase 6 | GPR 元数据库 + 检索引擎 |
| Phase 7 | LLM 精修 Agent（tool calling） |
| Phase 8 | Benchmark（vs CarveMe/gapseq） |
