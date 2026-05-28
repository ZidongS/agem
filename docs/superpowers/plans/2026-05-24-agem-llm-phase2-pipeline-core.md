# agem_llm Phase 2: GPR 全量构建 + SBML 生成 + 拓扑评估 + LLM 精修

**Goal:** 实现完整 pipeline：蛋白质预过滤 → 多源 GPR 汇编 → SBML 构建 → 拓扑评估 + 模板对比 → GPR 元数据库 → LLM 精修。

**Architecture:** 确定性 GPR 构建 + SBML 生成在前，LLM 仅在精修阶段做决策。

**Tech Stack:** Python 3.10+, asyncio, cobrapy, libsbml, sqlite3, litellm

---

## Task 1: GPR 全量构建器
- Create: agem_llm/agem_llm/gpr_builder.py
- build_full_gpr(): 对每个 CDS 查询 KEGG+BRENDA+Rhea+ModelSEED2+EZSpecificity，合并多维证据

## Task 2: SBML 构建器（GPR JSON → SBML）
- Create: agem_llm/agem_llm/sbml_builder.py
- build_sbml_from_gpr(): Rhea equation 解析 stoichiometry → COBRApy → SBML L3V1 FBCv3
- 添加生物量方程（Gram 染色状态驱动）、交换反应（BacDive 碳源）、转运反应

## Task 3: 拓扑连通性评估
- Create: agem_llm/agem_llm/topology.py
- find_deadends(): 死胡同代谢物检测
- check_biomass_reachability(): 生物量前体可达性
- find_isolated_subnetworks(): 断连子网检测
- assess_topology(): 完整拓扑评估

## Task 4: CarveMe 通用模板对比
- Create: agem_llm/agem_llm/template_compare.py
- compare_with_template(): 加载 gramneg/grampos 模板，求差集，回溯 GPR 证据，按通路分组

## Task 5: GPR 元数据检索引擎
- Create: agem_llm/agem_llm/gpr_store.py
- GPRStore: SQLite-backed 元数据索引
- search(where_clause), search_by_category(), get_stats()

## Task 6: LLM 精修 Agent
- Create: agem_llm/agem_llm/agents/refinement_agent.py
- RefinementAgent: tool calling 检索 GPR store → 诊断 → 提案修复 → libsbml 执行
- 三类修复：模板找缺补全、拓扑死胡同修复、碳源表型修复

## Task 7: 完整 Pipeline 组装
- Create: agem_llm/agem_llm/pipeline.py
- Pipeline.run(): 7 阶段完整 pipeline + 审计日志
- 集成测试
