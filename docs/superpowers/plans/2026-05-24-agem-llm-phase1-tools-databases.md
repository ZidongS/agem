# agem_llm Phase 1: 工具层 + 数据库层实现计划

**Goal:** 搭建 agem_llm 项目骨架，搬运 gem_agent 工具层，构建 Rhea/KEGG/BRENDA/ModelSEED2 数据库客户端，实现 BioRegistry 统一 ID 转换。

**Architecture:** 纯确定性代码层。每个数据库封装为独立 Python 模块，统一返回 {success, data, error} 格式。

**Tech Stack:** Python 3.10+, asyncio, SQLite, requests, bioregistry, cobrapy, pyrodigal, litellm, libsbml, memote, pydantic, pyyaml

---

## Task 1: 项目骨架
- Create: agem_llm/__init__.py, agem_llm/config.yaml, agem_llm/schemas.py, agem_llm/cache_manager.py
- Create: setup.py, requirements.txt, .gitignore

## Task 2: 搬运 gem_agent 工具层
- 拷贝 pyrodigal.py, kofamscan.py, eggnog.py, mmseqs.py, cobrapy_utils.py, sbml_ops.py, memote.py, base.py
- 替换 import 前缀 gem_agent → agem_llm

## Task 3: BioRegistry 统一 ID 转换层
- Create: agem_llm/databases/bioregistry_client.py
- normalize_prefix(), expand_curie(), cross_map()

## Task 4: Rhea REST 客户端 + 本地 SQLite 缓存
- Create: agem_llm/databases/rhea_client.py, rhea_local.py, data/build_rhea_db.py
- query_rhea_by_ec(), get_rhea_reaction(), query_rhea_by_uniprot()
- RheaLocalDB: build_from_rhea_api(), get_reactions_for_ec()

## Task 5: KEGG 本地数据库
- Create: agem_llm/databases/kegg_local.py, data/build_kegg_db.py
- KEGGLocalDB: build_from_kegg_api(), get_reactions_for_ko(), get_reactions_for_ec()

## Task 6: BRENDA 客户端
- Create: agem_llm/databases/brenda_client.py
- search_ec_in_brenda(): 查 EC 在某物种/近缘物种是否有实验支持
- 支持 JSON 批量下载 + SOAP API 实时查询

## Task 7: ModelSEED2 客户端
- Create: agem_llm/databases/modelseed_client.py
- ModelSEEDClient: search_reactions_by_ec(), download_from_github()
- 解析 compounds.tsv + reactions.tsv

## Task 8: BacDive + ProTraits/细菌-古菌性状 客户端
- Create: agem_llm/databases/bacdive_client.py, traits_client.py
- search_bacdive(), extract_metabolic_traits()
- TraitsClient: search_by_taxon_id(), get_traits_summary()

## Task 9: EZSpecificity 推理接口（骨架）
- Create: agem_llm/databases/ezspecificity_client.py
- predict_substrate_specificity(), batch_predict()
- GPU 部署后填充实际推理
