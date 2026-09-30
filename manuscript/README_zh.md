# 研究说明（中文）：升主动脉瘤（ATAA）与急性 A 型主动脉夹层（TAAD）的分子机制差异

**目标期刊**：Journal of Translational Medicine（Springer Nature/BMC，JCR Q1，IF≈6）。全文英文稿见 `manuscript.md` 与 `manuscript_JTM.docx`，队列汇总表见 `Table1_cohorts.csv`。

## 一、研究设计

- **数据来源**：2026-09-29 系统检索 GEO，纳入所有人类升/胸主动脉组织"疾病 vs 非病变主动脉"的转录组。
  - ATAA：4 个 bulk 队列，n=115（GSE26155、GSE140947、GSE235161、GSE202267）
  - TAAD：6 个 bulk 队列，n=76（GSE52093、GSE153434、GSE98770、GSE294606、GSE267434、GSE147026）
  - 剔除 GSE190635：其差异表达谱与其他所有夹层队列呈强负相关（ρ −0.37～−0.68），提示 GEO 上传的分组标签可能颠倒。
  - 独立验证：GSE318877（同一平台内直接比较夹层 vs 扩张的中膜显微区域 RNA-seq，70 个样本）；单细胞 GSE155468（ATAA）+ GSE213740（TAAD），整合后 73,755 个细胞。
- **统一流程**：基因符号统一 → 每队列差异分析（芯片/FPKM 用 limma-trend，counts 用 DESeq2）→ Stouffer 加权 Z + DerSimonian-Laird 随机效应 meta 分析 → 留一队列（LOCO）稳健性 → 两病基因分类（共享/夹层特异/动脉瘤特异/方向相反，附 z_diff 差异检验）→ GSEA（Hallmark/KEGG/Reactome）与 ORA → 单细胞衍生标志物 ssGSEA 推断细胞组成 → STRING PPI 枢纽 → 对照参照 z 分数 + LOCO 弹性网分类器。
- **关键处理修正**：GSE98770 从原始 Agilent 文件重处理（GEO 矩阵无信号，重处理后 1323 个 DEG，MYOCD/CXCL14/PI16 下调）；GSE235161 限定蛋白编码基因并重新标准化（rRNA/7SL 主导 FPKM）。

## 二、主要发现

1. **可重复性不对称**：夹层签名在 6 个队列/4 种平台间高度一致（队列间 ρ 中位数 0.50；LOCO 保留 73–85%），动脉瘤签名弱且依赖队列（ρ 中位数 0.12；去掉 GSE26155 仅保留 9%）。提示单队列动脉瘤"枢纽基因"结论需谨慎。
2. **两病分子上大部分不同**：全基因组 meta log2FC 相关仅 ρ=0.11；共识 DEG 重叠 137 个（Jaccard 0.056）。
3. **TAAD 特异（1780 基因）**：SMC 收缩程序丢失（MYOCD −1.60、ACTA2、CNN1、LMOD1、MYH11）；内皮（VWF、PECAM1、PTPRB、TEK）与外膜成纤维（PI16、DCN、APOD、CXCL14）转录本丢失；缺氧/应激（HIF1A、VEGFA、HMOX1、MT2A、CDKN1A）；MYC/mTORC1/E2F、核糖体生物合成、DNA 修复、ROS、UPR 上调；固有炎症（IL6、CXCL8、CCL2、S100A8/A9、SPP1、SERPINE1、TIMP1、MMP14）。PPI 枢纽：AKT1、STAT3、MYC。
4. **ATAA 特异/富集（230 基因）**：I 型干扰素（OAS1-3、IFIT1-3、MX1/2、CXCL10）、MHC-II 抗原呈递（CD74、HLA-DRA/DPA1/DPB1、RFX5 调控子）、粘附与 ECM（VCAM1、ITGA4、COL1A2、MXRA5）；无 SMC 收缩基因丢失；生物钟基因（BMAL1、NPAS2、PER3、CRY2）下调为新发现。
5. **共享核心（107 基因 / 15 个 Hallmark）**：TNFα-NF-κB、IL6-JAK-STAT3、缺氧、EMT、TGF-β、补体、血管生成（VEGFA、APLN）、UPR；共同下调 SMC 钙处理基因（RYR2、CASQ1、DES、ACTN2）。MMP2/MMP9 在两病均无稳定改变。
6. **细胞组成**：TAAD = SMC/成纤维/内皮评分下降 + S100A8/9⁺ 炎性髓系浸润；ATAA = T/NK/浆细胞浸润、无 SMC 丢失；单细胞比例方向一致（ATAA 中 T 细胞 7.9%→46.4%，P=0.024；TAAD 中 S100A8/9⁺ 髓系 4.2%→18.6%）。
7. **验证**：直接对比队列通路层面 ρ=0.42（P=5e-48），夹层特异基因 60% 方向一致（P=1e-28），夹层特异评分在 4 个中膜区域均升高（P=5e-6）；LOCO 分类器汇总 AUC=0.92。动脉瘤特异的干扰素程序在直接对比队列中未被复现（夹层中膜干扰素活性反而更高），已在文中如实说明。

## 三、你需要补充/决定的内容

- 作者、单位、通讯作者、基金、作者贡献（稿中留有占位符）。
- 部分数据集无正式发表文献（GSE235161、GSE202267、GSE294606、GSE267434、GSE147026、GSE190635、GSE213740、GSE318877），目前按 GEO 编号引用；投稿前请核对是否已有对应论文并补充引用。
- 参考文献按 JTM（Vancouver 顺序编号）格式编排，建议投稿前用 EndNote/Zotero 复核每条文献的卷期页码。
- 建议增加湿实验验证以提高录用概率（JTM 审稿人常要求）：在 ATAA 与 TAAD 手术标本中用 qPCR/免疫组化验证 MYOCD、HMOX1、S100A8/A9、IL6（夹层特异）与 CXCL10、CD74/HLA-DR、OAS1（动脉瘤富集）的表达；或用血浆 IL-6/S100A8/A9/SPP1 做初步循环标志物验证。

## 四、文件位置

- 图：`figures/Fig1_study_design.png`、`Fig1b_cohort_concordance.png`、`Fig2_gene_level_comparison.png`、`Fig3_hallmark_NES_heatmap.png`、`Fig3b_pathway_scatter.png`、`Fig4_cell_composition.png`、`Fig5_signature_ml.png`、`Fig6_PPI_hubs.png`、`Fig7_scRNA_overview.png`、`Fig8_sc_validation.png`（均有 PDF 版本）；QC 图在 `figures/qc/`。
- 结果表：`results/meta/`（meta 分析、基因分类、LOCO、直接验证）、`results/pathways/`（GSEA/ORA）、`results/cellcomp/`、`results/ml/`、`results/ppi/`、`results/sc/`、`results/sc_validation/`。
- 代码：`scripts/00`–`13`，可按 README.md 顺序完整复现。
