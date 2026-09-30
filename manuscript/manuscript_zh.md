# 升主动脉瘤与 A 型主动脉夹层中的收缩程序与损伤应答程序：人主动脉转录组的队列水平与患者水平比较

**Journal of Translational Medicine – 研究论文（投稿草稿，中文版）**

作者：[第一作者]^1^，[合作者]^1,2^，[通讯作者]^1*^

^1^ [科室，单位，城市，国家]；^2^ [科室，单位，城市，国家]

\* 通讯作者：[姓名，电子邮件]

---

## 摘要

**背景。** 手术切除主动脉组织的转录组由壁细胞功能表型、组织组成以及组织的应激与修复反应共同构成。因此，升主动脉瘤（ATAA）与急性 A 型主动脉夹层（TAAD）的差异既可能是同一组程序的强度差异，也可能是收缩、重塑与损伤应答程序的不同组合。仅在队列间比较疾病标签无法区分这两种可能，而队列内的患者级结构可以。

**方法。** 对 10 个 GEO bulk 队列（ATAA 4 个队列，63 例患者对 52 例对照；TAAD 6 个队列，41 对 35）用同一流程重新处理并做随机效应 meta 分析。在每个数据集中以与标签无关的方式计算两个外部定义的程序分：平滑肌收缩分（C）和损伤应答综合分（I；缺氧、糖酵解、氧化应激/NFE2L2、MYC-核糖体、p53-DNA 损伤与 NF-κB/IL-6 六个模块），以及二者之和 R = I + C 作为相对状态量。分析计划在探索阶段之后锁定：疾病间差异的一个队列级检验家族（C、I、R；Holm 校正）；每个单细胞数据集一个精确的患者级检验（收缩型平滑肌细胞内的 R；GSE155468、GSE213740、GSE189795）；一个患者内区域检验（瘤体对瘤颈，GSE140947）；病程（GSE222318）和中膜分层（GSE318877）数据集仅报告估计；一个队列内有序扩张程度的探索性趋势分析（GSE26155）。

**结果。** 收缩程序的下降在两病之间存在差异（合并 Hedges' g：TAAD −2.11，ATAA −0.34；差值 −1.77，Holm P < 0.001）；损伤应答综合分在两病中均升高、差异仅为临界（+1.75 对 +0.56；Holm P = 0.10）；相对状态 R 在两病之间无差异（+0.20 对 +0.24；P = 0.92）。糖酵解、氧化应激和 MYC-核糖体模块在 TAAD 中更高（BH ≤ 0.02），且在 TAAD 中这些效应在条件化三个细胞签名分后仍保留。干扰素-α、MHC II 类、T 细胞和 ECM-胶原分在 ATAA 中升高，但与 TAAD 的差异不显著。在收缩型平滑肌细胞内，两病的 C 均低于对照（g −0.95 与 −0.98；区间不含 0），I 与 R 无明确变化，R 的患者级主要检验均未通过（精确 P = 0.37 与 0.79）。[[STAGE_ABSTRACT_ZH]] 在动脉瘤患者内，瘤体与瘤颈的 R 无差异（P = 0.19）。在一个队列的有序扩张程度中，C 随扩张下降（ρ = −0.35，BH = 0.016），干扰素、MHC II 类、淋巴细胞和 ECM 分随扩张升高（BH ≤ 0.03），I 无趋势。嵌套留一队列分类器可区分两病（合并 AUC 0.90）。

**结论。** 在组织水平，TAAD 与 ATAA 的差异主要在于收缩程序丢失与代谢/应激激活的幅度，而两类程序的相对状态并无差异；"两病以不同方式组织这些程序"的假设未获现有患者级数据支持，且这些数据的检验效能不足。动脉瘤扩张伴随收缩程序下降以及免疫与基质程序，但不伴随夹层所特有的损伤应答综合分。

**关键词：** 升主动脉瘤；A 型主动脉夹层；转录组；meta 分析；单细胞 RNA 测序；平滑肌细胞；损伤应答；干扰素

---

## 背景

升主动脉瘤（ATAA）与急性 Stanford A 型主动脉夹层（TAAD）在临床上被作为同一退行性过程的不同阶段处理，然而多数 A 型夹层发生在低于手术直径阈值的主动脉，而许多动脉瘤可数十年保持稳定 [1–3]。两者均表现为中膜退变，伴弹力纤维断裂、蛋白聚糖沉积和平滑肌细胞（SMC）丢失 [4,5]；夹层另有内膜撕裂、壁内血肿和固有免疫浸润 [6–9]，动脉瘤组织则以淋巴细胞为主的浸润和增生性重塑为特征 [10,11]。共同的遗传易感性涉及 SMC 收缩、细胞外基质（ECM）和 TGF-β 基因 [12–14]，但不能解释为何一条主动脉扩张而另一条夹层。

手术标本的转录组是一个复合信号。它反映常驻壁细胞的功能表型（首先是 SMC 的收缩状态）、标本的细胞组成（包括浸润细胞），以及手术时活跃的应激与修复反应（夹层中包括缺氧、氧化与代谢应激和急性炎症）。当通过疾病标签比较两种疾病时，差异可来自上述任一成分及其组合。同一个"疾病对对照"效应差异与两种截然不同的情形相容：两病可能以不同强度表达同一组程序，也可能以不同方式组合收缩、重塑与损伤应答程序。要区分二者，需要疾病标签之下的结构：在同一类壁细胞内、在同一患者的不同部位、在同一疾病的不同病程阶段以及在同一队列的扩张程度范围内测量同一组程序。

近期两项研究为该问题提供了框架。一项 80 例人升主动脉的单细胞图谱在夹层中识别出糖酵解 ENO1/MIF 轴和 FBN1 阳性成纤维细胞的耗竭，并含干预实验 [15]；一项动脉瘤与夹层细胞异质性比较被概括为"侵蚀或爆炸" [16]。二者均强调两病不同，但都未问差异是强度上的还是组织方式上的。主动脉转录组的 meta 分析迄今只考察单一队列或单一疾病，已发表主动脉签名的跨队列可重复性尚不清楚 [17,18]。

本研究重新处理了所有符合条件的人升/胸主动脉组织公共 bulk 转录组，以随机效应 meta 分析合并，然后用外部定义的程序分提出一个更窄的问题：ATAA 与 TAAD 在 SMC 收缩程序与损伤应答程序的相对状态上是否不同，该关系是否在壁细胞内、随夹层病程、在患者内跨区域与中膜分层而变化？分析计划（包括唯一主统计量和每个患者级数据集一个推断检验）在探索阶段之后、患者级分析运行之前锁定。

## 方法

### 数据来源与纳入标准

于 2026 年 9 月 29 日在 GEO DataSets 以 "aortic dissection"、"ascending aortic aneurysm"、"thoracic aortic aneurysm" 与 "ascending aorta AND aneurysm" 检索人表达谱系列（芯片或高通量测序）（附加文件 1）。队列级分析的纳入条件：对散发性 ATAA 或急性 TAAD 患者的 bulk 升/胸主动脉组织与非病变主动脉组织（器官供体、心脏移植供体或心脏手术的非扩张主动脉）进行比较，并提供处理后的表达数据或可重新处理的原始数据。排除细胞培养实验、外周血、腹主动脉瘤、动物模型以及将动脉瘤与夹层合并而无法分离标签的系列（GSE314998）。10 个疾病对对照 bulk 队列符合条件（表 1）。具有疾病标签之下结构的数据集用于患者级分析：ATAA（GSE155468）和 TAAD（GSE213740）患者与对照升主动脉的单细胞 RNA 测序（scRNA-seq）；急性 TAAD 对对照的 scRNA-seq（GSE189795）；急性、亚急性和慢性夹层及对照的 scRNA-seq（GSE222318）；患者内瘤体与瘤颈样本（GSE140947）；夹层与扩张主动脉四层中膜的微区域 RNA-seq（GSE318877，3 例夹层与 4 例扩张）[19]；以及 GSE26155 的有序扩张分组，其中"临界"组仅保留用于趋势分析。GSE219204（培养 SMC）因无基因水平数据未分析。

### 预处理

所有队列以同一 Python 流程处理（附加文件 2）。基因符号依 NCBI gene_info 别名表统一为当前 NCBI Gene 命名；同一基因的多个探针取平均强度最高者（芯片）或求和（计数/FPKM）。平台特异步骤：GSE26155，Affymetrix Human Exon 1.0 ST 转录簇值经 GPL5175 注释，限于三叶瓣患者的内膜-中膜活检，扩张（> 45 mm）对非扩张（< 40 mm；临界样本不进入队列级分析）；GSE52093，Illumina HumanHT-12 v4 原始强度经 log2 转换与分位数标准化，保留在 ≥ 3 张芯片中检出（P < 0.05）的探针；GSE98770，因提交的 GeneSpring 矩阵不含组别相关信号，由 Agilent Feature Extraction 原始 gProcessedSignal 重新处理；GSE190635，MAS5 log2 值；GSE140947，基因计数，每受试者两个区域样本在队列级分析中求和、在患者内分析中分开；GSE153434、GSE202267 与 GSE294606，基因计数；GSE235161，总 RNA 文库 FPKM 限于蛋白编码基因、换算为 TPM 并做分位数标准化（原始值由 rRNA 与 7SL 转录本主导）；GSE267434，转录本级 FPKM 按基因求和（Ensembl 110）；GSE147026，提交的标准化值。GSE318877 的打孔样本计数在患者级估计中按患者求和，在分层分析中按患者与中膜层（腔侧 I、L、M、外膜侧 A）平均。GSE235161 提交的样本名代码无法映射到其论文所述亚组，因此该队列仅按动脉瘤对对照分析（表 1）。

### 差异表达与跨队列一致性

计数数据以 DESeq2（PyDESeq2 实现；Wald 检验，设计 ~ group）估计差异表达 [20,21]；芯片与标准化数据用 limma-trend 模型（inmoose 实现）[22]。假发现率（FDR）以 Benjamini–Hochberg 控制 [23]。队列间一致性为一对队列中任一队列显著（FDR < 0.05）的基因上调节统计量的 Spearman 相关。GSE190635 与其他所有 TAAD 队列负相关（ρ = −0.37 至 −0.68）；不用于主要分析，在一致性矩阵中显示，其纳入作为敏感性分析（附加文件 5）。其原始 CEL 文件存在但在本环境中无法重新处理。

### 随机效应 meta 分析与基因类别

每种疾病的 log2 倍数变化以 DerSimonian–Laird 随机效应模型合并（τ² 截断为零；主估计量）[24,25]；Hartung–Knapp–Sidik–Jonkman（HKSJ）校正为敏感性分析，并以完全相同的效应、方差与缺失处理与 statsmodels 逐基因核对（附加文件 5）。仅考虑在 ≥ 75 % 队列中存在的基因。共识差异表达基因（DEG）要求随机效应 FDR < 0.05、≥ 75 % 队列方向一致、≥ 50 % 队列同方向名义 P < 0.05。稳健性以留一队列（LOCO）重分析评估；对 TAAD 另以全部 15 个 4-of-6 队列子集评估以匹配 ATAA 的队列数。疾病间差异以 z_diff = (LFC_TAAD − LFC_ATAA)/√(SE_TAAD² + SE_ATAA²) 检验并做 FDR 校正。当 FDR_diff < 0.05 且该基因为该病共识 DEG 时归为该病富集；在另一病中的零等效性双单侧检验（TOST；|log2FC| < 0.5，90 % CI）增加一个描述性层级（"富集且在另一病中近零"）；0.5 为声明的操作性界值而非经验证的生物学阈值（附加文件 6）。疾病排序与差异排序的预排序 GSEA（gseapy；MSigDB v2024.1 Hallmark、KEGG legacy、Reactome）保留为描述性通路总结 [26–28]。

### 两层模块体系与程序分

模块分两层定义。**Tier 1** 模块在患者级分析之前由外部来源固定（附加文件 7）：SMC 收缩（24 个基因；KEGG 血管平滑肌收缩与 Reactome 平滑肌收缩的交集并入列出的核心基因 *MYH11*、*ACTA2*、*CNN1*、*TAGLN*、*MYOCD*、*LMOD1*、*SMTN*、*MYLK*、*DES*、*CALD1*、*TPM2*、*MYL9*、*ACTG2*、*SYNPO2*、*PLN*、*KCNMB1*、*SORBS1*、*CASQ2*）、钙处理（57）、细胞-基质黏附（43）、ECM-胶原（300）、基质降解（158）、缺氧（Hallmark；200）、糖酵解（25）、氧化应激/NFE2L2（146）、MYC-核糖体（387；去核糖体蛋白基因版本 298 为敏感性集）、p53-DNA 损伤（264）、NF-κB/IL-6 炎症（269）、干扰素-α 应答（Hallmark；97）、IFNG_specific（Hallmark 干扰素-γ 应答去除与干扰素-α 集共有的基因；127）和 MHC II 类抗原呈递（18 个明确列出的基因，含 *CD74*、*CIITA* 与 RFX 因子）。9 个细胞签名（各 ≤ 50 个基因）为整合 scRNA-seq 图谱中**仅对照细胞**的一对多标志，因此疾病相关表达不会进入标志定义。**Tier 2** 签名由发现数据构建（共识 DEG 与通路主题的交集），仅在留一队列折内使用，且在每折内不含被保留队列重建（附加文件 8）。

程序分在每个数据集中以相同方式、不参照标签计算：每个基因在数据集的所有样本上标准化（scRNA-seq 中在细胞类型内对所有患者的伪 bulk 谱标准化），模块分为其检出成员（≥ 5 个基因）标准化值的均值。C 为 SMC 收缩模块分。I 为六个损伤应答模块（缺氧、糖酵解、氧化应激/NFE2L2、MYC-核糖体、p53-DNA 损伤、NF-κB/IL-6）各自在数据集内标准化后的等权均值。R = I + C 表示损伤应答相对于收缩表达减弱的状态（令 L = −C，则 R = I − L）。R 回答的是：在固定尺度下，两类程序的相对状态是否与背景（疾病、病程、壁层、区域）有关；它不检验是否存在两个独立机制，且 R 不变与"两程序同步移动的单一强度轴"相容。因此每个 R 的检验均同时报告 C 与 I 的原始分数与效应。

### 分析计划与假设家族

计划在探索性队列级阶段之后、患者级分析之前锁定（附加文件 9）。定义三个推断家族。**F1（bulk，队列级）**：在每个队列内将疾病对对照对 C、I、R 的效应表示为 Hedges' g 及其标准误；按疾病以随机效应合并，合并 TAAD 效应与 ATAA 效应之差以 z 检验，三个检验做 Holm 校正；不跨队列置换疾病标签，因为疾病与队列绑定。**F2（单细胞，患者级）**：在 GSE155468、GSE213740 与 GSE189795 中各做一个收缩型 SMC 伪 bulk 内 R 的检验，精确枚举全部分组分配（各自 α = 0.05；不同数据集与不同假设不合并为一个家族）。**F3（患者内）**：在 GSE140947 中以精确符号置换（64 种分配）检验瘤体减瘤颈的 R 差值。结构上不允许推断的数据集仅报告估计：GSE222318 因各阶段取材位置（内膜-中膜对全层）不同，GSE318877 为 3 对 4 例患者（可达到的最小双侧精确 P = 0.057）。各模块、其他细胞类型以及 GSE26155 的有序扩张趋势（Jonckheere–Terpstra 统计量，置换 P）为探索性分析，在各数据块内做 BH 校正。效应量为 Hedges' g 及 bootstrap 95 % 置信区间（2,000 次重抽样）。表 2 对每个预设比较记录：方向、效应与区间、有效患者数、模块是否独立于模块构建（Tier 1）、推断栏（仅主要检验）、估计栏（区间不含/含 0）与可估计性栏（可信、有限精度或不可估计）。全文区分"现有数据不足以检出"与"无效应证据"。

### 单细胞处理与患者级伪 bulk

GSE155468（8 例 ATAA，3 例对照）与 GSE213740（6 例 TAAD，3 例供体）的细胞经质量过滤（200–6,000 个基因，线粒体读数 < 20 %）、每样本降采样至 ≤ 4,000 个细胞、标准化、log 转换，并在 2,000 个高变基因的 40 个主成分上用 Harmony 整合（Scanpy）[29,30]；Leiden 聚类以经典标志注释，去除一个双细胞样聚类后保留 12 个群体共 73,755 个细胞（附加文件 4）。GSE189795 与 GSE222318 以同一流程分别处理。每个数据集按患者与细胞类型形成伪 bulk 谱（log 标准化表达的均值；该类型细胞 < 20 个的患者被排除，这使 GSE155468 收缩型 SMC 分析中少了 1 例 ATAA 患者），在细胞类型内对全部患者计算程序分，并按细胞类型估计疾病对对照效应。全流程核对以对照参照标准化并将标签置换贯穿整个评分过程重复收缩型 SMC 的主要检验。

### 患者内、分层与有序扩张分析

GSE140947 中，在 24 个区域文库（6 例患者的瘤体与瘤颈；6 例供体的升主动脉中段与远段）上计算模块分并形成患者内差值；供体中段减远段差值作为背景。GSE318877 中，在 28 个患者×层谱上计算模块分，比较夹层与扩张患者各分数的患者内外膜侧减腔侧梯度（仅估计）。GSE26155 中，三叶瓣患者的全部内膜-中膜样本（无扩张 31、临界 6、扩张 22）一并评分，并以 Jonckheere–Terpstra 统计量（5,000 次置换）检验有序分组趋势。

### 组成签名条件化敏感性分析

为考察 bulk 程序效应是否由组织组成携带，在每个队列内以线性模型重新估计疾病对各 Tier 1 分的效应，并条件化三个预设的对照细胞签名分（收缩型 SMC、炎性髓系细胞、T 细胞）。评分前从结局模块中剔除与三个标志集中任一共有的基因；剩余 < 5 个基因的模块记为不可估计。非条件化与条件化的标准化系数按疾病以随机效应合并。这是敏感性分析，而非细胞内在成分与组成成分的分解。

### 支持性分析：跨队列可区分性与 Tier 2 签名

无泄漏的留一队列弹性网分类器（对照参照谱上的 TAAD 对 ATAA；L1 比 0.5，C = 0.05，类别平衡）在每折内重新计算训练队列的随机效应共识 DEG [31,32]；第二方案每次保留一个 ATAA 与一个 TAAD 队列（24 折），使每个测试折都含两个类别。Tier 2 签名在每折内由训练队列重建并在保留队列中评分；这评价的是构建过程而非固定基因列表。

### 软件与代码

分析使用 Python 3.11 及 pandas、SciPy、statsmodels、scikit-learn、PyDESeq2、inmoose、gseapy、Scanpy 与 harmonypy。全部代码、模块定义、结果表与图均已提供（见数据与材料可得性）。

## 结果

### 10 个队列与跨队列可重复性的不对称

10 个 bulk 队列符合条件：4 个 ATAA 队列（63 例患者，52 例对照）与 6 个 TAAD 队列（41 例患者，35 例对照；表 1；图 1a）。各队列差异表达在 ATAA 队列分别得到 1,966、1,186、7 和 81 个 FDR < 0.05 基因，在 TAAD 队列分别为 1,928、3,097、1,323、1,408、1,051 和 1,161 个（附加文件 3）。跨队列一致性不对称（图 1b）：TAAD 队列在四种平台上相互一致（Spearman ρ 中位 0.50，范围 0.23–0.71），ATAA 队列一致性弱（中位 0.12，范围 −0.09 至 0.22），ATAA 与 TAAD 队列之间一致性低（中位 0.10）。GSE190635 与其他所有 TAAD 队列负相关，处理方式见方法。

随机效应 meta 分析识别 1,917 个 TAAD 共识 DEG（上调 1,100，下调 817；I² 中位 59 %）和 171 个 ATAA 共识 DEG（上调 96，下调 75；图 1c）。HKSJ 下二者降至 354 和 0；LOCO 保留率 TAAD 为 0.74（中位），ATAA 为 0.28，去除 GSE26155 仅保留 11 % 的 ATAA 共识。与 ATAA 队列数匹配的 4-of-6 TAAD 子集给出 1,305–1,992 个共识基因、两两 ρ 中位 0.38–0.57、最小 LOCO 保留率 0.37–0.55（ATAA：0.11）。全基因组水平两病 meta 效应弱相关（ρ = 0.11）；38 个基因在两病中均为共识 DEG（期望 27）。差异检验将 778 个基因归为 TAAD 富集（其中 595 个经 TOST 在 ATAA 中近零），34 个归为 ATAA 富集（21 个在 TAAD 中近零），另有 27 个共享同向与 9 个反向基因（附加文件 6）。此不对称是关于基因级发现在纳入队列间统计可重复性的陈述，不代表生物学稳定性。meta 排序的 Hallmark GSEA 总结了基因级图景：27 个基因集在 TAAD 中富集、13 个在 ATAA 中富集，26 个在两排序间不同，以 MYC 靶基因、mTORC1、E2F、G2M、糖酵解与未折叠蛋白反应（TAAD 更高）以及干扰素-α 应答、UV-response-down 与顶端连接（ATAA 更高；肌生成 TAAD −2.24 对 ATAA −1.32）为首（附加文件 7）。

### 在队列水平，哪些程序区分两种组织状态？

采用 Tier 1 分数（图 2a–c；表 2），收缩程序 C 在每个 TAAD 队列均降低（g −1.33 至 −2.68；合并 −2.11，95 % CI −2.68 至 −1.54；I² 0 %），在 ATAA 中不一致（4 个队列中 2 个为负；合并 −0.34，−0.94 至 0.25）。两病之差为 −1.77（−2.59 至 −0.94；Holm P < 0.001；F1 主要检验）。损伤应答综合分 I 在 6 个 TAAD 队列中的 5 个升高（合并 +1.75，0.73 至 2.77），在全部 4 个 ATAA 队列中较小幅度升高（+0.56，−0.06 至 1.18）；差值（+1.19，0.00 至 2.39）为临界（P = 0.050；Holm P = 0.10）。相对状态 R = I + C 在每个队列中均接近零（−0.62 至 +0.96），两病无差异（+0.20 对 +0.24；差值 −0.04，−0.82 至 0.74；P = 0.92）。在 C 对 I 的平面上，TAAD 队列沿 I = −C 线（R 不变线）远离原点，ATAA 队列在同一条线上靠近原点（图 6a）。

在模块水平（图 2d），糖酵解（ATAA g +0.53 对 TAAD +1.81；差异 BH = 0.011）、氧化应激/NFE2L2（+0.19 对 +1.51；BH = 0.022）、MYC-核糖体（+0.18 对 +1.28；BH = 0.022）和钙处理（−0.25 对 −1.79；BH = 0.002）在两病间不同，而缺氧（+0.56 对 +1.10）、NF-κB/IL-6 炎症（+0.65 对 +1.06）和 p53-DNA 损伤（+0.38 对 +1.51；BH = 0.081）未达阈值。干扰素-α 应答（+0.80，0.40 至 1.20；4/4 队列）、IFNG_specific（+0.71）、MHC II 类（+0.91，0.02 至 1.80）、ECM-胶原（+0.60，0.21 至 0.99）与基质降解（+0.52）在 ATAA 中升高且合并区间不含零，但与 TAAD 的差异不显著（BH 0.095–0.65；MHC II 类在 TAAD 中为 −0.57）。对照细胞签名显示内皮（−2.01）、成纤维细胞（−1.10）和收缩型 SMC（−1.67）签名仅在 TAAD 中丢失（差异 BH ≤ 0.007），T 细胞（+0.84）、NK（+0.70）、B 细胞（+0.65）与浆细胞（+0.63）签名在 ATAA 中升高，与 TAAD 的差异为临界（BH 0.08–0.30）。近期图谱报告的 ENO1/MIF 轴基因 [15] 表现与糖酵解模块一致：*ENO1*（TAAD log2FC +0.90，ATAA +0.18；差异 FDR 0.002）、*MIF*（+0.74 对 +0.21；0.007）、*PKM*、*LDHA* 与 *HK2* 均为 TAAD 富集且经 TOST 在 ATAA 中近零。

条件化三个对照细胞签名分（图 2e；附加文件 10）使 TAAD 对 C（标准化系数 −1.47 → −0.44，−0.75 至 −0.13）和 I（+1.41 → +1.33，0.64 至 2.03）的效应减弱但保留，糖酵解、氧化应激、MYC-核糖体、p53 与 NF-κB/IL-6 效应亦保留。在 ATAA 中，C 的小幅效应保留（−0.30 → −0.26，−0.50 至 −0.03），而损伤综合分（+0.56 → −0.14）、干扰素-α（+0.78 → +0.24，−0.10 至 0.58）、MHC II 类（+0.90 → +0.39）与 ECM-胶原（+0.64 → +0.29）效应变得不确定。因此在此设计下，队列级的 ATAA 相关免疫与基质信号无法与组织组成分离，而 TAAD 相关的收缩丢失与代谢/应激激活可以。

### 差异是否存在于同一类壁细胞内？

在收缩型 SMC 伪 bulk 谱中（图 3；表 2），两数据集的 C 在疾病中均低于对照（GSE155468，ATAA 7 对 3 例，g −0.95，95 % CI −3.93 至 −0.16，精确 P = 0.167；GSE213740，TAAD 6 对 3 例，g −0.98，−4.16 至 −0.05，P = 0.167），I 无明确变化（−0.18 与 +0.35；区间含零）。R 的主要检验未通过（GSE155468，g −0.58，−1.65 至 0.33，120 种分配的精确 P = 0.367；GSE213740，g −0.35，−1.27 至 1.33，84 种分配的 P = 0.786）；以对照参照标准化的全流程核对给出 P = 0.150 与 0.048。各患者分数分布很宽（图 3a–d）：R 最低的 ATAA 患者兼有收缩程序的轻度下降与损伤综合分的明显下降；每例动脉瘤患者的收缩型 SMC 数（50–905）远少于对照（263–2,165），与此前注意到的解离偏倚一致。探索性分析中，ATAA SMC 的 MHC II 类分高于对照（g +1.27，0.66 至 3.33；精确 P = 0.067；BH = 0.90），TAAD SMC 的干扰素-α 分较低（g −1.00；P = 0.167）。在其他细胞类型中（图 3e–g）所有区间均很宽；唯一精确 P < 0.05 的探索性效应为 TAAD 患者 T 细胞中较低的损伤综合分（g −1.96，−8.66 至 −1.41；P = 0.024）。每个数据集仅 3 例对照，这些估计精度有限；数据未显示细胞内相对状态的差异，也不能排除中等大小的差异。

### 两类程序与夹层病程的关联是否不同？

[[STAGE_RESULTS_ZH]]

### 程序在患者内和沿扩张程度如何分布？

在动脉瘤患者内（GSE140947；图 5a–c），瘤体减瘤颈的 R 差值为 −0.97（−2.42 至 0.47；64 种分配的精确符号置换 P = 0.188；F3 主要检验未通过），两个成分在瘤体中均较低（ΔC −0.37，−0.89 至 0.15，P = 0.156；ΔI −0.60，−1.60 至 0.39，P = 0.188）；供体中段减远段差值为正且同样不确定（ΔR +0.95，−1.44 至 3.34）。在 GSE318877 的中膜各层（图 5d），C 与 I 在患者内各层间变化而无一致方向，外膜侧减腔侧梯度在夹层与扩张患者间无差异（C，g −0.23；I，−0.09；R，−0.17；区间均跨 −5 至 +3）。在同一队列的患者水平，夹层中膜的 C（g −1.11，−4.39 至 −0.02）、I（−0.47）与 R（−0.73）均低于扩张中膜，由 3 对 4 例患者估计（描述性精确 P 0.23–0.57；可达到的最小值 0.057）。

在 GSE26155 的有序扩张程度中（无 31、临界 6、扩张 22；图 5e；附加文件 11），C 随扩张下降（Spearman ρ −0.35；Jonckheere–Terpstra P = 0.006；BH = 0.016；组均值 +0.17、+0.14、−0.28），钙处理亦然（BH = 0.016），而损伤综合分无趋势（ρ +0.03；BH = 0.86），其六个模块均无趋势（BH ≥ 0.58）。因此 R 仅通过 C 下降（ρ −0.14；BH = 0.42）。相反，干扰素-α（ρ +0.38；BH = 0.012）、IFNG_specific（+0.32；0.030）、MHC II 类（+0.49；0.010）、ECM-胶原（+0.32；0.030）、基质降解（+0.30；0.030）、细胞-基质黏附（+0.32；0.030）以及 T 细胞、B 细胞、NK、巨噬细胞与炎性髓系签名（ρ +0.38 至 +0.46；BH ≤ 0.015）随扩张升高。因此，在单一内膜-中膜样本队列内，扩张伴随收缩程序下降以及免疫与基质程序，而不伴随夹层特有的损伤综合分。

### 程序组织方式与跨队列可区分性

图 6 在不含时间或因果顺序的前提下总结了状态对照。队列水平获支持的是收缩丢失以及糖酵解、氧化应激与 MYC-核糖体模块的差异；单一队列内获支持的是 C 随扩张下降及干扰素、MHC II 类、淋巴细胞与 ECM 程序随扩张升高。尚未确定的是：相对状态 R 是否在同一类壁细胞内不同，[[STAGE_FIG6_ZH]] 以及 ATAA 相关免疫与基质信号是否为细胞内在。两种组织状态跨队列仍可区分：在每折内重新计算共识特征的嵌套留一队列分类器达到合并队列外 AUC 0.90（准确率 83.7 %；每个 TAAD 队列 ≥ 71 %；GSE26155 91 %、GSE140947 100 %、GSE235161 68 %、GSE202267 50 %），每次保留两病各一个队列的 24 折 AUC 中位 1.00（四分位距 0.93–1.00；附加文件 12）。在每折内重建的 Tier 2 签名与 Tier 1 结果一致（附加文件 8）：TAAD 签名在保留的 TAAD 队列中 g 中位 1.7–2.9（SMC 收缩相关 2.77，最小 1.73，6/6 折 P < 0.05），在保留的 ATAA 队列中 g 中位 0.2–0.7；ATAA 签名在保留的 ATAA 队列中 g 中位 0.9–1.7，但 ECM、抗原呈递与淋巴细胞主题仅能在 4 折中的 3 折构建（保留 GSE26155 时剩余共识基因过少），其在 TAAD 队列中的效应范围为 −2.1 至 +2.8。

## 讨论

本研究以同一组外部定义的分数在 10 个 bulk 队列和 6 个患者级数据集中追踪两类程序。在组织水平，ATAA 与 TAAD 的主要差异在于收缩丢失的幅度（Holm P < 0.001），其次是损伤综合分（Holm P = 0.10），后者主要由糖酵解、氧化应激与 MYC-核糖体模块承载。两类程序的相对状态在两病间无差异，所有队列均靠近 R 不变的直线。在收缩型 SMC 内，两病 C 的方向相同，相对状态无明确变化；在动脉瘤患者内瘤体与瘤颈无差异；在单一队列内，扩张伴随收缩程序下降以及免疫与基质程序而无损伤综合分趋势。[[STAGE_DISCUSSION_1_ZH]]

### 相对状态检验能与不能说明什么

受检假设是两病以不同方式组织收缩与损伤应答程序，若成立，应表现为 R 在疾病间、病程间、壁层间或区域间的差异。现有数据中该关系无差异，故假设未获支持；但也未被否定，因为 R 不变与"两程序同步移动的单一强度轴"相容，且患者级检验可达到的最小 P 值为 0.017–0.057，置信区间宽到足以容纳中等效应。数据确实显示的是强度差异：同一组"收缩丢失 + 损伤应答"的组合在组织水平于夹层中强表达、于动脉瘤中弱表达，代谢与应激模块提供了最清晰的区分。这一强度解读与近期单细胞图谱描述的糖酵解 ENO1/MIF 轴 [15] 一致——其基因在本研究跨队列中为 TAAD 富集且在 ATAA 中近零；并补充说明该轴表现为更广泛损伤综合分的一部分，而非夹层特异模块。

### 动脉瘤相关程序

ATAA 队列与对照的差异在于干扰素、MHC II 类、淋巴细胞与 ECM 程序，在单一队列内这些程序随扩张程度升高而收缩分下降。有两点须谨慎。第一，在队列水平这些效应在条件化对照细胞签名分后变得不确定，因此用 bulk 数据无法将其与组织组成分离；动脉瘤 SMC 伪 bulk 中较高的 MHC II 类分（7 例患者，精确 P = 0.067）有提示性但效能不足。第二，有序扩张队列由三叶瓣患者的内膜-中膜样本组成，比较的是扩张程度而非疾病与对照组织，其结果为探索性。在此限制下，动脉瘤状态可表述为收缩程序下降与免疫、基质程序的组合，而不含伴随夹层的缺氧、糖酵解、氧化应激、MYC 与 p53 模块。这是组织水平程序组成的不同，而非同一类壁细胞以不同方式组织程序的证据。

### 可重复性不对称

TAAD 签名在 6 个队列、4 种平台间可重复，4-of-6 子集保持该性质，而 ATAA 共识依赖最大队列。此不对称涉及基因级发现在纳入队列间的统计可重复性——这些队列在对照组织、取材部位与病例定义上均不同——并不意味动脉瘤生物学更不稳定。它意味着已发表的单队列动脉瘤签名需谨慎解读，而程序级分数（干扰素与淋巴细胞模块在各 ATAA 队列中一致）比基因列表更稳健。

### 与相邻研究的关系

被概括为"侵蚀或爆炸"的细胞异质性比较 [16] 与 80 例主动脉图谱 [15] 均报告动脉瘤与夹层不同。本研究不重复"两种不同疾病"的对照，而是检验差异如何组织，发现现有患者级结构未显示两类程序的不同组织方式，只显示不同的强度与不同的伴随组成。直接的夹层对扩张队列 [19] 给出方向一致的患者级估计（夹层中膜的 C、I、R 均较低），但 3 对 4 例患者无法检验。

### 局限性

各队列的对照组织不同（器官供体、移植供体、非扩张手术主动脉），疾病与队列绑定，因此队列级差异在队列水平检验而不混合样本。标签无关标准化使分数在数据集内可比但不跨数据集可比，故比较的是效应量而非分数。单细胞数据集各仅 3 例对照，主要检验虽精确但分配数少。[[STAGE_LIMITATION_ZH]] Tier 1 模块为外部定义，不能涵盖所有相关程序；损伤综合分按声明对六个模块等权。GSE235161 的亚组代码无法映射到已发表分组，GSE190635 无法重新处理。性别、年龄、瓣膜形态与用药信息不统一可得且未建模。手术组织无法推断时间或因果顺序。

## 结论

在组织水平，升主动脉瘤与 A 型夹层的差异主要在于收缩丢失与代谢/应激激活的幅度，而非收缩与损伤应答程序的相对状态。"两病以不同方式组织这些程序"的假设未获现有患者级数据支持，而这些数据对该问题效能不足。动脉瘤扩张伴随收缩程序下降以及干扰素、抗原呈递、淋巴细胞与基质程序，但不伴随夹层的损伤综合分。以外部定义模块为基础的程序级分数为跨队列比较主动脉组织状态提供了可重复的单元，并可作为前瞻性、标注病程、分层取材设计的基础。

## 缩写

ATAA：升主动脉瘤；TAAD：Stanford A 型主动脉夹层；SMC：平滑肌细胞；ECM：细胞外基质；GEO：Gene Expression Omnibus；DEG：差异表达基因；FDR：假发现率；BH：Benjamini–Hochberg；REM：随机效应模型；HKSJ：Hartung–Knapp–Sidik–Jonkman；LOCO：留一队列；TOST：双单侧检验；GSEA：基因集富集分析；AUC：受试者工作特征曲线下面积；scRNA-seq：单细胞 RNA 测序；MHC：主要组织相容性复合体；CI：置信区间。

## 声明

**伦理批准与知情同意。** 不适用；所有数据均为公开的去标识数据。

**发表同意。** 不适用。

**数据与材料可得性。** 所有数据集可从 GEO 获取：GSE26155、GSE140947、GSE235161、GSE202267、GSE52093、GSE153434、GSE98770、GSE294606、GSE267434、GSE147026、GSE190635、GSE318877、GSE155468、GSE213740、GSE189795 与 GSE222318。分析代码、模块定义、锁定的分析计划、各队列与合并结果表、患者级分数表、表 2 与全部图见随附代码库（scripts/、results/、figures/、manuscript/ 目录）及附加文件。

**利益冲突。** 作者声明无利益冲突。

**基金。** [待补充。]

**作者贡献。** [待补充：构思与设计；数据获取与分析；解释；起草；批判性修改。所有作者阅读并批准终稿。]

**致谢。** 感谢将原始数据集存入 GEO 的研究者。

## 参考文献

1. Isselbacher EM, Preventza O, Hamilton Black J, et al. 2022 ACC/AHA guideline for the diagnosis and management of aortic disease. Circulation. 2022;146:e334–e482.
2. Evangelista A, Isselbacher EM, Bossone E, et al. Insights from the International Registry of Acute Aortic Dissection: a 20-year experience of collaborative clinical research. Circulation. 2018;137:1846–60.
3. Pape LA, Tsai TT, Isselbacher EM, et al. Aortic diameter ≥ 5.5 cm is not a good predictor of type A aortic dissection: observations from the International Registry of Acute Aortic Dissection (IRAD). Circulation. 2007;116:1120–7.
4. Nienaber CA, Clough RE, Sakalihasan N, et al. Aortic dissection. Nat Rev Dis Primers. 2016;2:16053.
5. Rombouts KB, van Merrienboer TAR, Ket JCF, Bogunovic N, van der Velden J, Yeung KK. The role of vascular smooth muscle cells in the development of aortic aneurysms and dissections. Eur J Clin Invest. 2022;52:e13697.
6. Kurihara T, Shimizu-Hirota R, Shimoda M, et al. Neutrophil-derived matrix metalloproteinase 9 triggers acute aortic dissection. Circulation. 2012;126:3070–80.
7. Anzai A, Shimoda M, Endo J, et al. Adventitial CXCL1/G-CSF expression in response to acute aortic dissection triggers local neutrophil recruitment and activation leading to aortic rupture. Circ Res. 2015;116:612–23.
8. Son BK, Sawaki D, Tomida S, et al. Granulocyte macrophage colony-stimulating factor is required for aortic dissection/intramural haematoma. Nat Commun. 2015;6:6994.
9. del Porto F, Proietta M, Tritapepe L, et al. Inflammation and immune response in acute aortic dissection. Ann Med. 2010;42:622–9.
10. He R, Guo DC, Estrera AL, et al. Characterization of the inflammatory and apoptotic cells in the aortas of patients with ascending thoracic aortic aneurysms and dissections. J Thorac Cardiovasc Surg. 2006;131:671–8.
11. Tang PC, Coady MA, Lovoulos C, et al. Hyperplastic cellular remodeling of the media in ascending thoracic aortic aneurysms. Circulation. 2005;112:1098–105.
12. Guo DC, Pannu H, Tran-Fadulu V, et al. Mutations in smooth muscle alpha-actin (ACTA2) lead to thoracic aortic aneurysms and dissections. Nat Genet. 2007;39:1488–93.
13. Milewicz DM, Trybus KM, Guo DC, et al. Altered smooth muscle cell force generation as a driver of thoracic aortic aneurysms and dissections. Arterioscler Thromb Vasc Biol. 2017;37:26–34.
14. Pinard A, Jones GT, Milewicz DM. Genetics of thoracic and abdominal aortic diseases. Circ Res. 2019;124:588–606.
15. [Authors]. Single-cell atlas of 80 human ascending aortas identifies a glycolytic ENO1/MIF axis in aortic dissection. Adv Sci. 2026 (PMID 42107066).
16. [Authors]. Erosion or explosion: cellular heterogeneity of aortic aneurysm and dissection. Inflammation. 2026 (PMID 42120767).
17. Ramasamy A, Mondry A, Holmes CC, Altman DG. Key issues in conducting a meta-analysis of gene expression microarray datasets. PLoS Med. 2008;5:e184.
18. Leek JT, Scharpf RB, Bravo HC, et al. Tackling the widespread and critical impact of batch effects in high-throughput data. Nat Rev Genet. 2010;11:733–9.
19. Kimura N, et al. Micro-regional transcriptomics of the ascending aortic media in dissection and dilatation (GSE318877). 2026 (PMID 42045332).
20. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15:550.
21. Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. Bioinformatics. 2023;39:btad547.
22. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43:e47.
23. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Ser B. 1995;57:289–300.
24. DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7:177–88.
25. IntHout J, Ioannidis JPA, Borm GF. The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. BMC Med Res Methodol. 2014;14:25.
26. Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A. 2005;102:15545–50.
27. Fang Z, Liu X, Peltz G. GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. Bioinformatics. 2023;39:btac757.
28. Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 2015;1:417–25.
29. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018;19:15.
30. Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods. 2019;16:1289–96.
31. Zou H, Hastie T. Regularization and variable selection via the elastic net. J R Stat Soc Ser B. 2005;67:301–20.
32. Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: machine learning in Python. J Mach Learn Res. 2011;12:2825–30.
33. Lakens D. Equivalence tests: a practical primer for t tests, correlations, and meta-analyses. Soc Psychol Personal Sci. 2017;8:355–62.
34. Jonckheere AR. A distribution-free k-sample test against ordered alternatives. Biometrika. 1954;41:133–45.
35. Hedges LV. Distribution theory for Glass's estimator of effect size and related estimators. J Educ Stat. 1981;6:107–28.
36. Folkersen L, Wågsäter D, Paloschi V, et al. Unraveling divergent gene expression profiles in bicuspid and tricuspid aortic valve patients with thoracic aortic dilatation: the ASAP study. Mol Med. 2011;17:1365–73.
37. Li Y, Ren P, Dawson A, et al. Single-cell transcriptome analysis reveals dynamic cell populations and differential gene expression patterns in control and aneurysmal human aortic tissue. Circulation. 2020;142:1374–88.
38. Chen PY, Qin L, Li G, et al. Smooth muscle cell reprogramming in aortic aneurysms. Cell Stem Cell. 2020;26:542–57.
39. Pan S, Wu D, Teschendorff AE, et al. JAK2-centered interactome hotspot identified by an integrative network algorithm in acute Stanford type A aortic dissection. PLoS One. 2014;9:e89406.
40. Kimura N, Futamura K, Arakawa M, et al. Gene expression profiling of acute type A aortic dissection combined with in vitro assessment. Eur J Cardiothorac Surg. 2017;52:810–7.
41. Zhou Z, Liu Y, Zhu X, et al. Exaggerated autophagy in Stanford type A aortic dissection: a transcriptome pilot analysis of human ascending aortic tissues. Genes (Basel). 2020;11:1187.
42. Barrett T, Wilhite SE, Ledoux P, et al. NCBI GEO: archive for functional genomics data sets—update. Nucleic Acids Res. 2013;41:D991–5.

*分析时尚无关联原始论文的数据集（GSE235161、GSE202267、GSE294606、GSE267434、GSE147026、GSE190635、GSE213740、GSE189795、GSE222318）以 GEO 登录号引用 [42]。文献 [15]、[16]、[19] 的细节待按最终发表版本补全。*

## 图注

**图 1 证据结构与队列可比性。** **a** 探索阶段之后锁定的分析计划下的数据结构、分析单位与推断状态；程序分定义见底部。**b** 队列间调节差异表达统计量在任一队列显著（FDR < 0.05）基因上的两两 Spearman 相关；GSE190635（灰色）与其他所有 TAAD 队列负相关。**c** 随机效应主分析下基因级发现的可重复性：ATAA、TAAD 与 15 个 4-of-6 TAAD 子集的共识 DEG、留一队列保留率、疾病内一致性、HKSJ 敏感性与疾病差异基因。

**图 2 bulk 队列中的 Tier 1 程序。** **a–c** 各队列（蓝，ATAA；红，TAAD）收缩分 C、损伤综合分 I 与相对状态 R = I + C 的 Hedges' g（疾病对对照，95 % CI）及随机效应合并估计（菱形），以及队列级 TAAD 减 ATAA 差值与 Holm 校正 P（F1 主要家族）。**d** 各 Tier 1 模块与对照细胞签名按疾病的合并 g 及其差值；星号，合并 95 % CI 不含 0；井号，差值 BH < 0.05（探索性）。**e** 组成签名条件化敏感性：条件化收缩型 SMC、炎性髓系与 T 细胞签名分前（空心）后（实心）的合并标准化疾病系数。

**图 3 细胞类型内的患者级程序分（scRNA-seq）。** **a, c** GSE155468（ATAA）与 GSE213740（TAAD）收缩型 SMC 伪 bulk 谱的各患者 C 与 I 分；点大小为细胞数；虚线为 R = 0。**b, d** 各组的 C、I、R 及 Hedges' g 与主要检验（R）的精确枚举 P。**e–g** 两数据集各细胞类型的 C、I、R 的 Hedges' g（95 % CI，截断）及患者数（病例 v 对照）。**h** 主要检验与全流程核对。

**图 4 夹层病程各阶段的程序分（GSE222318 仅估计；GSE189795 主要检验）。** **a–c** GSE222318 收缩型 SMC 中各患者按阶段的 C、I、R；标记形状为取材位置（内膜-中膜或全层），各阶段不同。**d** 各阶段对比的 Hedges' g（95 % CI），含急性对非急性及同取材位置的急性对亚急性；未做推断检验。**e** GSE189795 急性 TAAD 对对照的收缩型 SMC，含 R 的精确主要检验。

**图 5 患者内与有序扩张程度中的程序。** **a–c** GSE140947 患者内区域配对（动脉瘤瘤颈至瘤体；供体升主动脉远段至中段）的 C、I、R 及精确符号置换 P（R 为 F3 主要检验）。**d** GSE318877 各患者中膜分层谱（腔侧 I 至外膜侧 A）的 C（实线）与 I（虚线）（红，夹层；蓝，扩张）；仅估计。**e** GSE26155 有序扩张分组（无、临界、扩张）的 C、I、干扰素-α、MHC II 类、T 细胞签名与 ECM-胶原，附 Spearman ρ、Jonckheere–Terpstra 置换 P 与 BH 校正 P（探索性）。

**图 6 程序组织方式。** **a** 各数据集对 C（x）与 I（y）的效应：bulk 队列（圆）与合并估计（菱形）、收缩型 SMC 患者级效应（星）与直接夹层对扩张对比（十字）；虚线 I = −C 表示相对状态不变。**b** 按程序的状态对照：各病合并效应、两病是否不同、队列内或患者级证据；区分获支持与尚未确定的条目；不隐含时间或因果顺序。

## 表

**表 1** 纳入研究的数据集（manuscript/Table1_cohorts.csv）。

**表 2** 预设比较的证据矩阵：方向、效应与 95 % CI、有效 n、Tier 1 状态、推断（仅主要检验）、估计与可估计性（manuscript/Table2_evidence_matrix.csv；正文版本限于 C、I、R，见 Table2_evidence_matrix_main.csv）。

## 附加文件

**附加图。** 图 S1 在每个留一队列折内重建并在保留队列中评分的 Tier 2 签名。图 S2 嵌套留一队列分类器：各折 AUC 与准确率的 4 × 6 矩阵。图 S3 全部 Tier 1 模块的组成签名条件化敏感性。图 S4–S6 基因级比较、meta 排序的 Hallmark GSEA 与整合单细胞图谱（探索阶段）。

附加文件 1：GEO 检索策略与筛选。附加文件 2：分析代码（scripts/）。附加文件 3：各队列 QC 与差异表达表。附加文件 4：单细胞处理、注释与对照细胞标志集。附加文件 5：随机效应 meta 分析、HKSJ 核对、LOCO 与 4-of-6 子集表、GSE190635 敏感性。附加文件 6：含疾病间差异检验与 TOST 层级的基因类别。附加文件 7：Tier 1 模块成员表与 meta 排序的 Hallmark GSEA。附加文件 8：Tier 2 签名定义与折内评价。附加文件 9：锁定的分析计划与统计执行说明。附加文件 10：组成签名条件化敏感性表。附加文件 11：患者内、分层与有序扩张分数表。附加文件 12：嵌套留一队列分类器结果。
