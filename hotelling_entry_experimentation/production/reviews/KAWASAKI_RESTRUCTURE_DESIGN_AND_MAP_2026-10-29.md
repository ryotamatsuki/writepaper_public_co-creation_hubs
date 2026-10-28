# Frozen-result reconstruction: design and complete correspondence

Baseline: a463d3340d8e6c2307136a6d64db484a5fefc5fa (PR #69). Parent: cf896985295278136f033244d95b3a8015927f7e (PR #72). No theory extension. Design made from the full 37-entry ledger and main/supplement/code evidence before prose revision.

## Title decision

| Candidate | Assessment |
|---|---|
| Technical Complementarity, Private Screening, and the Value of Public Innovation Access | Accurate original; leads with one mechanism rather than the full decision chain. |
| Private Innovation Intermediation, Regional Incidence, and Public Access | Selected: names actual private/public decisions and existing incidence primitives without claiming a general regional equilibrium theorem. |
| Regional Characteristics and Optimal Public Innovation Hubs | Rejected: suggests exhaustive regional comparative statics, physical facilities and a broader optimum than analyzed. |

## Architecture

The model remains §2 verbatim. Private pricing and Lemma 2 move to §3. Regional incidence (Lemma 1), conditional public cost ceilings (Lemma 3), and Proposition 1 form §4. Section 5 compares local/national boundaries, entry incentives and the certified three-configuration table. Section 6 holds Corollary 1, the exact single-index identity, repriced equal-mean technology comparisons, finite-matching certificate, decomposition, economic scale and all negative portability records. Measurement is §7, related literature §8, conclusion §9; institutional motivation is a subsection of §8. Original numbered result IDs are retained although presentation order changes. The Online Supplement is byte-identical and labels are retained.

The national wedge equation remains inside the unchanged Lemma 3 block in §4, and is interpreted in §5 by reference. This avoids splitting a frozen numbered result. Every formula and exhibit is retained, including the redundant introduction display. No new numerical graph or simulation was needed.

## Exhaustive result map

The companion CSV preserves every analytical dimension from PR #72 and adds old/new location and presentation role. Rows below are evidence records, not 37 independent discoveries. Detailed assumptions, primitives, prior-art overlap and proof limits remain available in that CSV and the unchanged PR #72 review.

| ID | Result | Old location | New location | Role |
|---|---|---|---|---|
| F01 | 地域・全国事業化帰属 | 04_propositions: conditional mechanisms | §4 Public provision / regional characteristics (04_propositions) | Separate regional incidence from national surplus before valuing access. |
| F02 | 民間料金と技能選別 | 04_propositions: conditional mechanisms | §3 Private intermediation (03_continuation) | Explain private selection before the public counterfactual. |
| F03 | 同一参加数の能力構成 | 04_propositions: conditional mechanisms | §6 Technology / certification / robustness (04_diagnostic) | Separate fixed-mass quality from reoptimised technology comparisons. |
| F04 | 地方・全国公共追加閾値 | 04_propositions: local/national thresholds | §4 Public provision / regional characteristics (04_propositions) | Define the exact conditional local and national cost boundaries. |
| F05 | スクリーニングと地域損失の符号逆転 | 04_diagnostic: central proposition / equal-mean diagnostic | §4 Public provision / regional characteristics (04_propositions) | Identify the primitive condition for the opposing local-incidence/value-slope signs. |
| D01 | Pのみ最適料金 | 03_continuation: private pricing / entry | §3 Private intermediation (03_continuation) | 準備 — 内点f*=(U+T)/2;完全カバレッジ角解f*=U |
| D02 | B共存価格 | 03_continuation: private pricing / entry | §3 Private intermediation (03_continuation) | 準備 — fB=Δ/2,有料率1/2,売上NΔ/4 |
| D03 | P/P同時価格ゲーム | 03_continuation: private pricing / entry | §3 Private intermediation (03_continuation) | 技術証明 — カットオフ連立＋全域最適反応、基準3環境で全球純均衡唯一 |
| D04 | 公共基本の売上クラウディングアウト | 03_continuation: private pricing / entry | §3 Private intermediation (03_continuation) | 補助 — P料金収入>B料金収入 |
| D05 | 民間参入支配条件 | 03_continuation: private pricing / entry | §3 Private intermediation (03_continuation) | 補助 — FP<NΔ/4で民間が厳密参入採算 |
| D06 | 技能と参加量の単一指数 | 04_propositions: conditional mechanisms | §6 Technology / certification / robustness (04_diagnostic) | 新規性制約 — q=m−m²/2 |
| D07 | 同平均技術の民間参加閾値 | 04_diagnostic: central proposition / equal-mean diagnostic | §6 Technology / certification / robustness (04_diagnostic) | 機構識別 — 料金再最適化でellH>=ellL、内点は厳密差>0 |
| D08 | 技能Hとpooled Lの公共閾値差 | 04_diagnostic: central proposition / equal-mean diagnostic | §6 Technology / certification / robustness (04_diagnostic) | 機構識別 — TH−TLの正負は条件により両向き |
| D09 | 基本サービス便益の必要条件 | 04_diagnostic: central proposition / equal-mean diagnostic | §4 Public provision / regional characteristics (04_propositions) | 主命題境界 — γ→0でT=N u(2Δ−3u)/(8Δ);ω=0かu=0では逆説不可条件 |
| D10 | 16構成と限定プランナー | 03_continuation: private pricing / entry | §5 Decentralised provision / welfare (05_numerical) | 証拠枠組 — 4主体二値逸脱を各構成で確認、plannerは16構成最大 |
| D11 | 研究努力と新商品市場 | 02_model: Cournot and research | §2 The model (unchanged) | 基礎設定 — e*=ζΠ/kで技能・λに独立、Miは既存企業と消費者の符号不定効果 |
| N01 | 非対称A/Bパラメータ | 05_numerical; supplement / certification records | §5 Decentralised provision / welfare (05_numerical) | 環境 — NA=3 NB=1,nA=5 nB=2,ωA=.8 ωB=.2,γ=.1など |
| N02 | 3技術環境のNEとplanner | 05_numerical; supplement / certification records | §5 Decentralised provision / welfare (05_numerical) | 均衡例 — skill0 (B,P)/(B,B);skill.5 (B,B)/(B,B);pooled.5 (B,P)/(B,B) |
| N03 | 有限γ二次元長方形 | 05_numerical; supplement / certification records | §6 Technology / certification / robustness (04_diagnostic) | 中心証拠 — ωB∈[.13,.27],λ∈[0,.5]全域でT>0,Tλ>0,K<0,技術増分負 |
| N04 | Bのω範囲 | 05_numerical; supplement / certification records | §6 Technology / certification / robustness (04_diagnostic) | 主命題例 — .1167164<ωB<.2781020が弱γの閾値帯 |
| N05 | Bの地元事業化は負 | 05_numerical; supplement / certification records | §6 Technology / certification / robustness (04_diagnostic) | 背景 — MB≈−.030775,ωΠ≈.0221322,K≈−.00864284,S0≈−.000568455 |
| N06 | 技能λ0と.5の費用閾値 | 05_numerical; supplement / certification records | §6 Technology / certification / robustness (04_diagnostic) | 効果量 — TB:.0404250710→.0405316701（差.0001065991） |
| N07 | 直接サービス・通常マッチング・技能項 | 05_numerical; supplement / certification records | §6 Technology / certification / robustness (04_diagnostic) | 機構識別 — skill0 .0404864−.00006136; skill.5 .04063509−.00006140−.00004202; pooled.5 .04022101−.00006129−.00011649 |
| N08 | 直接サービス増分と技能損失の相殺 | 05_numerical; supplement / certification records | §6 Technology / certification / robustness (04_diagnostic) | 効果量 — カットオフ+.0002323,直接+.00014866,技術−.00004202,通常−.0000000396 |
| N09 | 異質対均質の費用上限差 | 05_numerical; supplement / certification records | §6 Technology / certification / robustness (04_diagnostic) | 機構識別 — TH−TL≈+.0004884409（アクセス+.0004139716,技能帰属+.0000744693） |
| N10 | 3ケースの均衡切替え固定費帯 | 05_numerical; supplement / certification records | §6 Technology / certification / robustness (04_diagnostic) | 限定 — .0404250710071<FG_B<.0405316701024、幅が費用の約.2632% |
| N11 | FG_A FG_B ωBの開区間 | 05_numerical; supplement / certification records | §6 Technology / certification / robustness (04_diagnostic) | 限定 — FG_A∈(.089999,.090001),FG_B∈(.040499,.040501),ωB∈(.199999,.200001) |
| N12 | 政府閾値λ星 | 05_numerical; supplement / certification records | §6 Technology / certification / robustness (04_diagnostic) | 参考 — λ*≈.351184の条件付き判断境界 |
| N13 | 近隣Aの外部便益とplannerの乖離 | 05_numerical; supplement / certification records | §5 Decentralised provision / welfare (05_numerical) | 政策帰結 — A外部便益≈+.000414,+.000768,+.001394；planner常に(B,B) |
| N14 | Bプレミアム採算 | 05_numerical; supplement / certification records | §3 Private intermediation (03_continuation) | 例 — 売上.25、固定費.24、純民間利益+.01 |
| N15 | 396セルと不成立環境 | 05_numerical; supplement / certification records | §6 Technology / certification / robustness (04_diagnostic) | 限界 — ΔB=.8/.9,NB=.8で複数均衡;χ=.95で(B,P)へ戻る |
| N16 | Beta(2,2)で元の反転消失 | 05_numerical; supplement / certification records | §6 Technology / certification / robustness (04_diagnostic) | 限界 — λ0,.5ともNE=planner=(P,P),B政府差≈−.03056 |
| N17 | Bertrand別仕様で不成立 | 05_numerical; supplement / certification records | §6 Technology / certification / robustness (04_diagnostic) | 限界 — 差別化Bertrand代替の負の記録 |
| N18 | 16状態・全逸脱・Lean18項目・100区間 | 05_numerical; supplement / certification records | §6 Technology / certification / robustness (04_diagnostic) | 信頼性 — 3環境全料金大域認証/部分Lean/条件付き連続区間を区別 |
| E01 | 参加数と参加者能力 | 06_empirical_implications | §7 Empirical implications (unchanged) | 測定 — ∂m/∂f=−1/T;∂q/∂f=−ell/T;平均能力導関数=1/(2T) |
| E02 | 公共と民間共存の売上 | 06_empirical_implications | §7 Empirical implications (unchanged) | 測定 — B民間料金Δ/2・率1/2、P比較で収入減・条件付き存続 |
| E03 | 研究補完性と越境厚生 | 06_empirical_implications | §7 Empirical implications (unchanged) | 測定 — ∂²p/(∂θA∂θB)=λ、FN−FL=ΔWj |

## Decision / evidence boundary

The motivating question is partly answered: conditional fee and access decisions, local/national accounting, and three fully checked entry examples. A complete parameter classification of four-party equilibrium, general market-size/competition comparative statics, corrective instruments and empirical regional identification remain outside the frozen evidence. Reframing must not imply they have been solved.

A newly visible inherited strictness issue is recorded separately: the proof uses positive baseline success for strict S0 negativity, although the printed Proposition 1 does not explicitly list it. No frozen statement or proof was repaired. This is a submission blocker requiring a later authorized mathematical scope correction, not a reason to hide the issue or add a model.
