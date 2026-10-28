# Writing and proof exposition audit based on original papers by Noriaki Matsushima

**Date:** 2026-10-28. **Task:** infer transferrable *economic-theory writing and proof design*, not copy phrases or treat AI-generated blog abstracts as original prose.

## Blog as bibliographic index only

https://nmatsush.blog.shinobi.jp/ announces research and identifies related journal papers. The September 2026 posts explicitly say some Japanese summaries were prepared using generative AI. Their style was **not** used as evidence of Matsushima's own English academic prose. We instead inspected published journal full texts where available or identified earlier university discussion-paper versions. Version status is recorded below.

## Five selected papers and observed exposition techniques

| Original work / access status | Precise sections/pages actually examined | Relevant technique, not copied language | Direct application to this manuscript |
|---|---|---|---|
| **Lu, Matsushima & Shekhar (2025)**, *International Journal of Industrial Organization* 103, 103200, OPEN PUBLISHED PDF | Published PDF Introduction pp.1–2, Related Literature p.3, Proposition 1 p.9, proof of Proposition 1 in Appendix p.16 | Introduce a live policy/market puzzle, position mechanism against the closest prior result, derive economic predictions in the model, and supply explicit sign inequalities in appendix; paper separately devotes Sections 5–7 to extensions, closest competitor and policy | Stage11 revised Intro now distinguishes UK/EDIH facts from the private fee mechanism, names the **closest** competing models, states the numerical vulnerability early; exact sign of foreign welfare recorded in table |
| **Matsushima, Nishikawa & Qiu (2026)**, *Contemporary Economic Policy*, open full journal HTML, DOI 10.1111/coep.70050 | Original article Introduction, model and conclusion (online HTML sections); final page numbering not stabilised as early view | Two geographic labour markets link through one product market; derive why a regulation can influence a rival through strategic substitution; interpret local gains/foreign losses separately before discussing extensions | Clarify which region finances basic service and where downstream profits/displacement fall. Do **not** confuse local and national incidence; contrast physical boundary with market integration |
| **Cong & Matsushima (2026)**, *Management Science*, open published journal HTML DOI 10.1287/mnsc.2024.08024 | Published online abstract, Sections 1–2 and model/results discussion available in journal HTML, supplemental source indicated | Sequential information/price mechanisms described in economic sequence and then decomposed into cross-market repercussions; scope of consumer-surplus and profit outcomes separated | Stage11 revised Intro explicitly sets timing and separates mass, moment, venture profits and regional incumbents. Avoid asserting correlation not in the data |
| **Matsushima & Shinohara (2014)**, *Journal of Economic Behavior & Organization* 106:428–441; FULL 2011 Osaka ISER preprint, final version **not** verified word-for-word | 2011 preprint Introduction pp.1–3, Analysis p.9, Proposition 1 p.13, Appendix. Published final bibliographic identity cross-checked separately | Early statement of tangible buyer–supplier choice, backward induction, explicit one-vs-two-buyer cases, then a short proposition with intuitive cost/bargaining interpretation | Retain P-only full-coverage kink, Both premium case and coupled P/P case as separate active sets; put full revenue-regime proof in appendix and avoid jumping from FOC to equilibrium |
| **Kawasaki & Matsushima (2026)**, Osaka OSIPP DP-2026-E-005, FULL discussion paper, **not** peer reviewed | Introduction pp.2–4, piecewise profit Lemma 1 pp.11–12, Proposition 1 pp.12–13, proof Appendix A.1 pp.27–28 | Piecewise equilibrium profits first, explicit tie-breaking, proof of each active-set optimum, figures only when they clarify thresholds; conditional comparative statics tied to private/public ownership | Keep finite deviations, kink, tie/multiple equilibria and exact cost band visible. New sensitivity plot is generated from Python, not invented smooth ownership boundaries |

## Design rules applied and ongoing limits

1. **Introduction**: first paragraph is a specific empirical policy puzzle; second says why endogenous participation composition matters; third states game/timing; fourth states conditional mechanisms; fifth explicitly identifies prior literature threats; sixth quantifies tiny entry interval and negative scope; seventh states diagnostic institutional conclusion.
2. **Model**: original economic actors and actions are consolidated in one section, including four binary decisions. No physics of co-location is insinuated.
3. **Proofs**: retain central bridge equations `q_i=(1-ell_i²)/2`, `Delta L_0=C ell_i`, `Delta L_1=C ell_i²/4`, fee kink and two regional accounts. Place all simple differentiations in a traceable appendix, not opaque code.
4. **Numerical**: show three scenarios, local/foreign welfare signs and cost width. Draw phase diagram only with actual discrete solver classifications and remember ties/multiplicity.
5. **Conclusion**: report model-specific diagnostic only, not public-building optimum or general innovation prescription.
6. **Caveat**: none of these techniques is exclusive to Matsushima or copied from him. The audit is a reader-centered checklist gleaned from actual articles and theoretical economic writing; no claim is made that a particular phrasing is uniquely his.

## Exact original article links

- Lu–Matsushima–Shekhar published full article: https://ir.library.osaka-u.ac.jp/repo/ouka/all/103287/IntJIndOrganiz_103_103200.pdf
- Matsushima–Nishikawa–Qiu published full HTML: https://onlinelibrary.wiley.com/doi/full/10.1111/coep.70050
- Cong–Matsushima published HTML: https://pubsonline.informs.org/doi/10.1287/mnsc.2024.08024
- Matsushima–Shinohara university preprint (2011): https://www.iser.osaka-u.ac.jp/static/resources/docs/dp/DP0808.pdf and final bibliographic record DOI 10.1016/j.jebo.2014.07.014
- Kawasaki–Matsushima (2026) working paper: https://www.osipp.osaka-u.ac.jp/archives/DP/2026/DP2026E005.pdf

**No paper wording was directly copied.** Unavailable final-paper contents are not represented as checked. This writing memo alone does not certify novelty or publication readiness.
