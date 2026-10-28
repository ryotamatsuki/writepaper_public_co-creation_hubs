# Stage 14 — Reference and contribution audit

**Frozen source:** 7308052c0f6299597eba4605ad9b4906966eb8f8, `manuscript/references.bib`. Existing prior-theorem mapping is retained rather than restarted.

## Reference key audit
The current BibTeX file contains **36** entries, of which **31** are cited in the active main+supplement sources. **No undefined citation keys** were found on static inspection. The five intentionally unused source records are `NishimuraOkamuro2026`, `KawasakiMatsushima2026`, `MatsushimaShinohara2014`, `MatsumuraMatsushimaIshibashi2009`, and `HacknerHerzing2026`; they are not expected to appear in the generated bibliography without `\\nocite`. They may remain in the source BibTeX file without violating the JE requirement that the *printed* list contain cited works.

## Primary-source checks made in this audit
| Key | Author/title/year/venue/page and DOI | Status |
| --- | --- | --- |
| `HuaChenZhang2016` | Publisher, *Production and Operations Management* 25(8), 1430–1448, 10.1111/poms.12555 | BIBLIOGRAPHY VERIFIED; final technical article/supplement **UNVERIFIED** |
| `KimLiXu2025` | Publisher, *Theoretical Economics* 20(3), 1081–1134, 10.3982/TE5651 | BIBLIOGRAPHY VERIFIED |
| `OstrizekSartori2023` | Publisher, *Games and Economic Behavior* 139, 26–55, 10.1016/j.geb.2023.01.011 | BIBLIOGRAPHY VERIFIED; exact author-preprint/final technical alignment **UNVERIFIED** |
| `DamianoLi2008` | Publisher, *J Eur Econ Assoc* 6(4), 789–818, 10.1162/JEEA.2008.6.4.789 | BIBLIOGRAPHY VERIFIED; final proof-version comparison **UNVERIFIED** |
| `DamianoLi2007` | Springer/RePEc, *Economic Theory* 30(2), 243–263, 10.1007/s00199-005-0058-2 | BIBLIOGRAPHY VERIFIED; final proof-version comparison **UNVERIFIED** |

Direct official sources: https://journals.sagepub.com/doi/10.1111/poms.12555 ; https://onlinelibrary.wiley.com/doi/10.3982/TE5651 ; https://www.sciencedirect.com/science/article/pii/S0899825623000143 ; https://academic.oup.com/jeea/article-abstract/6/4/789/2295772 ; https://link.springer.com/article/10.1007/s00199-005-0058-2 .

**Remaining 31 BibTeX records:** prior manuscript reports list/source-scope provenance but full item-by-item contemporary publisher verification was *not repeated*. Consequently, do not claim all 36 current metadata records were newly verified. Full DOI/title/date inspection remains an ordinary bibliography QA item; unknown details must be marked `UNVERIFIED` rather than inferred.

## Overlap / absorption
Confirmed inherited: posted-fee screening, complementary bilateral matching, mixed public/private service, screening with an externality, and the general local/national welfare wedge. The inspected nearest propositions do not directly supply the same `K_i<0`, `S_i^0<0`, `T'_i>0` **conditional regional public-cost boundary** under the maintained separate public alternative. Full technical packages of Hua (2016), final Damiano–Li articles and exact Ostrizek–Sartori version correspondence remain narrow source gaps. **No affirmative general originality certification**, no confirmed absorption requiring Stage 6 rollback.

## JE-specific reference warning
Official JE guide says name-year citations, alphabetized cited references, full DOI links where available, and only published/accepted items in the reference list. The file includes government reports, working papers, RIETI/OSIPP research papers and online institutional documents; their *publication/citation status and correct treatment* must be checked in the final rendered reference list. The existing `plainnat` style is not automatically a certification of full JE style compliance. Do not delete substantive citations without editorial assessment.

**Verdict:** Conditional on ordinary full bibliography/style preflight; no newly discovered fabricated reference or direct theorem absorption.

## Sources
- JE official instructions (accessed 2026-10-29 JST): https://link.springer.com/journal/712/submission-guidelines
- JE publication options (accessed 2026-10-29 JST): https://link.springer.com/journal/712/how-to-publish-with-us
- PR #73 head (read 2026-10-29 JST): 7308052c0f6299597eba4605ad9b4906966eb8f8
- Workflow v1.1 tag: 488e5ab06c207909296a7564eaf9066f7f94319c

