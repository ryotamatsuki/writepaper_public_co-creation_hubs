# Stage 12 Disclosure / Format Audit

Checked 2026-09-03. Only current public policies located on official pages are treated as requirements.

| Journal | Submission system | First-submission format | Data/code | AI / generative AI | COI/funding | Preprints | Other material issue |
|---|---|---|---|---|---|---|---|
| **JPET** | Wiley Research Exchange / Wiley Authors | Free Format; LaTeX main source + PDF accepted | Data sharing expected; Data Availability Statement required for Original Articles | Wiley: AIGC/LLM cannot be author; generative use beyond copyediting should be transparently disclosed; human accountability | Funding/standard declarations | Accepted | Current author page requests structured abstract; recheck live portal at Stage 13/14 |
| JRS | Wiley Authors | Free Format; double blind | Data sharing mandated; DAS required | Wiley publisher policy applies; no conflicting journal-specific rule located | Standard Wiley declarations | Accepted | Anonymize manuscript/title page separation |
| Journal of Economics | Editorial Manager | LaTeX encouraged/accepted; source + PDF; abstract 150–250; JEL | Research Data Policy/DAS in guidelines; code/material availability statements as applicable | LLMs not authors; generative LLM use documented in Methods/suitable section; AI-assisted copyediting alone need not be declared | Competing Interests required; funding acknowledgments | Springer policy applies | Source files required at submission/revision |
| RRS | Scholastica | Free format first; accepted papers LaTeX; double blind | Data/code statement optional | **AI Usage Disclosure required** | **Conflict of Interest required**; funding recommended | No contrary prohibition located | Cover letter must include spatial-element statement |
| REGION | Journal submission site | LaTeX accepted; normally <=10,000 words | Open-science/code encouraged | Journal-specific AI rule **NOT LOCATED** in checked pages | Standard scholarly declarations | No material blocker located | No submission/APC fees |

## Primary Stage-13 disclosure package

Prepare for JPET:

- Data Availability Statement. For this theory paper, state accurately that no external empirical dataset is required and identify the public repository/code used for symbolic/numerical verification if appropriate and stable.
- Funding declaration based on actual funding status.
- Conflict-of-interest declaration if required by live portal.
- ORCID in the submission system where requested.
- AI disclosure that truthfully describes generative-AI use in manuscript/research preparation if it went beyond AI-assisted copyediting. Do not under-disclose and do not invent tool use that did not occur.

## Reproducibility-policy fit

The repository already contains generated results, symbolic/numerical verification, tests and a full manuscript build. This exceeds the minimum evidence needed for a theory submission. Stage 13/14 should expose only the appropriate public/stable artifacts and should not alter generated economic outputs or automatically publish any private/sensitive material.