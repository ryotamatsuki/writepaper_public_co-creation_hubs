# Stage 15 release disposition — source, identity, bibliography, and journal portal
**Date:** 2026-11-03 (JST)
**Article:** *Private Innovation Intermediation, Regional Incidence, and Public Access*.
**Target:** Springer Nature, *Journal of Economics* (JE).
**Scientific evidence baseline:** `24e241b0e3991a6e2ebca9394dd6d848e75bdc9a`.
**Last manuscript/bibliography-change commit:** `84d1a8f09ea46745370c11793f3cb3dbfed6bdaa`.
**Audited Stage 15 report/checklist HEAD:** `230743ce315f30f46c1277f98edfa7bf60ee8611`.

## Completed

- **Historical author details, source-checked:** The name Ryota Matsuki, affiliation *Independent Researcher, Matsuyama, Ehime, Japan*, Gmail correspondence address and ORCID in the JE title page agree with prior actual journal submissions/title pages. The prior submitted materials also record *no external funding* and *no competing interests*, but those are paper-specific declarations and cannot certify the present article. The full account is in `reviews/STAGE15_FINAL_REFERENCE_AND_AUTHOR_AUDIT_2026-11-03.md`. Avoid copying any phone/private residential address from historical records.
- **Thirty-one cited references:** Each cited key is mapped in the cited-item audit to a publisher, official institutional repository or primary record. Five detected factual metadata issues were corrected in `manuscript/references.bib` (Bank of Canada author order; Carpentier year-volume/pages; Eldebo–Hjelm published issue/year; 1988 d'Aspremont–Jacquemin DOI; 1986 Mankiw–Whinston DOI). Additionally, a DOI link for the `@techreport` Bank of Canada paper is now explicitly displayed. This is bibliography metadata QA, **not** an absence-of-prior-theorem proof: the previously documented unavailable complete Hua / Damiano–Li / Ostrizek–Sartori technical comparisons remain qualified.
- **Compilation and audit on exact new scientific presentation:** At `84d1a8f09ea46745370c11793f3cb3dbfed6bdaa`, GitHub Actions **all SUCCESS**: JE build run [41655446431](https://github.com/ryotamatsuki/writepaper_public_co-creation_hubs/actions/runs/41655446431), reproducibility run [41655446453](https://github.com/ryotamatsuki/writepaper_public_co-creation_hubs/actions/runs/41655446453), JE scope/formal core run [41655446418](https://github.com/ryotamatsuki/writepaper_public_co-creation_hubs/actions/runs/41655446418). A CI-produced preview artifact (ID **13240530081**) was downloaded, unpacked and validated. The actual CI PDFs are 34-page `manuscript.pdf` (SHA256 `64561176ca91ef2b6612f55b1a0c6981d5bdee294041c9c777976de94bc90e97`) and 8-page `ESM_1.pdf` (SHA256 `b669949dc19720a1f4c60da398ee5136fbd564587f6fdfc40869c3288d83d196`). The nested PREVIEW archive SHA256 is `c0a9c9bb234171db5ae6b2db53612f572962225db8416fb572aa466253b26c72`. The bibliography source matches the corrected CI file. No CI LaTeX warnings.
- **Full current report/checklist HEAD check:** For `230743ce315f30f46c1277f98edfa7bf60ee8611`, all three follow-up workflows also **SUCCESS**: JE build run [41655997251](https://github.com/ryotamatsuki/writepaper_public_co-creation_hubs/actions/runs/41655997251), reproduction [41655997206](https://github.com/ryotamatsuki/writepaper_public_co-creation_hubs/actions/runs/41655997206), scope [41655997337](https://github.com/ryotamatsuki/writepaper_public_co-creation_hubs/actions/runs/41655997337). Difference after scientific-source commit affects audit/checklist only.

## Hard blocker A — specific published journal portal is NOT live

The journal's own [submission guidelines](https://link.springer.com/journal/712/submission-guidelines) link to [Editorial Manager JEZN](https://www.editorialmanager.com/jezn/), which displays exactly: **"Important Message: Site under development. Do not use for live manuscript submission."**

The **Editorial Manager vendor documentation** [Site Under Development Message and "Go Live" Functionality](https://emhelp.editorialmanager.com/robohelp/Site_Under_Development_Message_and_Go_Live_Functionality.htm) states explicitly that this banner indicates **implementation mode and not accepting live manuscript submissions**. This supersedes earlier uncertainty that the banner might merely be a display warning. An official working submission destination must be confirmed with the journal/publisher. Do not submit to this test/implementation site.

The accessible Springer [journal contact page](https://link.springer.com/journal/712/contact-the-journal) links to [Springer Nature journal support](https://support.springernature.com/en/support/tickets/new?ticket_form=journal_contact_form).

### Proposed author-to-publisher inquiry (not sent)
Subject: *Journal of Economics* — live submission portal and mandatory charges

> Dear Springer Nature Journal Support,\n\nI am preparing a LaTeX research article for submission to the *Journal of Economics*. The journal's official "Submit your manuscript" link directs me to Editorial Manager (https://www.editorialmanager.com/jezn/), which currently displays the warning "Site under development. Do not use for live manuscript submission." Could you please provide the journal's currently active submission URL, or confirm when live submission will be available?\n\nPlease also confirm whether initial submissions require an anonymized manuscript, which file categories should be used for the editable main/supplementary LaTeX sources and `ESM_1.pdf`, and whether any mandatory submission fee, page charge or overlength charge applies when opting for traditional subscription publication (without open-access APC).\n\nThank you.

This draft is ready for the author to send **through the official support page**; no support request was submitted or reply obtained in this workflow. Journal standard subscription publication has no required OA APC, but zero submission/length fee is not separately guaranteed without publisher/portal confirmation.

## Hard blocker B — paper-specific attestations must still be true

Although the historical identity/contact details are established, the current manuscript's assertions **"No external funding"**, **"no competing interests"**, the employer versus independent work/resource/permission question, sole-author accountability and exclusive submission/non-overlap, and the completeness of substantive generative-AI use disclosure are **specific to this project and cannot be verified from other articles' declarations**. The listed claims are present in the current paper, but the author must affirm them or revise truthfully. The presence of paid public employment does not automatically disqualify independent affiliation, but is relevant to resource provenance, permissions and nonfinancial conflicts.

## Final Stage 15 decision

**PREFLIGHT COMPLETE; FINAL RELEASE HOLD.** The reviewed paper is academically and technically ready for a publisher-portal preflight; citation metadata repairs and compilation are complete to the stated scope. **Not yet submission-authorised**, because JE's linked EM destination is explicitly in implementation mode and the present-article author disclosures/permissions are not independently attested. No merge; PR #73 remains Draft; no journal submission. Do not relabel the attached CI `PREVIEW` files as final until those two hard gates close and immutable end-package SHA is captured.
