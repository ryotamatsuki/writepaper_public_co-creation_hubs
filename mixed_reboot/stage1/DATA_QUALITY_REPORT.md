# Stage 1R Data Quality Report

## Diagnostic panel

- Prefectures: **15**
- Years: **2024–2026**
- Region-year observations: **45**
- Adjudicated architecture observations: **42**
- `UNCERTAIN`: **3** (Kagawa P boundary)
- Adjudicated distribution: `HP=28`, `H=14`, `P=0`, `0=0`
- H component among adjudicated observations: **42/42 = 100%**

The key data-quality problem is therefore not random missingness. It is **architecture saturation induced by the prefecture-year union rule**.

## Double coding

Pass-A / Pass-B exact architecture agreement: **42/45 = 93.3%**.

This statistic must not be interpreted as strong measurement validity. Both passes identify H in every adjudicated region-year. H has zero cross-sectional variation, so Cohen's kappa for H is uninformative/undefined as a meaningful reliability measure. The three disagreements concern whether Setouchi-i-Base's broad extra-regional mission is enough for strict P coding.

The coding exercise shows that the unit is reproducible enough to diagnose its failure: the system-level union makes H nearly automatic and absorbs program-level P-only variation into HP.

## Zero regime

Five deliberately chosen candidate-zero systems (Akita, Aomori, Iwate, Mie, Yamanashi) were inventoried. None survived even the weaker `0_P` test in the pilot. `0_S` is still less plausible because universities, chambers, banks, public support organizations and sectoral intermediaries provide material intermediation even where no famous hub exists.

This does not prove a zero never exists in Japan. It shows that the specified prefecture-year definition cannot generate a credible zero at acceptable search cost without arbitrary exclusions.

## Regional-state measurement

Strongly measurable families:

- actor density;
- industrial structure;
- research capacity;
- accessibility.

Critical weak families:

- local fragmentation;
- independent baseline external connectivity.

NISTEP provides prefecture-level same-prefecture versus outside-prefecture university-industry collaboration. This is useful institutional evidence, but the two shares partition the same collaboration total and therefore cannot by themselves identify two independent state dimensions. The same-prefecture share is also a university-industry linkage measure rather than a general local-network-fragmentation measure.

A second theoretically cleaner local-network measure (within-prefecture co-patent/project-network density) and a second external-connectivity measure were **not validated** as reproducible prefecture-year public-data variables in this stage.

## Funding/governance

Governance form and decision-right proxies are observable for a meaningful subset: concession/PFI, ordinary outsourcing, designated management, multi-stakeholder committee and modular program contracting all appear.

Exact annual funding shares, H/P earmarks and decision rights are much less complete. This is a tractable missing-data problem, not the decisive Stage 1R failure.

## Participant geography

Directly structured geography is strong for NEXs Tokyo and TIB, enrichable for E:N BASE, and sparse/non-harmonized for most conventional hubs. Missingness is non-random by program type.

## Transition quality

At least four dated institutional changes are observable, but the prefecture-year binary union often records them as `HP -> HP` even when a distinct P module is added. The proposed unit therefore erases precisely the portfolio change the future theory wants to explain.

## Kill tests

| Kill test | Result |
|---|---|
| H/P coding unstable | not the main failure |
| P-only disappears / architecture effectively H-HP | **FAIL TRIGGERED** |
| credible 0 identifiable | **FAIL TRIGGERED** |
| local fragmentation distinct from external connectivity | **FAIL TRIGGERED** |
| governance almost unobservable | no |
| source coverage wholly incomparable | no, but uneven |
| regimes only branding | no; functional distinctions exist at program level |

## Stage 1R internal verdict

**STAGE 1R FAIL**

Diagnosis: **MEASUREMENT-UNIT FAILURE / AGGREGATION-INDUCED ARCHITECTURE SATURATION**.

Under the combined protocol, Stage 2R-A and Stage 2R-B must not be executed after this verdict.
