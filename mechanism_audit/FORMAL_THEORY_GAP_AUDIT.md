# Formal-Theory Gap Audit

## Verdict

**B2 — STRUCTURAL FORMAL-THEORY GAP.**

This is not a claim that no related formal theory exists. Related formal theory is extensive. The finding is narrower: the screened literature does not provide a near-exact formal model combining all layers of the target architecture.

## What is already formalized

### Search and intermediation

Gehrig (1993) formalizes heterogeneous agents choosing between costly decentralized search and intermediary-facilitated exchange. Therefore a claim that intermediaries reduce search friction is standard.

### Competition among matchmakers

Caillaud & Jullien (2003) model imperfect price competition between two intermediation service providers with indirect network externalities, multihoming/nonexclusive services, and identity/usage-based pricing. This is the strongest generic threat to any claim about “multiple hubs + matching externalities.”

Damiano & Li (2008) model competing matchmakers that price to sort heterogeneous participants when match value exhibits complementarity. Therefore heterogeneous productive matching plus competing matchmakers is also not new in itself.

Sülzle (2009) studies duopolistic competition between independent and collaborative B2B matching marketplaces with indirect network externalities and exclusive/multihoming users.

Rochet & Tirole (2003) provide the broader two-sided-platform competition framework.

### Innovation intermediation / TTOs

Hoppe & Özdenören (2005) develop a formal intermediary between creators and users of inventions; intermediary screening expertise economizes on evaluating profitable versus unprofitable inventions.

Macho-Stadler, Pérez-Castrillo & Veugelers (2007) formalize TTO reputation and pooling as responses to asymmetric information about invention quality.

Calcagnini et al. (2019) embed technology transfer in an endogenous matching model between TTOs and innovative firms; outcomes depend on TTO search costs and firms' project advertising.

Banal-Estañol, Macho-Stadler & Pérez-Castrillo (2018) formalize two-sided matching between heterogeneous academics and firms, including complementarity in collaboration value.

### Intermediaries and downstream firms

Ichihashi (2021) formalizes competition between data intermediaries that collect data from consumers and sell it to downstream firms. This is an important precedent for an intermediary layer feeding a downstream industry, but the mediated object is data rather than productive collaboration, and public/regional provision is absent.

## Empirically observed cooperation–competition duality

Feser (2023) documents a distinct dualism of cooperation and competition among Berlin innovation/creativity labs. Competition is associated with a growing number of similar providers, lack of transparency and competition for clients/utilization; cooperation is motivated by specialization, learning, missing resources/specific knowledge and joint projects. This supplies a direct institutional motivation for the duality but not a formal equilibrium model of public hub provision.

## Exact target architecture versus prior art

The intended theory requires all of the following in one model:

1. **multiple intermediary providers**;
2. **public/decentralized regional choice of provision**;
3. **heterogeneous complementary agents and productive collaboration formation**;
4. **overlapping intermediation opportunities** that create cannibalization;
5. **cross-hub partner-pool/referral spillovers** that can create positive interaction;
6. **beneficiary firms subsequently compete strategically downstream**;
7. **firm profits are endogenous**;
8. **regional welfare feeds back into public provision incentives**;
9. the central object is whether public hub provision is strategically complementary or substitutable.

No screened paper contains this full architecture or yields it by only relabeling variables, imposing symmetry, or adding one fixed cost.

## Why this is B2, not B3

The closest building blocks are mature and close enough that a broader search could still uncover an absorption route. In particular:

- M2 matching is saturated by search/matching theory;
- M6 shared-pool externality is adjacent to indirect-network-externality platform theory;
- public provision is adjacent to local-public-input/infrastructure theory;
- downstream firm competition is standard IO.

The candidate contribution must therefore be the **interaction of these layers**, not any individual component.

## Strongest referee objection

> “A co-creation hub is just a matching intermediary. Caillaud–Jullien, Damiano–Li and platform theory already tell us how competing matchmakers behave; adding a government owner and downstream firms is cosmetic.”

### Response at this stage

The objection succeeds against any claim based on matching, network effects, multihoming, heterogeneous matching, or intermediary competition alone. It does **not yet absorb the target architecture** because the closest matching models have private intermediary pricing/market-structure objectives and terminate at match formation or platform participation. The target public-policy game makes regional governments choose productive-matching infrastructure based on the downstream market consequences of matches, while another jurisdiction's intermediary simultaneously expands the effective partner pool and substitutes for local intermediation. Those feedbacks are absent from the screened closest models.

This response is mathematical/architectural rather than institutional-label based, so B2 is permissible. It is not strong enough for B3.

## One-page specialization standard

Top threats are audited separately. None passed the one-page test for the **full** target model. Several do absorb individual submodules, which must be cited as prior art in any eventual paper.

## Routing

Because a defensible direct primitive was selected and the gap is structural rather than application-only:

**GO TO MINIMAL FORMAL MODEL CONSTRUCTION.**

The next model must be designed to expose exactly which architectural link creates any new theorem. If the full game collapses to a standard competing-matcher/public-input corollary, return NO-GO rather than adding features.