# Stage 4R-MC-G Model Contract

Use `SELECTED_STATIC_GAME.md` unchanged unless a mathematical inconsistency forces rollback to Stage 3.

Solve or recover B1, B2, B3 and G from the same primitives.

## Kill order
1. **Solvability:** derive provider regions, private optimum, public FOCs/SOCs/KKT and existence.
2. **Matching specificity:** test scalar/characteristic absorption using actual equilibrium results.
3. **Private strategic relevance:** compare B3 with G.
4. **Public strategic relevance:** show the second government changes a meaningful equilibrium/welfare object beyond mechanical replication.
5. **Full-model-only result:** at least one nontrivial result must be unavailable in every B1-B3 benchmark.
6. **Welfare wedge:** compare regional Nash with social planner.

Gate 5 failure implies `NO-GO`.

Use SymPy for exact algebra and then numerical counterexample search. Treat clipping, no-use, price corners and public investment boundaries explicitly.

Do not add dynamics, referral, Cournot, taste shocks, arbitrary cross terms, public pricing, endogenous metro quality or partner entry as rescue features.