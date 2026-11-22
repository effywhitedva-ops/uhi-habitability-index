# Use of AI agents in UHI v6.1

**Agents.** Two custom software agents based on large language models took part in this version: **Shoggothsito** (implementation and execution) and **Aetherios** (technical review).

**What they did.**
- Shoggothsito translated the author's model decisions into Python, ran the Monte Carlo variants and the Saltelli–Sobol computations, wrote the test suite, produced the tables and figures from scripts, kept the run logs, and drafted this documentation.
- Aetherios audited the methodology and code against the cited literature and proposed corrections.
- Both assisted with literature search. No agent made a modelling decision on its own: each change was proposed, reviewed and approved by the author before it entered the code, and the approvals are recorded with dates in `CHANGELOG.md`.

**How their output was checked.**
1. With every new module switched off, the code reproduces the released v5.2 results on all 6,112 planets (maximum difference 9×10⁻¹⁵, the same difference the v5.2 code shows when rerun in this environment), and bit for bit on an in-memory subset.
2. Synthetic Earth, Venus and Mars at 10 pc behave as expected (Earth survives; Venus and Mars are vetoed).
3. Key quantities for GJ 273 b, Proxima Cen b and TRAPPIST-1 d were recomputed by an independent script and agree with the pipeline within Monte Carlo noise (≤ 1.4%).
4. The atmospheric-retention metric reproduces hand-computed values from Pass et al. (2025), Eqs. 10–11, within 0.01 for 14 planets.
5. The Sobol report was validated on the Ishigami function before use.
6. The automated tests (`tests/test_v61.py`, `tests/test_v6.py`) pass in full on the released code.

Literature values were checked against the original sources (DOI records and NASA ADS). The author is responsible for the scientific content. The agents are tools, not authors.
