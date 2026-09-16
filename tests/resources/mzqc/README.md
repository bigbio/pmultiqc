# mzQC fixtures

The small JSON files in this directory are synthetic structural fixtures for unit and integration tests.

For the real prideQC validation, set `PRIDEQC_MZQC_FIXTURE_DIR` to the directory containing:

- `Natalia_TMT0_07_120m_1pt5.raw.mzQC`
- `Prosser_1004.raw.mzQC`
- `TDM_M1808_198.raw.mzQC`

The real fixtures are intentionally not copied into this source repository by the module implementation pass; they are local test artifacts and should remain canonical prideQC outputs.
