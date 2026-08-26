# CHANGELOG

<!-- version list -->

## v2.0.1 (2026-08-26)

### Bug Fixes

- **build**: Guard the second all-feature read for split variants
  ([`9bea8f6`](https://github.com/alba1735/tRNAgraph/commit/9bea8f6b3bad4af25694f46d8f42c87224d3d224))

- **build**: Preserve on-disk count file format when filtering splits
  ([`46fa1ad`](https://github.com/alba1735/tRNAgraph/commit/46fa1adf53237be23f51b327bfb4d9a5de248f5b))

### Continuous Integration

- **release**: Pin semantic-release and gitpython
  ([`57597b1`](https://github.com/alba1735/tRNAgraph/commit/57597b1f9dec6cbf88cbc7afa76198c3959ac893))

### Documentation

- **build**: Drop the stale norm_allfeatures layer from the merge docstring
  ([`5b74b14`](https://github.com/alba1735/tRNAgraph/commit/5b74b14bfc959dab76d91e478f1e5f1718762540))

- **roadmap**: Log the tools test delete footgun and console parity
  ([`1390925`](https://github.com/alba1735/tRNAgraph/commit/1390925d2f712f5dbe559ce907276bb1543a8c8f))

- **roadmap**: Various ideas added
  ([`a6e3a13`](https://github.com/alba1735/tRNAgraph/commit/a6e3a138a3a310a6028687cbb74bccb752eae879))


## v2.0.0 (2026-08-25)

### Bug Fixes

- **audit**: Route pretRNAcoverage to the normalization-matched reference
  ([`2b29f84`](https://github.com/alba1735/tRNAgraph/commit/2b29f8431924dfc14562ebfbbddb9d9f4506b34b))

- **cli**: Capture real errors in logs and stop rich from polluting them
  ([`a2fffd5`](https://github.com/alba1735/tRNAgraph/commit/a2fffd594f23985d55186ba74bb30439eb4bcec3))

- **coverage**: Always exclude multi-mapped reads to match tRAX
  ([`f3db434`](https://github.com/alba1735/tRNAgraph/commit/f3db434cc0a165f1d359797b2bb7ff246cd8ec42))

- **coverage**: Number the loci alignment before adding its edge margin
  ([`51ab27b`](https://github.com/alba1735/tRNAgraph/commit/51ab27b3ac9df7a21cbe7d17798e638d326a0f53))

- **docs**: Update improvements_from_trax and roadmap for bug fixes and clarifications
  ([`1b6c07e`](https://github.com/alba1735/tRNAgraph/commit/1b6c07e30525305a061349d9c2050f54a6f4d5c6))

- **docs**: Update roadmap with detailed bug fixes and validation results for remote pipeline
  failures
  ([`295bb5c`](https://github.com/alba1735/tRNAgraph/commit/295bb5ca2eda1f9cdda8a130d694991868fc5174))

- **makedb**: Validate organism mode and unify the Sprinzl position tables
  ([`5b57713`](https://github.com/alba1735/tRNAgraph/commit/5b5771365fc284272136a2997397bdf932f73312))

### Continuous Integration

- Derive the release version automatically instead of forcing a minor
  ([`d32b9e4`](https://github.com/alba1735/tRNAgraph/commit/d32b9e41da34e2fa3ed7ffe390d85c53bb0dfd12))

- Force the release bump level and fix dev pre-releases
  ([`31dbebd`](https://github.com/alba1735/tRNAgraph/commit/31dbebd01a2d88a2145f2cd4a542c7871bcc5f2d))

### Documentation

- Add citation metadata and keep its version in step with releases
  ([`9427781`](https://github.com/alba1735/tRNAgraph/commit/9427781d8095e7d3f6a8ca7b717fd61990f15881))

- Document the read basis and the coverage specificity partition
  ([`cca3ca0`](https://github.com/alba1735/tRNAgraph/commit/cca3ca047b8d2cefe6e82a147140f8087ad5bb5c))

- Rewrite the tRAX migration guide and drop internal-tool references
  ([`a981d10`](https://github.com/alba1735/tRNAgraph/commit/a981d100e69e177138b6e2028b2fbcc4952b4e3a))

- **roadmap**: Record hg38 validation outcome and narrow the open divergences
  ([`2b29f84`](https://github.com/alba1735/tRNAgraph/commit/2b29f8431924dfc14562ebfbbddb9d9f4506b34b))

### Features

- **coverage**: Make tRAX's read-specificity partition selectable and plot it
  ([`2b5cfbd`](https://github.com/alba1735/tRNAgraph/commit/2b5cfbd57fdd04d0fa5b6d9d1169dbece0e11f15))

- **release**: Enhance version bump options and improve semantic-release configuration
  ([`dfad2a6`](https://github.com/alba1735/tRNAgraph/commit/dfad2a6bb6ad2a889499b64ec4b19bdc1253c38b))

### Refactoring

- **cli**: Report the dev channel as prerelease rather than beta
  ([`67d0790`](https://github.com/alba1735/tRNAgraph/commit/67d0790d8f1ec6030b5e79afdb2f1b1cdcc1fd91))

- **graph**: Default every graph type to unique counts
  ([`e78fb38`](https://github.com/alba1735/tRNAgraph/commit/e78fb38183a20fd09bf1a699e0bd3e6f84d4df0a))


## v1.8.0 (2026-08-21)

### Bug Fixes

- **build**: Only delete split-BAM directories a run actually created
  ([`321c9e8`](https://github.com/alba1735/tRNAgraph-private/commit/321c9e8666698a632fb2262e043ec915c1c681b8))

### Chores

- **release**: 1.8.0
  ([`09ffb7d`](https://github.com/alba1735/tRNAgraph-private/commit/09ffb7d822c077d94291a8938d97623470d00079))

### Documentation

- Fix README/CLI reference inaccuracies; reword provenance terminology
  ([`137abc2`](https://github.com/alba1735/tRNAgraph-private/commit/137abc2b6764eccd94dcb6c5af1551d4ded1732c))

- **cli**: Clarify --maxmismatches flag description
  ([`86a275a`](https://github.com/alba1735/tRNAgraph-private/commit/86a275a90e62241b8c5408a15ac43aaff1486a5f))

### Refactoring

- Remove --maponly flag from tools test
  ([`eaae637`](https://github.com/alba1735/tRNAgraph-private/commit/eaae63738b2e1f8b844bd3b5bd1582ae8da87833))

- Remove --nofrag flag and its dead code paths
  ([`968f95a`](https://github.com/alba1735/tRNAgraph-private/commit/968f95abf65d60cee2234232345e9b752740a4a1))

- Remove broken --nosizefactors flag
  ([`5b2afd6`](https://github.com/alba1735/tRNAgraph-private/commit/5b2afd68f07ae09757bffa42ed91eb05c845b90c))

- Remove dead graphsdir creation from AnalysisPipeline and MapSamples
  ([`d9d302b`](https://github.com/alba1735/tRNAgraph-private/commit/d9d302b31e7b0c3495a63ce7111dbe10a1beb656))

- Rename --dumpother flag to --filterother
  ([`1f7cfc3`](https://github.com/alba1735/tRNAgraph-private/commit/1f7cfc3182cf894ef5d8e4f9a5bb5b383fdbd28b))

- Rename --lazy flag to --force-remap
  ([`0931f39`](https://github.com/alba1735/tRNAgraph-private/commit/0931f3961bca42ab5267d0749ab5a2ff17c2e3c8))

- Rename --mincoverage to --minfeaturereads, scoping it to the VST fit only
  ([`018e214`](https://github.com/alba1735/tRNAgraph-private/commit/018e2148cd6261bef7a0dffc6b9d7bbe974727a8))

- Rename --uniqueonly flag to --filtermultimapped
  ([`28d1fd5`](https://github.com/alba1735/tRNAgraph-private/commit/28d1fd56a05dea30aefbda7509140af034c73751))

### Testing

- Add reproducibility regression test for trnaloci classification ordering
  ([`3ee6189`](https://github.com/alba1735/tRNAgraph-private/commit/3ee6189111a5d84b70f41fd16b40bc15998eda90))


## v1.7.0 (2026-08-20)

### Chores

- **release**: 1.7.0
  ([`fffe004`](https://github.com/alba1735/tRNAgraph-private/commit/fffe004f10df87ca0db1d535a03157bb27ec1776))

### Features

- **build**: Add phase-tracked progress reporting to toolsCountReads and the build pipeline
  ([`d344199`](https://github.com/alba1735/tRNAgraph-private/commit/d344199bb3aa61fe1bf24b1ab5d1e9062cb8d86c))

- **graph**: Extend rich progress reporting to graph generation
  ([`9d64d11`](https://github.com/alba1735/tRNAgraph-private/commit/9d64d117bdadefda8dcb7e3c1acebafcdff96dc1))

- **logging**: Add centralized logging configuration
  ([`8fe1bc4`](https://github.com/alba1735/tRNAgraph-private/commit/8fe1bc4386eb407c023e93a14964d85935f96db8))

- **map**: Report mapping progress via toolsTG.progress_iterator()
  ([`8b1d33e`](https://github.com/alba1735/tRNAgraph-private/commit/8b1d33e61947cbc6c18f5d9fd04ccc1fa453d088))

### Refactoring

- **logging**: Apply centralized logging across all modules
  ([`ab91c35`](https://github.com/alba1735/tRNAgraph-private/commit/ab91c3508968dcaaea0663211f8050f66048f86b))


## v1.6.0 (2026-08-20)

### Chores

- **release**: 1.6.0
  ([`6125bd6`](https://github.com/alba1735/tRNAgraph-private/commit/6125bd63f37970a0238d59a0d78f7e74844494de))

### Features

- Resolve grouping column with fallback across visualizer functions
  ([`baa36f7`](https://github.com/alba1735/tRNAgraph-private/commit/baa36f75b23c590bb2a91d6bdb9ec4cf828a5630))

- **merge**: Add --force option and build-provenance conflict validation
  ([`da5c423`](https://github.com/alba1735/tRNAgraph-private/commit/da5c423fd63e364058532c6621353c3c80b129e4))

- **validation**: Add pydantic schemas for graph filter, colormap, metadata, and pairs files
  ([`2d08ddd`](https://github.com/alba1735/tRNAgraph-private/commit/2d08dddadf0775a2ea4e25c4e39dc52de4d491bb))


## v1.5.1 (2026-08-19)

### Chores

- **release**: 1.5.1
  ([`d8133b2`](https://github.com/alba1735/tRNAgraph-private/commit/d8133b2080610817334f189cf61d2b24e21bc8ef))

### Documentation

- Add improvements_from_trax.md documenting statistical methodology differences
  ([`0d4467b`](https://github.com/alba1735/tRNAgraph-private/commit/0d4467b3c379f443ec24d923c9a5b1696a941e1f))

### Features

- **build**: Validate metadata sample names against coverage data
  ([`7ed5497`](https://github.com/alba1735/tRNAgraph-private/commit/7ed5497c9754f1928b372c0146dbe9c2a696298a))

- **plots**: Support non-tRNA features in correlation matrix plots
  ([`8adc91f`](https://github.com/alba1735/tRNAgraph-private/commit/8adc91fe673563f7327b4f93634f47c84f43c4fb))

- **trim**: Add colormap support for trim-stats plot
  ([`b7cdae3`](https://github.com/alba1735/tRNAgraph-private/commit/b7cdae3640d829eb034970b8172798b2eab83755))


## v1.5.0 (2026-08-19)

### Bug Fixes

- **build**: Fix VST hang on larger sample counts
  ([`e3d39c2`](https://github.com/alba1735/tRNAgraph-private/commit/e3d39c255a69346762678ec67cdb4de81b7da641))

- **cluster**: Store sample cluster results in obsm instead of uns
  ([`7b10801`](https://github.com/alba1735/tRNAgraph-private/commit/7b1080160225663f0543587fdf90f418524294ab))

### Chores

- **release**: 1.5.0
  ([`f747a96`](https://github.com/alba1735/tRNAgraph-private/commit/f747a960718731a4c9f269a9f98537b52a4b23d0))

### Refactoring

- **log2fc**: Replace manual t-test with native PyDESeq2 GLM
  ([`5656026`](https://github.com/alba1735/tRNAgraph-private/commit/565602619af6ed763a18b4afb5ef2447e1f095b1))


## v1.4.0 (2026-08-19)

### Chores

- **release**: 1.4.0
  ([`01cf300`](https://github.com/alba1735/tRNAgraph-private/commit/01cf30098f8cda2860f3e0fc7dc47b819ca057ca))

### Documentation

- Fix notes in advanced_usage/data_structure/roadmap
  ([`d7c366d`](https://github.com/alba1735/tRNAgraph-private/commit/d7c366d0156226bcbc1375699f407986b99500dd))

- **readme**: Clarify preprocessing/analysis/graphing flowchart
  ([`4f370f0`](https://github.com/alba1735/tRNAgraph-private/commit/4f370f031242967835cb42ac52493e6d32a23bc6))

### Features

- **build**: Add all-feature-controlled normalization; improve PCA plots
  ([`f3e6e4b`](https://github.com/alba1735/tRNAgraph-private/commit/f3e6e4bc6cebc362b4df02f7d6f8b4e196ef16ab))

- **build**: Refactor AnnData handling for read-length split variants
  ([`aee76a5`](https://github.com/alba1735/tRNAgraph-private/commit/aee76a5f8dd19137704f468663f80c8492cb30a3))

- **graph**: Enhance volcano/PCA/heatmap plots and documentation
  ([`fc7dd25`](https://github.com/alba1735/tRNAgraph-private/commit/fc7dd254e0b886f3a787df74d2b0d50716541f08))

### Refactoring

- Streamline BAM splitting and analysis CLI/internal modules
  ([`724c73c`](https://github.com/alba1735/tRNAgraph-private/commit/724c73c0a5d7ae86fb34e1a83c34bc61d10b9723))


## v1.3.2 (2026-08-12)

### Chores

- **release**: 1.3.2
  ([`4355974`](https://github.com/alba1735/tRNAgraph-private/commit/4355974d999988958691187cf5eb07a9e5a5ec5d))

### Features

- **graph**: Expand CLI flags and add Venn diagram support prototype
  ([`48545a6`](https://github.com/alba1735/tRNAgraph-private/commit/48545a65f875a9374fff719656bfb45662ef63a2))

- **testsuite**: Add colormap.json fixture and migrationTool integration dependency
  ([`c597bd7`](https://github.com/alba1735/tRNAgraph-private/commit/c597bd79971caef05e1b36a5b4a119dd6243e97b))


## v1.3.1 (2026-03-25)

### Bug Fixes

- Add tqdm dependency and switch default VST method to log1p
  ([`8ca1989`](https://github.com/alba1735/tRNAgraph-private/commit/8ca1989b82fb24759f74c69a5c9a12a9a9e94ee4))

- **plots**: Split bar plot into its own count module
  ([`76c2460`](https://github.com/alba1735/tRNAgraph-private/commit/76c246027dfc0be2b9dd5f1617a3d0887d0fb342))

### Chores

- **release**: 1.3.1
  ([`a7fc21d`](https://github.com/alba1735/tRNAgraph-private/commit/a7fc21da2114acee05da603f2799720e2cb0671a))

### Documentation

- **readme**: Minor updates
  ([`308ad4a`](https://github.com/alba1735/tRNAgraph-private/commit/308ad4a677f5a755c69709587039c530596ab25f))


## v1.3.0 (2026-03-16)

### Chores

- **release**: 1.3.0
  ([`5de5fca`](https://github.com/alba1735/tRNAgraph-private/commit/5de5fcad7ff5280adbbe175c327cc157726aec4d))

### Features

- Standardize metadata format to include fastq paths; remove trim's -o flag
  ([`feebf93`](https://github.com/alba1735/tRNAgraph-private/commit/feebf93272fe86611931c038b4f2320f5ceca6d3))

- **build**: Add tRNA-specific size factor calculation and multiple VST methods
  ([`835b56d`](https://github.com/alba1735/tRNAgraph-private/commit/835b56d782051a6ad2763418495934329a5a2b2c))

- **split**: Refactor BAM splitting into a dedicated tool with configurable output dirs
  ([`82c69d2`](https://github.com/alba1735/tRNAgraph-private/commit/82c69d2be99a856ff4f7d6256f2ddee1aaf21fa1))


## v1.2.0 (2025-12-10)

### Bug Fixes

- Correct mapping and trimming issues to better match tRAX output
  ([`e10350c`](https://github.com/alba1735/tRNAgraph-private/commit/e10350c18d542cc1028e09f6bbb8208829b76fc9))

- **env-check**: Handle version-operator requirements; widen glibc compatibility window
  ([`fbda867`](https://github.com/alba1735/tRNAgraph-private/commit/fbda867f1bc7c27cd8b2f12f4fafb6c821d57829))

- **plots**: Capture original tRNA sequence positions in seqlogo generation
  ([`93b971a`](https://github.com/alba1735/tRNAgraph-private/commit/93b971a0605402e83edabb9e6b3e59c45f47afb9))

### Chores

- **release**: 1.2.0
  ([`b5ee928`](https://github.com/alba1735/tRNAgraph-private/commit/b5ee92847db8f36b2621e98abde90ae4b7212796))

### Documentation

- Split README into advanced_usage/cli_reference/data_structure guides
  ([`2043c89`](https://github.com/alba1735/tRNAgraph-private/commit/2043c8927d947d0696547e44431a60fb43f1d513))

- Update CLI docs and commands for the analyze workflow
  ([`eaf6cf5`](https://github.com/alba1735/tRNAgraph-private/commit/eaf6cf5f8b30675092169dc23a8644e4e04467e5))

### Features

- Add comprehensive tRAX documentation and align tool output with tRAX
  ([`82902e8`](https://github.com/alba1735/tRNAgraph-private/commit/82902e8455660acd97cedbda1e88e7c5432709c7))

- Introduce adataBuild/adataCluster/adataGraph/adataMerge modules
  ([`dca2673`](https://github.com/alba1735/tRNAgraph-private/commit/dca2673861ef33642acd80039722f51aa5b7e1f2))

- **build**: Add DESeq2 dispersion fit type option
  ([`a56bb79`](https://github.com/alba1735/tRNAgraph-private/commit/a56bb798442bb7e023e4c206daf4f5e2948a6ea5))

- **cli**: Add 'split' command for read-length BAM splitting
  ([`b026aae`](https://github.com/alba1735/tRNAgraph-private/commit/b026aae1c63e6c793ca8a3b9910702e17190100c))

- **cli**: Add environment validation checks with a skip option
  ([`078859d`](https://github.com/alba1735/tRNAgraph-private/commit/078859d51597c18c1e05459ee1e5df41e97226e6))

- **env-check**: Add umi_tools to requirements and validation
  ([`c1dde48`](https://github.com/alba1735/tRNAgraph-private/commit/c1dde48bf3ce3b01a1375108566fd14ad5e9fbd8))

- **testsuite**: Add tqdm progress bars
  ([`16ccae4`](https://github.com/alba1735/tRNAgraph-private/commit/16ccae45de52d5421fbf069681a67c5fb02d0498))

### Refactoring

- Reorganize CLI, docs, and coverage/mapping logic for clarity
  ([`2c59d8a`](https://github.com/alba1735/tRNAgraph-private/commit/2c59d8ab98dc8985883e65fe32c078cb66def7ec))

- Reorganize tRNAgraph scripts and assets
  ([`9345292`](https://github.com/alba1735/tRNAgraph-private/commit/934529202b8d269ab0c02ddd49ff924a241b6017))

- Restructure project as an installable Python CLI package
  ([`35bf81b`](https://github.com/alba1735/tRNAgraph-private/commit/35bf81bdc1fe7ab6986f12a092ce5fda0b8c8d7e))

- **cli**: Restructure CLI and module imports for clarity and performance
  ([`45c4bfe`](https://github.com/alba1735/tRNAgraph-private/commit/45c4bfe07912797547663837d05072ae739512fd))


## v1.1.0 (2025-11-21)

### Features

- Add experimental trimming, mapping, and track hub tooling
  ([`09da375`](https://github.com/alba1735/tRNAgraph-private/commit/09da3753ad4ad2e7201e96de0c5d462f7e9a2d2d))

- Implement lazy loading for plot and tool modules
  ([`6267c9d`](https://github.com/alba1735/tRNAgraph-private/commit/6267c9da4692792dd48aacca28dde0f40ea81e08))

- Improve MakeDB docs; refine MapReads/demoPipeline thread handling
  ([`680c257`](https://github.com/alba1735/tRNAgraph-private/commit/680c2572c7d8cb92e4479dca1a08e4fbbc4a17de))

- **plots**: Add experimental legacy plotting scripts
  ([`60a9291`](https://github.com/alba1735/tRNAgraph-private/commit/60a929180dbc3966cc5eac3b6e186925f5e7dc3c))

- **testsuite**: Add preflight testing functionality
  ([`1cb740a`](https://github.com/alba1735/tRNAgraph-private/commit/1cb740a30478e2f7bb7abe89a7b6461fbf25ae92))

### Refactoring

- **cli**: Migrate from argparse to Typer with type hints
  ([`02cca0b`](https://github.com/alba1735/tRNAgraph-private/commit/02cca0b3221483a00385f76ed04dd28e5e0482eb))


## v1.0.1 (2025-05-27)

### Chores

- **deps**: Consolidate dependencies into requirements.yaml
  ([`0ea7e8b`](https://github.com/alba1735/tRNAgraph-private/commit/0ea7e8b80af46ba7bea9efb2bbbaa3d8ee36f4ae))

### Documentation

- **readme**: Add Zenodo DOI citation badge
  ([`4a1a7e7`](https://github.com/alba1735/tRNAgraph-private/commit/4a1a7e76e000658956b0e535d94fd3a1dd20c50b))

- **readme**: Fix incorrect conda environment file name
  ([`4d1c7a1`](https://github.com/alba1735/tRNAgraph-private/commit/4d1c7a11678cf12e14b214b5d54471657f00a627))

- **readme**: Fix typo
  ([`b9fc9fb`](https://github.com/alba1735/tRNAgraph-private/commit/b9fc9fb08db55baffb42c2bebcc6433f61b7c3bf))

### Features

- **plots**: Drop ACT from radar plots, add RNA mode to seqlogo, improve cluster color mapping
  ([`ef7f993`](https://github.com/alba1735/tRNAgraph-private/commit/ef7f993a15454170f9860ed167a30cdd4c2cb4dc))


## v1.0.0 (2025-01-16)


## v0.9.2 (2024-02-17)


## v0.9.1 (2024-02-07)

- Initial Release
