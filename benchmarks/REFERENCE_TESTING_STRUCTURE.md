# Reference Testing Directory Structure

This document describes the directory structure for Python-based
numerical reference testing of the `batch_chem` driver. The directories and
files shown here will be added as the reference-testing infrastructure is
implemented.

The current per-benchmark coverage and remaining work are tracked in
[`REFERENCE_TEST_STATUS.md`](REFERENCE_TEST_STATUS.md).

## Source tree

```text
alquimia/
├── benchmarks/
│   ├── CMakeLists.txt
│   ├── tools/
│   │   ├── batch_results.py
│   │   ├── run_batch_reference_test.py
│   │   └── generate_batch_references.py
│   │
│   └── batch_chem/
│       ├── CMakeLists.txt
│       ├── *.cfg
│       ├── *-pflotran.in
│       ├── *-crunch.in
│       ├── *.dat
│       ├── *.dbs
│       │
│       └── references/
│           ├── manifest.json
│           ├── pflotran/
│           │   ├── general-reaction.csv
│           │   ├── calcite-vf.csv
│           │   ├── isotherms.csv
│           │   ├── ion-exchange-valocchi.csv
│           │   └── surface-complexation-2.csv
│           └── crunchflow/
│               ├── general-reaction.csv
│               ├── calcite-vf.csv
│               ├── isotherms.csv
│               ├── ion-exchange-valocchi.csv
│               └── surface-complexation-2.csv
│
├── unit_tests/
│   └── python/
│       ├── test_batch_results.py
│       └── fixtures/
│           ├── pflotran-observation.tec
│           ├── crunchflow-timeseries.tec
│           ├── expected.csv
│           └── malformed.csv
│
├── drivers/
│   ├── DriverOutput.c
│   ├── DriverOutput.h
│   └── input_util.c
│
└── cmake/
    └── Modules/
        └── add_alquimia_batch_chem_benchmark.cmake
```

## Responsibilities

- `benchmarks/tools/batch_results.py` provides the small reusable library for
  reading CSV files, comparing numerical results, and parsing native PFLOTRAN
  and CrunchFlow ASCII output.
- `benchmarks/tools/run_batch_reference_test.py` is the routine CTest entry
  point. It runs `batch_chem` and compares its CSV output with a committed
  reference.
- `benchmarks/tools/generate_batch_references.py` is a maintainer-only command
  that runs standalone chemistry engines and creates or deliberately updates
  normalized references.
- `benchmarks/batch_chem/references/manifest.json` defines reference cases,
  native-to-Alquimia column mappings, numerical tolerances, engine revisions,
  input checksums, and other provenance.
- `benchmarks/batch_chem/references/` contains reviewed, normalized CSV
  reference solutions. Raw engine logs and transient native output do not
  belong in this directory.
- `unit_tests/python/fixtures/` contains only small synthetic inputs used to
  test parsing and comparison behavior. Full simulation results belong with
  the benchmark references.

## Generated build tree

CTest output and logs should be written under the build tree rather than the
source checkout. Each test receives its own results directory.

```text
build/
└── benchmarks/
    └── batch_chem/
        └── results/
            ├── general-reaction-ac-pflotran/
            │   ├── actual.csv
            │   └── batch_chem.log
            ├── calcite-vf-pc-pflotran/
            │   ├── actual.csv
            │   └── batch_chem.log
            └── ...
```

Generated CSV files, native engine output, and logs remain untracked in the
build tree. Only reviewed normalized reference CSV files and their manifest are
stored in source control. Reference files must be updated through an explicit
maintainer command and must never be regenerated automatically by CTest or CI.

## Refreshing references

From the repository root, references can be refreshed explicitly with:

```bash
PYTHONPATH=benchmarks/tools python3 \
  benchmarks/tools/generate_batch_references.py \
  --manifest benchmarks/batch_chem/references/manifest.json \
  --pflotran /path/to/pflotran \
  --crunchflow /path/to/crunchflow \
  --update
```

Use one or more `--case CASE_NAME` arguments to update only selected cases.
The command runs each native engine in a temporary directory and records the
input checksum and executable name in the manifest. Review both the normalized
CSV changes and provenance before committing an update.

The initial CTest integration performs numerical comparisons for four
PFLOTRAN cases (general reaction, ion exchange, isotherms, and surface
complexation) and three CrunchFlow cases (general reaction, ion exchange, and
isotherms). Calcite evolution and CrunchFlow surface complexation remain smoke
tests because their current interface results do not agree sufficiently with
standalone output; the latter currently produces non-finite values.

CTest benchmark names end in `__smoke` or `__engine_reference`, making the
test type visible in the standard per-test summary. The same classification is
stored in CTest's `LABELS` property, so either group can be selected directly:

```bash
ctest -L smoke
ctest -L engine-reference
```

Unit tests use the `unit` label and do not carry either benchmark suffix.
