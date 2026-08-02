#!/usr/bin/env python3
"""Run one batch_chem benchmark and compare it with a CSV reference."""

import argparse
import json
import os
import shutil
import subprocess
import sys

from batch_results import ResultsError, compare_tables, read_csv_table


INPUT_SUFFIXES = (".cfg", ".in", ".dat", ".dbs", ".dbsx")


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--batch-chem", required=True)
  parser.add_argument("--benchmark-dir", required=True)
  parser.add_argument("--manifest", required=True)
  parser.add_argument("--case", required=True)
  parser.add_argument("--work-dir", required=True)
  args = parser.parse_args()

  with open(args.manifest, "r") as stream:
    manifest = json.load(stream)
  try:
    case = manifest["cases"][args.case]
  except KeyError:
    print("Unknown reference case: {}".format(args.case), file=sys.stderr)
    return 2

  os.makedirs(args.work_dir, exist_ok=True)
  for filename in os.listdir(args.benchmark_dir):
    source = os.path.join(args.benchmark_dir, filename)
    if os.path.isfile(source) and filename.endswith(INPUT_SUFFIXES):
      shutil.copy2(source, os.path.join(args.work_dir, filename))

  input_name = case["alquimia_input"]
  input_path = os.path.join(args.work_dir, input_name)
  output_path = os.path.join(
      args.work_dir, os.path.splitext(input_name)[0]+".csv")
  log_path = os.path.join(args.work_dir, "batch_chem.log")
  if os.path.exists(output_path):
    os.remove(output_path)
  with open(log_path, "w") as log:
    result = subprocess.run([args.batch_chem, input_name], cwd=args.work_dir,
                            stdout=log, stderr=subprocess.STDOUT)
  if result.returncode != 0:
    print("batch_chem failed with status {}; see {}".format(
        result.returncode, log_path), file=sys.stderr)
    return result.returncode or 1
  if not os.path.isfile(output_path):
    print("batch_chem did not create {}".format(output_path), file=sys.stderr)
    return 1

  reference_path = os.path.join(
      os.path.dirname(args.manifest), case["reference"])
  defaults = manifest.get("tolerances", {})
  try:
    actual = read_csv_table(output_path)
    reference = read_csv_table(reference_path)
    failures = compare_tables(
        actual, reference,
        rel_tol=float(defaults.get("rel_tol", 1.0e-6)),
        abs_tol=float(defaults.get("abs_tol", 1.0e-12)),
        time_abs_tol=float(defaults.get("time_abs_tol", 1.0e-9)),
        overrides=case.get("tolerances", {}))
  except (OSError, ResultsError, ValueError) as error:
    print("Result comparison error: {}".format(error), file=sys.stderr)
    return 1
  if failures:
    print("Reference comparison failed for {}:".format(args.case),
          file=sys.stderr)
    for failure in failures:
      print("  {}".format(failure), file=sys.stderr)
    return 1
  print("Reference comparison passed for {}".format(args.case))
  return 0


if __name__ == "__main__":
  sys.exit(main())
