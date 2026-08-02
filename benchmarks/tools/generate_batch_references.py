#!/usr/bin/env python3
"""Generate normalized batch references from standalone chemistry engines."""

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

from batch_results import (ResultsError, normalize_table,
                           read_crunchflow_tecplot,
                           read_pflotran_observation, write_csv_table)


INPUT_SUFFIXES = (".in", ".dat", ".dbs", ".dbsx")


def _sha256(path):
  digest = hashlib.sha256()
  with open(path, "rb") as stream:
    for block in iter(lambda: stream.read(65536), b""):
      digest.update(block)
  return digest.hexdigest()


def _wrap_legacy_pflotran_input(path):
  with open(path, "r") as stream:
    contents = stream.read()
  contents = "\n".join(
      "#"+line[1:] if line.lstrip().startswith(":") else line
      for line in contents.splitlines())+"\n"
  contents = re.sub(
      r"(?ms)^TIMESTEPPER(?:\s+TRANSPORT)?\s*\n(.*?)^/\s*$",
      r"NUMERICAL_METHODS TRANSPORT\n  TIMESTEPPER\n\1  /\nEND",
      contents, count=1)
  bounds = re.search(
      r"(?ms)^(\s*)BOUNDS\s*\n\s*(\S+)\s+(\S+)\s*\n"
      r"\s*(\S+)\s+(\S+)\s*\n\s*(\S+)\s+(\S+)\s*\n\s*/",
      contents)
  if bounds:
    indent = bounds.group(1)
    replacement = ("{}BOUNDS\n{}  {} {} {}\n{}  {} {} {}\n{}/".format(
        indent, indent, bounds.group(2), bounds.group(4), bounds.group(6),
        indent, bounds.group(3), bounds.group(5), bounds.group(7), indent))
    contents = contents[:bounds.start()]+replacement+contents[bounds.end():]
  if "ISOTHERM_REACTIONS" in contents and "ROCK_DENSITY" not in contents:
    contents = re.sub(r"(?m)^(\s*ID\s+\d+\s*)$",
                      r"\1\n  ROCK_DENSITY 2500.0", contents, count=1)
  wrapper = ("SIMULATION\n"
             "  SIMULATION_TYPE SUBSURFACE\n"
             "  PROCESS_MODELS\n"
             "    SUBSURFACE_TRANSPORT transport\n"
             "      MODE GIRT\n"
             "    /\n"
             "  /\n"
             "END\n\nSUBSURFACE\n\n")
  with open(path, "w") as stream:
    stream.write(wrapper)
    stream.write(contents)
    stream.write("\nEND_SUBSURFACE\n")


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--manifest", required=True)
  parser.add_argument("--case", action="append", dest="cases")
  parser.add_argument("--pflotran")
  parser.add_argument("--crunchflow")
  parser.add_argument("--update", action="store_true",
                      help="allow committed reference files to be replaced")
  args = parser.parse_args()
  if not args.update:
    parser.error("--update is required to write reference files")

  manifest_path = os.path.abspath(args.manifest)
  reference_dir = os.path.dirname(manifest_path)
  benchmark_dir = os.path.dirname(reference_dir)
  with open(manifest_path, "r") as stream:
    manifest = json.load(stream)
  selected = args.cases or sorted(manifest["cases"])

  for case_name in selected:
    case = manifest["cases"][case_name]
    native = case["native"]
    engine = native["engine"]
    executable = args.pflotran if engine == "pflotran" else args.crunchflow
    if not executable:
      parser.error("--{} is required for {}".format(engine, case_name))
    with tempfile.TemporaryDirectory(prefix="alquimia-reference-") as work_dir:
      for filename in os.listdir(benchmark_dir):
        source = os.path.join(benchmark_dir, filename)
        if os.path.isfile(source) and filename.endswith(INPUT_SUFFIXES):
          shutil.copy2(source, os.path.join(work_dir, filename))
      input_name = native["input"]
      if engine == "pflotran":
        if native.get("legacy_wrapper", False):
          _wrap_legacy_pflotran_input(os.path.join(work_dir, input_name))
        prefix = os.path.splitext(input_name)[0]
        command = [os.path.abspath(executable), "-input_prefix", prefix]
      else:
        command = [os.path.abspath(executable), input_name]
      log_path = os.path.join(work_dir, "native.log")
      with open(log_path, "w") as log:
        result = subprocess.run(command, cwd=work_dir, stdout=log,
                                stderr=subprocess.STDOUT)
      if result.returncode != 0:
        raise ResultsError("{} failed with status {} (log: {})".format(
            case_name, result.returncode, log_path))
      output_path = os.path.join(work_dir, native["output"])
      if engine == "pflotran":
        table = read_pflotran_observation(output_path)
      else:
        table = read_crunchflow_tecplot(output_path)
      normalized = normalize_table(
          table, native["columns"], float(native.get("time_scale", 1.0)))
      destination = os.path.join(reference_dir, case["reference"])
      os.makedirs(os.path.dirname(destination), exist_ok=True)
      write_csv_table(destination, normalized)
      case["provenance"] = {
          "executable": os.path.basename(executable),
          "input_sha256": _sha256(os.path.join(benchmark_dir, input_name)),
          "command": [os.path.basename(command[0])]+command[1:],
      }
      print("updated {}".format(destination))

  with open(manifest_path, "w") as stream:
    json.dump(manifest, stream, indent=2, sort_keys=True)
    stream.write("\n")
  return 0


if __name__ == "__main__":
  try:
    sys.exit(main())
  except (OSError, ResultsError, ValueError, KeyError) as error:
    print("Reference generation failed: {}".format(error), file=sys.stderr)
    sys.exit(1)
