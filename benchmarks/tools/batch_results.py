#!/usr/bin/env python3
"""Readers and numerical comparisons for Alquimia batch results."""

import csv
import math
import re


class ResultsError(Exception):
  """Raised when a results file is malformed or incompatible."""


class Table(object):
  def __init__(self, columns, rows):
    self.columns = columns
    self.rows = rows


def _validate_table(path, columns, rows):
  if not columns:
    raise ResultsError("{} has no header".format(path))
  duplicates = sorted(set(name for name in columns if columns.count(name) > 1))
  if duplicates:
    raise ResultsError("{} has duplicate columns: {}".format(
        path, ", ".join(duplicates)))
  for row_number, row in enumerate(rows, 2):
    if len(row) != len(columns):
      raise ResultsError("{}:{} has {} values; expected {}".format(
          path, row_number, len(row), len(columns)))
    for column, value in zip(columns, row):
      if not math.isfinite(value):
        raise ResultsError("{}:{} column '{}' is not finite".format(
            path, row_number, column))


def read_csv_table(path):
  with open(path, "r", newline="") as stream:
    reader = csv.reader(stream)
    try:
      columns = [name.strip() for name in next(reader)]
    except StopIteration:
      raise ResultsError("{} is empty".format(path))
    rows = []
    for row_number, fields in enumerate(reader, 2):
      if not fields or all(not field.strip() for field in fields):
        continue
      try:
        rows.append([float(field) for field in fields])
      except ValueError as error:
        raise ResultsError("{}:{} contains a non-numeric value: {}".format(
            path, row_number, error))
  _validate_table(path, columns, rows)
  return Table(columns, rows)


def write_csv_table(path, table):
  with open(path, "w", newline="") as stream:
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(table.columns)
    for row in table.rows:
      writer.writerow(["{:.17g}".format(value) for value in row])


def compare_tables(actual, reference, rel_tol=1.0e-6, abs_tol=1.0e-12,
                   time_abs_tol=1.0e-9, overrides=None, max_failures=20):
  overrides = overrides or {}
  failures = []
  actual_indices = {name: index for index, name in enumerate(actual.columns)}
  missing = [name for name in reference.columns if name not in actual_indices]
  if missing:
    failures.append("actual output is missing reference columns: {}".format(
        ", ".join(missing)))
  if failures:
    return failures

  if not reference.columns or reference.columns[0] != "time":
    return ["reference table must use 'time' as its first column"]
  actual_time_index = actual_indices["time"]
  row_pairs = []
  for reference_row_index, reference_row in enumerate(reference.rows):
    expected_time = reference_row[0]
    matches = [index for index, row in enumerate(actual.rows)
               if math.isclose(row[actual_time_index], expected_time,
                               rel_tol=0.0, abs_tol=time_abs_tol)]
    if len(matches) != 1:
      failures.append(
          "reference time {} has {} matching actual rows; expected 1".format(
              expected_time, len(matches)))
    else:
      row_pairs.append((matches[0], reference_row_index))
  if failures:
    return failures

  for reference_index, name in enumerate(reference.columns):
    actual_index = actual_indices[name]
    tolerance = overrides.get(name, {})
    column_rel_tol = float(tolerance.get("rel_tol", rel_tol))
    column_abs_tol = float(tolerance.get("abs_tol", abs_tol))
    if name == "time":
      column_rel_tol = 0.0
      column_abs_tol = time_abs_tol
    worst = None
    for actual_row_index, reference_row_index in row_pairs:
      expected = reference.rows[reference_row_index][reference_index]
      observed = actual.rows[actual_row_index][actual_index]
      if math.isclose(observed, expected, rel_tol=column_rel_tol,
                      abs_tol=column_abs_tol):
        continue
      absolute_error = abs(observed-expected)
      scale = max(abs(observed), abs(expected))
      relative_error = absolute_error/scale if scale else 0.0
      candidate = (absolute_error, reference_row_index, observed, expected,
                   relative_error)
      if worst is None or candidate > worst:
        worst = candidate
    if worst is not None:
      absolute_error, row_index, observed, expected, relative_error = worst
      time_value = reference.rows[row_index][0]
      failures.append(
          "{} at row {} (time={}): actual={:.17g}, reference={:.17g}, "
          "abs_error={:.6g}, rel_error={:.6g}".format(
              name, row_index+1, time_value, observed, expected,
              absolute_error, relative_error))
      if len(failures) >= max_failures:
        failures.append("additional failures omitted")
        break
  return failures


def _fortran_float(value):
  value = re.sub(r"(?<=\d)([+-]\d{3})$", r"e\1", value)
  return float(value.replace("D", "e").replace("d", "e"))


def read_pflotran_observation(path):
  with open(path, "r", newline="") as stream:
    header = stream.readline()
    if not header:
      raise ResultsError("{} is empty".format(path))
    columns = [name.strip().strip('"')
               for name in next(csv.reader([header]))]
    rows = []
    for line in stream:
      fields = line.split()
      if fields:
        rows.append([_fortran_float(field) for field in fields])
  _validate_table(path, columns, rows)
  return Table(columns, rows)


def read_crunchflow_tecplot(path):
  columns = None
  rows = []
  with open(path, "r") as stream:
    for line_number, line in enumerate(stream, 1):
      stripped = line.strip()
      if not stripped:
        continue
      if stripped.upper().startswith("VARIABLES"):
        header = stripped.split("=", 1)[1]
        columns = [name.strip().strip('"').strip()
                   for name in next(csv.reader([header]))]
        columns = [name for name in columns if name]
        continue
      if stripped.upper().startswith(("TITLE", "ZONE")):
        continue
      if columns is None:
        continue
      fields = stripped.split()
      try:
        rows.append([_fortran_float(field) for field in fields])
      except ValueError as error:
        raise ResultsError("{}:{} contains invalid Tecplot data: {}".format(
            path, line_number, error))
  if columns is None:
    raise ResultsError("{} has no VARIABLES header".format(path))
  _validate_table(path, columns, rows)
  return Table(columns, rows)


def normalize_table(table, mappings, time_scale=1.0):
  source_indices = {name: index for index, name in enumerate(table.columns)}
  columns = []
  selections = []
  for mapping in mappings:
    source = mapping.get("source")
    if source is None:
      prefix = mapping["source_prefix"]
      matches = [name for name in table.columns if name.startswith(prefix)]
      if len(matches) != 1:
        raise ResultsError(
            "native output has {} columns beginning with '{}'; expected 1".
            format(len(matches), prefix))
      source = matches[0]
    if source not in source_indices:
      raise ResultsError("native output is missing column '{}'".format(source))
    columns.append(mapping["target"])
    selections.append((source_indices[source], float(mapping.get("scale", 1.0))))
  rows = []
  for source_row in table.rows:
    row = [source_row[index]*scale for index, scale in selections]
    if columns and columns[0] == "time":
      row[0] *= time_scale
    rows.append(row)
  _validate_table("normalized table", columns, rows)
  return Table(columns, rows)
