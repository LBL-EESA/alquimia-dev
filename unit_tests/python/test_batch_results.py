#!/usr/bin/env python3

import os
import tempfile
import unittest

from batch_results import (ResultsError, Table, compare_tables,
                           normalize_table, read_crunchflow_tecplot,
                           read_csv_table, read_pflotran_observation)


FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")


class BatchResultsTest(unittest.TestCase):
  def test_csv_and_reordered_columns(self):
    reference = read_csv_table(os.path.join(FIXTURES, "expected.csv"))
    actual = Table(["value", "extra", "time"],
                   [[1.0, 99.0, 0.0], [2.000001, 99.0, 1.0]])
    self.assertEqual([], compare_tables(actual, reference))

  def test_reports_missing_column(self):
    reference = Table(["time", "value"], [[0.0, 1.0]])
    actual = Table(["time"], [[0.0]])
    failures = compare_tables(actual, reference)
    self.assertIn("missing reference columns", failures[0])

  def test_reports_row_and_time_mismatch(self):
    reference = Table(["time", "value"], [[0.0, 1.0], [1.0, 2.0]])
    actual = Table(["time", "value"], [[0.1, 1.0]])
    failures = compare_tables(actual, reference)
    self.assertIn("matching actual rows", failures[0])

  def test_rejects_nonfinite_and_duplicate_columns(self):
    with self.assertRaises(ResultsError):
      read_csv_table(os.path.join(FIXTURES, "malformed.csv"))
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv") as stream:
      stream.write("time,time\n0,0\n")
      stream.flush()
      with self.assertRaises(ResultsError):
        read_csv_table(stream.name)

  def test_native_ascii_readers_and_normalization(self):
    pflotran = read_pflotran_observation(
        os.path.join(FIXTURES, "pflotran-observation.tec"))
    crunchflow = read_crunchflow_tecplot(
        os.path.join(FIXTURES, "crunchflow-timeseries.tec"))
    self.assertEqual(2, len(pflotran.rows))
    self.assertEqual(2, len(crunchflow.rows))
    normalized = normalize_table(
        crunchflow,
        [{"source": "Time (day)", "target": "time"},
         {"source": "A", "target": "total_mobile[A]", "scale": 2.0}],
        time_scale=86400.0)
    self.assertEqual([86400.0, 4.0], normalized.rows[1])


if __name__ == "__main__":
  unittest.main()
