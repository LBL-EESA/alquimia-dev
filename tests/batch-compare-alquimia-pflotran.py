#!/usr/bin/env python
"""
Program to compare the results of alquimia driven and pflotran native simulations.

Alquimia Copyright (c) 2013, The Regents of the University of California, 
through Lawrence Berkeley National Laboratory (subject to receipt of any 
required approvals from the U.S. Dept. of Energy).  All rights reserved.

Alquimia is available under a BSD license. See LICENSE.txt for more
information.

If you have questions about your rights to use or distribute this software, 
please contact Berkeley Lab's Technology Transfer and Intellectual Property 
Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).

NOTICE.  This software was developed under funding from the U.S. Department 
of Energy.  As such, the U.S. Government has been granted for itself and 
others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
license in the Software to reproduce, prepare derivative works, and perform 
publicly and display publicly.  Beginning five (5) years after the date 
permission to assert copyright is obtained from the U.S. Department of Energy, 
and subject to any subsequent five (5) year renewals, the U.S. Government is 
granted for itself and others acting on its behalf a paid-up, nonexclusive, 
irrevocable, worldwide license in the Software to reproduce, prepare derivative
works, distribute copies to the public, perform publicly and display publicly, 
and to permit others to do so.

Authors: Benjamin Andre <bandre@lbl.gov>

"""


from __future__ import print_function
from __future__ import division

import argparse
from collections import deque
import datetime
import math
import os
import pprint
import re
import subprocess
import sys
import textwrap
import time
import traceback

if sys.version_info[0] == 2:
    import ConfigParser as config_parser
else:
    import configparser as config_parser


import numpy as np

pprinter = pprint.PrettyPrinter(indent=2)
txtwrap = textwrap.TextWrapper(width=78, subsequent_indent=4*" ")


class PFloTranObservationData(object):
    """
    Simple class to parse and store a pflotran observation file
    """

    _time_re = re.compile("^Time[\s]+\[([\w]+)\]$")

    def __init__(self):
        self.filename = None
        self.num_columns = None
        self.time_units = None
        self.column_names = []
        self.columns = None

    def __str__(self):
        message = []
        message.append("PFloTran data :")
        message.append("    Filename : {0}".format(self.filename))
        message.append("    Number of columns : {0}".format(self.num_columns))
        message.append("    Time units : {0}".format(self.time_units))
        message.append("    Column names ({0}) : {1}".format(len(self.column_names),
                                                         self.column_names))
        message.append("    Number of rows : {0}".format(self.columns.shape[0]))

        return "\n".join(message)

    def read(self, filebase):
        self.filename = filebase + "-obs-0.tec"
        # manually process the header, then read the data columns with numpy
        with open(self.filename, 'r') as obs_file:
            header_line = obs_file.readline()
            self._process_header(header_line)

        self.columns = np.loadtxt(self.filename, comments='"')
        #print(self.columns)


    def _process_header(self, header_line):
        header = header_line.split(',')
        self.num_columns = len(header)

        self.column_names.append("Time")
        time_field = header[0].strip(' "')
        match = self._time_re.match(time_field)
        if match:
            self.time_units = match.group(1)
        else:
            # default seconds or error...?
            raise Exception("ERROR: Could not determine time units in file '{0}'."
                            "\n  Time field : '{1}'".format(self.filename,
                                                            time_field))
        for column_header in header[1:]:
            column = column_header.strip(' "')
            fields = column.split("all")
            self.column_names.append(fields[0].strip())
        if len(self.column_names) != self.num_columns:
            raise Exception("ERROR: Could not determine all column names in "
                            "observation file.\n  Expected {0}, found {1}.\n"
                            "  Identified : {2}\n  Header: {3}".format(
                    self.num_columns, len(self.column_names), self.column_names,
                    header_line))

class AlquimiaObservationData(object):
    """
    Simple class to parse and store the data from an alquimia
    batch_chem output file.
    """

    _time_re = re.compile("^Time[\s]+\[([\w]+)\]$")

    def __init__(self):
        self.filename = None
        self.num_columns = None
        self.time_units = None
        self.column_names = []
        self.columns = None

    def __str__(self):
        message = []
        message.append("Alquimia data :")
        message.append("    Filename : {0}".format(self.filename))
        message.append("    Number of columns : {0}".format(self.num_columns))
        message.append("    Time units : {0}".format(self.time_units))
        message.append("    Column names ({0}) : {1}".format(len(self.column_names),
                                                         self.column_names))
        message.append("    Number of rows : {0}".format(self.columns.shape[0]))

        return "\n".join(message)

    def read(self, filebase):
        self.filename = filebase + ".txt"
        # manually process the header, then read the data columns with numpy
        with open(self.filename, 'r') as obs_file:
            header_line = obs_file.readline()
            self._process_header(header_line)

        self.columns = np.loadtxt(self.filename, comments='#')
        #print(self.columns)


    def _process_header(self, header_line):
        header = header_line.strip("#")
        header = header.strip()
        header = header.split(',')
        self.num_columns = len(header)

        self.column_names.append("Time")
        time_field = header[0].strip(' "')
        match = self._time_re.match(time_field)
        if match:
            self.time_units = match.group(1)
        else:
            # default seconds or error...?
            raise Exception("ERROR: Could not determine time units in file '{0}'."
                            "\n  Time field : '{1}'".format(self.filename,
                                                            time_field))
        for column_header in header[1:]:
            column = column_header.strip(' "')
            self.column_names.append(column)
        if len(self.column_names) != self.num_columns:
            raise Exception("ERROR: Could not determine all column names in "
                            "observation file.\n  Expected {0}, found {1}.\n"
                            "  Identified : {2}\n  Header: {3}".format(
                    self.num_columns, len(self.column_names), self.column_names,
                    header_line))


class TestManager(object):

    def __init__(self):
        self._config_filename = None
        self._alquimia = None
        self._pflotran = None
        self._test_info = None
        self._log_filename = None
        self._test_log = None
        self.report = {}

        self._setup_logfile()

    def __str__(self):
        message = ["Alquimia Batch Chemistry Test Manager:"]
        message.append("    log file : {0}".format(self._log_filename))
        message.append("    config file : {0}".format(self._config_filename))
        message.append("    alquimia : {0}".format(self._alquimia))
        message.append("    pflotran : {0}".format(self._pflotran))

        return "\n".join(message)

    def set_alquimia(self, alquimia):
        self._alquimia = self._check_executable(alquimia)


    def set_pflotran(self, pflotran):
        self._pflotran = self._check_executable(pflotran)


    def get_tests_from_config_file(self, config_filename):
        """
        Read a configuration file of the form:

        [test-name]
        alquimia = alquimia-inputfile-prefix
        pflotran = pflotran-inputfile-prefix


        NOTES:

            * the input filenames are prefix only. The suffix .cfg,
              .in, etc will be added automatically.
        """
        if config_filename == None:
            raise Exception("Error, must provide a test config filename")
        self._config_filename = config_filename
        #print("Reading config file : {0}".format(self._config_filename))
        config = config_parser.SafeConfigParser()
        config.read(self._config_filename)

        self._test_info = {}
        for name in config.sections():
            self._test_info[name] = {}
            for item in config.items(name):
                self._test_info[name][item[0]] = item[1]

        #pprinter.pprint(self._test_info)


    def run_all_tests(self, timeout):
        print("Running tests :")
        for t in self._test_info:
            status = self.run_test(t, timeout)
            if status == 0:
                print(".", end='', file=sys.stdout)
            else:
                print("F", end='', file=sys.stdout)
            sys.stdout.flush()

        print("", file=sys.stdout)


    def run_test(self, test_name, timeout):
        status = -1

        if test_name in self._test_info:
            print(70*'_', file=self._test_log)
            print("Running test : {0}".format(test_name), file=self._test_log)
            test = self._test_info[test_name]
            alquimia_status = self._run_alquimia(test["alquimia"], timeout)
            pflotran_status = self._run_pflotran(test["pflotran"], timeout)
            if alquimia_status == 0 and pflotran_status == 0:
                status = compare_results(self._test_log,
                                         test["alquimia"],
                                         test["pflotran"])
        else:
            print("ERROR: test name is not "
                  "valid: '{0}'".format(options.test))
            print("  Valid tests are: ")
            pprinter.pprint(tests.keys())
            status = 1

        self.report[test_name] = status

        return status

    def _setup_logfile(self):
        self._log_filename = "alquimia-tests-{0}.testlog".format(datetime.datetime.today().strftime("%Y-%m-%d_%H-%M-%S"))
        self._test_log = open(self._log_filename, 'w')
        print("Alquimia / PFloTran Batch Chemistry Test Log", file=self._test_log)
        print("System Info :", file=self._test_log)
        print("    platform : {0}".format(sys.platform), file=self._test_log)


    def _run_job(self, name, command, timeout):
        print("    {0}".format(" ".join(command)), file=self._test_log)

        run_stdout = open(name + ".stdout", 'w')
        start = time.time()
        proc = subprocess.Popen(command,
                                shell=False,
                                stdout=run_stdout,
                                stderr=subprocess.STDOUT)
        while proc.poll() is None:
            time.sleep(0.1)
            if time.time() - start > timeout:
                proc.terminate()
                time.sleep(0.1)
                message = txtwrap.fill(
                    "ERROR: job '{0}' has exceeded timeout limit of "
                    "{1} seconds.".format(name, timeout))
                print(''.join(['\n', message, '\n']), file=self._test_log)
        run_stdout.close()
        return abs(proc.returncode)


    def _run_alquimia(self, test_name, timeout):
        """
        Generate the run command, then run alquimia's batch_chem driver:

        batch_chem -d -i test_name.cfg

        """
        command = []

        command.append(self._alquimia)
        command.append("-d")
        command.append("-i")
        input_file_name = test_name + ".cfg"
        command.append(input_file_name)

        if os.path.isfile(test_name + ".txt"):
            os.rename(test_name + ".txt",
                      test_name + ".txt.old")

        if os.path.isfile(test_name + ".stdout"):
            os.rename(test_name + ".stdout",
                      test_name + ".stdout.old")

        status = self._run_job(test_name, command, timeout)
        if status != 0:
            message = txtwrap.fill(
                "ERROR : {name} : alquimia driver return an error "
                "code ({status}) indicating the simulation may have "
                "failed. Please check '{name}.stdout' for error "
                "messages.".format(
                    name=test_name, status=status))
            print("".join(['\n', message, '\n']), file=self._test_log)
        return status


    def _run_pflotran(self, test_name, timeout):
        PFLOTRAN_SUCCESS = 86

        command = []

        command.append(self._pflotran)
        command.append("-input_prefix")
        command.append(test_name)

        if os.path.isfile(test_name + "-obs-0.tec"):
            os.rename(test_name + "-obs-0.tec",
                      test_name + "-obs-0.tec.old")

        if os.path.isfile(test_name + ".out"):
            os.rename(test_name + ".out",
                      test_name + ".out.old")

        if os.path.isfile(test_name + ".stdout"):
            os.rename(test_name + ".stdout",
                      test_name + ".stdout.old")

        pflotran_status = self._run_job(test_name, command, timeout)

        # pflotran returns 0 on an error (e.g. can't find an input
        # file), 86 on success. 59 for timeout errors?
        if pflotran_status != PFLOTRAN_SUCCESS:
            status = 1
            message = txtwrap.fill(
                "ERROR : {name} : pflotran return an error "
                "code ({status}) indicating the simulation may have "
                "failed. Please check '{name}.out' and '{name}.stdout' "
                "for error messages.".format(
                    name=test_name, status=status))
            print("".join(['\n', message, '\n']), file=self._test_log)
        else:
            status = 0
        return status

    def _check_executable(self, executable):
        """
        Try to verify that we have something reasonable for the executable
        """
        # absolute path to the executable
        executable = os.path.abspath(executable)
        # is it a valid file?
        if not os.path.isfile(executable):
            raise Exception("ERROR: executable is not a valid file: "
                            "'{0}'".format(executable))
        return executable



def compare_columns(test_log, name,
                    alquimia_data, aindex,
                    pflotran_data, pindex):
    status = 0
    if not np.array_equal(pflotran_data.columns[:, pindex],
                          alquimia_data.columns[:, aindex]):
        status = 1
        different = (pflotran_data.columns[:, pindex] !=  
                     alquimia_data.columns[:, aindex])
        for i in range(len(different)):
            if different[i]:
                print("FAIL: observations are not equal for column '{0}':  row = "
                      "{1}   pflotran = {2}   alquimia = {3}".format(name,
                        i, pflotran_data.columns[i, pindex],
                        alquimia_data.columns[i, aindex]), file=test_log)

    return status
    
def compare_results(test_log, alquimia_name, pflotran_name):
    
    pflotran_data = PFloTranObservationData()
    pflotran_data.read(pflotran_name)
    print(pflotran_data, file=test_log)
    alquimia_data = AlquimiaObservationData()
    alquimia_data.read(alquimia_name)
    print(alquimia_data, file=test_log)

    status = 0
    if pflotran_data.time_units != alquimia_data.time_units:
        status += 1
        print("FAIL: time units differ : "
              "pflotran = {0}   alquimia = {1}".format(
                pflotran_data.time_units, alquimia_data.time_units), file=test_log)

    pshape = pflotran_data.columns.shape
    ashape = alquimia_data.columns.shape
    if (pshape[1] != pflotran_data.num_columns):
        status += 1
        print("ERROR: pflotran number of columns header does not equal number "
              "of data columns: header = {0}   data = {1}".format(
                pflotran_data.num_columns, pshape[1]), file=test_log)

    if (ashape[1] != alquimia_data.num_columns):
        status += 1
        print("ERROR: alquimia number of columns header does not equal number "
              "of data columns: header = {0}   data = {1}".format(
                alquimia_data.num_columns, ashape[1]), file=test_log)

    if pshape[1] != ashape[1]:
        status += 1
        print("FAIL: number of columns not equal : "
              "pflotran = {0}   alquimia = {1}".format(pshape[1], ashape[1]), file=test_log)

    if pshape[0] != ashape[0]:
        status += 1
        print("FAIL: number of rows not equal : "
              "pflotran = {0}   alquimia = {1}".format(pshape[0], ashape[0]), file=test_log)

    for (alquimia_index, name) in enumerate(alquimia_data.column_names):
        try:
            pflotran_index = pflotran_data.column_names.index(name)
            status += compare_columns(test_log, name,
                                      alquimia_data, alquimia_index,
                                      pflotran_data, pflotran_index)
        except ValueError as e:
            status += 1
            print("ERROR: Alquimia column '{0}' is not in pflotran data.".format(name), file=test_log)
        
    return status

def commandline_options():
    parser = argparse.ArgumentParser(
        description="Run alquimia driven and pflotran native batch chemistry "
        "simulations and compare the results.")

    parser.add_argument("--backtrace", action='store_true',
                        help="show exception backtraces as extra debugging "
                        "output")

    parser.add_argument('-a', "--alquimia-executable", required=True,
                        help="path to the alquimia batch_chem executable")

    parser.add_argument('-c', "--config-file", nargs=1, default=None,
                        required=True,
                        help="test configuration file to use")

    parser.add_argument('-p', "--pflotran-executable", required=True,
                        help="path to the pflotran executable")

    parser.add_argument('-t', "--test", default=None,
                        help="name of the test to run")
    parser.add_argument("--timeout", default=30,
                        help="timeout limit for simulations [seconds]")

    options = parser.parse_args()
    return options


def summary_report(run_time, report):
    status = 0
    print(70 * '-')
    print("Regression test summary :")
    print("    Total run time: {0:4g} [s]".format(run_time))
    for t in report:
        status += report[t]
        if report[t] > 0:
            print("    {0}... failed".format(t))
        elif report[t] == 0:
            print("    {0}... passed".format(t))
    print("\n\n")
    return status

def main(options):
    #pprinter.pprint(options)
    test_manager = TestManager()
    test_manager.get_tests_from_config_file(options.config_file[0])
    #pprinter.pprint(tests)

    test_manager.set_alquimia(options.alquimia_executable)
    test_manager.set_pflotran(options.pflotran_executable)

    if True:
        print(test_manager)

    start = time.time()
    if options.test is None:
        # run all tests in the config file
        test_manager.run_all_tests(options.timeout)
    else:
        test_manager.run_test(options.test, options.timeout)

    stop = time.time()

    summary_report(stop - start, test_manager.report)



if __name__ == "__main__":
    options = commandline_options()
    try:
        status = main(options)
        sys.exit(status)
    except Exception as e:
        print(str(e))
        if options.backtrace:
            traceback.print_exc()
        sys.exit(1)
