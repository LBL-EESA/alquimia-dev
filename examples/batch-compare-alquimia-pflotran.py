#!/usr/bin/env python

"""
Program to compare the results of alquimia driven and pflotran native simulations.
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


def read_config_file(config_filename):
    """
    Read a configuration file of the form:

    [test-name]
    alquimia = alquimia-inputfile.cfg
    pflotran = pflotran-inputfile.in
    """
    if config_filename == None:
        raise Exception("Error, must provide a test config filename")
    print("Reading config file : {0}".format(config_filename))
    config = config_parser.SafeConfigParser()
    config.read(config_filename)

    tests = {}
    for name in config.sections():
        tests[name] = {}
        for item in config.items(name):
            tests[name][item[0]] = item[1]

    #pprinter.pprint(tests)
    return tests

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

def check_executable(executable):
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

def run_job(name, command, timeout):
    if False :
        print("    {0}".format(" ".join(command)), file=sys.stdout)

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
            message = self._txtwrap.fill(
                "ERROR: job '{0}' has exceeded timeout limit of "
                "{1} seconds.".format(name, timeout))
            print(''.join(['\n', message, '\n']), file=sys.stdout)
    run_stdout.close()
    return abs(proc.returncode)


def run_alquimia(alquimia, test_name, timeout):
    """
    Generate the run command, then run alquimia's batch_chem driver:

    batch_chem -d -i test_name.cfg

    """
    command = []

    command.append(alquimia)
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

    status = run_job(test_name, command, timeout)
    if status != 0:
        message = txtwrap.fill(
            "WARNING : {name} : alquimia driver return an error "
            "code ({status}) indicating the simulation may have "
            "failed. Please check '{name}.stdout' for error "
            "messages.".format(
                name=test_name, status=status))
        print("".join(['\n', message, '\n']), file=sys.stdout)
    return status


def run_pflotran(pflotran, test_name, timeout):
    PFLOTRAN_SUCCESS = 86

    command = []

    command.append(pflotran)
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

    pflotran_status = run_job(test_name, command, timeout)

    # pflotran returns 0 on an error (e.g. can't find an input
    # file), 86 on success. 59 for timeout errors?
    if pflotran_status != PFLOTRAN_SUCCESS:
        status = 1
        message = txtwrap.fill(
            "WARNING : {name} : pflotran return an error "
            "code ({status}) indicating the simulation may have "
            "failed. Please check '{name}.out' and '{name}.stdout' "
            "for error messages.".format(
                name=test_name, status=status))
        print("".join(['\n', message, '\n']), file=sys.stdout)
    else:
        status = 0
    return status

def run_test(alquimia, pflotran, test):
    status = -1
    print("Running test : {0}".format(test))
    alquimia_status = run_alquimia(alquimia, test["alquimia"], options.timeout)
    pflotran_status = run_pflotran(pflotran, test["pflotran"], options.timeout)
    if alquimia_status == 0 and pflotran_status == 0:
        status = compare_results(test["alquimia"],
                                 test["pflotran"])
        
    return status

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
    Simple class to parse and store the data from an alquimia batch_chem output file.
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

    
def compare_results(alquimia_name, pflotran_name):
    
    pflotran_data = PFloTranObservationData()
    pflotran_data.read(pflotran_name)
    print(pflotran_data)
    alquimia_data = AlquimiaObservationData()
    alquimia_data.read(alquimia_name)
    print(alquimia_data)

    status = 0
    if pflotran_data.time_units != alquimia_data.time_units:
        status = 1
        print("FAIL: time units differ : pflotran = {0}   alquimia = {1}".format(pflotran_data.time_units, alquimia_data.time_units))

    pshape = pflotran_data.columns.shape
    ashape = alquimia_data.columns.shape
    if (pshape[1] != pflotran_data.num_columns):
        print("ERROR: pflotran number of columns header does not equal number of data columns: header = {0}   data = {1}".format(pflotran_data.num_columns, pshape[1]))

    if (ashape[1] != alquimia_data.num_columns):
        print("ERROR: alquimia number of columns header does not equal number of data columns: header = {0}   data = {1}".format(alquimia_data.num_columns, ashape[1]))

    if pshape[1] != ashape[1]:
        print("FAIL: number of columns not equal : pflotran = {0}   alquimia = {1}".format(pshape[1], ashape[1]))

    if pshape[0] != ashape[0]:
        print("FAIL: number of rows not equal : pflotran = {0}   alquimia = {1}".format(pshape[0], ashape[0]))

    if not np.array_equal(pflotran_data.columns[:, 0], alquimia_data.columns[:, 0]):
        different = pflotran_data.columns[:,0] != alquimia_data.columns[:,0]
        for i in range(len(different)):
            if different[i]:
                print("FAIL: observation times are not equal:  row = "
                      "{0}   pflotran = {1}   alquimia = {2}".format(
                        i, pflotran_data.columns[i,0], alquimia_data.columns[i, 0]))

    #if not np.allclose(pflotran_data.columns[:,0], alquimia_data.columns[:,0]):
        

def main(options):
    #pprinter.pprint(options)
    tests = read_config_file(options.config_file[0])
    #pprinter.pprint(tests)

    alquimia = check_executable(options.alquimia_executable)
    pflotran = check_executable(options.pflotran_executable)

    print("  Using alquimia executable : {0}".format(alquimia))
    print("  Using pflotran executable : {0}".format(pflotran))
    
    report = {}
    if options.test is None:
        # run all tests in the config file
        for t in tests:
            report[t] = run_test(alquimia, pflotran, tests[t])
    else:
        # run just the specified test
        if options.test in tests:
            report[options.test] = run_test(alquimia, pflotran, tests[options.test])
        else:
            print("ERROR: test name specified on the command line is not "
                  "valid: '{0}'".format(options.test))
            print("  Valid tests are: ")
            pprinter.pprint(tests.keys())

    pprinter.pprint(report)

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
