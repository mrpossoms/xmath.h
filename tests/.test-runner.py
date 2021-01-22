#!/usr/bin/env python

import os
import sys
import subprocess

COLOR_GREEN='\033[0;32m'
COLOR_RED='\033[0;31m'
COLOR_OFF='\033[0m'

def getTestList():
    tests = []
    for file_name in os.listdir("./bin"):
        if '.dSYM' in file_name:
            continue
        tests.append(os.path.join("./bin", file_name))
    return tests

testNames = getTestList()
testsPassed = 0
testsRan = 0

def runTest(file_path):
    args = [file_path, str(testsRan), str(len(testNames))]
    return subprocess.call(args)

print('Running xmath Tests')

for file_path in testNames:
    sys.stdout.write('{}...'.format(file_path))
    errored = runTest(file_path)

    if not errored:
        sys.stdout.write('{}OK{}\n'.format(COLOR_GREEN, COLOR_OFF))
        testsPassed += 1
    else:
        sys.stdout.write('{}FAIL{}\n'.format(COLOR_RED, COLOR_OFF))

    testsRan += 1

status = COLOR_GREEN

# check o see if any failed
if testsPassed != len(testNames):
    status = COLOR_RED

print(str(testsRan) + " of " + str(len(testNames)) + " ran")
print(status + str(testsPassed) + " of " + str(len(testNames)) + " passed\033[0m")

# some tests failed, do not return success
if status == COLOR_RED:
    exit(1)
