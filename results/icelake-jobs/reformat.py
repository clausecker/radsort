#!/usr/bin/env python3

"""Reformat radsort output files to a single data file for analysis."""

# usage: cat *.out | python3 reformat.py > results-all-reformatted.output
# or, if reformat.py is executable
# usage: cat *.out | ./reformat.py > results-all-reformatted.output

import sys

headers = {}
for line in sys.stdin:
    line = line.strip()
    if '' == line:
        continue
    if not 'Benchmark' in line:
        key, value = line.split(':')
        headers[key] = value.strip()
    else:
        parts = line.split()
        name = parts[0]
        count = parts[1]
        time = parts[2]
        timeunit = parts[3]
        throughput = parts[4]
        throughputunit = parts[5]
        new_line = name + "-" + headers["threads"]
        new_line += " " + headers["threads"]
        new_line += " " + count
        new_line += " " + time
        new_line += " " + timeunit
        new_line += " " + throughput
        new_line += " " + throughputunit
        new_line += " " + "arch: " + headers["arch"]
        new_line += " " + "os: " + headers["os"]
        new_line += " " + "items: " + headers["items"]
        print(new_line)
