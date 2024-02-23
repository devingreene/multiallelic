#!/usr/bin/env python3

import pyparsing as pp
import sys
import os
import struct
import ctypes
import global_vars
import oneline
import multiline
from oneline import elements
from error_print import eprint

assert global_vars.nloci is not None and global_vars.nalleles is not None

all_patterns = oneline.oneline_any | multiline.multiline_any

file_name_num = sys.stdin.fileno()
if os.isatty(0):
    file_name_num = 'input.txt'

# TODO Parsing error should happen if whole line (without comments) is
# not parsed.
with open(file_name_num, "r") as file:
    for line in file:
        line = line.strip()

        check = 0
        for ml in multiline.multiline_objects:
            check += ml.state
        assert check <= 1

        for ml in multiline.multiline_objects:
            if ml.state:
                break

        if ml.state:
            if line == '':
                ml.exit()
            else:
                # Mid-block full comments are allowed
                line = oneline.comment_filter.parse_string(line).asList()[0]
                ml.lines.append(line)
            continue

        # Drop comments from line
        # This uses regexps, so I don't think it can fail.
        line = oneline.comment_filter.parse_string(line).asList()[0]
        if line == '':
            continue

        try:
            all_patterns.parse_string(line)
        except pp.ParseException as e:
            eprint(e)
            eprint("Failed to parse line:\n---> {}".format(line))
            exit(1)
        except pp.ParseException as e:
            eprint(e)
            exit(1)

# Reasonable parameter checking
global_vars.ngtypes = global_vars.nalleles**global_vars.nloci
x = global_vars.ngtypes
if  x*x*(x + 1) * 4 > 0x8000000:
    eprint("Size of tables would exceed memory bound")
    exit(1)

cstruct = struct.pack("IIddId",
                      elements.get('number_of_loci', 0),
                      elements.get('number_of_alleles', 0),
                      elements.get('rate_of_mutation', 1e-9),
                      elements.get('rate_of_recombination', 0),
                      # elements.get('population_size', 10000),
                      elements.get('number_of_generations', 100),
                      elements.get('threshold', 0.99)
                      )
sys.stdout.buffer.write(cstruct)

for mll,mlt,*_ in multiline.multiline_labels_types:
    if mll not in multiline.multiline_labels_found:
        eprint("Parse error: Couldn't find line with: {}".
               format(mll))
        exit(1)
    mlo = multiline.Multiline.register[mll]
    assert mlo.type == mlt

    if mlo.type == 'dict':
        Type = 'double'
    elif mlo.type == 'list':
        Type = 'uint'
    else: assert False

    sys.stdout.buffer.write(
            (getattr(ctypes,'c_' + Type) * global_vars.ngtypes)(*mlo.data)
    )
