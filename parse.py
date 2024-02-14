#!/usr/bin/env python3

import pyparsing as pp
import sys
import struct
import ctypes
import global_vars
import oneline
import multiline
from oneline import elements
from error_print import eprint

eprint("Python program:")

assert global_vars.nloci is not None and global_vars.nalleles is not None

all_patterns = oneline.oneline_any | multiline.multiline_any

with open("input.txt", "r") as file:
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

global_vars.ngtypes = global_vars.nalleles**global_vars.nloci
x = global_vars.ngtypes
if  x*x*(x + 1) * 4 > 0x8000000:
    eprint("Size of tables would exceed memory bound")
    exit(1)

cstruct = struct.pack("IIddIxxxx",
                      elements.get('number_of_loci', 0),
                      elements.get('number_of_alleles', 0),
                      elements.get('rate_of_mutation', 1e-9),
                      elements.get('rate_of_recombination', 0),
                      # elements.get('population_size', 10000),
                      elements.get('number_of_generations', 100),
                      )
sys.stdout.buffer.write(cstruct)

for mll in multiline.multiline_labels:
    mlo = multiline.Multiline.register[mll]
    eprint("{}: {}".format(mlo.label, mlo.data))
    sys.stdout.buffer.write(
            (ctypes.c_double * global_vars.ngtypes)(*mlo.data)
    )
