import pyparsing as pp
import global_vars
from error_print import eprint
from functools import reduce

elements = {}

comment_filter = pp.Regex('[^#]*') + pp.Regex('#?.*')
number_of_loci = pp.Literal("number_of_loci") + pp.common.integer
number_of_alleles = pp.Literal('number_of_alleles') + pp.common.integer
rate_of_mutation = pp.Literal('rate_of_mutation') + pp.common.number
rate_of_recombination = pp.Literal('rate_of_recombination') + pp.common.number
# Not implemented
# population_size = pp.Literal('population_size') + pp.common.integer
number_of_generations = pp.Literal('number_of_generations') + pp.common.integer
threshold = pp.Literal('threshold') + pp.common.number

oneline_patterns = [number_of_loci,
            number_of_alleles,
            rate_of_mutation,
            rate_of_recombination,
            number_of_generations,
                    threshold]

oneline_any = reduce(lambda x,y: x | y, oneline_patterns)

def insert(toks):
    elements.update({toks[0]:toks[1]})

@comment_filter.add_parse_action
def drop(toks):
    return toks[0]

@number_of_loci.add_parse_action
def getnloci(toks):
    global_vars.nloci = int(toks[1])
    global_vars.ngtypes = global_vars.nalleles**global_vars.nloci

@number_of_alleles.add_parse_action
def getnalleles(toks):
    global_vars.nalleles = int(toks[1])

    # TODO Only permiting nallles <= 10 at this point
    if global_vars.nalleles > 10:
        eprint("No support for number of alleles exceeding 10")
        exit(1)

    global_vars.ngtypes = global_vars.nalleles**global_vars.nloci

@rate_of_mutation.add_parse_action
def assert_rate_clamp(toks):
    assert 0. <= toks[1] <= 1.

@rate_of_recombination.add_parse_action
def assert_rate_clamp(toks):
    assert 0. <= toks[1] <= 1.

@threshold.add_parse_action
def assert_threshold_clamp(toks):
    assert 0. <= toks[1] <= 1.

for p in oneline_patterns:
    if p is not comment_filter:
        p.add_parse_action(insert)
