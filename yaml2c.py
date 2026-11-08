#!/usr/bin/env python3
"""
yaml2c.py -- fill the `// yaml`-tagged declarations in evolve.c from input.yaml.

Usage:
    python3 yaml2c.py input.yaml evolve.c -o evolve_gen.c

The template (evolve.c) is never modified; a filled copy is written to -o
(or to stdout).  Genotype keys ("00", "12", ...) are read as base
`number_of_alleles` integers, most significant digit first, and used as array
indices.  Missing keys default to 0 in initial_state and 1 in fitness.
Fields of struct params are listed in SCALAR_PARAMS; stand-alone tagged
variables in main() go in STANDALONE.

The YAML is read with every scalar as a string, so unquoted keys such as 01,
08 or 10 are not mangled by YAML's int/octal rules.
"""

import argparse
import math
import re
import sys

import yaml

SENTINEL = "(uint)-1"

def flag(text):
    """0/1, or YAML-style true/false, yes/no, on/off -> 0 or 1."""
    t = str(text).strip().lower()
    if t in ("1", "true", "yes", "on"):
        return 1
    if t in ("0", "false", "no", "off"):
        return 0
    raise ValueError(text)


SCALAR_PARAMS = {             # name -> (C type, Python converter)
    "number_of_loci":        ("uint",   int),
    "number_of_alleles":     ("uint",   int),
    "rate_of_mutation":      ("double", float),
    "rate_of_recombination": ("double", float),
    "number_of_generations": ("uint",   int),
    "threshold":             ("double", float),
    "allow_double_mutations": ("uint",  flag),
}
SECTIONS = {"initial_state", "fitness", "target_genotypes", "report_points"}


# Stand-alone variables in main(), each filled from a top-level YAML key of the
# same name:  `<ctype> <name>; // yaml`  ->  `<ctype> <name> = <value>;`
# To add another one, add a line here and tag the C declaration.
STANDALONE = {                # name -> (C type, Python converter)
}


class ConfigError(Exception):
    pass


def warn(msg):
    print(f"yaml2c: warning: {msg}", file=sys.stderr)


# ---------------------------------------------------------------- YAML ----

class StrictStringLoader(yaml.BaseLoader):
    """All scalars as strings; duplicate mapping keys are an error."""


def _mapping_no_dupes(loader, node):
    seen = set()
    for key_node, _ in node.value:
        key = loader.construct_object(key_node)
        if key in seen:
            raise ConfigError(f"duplicate key {key!r} (line {key_node.start_mark.line + 1})")
        seen.add(key)
    return loader.construct_mapping(node)


StrictStringLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping_no_dupes)


def to_number(text, conv, what):
    try:
        val = conv(text)
    except (TypeError, ValueError):
        raise ConfigError(f"{what}: {text!r} is not a valid {conv.__name__}")
    if conv is int and val < 0:
        raise ConfigError(f"{what}: {text!r} must be >= 0")
    if conv is float and not math.isfinite(val):
        raise ConfigError(f"{what}: {text!r} must be finite")
    return val


def genotype_index(key, nloci, nalleles, what):
    key = str(key)
    if len(key) != nloci or not key.isdigit() or any(int(d) >= nalleles for d in key):
        raise ConfigError(
            f"{what}: genotype {key!r} must be {nloci} digits, each 0..{nalleles - 1}")
    return int(key, nalleles)


def load_config(path):
    with open(path) as fh:
        raw = yaml.load(fh, Loader=StrictStringLoader)
    if not isinstance(raw, dict):
        raise ConfigError("top level of the YAML file must be a mapping")

    unknown = set(raw) - set(SCALAR_PARAMS) - SECTIONS - set(STANDALONE)
    if unknown:
        raise ConfigError(f"unknown key(s): {', '.join(sorted(unknown))}")
    missing = (set(SCALAR_PARAMS) | set(STANDALONE)
               | {"initial_state", "target_genotypes"}) - set(raw)
    if missing:
        raise ConfigError(f"missing key(s): {', '.join(sorted(missing))}")

    cfg = {name: to_number(raw[name], conv, name)
           for name, (_, conv) in {**SCALAR_PARAMS, **STANDALONE}.items()}

    nloci, nalleles = cfg["number_of_loci"], cfg["number_of_alleles"]
    if nloci < 1:
        raise ConfigError("number_of_loci must be >= 1")
    if not 2 <= nalleles <= 10:   # mk_mutation_tbl and base_converter assert this
        raise ConfigError("number_of_alleles must be between 2 and 10")
    for name in ("rate_of_mutation", "rate_of_recombination"):
        if not 0.0 <= cfg[name] <= 1.0:
            raise ConfigError(f"{name} must be in [0, 1]")
    # With double mutations off, mk_mutation_tbl gives the unmutated genotype
    # 1 - nloci * rate_of_mutation, which must not go negative.
    if not cfg["allow_double_mutations"] and nloci * cfg["rate_of_mutation"] > 1.0:
        raise ConfigError(
            f"with allow_double_mutations off, number_of_loci * rate_of_mutation "
            f"must be <= 1 (got {nloci} * {cfg['rate_of_mutation']} = "
            f"{nloci * cfg['rate_of_mutation']:g})")
    if not 0.0 < cfg["threshold"] <= 1.0:
        raise ConfigError("threshold must be in (0, 1]")
    if cfg["number_of_generations"] < 1:
        raise ConfigError("number_of_generations must be >= 1")

    ngtypes = nalleles ** nloci
    if ngtypes ** 3 > 50_000_000:
        warn(f"{ngtypes} genotypes: the recombination table has "
             f"{ngtypes * ngtypes * (ngtypes + 1) // 2} doubles")
    cfg["ngtypes"] = ngtypes

    def genotype_array(section, default):
        values = [default] * ngtypes
        table = raw.get(section) or {}
        if not isinstance(table, dict):
            raise ConfigError(f"{section} must be a mapping of genotype: value")
        for key, text in table.items():
            idx = genotype_index(key, nloci, nalleles, section)
            val = to_number(text, float, f"{section}[{key}]")
            if val < 0:
                raise ConfigError(f"{section}[{key}] must be >= 0")
            values[idx] = val
        return values

    cfg["initial_state"] = genotype_array("initial_state", 0.0)
    cfg["fitness"] = genotype_array("fitness", 1.0)
    if sum(cfg["initial_state"]) <= 0:
        raise ConfigError("initial_state must have a positive total population")
    if not any(f > 0 and s > 0 for f, s in zip(cfg["fitness"], cfg["initial_state"])):
        raise ConfigError("every initially present genotype has fitness 0")

    targets = raw["target_genotypes"] or []
    if not isinstance(targets, list):
        raise ConfigError("target_genotypes must be a list")
    idxs = [genotype_index(t, nloci, nalleles, "target_genotypes") for t in targets]
    if len(set(idxs)) != len(idxs):
        raise ConfigError("target_genotypes contains duplicates")
    cfg["target_genotypes"] = idxs

    points = raw.get("report_points") or []
    if not isinstance(points, list):
        raise ConfigError("report_points must be a list")
    pts = [to_number(x, int, "report_points") for x in points]
    if any(x < 0 for x in pts):
        raise ConfigError("report_points must be >= 0")
    clean = sorted(set(pts))
    if clean != pts:
        warn(f"report_points sorted and de-duplicated: {clean}")
    late = [x for x in clean if x >= cfg["number_of_generations"]]
    if late:
        warn(f"report_points {late} are >= number_of_generations "
             f"({cfg['number_of_generations']}); the loop runs i = 0..n-1, so they never print")
    cfg["report_points"] = clean
    return cfg


# ----------------------------------------------------------- C output ----

def c_double(x):
    s = repr(float(x))
    return s if ("." in s or "e" in s) else s + "."


def c_array(values, fmt, per_line=8):
    items = [fmt(v) for v in values]
    rows = [", ".join(items[i:i + per_line]) for i in range(0, len(items), per_line)]
    return "{ " + (",\n          ".join(rows) if rows else "0") + " }"


def indent_block(text, ind):
    return "\n".join(ind + line if line else line for line in text.splitlines())


def gen_params(cfg, ind, _):
    fields = ",\n".join(
        f"    .{name} = {cfg[name] if ctype == 'uint' else c_double(cfg[name])}"
        + ("u" if ctype == "uint" else "")
        for name, (ctype, _) in SCALAR_PARAMS.items())
    return indent_block(f"struct params p = {{ /* yaml */\n{fields}\n}};", ind)


def gen_filled_array(name, ctype, values, fmt, size_expr, ind, original):
    """Keep the user's malloc line, then copy the YAML data into it."""
    body = (f"{original.strip()}\n"
            f"{{ /* yaml: {name} */\n"
            f"    static const {ctype} yaml_{name}[] =\n"
            f"        {c_array(values, fmt)};\n"
            f"    assert(sizeof yaml_{name} == {size_expr});\n"
            f"    memcpy({name}, yaml_{name}, sizeof yaml_{name});\n"
            f"}}")
    return indent_block(body, ind)


def gen_state(cfg, ind, line):
    return gen_filled_array("state", "double", cfg["initial_state"], c_double,
                            "state_size", ind, line)


def gen_fitness(cfg, ind, line):
    return gen_filled_array("fitness", "double", cfg["fitness"], c_double,
                            "state_size", ind, line)


def gen_targets(cfg, ind, line):
    # ngtypes slots; unused slots hold the (uint)-1 terminator check_threshold expects
    vals = cfg["target_genotypes"] + [None] * (cfg["ngtypes"] - len(cfg["target_genotypes"]))
    fmt = lambda v: SENTINEL if v is None else f"{v}u"
    return gen_filled_array("target_genotypes", "uint", vals, fmt,
                            "target_genotypes_size", ind, line)


def gen_report_size(cfg, ind, _):
    k = len(cfg["report_points"]) + 1          # +1 for the terminator
    return f"{ind}uint report_points_size = {k} * sizeof(uint); /* yaml: incl. terminator */"


def gen_report_points(cfg, ind, line):
    vals = cfg["report_points"] + [None]       # terminator: i never equals (uint)-1
    fmt = lambda v: SENTINEL if v is None else f"{v}u"
    return gen_filled_array("report_points", "uint", vals, fmt,
                            "report_points_size", ind, line)


def gen_standalone(name):
    def gen(cfg, ind, _):
        ctype, _conv = STANDALONE[name]
        val = c_double(cfg[name]) if ctype == "double" else f"{cfg[name]}u"
        return f"{ind}{ctype} {name} = {val}; /* yaml */"
    return gen


# Which tagged line gets which generator, keyed by a regex on the line.
RULES = [
    ("params",           re.compile(r"^\s*struct\s+params\s+p\s*;"),                    gen_params),
    ("state",            re.compile(r"^\s*double\s*\*\s*state\s*="),                    gen_state),
    ("fitness",          re.compile(r"^\s*double\s*\*\s*fitness\s*="),                  gen_fitness),
    ("target_genotypes", re.compile(r"^\s*uint\s*\*\s*target_genotypes\s*="),           gen_targets),
    ("report_size",      re.compile(r"^\s*uint\s+report_points_size\s*="),              gen_report_size),
    ("report_points",    re.compile(r"^\s*uint\s*\*\s*report_points\s*="),              gen_report_points),
] + [
    (name, re.compile(rf"^\s*{re.escape(ctype)}\s+{name}\s*;"), gen_standalone(name))
    for name, (ctype, _) in STANDALONE.items()
]
TAG = re.compile(r"//\s*(yaml|ymal|TODO)\b")   # 'ymal' is the typo on `struct params p`


def fill_template(src, cfg):
    out, done = [], set()
    for lineno, line in enumerate(src.splitlines(), 1):
        rule = next((r for r in RULES if r[1].search(line)), None)
        if rule and TAG.search(line):
            name, _, gen = rule
            if name in done:
                raise ConfigError(f"template line {lineno}: {name} appears twice")
            ind = re.match(r"\s*", line).group(0)
            out.append(gen(cfg, ind, line))
            done.add(name)
        else:
            if re.search(r"//\s*yaml\b", line) and not line.lstrip().startswith(("/*", "*")):
                raise ConfigError(f"template line {lineno}: tagged '// yaml' but not recognised:\n  {line.strip()}")
            out.append(line)
    missing = {r[0] for r in RULES} - done
    if missing:
        raise ConfigError(f"template is missing tagged line(s) for: {', '.join(sorted(missing))}")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("yaml_file")
    ap.add_argument("c_template")
    ap.add_argument("-o", "--output", help="output .c file (default: stdout)")
    args = ap.parse_args()

    try:
        cfg = load_config(args.yaml_file)
        with open(args.c_template) as fh:
            result = fill_template(fh.read(), cfg)
    except (ConfigError, yaml.YAMLError) as e:
        sys.exit(f"yaml2c: error: {e}")

    if args.output:
        with open(args.output, "w") as fh:
            fh.write(result)
    else:
        sys.stdout.write(result)


if __name__ == "__main__":
    main()
