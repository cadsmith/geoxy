#!/usr/bin/env python3
"""Checks that the proved computations are the code's: program.py's graphs against the
compiler's LLVM IR.

    python3 -I tests/proof/translation.py [clang++]

The inline first stages are compiled from a probe, Orient and InCircle on doubles; the later
stages from source/exact/circle.cpp, orient.cpp and batch.cpp themselves, at -O2 with geoxy's
flags. The IR's floating point values are read back into graphs: its fadd, fsub, fmul, fneg,
fabs and fma, its double constants to the bit, its inputs the probe's arguments or the loads of
the corners a function is given. A load from anywhere else is left opaque, and an opaque value
equals no graph's, so a doubt fails the check rather than passing it. Every value the predicates
compare against a bound, and every bound, must then equal program.py's, commutative operands in
either order: the first stages' tests and the signs they return, Compensated's estimate against
both its thresholds, the third stage's sum against its bound, CrossSign's two fused products
against each other with Safe's bounds and the roundoffs' tests for zero, and the batch loop's
first stage, point by point, with the corners it hands the later stages. No floating point
operation may carry a fast math flag.

The exit status is 1 when a check fails.
"""

import pathlib
import re
import shutil
import struct
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
import program as P

ROOT = pathlib.Path(__file__).resolve().parents[2]

PROBE = """#include "exact/circle.hpp"
#include "exact/orient.hpp"

extern "C" int probe_orient(double ax, double ay, double bx, double by, double cx, double cy)
{
  return Geoxy::Exact::Orient({ax, ay}, {bx, by}, {cx, cy});
}

extern "C" int probe_circle(double ax, double ay, double bx, double by, double cx, double cy,
                            double dx, double dy)
{
  return Geoxy::Exact::InCircle({{ax, ay}, {bx, by}, {cx, cy}}, {dx, dy});
}
"""

FAST = re.compile(r"\b(fast|nnan|ninf|nsz|arcp|contract|afn|reassoc)\b")
CORNER_OFFSETS = {0: "ax", 8: "ay", 16: "bx", 24: "by", 32: "cx", 40: "cy"}
LINE_OFFSETS = {"line": {0: "ax", 8: "ay", 16: "bx", 24: "by"}, "point": {0: "cx", 8: "cy"}}


def compile_ir(compiler, source, folder):
    target = pathlib.Path(folder) / (source.stem + ".ll")
    command = [compiler, "-std=c++23", "-O2", "-fno-slp-vectorize", "-fno-vectorize",
               "-ffp-contract=off", "-fno-discard-value-names",
               "-DGEOXY_STATIC", "-I", str(ROOT / "source"), "-I", str(ROOT / "include"),
               "-S", "-emit-llvm", "-o", str(target), str(source)]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return target.read_text()


def function(ir, pattern):
    """The text of the first function whose name matches."""
    match = re.search(r"^define [^\n]*@" + pattern + r"\(.*?^}", ir, re.S | re.M)
    if not match:
        raise SystemExit(f"no function {pattern} in the IR")
    return match.group(0)


def double(text):
    """An IR double constant, its bits exactly."""
    if text.startswith(("0x", "f0x")):
        digits = text.split("x", 1)[1]
        return struct.unpack(">d", bytes.fromhex(digits.rjust(16, "0")))[0]
    return float(text)


class Graph:
    """A function's IR read back as values: each SSA name its definition."""

    def __init__(self, body, inputs, loads=None):
        self.lines = {}
        self.inputs = inputs
        self.loads = loads if loads is not None else {"circle": CORNER_OFFSETS}
        self.fast = []
        self.compares = []
        self.memo = {}
        for line in body.splitlines():
            definition = re.match(r"\s*%([\w.]+) = (.*)$", line)
            if definition:
                self.lines[definition.group(1)] = definition.group(2)
                if re.match(r"(fadd|fsub|fmul|fneg|fdiv)\b", definition.group(2)) and \
                        FAST.search(definition.group(2).split("double")[0]):
                    self.fast.append(line.strip())
                compare = re.match(r"fcmp (\w+) double (\S+), (\S+)", definition.group(2))
                if compare:
                    self.compares.append((definition.group(1), compare.group(1),
                                          compare.group(2), compare.group(3)))

    def operand(self, text):
        text = text.strip().rstrip(",")
        if text.startswith("%"):
            return self.value(text[1:])
        return P.constant(double(text))

    def address(self, text):
        """A pointer as (base, offset), following constant geps."""
        name = text.strip().lstrip("%")
        definition = self.lines.get(name)
        if definition and re.match(r"getelementptr inbounds (?:nuw )?(?:\[16 x i8\]|"
                                   r"%\"struct\.Geoxy::Point\"), ptr %[\w.]+, i64 %", definition):
            return "point", 0
        if definition:
            gep = re.match(r"getelementptr inbounds (?:nuw )?i8, ptr %([\w.]+), i64 (\d+)",
                           definition)
            if gep:
                base, offset = self.address(gep.group(1))
                return base, offset + int(gep.group(2))
        return name, 0

    def value(self, name):
        if name in self.memo:
            return self.memo[name]
        if name in self.inputs:
            result = P.given(self.inputs[name])
        else:
            result = self.read(name, self.lines.get(name, ""))
        self.memo[name] = result
        return result

    def read(self, name, text):
        binary = re.match(r"(fadd|fsub|fmul) (?:[a-z ]+ )?double (\S+), (\S+)", text)
        if binary:
            op = {"fadd": "add", "fsub": "sub", "fmul": "mul"}[binary.group(1)]
            return P.Value(op, self.operand(binary.group(2)), self.operand(binary.group(3)))
        negation = re.match(r"fneg (?:[a-z ]+ )?double (\S+)", text)
        if negation:
            return P.neg(self.operand(negation.group(1)))
        call = re.match(r"(?:tail )?call (?:[a-z ]+ )?double @llvm\.(fabs|fma)\.f64\((.*)\)", text)
        if call:
            args = [self.operand(arg.split()[-1]) for arg in call.group(2).split(",")]
            if call.group(1) == "fabs":
                return P.magnitude(args[0])
            return P.fma(*args)
        cast = re.match(r"bitcast i64 %([\w.]+) to double", text)
        if cast:
            text = self.lines.get(cast.group(1), "").replace("load i64,", "load double,", 1)
        load = re.match(r"load double, ptr (%[\w.]+)", text)
        if load:
            base, offset = self.address(load.group(1))
            if offset in self.loads.get(base, {}):
                return P.given(self.loads[base][offset])
        extract = re.match(r"extractvalue \[2 x double\] %point\.coerce, (\d)", text)
        if extract:
            return P.given("dx" if extract.group(1) == "0" else "dy")
        return P.given(f"opaque {name}")

    def compared(self):
        """Every comparison as (predicate, value, value)."""
        return [(predicate, self.operand(left), self.operand(right))
                for _, predicate, left, right in self.compares]


def against(graph, value, bound):
    """Whether the IR compares the value, or its negation, with the bound, or its negation: the
    forms Decided and the strict test take after the compiler."""
    values = {value.key(), P.neg(value).key()}
    bounds = {bound.key(), P.neg(bound).key()}
    return any({left.key(), right.key()} & values and {left.key(), right.key()} & bounds
               for _, left, right in graph.compared())


def signs(graph, value):
    """Whether the IR takes the value's sign: compared with zero both ways."""
    zero = P.constant(0.0).key()
    found = {predicate for predicate, left, right in graph.compared()
             if left.key() == value.key() and right.key() == zero}
    return {"ogt", "olt"} <= found


def handed_corners(body, graph):
    """Whether the triangle OrientSigns hands OrientAdapt is the line's two corners, copied, then
    the point: memcpy of 16 bytes from the line and from 16 bytes in, and the point's two
    coordinates stored 32 and 40 bytes in."""
    call = re.search(r"call (?:noundef )?i32 @_ZN5Geoxy5Exact11OrientAdapt\w*\(ptr [^%]*%([\w.]+)\)",
                     body)
    if not call:
        return False
    triangle = call.group(1)
    copies = re.findall(r"llvm\.memcpy\S*\(ptr [^%]*%([\w.]+), ptr [^%]*%([\w.]+), i64 16", body)
    copied = {(graph.address(target), graph.address(source)) for target, source in copies}
    stores = re.findall(r"store double %([\w.]+), ptr %([\w.]+)", body)
    stored = {graph.address(target): graph.value(value) for value, target in stores}
    return ({((triangle, 0), ("line", 0)), ((triangle, 16), ("line", 16))} <= copied
            and stored.get((triangle, 32)) == P.given("cx")
            and stored.get((triangle, 40)) == P.given("cy"))


def main():
    compiler = sys.argv[1] if len(sys.argv) > 1 else shutil.which("clang++")
    failed = False

    def report(name, held, told=""):
        nonlocal failed
        failed = failed or not held
        print(f"{name:62} {'held' if held else 'BROKEN':7} {told}", flush=True)

    with tempfile.TemporaryDirectory() as folder:
        probe = pathlib.Path(folder) / "probe.cpp"
        probe.write_text(PROBE)
        probe_ir = compile_ir(compiler, probe, folder)
        circle_ir = compile_ir(compiler, ROOT / "source" / "exact" / "circle.cpp", folder)

        names = {n: n for n in ("ax", "ay", "bx", "by", "cx", "cy")}
        orient = Graph(function(probe_ir, "probe_orient"), names)
        model = P.orient_first()
        decided = model["decided"]
        report("orient: the first stage's test is the proved one",
               against(orient, P.magnitude(model["det"]), model["bound"])
               and any(p == "ogt" and P.greater(l, r) == decided
                       for p, l, r in orient.compared()))
        report("orient: its sign is the proved determinant's", signs(orient, model["det"]))
        report("orient: no fast math", not orient.fast, "; ".join(orient.fast[:2]))

        circle = Graph(function(probe_ir, "probe_circle"), {**names, "dx": "dx", "dy": "dy"})
        model = P.circle_first()
        report("circle: the first stage's test is the proved one",
               any(p == "ogt" and P.greater(l, r) == model["decided"]
                   for p, l, r in circle.compared()))
        report("circle: its sign is the proved determinant's", signs(circle, model["det"]))
        report("circle: no fast math", not circle.fast, "; ".join(circle.fast[:2]))
        handed = [line for line in circle.lines.values() if "CircleAdapt" in line]
        passed = [re.search(r"double noundef (%[\w.]+)\)", line) for line in handed]
        report("circle: the later stages get the first stage's permanent",
               bool(passed) and all(m and circle.operand(m.group(1)) == model["perm"]
                                    for m in passed))
        for name in ("lift", "m", "t", "pm", "q"):
            report(f"circle: the corners' {name} are one turned", P.corners_alike(model, name))

        adapt = Graph(function(circle_ir, r"_ZN5Geoxy5Exact11CircleAdapt\w*"),
                      {"permanent": "permanent"})
        model = P.compensated()
        report("compensated: the estimate meets the exact threshold",
               against(adapt, model["E"], model["exact bound"]))
        report("compensated: the estimate meets the second threshold",
               against(adapt, model["E"], model["second bound"]))
        report("compensated: the third stage meets its bound",
               against(adapt, model["C"], model["third bound"]))
        report("compensated: no fast math", not adapt.fast, "; ".join(adapt.fast[:2]))
        for name in ("hB", "Lo", "mB", "Mlo", "p", "lo", "corr"):
            report(f"compensated: the corners' {name} are one turned",
                   P.corners_alike(model, name))

        orient_ir = compile_ir(compiler, ROOT / "source" / "exact" / "orient.cpp", folder)
        cross = Graph(function(orient_ir, r"_ZN5Geoxy5Exact11OrientAdapt\w*"), {},
                      {"corners": CORNER_OFFSETS})
        model = P.orient_adapt()
        compared = cross.compared()
        rest, roundoff = model["rest"], model["roundoff"]
        report("cross sign: the fused products are compared, both ways",
               any(p == "ogt" and (l, r) == (rest, roundoff) for p, l, r in compared)
               and any((p == "olt" and (l, r) == (rest, roundoff))
                       or (p == "uge" and (l, r) == (rest, roundoff)) for p, l, r in compared))
        zero = P.constant(0.0)
        report("cross sign: every difference's roundoff is tested for zero",
               all(any(p in ("oeq", "une") and l == model[f"{d} tail"] and r == zero
                       for p, l, r in compared) for d in ("fx", "fy", "sx", "sy")))
        report("cross sign: Safe bounds both products below and above",
               all(any(p == "oge" and l == model[f"safe {which}"]
                       and r == P.constant(P.SAFE_LOW) for p, l, r in compared)
                   and any(p == "ole" and l == model[f"safe {which}"]
                           and r == P.constant(P.SAFE_HIGH) for p, l, r in compared)
                   for which in ("first", "second")))
        report("cross sign: no fast math", not cross.fast, "; ".join(cross.fast[:2]))

        batch_ir = compile_ir(compiler, ROOT / "source" / "exact" / "batch.cpp", folder)
        body = function(batch_ir, r"_ZN5Geoxy5Exact11OrientSigns\w*")
        batch = Graph(body, {}, LINE_OFFSETS)
        model = P.orient_first()
        report("batch: each point's first stage is the proved one",
               against(batch, P.magnitude(model["det"]), model["bound"]))
        report("batch: its sign is the proved determinant's", signs(batch, model["det"]))
        report("batch: no fast math", not batch.fast, "; ".join(batch.fast[:2]))
        report("batch: the later stages get the line's corners and the point",
               handed_corners(body, batch))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
