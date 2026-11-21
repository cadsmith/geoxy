#!/usr/bin/env python3
"""Reads the machine code compilers emit for the exact predicates, on arm64 and x86-64.

    python3 -I tests/proof/machine.py [clang++ ...]

The proofs take every operation to be rounded on its own, to a double, in the order written
(source/exact/strict.hpp). This compiles every source of source/exact, and a probe of the inline
first stages, with each compiler named (Apple's and LLVM's Clang and Homebrew's GCC by default):
Clang for arm64, x86-64, and x86-64 with FMA, GCC for the machine it runs on, each at -O0, -O2
and -O3, with geoxy's flags and with none; disassembles each object with llvm-objdump; and checks:

  no fused multiply-add, scalar or vector, where the source asks for none: everything built with
  GEOXY_EXACT_FMA 0, which takes std::fma out, and the first stages in every build;
  no x87 instruction on x86-64, whose 80 bit registers would round twice;
  that the reading sees a fusion where there is one: a * b + c fuses without the strict header,
  and does not with it;
  that fast math stops the build.

The exit status is 1 when a check fails.
"""

import concurrent.futures
import os
import pathlib
import platform
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCES = sorted((ROOT / "source" / "exact").glob("*.cpp"))
LLVM = pathlib.Path("/opt/homebrew/opt/llvm/bin")
HOMEBREW = pathlib.Path("/opt/homebrew/bin")
ARCHES = {"arm64": ["-arch", "arm64"],
          "x86-64": ["-arch", "x86_64"],
          "x86-64 fma": ["-arch", "x86_64", "-mfma"]}
LEVELS = ("-O0", "-O2", "-O3")
FLAGS = {"geoxy's flags": ["-ffp-contract=off"], "no flags": []}

PROBE = """#include "exact/circle.hpp"
#include "exact/orient.hpp"
#include "exact/triangle.hpp"

#include "geoxy/point.hpp"

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

CONTROL = "double fused(double a, double b, double c) { return a * b + c; }\n"

FUSED = re.compile(r"^(v?fn?m(add|sub)\w*|fml[as])$")
X87 = re.compile(r"^f[a-z0-9]*$")


def objdump():
    """LLVM's disassembler: Homebrew's, the path's, or the one Xcode ships."""
    found = LLVM / "llvm-objdump"
    if found.exists():
        return [str(found)]
    if shutil.which("llvm-objdump"):
        return ["llvm-objdump"]
    return ["xcrun", "objdump"]


def compilers():
    named = sys.argv[1:]
    if named:
        return named
    found = [path for path in ("/usr/bin/clang++", str(LLVM / "clang++"))
             if pathlib.Path(path).exists()]
    found += sorted(str(path) for path in HOMEBREW.glob("g++-[0-9]*"))[-1:]
    return found or [shutil.which("clang++")]


def gnu(compiler):
    """Whether the compiler is GCC, which builds for the machine it runs on alone."""
    result = subprocess.run([compiler, "--version"], capture_output=True, text=True, check=False)
    return "Free Software Foundation" in result.stdout


def arches(compiler):
    """The targets a compiler builds for: Clang's three; GCC's own, with no -arch."""
    if not gnu(compiler):
        return ARCHES
    native = "arm64" if platform.machine() == "arm64" else "x86-64"
    return {native: []}


def sysroot():
    result = subprocess.run(["xcrun", "--show-sdk-path"], capture_output=True, text=True,
                            check=False)
    return ["-isysroot", result.stdout.strip()] if result.returncode == 0 else []


def compile_one(compiler, source, target, flags):
    """Compiles a source to an object; returns the object, or the compiler's complaint."""
    command = [compiler, "-std=c++23", "-c", str(source), "-o", str(target),
               "-I", str(ROOT / "source"), "-I", str(ROOT / "include"), "-DGEOXY_STATIC",
               *sysroot(), *flags]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    return (target, "") if result.returncode == 0 else (None, result.stderr)


def instructions(target):
    """Each function's mnemonics, from the disassembly."""
    text = subprocess.run([*objdump(), "-d", "--no-show-raw-insn", str(target)],
                          capture_output=True, text=True, check=True).stdout
    functions = {}
    current = None
    for line in text.splitlines():
        label = re.match(r"^[0-9a-f]+ <(.+)>:$", line)
        if label:
            current = label.group(1)
            functions[current] = []
            continue
        code = re.match(r"^\s*[0-9a-f]+:\s+(\S+)", line)
        if code and current is not None:
            functions[current].append(code.group(1).split(".")[0].lower())
    return functions


def offending(functions, pattern, only=None):
    found = []
    for name, mnemonics in functions.items():
        if only is not None and not only(name):
            continue
        found += [f"{name}: {mnemonic}" for mnemonic in mnemonics if pattern.match(mnemonic)]
    return found


def audit(job):
    """Builds and reads one source in one configuration: its fusions, its first stages' fusions
    and its x87 instructions, or the compiler's complaint."""
    compiler, arch, common, source, forced, target = job
    extra = ["-DGEOXY_EXACT_FMA=0"] if forced else []
    built, complaint = compile_one(compiler, source, target, common + extra)
    if built is None:
        return None, complaint
    functions = instructions(built)
    fused = offending(functions, FUSED) if forced else []
    first = offending(functions, FUSED, lambda name: "probe" in name)
    x87 = offending(functions, X87) if arch.startswith("x86") else []
    return (fused, first, x87), ""


def main():
    failed = False

    def report(name, held, told=""):
        nonlocal failed
        failed = failed or not held
        print(f"{name:58} {'held' if held else 'BROKEN':7} {told}", flush=True)

    with tempfile.TemporaryDirectory() as folder:
        work = pathlib.Path(folder)
        probe = work / "probe.cpp"
        probe.write_text(PROBE)
        control = work / "control.cpp"
        control.write_text(CONTROL)
        strict = work / "strict.cpp"
        strict.write_text('#include "exact/strict.hpp"\n' + CONTROL)

        for compiler in compilers():
            tag = ("gcc" if gnu(compiler) else "llvm" if "homebrew" in compiler
                   else pathlib.Path(compiler).name)
            native = [] if gnu(compiler) else ["-arch", "arm64"]
            configurations = [(arch, level, flagged, [*arch_flags, level, *flags])
                              for arch, arch_flags in arches(compiler).items() for level in LEVELS
                              for flagged, flags in FLAGS.items()]
            jobs = {}
            for number, (arch, level, flagged, common) in enumerate(configurations):
                for index, source in enumerate([*SOURCES, probe]):
                    for forced in (True, False):
                        target = work / f"{number}-{index}-{int(forced)}.o"
                        jobs[(number, index, forced)] = (compiler, arch, common, source, forced,
                                                         target)
            with concurrent.futures.ThreadPoolExecutor(os.cpu_count()) as pool:
                results = dict(zip(jobs, pool.map(audit, jobs.values())))
            for number, (arch, level, flagged, _) in enumerate(configurations):
                name = f"{tag} {arch} {level} {flagged}"
                fused, first, x87 = [], [], []
                for key, (found, complaint) in results.items():
                    if key[0] != number:
                        continue
                    if found is None:
                        report(f"{name}: builds", False, complaint.splitlines()[0])
                        continue
                    fused += found[0]
                    first += found[1]
                    x87 += found[2]
                report(f"{name}: no fusion", not fused, "; ".join(fused[:3]))
                report(f"{name}: first stages unfused", not first, "; ".join(first[:3]))
                if arch.startswith("x86"):
                    report(f"{name}: no x87", not x87, "; ".join(x87[:3]))

            built, _ = compile_one(compiler, control, work / "control.o", [*native, "-O2"])
            seen = built is not None and offending(instructions(built), FUSED)
            report(f"{tag} the reading sees a * b + c fused", bool(seen))
            built, _ = compile_one(compiler, strict, work / "strict.o", [*native, "-O2"])
            kept = built is not None and not offending(instructions(built), FUSED)
            report(f"{tag} strict.hpp keeps it unfused", bool(kept))
            built, complaint = compile_one(compiler, probe, work / "fast.o", ["-ffast-math"])
            report(f"{tag} fast math stops the build", built is None and "fast" in complaint)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
