#!/usr/bin/env python3
"""Checks tests/proof's Rocq proofs against Flocq: rounding.v, the IEEE 754 facts, integers.v,
the whole-number arithmetic, and angle.v, the exact reductions of angles held in degrees.

    python3 -I tests/proof/rocq.py ROCQ FLOCQ_SOURCE FOLDER

Builds Flocq from its source, as fetched by vendor/flocq, then compiles each proof, whose every
theorem must be proved: the file holds no Admitted and no admit, and Print Assumptions names no
axiom but the four of Rocq's classical reals that Flocq and every real analysis in Rocq rest on.
The exit status is 1 when a check fails.
"""

import os
import pathlib
import re
import subprocess
import sys

ALLOWED = {"ClassicalDedekindReals.sig_forall_dec", "ClassicalDedekindReals.sig_not_dec",
           "Classical_Prop.classic", "FunctionalExtensionality.functional_extensionality_dep"}

HERE = pathlib.Path(__file__).resolve().parent
PROOFS = ("rounding", "integers", "angle")


def build_flocq(rocq, source):
    """Builds Flocq in its source folder, once; make is incremental after."""
    version = source / "src" / "Version.v"
    if not version.exists():
        template = (source / "src" / "Version.v.in").read_text()
        version.write_text(template.replace("@PACKAGE_VERSION_NUM@", "40202")
                           .replace("@PACKAGE_VERSION@", "4.2.2"))
    files = sorted(str(path.relative_to(source)) for path in (source / "src").rglob("*.v"))
    subprocess.run([rocq, "makefile", "-f", "_CoqProject", *files, "-o", "Makefile.coq"],
                   cwd=source, check=True, capture_output=True)
    jobs = str(os.cpu_count() or 1)
    for _ in range(3):
        result = subprocess.run(["make", "-f", "Makefile.coq", "-j", jobs], cwd=source,
                                capture_output=True, text=True, check=False)
        if result.returncode == 0:
            return True, ""
    return False, result.stderr.strip().splitlines()[-1] if result.stderr else "make failed"


def main():
    rocq, source, folder = sys.argv[1], pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
    failed = False

    def report(name, held, told=""):
        nonlocal failed
        failed = failed or not held
        print(f"{name:40} {'held' if held else 'BROKEN':7} {told}", flush=True)

    built, told = build_flocq(rocq, source)
    report("Flocq builds", built, told)
    if not built:
        return 1
    folder.mkdir(parents=True, exist_ok=True)
    for name in PROOFS:
        proof = HERE / f"{name}.v"
        text = proof.read_text()
        report(f"{name}: no theorem is admitted", not re.search(r"\b(Admitted|admit)\b", text))
        result = subprocess.run([rocq, "compile", "-noglob", "-R", str(source / "src"), "Flocq",
                                 "-o", str(folder / f"{name}.vo"), str(proof)],
                                cwd=folder, capture_output=True, text=True, check=False)
        report(f"Rocq accepts {name}.v", result.returncode == 0,
               next((line for line in result.stdout.splitlines() + result.stderr.splitlines()
                     if "Error" in line), ""))
        theorems = re.findall(r"^Print Assumptions (\w+)\.", text, re.M)
        blocks = re.split(r"^(?:Axioms:|Closed under the global context)\s*$", result.stdout,
                          flags=re.M)[1:]
        report(f"{name}: every theorem's assumptions are printed", len(blocks) == len(theorems),
               f"{len(blocks)} of {len(theorems)}")
        for theorem, block in zip(theorems, blocks):
            axioms = set(re.findall(r"^([\w.]+) :", block, re.M))
            report(f"{theorem}: only Rocq's real axioms", axioms <= ALLOWED,
                   ", ".join(sorted(axioms - ALLOWED)))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
