#!/usr/bin/env python3
"""Runs every check a change must pass, side by side, and says how long each took.

    python3 -I tools/gate/run.py [file ...]

The debug, release and sanitize presets are each built and tested; the coverage preset's `coverage`
target is built, which gates lines, functions, regions and branches at 100%; the `tidy` preset is
built, which lints the library through clang-tidy; the `bench` preset is built, so that
`geoxy_bench` builds; the `unity` preset is built, configured first where it is not, so that no two
files' helpers of one name meet when a unity build compiles them as one; and `geoxy_docs` and
`geoxy_api` are built in the debug preset once its tests are done. Each preset has a build folder of
its own, so they run at once: on ten cores the whole gate takes under twenty seconds where nothing
but a file or two changed, as every step is incremental, each test case runs on its own, and the
proofs of `tests/proof` and the Python checks of `tools/suites` run again only when what they check
changes, once for every preset.

Files named are also checked for their formatting, by tools/format/style.py, and for what an
editor running clangd shows, by tools/format/editor.py: name the C++ files a change touches.

Each step prints one line: how it ended, how long it took, and what it was. A step that fails
prints the lines that tell of a failure and the end of what it wrote. The exit status is 1 if any
step failed.
"""

import pathlib
import subprocess
import sys
import threading
import time

ROOT = pathlib.Path(__file__).resolve().parents[2]
TAIL = 25
TELLING = ("***", "(Failed)", "FAILED", "error:", "ERROR:", "runtime error")


def run(command):
    """A command's exit status and what it wrote, run in the repository."""
    done = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=False)
    return done.returncode, done.stdout + done.stderr


def chain(name, commands, results):
    """Runs commands one after another, stopping at the first that fails, and files how the
    chain ended under its name."""
    began = time.time()
    for command in commands:
        status, written = run(command)
        if status != 0:
            results[name] = (False, time.time() - began, " ".join(command), written)
            return
    results[name] = (True, time.time() - began, "", "")


def chains(files):
    """The gate's chains of commands, each with a build folder of its own."""
    preset = ["cmake", "--build", "--preset"]
    steps = {
        "debug": [preset + ["debug"], ["ctest", "--preset", "debug"],
                  preset + ["debug", "--target", "geoxy_docs"],
                  preset + ["debug", "--target", "geoxy_api"]],
        "release": [preset + ["release"], ["ctest", "--preset", "release"]],
        "sanitize": [preset + ["sanitize"], ["ctest", "--preset", "sanitize"]],
        "coverage": [preset + ["coverage", "--target", "coverage"]],
        "tidy": [preset + ["tidy"]],
        "bench": [preset + ["bench"]],
        "unity": [preset + ["unity"]],
    }
    if not (ROOT / "build" / "unity" / "CMakeCache.txt").exists():
        steps["unity"].insert(0, ["cmake", "--preset", "unity"])
    sources = [name for name in files if name.endswith((".cpp", ".hpp"))]
    if sources:
        python = [sys.executable, "-I"]
        steps["format"] = [python + ["tools/format/style.py", "--check"] + sources]
        steps["editor"] = [python + ["tools/format/editor.py"] + sources]
    return steps


def main():
    """Runs the chains at once and prints how each ended."""
    began = time.time()
    results = {}
    threads = [threading.Thread(target=chain, args=(name, commands, results))
               for name, commands in chains(sys.argv[1:]).items()]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    failed = False
    for name, (passed, took, command, written) in sorted(results.items()):
        print(f"{'ok    ' if passed else 'FAILED'} {took:6.1f} s  {name}")
        if not passed:
            failed = True
            print(f"       {command}")
            lines = written.strip().splitlines()
            told = [line for line in lines[:-TAIL] if any(word in line for word in TELLING)]
            for line in told[:TAIL] + lines[-TAIL:]:
                print(f"       {line}")
    print(f"{'FAILED' if failed else 'ok    '} {time.time() - began:6.1f} s  in all")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
