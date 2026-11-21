"""Build and run every C++ program in the docs, and hold each to the output written under it.

A program is a fence opened as ```cpp [main.cpp]; the fence right after it opened as
```text [Output] is what it must print. Each is compiled against the library of a build folder
and run, by the compiler CXX names, or c++.

usage: python3 -I tools/docs/examples.py [--build build/release] [--write] [page.md ...]
       --write puts what each program prints in its Output fence instead of checking it.
"""

import argparse
import os
import pathlib
import re
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
PROGRAM = re.compile(r"```cpp \[main\.cpp\]\n(.*?)```\n(\s*```text \[Output\]\n(.*?)```\n)?", re.S)


def run(source, build, folder):
    """What a program prints, or None with the reason on stderr where it does not build or run."""
    path = folder / "main.cpp"
    path.write_text(source)
    binary = folder / "main"
    made = subprocess.run(
        [os.environ.get("CXX", "c++"), "-std=c++23", "-O1", "-DGEOXY_STATIC", "-I",
         str(ROOT / "include"), str(path), str(build / "lib" / "libgeoxy.a"), "-o", str(binary)],
        capture_output=True, text=True)
    if made.returncode != 0:
        sys.stderr.write(made.stderr)
        return None
    ran = subprocess.run([str(binary)], capture_output=True, text=True, timeout=60)
    if ran.returncode != 0:
        sys.stderr.write(f"exit {ran.returncode}\n{ran.stderr}")
        return None
    return ran.stdout


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", default="build/release")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("pages", nargs="*")
    options = parser.parse_args()
    build = (ROOT / options.build).resolve()
    pages = [pathlib.Path(page) for page in options.pages] or sorted((ROOT / "docs").rglob("*.md"))
    programs = failures = 0
    for page in pages:
        text = page.read_text()
        pieces, last = [], 0
        for index, found in enumerate(PROGRAM.finditer(text)):
            programs += 1
            name = f"{page.relative_to(ROOT) if page.is_absolute() else page} program {index + 1}"
            with tempfile.TemporaryDirectory() as folder:
                printed = run(found.group(1), build, pathlib.Path(folder))
            if printed is None:
                failures += 1
                print(f"FAILED {name}: does not build or run")
                continue
            if options.write:
                pieces.append(text[last:found.start()])
                pieces.append(f"```cpp [main.cpp]\n{found.group(1)}```\n\n```text [Output]\n{printed}```\n")
                last = found.end()
            elif found.group(3) != printed:
                failures += 1
                print(f"FAILED {name}: prints\n{printed}--- but the page says\n{found.group(3)}")
        if options.write and pieces:
            page.write_text("".join(pieces) + text[last:])
    print(f"{programs} programs, {failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
