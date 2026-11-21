"""Lemmas for Z3 about the exact predicates' floating point arithmetic, and their runner.

A double is modelled as the real it holds. A rounded sum or difference r of an exact value z
obeys |r - z| <= u |r| and |r - z| <= u |z|, u = 2^-53, sums being exact where they would
underflow. A rounded product obeys them where |z| is at least the smallest normal double,
2^-1022; below it, where gradual underflow rounds to a multiple of 2^-1074, it is off by at most
2^-1075 and at most 2^-1022 itself; either way it keeps the exact value's sign, and zero stays
zero. An error-free transformation adds that the double and its roundoff sum to the exact value.
A rounding is monotone. Overflow is left to the lemmas in IEEE 754's own arithmetic, Z3's
floating point theory.

None of that is assumed: tests/proof/rounding.v proves each in Rocq with Flocq, for binary64
rounding to nearest, ties to even, and that IEEE 754's addition, subtraction, multiplication and
fused multiply-add are that rounding of the exact value: sum_relative, rounded_relative,
multiplied, monotone, symmetric; two_sum, two_diff_tail, fast_two_sum, two_product_fma and
two_product_split, each for source/exact/two.hpp's operations in its own order; ieee_plus,
ieee_minus, ieee_mult and ieee_fma.

Each lemma's conclusions are the next one's facts, by name. Z3 must find every conclusion's
negation unsatisfiable, and the facts alone satisfiable, so that no lemma holds for want of a
case. And every rounding a lemma states must be the code's: each lemma binds its variables to
program.py's graphs of the code, which translation.py checks against the compiler's IR, and the
runner checks each rounded variable is the graph's operation on the graph's operands, each
roundoff the graph's own.
"""

import fractions
import pathlib
import shutil
import subprocess
import sys
import tempfile
import time

EPSILON = 2.0 ** -53
NORMAL = 2.0 ** -1022
ORIENT_BOUND_A = (3.0 + (16.0 * EPSILON)) * EPSILON
CIRCLE_BOUND_A = (10.0 + (96.0 * EPSILON)) * EPSILON
CIRCLE_BOUND_B = (4.0 + (128.0 * EPSILON)) * EPSILON
CIRCLE_BOUND_C = (44.0 + (576.0 * EPSILON)) * EPSILON * EPSILON
CIRCLE_BOUND_COMPENSATED = 64.0 * EPSILON * EPSILON
CIRCLE_BOUND_UNDERFLOW = 2.0 ** -1069
RESULT_BOUND = (3.0 + (8.0 * EPSILON)) * EPSILON

PRELUDE = f"""(define-fun absr ((x Real)) Real (ite (>= x 0.0) x (- x)))
(define-fun u () Real (/ 1.0 9007199254740992.0))
(define-fun uN () Real (/ 1.0 {2 ** 1022}.0))
(define-fun h () Real (/ 1.0 {2 ** 1075}.0))
"""


def rational(value):
    """A double as an exact SMT-LIB rational."""
    exact = fractions.Fraction(value)
    return f"(/ {exact.numerator}.0 {exact.denominator}.0)"


def power(base, exponent):
    """base^exponent as an SMT-LIB product."""
    return "(* " + " ".join([base] * exponent) + ")" if exponent > 1 else base


def growth(exponent):
    """(1 + u)^exponent - 1, the error of that many roundings in a row."""
    return f"(- {power('(+ 1.0 u)', exponent)} 1.0)"


class Lemma:
    """Facts about reals and a conclusion, which Z3 must prove."""

    def __init__(self, name, *reals):
        self.name = name
        self.reals = list(reals)
        self.facts = []
        self.goals = []
        self.operations = []
        self.instances = []

    def fact(self, text):
        self.facts.append(text)

    def rounded(self, result, exact):
        """result is a sum or difference exact rounded to nearest, or a product clear of
        underflow."""
        self.operations.append(("rounded", result, exact))
        self.fact(f"(<= (absr (- {result} {exact})) (* u (absr {result})))")
        self.fact(f"(<= (absr (- {result} {exact})) (* u (absr {exact})))")

    def multiplied(self, result, exact):
        """result is a product exact rounded to nearest, gradual underflow allowed."""
        self.operations.append(("rounded", result, exact))
        relative = (f"(and (>= (absr {exact}) uN)"
                    f" (<= (absr (- {result} {exact})) (* u (absr {result})))"
                    f" (<= (absr (- {result} {exact})) (* u (absr {exact}))))")
        absolute = (f"(and (< (absr {exact}) uN) (<= (absr (- {result} {exact})) h)"
                    f" (<= (absr {result}) uN))")
        self.fact(f"(or {relative} {absolute})")
        self.fact(f"(=> (>= {exact} 0.0) (>= {result} 0.0))")
        self.fact(f"(=> (<= {exact} 0.0) (<= {result} 0.0))")

    def split(self, big, small, exact):
        """big + small is exact, big rounded to nearest: TwoSum, TwoDiff, TwoProduct, Square."""
        self.fact(f"(= (+ {big} {small}) {exact})")
        self.fact(f"(<= (absr (- {big} {exact})) (* u (absr {big})))")
        self.fact(f"(<= (absr (- {big} {exact})) (* u (absr {exact})))")
        self.operations.append(("split", big, small, exact))

    def decided(self, value, bound):
        """Decided(value, bound): value >= bound or -value >= bound."""
        self.fact(f"(or (>= {value} {bound}) (>= (- {value}) {bound}))")

    def exceeds(self, value, bound):
        """|value| > bound, as the first stages test it."""
        self.fact(f"(> (absr {value}) {bound})")

    def goal(self, text):
        self.goals.append(text)

    def code(self, binding):
        """Binds the lemma's variables to a place in the code: program.py's values. A rounded
        variable is bound to the code's rounded value; an exact one, a variable a lemma rounds,
        to the code's value whose exact result it stands for. A lemma used for several places is
        bound once for each."""
        self.instances.append(binding)

    def same_sign(self, truth, estimate):
        self.goal(f"(=> (> {estimate} 0.0) (> {truth} 0.0))")
        self.goal(f"(=> (< {estimate} 0.0) (< {truth} 0.0))")
        self.goal(f"(=> (= {estimate} 0.0) (= {truth} 0.0))")

    def text(self, negated=True):
        """The lemma in SMT-LIB: its facts and its conclusion's negation, or its facts alone."""
        lines = [PRELUDE]
        lines += [f"(declare-const {name} Real)" for name in self.reals]
        lines += [f"(assert {fact})" for fact in self.facts]
        if negated:
            lines.append(f"(assert (not (and {' '.join(self.goals)})))")
        lines.append("(check-sat)")
        return "\n".join(lines) + "\n"


class Floating:
    """A lemma in IEEE 754 binary64 itself, Z3's floating point theory: every value a double,
    every operation rounded to nearest with its infinities, NaNs and subnormals."""

    def __init__(self, name, *doubles):
        self.name = name
        self.doubles = list(doubles)
        self.facts = []
        self.goals = []
        self.operations = []
        self.instances = []

    def fact(self, text):
        self.facts.append(text)

    def goal(self, text):
        self.goals.append(text)

    def code(self, binding):
        """Binds the lemma's variables to a place in the code: program.py's values. A rounded
        variable is bound to the code's rounded value; an exact one, a variable a lemma rounds,
        to the code's value whose exact result it stands for. A lemma used for several places is
        bound once for each."""
        self.instances.append(binding)

    def text(self, negated=True):
        lines = ["(define-sort Double () (_ FloatingPoint 11 53))"]
        lines += [f"(declare-const {name} Double)" for name in self.doubles]
        lines += [f"(assert {fact})" for fact in self.facts]
        if negated:
            lines.append(f"(assert (not (and {' '.join(self.goals)})))")
        lines.append("(check-sat)")
        return "\n".join(lines) + "\n"


def parse(text):
    """An SMT-LIB expression as nested lists."""
    tokens = text.replace("(", " ( ").replace(")", " ) ").split()
    stack = [[]]
    for token in tokens:
        if token == "(":
            stack.append([])
        elif token == ")":
            done = stack.pop()
            stack[-1].append(done)
        else:
            stack[-1].append(token)
    return stack[0][0]


def evaluate(tree, binding, program):
    """An expression of a lemma as program.py's value, its variables bound."""
    if isinstance(tree, str):
        if tree in binding:
            return binding[tree]
        number = fractions.Fraction(tree)
        if fractions.Fraction(float(number)) != number:
            raise ValueError(f"{tree} is not a double")
        return program.constant(float(number))
    head, *rest = tree
    if head == "/":
        number = fractions.Fraction(rest[0][0] if isinstance(rest[0], list) else rest[0])
        number /= fractions.Fraction(rest[1])
        if fractions.Fraction(float(number)) != number:
            raise ValueError("a constant that is not a double")
        return program.constant(float(number))
    args = [evaluate(arg, binding, program) for arg in rest]
    if head == "absr":
        return program.magnitude(args[0])
    if head == "+" and len(args) == 2:
        return program.add(*args)
    if head == "-" and len(args) == 2:
        return program.sub(*args)
    if head == "-" and len(args) == 1:
        return program.neg(args[0])
    if head == "*" and len(args) == 2:
        return program.mul(*args)
    raise ValueError(f"no operation of the code is {head}")


def linked(lemma, program):
    """Whether every rounding the lemma states is the code's, at every place it is bound to:
    each rounded variable the code's operation on the code's operands, each roundoff of an
    error-free transformation the code's. Returns whether it is, and the first that is not."""
    if lemma.operations and not lemma.instances:
        return False, "rounds values bound to no code"
    for binding in lemma.instances:
        for operation in lemma.operations:
            kind, result = operation[0], operation[1]
            exact = operation[-1]
            if result not in binding:
                return False, f"{result} is bound to no code"
            try:
                tree = parse(exact) if exact.startswith("(") else exact
                value = evaluate(tree, binding, program)
            except (ValueError, KeyError) as error:
                return False, f"{result}: {error}"
            if binding[result] != value:
                return False, f"{result} is not the code's rounding of {exact}"
            if kind == "split":
                small = operation[2]
                if small not in binding:
                    return False, f"{small} is bound to no code"
                pair = {"mul": program.two_product, "add": program.two_sum,
                        "sub": program.two_diff}.get(value.op)
                if pair is None or binding[small] != pair(*value.args)[1]:
                    return False, f"{small} is not the code's roundoff of {result}"
    return True, ""


def program_module():
    """program.py, imported late: it imports this module's constants."""
    import program  # pylint: disable=import-outside-toplevel
    return program


def solve(solver, path, text):
    """Z3's verdict on an SMT-LIB file: sat, unsat, unknown, timeout or error."""
    path.write_text(text)
    result = subprocess.run([solver, "-T:600", str(path)], capture_output=True, text=True,
                            check=False)
    lines = result.stdout.strip().splitlines()
    return lines[0] if lines else "error"


def run(lemmas, checks=(), mutants=()):
    """Proves every lemma with the Z3 named by the first argument, z3 on the path by default.

    Each check is a function returning whether it held and what it found, run first. Each mutant
    is a lemma with one of the code's guards taken out or a bound shrunk below the proof's, which
    Z3 must refute with a counterexample: so the proofs need every guard they are given.

    Returns:
        The exit status: 1 when a check fails, a lemma is not proved or a mutant is.
    """
    solver = sys.argv[1] if len(sys.argv) > 1 else shutil.which("z3")
    if not solver:
        raise SystemExit("z3 not found")
    failed = False
    for check in checks:
        held, told = check()
        failed = failed or not held
        print(f"{check.__name__:22} {'held' if held else 'BROKEN':8} {told}", flush=True)
    with tempfile.TemporaryDirectory() as folder:
        for make in lemmas:
            lemma = make()
            held, told = linked(lemma, program_module())
            if not held:
                failed = True
                print(f"{lemma.name:22} UNLINKED {told}", flush=True)
            start = time.monotonic()
            proved = solve(solver, pathlib.Path(folder) / f"{make.__name__}.smt2", lemma.text())
            consistent = solve(solver, pathlib.Path(folder) / f"{make.__name__}.facts.smt2",
                               lemma.text(negated=False))
            verdict = "proved" if proved == "unsat" else proved.upper()
            if consistent != "sat":
                verdict = "VACUOUS"
            failed = failed or verdict != "proved"
            print(f"{lemma.name:22} {verdict:8} {time.monotonic() - start:6.1f} s", flush=True)
        for name, make in mutants:
            verdict = solve(solver, pathlib.Path(folder) / "mutant.smt2", make().text())
            refuted = verdict == "sat"
            failed = failed or not refuted
            print(f"{'mutant: ' + name:48} {'refuted' if refuted else verdict.upper()}",
                  flush=True)
    return 1 if failed else 0
