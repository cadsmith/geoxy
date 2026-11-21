# Security Policy

## Reporting a vulnerability

Report it through GitHub's private vulnerability reporting:
<https://github.com/cadsmith/geoxy/security/advisories/new>. The report is
visible to the maintainers alone until an advisory is published. Do not open
an issue or a pull request for it; both are public.

Include the curves and the call that trigger the problem, their numbers
written so that they read back to the bit (`%a`, or all 17 digits); the
geoxy release, compiler and platform; and the output as it was printed: the
sanitizer report, the stack trace or the fuzzer's artifact.

You get a first response within 14 days. Please allow 90 days from the
report before disclosing it publicly.

## What counts

A security vulnerability is a crash, a read or write outside bounds, a read
of uninitialised memory, undefined behaviour, a hang or unbounded memory use
that curves, points or numbers of any value cause through a public function,
numbers that are not finite among them.

It is not a vulnerability, and belongs in a public issue, when the calling
program misuses the API, such as reading a curve at an index a `Tree` no
longer holds; when an answer is inexact or wrong without any of the above;
or when it affects the documentation site or the build and test tooling.

## Supported versions

geoxy is before 1.0. Fixes land on `main` and ship in the next release, and
only the latest release receives them.

## Disclosure and credit

Once a fix is released, the advisory is published and the release notes name
the vulnerability. Both credit the reporter unless they ask not to be named.
The fix comes with a test that pins the input.

## Reports from tools

A report from a static analyser, a fuzzer or an AI tool goes through the
same private channel, and is most useful with an input that reproduces the
problem. Triage asks for one when it is missing, and a report that nothing
reproduces is closed, still in private.
