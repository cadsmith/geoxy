#!/usr/bin/env python3
"""Shows what an editor running clangd shows: every diagnostic of every file, and hovers.

    python3 -I tools/format/editor.py [--hover FILE:LINE:COLUMN] [file ...]

It talks to clangd over the language server protocol, as VS Code's clangd extension does,
so it reports what `clangd --check` leaves out: the include cleaner's missing and unused
includes, and every hint. With no file named, every tracked .hpp and .cpp is opened. The
exit status is 1 when any diagnostic is reported. With --hover, the hover text at a
1-based place is printed as the editor would show it.

CLANGD names the clangd to run, LLVM's own by default.
"""

import json
import os
import pathlib
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
CLANGD_HINTS = ("/opt/homebrew/opt/llvm/bin/clangd", "/usr/local/opt/llvm/bin/clangd")


def clangd():
    named = os.environ.get("CLANGD")
    if named:
        return named
    for hint in CLANGD_HINTS:
        if pathlib.Path(hint).exists():
            return hint
    found = shutil.which("clangd")
    if not found:
        raise SystemExit("clangd not found")
    return found


class Server:
    def __init__(self):
        self.process = subprocess.Popen([clangd(), "--log=error", "--background-index=false"],
                                        cwd=ROOT, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                        stderr=subprocess.DEVNULL)
        self.next_id = 0
        self.diagnostics = {}
        self.request("initialize", {"processId": os.getpid(), "rootUri": ROOT.as_uri(),
                                    "capabilities": {"textDocument": {"hover": {
                                        "contentFormat": ["markdown", "plaintext"]},
                                        "publishDiagnostics": {}}}})
        self.notify("initialized", {})

    def send(self, message):
        body = json.dumps(message).encode()
        self.process.stdin.write(b"Content-Length: %d\r\n\r\n" % len(body) + body)
        self.process.stdin.flush()

    def receive(self):
        length = 0
        while True:
            line = self.process.stdout.readline()
            if not line:
                raise SystemExit("clangd stopped")
            if line in (b"\r\n", b"\n"):
                break
            if line.lower().startswith(b"content-length:"):
                length = int(line.split(b":")[1])
        message = json.loads(self.process.stdout.read(length))
        if message.get("method") == "textDocument/publishDiagnostics":
            params = message["params"]
            self.diagnostics[params["uri"]] = params["diagnostics"]
        return message

    def notify(self, method, params):
        self.send({"jsonrpc": "2.0", "method": method, "params": params})

    def request(self, method, params):
        self.next_id += 1
        ident = self.next_id
        self.send({"jsonrpc": "2.0", "id": ident, "method": method, "params": params})
        while True:
            message = self.receive()
            if message.get("id") == ident and "method" not in message:
                return message.get("result")

    def open(self, path):
        self.notify("textDocument/didOpen", {"textDocument": {
            "uri": path.as_uri(), "languageId": "cpp", "version": 1, "text": path.read_text()}})
        # A request that needs the file built returns only once clangd has built it, and the
        # diagnostics are published before it.
        self.request("textDocument/documentSymbol", {"textDocument": {"uri": path.as_uri()}})
        while path.as_uri() not in self.diagnostics:
            self.receive()
        return self.diagnostics[path.as_uri()]

    def hover(self, path, line, column):
        result = self.request("textDocument/hover", {"textDocument": {"uri": path.as_uri()},
                                                     "position": {"line": line - 1,
                                                                  "character": column - 1}})
        return result["contents"]["value"] if result else ""

    def close(self):
        self.request("shutdown", None)
        self.notify("exit", None)
        self.process.wait(timeout=10)



# Built apart, against geoxy installed: no build folder of ours knows how it compiles.
APART = ":(exclude)tests/package"


def sources(named):
    if named:
        return [pathlib.Path(name).resolve() for name in named]
    listed = subprocess.run(["git", "ls-files", "--cached", "--others", "--exclude-standard",
                             "*.hpp", "*.cpp", APART], cwd=ROOT, check=True,
                            capture_output=True, text=True).stdout.split()
    return [ROOT / name for name in listed]


def main():
    arguments = sys.argv[1:]
    hover = None
    if "--hover" in arguments:
        at = arguments.index("--hover")
        hover = arguments[at + 1]
        del arguments[at:at + 2]
    server = Server()
    try:
        if hover:
            name, line, column = hover.rsplit(":", 2)
            path = pathlib.Path(name).resolve()
            server.open(path)
            print(server.hover(path, int(line), int(column)))
            return
        count = 0
        for path in sources(arguments):
            for diagnostic in server.open(path):
                start = diagnostic["range"]["start"]
                code = diagnostic.get("code", diagnostic.get("source", ""))
                print(f"{path.relative_to(ROOT)}:{start['line'] + 1}:{start['character'] + 1}: "
                      f"[{code}] {diagnostic['message'].splitlines()[0]}")
                count += 1
        sys.exit(1 if count else 0)
    finally:
        server.close()


if __name__ == "__main__":
    main()
