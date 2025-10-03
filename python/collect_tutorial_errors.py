#!/usr/bin/env python
"""
Collect and summarize error information from tutorial run directories.

Compatibility version (Python 2.7+ / 3.x). Enhanced to ignore 'Error' tokens
when they appear in table-like lines (Markdown or ASCII tables), unless
--include-table-errors is specified.

Usage examples:
  python scripts/collect_tutorial_errors.py --root .
  python scripts/collect_tutorial_errors.py --root . --json errors.json
  python scripts/collect_tutorial_errors.py --root . --include-logs
  python scripts/collect_tutorial_errors.py --root . --pattern "Segmentation fault"
  python scripts/collect_tutorial_errors.py --root . --include-table-errors   (to revert ignore behavior)

Exit codes:
  0 if no errors (or --always-zero)
  1 if errors detected (default)
"""

import argparse
import json
import os
import re
import sys

DEFAULT_FILE_NAMES = ["OUT.SCREEN"]
ERROR_TOKENS = [
    r"\bQBException\b",
    r"\bERROR\b",
    r"\bError\b",
    r"\bFAILED\b",
    r"\bException:\b",
]

TRACEBACK_START_RE = re.compile(r"^Traceback \(most recent call last\):")
QBEXCEPTION_RE = re.compile(r"QBException", re.IGNORECASE)

# Precompile table heuristics
MD_TABLE_ROW_RE = re.compile(r"^\s*\|.*\|\s*$")                     # | col | col |
MD_TABLE_SEP_RE = re.compile(r"^\s*\|?\s*:?-{3,}:?\s*(\|\s*:?-{3,}:?\s*)+\|?\s*$")  # |---|---|
ASCII_BORDER_RE = re.compile(r"^\s*\+[-+:=]{2,}\+\s*$")             # +-----+-----+
PIPE_HEAVY_LINE_RE = re.compile(r"^\s*(\|[^|]{0,120}){2,}\|?\s*$")  # multiple columns
BORDER_CHARS_RE = re.compile(r"^[\s\|\+\-:=]+$")                    # lines of just border chars

class ErrorBlock(object):
    def __init__(self, tutorial, file_path, error_type, summary,
                 lines, start_line, end_line):
        self.tutorial = tutorial
        self.file = file_path
        self.error_type = error_type
        self.summary = summary
        self.lines = lines
        self.start_line = start_line
        self.end_line = end_line


def color(s, code, enable):
    if not enable:
        return s
    return "\033[%sm%s\033[0m" % (code, s)


def iter_candidate_files(root, include_logs):
    for entry_name in sorted(os.listdir(root)):
        full = os.path.join(root, entry_name)
        if not os.path.isdir(full):
            continue
        tutorial = entry_name
        paths = []
        for fname in DEFAULT_FILE_NAMES:
            p = os.path.join(full, fname)
            if os.path.isfile(p):
                paths.append(p)
        if include_logs:
            for sub_name in os.listdir(full):
                if sub_name.endswith(".log"):
                    p = os.path.join(full, sub_name)
                    if os.path.isfile(p):
                        paths.append(p)
        for p in paths:
            yield tutorial, p


def capture_block(lines, start_index):
    block = []
    i = start_index
    n = len(lines)
    while i < n:
        line = lines[i]
        if i != start_index and line.strip() == "":
            break
        block.append(line)
        i += 1
    return block, i


def is_table_line(line):
    """Heuristic to detect if a line is part of a textual table, so we can
    suppress generic 'Error' token matches inside it."""
    stripped = line.rstrip("\n")
    if MD_TABLE_ROW_RE.match(stripped):
        return True
    if MD_TABLE_SEP_RE.match(stripped):
        return True
    if ASCII_BORDER_RE.match(stripped):
        return True
    # Lines that are just border chars
    if BORDER_CHARS_RE.match(stripped) and len(stripped.strip()) > 0 and ('-' in stripped or '|' in stripped or '+' in stripped):
        return True
    # Multi-column lines with many pipes
    if stripped.count('|') >= 2 and PIPE_HEAVY_LINE_RE.match(stripped):
        return True
    return False


def parse_file(tutorial, filepath, extra_patterns, generic_tokens, context,
               include_table_errors):
    try:
        with open(filepath, "r") as f:
            raw_lines = f.readlines()
    except Exception as e:
        return [ErrorBlock(
            tutorial, filepath, "FileReadError",
            "Could not read file: %s" % e, [], 0, 0
        )]

    results = []
    n = len(raw_lines)
    i = 0
    while i < n:
        line = raw_lines[i]

        # Traceback
        if TRACEBACK_START_RE.match(line):
            block, end = capture_block(raw_lines, i)
            results.append(ErrorBlock(
                tutorial, filepath, "Traceback",
                block[0].rstrip(), [l.rstrip("\n") for l in block],
                i + 1, end
            ))
            i = end
            continue

        # QBException
        if QBEXCEPTION_RE.search(line):
            block, end = capture_block(raw_lines, i)
            results.append(ErrorBlock(
                tutorial, filepath, "QBException",
                line.strip(), [l.rstrip("\n") for l in block],
                i + 1, end
            ))
            i = end
            continue

        # Extra patterns (with context)
        matched_extra = False
        for pat in extra_patterns:
            if pat.search(line):
                start_ctx = max(0, i - context)
                end_ctx = min(n, i + context + 1)
                block = raw_lines[start_ctx:end_ctx]
                results.append(ErrorBlock(
                    tutorial, filepath, "PatternMatch",
                    "Pattern %s matched: %s" % (pat.pattern, line.strip()),
                    [l.rstrip("\n") for l in block],
                    start_ctx + 1, end_ctx
                ))
                matched_extra = True
                break
        if matched_extra:
            i += 1
            continue

        # Generic tokens (skip if table line and user wants to ignore table errors)
        table_line = (not include_table_errors) and is_table_line(line)
        if not table_line:
            for tok in generic_tokens:
                if tok.search(line):
                    start_ctx = max(0, i - context)
                    end_ctx = min(n, i + context + 1)
                    block = raw_lines[start_ctx:end_ctx]
                    results.append(ErrorBlock(
                        tutorial, filepath, "GenericError",
                        line.strip(), [l.rstrip("\n") for l in block],
                        start_ctx + 1, end_ctx
                    ))
                    break

        i += 1

    return results


def build_regex_list(patterns):
    compiled = []
    for p in patterns:
        try:
            compiled.append(re.compile(p))
        except re.error as e:
            sys.stderr.write("Warning: invalid regex '%s': %s\n" % (p, e))
    return compiled


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Collect and summarize QB tutorial errors."
    )
    parser.add_argument("--root", default=".",
                        help="Root directory containing tutorial run subdirectories.")
    parser.add_argument("--include-logs", action="store_true",
                        help="Also scan *.log files in each tutorial directory.")
    parser.add_argument("--pattern", action="append", default=[],
                        help="Additional regex pattern to flag (repeatable).")
    parser.add_argument("--context", type=int, default=2,
                        help="Context lines before/after matches (default 2).")
    parser.add_argument("--json",
                        help="Write JSON summary to this file.")
    parser.add_argument("--no-color", action="store_true",
                        help="Disable colored output.")
    parser.add_argument("--always-zero", action="store_true",
                        help="Exit 0 even if errors found.")
    parser.add_argument("--include-table-errors", action="store_true",
                        help="Do NOT ignore lines that look like table rows.")
    args = parser.parse_args(argv)

    if not os.path.isdir(args.root):
        sys.stderr.write("ERROR: Root directory not found: %s\n" % args.root)
        return 1

    extra_patterns = build_regex_list(args.pattern)
    generic_tokens = [re.compile(t) for t in ERROR_TOKENS]

    tutorial_map = {}

    for tutorial, path in iter_candidate_files(args.root, args.include_logs):
        blocks = parse_file(
            tutorial, path,
            extra_patterns, generic_tokens,
            args.context, args.include_table_errors
        )
        if blocks:
            tutorial_map.setdefault(tutorial, []).extend(blocks)

    tutorials_with_errors = sorted(tutorial_map.keys())
    total_errors = sum(len(v) for v in tutorial_map.values())

    use_color = (not args.no_color) and sys.stdout.isatty()

    if not tutorials_with_errors:
        print(color("No errors detected in scanned tutorial outputs.", "32", use_color))
    else:
        print(color("Errors detected in %d tutorial(s), total error blocks: %d" %
                    (len(tutorials_with_errors), total_errors), "31", use_color))
        for tut in tutorials_with_errors:
            errs = tutorial_map[tut]
            print(color("\n=== %s (%d error block(s)) ===" %
                        (tut, len(errs)), "36", use_color))
            for idx, block in enumerate(errs, 1):
                label_color = "31" if block.error_type in ("QBException", "Traceback") else "33"
                header = "[%d] %s %s:%d-%d" % (idx, block.error_type, block.file,
                                               block.start_line, block.end_line)
                print(color(header, label_color, use_color))
                max_preview = 12
                for line in block.lines[:max_preview]:
                    print(color("    " + line, "2", use_color))
                if len(block.lines) > max_preview:
                    print(color("    ... (%d more lines)" %
                                (len(block.lines) - max_preview), "2", use_color))

    if args.json:
        out = []
        for tut in tutorials_with_errors:
            for block in tutorial_map[tut]:
                out.append({
                    "tutorial": block.tutorial,
                    "file": block.file,
                    "error_type": block.error_type,
                    "summary": block.summary,
                    "start_line": block.start_line,
                    "end_line": block.end_line,
                    "lines": block.lines,
                })
        try:
            with open(args.json, "w") as jf:
                json.dump({
                    "root": os.path.abspath(args.root),
                    "tutorial_count": len(tutorials_with_errors),
                    "error_block_count": total_errors,
                    "ignored_table_errors": (not args.include_table_errors),
                    "errors": out
                }, jf, indent=2)
            print(color("\nJSON report written to %s" % args.json, "32", use_color))
        except Exception as e:
            sys.stderr.write("Failed to write JSON file %s: %s\n" % (args.json, e))

    if tutorials_with_errors and not args.always_zero:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())