#!/usr/bin/env python3
"""
Collect and summarize error information from tutorial run directories.

Features:
- Recursively scans a root directory (default: current working dir).
- For each subdirectory, inspects:
    * OUT.SCREEN (default primary file)
    * Any *.log files (optional, enabled with --include-logs)
- Extracts:
    * QBException blocks (line containing 'QBException' + following context until blank line)
    * Python traceback blocks starting at 'Traceback (most recent call last):'
    * Lines containing generic error tokens (ERROR, FAILED, Exception:)
- Provides configurable regex patterns via --pattern (multiple allowed)
- Outputs a colored human-readable summary plus (optionally) a JSON file (--json)
- Exit code:
    * 0 if no errors found
    * 1 if errors found (override with --always-zero)

Usage:
    python scripts/collect_tutorial_errors.py --root ./ --json errors.json
"""

from __future__ import annotations
import argparse
import json
import os
import re
import sys
from dataclasses import dataclass, asdict
from typing import List, Optional, Dict, Iterable, Tuple

# ----------------------------- Configuration ---------------------------------

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

# ----------------------------- Data Classes ----------------------------------

@dataclass
class ErrorBlock:
    tutorial: str
    file: str
    error_type: str          # e.g. QBException, Traceback, GenericError
    summary: str             # First line / headline
    lines: List[str]         # Full captured lines
    start_line: int          # 1-based line number
    end_line: int            # inclusive

@dataclass
class TutorialResult:
    tutorial: str
    errors: List[ErrorBlock]

# ----------------------------- Utility Functions -----------------------------

def color(s: str, code: str, enable: bool) -> str:
    if not enable:
        return s
    return f"\033[{code}m{s}\033[0m"

def iter_candidate_files(root: str,
                         include_logs: bool) -> Iterable[Tuple[str, str]]:
    """
    Yields (tutorial_name, file_path).
    A tutorial directory is any directory containing at least one target file.
    """
    for entry in os.scandir(root):
        if not entry.is_dir():
            continue
        tutorial = entry.name
        paths = []
        for fname in DEFAULT_FILE_NAMES:
            p = os.path.join(entry.path, fname)
            if os.path.isfile(p):
                paths.append(p)
        if include_logs:
            for sub in os.scandir(entry.path):
                if sub.is_file() and sub.name.endswith(".log"):
                    paths.append(sub.path)
        for p in paths:
            yield tutorial, p

def capture_block(lines: List[str], start_index: int) -> Tuple[List[str], int]:
    """
    Capture lines from start_index forward until a blank line or end of file.
    Returns (captured_lines, end_index_exclusive).
    """
    block = []
    i = start_index
    while i < len(lines):
        line = lines[i]
        if i != start_index and line.strip() == "":
            break
        block.append(line)
        i += 1
    return block, i

def parse_file(tutorial: str, filepath: str,
               extra_patterns: List[re.Pattern],
               generic_tokens: List[re.Pattern],
               context: int) -> List[ErrorBlock]:
    """
    Parse a log file for:
      - QBException blocks
      - Python tracebacks
      - Additional user-specified patterns
      - Generic error lines
    """
    try:
        with open(filepath, "r", errors="replace") as f:
            raw_lines = f.readlines()
    except (OSError, UnicodeDecodeError) as e:
        return [ErrorBlock(
            tutorial=tutorial,
            file=filepath,
            error_type="FileReadError",
            summary=f"Could not read file: {e}",
            lines=[],
            start_line=0,
            end_line=0
        )]

    results: List[ErrorBlock] = []
    n = len(raw_lines)
    i = 0
    while i < n:
        line = raw_lines[i]
        # Traceback
        if TRACEBACK_START_RE.match(line):
            block, end = capture_block(raw_lines, i)
            # Optionally extend context lines after blank? Not needed.
            results.append(ErrorBlock(
                tutorial=tutorial,
                file=filepath,
                error_type="Traceback",
                summary=block[0].rstrip(),
                lines=[l.rstrip("\n") for l in block],
                start_line=i + 1,
                end_line=end
            ))
            i = end
            continue
        # QBException
        if QBEXCEPTION_RE.search(line):
            block, end = capture_block(raw_lines, i)
            results.append(ErrorBlock(
                tutorial=tutorial,
                file=filepath,
                error_type="QBException",
                summary=line.strip(),
                lines=[l.rstrip("\n") for l in block],
                start_line=i + 1,
                end_line=end
            ))
            i = end
            continue
        # Extra patterns (capture with context)
        matched_extra = False
        for pat in extra_patterns:
            if pat.search(line):
                start_ctx = max(0, i - context)
                end_ctx = min(n, i + context + 1)
                block = raw_lines[start_ctx:end_ctx]
                results.append(ErrorBlock(
                    tutorial=tutorial,
                    file=filepath,
                    error_type="PatternMatch",
                    summary=f"Pattern {pat.pattern} matched: {line.strip()}",
                    lines=[l.rstrip("\n") for l in block],
                    start_line=start_ctx + 1,
                    end_line=end_ctx
                ))
                matched_extra = True
                break
        if matched_extra:
            i += 1
            continue
        # Generic tokens (single line)
        for tok in generic_tokens:
            if tok.search(line):
                start_ctx = max(0, i - context)
                end_ctx = min(n, i + context + 1)
                block = raw_lines[start_ctx:end_ctx]
                results.append(ErrorBlock(
                    tutorial=tutorial,
                    file=filepath,
                    error_type="GenericError",
                    summary=line.strip(),
                    lines=[l.rstrip("\n") for l in block],
                    start_line=start_ctx + 1,
                    end_line=end_ctx
                ))
                break
        i += 1

    return results

# ----------------------------- Main Logic -------------------------------------

def build_regex_list(patterns: List[str]) -> List[re.Pattern]:
    compiled = []
    for p in patterns:
        try:
            compiled.append(re.compile(p))
        except re.error as e:
            print(f"Warning: invalid regex '{p}': {e}", file=sys.stderr)
    return compiled

def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Collect and summarize QB tutorial errors from OUT.SCREEN and log files."
    )
    parser.add_argument("--root", default=".", help="Root directory containing tutorial run subdirectories.")
    parser.add_argument("--include-logs", action="store_true", help="Also scan *.log files in each tutorial directory.")
    parser.add_argument("--pattern", action="append", default=[], help="Additional regex pattern to flag (can pass multiple).")
    parser.add_argument("--context", type=int, default=2, help="Context lines before/after pattern or generic errors (default: 2).")
    parser.add_argument("--json", help="Write JSON summary to this file.")
    parser.add_argument("--no-color", action="store_true", help="Disable colored output.")
    parser.add_argument("--always-zero", action="store_true", help="Exit 0 even if errors found.")
    args = parser.parse_args(argv)

    if not os.path.isdir(args.root):
        print(f"ERROR: Root directory not found: {args.root}", file=sys.stderr)
        return 1

    extra_patterns = build_regex_list(args.pattern)
    generic_tokens = [re.compile(t) for t in ERROR_TOKENS]

    tutorial_map: Dict[str, List[ErrorBlock]] = {}

    for tutorial, path in iter_candidate_files(args.root, args.include_logs):
        blocks = parse_file(tutorial, path, extra_patterns, generic_tokens, args.context)
        if blocks:
            tutorial_map.setdefault(tutorial, []).extend(blocks)

    # Summarize
    tutorials_with_errors = sorted(tutorial_map.keys())
    total_errors = sum(len(v) for v in tutorial_map.values())

    use_color = not args.no_color and sys.stdout.isatty()

    if not tutorials_with_errors:
        print(color("No errors detected in scanned tutorial outputs.", "32", use_color))
    else:
        print(color(f"Errors detected in {len(tutorials_with_errors)} tutorial(s), total error blocks: {total_errors}", "31", use_color))
        for tut in tutorials_with_errors:
            errs = tutorial_map[tut]
            print(color(f"\n=== {tut} ({len(errs)} error block(s)) ===", "36", use_color))
            for idx, block in enumerate(errs, 1):
                label_color = "31" if block.error_type in ("QBException", "Traceback") else "33"
                header = f"[{idx}] {block.error_type} {block.file}:{block.start_line}-{block.end_line}"
                print(color(header, label_color, use_color))
                # Show first few lines (limit)
                max_preview = 12
                for j, line in enumerate(block.lines[:max_preview], 1):
                    print(color(f"    {line}", "2", use_color))
                if len(block.lines) > max_preview:
                    print(color(f"    ... ({len(block.lines) - max_preview} more lines)", "2", use_color))

    # JSON output
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
                    "errors": out
                }, jf, indent=2)
            print(color(f"\nJSON report written to {args.json}", "32", use_color))
        except OSError as e:
            print(f"Failed to write JSON file {args.json}: {e}", file=sys.stderr)

    if tutorials_with_errors and not args.always_zero:
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())