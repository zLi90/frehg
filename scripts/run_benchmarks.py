#!/usr/bin/env python3
"""Run and validate Frehg2 benchmarks with one repeatable command.

P18 intentionally separates execution from comparison:

- b1/b2 are legacy-regression cases with numeric Frehg references.
- b3-b6 are converted SERGHEI-style fixtures. They are included in the
  validation report, but are not release-passing numerical comparisons until
  P18/P19 define comparable variables and runtime support for each case.
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


DEFAULT_CASES = [
    ("b1-sw", "benchmarks/b1-sw/b1-sw.production.yaml"),
    ("b2-gw", "benchmarks/b2-gw/b2-gw.production.yaml"),
    ("b3-kirkland", "benchmarks/b3-kirkland/b3-kirkland.yaml"),
    ("b4-govindaraju", "benchmarks/b4-govindaraju/b4-govindaraju.yaml"),
    ("b5-vcatchment", "benchmarks/b5-vcatchment/b5-vcatchment.yaml"),
    ("b6-superslab", "benchmarks/b6-superslab/b6-superslab.yaml"),
]


@dataclass
class VariableResult:
    name: str
    status: str
    compared_files: int = 0
    failed_files: int = 0
    max_relative_l2: float | None = None
    max_abs: float | None = None
    message: str = ""


@dataclass
class CaseResult:
    case_id: str
    config: str
    status: str
    executed: bool
    output_dir: str | None
    reference_type: str
    variables: list[VariableResult] = field(default_factory=list)
    messages: list[str] = field(default_factory=list)
    returncode: int = 0


def strip_comment(line: str) -> str:
    in_quote = False
    quote = ""
    for index, char in enumerate(line):
        if char in {"'", '"'} and (index == 0 or line[index - 1] != "\\"):
            if not in_quote:
                in_quote = True
                quote = char
            elif quote == char:
                in_quote = False
        if char == "#" and not in_quote:
            return line[:index]
    return line


def parse_scalar(value: str) -> Any:
    value = value.strip()
    if value in {"null", "~"}:
        return None
    if value in {"true", "false"}:
        return value == "true"
    if value.startswith('"') and value.endswith('"'):
        return value[1:-1]
    if value.startswith("'") and value.endswith("'"):
        return value[1:-1]
    if value.startswith("[") and value.endswith("]"):
        inner = value[1:-1].strip()
        if not inner:
            return []
        return [parse_scalar(item.strip()) for item in inner.split(",")]
    try:
        if any(char in value for char in [".", "e", "E"]):
            return float(value)
        return int(value)
    except ValueError:
        return value


def parse_simple_yaml(path: Path) -> dict[str, Any]:
    """Parse the YAML subset used by benchmark configs."""
    root: dict[str, Any] = {}
    stack: list[tuple[int, Any]] = [(-1, root)]
    pending_list_key: tuple[int, dict[str, Any], str] | None = None

    for raw in path.read_text().splitlines():
        line = strip_comment(raw).rstrip()
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip(" "))
        text = line.strip()

        while stack and stack[-1][0] >= indent:
            stack.pop()
        parent = stack[-1][1]

        if text.startswith("- "):
            item_text = text[2:].strip()
            if not isinstance(parent, list):
                if pending_list_key is None:
                    raise ValueError(f"list item without list parent in {path}: {text}")
                _, owner, key = pending_list_key
                owner[key] = []
                parent = owner[key]
                stack.append((indent - 2, parent))
            if item_text and ":" not in item_text:
                parent.append(parse_scalar(item_text))
                continue
            item: dict[str, Any] = {}
            parent.append(item)
            if item_text:
                key, value = item_text.split(":", 1)
                item[key.strip()] = parse_scalar(value) if value.strip() else {}
            stack.append((indent, item))
            continue

        if ":" not in text:
            continue
        key, value = text.split(":", 1)
        key = key.strip()
        value = value.strip()
        if value:
            parent[key] = parse_scalar(value)
            pending_list_key = None
        else:
            parent[key] = {}
            pending_list_key = (indent, parent, key)
            stack.append((indent, parent[key]))

    return root


def get_path(data: dict[str, Any], path: str, default: Any = None) -> Any:
    value: Any = data
    for part in path.split("."):
        if not isinstance(value, dict) or part not in value:
            return default
        value = value[part]
    return value


def read_numbers(path: Path) -> list[float]:
    values: list[float] = []
    with path.open() as handle:
        for line in handle:
            for item in line.split():
                try:
                    values.append(float(item))
                except ValueError:
                    continue
    return values


def metrics(reference: list[float], output: list[float]) -> tuple[float, float]:
    count = min(len(reference), len(output))
    if count == 0:
        return math.inf, math.inf
    ref = reference[:count]
    new = output[:count]
    diff = [left - right for left, right in zip(ref, new)]
    norm = math.sqrt(sum(value * value for value in ref))
    rel_l2 = math.sqrt(sum(value * value for value in diff)) / norm if norm > 0.0 else 0.0
    max_abs = max(abs(value) for value in diff)
    return rel_l2, max_abs


def resolve_config_path(repo_root: Path, case_id: str) -> Path:
    for known_id, relative in DEFAULT_CASES:
        if known_id == case_id:
            return repo_root / relative
    candidate = repo_root / case_id
    if candidate.exists():
        return candidate
    raise FileNotFoundError(f"unknown benchmark case: {case_id}")


def executable_path(repo_root: Path, override: str | None) -> Path:
    if override:
        return Path(override)
    candidates = [repo_root / "build/src/frehg2", repo_root / "build/frehg2"]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def output_dir_for(config_path: Path, data: dict[str, Any]) -> Path:
    filename = Path(str(get_path(data, "output.filename", "out/output.h5")))
    if not filename.is_absolute():
        filename = config_path.parent / filename
    return filename.parent


def reference_dir_for(config_path: Path, validation: dict[str, Any]) -> Path:
    raw = Path(str(get_path(validation, "reference.path", "")))
    if raw.is_absolute():
        return raw
    return (config_path.parent / raw).resolve()


def should_execute(case_id: str, reference_type: str, execute_converted: bool) -> bool:
    if reference_type == "frehg_legacy":
        return True
    return execute_converted and case_id.startswith(("b3", "b4", "b5", "b6"))


def compare_variables(
    output_dir: Path,
    reference_dir: Path,
    variables: list[dict[str, Any]],
    tolerance: float | None,
    max_abs_tolerance: float | None,
) -> list[VariableResult]:
    results: list[VariableResult] = []
    for variable in variables:
        name = str(variable.get("name", "unnamed"))
        output_pattern = str(variable.get("output", ""))
        reference_pattern = str(variable.get("reference", output_pattern))
        reference_files = sorted(path for path in reference_dir.iterdir()
                                 if fnmatch.fnmatch(path.name, reference_pattern))
        if not reference_files:
            results.append(VariableResult(name=name, status="skip", message="no reference files"))
            continue

        failed = 0
        max_rel = 0.0
        max_abs = 0.0
        compared = 0
        for reference_file in reference_files:
            output_name = reference_file.name
            if "*" in output_pattern:
                prefix = output_pattern.split("*", 1)[0]
                ref_prefix = reference_pattern.split("*", 1)[0]
                output_name = prefix + reference_file.name[len(ref_prefix):]
            output_file = output_dir / output_name
            if not output_file.exists():
                failed += 1
                continue
            rel_l2, abs_err = metrics(read_numbers(reference_file), read_numbers(output_file))
            max_rel = max(max_rel, rel_l2)
            max_abs = max(max_abs, abs_err)
            compared += 1
            if tolerance is not None and rel_l2 > tolerance:
                failed += 1
            if max_abs_tolerance is not None and abs_err > max_abs_tolerance:
                failed += 1

        status = "pass" if failed == 0 and compared > 0 else "fail"
        results.append(
            VariableResult(
                name=name,
                status=status,
                compared_files=compared,
                failed_files=failed,
                max_relative_l2=max_rel,
                max_abs=max_abs,
            )
        )
    return results


def run_case(
    repo_root: Path,
    case_id: str,
    config_path: Path,
    exe: Path,
    execute_converted: bool,
    fresh: bool,
) -> CaseResult:
    data = parse_simple_yaml(config_path)
    validation = get_path(data, "validation", {}) or {}
    reference_type = str(get_path(validation, "reference.type", "none"))
    variables = get_path(validation, "variables", []) or []
    tolerance = get_path(validation, "tolerances.relative_l2", None)
    max_abs_tolerance = get_path(validation, "tolerances.max_abs", None)
    output_dir = output_dir_for(config_path, data)
    reference_dir = reference_dir_for(config_path, validation)
    execute = should_execute(case_id, reference_type, execute_converted)

    result = CaseResult(
        case_id=case_id,
        config=str(config_path.relative_to(repo_root)),
        status="skip",
        executed=False,
        output_dir=str(output_dir.relative_to(repo_root)) if output_dir.is_relative_to(repo_root) else str(output_dir),
        reference_type=reference_type,
    )

    if not execute:
        result.messages.append(
            "converted SERGHEI fixture recorded; use --execute-converted to run experimental driver path"
        )
        result.status = "skip"
        return result

    if fresh and output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    completed = subprocess.run([str(exe), str(config_path)], cwd=repo_root, check=False)
    result.executed = True
    result.returncode = completed.returncode
    if completed.returncode != 0:
        result.status = "fail"
        result.messages.append(f"execution failed with return code {completed.returncode}")
        return result

    if not variables:
        result.status = "pass"
        result.messages.append("execution completed; no validation variables configured")
        return result

    result.variables = compare_variables(
        output_dir, reference_dir, variables, tolerance, max_abs_tolerance)
    result.status = "pass" if all(variable.status in {"pass", "skip"} for variable in result.variables) else "fail"
    return result


def result_to_dict(result: CaseResult) -> dict[str, Any]:
    return {
        "case_id": result.case_id,
        "config": result.config,
        "status": result.status,
        "executed": result.executed,
        "output_dir": result.output_dir,
        "reference_type": result.reference_type,
        "returncode": result.returncode,
        "messages": result.messages,
        "variables": [variable.__dict__ for variable in result.variables],
    }


def write_reports(report_dir: Path, results: list[CaseResult]) -> None:
    report_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "overall_status": "pass" if all(result.status != "fail" for result in results) else "fail",
        "cases": [result_to_dict(result) for result in results],
    }
    (report_dir / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    lines = ["# Frehg2 Validation Summary", ""]
    lines.append(f"Overall status: `{summary['overall_status']}`")
    lines.append("")
    lines.append("| Case | Status | Executed | Reference | Notes |")
    lines.append("| --- | --- | --- | --- | --- |")
    for result in results:
        notes = "; ".join(result.messages)
        if result.variables:
            variable_notes = [
                f"{variable.name}:{variable.status}"
                for variable in result.variables
            ]
            notes = "; ".join(filter(None, [notes, ", ".join(variable_notes)]))
        lines.append(
            f"| `{result.case_id}` | `{result.status}` | `{result.executed}` | "
            f"`{result.reference_type}` | {notes} |"
        )
    lines.append("")
    (report_dir / "validation_summary.md").write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--executable", help="Path to frehg2 executable")
    parser.add_argument("--cases", nargs="*", default=[case_id for case_id, _ in DEFAULT_CASES])
    parser.add_argument("--report-dir", type=Path, default=Path("build/validation"))
    parser.add_argument("--execute-converted", action="store_true",
                        help="run b3-b6 through current experimental driver paths")
    parser.add_argument("--no-fresh", action="store_true", help="do not remove existing output dirs")
    parser.add_argument("--strict-skips", action="store_true",
                        help="treat skipped converted cases as failures")
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    exe = executable_path(repo_root, args.executable)
    results = [
        run_case(
            repo_root,
            case_id,
            resolve_config_path(repo_root, case_id),
            exe,
            args.execute_converted,
            not args.no_fresh,
        )
        for case_id in args.cases
    ]
    report_dir = args.report_dir if args.report_dir.is_absolute() else repo_root / args.report_dir
    write_reports(report_dir, results)

    has_failure = any(result.status == "fail" for result in results)
    has_skip = any(result.status == "skip" for result in results)
    if has_failure or (args.strict_skips and has_skip):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
