"""CLI/deployment regressions; all summary fixtures here are synthetic."""
from pathlib import Path
import copy
import json
import os
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'python'))
import summarize_MD_comparisons as summary
import openmm_xml_precision as precision
from test_openmm_xml_precision import legacy_summary, setup_summary

REPO = Path(__file__).resolve().parents[1]


def isolated_script(tmp_path, helper=None):
    app = tmp_path / 'app'
    app.mkdir()
    shutil.copy2(REPO/'python/summarize_MD_comparisons.py', app)
    if helper is True:
        shutil.copy2(REPO/'python/openmm_xml_precision.py', app)
    elif isinstance(helper, str):
        (app/'openmm_xml_precision.py').write_text(helper)
    return app/'summarize_MD_comparisons.py'


def run_cli(script, campaign, prefix, precision_enabled=True):
    cmd = [sys.executable, str(script), '--physics-only', '--root', str(campaign), '--prefix', str(prefix)]
    if precision_enabled:
        cmd.append('--recheck-xml-precision')
    env = os.environ.copy()
    env.pop('PYTHONPATH', None)
    return subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=30)


@pytest.mark.parametrize('helper, expected', [
    (None, 'required companion module is missing'),
    ('def broken(:\n', 'SyntaxError'),
    ('POLICY = {}\n', 'no callable audit_system_pair'),
    ('raise RuntimeError("broken helper import")\n', 'broken helper import'),
])
def test_setup_failure_is_single_error_and_preserves_reports(tmp_path, helper, expected):
    script = isolated_script(tmp_path, helper)
    campaign = tmp_path/'campaign'
    campaign.mkdir()
    case = campaign/'0000'
    case.mkdir()
    setup_summary(case, legacy_summary()).rename(case/'0000_summary.json')
    prefix = tmp_path/'report'
    sentinel = tmp_path/'report_summary.csv'
    sentinel.write_text('existing-results-must-not-be-overwritten\n')
    proc = run_cli(script, campaign, prefix)
    assert proc.returncode == 2
    assert expected in proc.stderr
    assert 'No structures were classified' in proc.stderr
    assert 'Structures:' not in proc.stdout
    assert sentinel.read_text() == 'existing-results-must-not-be-overwritten\n'
    assert not (tmp_path/'report_report.txt').exists()


def test_correct_helper_allows_regular_cli(tmp_path):
    script = isolated_script(tmp_path, True)
    case = tmp_path/'campaign'/'0000'
    case.mkdir(parents=True)
    setup_summary(case, legacy_summary()).rename(case/'0000_summary.json')
    proc = run_cli(script, case.parent, tmp_path/'report')
    assert proc.returncode == 0, proc.stderr
    assert 'REVIEW=0  GOOD=1' in proc.stdout
    assert 'MATCH_WITHIN_TOLERANCE=1' in proc.stdout


def test_no_option_needs_no_helper(tmp_path):
    script = isolated_script(tmp_path)
    case = tmp_path/'campaign'/'0000'
    case.mkdir(parents=True)
    setup_summary(case, legacy_summary()).rename(case/'0000_summary.json')
    proc = run_cli(script, case.parent, tmp_path/'report', False)
    assert proc.returncode == 0, proc.stderr
    assert 'REVIEW=1  GOOD=0' in proc.stdout


def test_unexpected_audit_exception_is_visible(tmp_path, monkeypatch, capsys):
    path = setup_summary(tmp_path, legacy_summary())
    def broken(*args):
        raise RuntimeError('synthetic-audit-error')
    monkeypatch.setattr(precision, 'audit_system_pair', broken)
    row = summary.analyze_one(path, physics_only=True, recheck_xml_precision=True)
    assert row['status'] == 'ERROR'
    assert 'RuntimeError' in row['review_reasons']
    summary.print_console_summary([row], {}, physics_only=True)
    out = capsys.readouterr().out
    assert 'Errors / failed runs:' in out
    assert 'synthetic-audit-error' in out


def test_malformed_summary_error_visible_without_precision(tmp_path, capsys):
    path = tmp_path/'0000_summary.json'
    path.write_text('{broken')
    row = summary.analyze_one(path, physics_only=True)
    assert row['status'] == 'ERROR'
    summary.print_console_summary([row], {}, physics_only=True)
    assert 'cannot parse JSON' in capsys.readouterr().out
