#!/usr/bin/env python3
"""Read-only, fail-closed comparison of standard OpenMM System XML files.

This is a numerical-equivalence audit, not bitwise equality and not an MD test.
The narrow LJ-exception tolerances are calibrated to the supplied 1fd0 AMBER
E16.8 reconstruction case. They are explicit policy, not OpenMM defaults.
Unsupported forces, schemas, virtual sites, or nonempty parameter offsets never
produce MATCH. No System, State, topology, or original summary is modified.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
from typing import Any

POLICY = {
    "id": "amber-exception-reconstruction-v1",
    "other_numeric": {"atol": 1e-12, "rtol": 1e-12},
    "exception_charge_product_e2": {"atol": 1e-12, "rtol": 1e-12},
    "exception_sigma_nm": {"atol": 1e-12, "rtol": 2e-9},
    "exception_epsilon_kj_mol": {"atol": 1e-12, "rtol": 2e-8},
    "legacy_diagnostic": {"atol": 1e-10, "rtol": 1e-10, "relative_to": "input2"},
    "activity_rule": "Exact zero/nonzero and charge-sign transitions always fail; no activity epsilon.",
}
FORCE_SCHEMAS = {
    "HarmonicBondForce": ("2", "Bonds", "Bond", ("p1", "p2"), ("d", "k")),
    "HarmonicAngleForce": ("2", "Angles", "Angle", ("p1", "p2", "p3"), ("a", "k")),
    "PeriodicTorsionForce": ("2", "Torsions", "Torsion", ("p1", "p2", "p3", "p4", "periodicity"), ("phase", "k")),
}
NB_ATTRS = {
    "alpha", "cutoff", "dispersionCorrection", "ewaldTolerance", "exceptionsUsePeriodic",
    "forceGroup", "includeDirectSpace", "ljAlpha", "ljnx", "ljny", "ljnz", "method",
    "name", "nx", "ny", "nz", "recipForceGroup", "rfDielectric", "switchingDistance",
    "type", "useSwitchingFunction", "version",
}
NB_FLOAT_ATTRS = {"alpha", "cutoff", "ewaldTolerance", "ljAlpha", "rfDielectric", "switchingDistance"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def number(value: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("non-finite XML parameter")
    return result


def close(a: float, b: float, rule: dict) -> bool:
    return abs(a-b) <= rule["atol"] + rule["rtol"] * max(abs(a), abs(b))


def changed_zero(a: float, b: float) -> bool:
    return (a == 0) != (b == 0)


def legacy_different(a: float, b: float) -> bool:
    rule = POLICY["legacy_diagnostic"]
    return abs(a-b) > rule["atol"] + rule["rtol"] * abs(b)


def parse_xml(path: Path) -> ET.Element:
    raw = path.read_bytes()
    if b"<!DOCTYPE" in raw.upper() or b"<!ENTITY" in raw.upper():
        raise ValueError("DTD/entity declarations are not supported")
    return ET.fromstring(raw)


class Audit:
    def __init__(self) -> None:
        self.problems: list[str] = []
        self.unsupported: list[str] = []
        self.sections: dict[str, Any] = {}
        self.exceptions: list[dict[str, Any]] = []
        self.particles = 0

    def schema(self, node: ET.Element, attrs: set, children: list[str], label: str) -> None:
        # Explicit schemas prevent unknown fields from being silently skipped.
        if set(node.attrib) != attrs or Counter(c.tag for c in node) != Counter(children):
            self.unsupported.append(f"{label}: unsupported attributes or children")

    def attrs(self, a: ET.Element, b: ET.Element, label: str, numeric: set = frozenset(), ignore: set = frozenset()) -> None:
        for k in sorted((set(a.attrib) | set(b.attrib)) - set(ignore)):
            x, y = a.get(k), b.get(k)
            if x is None or y is None:
                self.problems.append(f"{label}: attribute {k} missing from one input")
            elif k in numeric:
                if not close(number(x), number(y), POLICY["other_numeric"]):
                    self.problems.append(f"{label}: {k} differs")
            elif x != y:
                self.problems.append(f"{label}: {k} differs ({x} vs {y})")

    def table(self, a: ET.Element, b: ET.Element, label: str, tag: str, ids: tuple, vals: tuple, unordered_pair: bool = False, guard: tuple = ()) -> None:
        maps = []
        for container in (a, b):
            if container.attrib:
                self.unsupported.append(f"{label}: unexpected table attributes")
            records = defaultdict(list)
            for n, row in enumerate(container):
                if row.tag != tag or set(row.attrib) != set(ids+vals) or len(row):
                    self.unsupported.append(f"{label}: unsupported row schema")
                    continue
                key = tuple(int(row.get(i)) for i in ids) if ids else (n,)
                if unordered_pair:
                    key = tuple(sorted(key))
                for field in ids:
                    if field.startswith('p') and field != 'periodicity':
                        index = int(row.get(field))
                        if not 0 <= index < self.particles:
                            raise ValueError(f"{label}: particle index out of range")
                records[key].append(tuple(number(row.get(v)) for v in vals))
            maps.append(records)
        x, y = maps
        only1, only2 = set(x)-set(y), set(y)-set(x)
        multiplicity = sum(len(x[k]) != len(y[k]) for k in set(x)&set(y))
        count = 0
        max_delta = dict.fromkeys(vals, 0.0)
        for k in set(x)&set(y):
            for va, vb in zip(sorted(x[k]), sorted(y[k])):
                bad = False
                for field, aa, bb in zip(vals, va, vb):
                    max_delta[field] = max(max_delta[field], abs(aa-bb))
                    if not close(aa, bb, POLICY['other_numeric']) or (field in guard and (changed_zero(aa, bb) or ((aa < 0 < bb) or (bb < 0 < aa)))):
                        bad = True
                count += bad
        self.sections[label] = dict(count1=len(a), count2=len(b), only1=len(only1), only2=len(only2), multiplicity_mismatches=multiplicity, parameter_mismatches=count, max_absolute_deltas=max_delta)
        if only1 or only2 or multiplicity or count:
            self.problems.append(f"{label}: inventory/multiplicity/parameter mismatch")

    def nonbonded(self, a: ET.Element, b: ET.Element, label: str) -> None:
        children = ['GlobalParameters', 'ParticleOffsets', 'ExceptionOffsets', 'Particles', 'Exceptions']
        for f in (a, b):
            self.schema(f, NB_ATTRS, children, label)
            if f.get('version') != '4':
                self.unsupported.append(f"{label}: unsupported serializer version")
            for name in children[:3]:
                c = f.find(name)
                if c is None or len(c) or c.attrib:
                    self.unsupported.append(f"{label}: nonempty or missing {name}; cannot certify effective parameters")
        self.attrs(a, b, label, NB_FLOAT_ATTRS, {'name'})
        self.table(a.find('Particles'), b.find('Particles'), label+'.particles', 'Particle', (), ('q', 'sig', 'eps'), guard=('q', 'eps'))
        maps = []
        for f in (a, b):
            records = {}
            c = f.find('Exceptions')
            if c is None:
                raise ValueError('missing Exceptions')
            if c.attrib:
                self.unsupported.append(label+': unexpected Exceptions attributes')
            for e in c:
                if e.tag != 'Exception' or set(e.attrib) != {'p1','p2','q','sig','eps'} or len(e):
                    self.unsupported.append(label+': unsupported Exception row')
                    continue
                key = tuple(sorted((int(e.get('p1')), int(e.get('p2')))))
                if key in records or key[0] == key[1] or not (0 <= key[0] < key[1] < self.particles):
                    raise ValueError('duplicate, self, or invalid exception pair')
                values = tuple(number(e.get(k)) for k in ('q','sig','eps'))
                if values[1] < 0 or values[2] < 0:
                    raise ValueError('negative sigma/epsilon')
                records[key] = values
            maps.append(records)
        x,y = maps
        d = dict(label=label, count1=len(x), count2=len(y), only1=len(set(x)-set(y)), only2=len(set(y)-set(x)),
                 material_pair_count=0, legacy_active_mismatch_count=0, legacy_strict_mismatch_count=0,
                 legacy_charge_product_mismatch_count=0, legacy_sigma_mismatch_count=0, legacy_epsilon_mismatch_count=0,
                 numeric_active_pair_count=0, raw_active_difference_count=0, inert_sigma_difference_count=0,
                 coulomb_activity_changes=0, lj_activity_changes=0, charge_sign_changes=0,
                 channel_activity_transitions={}, max_absolute_active_sigma_nm=0.0, max_relative_active_sigma=0.0,
                 max_absolute_active_epsilon_kj_mol=0.0, max_relative_active_epsilon=0.0,
                 max_absolute_charge_product_e2=0.0, material_examples=[])
        transitions = Counter()
        for k in sorted(set(x)&set(y)):
            qa,sa,ea=x[k]; qb,sb,eb=y[k]
            lj = ea != 0 or eb != 0
            active = qa != 0 or qb != 0 or lj
            qflip=changed_zero(qa,qb); eflip=changed_zero(ea,eb); sign=(qa < 0 < qb) or (qb < 0 < qa)
            d['coulomb_activity_changes']+=qflip; d['lj_activity_changes']+=eflip; d['charge_sign_changes']+=sign
            transitions[f"{int(qa!=0)},{int(ea!=0)}->{int(qb!=0)},{int(eb!=0)}"]+=1
            legacy=[legacy_different(aa,bb) for aa,bb in zip(x[k],y[k])]
            d['legacy_strict_mismatch_count']+=any(legacy)
            d['legacy_charge_product_mismatch_count']+=legacy[0]
            d['legacy_sigma_mismatch_count']+=legacy[1]
            d['legacy_epsilon_mismatch_count']+=legacy[2]
            d['legacy_active_mismatch_count']+=(legacy[0] or legacy[2] or (lj and legacy[1]))
            if not lj and sa != sb:
                d['inert_sigma_difference_count']+=1
            any_active_change=(qa!=qb or ea!=eb or (lj and sa!=sb))
            d['raw_active_difference_count']+=any_active_change
            d['max_absolute_charge_product_e2']=max(d['max_absolute_charge_product_e2'],abs(qa-qb))
            if lj:
                d['max_absolute_active_sigma_nm']=max(d['max_absolute_active_sigma_nm'],abs(sa-sb))
                d['max_relative_active_sigma']=max(d['max_relative_active_sigma'], abs(sa-sb)/max(abs(sa),abs(sb)) if sa or sb else 0.0)
                d['max_absolute_active_epsilon_kj_mol']=max(d['max_absolute_active_epsilon_kj_mol'],abs(ea-eb))
                d['max_relative_active_epsilon']=max(d['max_relative_active_epsilon'],abs(ea-eb)/max(ea,eb))
            bad=(qflip or eflip or sign or not close(qa,qb,POLICY['exception_charge_product_e2'])
                 or not close(ea,eb,POLICY['exception_epsilon_kj_mol'])
                 or (lj and not close(sa,sb,POLICY['exception_sigma_nm'])))
            if bad:
                d['material_pair_count']+=1
                if len(d['material_examples']) < 25:
                    d['material_examples'].append(dict(pair=k,input1=x[k],input2=y[k]))
            elif any_active_change:
                d['numeric_active_pair_count']+=1
        d['channel_activity_transitions']=dict(transitions)
        if d['only1'] or d['only2'] or d['material_pair_count']:
            self.problems.append(label+': material exception mismatch under precision policy')
        self.exceptions.append(d)

    def system(self, a: ET.Element, b: ET.Element) -> None:
        for r in (a,b):
            if r.tag != 'System' or r.get('version') != '1':
                self.unsupported.append('unsupported System schema')
            self.schema(r, {'openmmVersion','type','version'}, ['PeriodicBoxVectors','Particles','Constraints','Forces'], 'System')
        self.attrs(a,b,'System',ignore={'openmmVersion'})
        self.particles=len(a.find('Particles'))
        if self.particles != len(b.find('Particles')):
            self.problems.append('particle count differs')
            return
        self.table(a.find('Particles'),b.find('Particles'),'particle_masses','Particle',(),('mass',),guard=('mass',))
        self.table(a.find('Constraints'),b.find('Constraints'),'constraints','Constraint',('p1','p2'),('d',),True)
        for v in (a.find('PeriodicBoxVectors'),b.find('PeriodicBoxVectors')):
            self.schema(v,set(),['A','B','C'],'periodic_box')
            for c in v:
                self.schema(c,{'x','y','z'},[], 'box_vector')
        for tag in ('A','B','C'):
            self.attrs(a.find('PeriodicBoxVectors/'+tag),b.find('PeriodicBoxVectors/'+tag),'box.'+tag,{'x','y','z'})
        groups=[]
        for r in (a,b):
            d=defaultdict(list)
            f=r.find('Forces')
            if f.attrib:
                self.unsupported.append('unexpected Forces attributes')
            for force in f:
                if force.tag != 'Force':
                    self.unsupported.append('unexpected force tag')
                d[force.get('type')].append(force)
            groups.append(d)
        inventory=[{k:len(v) for k,v in g.items()} for g in groups]
        self.sections['force_inventory']=dict(input1=inventory[0],input2=inventory[1])
        if inventory[0] != inventory[1]:
            self.problems.append('force inventory differs')
        for typ in sorted(set(groups[0])|set(groups[1])):
            for i,(fa,fb) in enumerate(zip(groups[0][typ],groups[1][typ])):
                label=f'{typ}[{i}]'
                if typ=='NonbondedForce':
                    self.nonbonded(fa,fb,label)
                elif typ in FORCE_SCHEMAS:
                    version,child,tag,ids,vals=FORCE_SCHEMAS[typ]
                    for f in (fa,fb):
                        self.schema(f,{'forceGroup','name','type','usesPeriodic','version'},[child],label)
                        if f.get('version')!=version:
                            self.unsupported.append(label+': unsupported serializer version')
                    self.attrs(fa,fb,label,ignore={'name'})
                    self.table(fa.find(child),fb.find(child),label,tag,ids,vals,typ=='HarmonicBondForce',guard=('k',))
                elif typ=='CMMotionRemover':
                    for f in (fa,fb):
                        self.schema(f,{'forceGroup','frequency','name','type','version'},[],label)
                        if f.get('version')!='1':
                            self.unsupported.append(label+': unsupported serializer version')
                    self.attrs(fa,fb,label,ignore={'name'})
                else:
                    self.unsupported.append(label+': unsupported force type')


def audit_system_pair(input1: Path | str, input2: Path | str) -> dict[str, Any]:
    paths=[Path(input1),Path(input2)]
    result: dict[str, Any]=dict(policy=POLICY,status='INCOMPLETE',fully_compared=False,problems=[],unsupported=[],inputs=[])
    audit=Audit()
    try:
        for p in paths:
            result['inputs'].append(dict(path=str(p),sha256=sha256(p)))
        roots=[parse_xml(p) for p in paths]
        result['openmm_versions']=[r.get('openmmVersion') for r in roots]
        audit.system(*roots)
        result.update(problems=audit.problems,unsupported=audit.unsupported,sections=audit.sections,exceptions=audit.exceptions)
        result['fully_compared']=not audit.unsupported
        result['status']='MISMATCH' if audit.problems else 'INCOMPLETE' if audit.unsupported else 'MATCH_WITHIN_TOLERANCE'
    except (OSError, ValueError, TypeError, AttributeError, KeyError, ET.ParseError) as exc:
        result['error']=f'{type(exc).__name__}: {exc}'
    return result


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input1',type=Path); p.add_argument('input2',type=Path)
    p.add_argument('--json',type=Path,help='Write audit JSON; input XML files remain unchanged.')
    args=p.parse_args()
    result=audit_system_pair(args.input1,args.input2)
    text=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if args.json:
        args.json.parent.mkdir(parents=True,exist_ok=True); args.json.write_text(text)
    else:
        print(text,end='')
    print(result['status'])
    return 0 if result['status']=='MATCH_WITHIN_TOLERANCE' else 1 if result['status']=='MISMATCH' else 2

if __name__=='__main__':
    raise SystemExit(main())
