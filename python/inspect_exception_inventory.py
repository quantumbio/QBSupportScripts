#!/usr/bin/env python3
"""Inspect OpenMM NonbondedForce exception keys present on only one side.

Usage, from the campaign root:
    python inspect_exception_inventory.py --root . 4ibf 5yz2 6d9x

Writes one CSV per case to --out. Pure Python standard library; reads existing
OpenMM XML and optional cpp companion AMBER parm7 without changing them.

Interpretation: This is a structural/parameter inventory audit, NOT a substitute
for evaluation of the two Hamiltonians at identical coordinates in OpenMM.
"""
from __future__ import annotations

import argparse
import bisect
from collections import Counter, deque
import csv
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET

FIELDS = (
    'case','side','i','j','atom_i','atom_j','exception_qprod_e2',
    'exception_sigma_nm','exception_epsilon_kj_mol',
    'ordinary_qprod_e2','ordinary_sigma_nm','ordinary_epsilon_kj_mol',
    'exception_coulomb_active','exception_lj_active',
    'ordinary_coulomb_active','ordinary_lj_active',
    'both_paths_inactive','in_raw_amber_exclusion_list','bond_graph_distance',
    'dihedral_endpoint_records','dihedral_endpoint_allowed','dihedral_endpoint_suppressed',
    'dihedral_endpoint_improper','possible_physics_consequence',
)


def read_nonbonded(path: Path):
    # Keep consistent with the fail-closed XML comparator's DTD restrictions.
    raw = path.read_bytes()
    if b'<!DOCTYPE' in raw.upper() or b'<!ENTITY' in raw.upper():
        raise ValueError(f'{path}: XML entity declarations not supported')
    root = ET.fromstring(raw)
    if root.tag != 'System':
        raise ValueError(f'{path}: expected serialized System')
    forces = [f for f in root.findall('./Forces/Force') if f.get('type') == 'NonbondedForce']
    if len(forces) != 1:
        raise ValueError(f'{path}: expected exactly one NonbondedForce; found {len(forces)}')
    force = forces[0]
    particles = []
    for node in force.findall('./Particles/Particle'):
        particles.append(tuple(float(node.get(k)) for k in ('q','sig','eps')))
    exceptions = {}
    for node in force.findall('./Exceptions/Exception'):
        pair = tuple(sorted((int(node.get('p1')),int(node.get('p2')))))
        if pair in exceptions or pair[0] == pair[1] or pair[0] < 0 or pair[1] >= len(particles):
            raise ValueError(f'{path}: duplicate or invalid exception pair {pair}')
        exceptions[pair] = tuple(float(node.get(k)) for k in ('q','sig','eps'))
    return particles, exceptions


_FORMAT = re.compile(r'\(\s*(\d+)\s*([AaIiDdEeFf])\s*(\d+)(?:\.\d+)?\s*\)')


def read_amber_topology(path: Path, natoms: int):
    """Read basic AMBER atom/residue/connectivity/exclusion fields; optional."""
    blocks = {}
    flag = None
    width = None
    typ = None
    with path.open(encoding='ascii', errors='replace') as handle:
        for line in handle:
            if line.startswith('%FLAG '):
                flag = line[6:].strip()
                width, typ = None, None
                if flag in {
                    'ATOM_NAME','RESIDUE_LABEL','RESIDUE_POINTER',
                    'BONDS_INC_HYDROGEN','BONDS_WITHOUT_HYDROGEN',
                    'DIHEDRALS_INC_HYDROGEN','DIHEDRALS_WITHOUT_HYDROGEN',
                    'NUMBER_EXCLUDED_ATOMS','EXCLUDED_ATOMS_LIST',
                }:
                    blocks[flag] = []
            elif line.startswith('%FORMAT('):
                match = _FORMAT.search(line)
                if match:
                    typ, width = match.group(2).upper(), int(match.group(3))
            elif line.startswith('%') or flag not in blocks:
                continue
            elif width:
                for start in range(0, len(line.rstrip('\r\n')), width):
                    token = line[start:start+width].strip()
                    if not token:
                        continue
                    if typ == 'I':
                        blocks[flag].append(int(token))
                    elif typ == 'A':
                        blocks[flag].append(token)
    names = blocks.get('ATOM_NAME', [])
    residue_labels = blocks.get('RESIDUE_LABEL', [])
    residue_pointers = blocks.get('RESIDUE_POINTER', [])
    if len(names) != natoms:
        raise ValueError(f'{path}: atom names {len(names)} != XML particles {natoms}')

    def atom_name(idx):
        ridx = bisect.bisect_right(residue_pointers,idx+1)-1
        residue = f'{residue_labels[ridx]}{ridx+1}' if 0 <= ridx < len(residue_labels) else '?'
        return f'{residue}:{names[idx]}'

    graph = [set() for _ in range(natoms)]
    for field in ('BONDS_INC_HYDROGEN','BONDS_WITHOUT_HYDROGEN'):
        values = blocks.get(field,[])
        if len(values) % 3:
            raise ValueError(f'{path}: malformed {field}')
        for n in range(0,len(values),3):
            i,j = abs(values[n])//3, abs(values[n+1])//3
            if not 0 <= i < natoms or not 0 <= j < natoms:
                raise ValueError(f'{path}: invalid bond indices')
            graph[i].add(j); graph[j].add(i)
    excluded = set()
    ex_counts = blocks.get('NUMBER_EXCLUDED_ATOMS', [])
    ex_items = blocks.get('EXCLUDED_ATOMS_LIST', [])
    if ex_counts and len(ex_counts) != natoms:
        raise ValueError(f'{path}: malformed NUMBER_EXCLUDED_ATOMS')
    if ex_counts:
        cursor = 0
        for i, n in enumerate(ex_counts):
            for val in ex_items[cursor:cursor+n]:
                if val > 0:
                    j = val - 1
                    excluded.add(tuple(sorted((i,j))))
            cursor += n
        if cursor != len(ex_items):
            raise ValueError(f'{path}: malformed excluded atoms list')
    dihedrals = {}
    for field in ('DIHEDRALS_INC_HYDROGEN','DIHEDRALS_WITHOUT_HYDROGEN'):
        values = blocks.get(field,[])
        if len(values)%5:
            raise ValueError(f'{path}: malformed {field}')
        for n in range(0,len(values),5):
            a,b,c,d,term=values[n:n+5]
            key=tuple(sorted((abs(a)//3,abs(d)//3)))
            flags=dihedrals.setdefault(key,Counter())
            flags['count']+=1
            # AMBER signed third atom denotes suppress-end interactions;
            # negative fourth marks an improper term.
            if c<0: flags['suppressed']+=1
            elif d<0: flags['improper']+=1
            else: flags['allowed']+=1
    return atom_name, graph, excluded, dihedrals


def distance(graph, src, dst, max_edges=4):
    if src == dst:
        return 0
    seen={src}
    queue=deque([(src,0)])
    while queue:
        current,n=queue.popleft()
        if n>=max_edges:
            continue
        for nxt in graph[current]:
            if nxt == dst:
                return n+1
            if nxt not in seen:
                seen.add(nxt)
                queue.append((nxt,n+1))
    return '>4 or disconnected'


def make_rows(case, side, pairset, one_particles, one_exceptions, other_particles, topology):
    atom_name, graph, excluded, dihedrals = topology if topology else (None,None,None,None)
    for i,j in sorted(pairset):
        q,s,e = one_exceptions[(i,j)]
        qi,si,ei = other_particles[i]
        qj,sj,ej = other_particles[j]
        ordinary_q = qi*qj
        ordinary_s = (si+sj)/2
        ordinary_e = math.sqrt(ei*ej)
        both_inactive = (q == 0 and e == 0 and ordinary_q == 0 and ordinary_e == 0)
        dflags = dihedrals.get((i,j),Counter()) if dihedrals is not None else Counter()
        yield dict(
            case=case,side=side,i=i,j=j,
            atom_i=atom_name(i) if atom_name else '',
            atom_j=atom_name(j) if atom_name else '',
            exception_qprod_e2=repr(q),exception_sigma_nm=repr(s),exception_epsilon_kj_mol=repr(e),
            ordinary_qprod_e2=repr(ordinary_q),ordinary_sigma_nm=repr(ordinary_s),ordinary_epsilon_kj_mol=repr(ordinary_e),
            exception_coulomb_active=bool(q),exception_lj_active=bool(e),
            ordinary_coulomb_active=bool(ordinary_q),ordinary_lj_active=bool(ordinary_e),
            both_paths_inactive=both_inactive,
            in_raw_amber_exclusion_list=((i,j) in excluded) if excluded is not None else '',
            bond_graph_distance=distance(graph,i,j) if graph is not None else '',
            dihedral_endpoint_records=dflags['count'],
            dihedral_endpoint_allowed=dflags['allowed'],
            dihedral_endpoint_suppressed=dflags['suppressed'],
            dihedral_endpoint_improper=dflags['improper'],
            possible_physics_consequence='likely inert on this nonbonded pair' if both_inactive else 'requires Hamiltonian/energy verification',
        )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cases',nargs='+',help='Structure directories, e.g. 4ibf 5yz2 6d9x')
    p.add_argument('--root',type=Path,default=Path('.'))
    p.add_argument('--out',type=Path,default=Path('precision_details/exception_inventory'))
    p.add_argument('--no-topology',action='store_true',help='Skip optional AMBER topology annotation')
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    for case in args.cases:
        parent = args.root/case
        x_part,x_ex=read_nonbonded(parent/'cpp/openmm_system.xml')
        y_part,y_ex=read_nonbonded(parent/'python/openmm_system.xml')
        if len(x_part)!=len(y_part):
            raise ValueError(f'{case}: particle counts disagree')
        onlycpp=set(x_ex)-set(y_ex)
        onlypy=set(y_ex)-set(x_ex)
        topo=None
        if not args.no_topology:
            candidates=sorted((parent/'cpp').glob('*.parm7'))
            if not candidates:
                candidates=sorted(parent.glob('*.parm7'))
            if len(candidates)==1:
                topo=read_amber_topology(candidates[0],len(x_part))
                print(f'  using topology {candidates[0]}')
            elif candidates:
                print(f'  {case}: multiple parm7 candidates; topology annotation omitted')
            else:
                print(f'  {case}: no cpp/root parm7 available; topology annotation omitted')
        rows=list(make_rows(case,'C++ only',onlycpp,x_part,x_ex,y_part,topo))
        rows+=list(make_rows(case,'Python only',onlypy,y_part,y_ex,x_part,topo))
        dest=args.out/f'{case}_unique_exceptions.csv'
        with dest.open('w',newline='',encoding='utf-8') as handle:
            writer=csv.DictWriter(handle,fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        nontrivial=sum(row['both_paths_inactive'] is False for row in rows)
        raw_excluded=sum(row['in_raw_amber_exclusion_list'] is True for row in rows)
        print(f'{case}: C++ only={len(onlycpp)}, Python only={len(onlypy)}, '
              f'potentially interacting={nontrivial}, AMBER excluded={raw_excluded}; {dest}')
        for row in rows[:8]:
            print(f"  {row['side']} ({row['i']}, {row['j']}) {row['atom_i']} / {row['atom_j']} "
                  f"bond_distance={row['bond_graph_distance']} "
                  f"raw_excluded={row['in_raw_amber_exclusion_list']} "
                  f"C++/Python exception={row['exception_qprod_e2']}, {row['exception_epsilon_kj_mol']} "
                  f"ordinary={row['ordinary_qprod_e2']}, {row['ordinary_epsilon_kj_mol']}")
    return 0

if __name__=='__main__':
    raise SystemExit(main())
