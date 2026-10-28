#!/usr/bin/env python3
"""Diagnose SystemModelSnapshotWriter v1-v3 JSON or gzip snapshots.

Usage: python3 model_comparisons.py A.json.gz B.json.gz [--details]
Indices identify flat atoms/residues at these checkpoints; they do not establish
identity across independent model reconstruction. Numeric enums stay numeric.
Only Python's standard library is required. Exit 0 for a successful comparison
(including differences), 2 for invalid input or arguments.

Schema 1/2 fields are compared when present; absent features remain unavailable.
The current writer emits v3 policies: ordinary exclusions and independent special
1-4 energy multipliers with optional unscaled LJ [A,B] coefficients. Pair identity
uses flat atom indices, not collection position. Ring identity uses the ordered
perimeter modulo rotation/reversal; raw ID/order differences remain diagnostic.
"""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import gzip
import json
import math
import reprlib
import sys

SCHEMA = 'quantumbio.system_model_snapshot'
MISSING = object()
POLICIES = ('nonbonded_exclusions', 'special_pairs')
PAIR_FIELDS = ('atomA_index', 'atomB_index')


def cycle_key(atoms):
    """A perimeter, modulo rotation and reversal, never an unordered atom set."""
    if (not isinstance(atoms, list) or len(atoms) < 3 or
            any(type(a) is not int or a < 0 for a in atoms) or len(set(atoms)) != len(atoms)):
        raise ValueError('invalid ring perimeter')
    start = atoms.index(min(atoms))
    forward = atoms[start:] + atoms[:start]
    reverse = [forward[0]] + forward[:0:-1]
    return min(tuple(forward), tuple(reverse))


def display(value):
    if value is MISSING:
        return '<missing>'
    # Bound the traversal/formatting too: a membership list may have 30k atoms.
    result = reprlib.repr(value)
    if isinstance(value, list) and len(value) > 6:
        result += f' (length {len(value)})'
    return result


def equal(a, b):
    """Float noise tolerance is separate from coordinate displacement tolerance."""
    if type(a) is not type(b):
        return False
    if isinstance(a, float):
        return a == b or math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)
    if isinstance(a, list):
        return len(a) == len(b) and all(equal(x, y) for x, y in zip(a, b))
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(equal(v, b[k]) for k, v in a.items())
    return a == b


def validate(data, label):
    if not isinstance(data, dict):
        raise ValueError(f'{label}: expected a JSON object')
    header = data.get('snapshot')
    if not isinstance(header, dict) or header.get('schema') != SCHEMA:
        raise ValueError(f'{label}: missing or incompatible $.snapshot.schema (expected {SCHEMA})')
    version = header.get('schema_version')
    if type(version) is not int or version not in (1, 2, 3):
        raise ValueError(f'{label}: unsupported $.snapshot.schema_version: {version!r}')
    if not isinstance(data.get('system_model'), dict):
        raise ValueError(f'{label}: missing or malformed $.system_model')
    for name in POLICIES:
        policy = data['system_model'].get(name, MISSING)
        if policy is MISSING:
            continue  # Feature availability, rather than version alone, is authoritative.
        path = '$.system_model.' + name
        if not isinstance(policy, dict):
            raise ValueError(f'{label}: {path} must be an object')
        if 'pairs' in policy and not isinstance(policy['pairs'], list):
            raise ValueError(f'{label}: {path}.pairs must be an array')
        if 'complete' in policy and type(policy['complete']) is not bool:
            raise ValueError(f'{label}: {path}.complete must be a boolean')
        if 'provenance' in policy and policy['provenance'] not in (
                'unavailable', 'imported', 'generated', 'mixed'):
            raise ValueError(f'{label}: unsupported {path}.provenance: {policy["provenance"]!r}')
    # Collections can be absent in older checkpoints, but cannot have wrong types.
    def visit(obj, path):
        if isinstance(obj, dict):
            if 'counts' in obj:
                counts = obj['counts']
                if not isinstance(counts, dict) or any(
                        type(v) is not int or v < 0 for v in counts.values()):
                    raise ValueError(f'{label}: {path}.counts must contain nonnegative integers')
            for key, val in obj.items():
                child_path = path + '.' + key
                if key == 'partial_charge_imported' and type(val) is not bool:
                    raise ValueError(f'{label}: {child_path} must be a boolean')
                if identity_fields(child_path) and not isinstance(val, list):
                    raise ValueError(f'{label}: {child_path} must be an array')
                visit(val, child_path)
        elif isinstance(obj, list):
            for val in obj:
                if isinstance(val, (dict, list)):
                    visit(val, path + '[]')
    for key in ('strands', 'atoms', 'residues', 'bonds', 'residue_connections',
                'rings', 'fragments', 'alternate_groups', 'alternate_group_bridge'):
        if key in data['system_model'] and not isinstance(data['system_model'][key], list):
            raise ValueError(f'{label}: $.system_model.{key} must be an array')
    visit(data, '$')


def load_snapshot(path):
    try:
        opener = gzip.open if str(path).endswith('.gz') else open
        with opener(path, 'rt', encoding='utf-8') as stream:
            data = json.load(stream, parse_constant=lambda s: (_ for _ in ()).throw(
                ValueError(f'nonstandard JSON number {s}; writer uses strings')))
        validate(data, str(path))
        return data
    except (OSError, ValueError, EOFError, UnicodeError, RecursionError) as exc:
        raise ValueError(f'Cannot read snapshot {path}: {exc}') from exc


@lru_cache(maxsize=1024)
def identity_fields(path):
    """Use writer fields, never a generic 'index' for unknown collections."""
    flat = '$.system_model.'
    if path in (flat + 'atoms', flat + 'residues', flat + 'strands'):
        return ('index',)
    if path == flat + 'bonds':
        return ('atomA_index', 'atomB_index')
    if path in (flat + name + '.pairs' for name in POLICIES):
        return PAIR_FIELDS
    if path == flat + 'residue_connections':
        return ('residue_a_index', 'residue_b_index')
    if path in (flat + 'rings', flat + 'fragments'):
        return ('atom_indices',)
    if path == flat + 'alternate_group_bridge':
        return ('atom_index',)
    if path == flat + 'alternate_groups':
        return ('group_id',)
    if path.endswith('.alternate_groups[].alternates'):
        return ('alternate',)
    if path == '$.constitutional_graph.atomic_graph.nodes':
        return ('atom_index',)
    if path == '$.constitutional_graph.atomic_graph.edges':
        return ('atom_a_index', 'atom_b_index')
    if path == '$.constitutional_graph.residue_connectivity_graph.nodes':
        return ('residue_index',)
    if path == '$.constitutional_graph.residue_connectivity_graph.edges':
        return ('residue_a_index', 'residue_b_index')
    if '.gemmi.' in path:
        last = path.rsplit('.', 1)[-1]
        return {'models': ('num',), 'chains': ('name',),
                'residues': ('sequence_number', 'insertion_code', 'name', 'segment'),
                'atoms': ('name', 'alternate_location'), 'conect_map': ('atom_serial',),
                'assemblies': ('name',)}.get(last)
    return None


def record_key(record, fields, path=''):
    if not isinstance(record, dict) or any(k not in record for k in fields):
        raise ValueError('missing identity fields')
    values = tuple(record[k] for k in fields)
    if fields == ('atom_indices',):
        if path == '$.system_model.rings':
            return cycle_key(values[0])
        return tuple(sorted(values[0]))
    if len(fields) == 2 and fields[0].endswith(('_a_index', 'A_index')):
        return tuple(sorted(values))
    hash(values)
    return values


def inventory(data, label):
    """Aggregate nested array lengths by full wildcard path, without copying data."""
    sizes = defaultdict(lambda: [0, 0])  # total entries, number of arrays
    counts = {}
    warnings = Counter()
    def visit(obj, path):
        if isinstance(obj, dict):
            for key, val in obj.items():
                if key == 'counts' and isinstance(val, dict):
                    for name, count in val.items():
                        cp = path + '.counts.' + name
                        counts[cp] = count
                        actual = obj.get(name, MISSING)
                        if name == 'gemmi_models':
                            actual = obj.get('gemmi', {}).get('models', MISSING)
                        if isinstance(actual, list) and count != len(actual):
                            warnings[f'{label}: {cp}={count}, actual length={len(actual)}'] += 1
                visit(val, path + '.' + key)
        elif isinstance(obj, list):
            sizes[path][0] += len(obj)
            sizes[path][1] += 1
            for val in obj:
                if isinstance(val, (dict, list)):
                    visit(val, path + '[]')
    visit(data, '$')
    return sizes, counts, warnings


class Comparison:
    def __init__(self, a, b, max_examples=3, coordinate_tolerance=1e-6):
        validate(a, 'A')
        validate(b, 'B')
        self.a, self.b = a, b
        self.max_examples = max_examples
        self.tolerance = coordinate_tolerance
        self.fields = Counter()
        self.examples = defaultdict(list)
        self.metadata = []
        self.collections = defaultdict(Counter)
        self.shapes = Counter()
        self.warnings = Counter()
        self.raw_fields = Counter()
        self.raw_examples = defaultdict(list)
        self.noise = Counter()
        self.invalid_collections = set()
        self.coords = defaultdict(lambda: [0, 0.0, 0.0, 0, 0])
        self.sizes_a, self.counts_a, wa = inventory(a, 'A')
        self.sizes_b, self.counts_b, wb = inventory(b, 'B')
        self.warnings.update(wa)
        self.warnings.update(wb)
        self.rings_a = self.audit(a, 'A')
        self.rings_b = self.audit(b, 'B')
        self.compare(a, b, '$', '$')

    def audit(self, data, label):
        """Keep malformed/duplicate raw records visible before identity matching."""
        model = data['system_model']
        atom_ids = {a.get('index') for a in model.get('atoms', []) if isinstance(a, dict)
                    and type(a.get('index')) is int}
        for name in POLICIES:
            policy = model.get(name, {})
            path = '$.system_model.' + name + '.pairs'
            pairs = policy.get('pairs', [])
            seen_pairs = set()
            if policy.get('complete') is True and policy.get('provenance') == 'unavailable':
                self.warnings[f'{label}: {name}: complete policy has unavailable provenance'] += 1
            if pairs and policy.get('provenance') == 'unavailable':
                self.warnings[f'{label}: {name}: pairs stored with unavailable provenance'] += 1
            for i, pair in enumerate(pairs):
                valid = isinstance(pair, dict) and all(type(pair.get(k)) is int for k in PAIR_FIELDS)
                if valid:
                    endpoints = [pair[k] for k in PAIR_FIELDS]
                    valid = endpoints[0] != endpoints[1] and all(a >= 0 for a in endpoints)
                    if 'atoms' in model:
                        valid &= all(a in atom_ids for a in endpoints)
                    key = tuple(sorted(endpoints))
                    if key in seen_pairs:
                        self.warnings[f'{label}: {path}: duplicate identity keys; positional comparison'] += 1
                        self.invalid_collections.add(path)
                    seen_pairs.add(key)
                if name == 'special_pairs' and isinstance(pair, dict):
                    for field in ('coulombMultiplier', 'ljMultiplier'):
                        number = pair.get(field)
                        valid &= type(number) in (int, float) and math.isfinite(number) and number >= 0
                    coefficients = pair.get('ljCoefficients', MISSING)
                    valid &= coefficients is None or (isinstance(coefficients, list) and
                        len(coefficients) == 2 and all(type(v) in (int, float) and
                            math.isfinite(v) and v >= 0 for v in coefficients))
                if not valid:
                    self.warnings[f'{label}: {path}: invalid pair record at position {i}; positional comparison'] += 1
                    self.invalid_collections.add(path)
        rings = {}
        perimeters = set()
        expected_atoms, expected_edges, expected_residues = defaultdict(set), defaultdict(set), defaultdict(set)
        for i, ring in enumerate(model.get('rings', [])):
            try:
                ring_id = ring.get('index', i)
                if type(ring_id) is not int or ring_id < 0 or ring_id in rings:
                    raise ValueError('invalid or duplicate ring index')
                perimeter = ring['atom_indices']
                rings[ring_id] = cycle_key(perimeter)
                if rings[ring_id] in perimeters:
                    self.warnings[f'{label}: rings: duplicate perimeter identity; positional comparison'] += 1
                    self.invalid_collections.add('$.system_model.rings')
                perimeters.add(rings[ring_id])
                if type(ring.get('residue_index')) is int:
                    expected_residues[ring['residue_index']].add(ring_id)
                if 'atoms' in model and any(a not in atom_ids for a in perimeter):
                    self.warnings[f'{label}: rings: invalid perimeter atom references'] += 1
                for a, b in zip(perimeter, perimeter[1:] + perimeter[:1]):
                    expected_atoms[a].add(ring_id)
                    expected_edges[tuple(sorted((a, b)))].add(ring_id)
            except (ValueError, TypeError, KeyError, AttributeError):
                self.warnings[f'{label}: rings: invalid perimeter/index at position {i}'] += 1
                self.invalid_collections.add('$.system_model.rings')
        for name in ('atoms', 'bonds', 'residues'):
            for record in model.get(name, []):
                if not isinstance(record, dict):
                    continue
                for field in ('ringParticipation', 'systemRingIndices'):
                    if field not in record:
                        continue
                    refs = record[field]
                    if (not isinstance(refs, list) or any(type(r) is not int for r in refs)):
                        self.warnings[f'{label}: {name}.{field}: malformed ring references'] += 1
                        continue
                    if len(refs) != len(set(refs)):
                        self.warnings[f'{label}: {name}.{field}: duplicate ring references'] += 1
                    if 'rings' in model and any(r not in rings for r in refs):
                        self.warnings[f'{label}: {name}.{field}: invalid ring references'] += 1
                    if field == 'ringParticipation' and 'ring_count' in record and record['ring_count'] != len(refs):
                        self.warnings[f'{label}: {name}: ring_count disagrees with ringParticipation'] += 1
                    if 'rings' in model and name in ('atoms', 'bonds', 'residues'):
                        if name == 'atoms':
                            expected = expected_atoms.get(record.get('index'), set())
                        elif name == 'bonds':
                            try:
                                expected = expected_edges.get(record_key(record, PAIR_FIELDS), set())
                            except (ValueError, TypeError):
                                continue
                        else:
                            expected = expected_residues.get(record.get('index'), set())
                        if set(refs) != expected:
                            self.warnings[f'{label}: {name}.{field}: disagrees with registered perimeters'] += 1
                        if name == 'atoms' and 'smallest_ring' in record:
                            smallest = min((len(rings[r]) for r in expected), default=0)
                            if record['smallest_ring'] != smallest:
                                self.warnings[f'{label}: atoms.smallest_ring: disagrees with registered perimeters'] += 1
        if 'bonds' in model and 'rings' in model:
            edges = set()
            for bond in model['bonds']:
                try:
                    edges.add(record_key(bond, PAIR_FIELDS))
                except (ValueError, TypeError):
                    pass
            missing = expected_edges.keys() - edges
            if missing:
                self.warnings[f'{label}: rings: {len(missing)} perimeter edges absent from bonds'] += 1
        return rings

    def raw_change(self, path, identifier, a, b):
        self.raw_fields[path] += 1
        if len(self.raw_examples[path]) < self.max_examples:
            self.raw_examples[path].append(f'{identifier}: A={display(a)}, B={display(b)}')

    def ring_references(self, a, b, path, identifier):
        try:
            if not isinstance(a, list) or not isinstance(b, list):
                return None
            if any(type(r) is not int for r in a + b):
                return None
            left = Counter(self.rings_a[r] for r in a)
            right = Counter(self.rings_b[r] for r in b)
        except (KeyError, TypeError):
            return None  # Undefined IDs cannot be interpreted; compare raw values.
        if left != right:
            return self.changed(path, identifier,
                                {'ring_ids': a, 'perimeters': list(left.elements())},
                                {'ring_ids': b, 'perimeters': list(right.elements())})
        if a != b:
            self.raw_change(path, identifier, a, b)
        return False

    def changed(self, path, identifier, a, b):
        self.fields[path] += 1
        if len(self.examples[path]) < self.max_examples:
            self.examples[path].append(f'{identifier}: A={display(a)}, B={display(b)}')
        if '[]' not in path and '.counts.' not in path:
            self.metadata.append(f'{path}: A={display(a)}, B={display(b)}')
        return True

    def coordinates(self, a, b, path, identifier):
        stats = self.coords[path]
        if (not isinstance(a, list) or not isinstance(b, list) or len(a) != 3 or len(b) != 3
                or any(type(v) not in (int, float) or not math.isfinite(v) for v in a + b)):
            stats[4] += 1
            if not equal(a, b):
                return self.changed(path, identifier, a, b)
            return False
        displacement = math.dist(a, b)
        stats[0] += 1
        stats[1] += displacement * displacement
        stats[2] = max(stats[2], displacement)
        stats[3] += displacement > self.tolerance
        return displacement > self.tolerance

    def compare(self, a, b, path, identifier):
        if path.startswith('$.system_model.') and path.endswith(('.ringParticipation', '.systemRingIndices')):
            result = self.ring_references(a, b, path, identifier)
            if result is not None:
                return result
        if path.endswith(('.atoms[].coordinates', '.atoms[].position')):
            return self.coordinates(a, b, path, identifier)
        if isinstance(a, dict) and isinstance(b, dict):
            changed = False
            for key in sorted(a.keys() | b.keys()):
                changed |= self.compare(a.get(key, MISSING), b.get(key, MISSING),
                                        path + '.' + key, identifier)
            return changed
        if isinstance(a, list) and isinstance(b, list):
            if len(a) != len(b):
                self.shapes[path] += 1
            # Scalar arrays and Eigen matrices are one field per affected record.
            if not identity_fields(path) and not any(isinstance(x, dict) for x in a) and not any(isinstance(x, dict) for x in b):
                if not equal(a, b):
                    return self.changed(path, identifier, a, b)
                return False
            return self.collection(a, b, path, identifier)
        if not equal(a, b):
            return self.changed(path, identifier, a, b)
        if isinstance(a, float) and isinstance(b, float) and a != b:
            self.noise[path] += 1
        return False

    def collection(self, a, b, path, parent):
        summary = self.collections[path]
        fields = identity_fields(path)
        maps = []
        if path in self.invalid_collections:
            fields = None
        if fields:
            requested_fields = fields
            for label, records in (('A', a), ('B', b)):
                mapping = {}
                duplicates = 0
                try:
                    for i, record in enumerate(records):
                        key = record_key(record, requested_fields, path)
                        if key in mapping:
                            duplicates += 1
                        else:
                            mapping[key] = (i, record)
                except (ValueError, TypeError, KeyError, AttributeError):
                    self.warnings[f'{label}: {path}: unusable identity; positional comparison'] += 1
                    fields = None
                    break
                if duplicates:
                    self.warnings[f'{label}: {path}: duplicate identity keys ({duplicates} extra records); positional comparison'] += 1
                    fields = None
                maps.append(mapping)
        if not fields:
            summary['positional'] += 1
            summary['extra_A'] += max(0, len(a) - len(b))
            summary['extra_B'] += max(0, len(b) - len(a))
            modified = 0
            for i, (x, y) in enumerate(zip(a, b)):
                modified += self.compare(x, y, path + '[]', f'{parent}/{path}[position={i}]')
            summary['modified'] += modified
            return len(a) != len(b) or bool(modified)
        ma, mb = maps
        summary['matched'] += 1
        summary['removed'] += len(ma.keys() - mb.keys())
        summary['added'] += len(mb.keys() - ma.keys())
        for kind, keys, records in (('removed', ma.keys() - mb.keys(), ma),
                                    ('added', mb.keys() - ma.keys(), mb)):
            example_path = path + '[].' + kind
            for key in sorted(keys, key=repr)[:self.max_examples]:
                self.examples[example_path].append(f'{path}[identity={display(key)}]: {display(records[key][1])}')
        modified = 0
        # Insertion order of A makes capped examples reproducible and avoids sorting huge maps.
        for key, (i, x) in ma.items():
            if key not in mb:
                continue
            j, y = mb[key]
            summary['moved'] += i != j
            identifier = f'{path}[{",".join(fields)}={display(key)}]'
            if '[]' in path:
                identifier = f'{parent} -> {identifier}'
            # Vector positions are bookkeeping for endpoint/membership identities.
            # Endpoint orientation is recorded separately, not a chemical change.
            ignore = {'index'} if fields != ('index',) else set()
            if 'index' in ignore and x.get('index') != y.get('index'):
                self.raw_change(path + '[].index', identifier, x.get('index', MISSING), y.get('index', MISSING))
            if len(fields) == 2 and fields[0].endswith(('_a_index', 'A_index')):
                summary['orientation'] += tuple(x[k] for k in fields) != tuple(y[k] for k in fields)
                # Stereo, image offsets and signed rotations retain endpoint orientation.
                directional = path == '$.system_model.bonds' and any(
                    r.get(k, 0) != 0 for r in (x, y) for k in
                    ('stereo', 'symmetry_index', 'symmetry_x', 'symmetry_y', 'symmetry_z', 'rotation_angle'))
                if not directional:
                    ignore.update(fields)
            record_changed = False
            for field in sorted(x.keys() | y.keys()):
                if field not in ignore:
                    if fields == ('atom_indices',) and field == 'atom_indices':
                        summary['membership_order'] += x[field] != y[field]
                        if path == '$.system_model.rings':
                            if x[field] != y[field]:
                                self.raw_change(path + '[].' + field, identifier, x[field], y[field])
                            continue
                        # Fragment membership order has no documented cycle semantics.
                    record_changed |= self.compare(x.get(field, MISSING), y.get(field, MISSING),
                                                   path + '[].' + field, identifier)
            modified += record_changed
        summary['modified'] += modified
        return bool(summary['added'] or summary['removed'] or modified)

    def report(self, filenames=('A', 'B'), details=False):
        lines = ['Summary']
        for label, name, data in zip(('A', 'B'), filenames, (self.a, self.b)):
            model = data['system_model']
            lines.append(f'  {label}: {name}')
            lines.append(f'     {SCHEMA} v{data["snapshot"]["schema_version"]}; '
                         f'name={display(model.get("name", MISSING))}; entryID={display(model.get("entryID", MISSING))}')
        lines.append(f'  {len(self.fields)} changed field paths; '
                     f'{sum(v["modified"] for v in self.collections.values())} modified records '
                     '(nested collections counted separately).')
        lines.append('  Flat indices are checkpoint identities; enums are stored numeric values.')
        lines.append('  Missing fields mean unavailable, not empty/false/zero. Atom correspondence uses flat index; no inferred renumbering.')
        def section(title, rows, limit=12):
            lines.extend(['', title])
            if not rows:
                lines.append('  None.')
            else:
                shown = rows if details else rows[:limit]
                lines.extend('  ' + row for row in shown)
                if len(shown) < len(rows):
                    lines.append(f'  ... {len(rows)-len(shown)} more; use --details.')
        structural = []
        for path, s in sorted(self.collections.items()):
            if s['positional']:
                if details or s['modified'] or s['extra_A'] or s['extra_B']:
                    structural.append(f'{path}: positional comparison (no reliable unique identity); '
                                      f'modified={s["modified"]}, unmatched A={s["extra_A"]}, B={s["extra_B"]}')
            elif details or s['modified'] or s['added'] or s['removed'] or s['moved'] or s['orientation'] or s['membership_order'] or path == '$.system_model.bonds':
                structural.append(f'{path}: added={s["added"]}, removed={s["removed"]}, modified={s["modified"]}; '
                                  f'reordered positions={s["moved"]}, reversed orientation={s["orientation"]}')
                if s['membership_order']:
                    structural.append(f'{path}: {s["membership_order"]} records with reordered atom membership')
                if path == '$.system_model.bonds' and not (s['added'] or s['removed']):
                    structural.append('Bond endpoint-pair sets are identical.')
                if s['added'] or s['removed']:
                    for kind in ('removed', 'added'):
                        structural.extend(f'{kind}: {ex}' for ex in self.examples[path + '[].' + kind])
        for path, (n, ss, maximum, exceed, invalid) in sorted(self.coords.items()):
            structural.append(f'{path}: {n} matched atoms, max displacement={maximum:.6g}, '
                              f'RMS={math.sqrt(ss/n) if n else 0:.6g}; '
                              f'{exceed} exceeding tolerance {self.tolerance:g}; invalid/missing={invalid}')
        # Always show coordinates even when structural summaries are capped.
        section('Important structural changes', structural, 20)
        section('Metadata changes', self.metadata)
        availability = []
        for label, data in (('A', self.a), ('B', self.b)):
            model = data['system_model']
            for name in POLICIES:
                policy = model.get(name)
                if policy is None:
                    availability.append(f'{label}: {name}: unavailable (not serialized)')
                else:
                    count = len(policy['pairs']) if 'pairs' in policy else 'unavailable'
                    availability.append(f'{label}: {name}: provenance={display(policy.get("provenance", MISSING))}, '
                                        f'complete={display(policy.get("complete", MISSING))}, pairs={count}')
            atoms = model.get('atoms', [])
            available = sum('partial_charge_imported' in atom for atom in atoms if isinstance(atom, dict))
            availability.append(f'{label}: imported-charge provenance: {available}/{len(atoms)} atoms serialized; '
                                'absent flags are unavailable')
        section('Policy and imported-charge availability', availability)
        section('Special-pair semantics', ['Independent dimensionless Coulomb/LJ energy multipliers; '
                'ljCoefficients=[A,B] are unscaled (Angstrom, kcal/mol), never reordered or inverted; '
                'null requires finalized MMData coefficients.'])
        rows = []
        unchanged = []
        def row(path, av, bv):
            delta = f'{bv-av:+d}' if type(av) is int and type(bv) is int else 'n/a'
            return f'{path:<76} {str(av):>10} {str(bv):>10} {delta:>10}'
        rows.append(f'{"Path (reported count / actual array length)":<76} {"A":>10} {"B":>10} {"B - A":>10}')
        for path in sorted(self.counts_a.keys() | self.counts_b.keys()):
            av, bv = self.counts_a.get(path, '<missing>'), self.counts_b.get(path, '<missing>')
            if av != bv or details:
                rows.append(row(path + ' [reported]', av, bv))
            if av == bv and path in ('$.system_model.counts.' + k for k in ('atoms', 'residues', 'bonds')):
                unchanged.append(row(path + ' [reported]', av, bv))
        for path in sorted(self.sizes_a.keys() | self.sizes_b.keys()):
            av, bv = self.sizes_a.get(path), self.sizes_b.get(path)
            if av != bv or details:
                rows.append(row(path + ' [length total]', av[0] if av else '<missing>', bv[0] if bv else '<missing>'))
                if (av and av[1] != 1) or (bv and bv[1] != 1):
                    rows.append(row(path + ' [number of arrays]', av[1] if av else '<missing>', bv[1] if bv else '<missing>'))
            if av == bv and path in ('$.system_model.' + k for k in (
                    'atoms', 'residues', 'bonds', 'nonbonded_exclusions.pairs', 'special_pairs.pairs')):
                unchanged.append(row(path + ' [actual length]', av[0], bv[0]))
        for path, n in sorted(self.shapes.items()):
            if '[]' in path:
                rows.append(f'{path}: {n} matched arrays with changed lengths (totals may cancel)')
        section('Collection-size differences', rows if len(rows) > 1 else [], 16)
        section('Important unchanged counts', unchanged)
        raw = []
        for path, count in sorted(self.raw_fields.items()):
            raw.append(f'{path}: {count} raw bookkeeping/order differences; semantic identities unchanged')
            raw.extend('  ' + example for example in self.raw_examples[path])
        section('Raw ring/identity metadata differences', raw)
        section('Numerical noise below rel/abs 1e-12 (not counted as changes)',
                [f'{path}: {count} values' for path, count in sorted(self.noise.items())])
        section('Internal consistency warnings (independent of cross-snapshot differences)',
                [f'{msg}' + (f' ({n} occurrences)' if n > 1 else '') for msg, n in sorted(self.warnings.items())])
        priorities = ('formal_charge', 'partial_charge', 'order', 'type', 'ringParticipation', 'enabled', 'mapping_valid')
        def rank(path):
            field = path.rsplit('.', 1)[-1]
            return (0 if field in priorities else 1 if path.startswith('$.system_model.') else 2, path)
        field_rows = []
        changed = sorted((p for p in self.fields if '[]' in p), key=rank)
        for path in (changed if details else changed[:12]):
            field_rows.append(f'{path}: {self.fields[path]} affected records')
            field_rows.extend('  ' + ex for ex in self.examples[path])
        if not details and len(changed) > 12:
            field_rows.append(f'... {len(changed)-12} more field paths; use --details.')
        section('Detailed field differences (float noise threshold: rel/abs 1e-12)', field_rows, len(field_rows))
        return '\n'.join(lines) + '\n'


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('snapshot_A')
    parser.add_argument('snapshot_B')
    parser.add_argument('--details', action='store_true')
    parser.add_argument('--max-examples', type=int, default=3)
    parser.add_argument('--coordinate-tolerance', type=float, default=1e-6)
    args = parser.parse_args(argv)
    if args.max_examples < 0:
        parser.error('--max-examples must be nonnegative')
    if not math.isfinite(args.coordinate_tolerance) or args.coordinate_tolerance < 0:
        parser.error('--coordinate-tolerance must be finite and nonnegative')
    try:
        comparison = Comparison(load_snapshot(args.snapshot_A), load_snapshot(args.snapshot_B),
                                args.max_examples, args.coordinate_tolerance)
        print(comparison.report((args.snapshot_A, args.snapshot_B), args.details), end='')
    except (ValueError, RecursionError) as exc:
        parser.exit(2, f'error: {exc}\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
