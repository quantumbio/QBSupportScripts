"""Focused stdlib tests; run: python3 -m unittest discover -s python -p 'test_model_comparisons.py'."""
import copy
import gzip
import json
from pathlib import Path
import tempfile
import unittest
import subprocess
import sys

from model_comparisons import Comparison, SCHEMA, load_snapshot


def snapshot():
    return {'snapshot': {'schema': SCHEMA, 'schema_version': 1},
            'system_model': {'name': 'synthetic', 'entryID': '', 'role': 0,
                             'counts': {'atoms': 2, 'residues': 1, 'bonds': 2},
                             'atoms': [{'index': 0, 'coordinates': [0., 0., 0.], 'formal_charge': 0},
                                       {'index': 1, 'coordinates': [1., 0., 0.], 'formal_charge': 0}],
                             'residues': [{'index': 0}],
                             'bonds': [{'index': 0, 'atomA_index': 0, 'atomB_index': 1, 'order': 1, 'type': 2},
                                       {'index': 1, 'atomA_index': 1, 'atomB_index': 2, 'order': 1, 'type': 2}]}}


class ComparisonTests(unittest.TestCase):
    def test_identical(self):
        d = snapshot()
        c = Comparison(d, copy.deepcopy(d))
        self.assertFalse(c.fields)
        self.assertFalse(c.warnings)
        self.assertIn('0 changed field paths', c.report())

    def test_counts_and_inconsistency(self):
        a, b = snapshot(), snapshot()
        b['system_model']['counts']['bonds'] = 4
        c = Comparison(a, b)
        self.assertEqual(c.fields['$.system_model.counts.bonds'], 1)
        self.assertTrue(any('actual length=2' in w for w in c.warnings))
        self.assertIn('+2', c.report())

    def test_missing_collection(self):
        a, b = snapshot(), snapshot()
        del b['system_model']['bonds']
        c = Comparison(a, b)
        self.assertIn('<missing>', c.report())
        self.assertFalse(c.warnings)  # Absent count-only output is not a proven inconsistency.

    def test_reordered_bonds_and_orientation(self):
        a, b = snapshot(), snapshot()
        b['system_model']['bonds'].reverse()
        for i, r in enumerate(b['system_model']['bonds']):
            r['index'] = i
            r['atomA_index'], r['atomB_index'] = r['atomB_index'], r['atomA_index']
        c = Comparison(a, b)
        self.assertFalse(c.fields)
        s = c.collections['$.system_model.bonds']
        self.assertEqual(s['moved'], 2)
        self.assertEqual(s['orientation'], 2)
        self.assertEqual(s['modified'], 0)
        self.assertIn('endpoint-pair sets are identical', c.report())

    def test_bond_properties_by_endpoints(self):
        a, b = snapshot(), snapshot()
        b['system_model']['bonds'].reverse()
        b['system_model']['bonds'][0]['order'] = 2
        b['system_model']['bonds'][0]['type'] = 3
        c = Comparison(a, b, max_examples=1)
        self.assertEqual(c.fields['$.system_model.bonds[].order'], 1)
        self.assertEqual(c.fields['$.system_model.bonds[].type'], 1)
        self.assertIn('(1, 2)', c.examples['$.system_model.bonds[].order'][0])

    def test_float_and_coordinate_tolerances(self):
        a, b = snapshot(), snapshot()
        a['system_model']['atoms'][0]['partial_charge'] = 0.123
        b['system_model']['atoms'][0]['partial_charge'] = 0.123 + 1e-14
        b['system_model']['atoms'][0]['coordinates'][0] = 5e-7
        c = Comparison(a, b)
        self.assertFalse(c.fields)
        stats = c.coords['$.system_model.atoms[].coordinates']
        self.assertEqual(stats[3], 0)
        self.assertAlmostEqual(stats[2], 5e-7)
        self.assertAlmostEqual((stats[1]/stats[0])**0.5, 5e-7 / 2**0.5)
        self.assertEqual(Comparison(a, b, coordinate_tolerance=1e-8).coords[
            '$.system_model.atoms[].coordinates'][3], 1)
        b['system_model']['atoms'][0]['formal_charge'] = True
        self.assertEqual(Comparison(a, b).fields['$.system_model.atoms[].formal_charge'], 1)

    def test_duplicate_keys_are_not_overwritten(self):
        a, b = snapshot(), snapshot()
        b['system_model']['bonds'][1].update(atomA_index=0, atomB_index=1)
        c = Comparison(a, b)
        self.assertTrue(any('duplicate identity keys' in w for w in c.warnings))
        self.assertTrue(c.collections['$.system_model.bonds']['positional'])

    def test_ring_membership_matching(self):
        a, b = snapshot(), snapshot()
        a['system_model']['rings'] = [{'index': 0, 'atom_indices': [1, 2, 3], 'exocyclic_double_bond': 0},
                                      {'index': 1, 'atom_indices': [4, 5, 6], 'exocyclic_double_bond': 0}]
        b['system_model']['rings'] = copy.deepcopy(a['system_model']['rings'][::-1])
        b['system_model']['rings'][0]['atom_indices'].reverse()
        c = Comparison(a, b)
        self.assertEqual(c.collections['$.system_model.rings']['moved'], 2)
        self.assertEqual(c.collections['$.system_model.rings']['modified'], 0)
        self.assertEqual(c.collections['$.system_model.rings']['membership_order'], 1)

    def test_positional_and_nested_sizes(self):
        a, b = snapshot(), snapshot()
        a['system_model']['unknown'] = [{'items': [1, 2]}, {'items': []}]
        b['system_model']['unknown'] = [{'items': [1]}, {'items': [2]}]
        c = Comparison(a, b)
        self.assertEqual(c.fields['$.system_model.unknown[].items'], 2)
        self.assertIn('positional comparison', c.report())
        self.assertIn('2 matched arrays with changed lengths', c.report(details=True))

    def test_graph_matching_and_exact_identifiers(self):
        a, b = snapshot(), snapshot()
        graph = {'atomic_graph': {'nodes': [{'atom_index': 0, 'enabled': True},
                                           {'atom_index': 1, 'enabled': True}],
                                 'edges': [{'index': 0, 'atom_a_index': 0, 'atom_b_index': 1,
                                            'bond_index': 0, 'enabled': True}]}}
        a['constitutional_graph'] = graph
        b['constitutional_graph'] = copy.deepcopy(graph)
        b['constitutional_graph']['atomic_graph']['nodes'].reverse()
        edge = b['constitutional_graph']['atomic_graph']['edges'][0]
        edge.update(atom_a_index=1, atom_b_index=0, bond_index=1, enabled=False)
        c = Comparison(a, b)
        self.assertEqual(c.collections['$.constitutional_graph.atomic_graph.nodes']['moved'], 2)
        self.assertEqual(c.fields['$.constitutional_graph.atomic_graph.edges[].bond_index'], 1)
        self.assertEqual(c.fields['$.constitutional_graph.atomic_graph.edges[].enabled'], 1)

    def test_added_removed_bonds_and_capped_examples(self):
        a, b = snapshot(), snapshot()
        b['system_model']['bonds'][0].update(atomA_index=7, atomB_index=8)
        c = Comparison(a, b, max_examples=0)
        self.assertEqual(c.collections['$.system_model.bonds']['added'], 1)
        self.assertEqual(c.collections['$.system_model.bonds']['removed'], 1)
        self.assertNotIn('endpoint-pair sets are identical', c.report())
        self.assertTrue(all(not examples for examples in c.examples.values()))

    def test_io_and_invalid_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name in ('a.json', 'a.json.gz'):
                path = Path(tmp) / name
                opener = gzip.open if name.endswith('.gz') else open
                with opener(path, 'wt') as stream:
                    json.dump(snapshot(), stream)
                self.assertEqual(load_snapshot(path), snapshot())
            path = Path(tmp) / 'bad.json'
            path.write_text('{')
            with self.assertRaisesRegex(ValueError, 'Cannot read snapshot'):
                load_snapshot(path)
            bad = snapshot()
            bad['snapshot']['schema_version'] = 99
            path.write_text(json.dumps(bad))
            with self.assertRaisesRegex(ValueError, 'unsupported'):
                load_snapshot(path)
            bad = snapshot()
            bad['system_model']['bonds'] = {}
            with self.assertRaisesRegex(ValueError, 'must be an array'):
                Comparison(snapshot(), bad)


def current_snapshot():
    """Writer-shaped v3, with a four-edge perimeter and independent policies."""
    model = {'atoms': [{'index': i, 'name': 'C' + str(i), 'coordinates': [float(i), 0., 0.],
                        'partial_charge': 0., 'ringParticipation': [0], 'ring_count': 1}
                       for i in range(4)],
             'residues': [{'index': 0, 'systemAtomIndices': [0, 1, 2, 3], 'systemRingIndices': [0]}],
             'bonds': [{'index': i, 'atomA_index': a, 'atomB_index': b, 'order': 1,
                        'ringParticipation': [0]} for i, (a, b) in enumerate([(0, 1), (1, 2), (2, 3), (3, 0)])],
             'rings': [{'index': 0, 'residue_index': 0, 'atom_indices': [0, 1, 2, 3],
                        'exocyclic_double_bond': 0}],
             'counts': {'atoms': 4, 'residues': 1, 'bonds': 4, 'rings': 1},
             'nonbonded_exclusions': {'provenance': 'imported', 'complete': True,
                                      'pairs': [{'atomA_index': 0, 'atomB_index': 1},
                                                {'atomA_index': 0, 'atomB_index': 2}]},
             'special_pairs': {'provenance': 'imported', 'complete': True,
                               'scaling': 'independent dimensionless energy multipliers',
                               'lj_coefficients': 'unscaled [A,B]',
                               'pairs': [{'atomA_index': 0, 'atomB_index': 2,
                                          'coulombMultiplier': 1./1.2, 'ljMultiplier': 0.5,
                                          'ljCoefficients': [12., 6.]}]}}
    return {'snapshot': {'schema': SCHEMA, 'schema_version': 3}, 'system_model': model}


class SchemaEvolutionTests(unittest.TestCase):
    def test_current_self_and_collection_counts(self):
        a = current_snapshot()
        c = Comparison(a, copy.deepcopy(a))
        self.assertFalse(c.fields)
        self.assertFalse(c.warnings)
        self.assertIn('special_pairs.pairs', c.report())
        self.assertIn('nonbonded_exclusions.pairs', c.report())
        self.assertIn('0/4 atoms serialized', c.report())

    def test_pair_record_order_and_endpoint_orientation(self):
        a, b = current_snapshot(), current_snapshot()
        for name in ('nonbonded_exclusions', 'special_pairs'):
            records = b['system_model'][name]['pairs']
            records.reverse()
            for pair in records:
                pair['atomA_index'], pair['atomB_index'] = pair['atomB_index'], pair['atomA_index']
        c = Comparison(a, b)
        self.assertFalse(c.fields)
        self.assertTrue(all(not s['modified'] for s in c.collections.values()))
        self.assertEqual(c.collections['$.system_model.nonbonded_exclusions.pairs']['orientation'], 2)

    def test_pair_add_remove_and_details(self):
        a, b = current_snapshot(), current_snapshot()
        for name in ('nonbonded_exclusions', 'special_pairs'):
            b['system_model'][name]['pairs'][0].update(atomA_index=1, atomB_index=3)
        c = Comparison(a, b)
        for name in ('nonbonded_exclusions', 'special_pairs'):
            path = '$.system_model.' + name + '.pairs'
            self.assertEqual(c.collections[path]['added'], 1)
            self.assertEqual(c.collections[path]['removed'], 1)
            self.assertIn('identity=(1, 3)', c.report(details=True))

    def test_special_scaling_coefficients_null_and_noise(self):
        for field, value in [('coulombMultiplier', 0.6), ('ljMultiplier', 0.25),
                             ('ljCoefficients', [6., 12.]), ('ljCoefficients', None)]:
            with self.subTest(field=field, value=value):
                a, b = current_snapshot(), current_snapshot()
                b['system_model']['special_pairs']['pairs'][0][field] = value
                c = Comparison(a, b)
                self.assertEqual(c.fields['$.system_model.special_pairs.pairs[].' + field], 1)
        a, b = current_snapshot(), current_snapshot()
        b['system_model']['special_pairs']['pairs'][0]['coulombMultiplier'] += 1e-14
        c = Comparison(a, b)
        self.assertFalse(c.fields)
        self.assertTrue(c.noise)
        b['system_model']['special_pairs']['pairs'][0]['ljCoefficients'][0] += 1e-14
        self.assertFalse(Comparison(a, b).fields)

    def test_policy_provenance_complete_and_authoritative_empty(self):
        for name in ('nonbonded_exclusions', 'special_pairs'):
            a, b = current_snapshot(), current_snapshot()
            a['system_model'][name].update(pairs=[], provenance='unavailable', complete=False)
            b['system_model'][name].update(pairs=[], provenance='imported', complete=True)
            c = Comparison(a, b)
            self.assertEqual(c.fields['$.system_model.' + name + '.provenance'], 1)
            self.assertEqual(c.fields['$.system_model.' + name + '.complete'], 1)
            self.assertEqual(c.collections['$.system_model.' + name + '.pairs']['modified'], 0)
            self.assertIn('pairs=0', c.report())

    def test_old_versions_and_charge_provenance_are_unavailable(self):
        for version in (1, 2):
            a, b = current_snapshot(), current_snapshot()
            a['snapshot']['schema_version'] = version
            for name in ('nonbonded_exclusions', 'special_pairs'):
                del a['system_model'][name]
            b['system_model']['atoms'][0]['partial_charge_imported'] = False
            c = Comparison(a, b)
            self.assertEqual(c.fields['$.snapshot.schema_version'], 1)
            self.assertEqual(c.fields['$.system_model.atoms[].partial_charge_imported'], 1)
            self.assertIn('unavailable (not serialized)', c.report())
            self.assertIn('1/4 atoms serialized', c.report())
            self.assertIn('<missing>', c.report(details=True))

    def test_duplicate_and_invalid_pairs_stay_visible(self):
        for name in ('nonbonded_exclusions', 'special_pairs'):
            a = current_snapshot()
            pairs = a['system_model'][name]['pairs']
            duplicate = copy.deepcopy(pairs[0])
            duplicate['atomA_index'], duplicate['atomB_index'] = duplicate['atomB_index'], duplicate['atomA_index']
            pairs.append(duplicate)
            c = Comparison(a, copy.deepcopy(a))
            self.assertTrue(any('duplicate identity keys' in w for w in c.warnings))
            older = copy.deepcopy(a)
            del older['system_model'][name]
            self.assertTrue(any('duplicate identity keys' in w for w in Comparison(older, a).warnings))
            for bad in (-1, 8, True, '1', 0):
                b = current_snapshot()
                b['system_model'][name]['pairs'][0]['atomB_index'] = bad
                c = Comparison(b, copy.deepcopy(b))
                self.assertTrue(any('invalid pair record' in w for w in c.warnings))
        for field, value in [('ljMultiplier', -1), ('coulombMultiplier', 'NaN'),
                             ('ljCoefficients', [1]), ('ljCoefficients', [True, 2.])]:
            a = current_snapshot()
            a['system_model']['special_pairs']['pairs'][0][field] = value
            self.assertTrue(any('invalid pair record' in w for w in Comparison(a, a).warnings))

    def test_unsupported_structures_are_explicit(self):
        for field, value, message in [('pairs', {}, 'must be an array'),
                                      ('complete', 1, 'must be a boolean'),
                                      ('provenance', 'mystery', 'unsupported')]:
            a = current_snapshot()
            a['system_model']['special_pairs'][field] = value
            with self.assertRaisesRegex(ValueError, message):
                Comparison(a, a)
        a = current_snapshot()
        a['system_model']['atoms'][0]['partial_charge_imported'] = 0
        with self.assertRaisesRegex(ValueError, 'partial_charge_imported must be a boolean'):
            Comparison(a, a)

    def test_ring_rotation_reversal_ids_and_membership_order(self):
        a, b = current_snapshot(), current_snapshot()
        b['system_model']['rings'][0].update(index=7, atom_indices=[2, 1, 0, 3])
        for name in ('atoms', 'bonds'):
            for record in b['system_model'][name]:
                record['ringParticipation'] = [7]
        b['system_model']['residues'][0]['systemRingIndices'] = [7]
        c = Comparison(a, b)
        self.assertFalse(c.fields)
        self.assertFalse(c.warnings)
        self.assertTrue(c.raw_fields)
        self.assertIn('Raw ring/identity metadata differences', c.report())

    def test_different_perimeters_with_same_atoms(self):
        a, b = current_snapshot(), current_snapshot()
        b['system_model']['rings'][0]['atom_indices'] = [0, 2, 1, 3]
        c = Comparison(a, b)
        self.assertEqual(c.collections['$.system_model.rings']['added'], 1)
        self.assertEqual(c.collections['$.system_model.rings']['removed'], 1)
        self.assertTrue(c.fields['$.system_model.atoms[].ringParticipation'])
        self.assertTrue(any('perimeter' in warning for warning in c.warnings))
        a['system_model']['rings'].append(copy.deepcopy(b['system_model']['rings'][0]))
        a['system_model']['rings'][1]['index'] = 1
        c = Comparison(a, a)
        self.assertFalse(any('duplicate identity keys' in warning for warning in c.warnings))

    def test_raw_ring_defects_and_directional_bonds(self):
        a = current_snapshot()
        a['system_model']['atoms'][0]['ringParticipation'] = [0, 0, 99]
        c = Comparison(a, copy.deepcopy(a))
        for defect in ('duplicate ring references', 'invalid ring references', 'ring_count disagrees'):
            self.assertTrue(any(defect in warning for warning in c.warnings))
        a, b = current_snapshot(), current_snapshot()
        for model in (a['system_model'], b['system_model']):
            model['bonds'][0]['symmetry_x'] = 1
        b['system_model']['bonds'][0].update(atomA_index=1, atomB_index=0)
        self.assertTrue(Comparison(a, b).fields['$.system_model.bonds[].atomA_index'])

    def test_fused_ring_remapping_and_shared_edge_membership_order(self):
        a = current_snapshot()
        model = a['system_model']
        model['rings'] = [{'index': i, 'atom_indices': ring, 'residue_index': 0}
                          for i, ring in enumerate(([0, 1, 2], [0, 1, 3]))]
        model['residues'][0]['systemRingIndices'] = [0, 1]
        for i, refs in enumerate(([0, 1], [0, 1], [0], [1])):
            model['atoms'][i].update(ringParticipation=refs, ring_count=len(refs))
        model['bonds'] = [{'index': i, 'atomA_index': x, 'atomB_index': y,
                           'ringParticipation': refs} for i, (x, y, refs) in enumerate(
                               [(0, 1, [0, 1]), (1, 2, [0]), (2, 0, [0]), (1, 3, [1]), (3, 0, [1])])]
        model['counts'].update(bonds=5, rings=2)
        b = copy.deepcopy(a)
        b['system_model']['rings'].reverse()
        for i, ring in enumerate(b['system_model']['rings']):
            ring['index'] = i
            ring['atom_indices'] = ring['atom_indices'][1:] + ring['atom_indices'][:1]
        for name in ('atoms', 'bonds', 'residues'):
            field = 'systemRingIndices' if name == 'residues' else 'ringParticipation'
            for record in b['system_model'][name]:
                record[field] = [1 - r for r in record[field]]
        c = Comparison(a, b)
        self.assertFalse(c.fields)
        self.assertFalse(c.warnings)
        self.assertEqual(c.collections['$.system_model.rings']['moved'], 2)
        self.assertTrue(c.raw_fields['$.system_model.bonds[].ringParticipation'])
        # Reordering memberships alone must not change either ring or shared edge.
        b = copy.deepcopy(a)
        b['system_model']['bonds'][0]['ringParticipation'].reverse()
        c = Comparison(a, b)
        self.assertFalse(c.fields)
        self.assertFalse(c.warnings)
        self.assertEqual(c.raw_fields['$.system_model.bonds[].ringParticipation'], 1)

    def test_new_fields_and_incomplete_older_policy_are_not_skipped(self):
        a, b = current_snapshot(), current_snapshot()
        del a['system_model']['special_pairs']['complete']
        b['system_model']['special_pairs']['future_policy'] = {'enabled': True}
        b['system_model']['special_pairs']['pairs'][0]['future_value'] = 2.
        c = Comparison(a, b)
        for field in ('complete', 'future_policy', 'pairs[].future_value'):
            self.assertTrue(c.fields['$.system_model.special_pairs.' + field])
        self.assertIn('complete=<missing>', c.report())

    def test_atom_index_correspondence_and_gzip_cli(self):
        a, b = current_snapshot(), current_snapshot()
        b['system_model']['atoms'].reverse()
        self.assertFalse(Comparison(a, b).fields)
        with tempfile.TemporaryDirectory() as tmp:
            left, right = Path(tmp) / 'a.json', Path(tmp) / 'b.json.gz'
            left.write_text(json.dumps(a))
            with gzip.open(right, 'wt') as stream:
                json.dump(b, stream)
            self.assertEqual(load_snapshot(right), b)
            result = subprocess.run([sys.executable, str(Path(__file__).with_name('model_comparisons.py')),
                                     str(left), str(right), '--details', '--max-examples', '1',
                                     '--coordinate-tolerance', '1e-8'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('0 changed field paths', result.stdout)
            self.assertIn('special_pairs.pairs', result.stdout)


if __name__ == '__main__':
    unittest.main()
