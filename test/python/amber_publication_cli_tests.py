"""Selected force-field publication through real preparation, including waterbox."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import shutil
import tempfile
import unittest


FIXTURE = Path(__file__).resolve().parent / 'fixtures/3fad-asn-fragment.pdb'
BINARY = os.environ.get('QMECHANIC_TEST_BINARY')
if BINARY:
    # Resolve before subprocess changes cwd to its isolated temporary directory.
    BINARY = str(Path(shutil.which(BINARY) or BINARY).resolve())


@unittest.skipUnless(os.environ.get('QMECHANIC_TEST_BINARY'), 'requires built qmechanic and force-field data')
class AmberPublicationCLI(unittest.TestCase):
    def run_publication(self, directory, inputs, extra, preload=None, disable_opt=True, expected_error=None):
        env = os.environ.copy()
        for name in list(env):
            if name.startswith(('QB_', 'DEV_')) or name == 'LD_PRELOAD':
                env.pop(name)
        env.update(QB_BIOASSEMBLY='0', QB_SET_BOX_WATER_NUMBER='32', QB_SET_IONIC_STRENGTH='0')
        if preload:
            env['LD_PRELOAD'] = preload
        directory.mkdir(exist_ok=True)
        options = ['--opt', 'off'] if disable_opt else []
        result = subprocess.run(
            [BINARY, *map(str, inputs), '--np', '1', '-v', '1',
             *options, '-d', str(directory), '-O', *extra],
            cwd=directory, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120)
        if expected_error:
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(expected_error, result.stdout)
            self.assertNotIn('PARMTOP publication requires an explicit force field', result.stdout)
        else:
            self.assertEqual(result.returncode, 0, result.stdout[-6000:])
        return result.stdout

    def test_existing_no_h_prepare_gate_precedes_default_resolution(self):
        source = FIXTURE
        with tempfile.TemporaryDirectory(prefix='amber-no-h-') as tmp:
            self.run_publication(Path(tmp), [source], ['--prepare', '-p', 'parm7'], disable_opt=False,
                                 expected_error='--opt requires that at least one functional has been chosen')

    def test_implicit_mm_selection_and_automatic_companion_topology(self):
        source = FIXTURE
        for output_format in ('parm7', 'inpcrd'):
            with self.subTest(output_format=output_format), tempfile.TemporaryDirectory(prefix='amber-default-') as tmp:
                directory = Path(tmp)
                # Only the QM Hamiltonian is explicit; main selects the default MM 19SB.
                self.run_publication(directory, [source], ['--prepare', '-h', 'pm6', '-p', output_format])
                snapshot = json.loads(gzip.decompress((directory / 'before_final_io.json.gz').read_bytes()))['system_model']
                asn = next(r for r in snapshot['residues'] if r['name'] == 'ASN')
                ca = next(snapshot['atoms'][i] for i in asn['systemAtomIndices'] if snapshot['atoms'][i]['name'] == 'CA')
                self.assertEqual(ca['type'], 'XC')
                topology = next(directory.glob('*.parm7' if output_format == 'parm7' else '*.prmtop'))
                self.assertIn('XC', topology.read_text().split('%FLAG AMBER_ATOM_TYPE\n')[1].split('%FLAG')[0])
                self.assertEqual(len(list(directory.glob('*.inpcrd'))), 1)

    def test_restart_import_publication_uses_main_selection(self):
        source = FIXTURE
        with tempfile.TemporaryDirectory(prefix='amber-republish-') as tmp:
            prepared, restarted = Path(tmp) / 'prepared', Path(tmp) / 'restarted'
            self.run_publication(prepared, [source], ['--prepare', '--waterbox', '-h', 'amberff19SB', '-p', 'parm7'])
            original = next(prepared.glob('*.parm7'))
            coordinates = next(prepared.glob('*.inpcrd'))
            self.run_publication(restarted, [coordinates], ['-h', str(original), '-p', 'parm7'])
            republished = next(restarted.glob('*.parm7')).read_text()
            self.assertIn('XC', republished.split('%FLAG AMBER_ATOM_TYPE\n')[1].split('%FLAG')[0])
            pointers = republished.split('%FLAG POINTERS\n')[1].split('%FLAG')[0].splitlines()[1:]
            self.assertEqual([int(value) for line in pointers for value in line.split()][27], 1)  # IFBOX
            for section in ('SOLVENT_POINTERS', 'ATOMS_PER_MOLECULE', 'BOX_DIMENSIONS'):
                self.assertEqual(original.read_text().split('%FLAG ' + section + '\n')[1].split('%FLAG')[0],
                                 republished.split('%FLAG ' + section + '\n')[1].split('%FLAG')[0])

    @unittest.skipUnless(os.environ.get('QMECHANIC_EDIT_PUBLICATION_CAPTURE'),
                         'requires prepared-model handoff probe; real edit reconstruction loses standard names')
    def test_edit_callback_publication_uses_outer_selection(self):
        source = FIXTURE
        with tempfile.TemporaryDirectory(prefix='amber-edit-') as tmp:
            directory = Path(tmp)
            replacement = directory / 'asn.pdb'
            replacement.write_text(''.join(line + '\n' for line in source.read_text().splitlines()
                                          if line.startswith('ATOM') and line[17:20] == 'ASN') + 'END\n')
            reference = directory / 'reference.pdb'
            reference.write_text(source.read_text())
            output = directory / 'output'
            stdout = self.run_publication(output, [reference], ['--prepare', '-h', 'amberff19SB', '-p', 'parm7',
                                                 '--edit', '/A/ASN/2/,' + str(replacement)],
                                         preload=str(Path(os.environ['QMECHANIC_EDIT_PUBLICATION_CAPTURE']).resolve()))
            self.assertIn('Prepared edit-publication fixture callback', stdout)
            self.assertIn('Prepared edit-publication callback returned successfully', stdout)
            # The callback publishes; the ordinary final-publication branch is not needed.
            callback_outputs = list(output.glob('*.parm7'))
            self.assertEqual(len(callback_outputs), 1)
            topology = callback_outputs[0].read_text()
            self.assertIn('XC', topology.split('%FLAG AMBER_ATOM_TYPE\n')[1].split('%FLAG')[0])

    def test_selected_table_after_preparation_and_waterbox(self):
        source = FIXTURE
        for field, ca_type in (('amberff14SB', 'CX'), ('amberff19SB', 'XC')):
            for waterbox in (False, True):
                with self.subTest(field=field, waterbox=waterbox), tempfile.TemporaryDirectory(prefix='amber-publication-cli-') as tmp:
                    directory = Path(tmp)
                    env = os.environ.copy()
                    for name in list(env):
                        if name.startswith('QB_') or name.startswith('DEV_') or name == 'LD_PRELOAD':
                            env.pop(name)
                    env.update(QB_BIOASSEMBLY='0', QB_SET_BOX_WATER_NUMBER='32', QB_SET_IONIC_STRENGTH='0')
                    command = [BINARY, str(source.resolve()), '--np', '1', '-v', '1',
                               '--prepare', '--opt', 'off', '-h', field, '-p', 'parm7', '-d', tmp, '-O']
                    if waterbox:
                        command.append('--waterbox')
                    result = subprocess.run(command, cwd=tmp, env=env, stdout=subprocess.PIPE,
                                            stderr=subprocess.STDOUT, text=True, timeout=120)
                    self.assertEqual(result.returncode, 0, result.stdout[-4000:])
                    snapshot = json.loads(gzip.decompress((directory / 'before_final_io.json.gz').read_bytes()))['system_model']
                    residue = next(r for r in snapshot['residues'] if r['name'] == 'ASN' and r['sequence'] == 2)
                    ca = next(snapshot['atoms'][i] for i in residue['systemAtomIndices'] if snapshot['atoms'][i]['name'] == 'CA')
                    self.assertEqual(ca['type'], ca_type)
                    topology = next(directory.glob('*.parm7')).read_text()
                    self.assertIn(ca_type, topology.split('%FLAG AMBER_ATOM_TYPE\n')[1].split('%FLAG')[0])
                    if waterbox:
                        self.assertEqual(snapshot['gemmi']['crystallography']['spacegroup_hm'], 'P 1')
                        self.assertIn('%FLAG BOX_DIMENSIONS', topology)
                        self.assertIn('%FLAG ATOMS_PER_MOLECULE', topology)


if __name__ == '__main__':
    unittest.main()
