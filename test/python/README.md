# qmechanic CLI regressions

This directory contains maintained regression/integration tests. Operational
scripts belong in `QBSupportScripts/python/`.

`amber_publication_cli_tests.py` uses Python's standard-library `unittest` and
launches the actual executable in isolated temporary directories. It requires a
built qmechanic, installed force-field data, and the executable's normal runtime
library environment. No developer checkout paths are embedded in the test.

From any directory:

```sh
export QMECHANIC_TEST_BINARY=/path/to/qmechanic.bin
export QMECHANIC_HOME=/path/to/installed/qmechanic/data
python3 /path/to/QBSupportScripts/test/python/amber_publication_cli_tests.py -v
```

Or from this directory, use discovery:

```sh
python3 -m unittest discover -s . -p '*_tests.py' -v
```

With no `QMECHANIC_TEST_BINARY`, executable tests are explicitly skipped.
A specified relative executable path is resolved before changing directories;
a binary name on PATH is also accepted. Fixtures resolve relative to this file.
`QB_*`, `DEV_*`, and inherited `LD_PRELOAD` settings are cleared in child processes
so unrelated developer settings do not change the tested behavior.

The suite covers explicit 14SB/19SB preparation and waterboxing, implicitly
selected MM 19SB with explicit QM `-h pm6`, automatic `inpcrd` companion topology,
and periodic restart/import → publication. The no-`-h` preparation case records
its existing earlier CLI rejection; it does not claim no-`-h` preparation succeeds.

The optional edit callback test requires Linux and
`QMECHANIC_EDIT_PUBLICATION_CAPTURE=/path/to/amber_edit_publication_capture.so`,
built with `make CONF=Release amber-edit-publication-capture` in QMechanic. This fixture
bypasses ComplexingSeries reconstruction; it exercises main's actual callback
and actual writer. Without the helper that test is skipped. Full edit
reconstruction currently loses standard residue names and is **not covered**
by this prepared-model handoff test.

`fixtures/3fad-asn-fragment.pdb` is the deposited heavy-atom MET1/ASN2/ILE3 fragment
(24 atoms, alternate A). ASN CA distinguishes 14SB CX from 19SB XC. Adjacent
residues preserve the preparation context. It contains no BME; no deferred
sulfur/protonation chemistry is asserted. Only this CLI suite needs the fixture.

Direct writer/API empty-identity and selected-table coefficient/diagnostic
coverage remains in QMechanic's C++ doctests (`make CONF=Release amber-tests`). The CLI suite
cannot distinguish a main variable from a getter returning the same string;
source review establishes that main explicitly forwards `forcefieldName`.
It also does not provide a full historical 14SB serialized golden baseline.
