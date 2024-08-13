#!/bin/bash

# MT-threaded spot test script


if [[ -z "${QBHOME}" ]]; then
  echo "Set QBHOME or source /path/to/DivConSuite/etc/qbenv.sh before running this script."
  exit 1
fi

DIVCON=${QBHOME}/bin/qmechanic
WORKDIR=${PWD}

# Density-driven docking: known pocket location 
#   * This is effectively re-docking of the ligand in density.
#   * We use the same appraoch we have always used in which the placed ligand defines the pocket.
#   * Since we know where we are docking, the "blobs" keyword is unnecessary.

echo "Tutorial #1a-pdb: Automated DivCon-based Docking – using qmechanic and XModeScore"
dir4=XModeScore-dock_ONIOM
rm -rf ${WORKDIR}/$dir4 ; mkdir -p ${WORKDIR}/$dir4 ; cd ${WORKDIR}/$dir4
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.mtz
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4XF_A_402.pdb

${DIVCON} 5C3K-H_refine_001.pdb 5C3K-H_refine_001.mtz   `# input files. Since MTZ file is provided, it is understood that density docking will be used`                                                                         \
    -h amberff14sb pm6                                  `# Use AMBERFF14SB for all MM and PM6 for all QM. For v1: all steps up to XModeScore should be done with MM and final XModeScore should be done with QM/MM`             \
    --ligand 4XF_A_402.pdb                              `# This is a placed ligand so no novel compound required. The novel compound (if provided) is the LIGAND in any subsequent steps which uses the LIGAND keyword.`        \
    --qm-region LIGAND,/A/ZN/401// 3.0 0                `# The QM region in ONIOM. Note: for v1, all steps up to XModeScore should be done with MM and final XModeScore should be done with QM/MM. `                            \
    --confSearch LIGAND opt off                         `# Use the input conformer (LIGAND keyword which references --ligand input), keep poses/enumerations w/in top 1 standard deviation, perform no opt on conformer`                \
    --flip all                                          `# Enumerate all rotomers (flip states) available in LIGAND`                                                                                                            \
    --chirality all                                     `# Enumerate all stereoisomers (chiral centers) available in LIGAND`                                                                                                    \
    --dock 5,1SD opt torsion POCKET 100 0.01            `# Dock the top 5 poses/enumerations from confSearch and refine/optimize the entire site after placement (using MM), keep the top 1 standard deviation of refined sets.`    \
    --protomers all [-1..1]                             `# Enumerate all protomers available in LIGAND w/in -1..1 and be sure to mirror protomer modifications to the target.`                                                  \
    --xmodescore opt all                                `# Use XModeScore to decide the final "winners" and "losers." If --qm-region and -h amberff14sb pm6 set, use QM/MM refinement of the target+ligand for XModeScore.`     \
    --np 8  -v 2                                        `# All refinements and all scoring happen in parallel. Verbosity=2 should be clean. Verbosity=3 should just have a bit more debug information but should still be fairly clean.`    \
    -p 5C3K-out.pdb 5C3K-out.mtz                        `# There should be a PDB and MTZ file for each refined/scored case from XModeScore. Should be called: 5C3K-out_pose-1.pdb 5C3K-out_pose-1.mtz, 5C3K-out_pose-2.pdb 5C3K-out_pose-2.mtz etc (note: the pose-# should equal the number in the Pose Scoring Table in the screen output)`   \
        >& OUT.SCREEN

echo "Tutorial #1a-mol2: Automated DivCon-based Docking – using qmechanic and XModeScore"
dir4=XModeScore-dock_MM
rm -rf ${WORKDIR}/$dir4 ; mkdir -p ${WORKDIR}/$dir4 ; cd ${WORKDIR}/$dir4
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.mtz
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4XF_A_402.pdb

${DIVCON} 5C3K-H_refine_001.pdb 5C3K-H_refine_001.mtz   `# input files. Since MTZ file is provided, it is understood that density docking will be used`     \
    -h amberff14sb                                      `# Use AMBERFF14SB for all MM. Since no PM6 is requested & no region is set, all calculations (including XModeScore) should be done with MM only.`      \
    --ligand 4XF_A_402.pdb                              `# This is a placed ligand so no novel compound required. `     \
    --confSearch LIGAND opt off                         `# Use the input conformer (LIGAND keyword which references --ligand input), keep poses/enumerations w/in top 1 standard deviation, perform no opt on conformer`        \
    --flip off                                          `# DONT Enumerate all rotomers (flip states) available in LIGAND`       \
    --chirality off                                     `# DONE Enumerate all stereoisomers (chiral centers) available in LIGAND`   \
    --dock 5,1SD opt torsion LIGAND 100 0.01            `# Dock the top 5 poses/enumerations from confSearch and refine/optimize the ligand after placement, keep the top 1 standard deviation of refined sets to pass to XModeScore.` \
    --protomers all [-1..1]                             `# Enumerate all protomers available in LIGAND w/in -1..1 and be sure to mirror protomer modifications to the target`       \
    --xmodescore opt POCKET                             `# Use XModeScore to decide the final "winners" and "losers." Only optimize the POCKET+LIGAND and not the whole target+ligand.`     \
    --np 8  -v 2                                        `# All refinements and all scoring happen in parallel`      \
    -p 5C3K-out.mol2 5C3K-out.mtz                       `# There should be a MOL2 and MTZ file for each refined/scored case from XModeScore. Should be called: 5C3K-out_pose-1.mol2 5C3K-out_pose-1.mtz, 5C3K-out_pose-2.mol2 5C3K-out_pose-2.mtz etc (note: the pose-# should equal the number in the Pose Scoring Table in the screen output and the XModeScore terms should be stored in the mol2 file)` \
        >& OUT.SCREEN

echo "Tutorial #1b: Automated DivCon-based Docking – using qmechanic and EndState MTScore"
dir4=MTDock-MTScore_garf
rm -rf ${WORKDIR}/$dir4 ; mkdir -p ${WORKDIR}/$dir4 ; cd ${WORKDIR}/$dir4
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.mtz
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4XF_A_402.pdb

${DIVCON} 5C3K-H_refine_001.pdb 5C3K-H_refine_001.mtz   `# input files. Since MTZ file is provided, it is understood that density docking will be used`     \
    -h garf                                             `# This is going to be an MT calculation so as always: MM-amberff14sb used for all optimizations and GARF for all MT calculations.`     \
    --ligand 4XF_A_402.pdb                              `# This is a placed ligand so no novel compound required. `     \
    --mtcs LIGAND opt off                               `# Use the input conformer (LIGAND keyword which references --ligand input), keep poses/enumerations w/in top 1 standard deviation, perform no opt on conformer`    \
    --flip all                                          `# Enumerate all rotomers (flip states) available in LIGAND`        \
    --chirality all                                     `# Enumerate all stereoisomers (chiral centers) available in LIGAND`        \
    --mtdock 5,1SD opt torsion POCKET 100 0.01          `# Dock the top 5 poses/enumerations from above and refine/optimize the entire site after placement, keep the top 1 standard deviation of refined sets to pass to MTScore.`        \
    --protomers off                                     `# DONT Enumerate all protomers available in LIGAND w/in -1..1 and be sure to mirror protomer modifications to the target`      \
    --mtscore endstate                                  `# Use MTScore EndState to decide the final "winners" and "losers"`     \
    --np 8  -v 2                                        `# All refinements and all scoring happen in parallel`              \
    -p 5C3K-out.mol2 5C3K-out.mtz                       `# There should be a MOL2 and MTZ file for each refined/scored case from MTScore. Should be called: 5C3K-out_pose-1.mol2 5C3K-out_pose-1.mtz, 5C3K-out_pose-2.mol2 5C3K-out_pose-2.mtz etc (note: the pose-# should equal the number in the Pose Scoring Table in the screen output and the MTScore terms should be stored in the mol2 file)`       \
        >& OUT.SCREEN

echo "Tutorial #1c: Automated DivCon-based Docking – using qmechanic and Ensemble MTScore"
dir4=MTDock-MTScore_amberff14sb
rm -rf ${WORKDIR}/$dir4 ; mkdir -p ${WORKDIR}/$dir4 ; cd ${WORKDIR}/$dir4
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.mtz
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4XF_A_402.pdb

${DIVCON} 5C3K-H_refine_001.pdb 5C3K-H_refine_001.mtz   `# input files. Since MTZ file is provided, it is understood that density docking will be used`     \
    -h amberff14sb                                      `# This is going to be an MT calculation so as always: MM-amberff14sb used for all optimizations and AMBERff14sb for all MT calculations.`  \
    --ligand 4XF_A_402.pdb                              `# This is a placed ligand so no novel compound required. ` \
    --confSearch LIGAND opt off                         `# Use the input conformer (LIGAND keyword which references --ligand input), keep poses/enumerations w/in top 1 standard deviation, perform no opt on conformer`    \
    --flip all                                          `# Enumerate all rotomers (flip states) available in LIGAND`        \
    --chirality all                                     `# Enumerate all stereoisomers (chiral centers) available in LIGAND`    \
    --mtdock 5,1SD opt torsion POCKET 100 0.01          `# Dock the top 5 poses/enumerations from above and refine/optimize the entire site after placement, keep the top 1 standard deviation of refined sets.`   \
    --protomers all [-1..1]                             `# Enumerate all protomers available in LIGAND w/in -1..1 and be sure to mirror protomer modifications to the target`       \
    --mtscore ensemble                                  `# Use MTScore Ensemble across the entire ensemble of protein/ligand complexes generated.`      \
    --np 8  -v 2                                        `# All refinements and all scoring happen in parallel`      \
    -p 5C3K-out.mol2 5C3K-out.mtz                       `# There should be a MOL2 and MTZ file for each refined/scored case from MTScore. Should be called: 5C3K-out_pose-1.mol2 5C3K-out_pose-1.mtz, 5C3K-out_pose-2.mol2 5C3K-out_pose-2.mtz etc (note: the pose-# should equal the number in the Pose Scoring Table in the screen output and the MTScore terms should be stored in the mol2 file)`       \
        >& OUT.SCREEN

echo "Tutorial #1d: DivCon --prepare followed by DivCon-based Docking – using qmechanic and XModeScore"
dir4=prepare_XModeScore-dock_ONIOM
rm -rf ${WORKDIR}/$dir4 ; mkdir -p ${WORKDIR}/$dir4 ; cd ${WORKDIR}/$dir4

${DIVCON} 5C3K                                          `# Run on the original PDB file downloaded automatically from the RCSB`     \
    -h amberff14sb                                      `# Use amberff14sb throughout the calculation`      \
    --prepare all                                       `# Prepare all missing residues/sidechains/etc`     \
    -p  5C3K-refined.pdb  5C3K-refined.mtz              `# Output the final PDB and MTZ files for use in the next step`     \
    -v 2 --np 8  >& OUT.SCREEN
    
${DIVCON} 5C3K-refined.pdb  5C3K-refined.mtz -O         `# Use the prepared files for this analysis. Overwrite the previous h5 file.`       \
    -h amberff14sb pm6                                  `# Use AMBERFF14SB for all MM and PM6 for all QM. For v1: all steps up to XModeScore should be done with MM and final XModeScore should be done with QM/MM`     \
    --ligand /A/4XF/402//                               `# Run on the selection (not a supplied ligand file).`      \
    --qm-region LIGAND,/A/ZN/401// 3.0 0                `# The QM region in ONIOM. Note: for v1, all steps up to XModeScore should be done with MM and final XModeScore should be done with QM/MM. `        \
    --confSearch LIGAND opt off                         `# Use the input conformer (LIGAND keyword which references --ligand input), keep poses/enumerations w/in top 1 standard deviation, perform no opt on conformer`        \
    --flip all                                          `# Enumerate all rotomers (flip states) available in LIGAND`        \
    --chirality all                                     `# Enumerate all stereoisomers (chiral centers) available in LIGAND`        \
    --dock 5,1SD opt torsion POCKET 100 0.01            `# Dock the top 5 poses/enumerations from above and refine/optimize the entire site after placement (using MM), keep the top 1 standard deviation of refined sets.`        \
    --protomers all [-1..1]                             `# Enumerate all protomers available in LIGAND w/in -1..1 and be sure to mirror protomer modifications to the target.`      \
    --xmodescore opt all                                `# Use XModeScore to decide the final "winners" and "losers." If --qm-region and -h amberff14sb pm6 set, use QM/MM refinement of the target+ligand for XModeScore.`     \
    --np 8  -v 2                                        `# All refinements and all scoring happen in parallel. Verbosity=2 should be clean. Verbosity=3 should just have a bit more debug information but should still be fairly clean.`        \
    -p 5C3K-out.pdb 5C3K-out.mtz                        `# There should be a PDB and MTZ file for each refined/scored case from XModeScore. Should be called: 5C3K-out_pose-1.pdb 5C3K-out_pose-1.mtz, 5C3K-out_pose-2.pdb 5C3K-out_pose-2.mtz etc (note: the pose-# should equal the number in the Pose Scoring Table in the screen output)`       \
        >& OUT.SCREEN


# Density Driven Docking: blob search (unknonw location)
#   * The "blobs" keyword is used to search for 
