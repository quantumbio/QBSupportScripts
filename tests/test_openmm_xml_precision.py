import copy
import importlib.util
import json
import sys
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'python'))
from openmm_xml_precision import audit_system_pair
import summarize_MD_comparisons as summary

# Small synthetic schema fixtures. No fixture is represented as a real MD run.
def fixture_root():
    return ET.fromstring('''<System openmmVersion="8.5.2" type="System" version="1">
<PeriodicBoxVectors><A x="3" y="0" z="0"/><B x="0" y="3" z="0"/><C x="0" y="0" z="3"/></PeriodicBoxVectors>
<Particles><Particle mass="12"/><Particle mass="12"/><Particle mass="12"/><Particle mass="12"/></Particles>
<Constraints><Constraint p1="0" p2="1" d="0.15"/></Constraints>
<Forces>
<Force forceGroup="0" name="bond" type="HarmonicBondForce" usesPeriodic="0" version="2"><Bonds><Bond p1="1" p2="2" d="0.15" k="100"/></Bonds></Force>
<Force forceGroup="1" name="angle" type="HarmonicAngleForce" usesPeriodic="0" version="2"><Angles><Angle p1="0" p2="1" p3="2" a="1.9" k="100"/></Angles></Force>
<Force forceGroup="2" name="torsion" type="PeriodicTorsionForce" usesPeriodic="0" version="2"><Torsions><Torsion p1="0" p2="1" p3="2" p4="3" periodicity="3" phase="0" k="1"/></Torsions></Force>
<Force alpha="0" cutoff="1" dispersionCorrection="1" ewaldTolerance=".0001" exceptionsUsePeriodic="0" forceGroup="3" includeDirectSpace="1" ljAlpha="0" ljnx="0" ljny="0" ljnz="0" method="4" name="nb" nx="0" ny="0" nz="0" recipForceGroup="-1" rfDielectric="78.5" switchingDistance="-1" type="NonbondedForce" useSwitchingFunction="0" version="4">
<GlobalParameters/><ParticleOffsets/><ExceptionOffsets/>
<Particles><Particle q="0.1" sig="0.3" eps="0.5"/><Particle q="0.1" sig="0.3" eps="0.5"/><Particle q="0.1" sig="0.3" eps="0.5"/><Particle q="0.1" sig="0.3" eps="0.5"/></Particles>
<Exceptions><Exception p1="0" p2="1" q="0" sig="1" eps="0"/><Exception p1="0" p2="3" q="0.02" sig="0.3" eps="0.3"/></Exceptions>
</Force>
<Force forceGroup="4" frequency="1" name="cmm" type="CMMotionRemover" version="1"/>
</Forces></System>''')

def nb(root):
    return root.find("./Forces/Force[@type='NonbondedForce']")

def write_case(tmp_path, a, b):
    for label,r in [('cpp',a),('python',b)]:
        d=tmp_path/label; d.mkdir(exist_ok=True)
        ET.ElementTree(r).write(d/'openmm_system.xml')
    return tmp_path/'cpp'/'openmm_system.xml',tmp_path/'python'/'openmm_system.xml'

@pytest.mark.parametrize('field,value,index,status', [
    ('sig','0.5',0,'MATCH_WITHIN_TOLERANCE'),
    ('eps','0.300000003',1,'MATCH_WITHIN_TOLERANCE'),
    ('sig','0.3000000003',1,'MATCH_WITHIN_TOLERANCE'),
    ('q','1e-25',0,'MISMATCH'),
    ('eps','1e-25',0,'MISMATCH'),
    ('eps','0',1,'MISMATCH'),
    ('q','0',1,'MISMATCH'),
    ('q','-0.02',1,'MISMATCH'),
    ('sig','0.30001',1,'MISMATCH'),
    ('eps','0.301',1,'MISMATCH'),
    ('eps','nan',1,'INCOMPLETE'),
    ('eps','-0.3',1,'INCOMPLETE'),
])
def test_exception_policy(tmp_path,field,value,index,status):
    a=fixture_root(); b=copy.deepcopy(a); nb(b).find('Exceptions')[index].set(field,value)
    assert audit_system_pair(*write_case(tmp_path,a,b))['status']==status

@pytest.mark.parametrize('path,field,value', [
    ('./Particles/Particle','mass','12.01'),
    ('./Constraints/Constraint','d','0.151'),
    ("./Forces/Force[@type='NonbondedForce']",'cutoff','1.1'),
    ("./Forces/Force[@type='NonbondedForce']/Particles/Particle",'q','0.101'),
    ("./Forces/Force[@type='HarmonicAngleForce']/Angles/Angle",'a','1.901'),
    ("./Forces/Force[@type='PeriodicTorsionForce']/Torsions/Torsion",'k','1.01'),
])
def test_other_physics_not_masked(tmp_path,path,field,value):
    a=fixture_root();b=copy.deepcopy(a);b.find(path).set(field,value)
    assert audit_system_pair(*write_case(tmp_path,a,b))['status']=='MISMATCH'

@pytest.mark.parametrize('kind',['custom','offset','unknown_attribute','virtual_site','duplicate_pair'])
def test_unsupported_not_good(tmp_path,kind):
    a=fixture_root();b=copy.deepcopy(a)
    for root in (a,b):
        if kind=='custom': ET.SubElement(root.find('Forces'),'Force',type='CustomBondForce')
        elif kind=='offset': ET.SubElement(nb(root).find('ExceptionOffsets'),'Offset',exception='0',q='1')
        elif kind=='unknown_attribute': nb(root).set('newPhysicsOption','1')
        elif kind=='virtual_site': ET.SubElement(root.find('Particles')[0],'VirtualSite')
        elif kind=='duplicate_pair': nb(root).find('Exceptions').append(copy.deepcopy(nb(root).find('Exceptions')[0]))
    assert audit_system_pair(*write_case(tmp_path,a,b))['status']=='INCOMPLETE'

def test_missing_pair(tmp_path):
    a=fixture_root();b=copy.deepcopy(a);nb(b).find('Exceptions').remove(nb(b).find('Exceptions')[0])
    assert audit_system_pair(*write_case(tmp_path,a,b))['status']=='MISMATCH'

def test_duplicate_torsion_multiplicity(tmp_path):
    a=fixture_root();b=copy.deepcopy(a);t=b.find("./Forces/Force[@type='PeriodicTorsionForce']/Torsions");t.append(copy.deepcopy(t[0]))
    assert audit_system_pair(*write_case(tmp_path,a,b))['status']=='MISMATCH'

def test_pair_order_and_metadata(tmp_path):
    a=fixture_root();b=copy.deepcopy(a);b.set('openmmVersion','8.1.2');e=nb(b).find('Exceptions')[0];e.set('p1','1');e.set('p2','0')
    assert audit_system_pair(*write_case(tmp_path,a,b))['status']=='MATCH_WITHIN_TOLERANCE'

def legacy_summary(count=1,version=4,problems=None):
    if problems is None: problems=['NonbondedForce[0] differs in 1 physics-affecting section(s)']
    return {'label1':'C++','label2':'Python','openmm_validation':{'serialized_xml':{
        'comparison_version':version, 'system':{'fully_compared':True,'physics_match':False,'problems':problems,
        'forces':{'NonbondedForce':[{'exception_diagnostic':{'physics_affecting_parameter_mismatch_count':count}}]}}}}}

def setup_summary(tmp_path, data):
    a=fixture_root();b=copy.deepcopy(a);nb(b).find('Exceptions')[1].set('eps','0.300000003')
    write_case(tmp_path,a,b)
    p=tmp_path/'test_summary.json';p.write_text(json.dumps(data));return p

def test_opt_in_default_unchanged_and_keeps_original(tmp_path):
    p=setup_summary(tmp_path,legacy_summary());before={f:f.read_bytes() for f in tmp_path.rglob('*') if f.is_file()}
    assert summary.analyze_one(p,physics_only=True)['status']=='REVIEW'
    r=summary.analyze_one(p,physics_only=True,recheck_xml_precision=True)
    assert r['status']=='GOOD'
    assert r['physics_affecting_exception_mismatch_count']==1
    assert r['xml_precision_material_exception_count']==0
    assert r['xml_precision_original_model_status']=='MISMATCH'
    assert 'MATCH_WITHIN_TOLERANCE' in r['notes']
    assert all(f.read_bytes()==v for f,v in before.items())

@pytest.mark.parametrize('data',[
    legacy_summary(count=2),
    legacy_summary(version=3),
    legacy_summary(problems=['other physics problem']),
    legacy_summary(problems=['NonbondedForce[0] differs in 2 physics-affecting section(s)']),
])
def test_no_unverified_downgrade(tmp_path,data):
    p=setup_summary(tmp_path,data)
    old=summary.analyze_one(p,physics_only=True)
    new=summary.analyze_one(p,physics_only=True,recheck_xml_precision=True)
    assert old['status']==new['status']

def test_material_stays_review(tmp_path):
    p=setup_summary(tmp_path,legacy_summary());f=tmp_path/'python'/'openmm_system.xml';r=ET.parse(f).getroot();nb(r).find('Exceptions')[1].set('eps','0');ET.ElementTree(r).write(f)
    assert summary.analyze_one(p,physics_only=True,recheck_xml_precision=True)['status']=='REVIEW'

def test_missing_xml_cannot_clear_review(tmp_path):
    p=setup_summary(tmp_path,legacy_summary());(tmp_path/'python'/'openmm_system.xml').unlink()
    assert summary.analyze_one(p,physics_only=True,recheck_xml_precision=True)['status']=='REVIEW'

def test_missing_xml_old_good_becomes_limited(tmp_path):
    data=legacy_summary(count=0,problems=[]);data['openmm_validation']['serialized_xml']['system']['physics_match']=True
    p=setup_summary(tmp_path,data);(tmp_path/'python'/'openmm_system.xml').unlink()
    assert summary.analyze_one(p,physics_only=True,recheck_xml_precision=True)['status']=='GOOD_LIMITED'
