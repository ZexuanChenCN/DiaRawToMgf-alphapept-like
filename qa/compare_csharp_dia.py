"""Development QA only; the shipping CLI is pure .NET."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from clean_pfind_mgf import records

p=argparse.ArgumentParser()
p.add_argument('--reference',required=True)
p.add_argument('--candidate',required=True)
p.add_argument('--mode',choices=['fixed-features','native','alphapept-port'],required=True)
p.add_argument('--reference-features')
p.add_argument('--reference-curves')
p.add_argument('--output',required=True)
a=p.parse_args()
ref,cand=Path(a.reference),Path(a.candidate)
r=pd.read_csv(ref/'spectrum_feature_mapping.tsv',sep='\t')
c=pd.read_csv(cand/'mapping.tsv',sep='\t')
result=dict(mode=a.mode,reference_records=len(r),candidate_records=len(c),
            reference_ms2=int(r.raw_scan.nunique()),candidate_ms2=int(c.raw_scan.nunique()),
            common_ms2=len(set(r.raw_scan)&set(c.raw_scan)))
if a.mode!='native':
    x=r.merge(c,on=['raw_scan','feature_id'],suffixes=('_ref','_cs'))
    result.update(common_pairs=len(x),reference_only_pairs=len(r)-len(x),candidate_only_pairs=len(c)-len(x),
                  identical_charge=bool((x.charge_ref==x.charge_cs).all()),
                  identical_rank=bool((x.candidate_rank_ref==x.candidate_rank_cs).all()),
                  identical_title=bool((x.title_ref==x.title_cs).all()),
                  max_precursor_difference=float(abs(x.precursor_mz_ref-x.precursor_mz_cs).max()),
                  max_local_intensity_difference=float(abs(x.local_intensity_ref-x.local_intensity_cs).max()))
else:
    rg={(int(scan),int(z)):np.sort(d.precursor_mz.to_numpy()) for (scan,z),d in r.groupby(['raw_scan','charge'])}
    matched=0
    for row in c.itertuples():
        values=rg.get((int(row.raw_scan),int(row.charge)))
        if values is None:continue
        i=np.searchsorted(values,row.precursor_mz)
        if any(abs(values[j]-row.precursor_mz)<=row.precursor_mz*10e-6 for j in (i-1,i) if 0<=j<len(values)):
            matched+=1
    result['candidate_pairs_with_same_charge_reference_within_10ppm']=matched
    result['candidate_pair_match_fraction']=matched/len(c)
    result['native_feature_algorithm_is_not_equivalent']=True

def peak_signatures(directory,ms2_column):
    path=next(directory.glob('*.mgf'))
    seen={}
    for h,peaks in records(path):
        scan=int(h.get('SCANS',h['TITLE'].rsplit('.',5)[1]))
        # Compare formatted peak values, independent of ordering of precursor candidates.
        signature=hashlib.sha256(np.round(peaks,6).astype('<f8').tobytes()).hexdigest()
        if scan in seen:
            assert seen[scan]==signature
        seen[scan]=signature
    return seen

if a.mode!='native':
    print('Comparing full MGF peak arrays...',flush=True)
    rs=peak_signatures(ref,'raw_scan');cs=peak_signatures(cand,'raw_scan')
    result['common_scan_peak_arrays_identical']=all(rs[k]==cs[k] for k in rs.keys()&cs.keys())
    result['identical_mgf_scan_set']=rs.keys()==cs.keys()
    result['association_alignment_pass']=bool(result['reference_only_pairs']==0 and result['candidate_only_pairs']==0 and
        result['identical_charge'] and result['identical_rank'] and result['identical_title'] and result['max_precursor_difference']<1e-7 and
        result['common_scan_peak_arrays_identical'] and result['identical_mgf_scan_set'])
if a.reference_features:
    rf=pd.read_csv(a.reference_features,sep='\t').set_index('feature_id')
    cf=pd.read_csv(cand/'ms1_features.tsv',sep='\t').set_index('feature_id').rename(columns={'intensity_apex':'ms1_int_sum_apex','area_min':'ms1_int_sum_area'})
    common=rf.index.intersection(cf.index)
    result['features']=dict(reference=len(rf),candidate=len(cf),same_ids=rf.index.equals(cf.index),columns={})
    for name in cf.columns.intersection(rf.columns):
        expected=rf.loc[common,name].to_numpy();actual=cf.loc[common,name].to_numpy()
        result['features']['columns'][name]=dict(max_abs_error=float(np.max(np.abs(expected-actual))),
            agreement_at_rtol_1e_12=float(np.mean(np.isclose(expected,actual,rtol=1e-12,atol=1e-8))))
if a.reference_curves:
    import h5py
    print('Comparing full MS1 chromatograms...',flush=True)
    with h5py.File(a.reference_curves) as h:
        offsets=h['offsets'][:];scan=h['ms1_scan_num'][:];fid=h['feature_id'][:];intensity=h['isotope_sum_intensity'][:]
    expected=pd.DataFrame(dict(feature_id=fid,ms1_scan=np.repeat(scan,np.diff(offsets)),expected=intensity))
    actual=pd.read_csv(cand/'precursor_curves.tsv',sep='\t',usecols=['feature_id','ms1_scan','isotope_sum_intensity'])
    merged=expected.merge(actual,on=['feature_id','ms1_scan'])
    result['curves']=dict(reference_points=len(expected),candidate_points=len(actual),common_points=len(merged),
        max_abs_error=float(np.max(np.abs(merged.expected-merged.isotope_sum_intensity))),
        allclose=bool(np.allclose(merged.expected,merged.isotope_sum_intensity,rtol=1e-12,atol=1e-8)))
if a.mode=='alphapept-port':
    result['port_alignment_pass']=result['association_alignment_pass']
    if 'features' in result:
        f=result['features'];result['port_alignment_pass'] &= f['same_ids'] and all(x['agreement_at_rtol_1e_12']==1 for x in f['columns'].values())
    if 'curves' in result:
        c=result['curves'];result['port_alignment_pass'] &= c['reference_points']==c['candidate_points']==c['common_points'] and c['allclose']
Path(a.output).write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2),flush=True)
if a.mode=='alphapept-port' and not result['port_alignment_pass']:raise SystemExit(1)
