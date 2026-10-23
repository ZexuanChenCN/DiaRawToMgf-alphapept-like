"""Development reference snapshots from the pinned AlphaPept 0.5.3 implementation.
Not a runtime dependency of the C# port.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time
import h5py
import numpy as np
import alphapept.performance as perf
perf.set_worker_count(8)
perf.set_compilation_mode('numba-multithread')
import alphapept.feature_finding as ff

p=argparse.ArgumentParser()
p.add_argument('--hdf',required=True)
p.add_argument('--out',required=True)
p.add_argument('--rt-start',type=float,default=10)
p.add_argument('--rt-end',type=float,default=12)
p.add_argument('--through',choices=['hills','stats','features'],default='features')
a=p.parse_args()
out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
manifest=dict(alphapept=importlib.metadata.version('alphapept'),numpy=np.__version__,arrays={},stages={})
assert manifest['alphapept']=='0.5.3'
def save(name,value,dtype=None):
    data=np.ascontiguousarray(value,dtype=dtype)
    path=out/(name+'.bin');data.tofile(path)
    manifest['arrays'][name]=dict(dtype=str(data.dtype),shape=list(data.shape),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
def mark(name):
    manifest['stages'][name]=time.time()
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print(name,flush=True)
with h5py.File(a.hdf) as h:
    g=h['Raw/MS1_scans']; rt=g['rt_list_ms1'][:];bounds=g['indices_ms1'][:]
    selected=np.flatnonzero((rt>=a.rt_start)&(rt<=a.rt_end));lo,hi=int(selected[0]),int(selected[-1])+1
    begin,end=int(bounds[lo]),int(bounds[hi])
    q=dict(indices_ms1=(bounds[lo:hi+1]-begin).astype(np.int64),rt_list_ms1=rt[lo:hi],
           mass_list_ms1=g['mass_list_ms1'][begin:end],int_list_ms1=g['int_list_ms1'][begin:end])
    save('scan_numbers',g['scan_list_ms1'][lo:hi],'<i4')
for k,v in q.items():save(k,v,'<i8' if k=='indices_ms1' else '<f8')
ptr,data,counts,median,std=ff.extract_hills(q,2,8)
for k,v in [('first_ptr',ptr),('first_data',data),('first_counts',counts)]:save(k,v,'<i4')
manifest['first_score_median']=float(median);manifest['first_score_std']=float(std)
tol=float(median+3*std);manifest['adaptive_tolerance']=tol
ptr,data,counts,_,_=ff.extract_hills(q,2,tol)
for k,v in [('second_ptr',ptr),('second_data',data),('second_counts',counts)]:save(k,v,'<i4')
ptr,data=ff.remove_duplicate_hills(ptr,data,counts)
save('dedup_ptr',ptr,'<i4');save('dedup_data',data,'<i4')
ptr=ff.split_hills(ptr,data,q['int_list_ms1'],1.3,1)
save('split_ptr',ptr,'<i4')
data,ptr=ff.filter_hills(data,ptr,q['int_list_ms1'],40,1)
save('filtered_ptr',ptr,'<i4');save('filtered_data',data,'<i4');mark('hills')
if a.through=='hills':raise SystemExit
stats,sortindex,upper,scanidx,data,ptr=ff.get_hill_data(q,ptr,data,300,150)
for k,v in [('stats',stats),('sortindex',sortindex),('upper',upper),('scanidx',scanidx)]:save(k,v,'<f8' if k=='stats' else '<i8')
mark('stats')
if a.through=='stats':raise SystemExit
pre=ff.get_pre_isotope_patterns(stats,upper,sortindex,ptr,data,q['int_list_ms1'],scanidx,ff.maximum_offset)
save('pre_ptr',np.r_[0,np.cumsum([len(x) for x in pre])],'<i8');save('pre_data',np.concatenate(pre),'<i8');mark('pre_isotopes')
from alphapept.constants import averagine_aa,isotopes
patterns,isoidx,charges=ff.get_isotope_patterns(pre,ptr,data,q['int_list_ms1'],scanidx,stats,sortindex,averagine_aa,isotopes)
for k,v in [('patterns',patterns),('isoidx',isoidx),('charges',charges)]:save(k,v,'<i8')
report,lookup=ff.feature_finder_report(q,patterns,charges,isoidx,stats,sortindex,ptr,data)
report.to_csv(out/'feature_table.tsv',sep='\t',index=False);save('lookup',lookup,'<i4');mark('features')
