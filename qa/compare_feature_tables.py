"""QA only: compare feature IDs and every shared AlphaPept numerical column."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

p=argparse.ArgumentParser()
p.add_argument('reference')
p.add_argument('candidate')
p.add_argument('output')
a=p.parse_args()
ref=pd.read_csv(a.reference,sep='\t')
candidate=pd.read_csv(a.candidate,sep='\t').rename(columns={'intensity_apex':'ms1_int_sum_apex','area_min':'ms1_int_sum_area'})
if 'feature_id' not in ref:ref.insert(0,'feature_id',np.arange(len(ref)))
ref=ref.set_index('feature_id');candidate=candidate.set_index('feature_id')
result=dict(reference_features=len(ref),candidate_features=len(candidate),same_ids=ref.index.equals(candidate.index),columns={})
common=ref.index.intersection(candidate.index)
for name in ref.columns.intersection(candidate.columns):
    x=ref.loc[common,name].to_numpy();y=candidate.loc[common,name].to_numpy()
    result['columns'][name]=dict(max_abs_error=float(np.max(np.abs(x-y))),
        exact_fraction=float(np.mean(x==y)),agreement_rtol_1e_12=float(np.mean(np.isclose(x,y,rtol=1e-12,atol=1e-8))))
result['pass']=bool(result['same_ids'] and all(x['agreement_rtol_1e_12']==1 for x in result['columns'].values()))
Path(a.output).write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2),flush=True)
raise SystemExit(0 if result['pass'] else 1)
