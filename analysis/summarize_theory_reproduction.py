from pathlib import Path
import json
R=Path(__file__).resolve().parents[1]
old=json.loads((R/'analysis/theory/curvi_test.json').read_text());new=json.loads((R/'data/checks/curvi_test.json').read_text())
keys=['cfd_over_lub','theory_over_lub','theory_over_cfd','no_curv_over_cfd'];errors=[];metric={}
for g in old:
 metric[g]={}
 for sig in old[g]:
  metric[g][sig]={}
  for cut,vals in old[g][sig].items():
   if cut=='terms_6_12':continue
   for k in keys:errors.append(abs(new[g][sig][cut][k]-vals[k]))
   metric[g][sig][cut]={k:new[g][sig][cut][k] for k in ['raw_signed_metric_min','filtered_signed_metric_min','raw_metric_nonpositive_fraction']}
out={'comparisons':80,'maximum_absolute_ratio_difference':max(errors),'metric_checks':metric}
(R/'data/checks/theory_reproduction_summary.json').write_text(json.dumps(out,indent=2))
print('Maximum ratio difference:',max(errors))
for g in metric:
 for sig in ['sig2','sig3','sig4']:
  print(g,sig,'raw min',min(v['raw_signed_metric_min'] for v in metric[g][sig].values()),'raw fraction 6/12',metric[g][sig]['6_12']['raw_metric_nonpositive_fraction'])
assert max(errors)<1e-7
