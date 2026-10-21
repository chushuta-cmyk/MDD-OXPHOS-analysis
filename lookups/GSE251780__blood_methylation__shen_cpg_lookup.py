"""Correct Table 2 extraction; re-extract all candidates, never overwrite sources."""
import csv, json, hashlib, gzip, sys
from pathlib import Path
# --- Path configuration: EDIT these placeholder roots for your machine ----------
DATA_ROOT    = "<DATA_ROOT>"
WORK_ROOT    = "<WORK_ROOT>"
ARCHIVE_ROOT = "<ARCHIVE_ROOT>"
OUTPUT_ROOT  = "./output"
OUT=Path(OUTPUT_ROOT)/'lookups'
OUT.mkdir(parents=True,exist_ok=True)
ROOT=Path(WORK_ROOT)/'analysis_framework/GSE251780'
SOURCE='https://www.nature.com/articles/s44220-025-00486-4/tables/2'
# Complete Table 2, checked against publisher. Source coefficients are standardized,
# local logFC is M-value difference: direction only; magnitudes must not be pooled.
raw='''cg19432791 LZTFL1 0.025 0.004 1.70e-9
cg11101030 NA -0.035 0.006 4.99e-9
cg10145533 NA 0.019 0.003 5.34e-9
cg12187985 SS18 0.021 0.004 5.43e-9
cg09726355 GPATCH8 0.019 0.003 8.53e-9
cg03671627 PIKFYVE 0.026 0.005 9.85e-9
cg17494199 NA -0.011 0.002 1.22e-8
cg17002328 CCDC88C 0.016 0.003 2.32e-8
cg12805451 CCDC14 0.024 0.004 3.11e-8
cg20793665 NA 0.019 0.003 3.79e-8
cg01150646 MYO1C 0.013 0.002 3.80e-8
cg02948555 ZNF106 0.036 0.007 4.72e-8
cg10604091 DLEU2 0.023 0.004 5.27e-8
cg04173586 DOT1L 0.025 0.005 6.30e-8
cg07063853 API5 0.020 0.004 6.42e-8'''
source=[]
for l in raw.splitlines():
    c,g,b,se,p=l.split()
    source.append(dict(CpG=c,Gene=g,Shen_standardized_beta=float(b),Shen_SE_rounded=float(se),Shen_P=float(p),Source=SOURCE))
assert len(source)==len({r['CpG'] for r in source})==15
def write(name,rows):
    with (OUT/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader();w.writerows(rows)
write('Shen2025_Table2_complete15.tsv',source)
targets={r['CpG'] for r in source}
results={}; hashes={}
for label,fn in [('sex_adjusted','GSE251780_MDD_DMP_all.tsv'),('cell_adjusted','GSE251780_MDD_DMP_cell_adjusted_all.tsv')]:
    p=ROOT/'results'/fn
    hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    with p.open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==723975
    sel=[r for r in rows if r['Probe'] in targets]
    assert len(sel)==len({r['Probe'] for r in sel})
    results[label]={r['Probe']:r for r in sel}
    del rows
out=[]
for s in source:
    for label in results:
        r=results[label].get(s['CpG']); q=dict(s,Model=label,Covered=r is not None)
        keys=['Local_M_difference','Local_moderated_SE','Local_moderated_t','Local_P','Local_genomewide_BH','Candidate_Bonferroni15','Sensitivity_Bonferroni30','Direction_match']
        q.update({k:None for k in keys})
        if r:
            b=float(r['logFC']);t=float(r['t']);p=float(r['P.Value'])
            q.update(Local_M_difference=b,Local_moderated_SE=abs(b/t) if t else None,Local_moderated_t=t,Local_P=p,Local_genomewide_BH=float(r['adj.P.Val']),Candidate_Bonferroni15=min(1,15*p),Sensitivity_Bonferroni30=min(1,30*p),Direction_match=('same' if b*s['Shen_standardized_beta']>0 else 'opposite' if b*s['Shen_standardized_beta']<0 else 'zero'))
        out.append(q)
write('Shen2025_crosscheck_corrected.tsv',out)
summary={}
for label in results:
    rr=[r for r in out if r['Model']==label and r['Covered']]
    summary[label]=dict(targets=15,covered=len(rr),direction_evaluable=len(rr),same_direction=sum(r['Direction_match']=='same' for r in rr),nominal=sum(r['Local_P']<.05 for r in rr),bonferroni15=sum(r['Candidate_Bonferroni15']<.05 for r in rr),bonferroni30=sum(r['Sensitivity_Bonferroni30']<.05 for r in rr),missing=[r['CpG'] for r in out if r['Model']==label and not r['Covered']])
(OUT/'step4_summary.json').write_text(json.dumps(summary,indent=2))
(OUT/'input_sha256.json').write_text(json.dumps(hashes,indent=2))
print(json.dumps(summary,indent=2))
for r in out:
    if r['CpG'] in ('cg04173586','cg20793665','cg01150646'):print(r)
# Extract only the 15 prespecified CpGs from submitter-processed matrix for
# target-level diagnostics; not a replacement for full-methylome eBayes.
matrix=ROOT/'data'/'GSE251780_Mokhtari_MatrixProcessed_AverageAndPvalue.txt.gz'
with gzip.open(matrix,'rt') as f, (OUT/'target15_processed_matrix.tsv').open('w') as g:
    g.write(next(f))
    seen=[]
    for line in f:
        if line.split('\t',1)[0].strip('"') in targets:
            g.write(line);seen.append(line.split('\t',1)[0])
print('Target probes present before local filtering:',len(seen))
