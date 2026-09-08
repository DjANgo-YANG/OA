from pathlib import Path
import hashlib,json,urllib.request,concurrent.futures,sys
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
base=Path(__file__).resolve().parent.parent
out=base/'results'
files=list(base.glob('GSE*.txt'))+[base/'GSE117999_normalized_matrix.xlsx']
def hashes():return {f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
before=hashes()
def fetch(acc):
 p=out/(acc+'_source.soft.txt')
 if p.exists():return acc,p.read_text(encoding='utf-8')
 url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+acc+'&targ=self&form=text&view=full'
 t=urllib.request.urlopen(url,timeout=45).read().decode('utf-8')
 assert '= '+acc in t,'Invalid GEO response'
 p.write_text(t,encoding='utf-8');return acc,t
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 sources=dict(pool.map(fetch,['GSE114007','GSM3130531','GSM3130549','GSE286154','GSM8718799','GSM8718807','GSE117999','GSM3316962']))
def read(name):
 p=base/name
 with p.open(encoding='utf-8-sig') as f:h=f.readline().strip().split('\t')
 assert len(h)==len(set(h)),'Duplicate headers: stop'
 d=pd.read_csv(p,sep='\t',keep_default_na=False)
 assert d.iloc[:,0].is_unique and d.iloc[:,0].ne('').all(),'Ambiguous IDs: stop'
 return d
summary={}
def qc(x,key,lines):
 a=x.to_numpy(dtype=float);assert np.isfinite(a).all() and x.columns.is_unique
 q=np.quantile(a,[0,.01,.25,.5,.75,.99,1])
 stat=dict(expression_shape=list(a.shape),missing_nonfinite=0,all_integer=bool((a==np.floor(a)).all()),zero_fraction=float((a==0).mean()),negative_fraction=float((a<0).mean()),quantiles_0_1_25_50_75_99_100=q.tolist(),constant_rows=int((np.ptp(a,axis=1)==0).sum()),all_zero_rows=int((a==0).all(axis=1).sum()))
 pd.DataFrame(dict(sample_id=x.columns,sum=a.sum(axis=0),minimum=a.min(axis=0),q25=np.quantile(a,.25,axis=0),median=np.median(a,axis=0),q75=np.quantile(a,.75,axis=0),maximum=a.max(axis=0),zero_fraction=(a==0).mean(axis=0))).to_csv(out/(key+'_sample_QC.csv'),index=False)
 x.corr(method='spearman').to_csv(out/(key+'_sample_spearman.csv'))
 v=np.log2(a+1) if key=='GSE286154' else a
 fig,ax=plt.subplots(figsize=(max(10,x.shape[1]*.34),5));ax.boxplot(v,labels=x.columns,showfliers=False);ax.tick_params(axis='x',rotation=90);ax.set_ylabel('log2(raw counts+1), display only' if key=='GSE286154' else 'Supplied expression scale');ax.set_title(key);fig.tight_layout();fig.savefig(out/(key+'_QC_boxplot.png'),dpi=150);plt.close(fig)
 z=v.T-v.T.mean(axis=0);u,s,_=np.linalg.svd(z,full_matrices=False);pc=u[:,:2]*s[:2];var=s*s/(s*s).sum();stat['PCA_variance_first2']=var[:2].tolist()
 pd.DataFrame(dict(sample_id=x.columns,PC1=pc[:,0],PC2=pc[:,1])).to_csv(out/(key+'_QC_PCA.csv'),index=False)
 fig,ax=plt.subplots(figsize=(8,6));ax.scatter(pc[:,0],pc[:,1])
 for i,c in enumerate(x.columns):ax.annotate(c,pc[i],fontsize=6)
 ax.set_xlabel('PC1 %.1f%%'%(100*var[0]));ax.set_ylabel('PC2 %.1f%%'%(100*var[1]));ax.set_title(key+' descriptive PCA');fig.tight_layout();fig.savefig(out/(key+'_QC_PCA.png'),dpi=150);plt.close(fig)
 lines+=['QC: no sample removal. PCA uses all features, gene centering without variance scaling; raw-count PCA uses log2(count+1) for display only.','Numeric QC:',json.dumps(stat,indent=2)]
 (out/(key+'_QC_report.txt')).write_text('\n'.join(lines)+'\n',encoding='utf-8');summary[key]=stat;print(key,stat,flush=True)
n=read('GSE114007_normal_normalized.counts.txt');o=read('GSE114007_OA_normalized.counts.txt')
assert n.columns[0]==o.columns[0]=='symbol' and set(n.symbol)==set(o.symbol),'Different gene sets: stop'
ns=[c for c in n if c not in ['symbol','Average Normal','Max']];os=[c for c in o if c not in ['symbol','Average OA','Max']]
assert len(ns)==18 and len(os)==20 and not set(ns)&set(os)
x=n.set_index('symbol')[ns].join(o.set_index('symbol')[os],validate='one_to_one');assert np.isfinite(x.to_numpy()).all()
x.to_csv(out/'GSE114007_merged_expression.tsv',sep='\t',index_label='gene_id')
pd.DataFrame(dict(sample_id=ns+os,condition=['Normal']*18+['OA']*20,source_file=['GSE114007_normal_normalized.counts.txt']*18+['GSE114007_OA_normalized.counts.txt']*20,platform=['GPL11154' if '_Cart_' in c else 'GPL18573' for c in ns+os])).to_csv(out/'GSE114007_metadata.csv',index=False)
qc(x,'GSE114007',['COMPLETE: merge and basic QC.','Input dimensions: Normal 23710 x 21, OA 23710 x 23 including ID and summary columns.','Output: 23710 genes x 38 samples; 18 Normal and 20 OA. Gene ID: symbols, no blanks or duplicates.','Gene sets identical, original order differs. Joined by exact symbol, preserving Normal gene order. Average Normal/Average OA/Max are excluded. No values changed.','Normalized: YES. GEO describes TMM-normalized CPM and limma-voom. Negative fractional values are consistent with log2-CPM; exact deposited log formula/prior count not stated. No additional normalization or log transformation.','Future DE: preferably raw counts with TMM + limma-voom and verified platform/covariates. Existing log matrix can support limma, but original voom precision weights are unavailable. Do not input supplied log values into DESeq2 or rerun voom on them.','Two platforms: GPL11154 for Cart samples, GPL18573 for other samples; verify platform effects before DE. No DE performed.','Sources: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3130531 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3130549'])
d=read('GSE286154_gene_count.txt');x=d.set_index('gene');a=x.to_numpy();assert np.isfinite(a).all() and (a>=0).all() and (a==np.floor(a)).all() and a.max()<np.iinfo(np.int32).max
assert 'Gene expression counts' in sources['GSM8718799']
pd.DataFrame(dict(sample_id=x.columns,label_from_column=[c.rstrip('0123456789') for c in x.columns],condition='unassigned',patient='unassigned')).to_csv(out/'GSE286154_metadata.csv',index=False)
qc(x,'GSE286154',['INPUT QC COMPLETE; normalization performed by normalize_GSE286154.R.','14927 gene symbols x 16 samples: MC1-8 and LC1-8. Duplicate gene IDs: 0; blank IDs: 0.','Raw nonnegative integer counts supported by complete numeric checks and GEO supplementary gene-count description. No normalization reported. GEO quantifier is imprecisely described as featureCounts or StringTie; exact tool cannot be determined.','Apply DESeq2 default ratio size-factor normalization only, design ~1; retain every gene and raw counts. No DE fitting or testing.','Study: 8 OA patients with paired smooth/damaged cartilage, NOT healthy Normal vs OA. MC/LC direction and individual patient mapping require verification before DE; no biological groups or pair IDs guessed.','Future DE: DESeq2 on original counts, patient + condition design after mapping is verified. Do not feed normalized matrix back as raw counts.','Sources: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286154 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8718799 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8718807'])
f=base/'GSE117999_normalized_matrix.xlsx';assert pd.ExcelFile(f).sheet_names==['Matrix Template']
d=pd.read_excel(f,keep_default_na=False);assert d.columns[:5].tolist()==['FeatureNum','ProbeUID','ControlType','ProbeName','GeneName'] and d.FeatureNum.is_unique
x=d.iloc[:,5:];d.iloc[:,:5].to_csv(out/'GSE117999_feature_annotation.csv',index=False)
qc(x,'GSE117999',['BASIC QC COMPLETE; additional transformations STOPPED pending confirmation of exact normalization/log2 provenance.','Agilent GPL20844 microarray, NOT RNA-seq. One worksheet: Matrix Template. 62976 features x 25 columns = 5 annotations + 20 samples.','First column FeatureNum is chip feature index, NOT gene ID. ProbeName: Agilent probes including controls. GeneName: mixed symbols, transcript accessions, lncRNA labels and unmapped probes. No automatic gene mapping or aggregation.','ControlType counts: '+str(d.ControlType.value_counts().to_dict()),'Duplicate ProbeName rows: '+str(d.ProbeName.duplicated().sum())+'; duplicate GeneName rows: '+str(d.GeneName.duplicated().sum())+'; blank GeneName rows: '+str(d.GeneName.eq('').sum()),'Normalized: supplied by GEO as normalized_matrix. Exact algorithm unconfirmed. Range 3.6292 to 19.0466 is consistent with log2 intensities but does not prove log base. Log2: likely, unconfirmed. No renormalization or log transformation.','QC retains controls; no removal of probes, duplicates or samples. 60901 ControlType=0 biological feature rows, 2075 controls.','GEO describes 24 initial samples, excluding P4-010/P4-012/P4-104/P4-108, leaving 20 supplied columns.','Future DE: microarray limma after confirming scale, group metadata and probe/control handling. Do not use DESeq2, count-based edgeR or voom on these intensities.','Sources: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117999 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3316962'])
after=hashes();assert before==after
(out/'input_integrity_SHA256.json').write_text(json.dumps(dict(before=before,after=after,unchanged=True),indent=2),encoding='utf-8')
(out/'QC_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
(out/'python_versions.txt').write_text('Python '+sys.version+'\npandas '+pd.__version__+'\nnumpy '+np.__version__+'\nmatplotlib '+matplotlib.__version__,encoding='utf-8')
