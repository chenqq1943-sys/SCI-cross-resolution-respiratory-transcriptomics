# =====================================================================
# Rebuild FULL-transcriptome pseudobulk from GSE234774 snRNA raw counts
# Same cohort/celltype/nucleus rules as original Phase 2C slicing.
# Output: gene x (library x celltype) raw UMI counts for ALL genes.
# =====================================================================
import gzip, csv, io, time, collections, os
import numpy as np
base = os.environ.get("SCI_RAW_DIR", os.getcwd())
outdir = os.path.join(base, "Phase2C_corrected_fulltranscriptome_pseudobulk")
os.makedirs(outdir, exist_ok=True)
t0=time.time()
log=[]
def note(s): log.append(s); print(s, flush=True)

mtx_p=os.path.join(base,"GSE234774_rnaseq_filtered_scRNA.mtx.gz")
bc_p =os.path.join(base,"GSE234774_rnaseq_barcodes.txt.gz")
ft_p =os.path.join(base,"GSE234774_rnaseq_features.txt.gz")
meta_p=os.path.join(base,"GSE234774_rnaseq_meta.txt.gz")

# ---- load features, barcodes, meta ----
with gzip.open(ft_p,"rt") as f: features=[l.rstrip("\n") for l in f]
with gzip.open(bc_p,"rt") as f: barcodes=[l.rstrip("\n") for l in f]
with gzip.open(meta_p,"rt") as f:
    hdr=f.readline().rstrip("\n").split("\t"); ci={h:i for i,h in enumerate(hdr)}
    meta_rows=[l.rstrip("\n").split("\t") for l in f]
Nfeat=len(features); Nbc=len(barcodes); Nmeta=len(meta_rows)
note(f"features={Nfeat} barcodes={Nbc} meta={Nmeta}")

# ---- selection + celltype (exact original rules) ----
def assign(l1,l3):
    if l1=="Neurons": return "Neurons"
    if l3=="Astrocytes": return "Astrocytes"
    if l3 in ("Mature oligodendrocytes (MOL)","Mature oligodendrocytes (MOL), ischemic","Myelin forming oligodendrocytes (MFOL)"): return "Oligodendrocyte_combined"
    if l3=="Microglia": return "Microglia"
    return None
sel_mask=np.zeros(Nbc,dtype=bool)
sel_sample=np.full(Nbc,-1,dtype=np.int64)
sample_keys=[]
for i,r in enumerate(meta_rows):
    if r[ci["experiment"]] not in ("uninjured","timecourse"): continue
    ct=assign(r[ci["layer1"]],r[ci["layer3"]])
    if ct is None: continue
    key=(r[ci["library"]],ct)
    if key not in sample_keys: sample_keys.append(key)
    j=sample_keys.index(key)
    sel_mask[i]=True; sel_sample[i]=j
nSamp=len(sample_keys)
sel_idx=np.where(sel_mask)[0]
note(f"selected nuclei={sel_idx.size}  pseudobulk samples={nSamp}")

# nuclei counts per sample (from selection)
n_nuclei=np.bincount(sel_sample[sel_idx], minlength=nSamp).astype(np.int64)

# nCount_RNA sum per sample (independent cross-check)
ncount=np.zeros(nSamp,dtype=np.int64)
with gzip.open(meta_p,"rb") as fh:
    rdr=csv.DictReader(io.TextIOWrapper(fh,encoding="utf-8"),delimiter="\t")
    for r in rdr:
        if r["experiment"] not in ("uninjured","timecourse"): continue
        ct=assign(r["layer1"],r["layer3"])
        if ct is None: continue
        key=(r["library"],ct)
        ncount[sample_keys.index(key)]+=int(float(r["nCount_RNA"]))
note("meta nCount_RNA summed per sample")

# ---- stream mtx, accumulate ALL genes x sample ----
f=gzip.open(mtx_p,"rb")
# skip header: read until dims line
while True:
    line=f.readline(); p=line.split()
    if len(p)==3 and all(t.isdigit() for t in p):
        M,N,NNZ=map(int,p); break
note(f"mtx header: M={M} N={N} NNZ={NNZ} ; rows==features({M==Nfeat}) cols==barcodes({N==Nbc})")
acc=np.zeros(Nfeat*nSamp,dtype=np.int64)
buf=b""
while True:
    chunk=f.read(1<<23)
    if not chunk: break
    buf+=chunk
    lines=buf.split(b"\n"); buf=lines.pop()
    if not lines: continue
    block=b"\n".join(lines)
    arr=np.fromstring(block.decode("latin-1"), dtype=np.int64, sep=" ").reshape(-1,3)
    r=arr[:,0]-1; c=arr[:,1]-1; v=arr[:,2]
    keep=sel_mask[c]
    if keep.any():
        rr=r[keep]; cc=sel_sample[c[keep]]; vv=v[keep]
        np.add.at(acc, rr*nSamp+cc, vv)
if buf.strip():
    arr=np.fromstring(buf.decode("latin-1"),dtype=np.int64,sep=" ").reshape(-1,3)
    r=arr[:,0]-1;c=arr[:,1]-1;v=arr[:,2];keep=sel_mask[c]
    if keep.any():
        np.add.at(acc, r[keep]*nSamp+sel_sample[c[keep]], v[keep])
f.close()
note("mtx streamed; accumulation done in %.0fs"%(time.time()-t0))
acc=acc.reshape(Nfeat,nSamp)
full_umi=acc.sum(axis=0).astype(np.int64)
detected_genes=np.array((acc>0).sum(axis=0),dtype=np.int64)

# ---- cross-check: colSums vs meta nCount ----
diff=full_umi-ncount
maxabs=int(np.abs(diff).max())
rel=max(np.abs(diff)/(ncount+1e-9)).max()
note(f"colSums vs nCount_RNA: max|diff|={maxabs}  max rel diff={rel:.6f}")
mism=sum(np.abs(diff)>1)
note(f"samples with |diff|>1: {mism}/{nSamp}")
if mism>0:
    for i in range(nSamp):
        if np.abs(diff[i])>1:
            note("   mismatch: %s full=%d nCount=%d"%(sample_keys[i],full_umi[i],ncount[i]))

# ---- write outputs ----
def sname(a,b): return f"{a}__{b}"
cols=[sname(a,b) for a,b in sample_keys]
with gzip.open(os.path.join(outdir,"fulltranscriptome_pseudobulk_counts.csv.gz"),"wt",encoding="utf-8",newline="") as fh:
    w=csv.writer(fh)
    w.writerow(["GeneSymbol"]+cols)
    for g in range(Nfeat):
        w.writerow([features[g]]+[int(x) for x in acc[g]])
note("fulltranscriptome_pseudobulk_counts written")
with open(os.path.join(outdir,"fulltranscriptome_pseudobulk_coldata.csv"),"w",encoding="utf-8",newline="") as fh:
    w=csv.writer(fh)
    w.writerow(["pseudobulk_sample","library","time","celltype","n_nuclei","full_transcriptome_total_UMI","n_detected_genes","meta_nCount_RNA_sum"])
    times=[meta_rows[i][ci["label"]] for i in sel_idx]  # not per-sample; fix below
    # time per sample = label of first nucleus in that sample
    first_time={}
    for i in sel_idx:
        key=(meta_rows[i][ci["library"]], assign(meta_rows[i][ci["layer1"]],meta_rows[i][ci["layer3"]]))
        if key not in first_time: first_time[key]=meta_rows[i][ci["label"]]
    for oi,(a,b) in enumerate(sample_keys):
        w.writerow([cols[oi],a,first_time[(a,b)],b,int(n_nuclei[oi]),int(full_umi[oi]),int(detected_genes[oi]),int(ncount[oi])])
note("fulltranscriptome_pseudobulk_coldata written")

# ---- why 77 not 80 ----
libs=sorted(set(x[0] for x in sample_keys))
present=collections.defaultdict(set)
for (a,b) in sample_keys: present[a].add(b)
allct=sorted(set(x[1] for x in sample_keys))
missing=[]
for lb in libs:
    for ct in allct:
        if ct not in present[lb]: missing.append((lb,ct))
note(f"libraries={len(libs)} celltypes={len(allct)} expected={len(libs)*len(allct)} actual={nSamp} missing={len(missing)}")
for m in missing: note("   missing combo: %s / %s"%(m[0],m[1]))

log.append("finished %.0fs"%(time.time()-t0))
open(os.path.join(outdir,"_rebuild_log.txt"),"w",encoding="utf-8").write("\n".join(log))
print("DONE", flush=True)
