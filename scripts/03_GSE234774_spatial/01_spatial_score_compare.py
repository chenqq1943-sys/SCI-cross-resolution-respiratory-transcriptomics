# =====================================================================
# Phase 2D spatial pilot - GSE234774 2D Visium
# Predefined programs: ComplexI_assembly, Aerobic_respiration
# Units: spot (within-section), section (spatial subsample), biological replicate (slide)
# No new pathways, no SIRT6, no deconvolution, no clustering, no region definition.
# =====================================================================
import gzip, csv, io, os, time, math
import numpy as np
from scipy import io as sio, sparse
from scipy import stats

base=os.environ.get("SCI_RAW_DIR", os.getcwd()); repo_root=os.environ.get("SCI_REPOSITORY_ROOT", os.getcwd()); recon=os.path.join(base,"recon"); out=os.path.join(base,"Phase2D_spatial_pilot")
os.makedirs(out, exist_ok=True)
LOG=[]
def note(s): print(s, flush=True); LOG.append(s)
t0=time.time()
note("Phase2D spatial pilot started: "+time.strftime("%Y-%m-%d %H:%M:%S"))

# ---------------- 1. load data ----------------
mtx = sio.mmread(recon+"/GSE234774_spatial_2d_filtered_spatial_2d.mtx.gz").tocsr()
M,N = mtx.shape
note(f"mtx loaded: {M} x {N}, nnz={mtx.nnz}")
with gzip.open(recon+"/GSE234774_spatial_2d_features.txt.gz","rt") as fh:
    feat=[l.rstrip("\n") for l in fh]
with gzip.open(recon+"/GSE234774_spatial_2d_barcodes.txt.gz","rt") as fh:
    bc=[l.rstrip("\n") for l in fh]
with gzip.open(recon+"/GSE234774_spatial_2d_meta.txt.gz","rb") as fh:
    meta=list(csv.DictReader(io.TextIOWrapper(fh,encoding="utf-8"), delimiter="\t"))
Nfeat=len(feat); Nbc=len(bc); Nmeta=len(meta)
note(f"features={Nfeat} barcodes={Nbc} meta={Nmeta}")

# meta order vs barcodes order
meta_bc=[x["barcode"] for x in meta]
order_ok = (meta_bc==bc)
note(f"meta barcode order == barcodes file order: {order_ok}")
if not order_ok:
    # reorder meta to match bc
    mpos={x["barcode"]:i for i,x in enumerate(meta)}
    meta=[meta[mpos[b]] for b in bc]
    note("meta reordered to matrix column order")

# ---------------- 2. technical QC ----------------
colsum=np.asarray(mtx.sum(axis=0)).ravel()
meta_ncount=np.array([float(x["nCount_Spatial"]) for x in meta])
qc_rows=[]
def qcrow(check,status,detail): qc_rows.append([check,status,detail])
qcrow("matrix dims", "PASS", f"{M}x{N} nnz={mtx.nnz}")
qcrow("rows==features","PASS" if M==Nfeat else "FAIL", str(M==Nfeat))
qcrow("cols==barcodes","PASS" if N==Nbc else "FAIL", str(N==Nbc))
qcrow("cols==meta","PASS" if N==Nmeta else "FAIL", str(N==Nmeta))
qcrow("meta order==barcodes","PASS" if order_ok else "REORDERED", str(order_ok))
mismatch=np.sum(np.abs(colsum-meta_ncount)>1)
qcrow("colsum vs meta nCount mismatch", "PASS" if mismatch==0 else "WARN", f"{mismatch}/{N}")
# low-UMI threshold (pre-specified from tech QC distribution): 10th percentile of total UMI
umi_low=np.percentile(colsum,10)
qcrow("pre-specified low-UMI threshold","INFO",f"10th pct total UMI = {umi_low:.1f}")
note(f"total UMI range {colsum.min():.0f}-{colsum.max():.0f}; p10={umi_low:.0f}")
with open(out+"/spatial_technical_QC.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["check","status","detail"]); w.writerows(qc_rows)
note("spatial_technical_QC written")

# ---------------- 3. normalization: log1p(CP10K), keep sparse ----------------
scale=10000.0/np.maximum(colsum,1.0)
dat=mtx.data.astype(np.float64)*scale[mtx.indices]
norm=sparse.csr_matrix((dat, mtx.indices, mtx.indptr), shape=mtx.shape)
norm.data=np.log1p(norm.data)
note("normalized log1p(CP10K) built (sparse); nnz="+str(norm.nnz))

# ---------------- 4. target gene mapping ----------------
with open(os.path.join(repo_root,"gene_programs","locked_respiratory_gene_programs.csv")) as fh:
    tl=list(csv.DictReader(fh))
feat_set={g:i for i,g in enumerate(feat)}
complexi_syms=[r["GeneSymbol"] for r in tl if r["ComplexI_assembly"]=="TRUE"]
aerobic_syms=[r["GeneSymbol"] for r in tl if r["Aerobic_respiration"]=="TRUE"]
def map_genes(syms, tag):
    idx=[]; matched=[]; unmatched=[]
    for g in syms:
        if g in feat_set: idx.append(feat_set[g]); matched.append(g)
        else: unmatched.append(g)
    note(f"{tag}: original={len(syms)} matched={len(matched)} unmatched={len(unmatched)} {unmatched}")
    return np.array(idx,dtype=np.int64), matched, unmatched
ci_idx,ci_m,ci_u=map_genes(complexi_syms,"ComplexI_assembly")
ae_idx,ae_m,ae_u=map_genes(aerobic_syms,"Aerobic_respiration")
with open(out+"/spatial_target_gene_mapping.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["program","gene","in_visium_features"])
    for g in complexi_syms: w.writerow(["ComplexI_assembly",g, str(g in feat_set)])
    for g in aerobic_syms: w.writerow(["Aerobic_respiration",g, str(g in feat_set)])
note("spatial_target_gene_mapping written")

progs={"ComplexI_assembly":ci_idx,"Aerobic_respiration":ae_idx}
prog_syms={"ComplexI_assembly":ci_m,"Aerobic_respiration":ae_m}

# ---------------- 5. program scoring ----------------
def score_methodA(prog_idx):
    nprog=len(prog_idx); nz=norm.data; ridx=norm.indices
    gsum=np.bincount(ridx,weights=nz,minlength=Nfeat)
    gsumsq=np.bincount(ridx,weights=nz*nz,minlength=Nfeat)
    mean_g=gsum/N; var_g=gsumsq/N-mean_g**2; var_g[var_g<0]=0; sd_g=np.sqrt(var_g)
    px=norm[prog_idx].toarray()  # (nprog,N) small
    gm=mean_g[prog_idx]; gs=sd_g[prog_idx]
    gs_ok=gs>1e-9
    z=np.zeros_like(px)
    if gs_ok.any():
        z[gs_ok]=(px[gs_ok]-gm[gs_ok][:,None])/gs[gs_ok][:,None]
    score=z.mean(axis=0)
    return score, nprog

def score_methodB(prog_idx):
    nprog=len(prog_idx); nz=norm.data; ridx=norm.indices
    nz=norm.data.astype(np.float64); csc=norm.tocsc()
    d=csc.data; p=csc.indptr
    px=norm[prog_idx].toarray()  # (nprog,N)
    nfeat=N
    score=np.zeros(N,dtype=np.float64)
    for j in range(N):
        sv=np.sort(d[p[j]:p[j+1]])
        nz_j=len(sv); n_zero=nfeat-nz_j
        v=px[:,j]
        fr=np.empty(nprog)
        for g in range(nprog):
            val=v[g]
            if val<=0:
                fr[g]=(0.0+0.5*n_zero)/nfeat
            else:
                lo=np.searchsorted(sv,val,'left'); hi=np.searchsorted(sv,val,'right')
                n_less=n_zero+lo; n_eq=hi-lo
                fr[g]=(n_less+0.5*n_eq)/nfeat
        score[j]=fr.mean()
    return score, nprog

scores={}  # prog -> {'A':array,'B':array}
for pn,pix in progs.items():
    A,npa=score_methodA(pix); B,npb=score_methodB(pix)
    scores[pn]={'A':A,'B':B}
    note(f"{pn}: methodA(A)={A.min():.3f}/{A.max():.3f} methodB(B)={B.min():.4f}/{B.max():.4f}")
note("scoring done in %.0fs"%(time.time()-t0))

# ---------------- save per-spot scores ----------------
gmap={"uninjured":"uninjured","7days":"7d","2months":"2m"}
with gzip.open(out+"/spatial_program_scores.csv.gz","wt",newline="") as fh:
    w=csv.writer(fh)
    w.writerow(["barcode","section_id","slide","spot","condition","time_label","replicate","sections","x","y","program","method","score"])
    for j in range(N):
        mj=meta[j]
        for pn in progs:
            for met in ("A","B"):
                w.writerow([mj["barcode"], mj["id"], mj["slide"], mj["spot"], mj["group"],
                            gmap.get(mj["group"],mj["group"]), mj["replicate"], mj["sections"],
                            mj["x"], mj["y"], pn, met, "%.6f"%scores[pn][met][j]])
note("spatial_program_scores written")

# meta arrays
grp=np.array([x["group"] for x in meta]); rep=np.array([x["replicate"] for x in meta])
sec_id=np.array([x["id"] for x in meta]); slide=np.array([x["slide"] for x in meta])
xs=np.array([float(x["x"]) for x in meta]); ys=np.array([float(x["y"]) for x in meta])
ncount=np.array([float(x["nCount_Spatial"]) for x in meta])
sections=np.array([x["sections"] for x in meta])
condition=grp  # uninjured/7days/2months

# ---------------- 6. section-level summary ----------------
sec_ids=sorted(set(sec_id))
sec_rows=[]
for sid in sec_ids:
    sel=np.where(sec_id==sid)[0]
    cd=condition[sel][0]; rp=rep[sel][0]; sl=slide[sel][0]
    for pn in progs:
        for met in ("A","B"):
            s=scores[pn][met][sel]
            med=float(np.median(s)); iqr=float(np.percentile(s,75)-np.percentile(s,25))
            p10=float(np.percentile(s,10)); p90=float(np.percentile(s,90))
            cv=float(np.std(s)/ (np.mean(s) if np.mean(s)!=0 else np.nan))
            sec_rows.append([sid, gmap.get(cd,cd), rp, sl, pn, met, len(sel), round(med,5), round(iqr,5), round(p10,5), round(p90,5), round(cv,5)])
with open(out+"/section_level_summary.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["section_id","condition","replicate","slide","program","method","n_spots","median","IQR","p10","p90","spatial_CV"]); w.writerows(sec_rows)
note("section_level_summary written")

# ---------------- 7. thresholds + replicate-level summary ----------------
# thresholds pre-specified from uninjured spots (per program x method)
lowq,highq=0.10,0.90
uni_mask=condition=="uninjured"
thr={}
for pn in progs:
    for met in ("A","B"):
        su=scores[pn][met][uni_mask]
        thr[(pn,met)]=(float(np.percentile(su,lowq*100)), float(np.percentile(su,highq*100)))
note("thresholds (low,high) from uninjured spots: "+str({f"{k[0]}/{k[1]}":(round(v[0],4),round(v[1],4)) for k,v in thr.items()}))

bio_ids=sorted(set(zip(condition,rep)))
rep_rows=[]
for (cd,rp) in bio_ids:
    sel=np.where((condition==cd)&(rep==rp))[0]
    sl=slide[sel][0]
    for pn in progs:
        for met in ("A","B"):
            s=scores[pn][met][sel]
            lo,hi=thr[(pn,met)]
            low_frac=float(np.mean(s<lo)); high_frac=float(np.mean(s>hi))
            rep_rows.append([gmap.get(cd,cd), rp, sl, pn, met, len(sel),
                round(float(np.median(s)),5), round(float(np.percentile(s,75)-np.percentile(s,25)),5),
                round(low_frac,5), round(high_frac,5)])
with open(out+"/replicate_level_summary.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["condition","replicate","slide","program","method","n_spots","median","IQR","low_fraction","high_fraction"]); w.writerows(rep_rows)
note("replicate_level_summary written")

# ---------------- 8. primary comparisons ----------------
def cohend(x,y):
    nx,ny=len(x),len(y)
    s2=((nx-1)*np.var(x,ddof=1)+(ny-1)*np.var(y,ddof=1))/(nx+ny-2)
    if s2<=0 or math.isnan(s2): return None
    return (np.mean(y)-np.mean(x))/math.sqrt(s2)

def rep_median(cd,pn,met):
    vals=[]
    for (c,rp) in bio_ids:
        if c==cd:
            sel=np.where((condition==c)&(rep==rp))[0]
            vals.append(np.median(scores[pn][met][sel]))
    return np.array(vals)

prim_rows=[]
for pn in progs:
    for met in ("A","B"):
        u=rep_median("uninjured",pn,met)
        for cd_t in ("7days","2months"):
            x=rep_median(cd_t,pn,met)
            d=float(np.mean(x)-np.mean(u))
            dirn="down" if d<0 else "up"
            es=cohend(u,x)
            t,p=stats.ttest_ind(x,u,equal_var=False)
            md=float(np.median(u)); below=int(np.sum(x<md))
            prim_rows.append([pn,met,gmap[cd_t],3,3, round(float(np.mean(u)),5),round(float(np.mean(x)),5),
                round(d,5),dirn, ("" if es is None else round(es,4)), round(float(p),5), f"{below}/3"])
with open(out+"/primary_spatial_comparisons.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["program","method","comparison","n_control","n_injured","control_mean","injured_mean","delta","direction","effect_size","nominal_p","replicate_consistency"]); w.writerows(prim_rows)
note("primary_spatial_comparisons written")

# ---------------- 9. sensitivity ----------------
# A. low-UMI exclusion (pre-specified threshold p10 of total UMI), recompute scores on kept spots
keep_umi=colsum>=umi_low
raw_sub=mtx[:,keep_umi].tocsr()
cs_sub=np.asarray(raw_sub.sum(axis=0)).ravel()
scale_sub=10000.0/np.maximum(cs_sub,1.0)
dd=raw_sub.data.astype(np.float64)*scale_sub[raw_sub.indices]
nsub=sparse.csr_matrix((dd,raw_sub.indices,raw_sub.indptr),shape=raw_sub.shape)
nsub.data=np.log1p(nsub.data)
def scoreA_sub(nsub,pix):
    nprog=len(pix); nz=nsub.data; ridx=nsub.indices
    gsum=np.bincount(ridx,weights=nz,minlength=Nfeat)
    gss=np.bincount(ridx,weights=nz*nz,minlength=Nfeat)
    mn=gsum/raw_sub.shape[1]; vr=gss/raw_sub.shape[1]-mn**2; vr[vr<0]=0; sd=np.sqrt(vr)
    px=nsub[pix].toarray(); gm=mn[pix]; gs=sd[pix]; ok=gs>1e-9; z=np.zeros_like(px)
    if ok.any(): z[ok]=(px[ok]-gm[ok][:,None])/gs[ok][:,None]
    return z.mean(axis=0)
def scoreB_sub(nsub,pix):
    nprog=len(pix); nfeat=Nfeat; csc=nsub.tocsc(); d=csc.data; p=csc.indptr; px=nsub[pix].toarray()
    ncol=raw_sub.shape[1]; sc=np.zeros(ncol)
    for j in range(ncol):
        sv=np.sort(d[p[j]:p[j+1]]); nz_j=len(sv); n_zero=nfeat-nz_j
        v=px[:,j]; fr=np.empty(nprog)
        for g in range(nprog):
            val=v[g]
            if val<=0: fr[g]=(0.5*n_zero)/nfeat
            else:
                lo=np.searchsorted(sv,val,'left'); hi=np.searchsorted(sv,val,'right')
                fr[g]=(n_zero+lo+0.5*(hi-lo))/nfeat
        sc[j]=fr.mean()
    return sc

sensA=[]
for pn in progs:
    for met in ("A","B"):
        if met=="A": sc=scoreA_sub(nsub,progs[pn])
        else: sc=scoreB_sub(nsub,progs[pn])
        keep_idx=np.where(keep_umi)[0]
        def rmed(cd):
            o=[]
            for (c,rp) in bio_ids:
                if c==cd:
                    ss=[jj for jj,j in enumerate(keep_idx) if condition[j]==c and rep[j]==rp]
                    if ss: o.append(np.median(sc[ss]))
            return np.array(o)
        u=rmed("uninjured")
        for cd_t in ("7days","2months"):
            x=rmed(cd_t); d=np.mean(x)-np.mean(u)
            sensA.append([pn,met,"exclude_lowUMI",gmap[cd_t],round(float(d),5),"down" if d<0 else "up"])
note("sensitivity A (low-UMI) computed; n kept spots=%d"%(np.sum(keep_umi)))

# B. leave-one-section-out: replicate median from 3 sections
def rep_median_drop(cd,pn,met,drop_section=None):
    vals=[]
    for (c,rp) in bio_ids:
        if c==cd:
            sel=np.where((condition==c)&(rep==rp))[0]
            if drop_section is not None:
                sel=sel[sec_id[sel]!=drop_section]
            vals.append(np.median(scores[pn][met][sel]))
    return np.array(vals)
sensB=[]
for pn in progs:
    for met in ("A","B"):
        u0=rep_median("uninjured",pn,met)
        for cd_t in ("7days","2months"):
            x0=rep_median(cd_t,pn,met)
            d0=np.mean(x0)-np.mean(u0)
            flips=0; details=[]
            for sid in sec_ids:
                if condition[sec_id==sid][0]!=cd_t: continue
                # dropping this section from the injured replicate
                x1=rep_median_drop(cd_t,pn,met,sid)
                u1=rep_median_drop("uninjured",pn,met,sid)  # sid belongs to cd_t, so uninjured unaffected
                d1=np.mean(x1)-np.mean(u1)
                if (d0<0)!=(d1<0): flips+=1; details.append(sid)
            sensB.append([pn,met,gmap[cd_t],str(flips),";".join(details) if details else ""])
note("sensitivity B (leave-one-section-out) done")

# C. method agreement
sensC=[]
for cd_t in ("uninjured","7days","2months"):
    mask=condition==cd_t
    for pn in progs:
        a=scores[pn]["A"][mask]; b=scores[pn]["B"][mask]
        r=float(np.corrcoef(a,b)[0,1])
        # high/low region agreement: Jaccard of high regions
        lo,hi=thr[(pn,"A")]; lo2,hi2=thr[(pn,"B")]
        ha=a>hi; hb=b>hi2
        inter=np.sum(ha&hb); union=np.sum(ha|hb)
        sensC.append([pn,cd_t,round(r,4), inter, union, round(inter/union,4) if union else None])
note("sensitivity C (method agreement) done")

with open(out+"/spatial_sensitivity.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["sensitivity","program","method","comparison","detail","direction_delta","direction"])
    for r in sensA: w.writerow(["A_lowUMI",r[0],r[1],r[3],r[2],r[4],r[5]])
    for r in sensB: w.writerow(["B_leave_one_section",r[0],r[1],r[2],r[4],r[3],"n_flips"])
    for r in sensC: w.writerow(["C_method_agreement",r[0],r[1],r[2],f"jac_high={r[3]}/{r[4]}={r[5]}",r[2],"corr"])
note("spatial_sensitivity written")

# ---------------- 10. figures ----------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

cond2label={"uninjured":"Uninjured","7days":"7d","2months":"2m"}
sec_ids_sorted=sorted(sec_ids, key=lambda s:(slide[sec_id==s][0], condition[sec_id==s][0], sections[sec_id==s][0]))

def draw_maps(pn):
    fig,axes=plt.subplots(6,6,figsize=(24,22))
    vmin=np.percentile(scores[pn]["A"],1); vmax=np.percentile(scores[pn]["A"],99)
    for ax,sid in zip(axes.ravel(),sec_ids_sorted):
        sel=np.where(sec_id==sid)[0]
        sc=scores[pn]["A"][sel]
        im=ax.scatter(xs[sel],ys[sel],c=sc,s=3,cmap="viridis",vmin=vmin,vmax=vmax)
        ax.set_title(f"{sid}\n{cond2label[condition[sel][0]]} rep{rep[sel][0]}",fontsize=7)
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"GSE234774 2D Visium - {pn} (Method A score), 36 sections, x-y spatial map",fontsize=13)
    fig.tight_layout(); fig.savefig(f"{out}/{pn.replace('_assembly','').replace('_respiration','Respiration')}_36section_spatial_maps.pdf")
    plt.close(fig); note(f"figure: {pn} 36-section maps written")

draw_maps("ComplexI_assembly")
draw_maps("Aerobic_respiration")

# replicate-level summary figure: rows=program, cols=condition; per panel 3 replicates, 2 boxes (A,B)
fig,axes=plt.subplots(2,3,figsize=(15,9))
cond_order=["uninjured","7days","2months"]
for ri,pn in enumerate(progs):
    for ci,cd_t in enumerate(cond_order):
        ax=axes[ri,ci]
        reps=sorted([rp for (c,rp) in bio_ids if c==cd_t])
        xpos=[]; labs=[]
        for k,rp in enumerate(reps):
            sel=np.where((condition==cd_t)&(rep==rp))[0]
            for mi,met in enumerate(("A","B")):
                bp=ax.boxplot([scores[pn][met][sel]],positions=[k*3+mi],widths=0.8,patch_artist=True)
                bp["boxes"][0].set_facecolor("#e15759" if met=="A" else "#4e79a7")
            xpos.append(k*3+0.5); labs.append(f"R{rp}")
        ax.set_xticks(xpos); ax.set_xticklabels(labs,fontsize=8)
        ax.set_title(f"{cond2label[cd_t]} - {pn}",fontsize=9)
        ax.set_ylabel("score",fontsize=8)
        if ri==0 and ci==0:
            from matplotlib.patches import Patch
            ax.legend(handles=[Patch(color="#e15759",label="Method A"),Patch(color="#4e79a7",label="Method B")],fontsize=7,loc="best")
fig.suptitle("GSE234774 2D Visium - replicate-level program score distribution (spots) by condition",fontsize=12)
fig.tight_layout(); fig.savefig(out+"/replicate_level_summary.pdf"); plt.close(fig)
note("replicate_level_summary figure written")

# ---------------- 11. sessionInfo / log ----------------
import sys, subprocess
env=[]
env.append("python "+sys.version)
for m in ["numpy","scipy","matplotlib"]:
    try: env.append(f"{m} {__import__(m).__version__}")
    except Exception as e: env.append(f"{m} MISSING")
with open(out+"/sessionInfo.txt","w") as fh:
    fh.write("GSE234774 2D Visium Phase 2D spatial pilot - environment\n")
    fh.write("\n".join(env)+"\n")
    fh.write("methodA: mean standardized expression (per-gene z across all spots, averaged over program genes)\n")
    fh.write("methodB: rank-based fractional-rank gene-set score (per spot, rank relative to all genes, ties=average)\n")
    fh.write("normalization: log1p(counts / spot_total_UMI * 10000)\n")
    fh.write("thresholds: low=10th, high=90th percentile of uninjured spot scores (per program x method)\n")
with open(out+"/analysis_log.txt","w") as fh:
    fh.write("Phase 2D spatial pilot log\nStarted: "+time.strftime("%Y-%m-%d %H:%M:%S")+"\n")
    fh.write("\n".join(LOG)+"\nFinished: "+time.strftime("%Y-%m-%d %H:%M:%S"))
note("sessionInfo + analysis_log written")
note("DONE in %.0fs"%(time.time()-t0))
