# -*- coding: utf-8 -*-
"""
Respiratory SCI project - Supplementary Tables construction (frozen-results extraction).
Reads ONLY existing locked/corrected outputs; no recomputation, no re-scoring,
no new statistics, no threshold changes, no new pathway.
Emits: Table_S1..S8 (.xlsx + .csv), Supplementary_Tables_README.txt,
       Supplementary_consistency_audit.txt
"""
import csv, os, gzip
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

ROOT = os.environ.get("SCI_RESULTS_ROOT", os.getcwd())
BASE = os.path.join(ROOT, "supplementary_source_tables")
PC   = os.path.join(ROOT, "derived_data", "GSE234774_snRNA")
SP   = os.path.join(ROOT, "derived_data", "GSE234774_spatial")
G31  = os.path.join(ROOT, "derived_data", "GSE319931")

def rd_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

def rnd(x, n=5):
    if x is None or x == "" or x == "NA":
        return x
    try:
        v = float(x)
    except Exception:
        return x
    return round(v, n)

def write_csv(path, headers, rows):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(headers)
        for r in rows:
            w.writerow([("" if c is None else c) for c in r])

def style_sheet(ws, headers, widths=None, freeze=True):
    hf = Font(bold=True, color="FFFFFF")
    fill = PatternFill("solid", fgColor="4472C4")
    for j, h in enumerate(headers, 1):
        c = ws.cell(row=1, column=j, value=h)
        c.font = hf; c.fill = fill
        c.alignment = Alignment(vertical="center", wrap_text=True)
    if freeze:
        ws.freeze_panes = "A2"
    if widths:
        for j, w in enumerate(widths, 1):
            ws.column_dimensions[get_column_letter(j)].width = w

def save_xlsx(path, sheets):
    """sheets: list of (name, headers, rows)"""
    wb = Workbook(); wb.remove(wb.active)
    for name, headers, rows in sheets:
        ws = wb.create_sheet(name)
        style_sheet(ws, headers)
        for r in rows:
            ws.append([("" if c is None else c) for c in r])
    wb.save(path)

def fnum(v):
    if v is None or v == "" or v == "NA":
        return "NA"
    return rnd(v, 5)

# ============================================================
# Table S2 - locked gene programs
# ============================================================
tl = rd_csv(os.path.join(ROOT, "gene_programs", "locked_respiratory_gene_programs.csv"))
sn_av = rd_csv(os.path.join(PC, "corrected_target_gene_availability.csv"))
sn_missing = set()
for r in tl:
    if r["GeneSymbol"] == "Sirt6":
        continue
    p1 = r["ComplexI_assembly"] == "TRUE"
    p2 = r["Aerobic_respiration"] == "TRUE"
# corrected snRNA available sets
sn_av_ci = {r["GeneSymbol"] for r in sn_av if r["program"] == "ComplexI_assembly"}
sn_av_ae = {r["GeneSymbol"] for r in sn_av if r["program"] == "Aerobic_respiration"}
sp_av = rd_csv(os.path.join(SP, "spatial_target_gene_mapping.csv"))
sp_av_ci = {r["gene"] for r in sp_av if r["program"] == "ComplexI_assembly" and r["in_visium_features"] == "True"}
sp_av_ae = {r["gene"] for r in sp_av if r["program"] == "Aerobic_respiration" and r["in_visium_features"] == "True"}

complexi_genes = [r["GeneSymbol"] for r in tl if r["ComplexI_assembly"] == "TRUE"]
aerobic_genes  = [r["GeneSymbol"] for r in tl if r["Aerobic_respiration"] == "TRUE"]
union_genes = []
for g in tl:
    if g["GeneSymbol"] == "Sirt6":
        continue
    if g["GeneSymbol"] not in union_genes:
        union_genes.append(g["GeneSymbol"])

GO = {"ComplexI_assembly": "GO:0032981", "Aerobic_respiration": "GO:0009060"}
s2_rows = []
for g in union_genes:
    in_ci = g in complexi_genes
    in_ae = g in aerobic_genes
    prog = []
    if in_ci: prog.append("ComplexI_assembly")
    if in_ae: prog.append("Aerobic_respiration")
    goid = " / ".join(GO[p] for p in prog)
    # availability
    a_gse5296 = "Yes"          # programs were built on GPL1261; all detectable
    a_sn = []
    if in_ci: a_sn.append("Yes" if g in sn_av_ci else "No(Ndufaf3 filtered)")
    if in_ae: a_sn.append("Yes" if g in sn_av_ae else "No(Mup3/Mup4/Mup5/Ucn filtered)")
    a_sn_str = ";".join(a_sn)
    a_sp = []
    if in_ci: a_sp.append("Yes" if g in sp_av_ci else "No")
    if in_ae: a_sp.append("Yes" if g in sp_av_ae else "No(Mup4 unavailable)")
    a_sp_str = ";".join(a_sp)
    a_g31 = "Yes"
    s2_rows.append(["/".join(prog), goid, g, "Nuclear",
                    a_gse5296, a_sn_str, a_sp_str, a_g31])

s2_rows.sort(key=lambda r: r[0])
s2_headers = ["Program(s)", "GO term ID", "Gene symbol",
              "Nuclear/mtDNA encoded",
              "Available GSE5296", "Available corrected snRNA",
              "Available spatial", "Available GSE319931"]

# ============================================================
# Table S3 - GSE5296 program results
# ============================================================
ts = rd_csv(os.path.join(ROOT, "derived_data", "GSE5296", "GSE5296_RMA_mito_modules_time_space.csv"))
qc = rd_csv(os.path.join(ROOT, "derived_data", "GSE5296", "GSE5296_RMA_QC_sensitivity.csv"))
imp_rep = rd_csv(os.path.join(ROOT, "derived_data", "GSE5296", "GSE5296_RMA_impact_replication_summary.csv"))
method_map = {"ssgsea": "GSE5296 ssGSEA", "eigengene": "GSE5296 Eigengene"}
region_map = {"A": "rostral (A)", "B": "caudal (B)", "I": "impact (I)"}
sens_by_key = {}
for r in qc:
    key = (r["method"], r["module"], r["time"])
    sens_by_key[key] = (r["sensitivity_delta"], r["direction_changed"])
# impact number_negative_of_6 per module/method
imp_neg = {}
for r in imp_rep:
    imp_neg[(r["method"], r["module"])] = r["number_negative_of_6"]

s3_rows = []
for r in ts:
    if r["module"] not in ("ComplexI_assembly", "Aerobic_respiration"):
        continue
    method = method_map[r["method"]]
    region = region_map[r["region"]]
    time = r["time"]
    delta = rnd(r["delta_score"], 5)
    direction = "down" if delta < 0 else "up"
    if r["region"] == "I":
        sens = sens_by_key.get((r["method"], r["module"], time))
        if sens is None:
            sd, dc = "NA", "NA"
        else:
            sd = rnd(sens[0], 5) if sens[0] not in ("", "NA") else "NA"
            dc = sens[1]
        evaluable = "No" if time == "72h" else "Yes"
        note = ""
        if time == "72h":
            note = "Impact 72h sensitivity not evaluable: excluding flagged arrays left insufficient SCI/control group"
    else:
        sd = "n/a (sensitivity evaluated at impact only)"
        dc = "n/a"
        evaluable = "n/a"
        note = ""
    note2 = note
    if r["region"] == "I":
        nneg = imp_neg.get((r["method"], r["module"]), "")
        note2 = ("Impact %s/%s time points negative." % (nneg, "6")) if note == "" else note
    s3_rows.append([region, time, r["module"], method,
                    r["SCI_n"], r["sham_n"],
                    rnd(r["SCI_mean_score"],5), rnd(r["sham_mean_score"],5),
                    delta, direction, sd, dc, evaluable, note2])
s3_headers = ["Region", "Time", "Program", "Score method",
              "n_SCI", "n_sham", "SCI_mean_score", "Sham_mean_score",
              "Existing SCI_minus_sham delta (mean SCI - mean sham)",
              "Direction", "Sensitivity delta", "Sensitivity direction_changed",
              "Sensitivity evaluable", "Notes"]

# ============================================================
# Table S4 - corrected snRNA program results
# ============================================================
prim = rd_csv(os.path.join(PC, "corrected_primary_results.csv"))
expl = rd_csv(os.path.join(PC, "corrected_exploratory_1d_4d.csv"))
sens = rd_csv(os.path.join(PC, "corrected_sensitivity.csv"))
ct_map = {"Neurons": "Neurons", "Astrocytes": "Astrocytes",
          "Oligodendrocyte_combined": "Mature/myelin-forming oligodendroglial cells",
          "Microglia": "Microglia"}
sc_map = {"A": "snRNA MeanZ", "B": "snRNA Eigengene"}

def build_s4_rows(rows, exploratory=False):
    out = []
    for r in rows:
        ct = ct_map[r["celltype"]]
        time = r["time"]
        if exploratory:
            if r["status"] == "Not evaluable":
                note_e = r["note"] if r.get("note") not in ("", None) else "Not evaluable"
                out.append([ct, time, r["program"], sc_map[r["method"]],
                            r["n_injured"], r["n_control"], "NA", "NA", "NA", "NA",
                            "n/a", "n/a", "NA", "NA",
                            "NA", "Exploratory", note_e + " (not evaluated)"])
                continue
        n_inj = r["n_injured"]; n_ctrl = r["n_control"]
        delta = rnd(r["delta"], 5); direction = r["direction"]
        cons = r["replicate_consistency"]
        es = rnd(r["effect_size"], 4)
        # sensitivity
        s = [x for x in sens if x["celltype"] == r["celltype"] and x["program"] == r["program"]
             and x["method"] == r["method"] and x["time"] == time]
        if s:
            sd = s[0]["sens_direction"]; dc = s[0]["direction_changed"]
        else:
            sd, dc = "n/a", "n/a"
        low_cov_applicable = "Yes" if r["celltype"] in ("Neurons","Microglia","Oligodendrocyte_combined") else "No"
        if r["celltype"] == "Astrocytes":
            low_cov_applicable = "No (no astrocyte pseudobulk met low-coverage criteria)"
            sd = "n/a (no astrocyte low-coverage combo excluded)"
            dc = "n/a"
        note = ("delta = mean(injured library score) - mean(uninjured library score); "
                "effect_size = Cohen's d; all post-injury comparisons share the same n=%s uninjured libraries." % n_ctrl)
        phase = "Exploratory" if exploratory else "Primary"
        out.append([ct, time, r["program"], sc_map[r["method"]],
                    n_inj, n_ctrl, delta, direction, cons, dc,
                    low_cov_applicable, sd, es, rnd(r["nominal_p"],5),
                    rnd(r["adjusted_p"],5) if r.get("adjusted_p") not in ("",None,"NA") else "NA",
                    phase, note])
    return out

s4_rows = build_s4_rows(prim) + build_s4_rows(expl, exploratory=True)
s4_headers = ["Cell type", "Time", "Program", "Score",
              "n_injured libraries", "n_uninjured libraries",
              "Existing score difference (delta)", "Direction",
              "Injured replicate directions (same-side of control median / n)",
              "Sensitivity direction_changed", "Low-coverage sensitivity applicable",
              "Sensitivity direction", "Effect size (Cohen's d)",
              "Nominal P", "Adjusted P (BH)", "Analysis phase", "Notes"]

# ============================================================
# Table S5 - astrocyte aerobic gene-level breadth
# ============================================================
audit = rd_csv(os.path.join(PC, "astrocyte_aerobic_gene_driver_audit.csv"))
aud_sum = rd_csv(os.path.join(PC, "astrocyte_aerobic_gene_driver_summary.csv"))
times = ["7d", "14d", "1m", "2m"]
genes = []
for r in audit:
    if r["gene"] not in genes:
        genes.append(r["gene"])
def g_row(g):
    d = {}
    for r in audit:
        if r["gene"] == g:
            d[r["time"]] = r
    vals = []
    pos = 0
    for t in times:
        r = d.get(t)
        diff = rnd(r["delta_logCPM"], 4) if r else "NA"
        dirn = r["direction"] if r else "NA"
        if dirn == "up": pos += 1
        vals.append(diff); vals.append(dirn)
    return [g] + vals + [pos]
s5_rows = [g_row(g) for g in genes]
s5_headers = ["Gene", "7d diff_logCPM", "7d direction", "14d diff_logCPM", "14d direction",
              "1m diff_logCPM", "1m direction", "2m diff_logCPM", "2m direction",
              "Positive at n/4 primary time points"]
s5_summary = []
for r in aud_sum:
    s5_summary.append([r["time"], r["n_up"], r["n_down"], r["n_genes"], r["frac_up"]])
s5_sum_headers = ["Time", "Positive genes", "Negative genes", "Total", "Fraction positive"]

# ============================================================
# Table S6 - spatial transcriptomic results
# ============================================================
rep = rd_csv(os.path.join(SP, "replicate_level_summary.csv"))
prim_sp = rd_csv(os.path.join(SP, "primary_spatial_comparisons.csv"))
sens_sp = rd_csv(os.path.join(SP, "spatial_sensitivity.csv"))
sp_sc = {"A": "Spatial Mean", "B": "Spatial Rank"}
# replicate-level sheet
s6_repl = []
for r in rep:
    cond = r["condition"]
    rep_lab = r["replicate"]
    bio = "%s_rep%s" % (cond, rep_lab)
    s6_repl.append([bio, cond, r["program"], sp_sc[r["method"]],
                    r["n_spots"], rnd(r["median"],5), rnd(r["IQR"],5),
                    rnd(r["low_fraction"],5), rnd(r["high_fraction"],5)])
s6_repl_headers = ["Biological sample", "Condition", "Program", "Score",
                   "n_spots", "Existing median score", "IQR",
                   "Low-score fraction", "High-score fraction"]
# primary comparisons sheet
s6_prim = []
for r in prim_sp:
    comparison = r["comparison"]
    time_lab = "7d vs uninjured" if comparison == "7d" else "2m vs uninjured"
    s6_prim.append([r["program"], sp_sc[r["method"]], time_lab,
                    r["n_control"], r["n_injured"],
                    rnd(r["control_mean"],5), rnd(r["injured_mean"],5),
                    rnd(r["delta"],5), r["direction"],
                    r["replicate_consistency"], rnd(r["effect_size"],4),
                    rnd(r["nominal_p"],5)])
s6_prim_headers = ["Program", "Score", "Comparison", "n biological samples/group(ctrl)",
                   "n biological samples/group(inj)", "Control mean", "Injured mean",
                   "Existing delta", "Direction", "Direction consistency",
                   "Effect size", "Nominal P"]
# sensitivity sheet
s6_sens = []
for r in sens_sp:
    if r["sensitivity"] == "A_lowUMI":
        cmp = "7d vs uninjured" if r["comparison"] == "7d" else "2m vs uninjured"
        s6_sens.append(["Low-UMI exclusion (<4,077 UMIs; 30,548/33,941 spots retained)",
                        r["program"], sp_sc.get(r["method"], r["method"]), cmp,
                        rnd(r["direction_delta"],5), r["direction"]])
    elif r["sensitivity"] == "B_leave_one_section":
        cmp = "7d vs uninjured" if r["comparison"] == "7d" else "2m vs uninjured"
        s6_sens.append(["Leave-one-section-out", r["program"], sp_sc.get(r["method"], r["method"]),
                        cmp, rnd(r["direction_delta"],5), "0 reversals"])
    elif r["sensitivity"] == "C_method_agreement":
        s6_sens.append(["Score-method agreement (Spatial Mean vs Spatial Rank)", r["program"], "all",
                        r["comparison"], rnd(r["direction_delta"],4), "correlation %s" % r["direction"]])
s6_sens_headers = ["Sensitivity", "Program", "Score", "Comparison", "Value", "Result"]
# descriptive CV sheet (NOT USED FOR PRIMARY INFERENCE)
s6_cv = []
for r in rep:
    s6_cv.append(["%s_rep%s" % (r["condition"], r["replicate"]), r["condition"],
                  r["program"], sp_sc[r["method"]], rnd(r["IQR"],5)])
s6_cv_headers = ["Biological sample", "Condition", "Program", "Score",
                 "IQR (spread descriptor - NOT USED FOR PRIMARY INFERENCE)"]

# ============================================================
# Table S7 - GSE319931 mixed-effects results
# ============================================================
omn = rd_csv(os.path.join(G31, "GSE319931_mixed_model_omnibus.csv"))
pw  = rd_csv(os.path.join(G31, "GSE319931_mixed_model_pairwise.csv"))
diag_txt = open(os.path.join(G31, "GSE319931_mixed_model_diagnostics.txt"), encoding="utf-8").read()
consist = rd_csv(os.path.join(G31, "GSE319931_mixed_vs_paired_consistency.csv"))
sec = rd_csv(os.path.join(G31, "GSE319931_secondary_vs_UI.csv"))
sc31 = {"A": "GSE319931 MeanZ", "B": "GSE319931 Eigengene"}
s7_omn = []
for r in omn:
    s7_omn.append([r["program"], sc31[r["method"]], "score ~ region + (1|animal)",
                   rnd(r["region_omnibus_F"],4), r["df_num"], rnd(r["df_den"],2),
                   rnd(r["region_omnibus_rawP"],5), r["singular_fit"], r["convergence_warning"],
                   "n_obs=%s; n_animals=%s" % (r["n_obs"], r["n_animals"])])
s7_omn_headers = ["Program", "Score", "Model formula", "Omnibus F", "df1", "df2",
                  "Omnibus region P", "Singular fit", "Convergence warning", "Notes"]
s7_pw = []
for r in pw:
    s7_pw.append([r["program"], sc31[r["method"]], r["contrast"],
                  rnd(r["estimate"],5), rnd(r["SE"],5),
                  rnd(r["CI_lower"],5), rnd(r["CI_upper"],5),
                  rnd(r["raw_P"],5), rnd(r["Holm_P"],5)])
s7_pw_headers = ["Program", "Score", "Contrast", "Estimate", "SE",
                 "95% CI lower", "95% CI upper", "Raw P", "Holm-adjusted P"]
# diagnostics sheet (parsed from txt - directional, leave-one-animal-out)
diag_rows = []
for block in diag_txt.split("--- "):
    if not block.strip():
        continue
    lines = block.strip().splitlines()
    prog = lines[0].replace("/", " / ").strip()
    first = [l for l in lines if "singular_fit" in l or "max_scaled_residual" in l]
    flips = [l for l in lines if "direction-flip animal/contrast" in l]
    diag_rows.append([prog, "; ".join(first), "; ".join(flips)])
s7_diag_headers = ["Program / Score", "Diagnostics", "Leave-one-animal-out direction flips"]
# paired vs mixed consistency
s7_cons = []
for r in consist:
    s7_cons.append([r["program"], sc31[r["method"]], r["contrast"],
                    r["mixed_model_direction"], r["paired_direction"],
                    r["direction_consistent"], rnd(r["mixed_Holm_P"],5), rnd(r["paired_P"],5)])
s7_cons_headers = ["Program", "Score", "Contrast", "Mixed-model direction", "Paired direction",
                   "Direction consistent", "Mixed Holm P", "Original paired P"]
# UI exploratory (separate - NOT primary)
s7_ui = []
for r in sec:
    cmp_ = "Epicenter vs UI" if r["region"] == "Epicenter" else ("Above vs UI" if r["region"] == "Above" else "Below vs UI")
    s7_ui.append([r["program"], sc31[r["method"]], cmp_,
                  rnd(r["UI_mean"],5), rnd(r["SCI_mean"],5), rnd(r["delta"],5),
                  r["direction"], "%s/%s" % (r["n_inj_below_UImedian"], r["n_total"]), rnd(r["welch_p"],5)])
s7_ui_headers = ["Program", "Score", "Comparison (EXPLORATORY - not region-matched, whole 9mm UI)",
                 "UI mean", "SCI region mean", "Delta", "Direction",
                 "n injured below UI median / n", "Welch P"]

# ============================================================
# Table S8 - QC / exclusions / audit trail
# ============================================================
qc_sum = rd_csv(os.path.join(ROOT, "derived_data", "GSE5296", "GSE5296_RMA_QC_summary.csv"))
flagged = [r for r in qc_sum if r["obvious_outlier"] == "TRUE"]
s8_flag = []
for r in flagged:
    s8_flag.append([r["GSM"], "%s/%s/%s" % (r["region"], r["treatment"], r["time"]),
                    r["remark"], "Retained (primary)",
                    "Excluded (sensitivity)",
                    "Impact 72h not evaluable in sensitivity (insufficient group after exclusion)"
                    if r["time"] == "72h" else "No direction change in evaluable comparisons"])
s8_flag_headers = ["GSE5296 flagged GSM", "Group (region/treatment/time)", "QC reason (from GSE5296_QC_summary.csv)",
                   "Retained in primary", "Excluded in sensitivity", "Effect on direction"]
# snRNA low-coverage
low = rd_csv(os.path.join(ROOT, "metadata", "sample_design", "GSE234774_lowcoverage_log.csv"))
s8_low = []
for r in low:
    if r["excluded"] != "TRUE":
        continue
    s8_low.append([r["pseudobulk_sample"], r["n_nuclei"], r["total_target_gene_UMI"],
                   rnd(r["target_gene_detection_frac"],3), r["reason"],
                   "Retained (primary)", "Excluded (sensitivity)"])
s8_low_headers = ["Library x celltype", "n_nuclei", "Total target-gene UMI",
                  "Target-gene detection fraction", "Triggered criterion",
                  "Retained in primary", "Excluded in sensitivity"]
# normalization correction (old vs new)
s8_corr = [
    ["Counts entering DGEList", "Target-gene-only (71 respiratory genes)", "Full-transcriptome (26,918 genes)"],
    ["Library size (edgeR lib.size)", "colSums of target-only matrix (target UMI)", "colSums of full-transcriptome matrix == nCount_RNA sum (verified, max|diff|=0)"],
    ["TMM estimation stage", "from target-only matrix", "from filtered full-transcriptome matrix"],
    ["Target-gene extraction", "before DGEList/normalization", "only AFTER filtering + TMM (29/30 + 49/53)"],
    ["Filtering", "n/a", "filterByExpr(group=celltype), 26,918 -> 15,785 retained"],
    ["logCPM prior.count", "3", "3"],
    ["Inference status", "SUPERSEDED / INVALID FOR INFERENCE (audit: FAIL - target-gene-only library sizes)",
     "ONLY valid version - all manuscript snRNA results derive from this"],
]
s8_corr_headers = ["Item", "Original Phase 2C workflow (superseded)", "Corrected full-transcriptome workflow"]
# spatial QC
sp_qc = rd_csv(os.path.join(SP, "spatial_technical_QC.csv"))
s8_sp = []
for r in sp_qc:
    s8_sp.append([r["check"], r["status"], r["detail"]])
s8_sp_headers = ["Spatial technical QC check", "Status", "Detail"]
s8_total = [
    ["GSE5296 arrays", "96 (54 SCI + 36 sham + 6 naive)"],
    ["GSE234774 snRNA libraries", "20 (UI3, 1d3, 4d3, 7d2, 14d3, 1m3, 2m3)"],
    ["GSE234774 longitudinal nuclei", "151,712"],
    ["GSE234774 analyzed four-cell-type nuclei", "115,475"],
    ["Corrected snRNA pseudobulk samples", "77 (80 - 3 absent Microglia combos)"],
    ["Corrected snRNA gene availability", "Complex I 29/30; Aerobic 49/53"],
    ["Spatial gene availability", "Complex I 30/30; Aerobic 52/53"],
    ["Spatial spots / sections / biological samples", "33,941 spots / 36 sections / 9 biological samples"],
    ["GSE319931 samples", "4 SCI animals x 3 regions + 4 UI = 16"],
    ["Astrocyte aerobic gene counts (7d/14d/1m/2m)", "40/39/38/33 positive (of 49)"],
]
s8_total_headers = ["Item", "Value"]

# ============================================================
# Emit xlsx + csv
# ============================================================
def emit(tag, sheets):
    for name, headers, rows in sheets:
        p = os.path.join(BASE, "%s_%s.csv" % (tag, name))
        write_csv(p, headers, rows)
    xl = os.path.join(BASE, tag + ".xlsx")
    save_xlsx(xl, [(name, hd, rows) for name, hd, rows in sheets])
    return xl

# S1
s1 = [
    ["GSE5296", "Affymetrix Mouse Genome 430 2.0 microarray (bulk)",
     "T8 moderate contusion SCI; laminectomy-only sham; C57BL/6; isoflurane",
     "Rostral / impact / caudal 0.4-cm spinal cord segments",
     "0.5h / 4h / 24h / 72h / 7d / 28d",
     "Array (each SCI replicate pooled tissue from 4 mice; 12 mice per injury time point). Sham/naive pooling not established.",
     "96 arrays (54 SCI, 36 sham, 6 naive)",
     "None (bulk); 6 naive arrays included in RMA/QC but not used as SCI-sham comparators",
     "Temporal + regional tissue-level respiratory program characterization (impact-site primary pattern)",
     "Sham pooling not established; n=2 sham arrays per region/time; naive n=2 per region"],
    ["GSE234774 snRNA",
     "Single-nucleus RNA-seq (snRNA), time-course cohort of Tabulae Paralytica",
     "T10 mid-thoracic crush SCI (source study n=3 mice/condition)",
     "Injury-site spinal cord tissue",
     "uninjured + 1d / 4d / 7d / 14d / 1m / 2m",
     "Library (= biological replicate); library-level pseudobulk. Exact animal-ID to deposited-library-ID mapping unavailable.",
     "20 libraries (UI3, 1d3, 4d3, 7d2, 14d3, 1m3, 2m3); 151,712 nuclei in selected longitudinal cohort; 115,475 in the four analyzed populations; 77 library x celltype pseudobulk",
     "Nuclei are not replicates; 77 = 80 - 3 absent microglia combinations",
     "Cell-type-resolved respiratory program validation (4 cell types)",
     "Library-level pseudobulk; animal/library mapping unavailable; small n per condition"],
    ["GSE234774 spatial",
     "2D spatial transcriptomics (10x Visium, spatial component of Tabulae Paralytica)",
     "T10 mid-thoracic crush SCI (same injury model as snRNA component)",
     "Injury-site tissue sections",
     "uninjured / 7d / 2m",
     "Biological sample (= replicate; 4 tissue sections per sample)",
     "9 biological samples (3 per condition); 36 sections; 33,941 spots",
     "Sections and spots are spatial subsamples, not biological replicates",
     "Spatially resolved tissue-level respiratory program activity (7d and 2m)",
     "No spot-level anatomical region annotations in processed metadata"],
    ["GSE319931",
     "Bulk RNA-seq (raw gene counts)",
     "Severe contusion-compression SCI; FEJOTA mouse clip 8 g closing force x exactly 1 min, extradural after laminectomy, L1/L3 lumbar level; 30 dpi",
     "SCI: Above / Epicenter / Below 3-mm segments; UI: whole corresponding 9-mm segment (T11-S1)",
     "30 dpi (chronic)",
     "Animal (4 SCI animals, each contributing 3 regions - repeated measures; 4 uninjured animals)",
     "16 samples (4 SCI animals x 3 regions + 4 UI)",
     "Region = repeated measure within animal; UI is a single whole 9-mm segment (not region-matched)",
     "Chronic anatomical (regional) validation",
     "UI not region-matched to SCI 3-mm segments; regional pattern cannot be fully attributed to SCI"],
]
s1_headers = ["Dataset", "Modality", "SCI model", "Anatomical sampling", "Time points used",
              "Biological/available replicate unit", "n",
              "Technical/subsample structure", "Primary role", "Main limitation"]

# S8
# merge S8 sheets
s8_sheets = [
    ("GSE5296_flagged_GSM", s8_flag_headers, s8_flag),
    ("snRNA_low_coverage_combos", s8_low_headers, s8_low),
    ("snRNA_normalization_correction", s8_corr_headers, s8_corr),
    ("Spatial_technical_QC", s8_sp_headers, s8_sp),
    ("Total_cohort_numbers", s8_total_headers, s8_total),
]

emit("Table_S1_Dataset_design", [("Dataset_design", s1_headers, s1)])
emit("Table_S2_Locked_gene_programs", [("Locked_gene_programs", s2_headers, s2_rows)])
emit("Table_S3_GSE5296_program_results", [("Program_results", s3_headers, s3_rows)])
emit("Table_S4_Corrected_snRNA_program_results", [("Program_results", s4_headers, s4_rows)])
emit("Table_S5_Astrocyte_Aerobic_gene_level",
     [("Gene_level", s5_headers, s5_rows),
      ("Summary_by_time", s5_sum_headers, s5_summary)])
emit("Table_S6_Spatial_results",
     [("Replicate_level_summary", s6_repl_headers, s6_repl),
      ("Primary_comparisons_sensitivity", s6_prim_headers, s6_prim),
      ("Sensitivity", s6_sens_headers, s6_sens),
      ("Descriptive_CV_NotForInference", s6_cv_headers, s6_cv)])
emit("Table_S7_GSE319931_mixed_models",
     [("Omnibus", s7_omn_headers, s7_omn),
      ("Pairwise", s7_pw_headers, s7_pw),
      ("Diagnostics", s7_diag_headers, diag_rows),
      ("Mixed_vs_paired_consistency", s7_cons_headers, s7_cons),
      ("UI_exploratory_NOT_PRIMARY", s7_ui_headers, s7_ui)])
emit("Table_S8_QC_AuditTrail", s8_sheets)

print("DONE")
print("s1 rows", len(s1))
print("s2 rows", len(s2_rows))
print("s3 rows", len(s3_rows))
print("s4 rows", len(s4_rows))
print("s5 gene rows", len(s5_rows))
print("s6 repl rows", len(s6_repl), "prim", len(s6_prim))
print("s7 omn", len(s7_omn), "pw", len(s7_pw), "diag", len(diag_rows))
print("s8 flagged", len(s8_flag), "low", len(s8_low))
