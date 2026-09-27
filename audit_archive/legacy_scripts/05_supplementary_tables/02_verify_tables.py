# -*- coding: utf-8 -*-
"""Mechanical consistency audit of generated Supplementary tables against
Manuscript v0.4 locked numbers. Extraction verification only - no recomputation
of the underlying analysis."""
import csv, os
BASE = os.path.join(os.environ.get("SCI_RESULTS_ROOT", os.getcwd()), "supplementary_source_tables")
def rd(p):
    with open(p, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

results = []
def check(name, cond, detail=""):
    results.append((name, "PASS" if cond else "FAIL", detail))

# --- S2 ---
s2 = rd(os.path.join(BASE, "Table_S2_Locked_gene_programs_Locked_gene_programs.csv"))
genes = [r["Gene symbol"] for r in s2]
ci = [r for r in s2 if "ComplexI_assembly" in r["Program(s)"]]
ae = [r for r in s2 if "Aerobic_respiration" in r["Program(s)"]]
both = [r for r in s2 if "ComplexI_assembly" in r["Program(s)"] and "Aerobic_respiration" in r["Program(s)"]]
check("S2 30 ComplexI genes", len(ci) == 30, "n=%d" % len(ci))
check("S2 53 Aerobic genes", len(ae) == 53, "n=%d" % len(ae))
check("S2 11 shared", len(both) == 11, "n=%d" % len(both))
check("S2 72 union", len(genes) == 72, "n=%d" % len(genes))
check("S2 0 mtDNA-encoded", all(r["Nuclear/mtDNA encoded"] == "Nuclear" for r in s2), "all nuclear")
# snRNA availability
sn_missing = {r["Gene symbol"] for r in s2 if "No(" in r["Available corrected snRNA"]}
check("S2 snRNA ComplexI 29/30 (Ndufaf3 filtered)", "Ndufaf3" in sn_missing and len([r for r in s2 if "No(" in r["Available corrected snRNA"] and "ComplexI" in r["Program(s)"]]) == 1,
      "missing=%s" % sorted(sn_missing))
sp_missing = {r["Gene symbol"] for r in s2 if "No(" in r["Available spatial"]}
check("S2 spatial 52/53 (Mup4 unavailable)", sp_missing == {"Mup4"}, "missing=%s" % sorted(sp_missing))

# --- S3 ---
s3 = rd(os.path.join(BASE, "Table_S3_GSE5296_program_results_Program_results.csv"))
imp = [r for r in s3 if r["Region"] == "impact (I)"]
for meth in ("GSE5296 ssGSEA", "GSE5296 Eigengene"):
    for prog in ("ComplexI_assembly", "Aerobic_respiration"):
        rows = [r for r in imp if r["Program"] == prog and r["Score method"] == meth]
        neg = sum(1 for r in rows if r["Direction"] == "down")
        check("S3 impact %s %s 6/6 down" % (prog, meth), neg == 6, "n_down=%d" % neg)
s3_72h = [r for r in imp if r["Time"] == "72h"]
check("S3 impact 72h sensitivity not evaluable", all(r["Sensitivity evaluable"] == "No" for r in s3_72h), "n=%d" % len(s3_72h))

# --- S4 ---
s4 = rd(os.path.join(BASE, "Table_S4_Corrected_snRNA_program_results_Program_results.csv"))
ast = [r for r in s4 if r["Cell type"] == "Astrocytes" and r["Analysis phase"] == "Primary" and r["Program"] == "Aerobic_respiration"]
for sc in ("snRNA MeanZ", "snRNA Eigengene"):
    rows = [r for r in ast if r["Score"] == sc]
    up = sum(1 for r in rows if r["Direction"] == "up")
    check("S4 Astro/Aerobic %s 4/4 primary up" % sc, up == 4, "n_up=%d" % up)
s4_primary = [r for r in s4 if r["Analysis phase"] == "Primary"]
check("S4 64 primary rows", len(s4_primary) == 64, "n=%d" % len(s4_primary))
s4_expl = [r for r in s4 if r["Analysis phase"] == "Exploratory"]
check("S4 32 exploratory rows", len(s4_expl) == 32, "n=%d" % len(s4_expl))
check("S4 Microglia 1d not evaluable x4", len([r for r in s4_expl if r["Cell type"] == "Microglia" and r["Time"] == "1d"]) == 4, "n=%d" % len([r for r in s4_expl if r["Cell type"] == "Microglia" and r["Time"] == "1d"]))

# --- S5 ---
s5s = rd(os.path.join(BASE, "Table_S5_Astrocyte_Aerobic_gene_level_Summary_by_time.csv"))
expect = {"7d": 40, "14d": 39, "1m": 38, "2m": 33}
for r in s5s:
    check("S5 %s positive count %d" % (r["Time"], expect[r["Time"]]), int(r["Positive genes"]) == expect[r["Time"]],
          "n_up=%s" % r["Positive genes"])
# verify Ndufs4/Cox10/Ndufs1/Fxn all negative
s5g = rd(os.path.join(BASE, "Table_S5_Astrocyte_Aerobic_gene_level_Gene_level.csv"))
four = {g: s5g for g in ("Ndufs4","Cox10","Ndufs1","Fxn")}
for g in ("Ndufs4","Cox10","Ndufs1","Fxn"):
    row = [r for r in s5g if r["Gene"] == g][0]
    check("S5 %s negative at all 4 primary times" % g, row["Positive at n/4 primary time points"] == "0",
          "pos=%s" % row["Positive at n/4 primary time points"])

# --- S6 ---
s6p = rd(os.path.join(BASE, "Table_S6_Spatial_results_Primary_comparisons_sensitivity.csv"))
check("S6 8 primary comparisons", len(s6p) == 8, "n=%d" % len(s6p))
check("S6 all 8 down", all(r["Direction"] == "down" for r in s6p), "n_down=%d" % sum(1 for r in s6p if r["Direction"]=="down"))
check("S6 3/3 consistency all", all(r["Direction consistency"].strip() == "3/3" for r in s6p),
      "consistency values=%s" % {r["Program"]+"_"+r["Score"]+"_"+r["Comparison"]: r["Direction consistency"] for r in s6p})
s6sens = rd(os.path.join(BASE, "Table_S6_Spatial_results_Sensitivity.csv"))
loos = [r for r in s6sens if r["Sensitivity"] == "Leave-one-section-out"]
check("S6 leave-one-section-out 0 reversals", all(r["Result"] == "0 reversals" for r in loos), "n=%d" % len(loos))

# --- S7 ---
s7o = rd(os.path.join(BASE, "Table_S7_GSE319931_mixed_models_Omnibus.csv"))
omnexp = {("ComplexI_assembly","GSE319931 MeanZ"): 0.0151,
          ("ComplexI_assembly","GSE319931 Eigengene"): 0.0071,
          ("Aerobic_respiration","GSE319931 MeanZ"): 0.0197,
          ("Aerobic_respiration","GSE319931 Eigengene"): 0.0183}
for (p, sc_), exp in omnexp.items():
    row = [r for r in s7o if r["Program"] == p and r["Score"] == sc_][0]
    check("S7 omnibus %s %s P" % (p, sc_), abs(float(row["Omnibus region P"]) - exp) < 0.0005, "P=%s" % row["Omnibus region P"])
s7p = rd(os.path.join(BASE, "Table_S7_GSE319931_mixed_models_Pairwise.csv"))
holmab = {("ComplexI_assembly","GSE319931 MeanZ"): 0.0157,
          ("ComplexI_assembly","GSE319931 Eigengene"): 0.0089,
          ("Aerobic_respiration","GSE319931 MeanZ"): 0.0279,
          ("Aerobic_respiration","GSE319931 Eigengene"): 0.0194}
for (p, sc_), exp in holmab.items():
    row = [r for r in s7p if r["Program"] == p and r["Score"] == sc_ and r["Contrast"] == "A-B"][0]
    check("S7 A-B Holm %s %s" % (p, sc_), abs(float(row["Holm-adjusted P"]) - exp) < 0.0005, "P=%s" % row["Holm-adjusted P"])

# --- cohort numbers (S1 + S8) ---
s8t = rd(os.path.join(BASE, "Table_S8_QC_AuditTrail_Total_cohort_numbers.csv"))
d = {r["Item"]: r["Value"] for r in s8t}
check("Cohort 96 GSE5296 arrays", "96 (54 SCI + 36 sham + 6 naive)" in d["GSE5296 arrays"])
check("Cohort 20 snRNA libraries", d["GSE234774 snRNA libraries"].startswith("20 "))
check("Cohort 151,712 nuclei", d["GSE234774 longitudinal nuclei"] == "151,712")
check("Cohort 115,475 four-cell-type nuclei", d["GSE234774 analyzed four-cell-type nuclei"] == "115,475")
check("Cohort 77 pseudobulk", d["Corrected snRNA pseudobulk samples"].startswith("77"))
check("Cohort spatial 33941/36/9", d["Spatial spots / sections / biological samples"].startswith("33,941 spots / 36 sections / 9"))
check("Cohort GSE319931 4x3+4", d["GSE319931 samples"].startswith("4 SCI animals x 3 regions + 4 UI"))

print("")
print("CONSISTENCY AUDIT RESULTS")
nf = 0
for name, status, detail in results:
    print("[%s] %s  %s" % ("PASS" if status else "FAIL", name, detail))
    if status == "FAIL": nf += 1
print("")
print("TOTAL: %d checks, %d FAIL" % (len(results), nf))
