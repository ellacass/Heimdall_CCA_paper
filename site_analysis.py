import math
import pandas as pd
from collections import Counter
from Bio import SeqIO

###########################################################
# -------------------- USER CONFIG -----------------------
###########################################################

FASTA         = "aln_0.5_no_termini.fasta"   # your trimmed alignment
CORE_START    = 1                               # first site of catalytic core
CORE_END      = 205                             # last site (full trimmed alignment)
OUTPUT_ALL    = "residue_pattern_audit.csv"
OUTPUT_SYN_E  = "synapomorphy_candidates_euks.csv"
OUTPUT_SYN_B  = "synapomorphy_candidates_bact.csv"
OUTPUT_SUMM   = "residue_pattern_summary.txt"

# Taxon filters — adjust prefixes to match your sequence IDs
HEIM_PREFIX   = "Heim12"
EB_PREFIX     = "EB_"
BACT_TARGETS  = ["Bacteroidetes", "Rhodothermales", "Salinibacter", "Chlorobi"]

# Hardcoded coordinates mapped exactly to your 205-site trimmed core format
# (Note: Motif B was filtered out during ClipKIT trimming sequence steps)
MOTIFS_TRIMMED = {
    "A": (4, 8),
    "FL": (17, 19),
    "C": (86, 92),
    "D": (133, 137),
    "E": (172, 177)
}

# Distance threshold (in amino acids) to classify a flanking position as "Proximal"
PROXIMITY_THRESHOLD = 5

###########################################################
# -------------------- UTILITIES -------------------------
###########################################################

def load_seqs(fasta):
    return {r.id: str(r.seq) for r in SeqIO.parse(fasta, "fasta")}


def dominant_residue(seq_dict, ids, site):
    aas = [seq_dict[t][site - 1] for t in ids if seq_dict[t][site - 1] != "-"]
    if not aas:
        return None, 0.0
    counter = Counter(aas)
    aa, count = counter.most_common(1)[0]
    return aa, count / len(aas)


def keff(seq_dict, ids, site):
    aas = [seq_dict[t][site - 1] for t in ids if seq_dict[t][site - 1] != "-"]
    if not aas:
        return 1.0
    total = len(aas)
    counts = Counter(aas)
    entropy = -sum((c / total) * math.log(c / total) for c in counts.values())
    return math.exp(entropy)


ALIPHATIC   = set("GAVLIP")
AROMATIC    = set("FYWH")
POLAR       = set("STNQ")
CHARGED_POS = set("RK")
CHARGED_NEG = set("DE")
SPECIAL     = set("CM")

def biochemical_class(aa):
    if aa in ALIPHATIC:   return "aliphatic"
    if aa in AROMATIC:    return "aromatic"
    if aa in POLAR:       return "polar"
    if aa in CHARGED_POS: return "pos_charged"
    if aa in CHARGED_NEG: return "neg_charged"
    if aa in SPECIAL:     return "special"
    return "other"


def is_conservative(aa1, aa2):
    return biochemical_class(aa1) == biochemical_class(aa2)

###########################################################
# -------------------- MAIN AUDIT ------------------------
###########################################################

def run_audit(fasta, core_start, core_end):
    seqs = load_seqs(fasta)
    n_sites = len(next(iter(seqs.values())))

    heim_ids = [x for x in seqs if HEIM_PREFIX in x]
    eb_ids   = [x for x in seqs if x.startswith(EB_PREFIX)]
    bact_ids = [x for x in seqs if any(t in x for t in BACT_TARGETS)]

    print(f"Loaded {len(seqs)} sequences, {n_sites} alignment sites")
    print(f"  Heim12 taxa : {len(heim_ids)}")
    print(f"  Eukaryote   : {len(eb_ids)}")
    print(f"  Bacteria    : {len(bact_ids)}")
    print()

    records = []

    for site in range(1, n_sites + 1):
        region = "CORE" if core_start <= site <= core_end else "TAIL"

        # Calculate structural motif metrics
        motif_status = "Outside"
        closest_motif = None
        min_dist = float('inf')
        
        for name, (start, end) in MOTIFS_TRIMMED.items():
            if start <= site <= end:
                motif_status = "Inside"
                closest_motif = name
                min_dist = 0
                break
            else:
                if site < start:
                    dist = start - site
                else:
                    dist = site - end
                if dist < min_dist:
                    min_dist = dist
                    closest_motif = name
                    
        if motif_status != "Inside" and min_dist <= PROXIMITY_THRESHOLD:
            motif_status = "Proximal"

        h_aa, h_freq = dominant_residue(seqs, heim_ids, site)
        e_aa, e_freq = dominant_residue(seqs, eb_ids,   site)
        b_aa, b_freq = dominant_residue(seqs, bact_ids, site)

        if None in (h_aa, e_aa, b_aa):
            interp = "GAPPY"
            conservative = None
            biochem_note = ""
        elif h_aa == e_aa == b_aa:
            interp = "CONSERVED"
            conservative = True
            biochem_note = "fully conserved"
        elif h_aa == e_aa and h_aa != b_aa:
            interp = "SUPPORT_EUKS"
            conservative = is_conservative(h_aa, b_aa)
            biochem_note = (
                "conservative substitution (weak synapomorphy)"
                if conservative else
                "biochemically distinct (strong synapomorphy)"
            )
        elif h_aa == b_aa and h_aa != e_aa:
            interp = "SUPPORT_BACT"
            conservative = is_conservative(h_aa, e_aa)
            biochem_note = (
                "conservative substitution (weak synapomorphy)"
                if conservative else
                "biochemically distinct (strong synapomorphy)"
            )
        else:
            interp = "CONFLICT"
            conservative = None
            biochem_note = f"Heim={h_aa}, EB={e_aa}, Bact={b_aa}"

        records.append({
            "Site"           : site,
            "Region"         : region,
            "Motif_Status"   : motif_status,
            "Closest_Motif"  : closest_motif,
            "Motif_Distance" : min_dist,
            "Heim_Residue"   : h_aa,
            "Heim_Freq"      : round(h_freq, 3),
            "EB_Residue"     : e_aa,
            "EB_Freq"        : round(e_freq, 3),
            "Bact_Residue"   : b_aa,
            "Bact_Freq"      : round(b_freq, 3),
            "Heim_Keff"      : round(keff(seqs, heim_ids, site), 3),
            "Interpretation" : interp,
            "Conservative"   : conservative,
            "Biochem_Note"   : biochem_note,
        })

    return pd.DataFrame(records)


###########################################################
# -------------------- REPORTING -------------------------
###########################################################

def write_summary(df, path):
    lines = []
    lines.append("=" * 70)
    lines.append("RESIDUE PATTERN AUDIT — SUMMARY")
    lines.append("=" * 70)

    for region in ["CORE", "TAIL"]:
        sub = df[df["Region"] == region]
        lines.append(f"\n>>> {region} ({len(sub)} sites)")
        lines.append("-" * 70)

        counts = sub["Interpretation"].value_counts()
        lines.append(f"{'Category':<30} {'Count':>6}  {'%':>6}")
        for cat, n in counts.items():
            lines.append(f"  {cat:<28} {n:>6}  {100*n/len(sub):>5.1f}%")

        euk_sites = sub[sub["Interpretation"] == "SUPPORT_EUKS"]
        if not euk_sites.empty:
            strong = euk_sites[euk_sites["Conservative"] == False]
            weak   = euk_sites[euk_sites["Conservative"] == True]
            lines.append(f"\n  SUPPORT_EUKS breakdown:")
            lines.append(f"    Strong synapomorphies (biochemically distinct) : {len(strong)}")
            lines.append(f"    Weak synapomorphies   (conservative change)     : {len(weak)}")

        bact_sites = sub[sub["Interpretation"] == "SUPPORT_BACT"]
        if not bact_sites.empty:
            strong = bact_sites[bact_sites["Conservative"] == False]
            weak   = bact_sites[bact_sites["Conservative"] == True]
            lines.append(f"\n  SUPPORT_BACT breakdown:")
            lines.append(f"    Strong synapomorphies (biochemically distinct) : {len(strong)}")
            lines.append(f"    Weak synapomorphies   (conservative change)     : {len(weak)}")

    # NEW SECTION: Motif Spatial Partitioning Audit Summary
    lines.append("\n" + "=" * 70)
    lines.append("MOTIF STRUCTURAL ENRICHMENT AUDIT (CORE MATRIX)")
    lines.append("=" * 70)
    
    core_df = df[df["Region"] == "CORE"]
    for target in ["SUPPORT_EUKS", "SUPPORT_BACT"]:
        sub_target = core_df[core_df["Interpretation"] == target]
        total_target = len(sub_target)
        lines.append(f"\n>>> {target} Positional Distribution (Total: {total_target} sites)")
        lines.append("-" * 70)
        
        if total_target > 0:
            status_counts = sub_target["Motif_Status"].value_counts()
            for status in ["Inside", "Proximal", "Outside"]:
                cnt = status_counts.get(status, 0)
                pct = (cnt / total_target) * 100
                note = f" (within {PROXIMITY_THRESHOLD} AA of motif boundaries)" if status == "Proximal" else ""
                lines.append(f"  {status:<8}{note:<42} : {cnt:>3} sites ({pct:>5.1f}%)")
        else:
            lines.append("  No sites found matching this classification category.")

    # Eukaryotic Candidate Details
    lines.append("\n" + "=" * 70)
    lines.append("STRONG EUKARYOTIC SYNAPOMORPHY CANDIDATES (CORE, SUPPORT_EUKS, biochemically distinct)")
    lines.append("=" * 70)
    candidates_euk = df[
        (df["Region"] == "CORE") &
        (df["Interpretation"] == "SUPPORT_EUKS") &
        (df["Conservative"] == False)
    ]
    lines.append(f"{'Site':<8} {'Compartment':<12} {'Heim':<6} {'EB':<6} {'Bact':<6} {'Heim class':<15} {'Note'}")
    for _, r in candidates_euk.iterrows():
        loc_str = f"{r.Motif_Status}_{r.Closest_Motif}" if r.Motif_Status != "Outside" else "Scaffold"
        lines.append(
            f"  {int(r.Site):<6} {loc_str:<12} {r.Heim_Residue:<6} {r.EB_Residue:<6} "
            f"{r.Bact_Residue:<6} {biochemical_class(r.Heim_Residue):<15} {r.Biochem_Note}"
        )

    # Bacterial Candidate Details
    lines.append("\n" + "=" * 70)
    lines.append("STRONG BACTERIAL SYNAPOMORPHY CANDIDATES (CORE, SUPPORT_BACT, biochemically distinct)")
    lines.append("=" * 70)
    candidates_bact = df[
        (df["Region"] == "CORE") &
        (df["Interpretation"] == "SUPPORT_BACT") &
        (df["Conservative"] == False)
    ]
    lines.append(f"{'Site':<8} {'Compartment':<12} {'Heim':<6} {'EB':<6} {'Bact':<6} {'Heim class':<15} {'Note'}")
    for _, r in candidates_bact.iterrows():
        loc_str = f"{r.Motif_Status}_{r.Closest_Motif}" if r.Motif_Status != "Outside" else "Scaffold"
        lines.append(
            f"  {int(r.Site):<6} {loc_str:<12} {r.Heim_Residue:<6} {r.EB_Residue:<6} "
            f"{r.Bact_Residue:<6} {biochemical_class(r.Heim_Residue):<15} {r.Biochem_Note}"
        )

    text = "\n".join(lines)
    print(text)
    with open(path, "w") as f:
        f.write(text)
    print(f"\nSummary report written to: {path}")


###########################################################
# ----------------------- MAIN ---------------------------
###########################################################

if __name__ == "__main__":
    df = run_audit(FASTA, CORE_START, CORE_END)

    # Full CSV Output
    df.to_csv(OUTPUT_ALL, index=False)
    print(f"Full audit matrix saved to: {OUTPUT_ALL}")

    # Euk Synapomorphy candidates only
    syn_euk = df[
        (df["Region"] == "CORE") &
        (df["Interpretation"] == "SUPPORT_EUKS") &
        (df["Conservative"] == False)
    ]
    syn_euk.to_csv(OUTPUT_SYN_E, index=False)

    # Bact Synapomorphy candidates only
    syn_bact = df[
        (df["Region"] == "CORE") &
        (df["Interpretation"] == "SUPPORT_BACT") &
        (df["Conservative"] == False)
    ]
    syn_bact.to_csv(OUTPUT_SYN_B, index=False)

    # Compile Summary Report
    write_summary(df, OUTPUT_SUMM)