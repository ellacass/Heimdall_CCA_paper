import re
import pandas as pd
from collections import defaultdict

# Define motif patterns (regular expressions)
motifs = {
    "PAP Upstream": re.compile(r"[LIV][LIV]G[RK][RK]F.[LIV][AILMVF][HQL][LIV]"),
    "Motif A": re.compile(r"GG..R."),
    "Motif B": re.compile(r"RRD"),
    "Motif C": re.compile(r"D...G"),
    "Motif D": re.compile(r"D..R..R"),
    "CCA Motif E": re.compile(r"ER...E"),
    "A motif E": re.compile(r"[ST]R.{3}E.{3}[AVLFIMW]{2}"),
    "PAP Motif E": re.compile(r"ARL.[ED]E..K.L")
}

def find_motifs(sequence):
    """Check motif positions in a given sequence and store them."""
    motif_positions = {}
    
    for motif, pattern in motifs.items():
        match = pattern.search(sequence)
        if match:
            # Store the start position of the first match
            motif_positions[motif] = match.start()
        else:
            # Store None if no match is found
            motif_positions[motif] = None
    return motif_positions

def process_fasta(file_path, output_csv):
    """Reads a FASTA file, checks motif order, and saves results to a CSV file."""
    sequences = defaultdict(str)
    header = None

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
            elif header:
                sequences[header] += line
    
    # --- NEW STEP: Create a DataFrame and remove duplicates ---
    # Create a DataFrame from the sequences dictionary
    df_raw = pd.DataFrame(sequences.items(), columns=['Protein', 'Sequence'])
    
    # Remove rows where the 'Sequence' column is a duplicate
    df_unique = df_raw.drop_duplicates(subset=['Sequence']).copy()
    print(f"Removed {len(df_raw) - len(df_unique)} duplicate sequences.")
    
    results = []

    # Iterate over the unique sequences DataFrame
    for _, row in df_unique.iterrows():
        header = row['Protein']
        sequence = row['Sequence']
        
        motif_positions = find_motifs(sequence)

        core_motifs = ["Motif A", "Motif B", "Motif C", "Motif D"]
        all_core_present = all(motif_positions[m] is not None for m in core_motifs)

        has_cca_e = motif_positions["CCA Motif E"] is not None
        has_a_e = motif_positions["A motif E"] is not None
        has_pap_e = motif_positions["PAP Motif E"] is not None
        has_pap_upstream = motif_positions["PAP Upstream"] is not None

        if all_core_present:
            if has_cca_e and not has_pap_e and not has_a_e and not has_pap_upstream:
                enzyme_type = "CCA"
            elif (has_pap_e or has_pap_upstream) and not has_a_e:
                enzyme_type = "PAP"
            elif has_a_e and not has_cca_e and not has_pap_e and not has_pap_upstream:
                enzyme_type = "A-adding"
            else:
                enzyme_type = "Unknown (Motif E type unclear)"
        else:
            enzyme_type = "Incomplete (Missing core motifs A-D)"

        missing_motifs = [motif for motif, pos in motif_positions.items() if pos is None]
        missing_str = ", ".join(missing_motifs) if missing_motifs else "None"

        detected_positions = [pos for pos in motif_positions.values() if pos is not None]
        first_motif_position = min(detected_positions) if detected_positions else None
        insufficient_upstream = "Yes" if first_motif_position is not None and first_motif_position < 20 else "No"
        final_classification = enzyme_type if "Unknown" not in enzyme_type and "Incomplete" not in enzyme_type else "Unclassified"

        results.append([header, enzyme_type, missing_str, insufficient_upstream, final_classification])

    df_final = pd.DataFrame(results, columns=["Protein", "Enzyme Type", "Missing Motifs", "Insufficient Upstream Space", "Final Classification"])
    df_final.to_csv(output_csv, index=False)
    print(f"Results saved to {output_csv}")

fasta_file = "/Users/pat/Documents/PhD/Heimdalls/origin_stories/bacteria/CCA2_ecoli_REF/subtilis_ref/missing_genomes/filter_blastoutputs/prodigal_output.fasta"
output_csv = "/Users/pat/Documents/PhD/Heimdalls/origin_stories/bacteria/CCA2_ecoli_REF/subtilis_ref/missing_genomes/filter_blastoutputs/PAP_vs_CCA_analysis.csv"
process_fasta(fasta_file, output_csv)