"""
Introduce in-frame stop codons into HIV consensus FASTA for CFEIntact testing.
Introducing in-frame stop codons (TAA):
  gag: position 892 | original codon: GTA -> TAA (stop)
  pol: position 2502 | original codon: ATA -> TAA (stop)
Stop codon used: TAA
"""

import sys

def read_fasta(filepath):
    with open(filepath) as f:
        lines = f.readlines()
    header = lines[0].strip()
    genome = list("".join(l.strip() for l in lines[1:]))
    return header, genome

def find_inframe_position(gene_start, target_pos):
    """Find nearest in-frame codon start at or after target_pos (1-based GTF coords)."""
    offset = (target_pos - gene_start) % 3
    if offset == 0:
        return target_pos
    return target_pos + (3 - offset)

def introduce_stop(genome, pos_1based, label):
    """Replace codon at pos_1based (1-based) with TAA stop codon."""
    idx = pos_1based - 1  # convert to 0-based
    original_codon = "".join(genome[idx:idx+3])
    genome[idx]   = "T"
    genome[idx+1] = "A"
    genome[idx+2] = "A"
    print(f"  {label}: position {pos_1based} | original codon: {original_codon} -> TAA (stop)")
    return genome

# Gene start positions from GTF (1-based)
genes = {
    "gag":    790,
    "pol":    2358
}

# Target roughly 100bp into each gene for a clean test
targets = {
    "gag":    890,
    "pol":    2500
}

if len(sys.argv) < 2:
    print("Usage: python introduce_stop_codons.py <consensus.fasta>")
    sys.exit(1)

input_fasta = sys.argv[1]
output_fasta = input_fasta.replace(".fasta", "_stoptest.fasta")

header, genome = read_fasta(input_fasta)
print(f"Loaded: {header} ({len(genome)} bp)\n")
print("Introducing in-frame stop codons (TAA):")

for gene, gene_start in genes.items():
    target = targets[gene]
    inframe_pos = find_inframe_position(gene_start, target)
    genome = introduce_stop(genome, inframe_pos, gene)

# Write output
with open(output_fasta, "w") as f:
    f.write(header + "_stoptest\n")
    # Write sequence in 60bp lines
    seq = "".join(genome)
    for i in range(0, len(seq), 60):
        f.write(seq[i:i+60] + "\n")

print(f"\nOutput written to: {output_fasta}")
print("Run CFEIntact on this file to verify stop codon detection.")
