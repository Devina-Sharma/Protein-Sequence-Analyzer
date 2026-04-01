from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def analyze_protein(sequence):
    analysis = ProteinAnalysis(sequence)

    return {
        "length": len(sequence),
        "molecular_weight": analysis.molecular_weight(),
        "aromaticity": analysis.aromaticity(),
        "instability_index": analysis.instability_index(),
        "isoelectric_point": analysis.isoelectric_point(),
        "gravy": analysis.gravy()  
    }


def process_fasta(file):
    results = []

    for record in SeqIO.parse(file, "fasta"):
        seq = str(record.seq)
        data = analyze_protein(seq)

        results.append({
            "id": record.id,
            **data
        })

    return results


def save_results(results, filename):
    with open(filename, "w") as f:
        f.write("ID\tLength\tMolWt\tAromaticity\tInstability\tpI\tGRAVY\n")
        for r in results:
            f.write(
                f"{r['id']}\t{r['length']}\t"
                f"{r['molecular_weight']:.2f}\t"
                f"{r['aromaticity']:.3f}\t"
                f"{r['instability_index']:.2f}\t"
                f"{r['isoelectric_point']:.2f}\t"
                f"{r['gravy']:.3f}\n"
            )


if __name__ == "__main__":
    fasta_file = "protein.fasta"

    results = process_fasta(fasta_file)

    print("PROTEIN ANALYSIS RESULTS:\n")
    for r in results:
        print(f"ID: {r['id']}")
        print(f" Length: {r['length']}")
        print(f" Molecular Weight: {r['molecular_weight']:.2f}")
        print(f" Aromaticity: {r['aromaticity']:.3f}")
        print(f" Instability Index: {r['instability_index']:.2f}")
        print(f" Isoelectric Point (pI): {r['isoelectric_point']:.2f}")
        print(f" GRAVY: {r['gravy']:.3f}")
        print("-" * 40)

    save_results(results, "protein_results.txt")

    print("\nResults saved to protein_results.txt")
