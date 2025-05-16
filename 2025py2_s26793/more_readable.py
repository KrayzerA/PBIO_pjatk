from Bio import Entrez, SeqIO
import pandas as pd, matplotlib.pyplot as plt
from io import StringIO


class NCBIRetriever:
    def __init__(self, email, key):
        Entrez.email, Entrez.api_key = email, key

    def search(self, taxid):
        r = Entrez.read(Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism]", usehistory="y"))
        self.webenv, self.qkey, self.count = r["WebEnv"], r["QueryKey"], int(r["Count"])

    def fetch(self, start, maxr):
        h = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                          retstart=start, retmax=maxr, webenv=self.webenv, query_key=self.qkey)
        return list(SeqIO.parse(StringIO(h.read()), "gb"))



def main():
    email = input("Enter your email: ")
    apikey = input("Enter your NCBI API key: ")
    taxid = input("Enter taxonomic ID (taxid): ")
    minlen = int(input("Minimum sequence length: "))
    maxlen = int(input("Maximum sequence length: "))
    limit = int(input("Max records to fetch (e.g. 500): "))

    ncbi = NCBIRetriever(email, apikey)
    ncbi.search(taxid)
    print(f"Found {ncbi.count} records.")

    seqs, step = [], 500
    for i in range(0, ncbi.count, step):
        if len(seqs) >= limit:
            break
        print(f"Fetching {i} - {i + step}...")
        for rec in ncbi.fetch(i, step):
            l = len(rec.seq)
            if minlen <= l <= maxlen:
                seqs.append((rec.id, l, rec.description))
            if l >= limit:
                break

    if not seqs:
        print("No records in specified length range.");
        return

    df = pd.DataFrame(seqs, columns=["Accession", "Length", "Description"])
    df.sort_values("Length", ascending=False, inplace=True)
    df.to_csv(f"taxid_{taxid}_report.csv", index=False)
    print(f"Saved CSV as taxid_{taxid}_report.csv")
    
    plt.figure(figsize=(10, 5))
    plt.plot(df["Accession"], df["Length"], marker='o')
    plt.xticks(rotation=90, fontsize=6)
    plt.xlabel("Accession")
    plt.ylabel("Sequence Length")
    plt.tight_layout()
    plt.savefig(f"taxid_{taxid}_plot.png", dpi=300)
    print(f"Saved plot as taxid_{taxid}_plot.png")


if __name__ == "__main__":
    main()
