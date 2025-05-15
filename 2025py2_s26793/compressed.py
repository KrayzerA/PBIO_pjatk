from Bio import Entrez,SeqIO;import pandas as pd,matplotlib.pyplot as plt;from io import StringIO
class NCBI:
 def __init__(s,e,k):Entrez.email,Entrez.api_key=e,k
 def s(s,t):r=Entrez.read(Entrez.esearch(db="nucleotide",term=f"txid{t}[Organism]",usehistory="y"));s.w,s.q,s.c=r["WebEnv"],r["QueryKey"],int(r["Count"])
 def f(s,i,m):h=Entrez.efetch(db="nucleotide",rettype="gb",retmode="text",retstart=i,retmax=m,webenv=s.w,query_key=s.q);return list(SeqIO.parse(StringIO(h.read()),"gb"))
def main():
 # e,k='s26793@pjwstk.edu.pl','d7af0ba424e4e31cd83a054c5bd530d13708'
 e=input("email:");k=input("apikey:");t=input("taxid: ");a=int(input("min len: "));b=int(input("max len: "));l=int(input("limit: "))
 n=NCBI(e,k);n.s(t);print(f"Found {n.c} records, fetching up to {l}...")
 x=[];s=500
 for i in range(0,min(n.c,l),s):
  print(f"{i}-{i+min(s,l-i)}");x+=[(r.id,len(r.seq),r.description)for r in n.f(i,min(s,l-i))if a<=len(r.seq)<=b]
 if not x:print("None in range.");return
 d=pd.DataFrame(x,columns=["Accession","Length","Description"]).sort_values("Length",ascending=0)
 d.to_csv(f"taxid_{t}_report.csv",index=0);print("CSV saved.")
 plt.figure(figsize=(10,5));plt.plot(d["Accession"],d["Length"],'o');plt.xticks(rotation=90,fontsize=6)
 plt.xlabel("Accession");plt.ylabel("Length");plt.tight_layout();plt.savefig(f"taxid_{t}_plot.png",dpi=300);print("Plot saved.")
main()
