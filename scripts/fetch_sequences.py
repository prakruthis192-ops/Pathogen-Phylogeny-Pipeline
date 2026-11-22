import os
import time
from Bio import Entrez, SeqIO

Entrez.email = "prakruthis192@gmail.com"

SEARCH_TERM = (
    '(Influenza A virus OR Alphainfluenzavirus influenzae) AND H3N2 AND '
    'hemagglutinin AND "complete cds" AND 1650:1750[SLEN]'
)
MAX_RECORDS = 25
OUT_PATH = "data/raw/h3n2_ha_raw.fasta"

def fetch_ncbi_sequences():
    print(f"Searching NCBI for: {SEARCH_TERM}")
    handle = Entrez.esearch(db="nuccore", term=SEARCH_TERM, retmax=MAX_RECORDS, sort="relevance")
    record = Entrez.read(handle)
    handle.close()
    
    id_list = record["IdList"]
    if not id_list:
        raise RuntimeError("No sequences found! Check your internet connection or search query.")
        
    print(f"Found {len(id_list)} matching records. Downloading FASTA...")
    
    time.sleep(0.5)
    fetch_handle = Entrez.efetch(db="nuccore", id=",".join(id_list), rettype="fasta", retmode="text")
    raw_data = fetch_handle.read()
    fetch_handle.close()
    
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w") as f:
        f.write(raw_data)
        
    seqs = list(SeqIO.parse(OUT_PATH, "fasta"))
    print(f"Saved {len(seqs)} sequences to {OUT_PATH}")

if __name__ == "__main__":
    fetch_ncbi_sequences()

