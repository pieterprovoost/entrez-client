from Bio import Entrez, SeqIO
import time
from rich.progress import Progress, TextColumn, BarColumn, TimeRemainingColumn, SpinnerColumn, MofNCompleteColumn


class EntrezClient:
    def __init__(self, email: str, api_key: str, db: str = "nucleotide"):
        self.email = email
        self.api_key = api_key
        self.db = db
        self.batch_size = 1000
        Entrez.email = email
        Entrez.api_key = api_key

    def search(self, term: str, retmax: int = 1000000000):
        handle = Entrez.esearch(db=self.db, term=term, retmax=retmax)
        results = Entrez.read(handle)
        handle.close()
        ids = results["IdList"]
        return ids

    def download_ids(self, term: str, path: str, retmax: int = 1000000000):
        ids = self.search(term, retmax=retmax)
        with open(path, "w") as f:
            for id in ids:
                f.write(f"{id}\n")

    def download(self, term: str, path: str, retmax: int = 1000000000, rettype: str = "fasta", retmode: str = "text"):
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            "[progress.percentage]{task.percentage:>3.0f}%",
            TimeRemainingColumn(),
        ) as progress:
            search_task = progress.add_task("[red]Searching", total=1)
            ids = self.search(term, retmax=retmax)
            progress.update(search_task, completed=1)

            open(path, "w").close()
            
            offsets = range(0, len(ids), self.batch_size)
            task = progress.add_task("[green]Downloading", total=len(ids))

            for start in offsets:
                end = min(start + self.batch_size, len(ids))
                batch_ids = ids[start:end]
                while True:
                    try:
                        handle = Entrez.efetch(db=self.db, id=",".join(batch_ids), rettype=rettype, retmode=retmode)
                        with open(path, "a") as out_handle:
                            records = SeqIO.parse(handle, "fasta")
                            count = SeqIO.write(records, out_handle, "fasta")
                        progress.update(task, advance=count)
                        handle.close()
                        break
                    except Exception as e:
                        time.sleep(10)
