from Bio import Entrez, SeqIO
import time
from rich.progress import Progress, TextColumn, BarColumn, TimeRemainingColumn, SpinnerColumn, MofNCompleteColumn
import diskcache as dc
from multiprocessing import Pool
import os
import logging
import shutil


class EntrezClient:
    def __init__(self, db: str = "nucleotide", batch_size: int = 1000, num_workers: int = 4, rettype: str = "fasta", retmode: str = "text"):
        self.db = db
        self.work_dir = "queue_work"
        self.pending_dir = "queue_pending"
        self.results_dir = "temp_results"
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.rettype = rettype
        self.retmode = retmode

        self.work = dc.Deque(directory=self.work_dir)
        self.pending = dc.Deque(directory=self.pending_dir)

    def search(self, term: str, retmax: int = 1000000000):
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            "[progress.percentage]{task.percentage:>3.0f}%",
            TimeRemainingColumn(),
        ) as progress:
            task = progress.add_task("[red]Searching", total=1)

            handle = Entrez.esearch(db=self.db, term=term, retmax=retmax)
            results = Entrez.read(handle)
            handle.close()
            ids = results["IdList"]
            progress.update(task, completed=1)
            return ids

    def recover_pending(self):
        logging.info(f"Recovering {len(self.pending)} pending items")
        while self.pending:
            self.work.append(self.pending.popleft())

    def initialize(self, term: str, retmax: int = 1000000000):
        self.work.clear()
        self.pending.clear()
        if os.path.exists(self.results_dir):
            shutil.rmtree(self.results_dir)
        os.mkdir(self.results_dir)
        ids = self.search(term, retmax=retmax)
        self.work.extend(ids)

    @staticmethod
    def worker(num: int, work_dir: str, pending_dir: str, results_dir: str, batch_size: int, db: str, rettype: str, retmode: str):
        work = dc.Deque(directory=work_dir)
        pending = dc.Deque(directory=pending_dir)

        while True:
            ids = []
            for _ in range(batch_size):
                try:
                    new_id = work.popleft()
                    ids.append(new_id)
                    pending.append(new_id)
                except IndexError:
                    break
            if not ids:
                break

            logging.debug(f"Worker {num} processing {len(ids)} items")

            while True:
                try:
                    handle = Entrez.efetch(db=db, id=",".join(ids), rettype=rettype, retmode=retmode)
                    filename = f"worker_{num}_{int(time.time())}.fasta"
                    with open(os.path.join(results_dir, filename), "a") as out_handle:
                        records = SeqIO.parse(handle, "fasta")
                        count = SeqIO.write(records, out_handle, "fasta")
                    handle.close()
                    logging.debug(f"Worker {num} wrote {count} items to results")

                    for id in ids:
                        pending.remove(id)

                    break
                except Exception as e:
                    if os.path.exists(os.path.join(results_dir, filename)):
                        os.remove(os.path.join(results_dir, filename))
                    logging.debug(f"Worker {num} encountered exception, waiting for retry: {e}")
                    time.sleep(10)

    def download(self, path: str):
        self.recover_pending()
        logging.info(f"Found {len(self.work)} items in work queue")

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            "[progress.percentage]{task.percentage:>3.0f}%",
            TimeRemainingColumn(),
        ) as progress:
            total_items = len(self.work) + len(self.pending)
            task = progress.add_task("[red]Downloading", total=total_items)

            with Pool(processes=self.num_workers) as pool:
                results = []
                for i in range(self.num_workers):
                    args = (i, self.work_dir, self.pending_dir, self.results_dir, self.batch_size, self.db, self.rettype, self.retmode)
                    results.append(pool.apply_async(self.worker, args))

                all_done = False
                while not all_done:
                    remaining = len(self.work) + len(self.pending)
                    completed = total_items - remaining
                    progress.update(task, completed=completed)
                    all_done = all(result.ready() for result in results)
                    time.sleep(10)
                
                progress.update(task, completed=total_items)

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            "[progress.percentage]{task.percentage:>3.0f}%",
            TimeRemainingColumn(),
        ) as progress:
            all_files = os.listdir(self.results_dir)
            task = progress.add_task("[red]Merging", total=len(all_files))

            with open(path, "w") as out_handle:
                for i, filename in enumerate(all_files):
                    with open(os.path.join(self.results_dir, filename), "r") as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
                        progress.update(task, completed=i+1)
