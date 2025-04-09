from Bio import Entrez, SeqIO
import time
from rich.progress import Progress, TextColumn, BarColumn, TimeRemainingColumn, SpinnerColumn, MofNCompleteColumn
import diskcache as dc
from multiprocessing import Pool, Lock, Manager
import os
import logging
from io import StringIO
import shutil


class EntrezClient:
    def __init__(self, db: str = "nucleotide"):
        self.db = db
        self.WORK_DIR = 'queue_work'
        self.PENDING_DIR = 'queue_pending'
        self.RESULTS_DIR = 'results'
        self.BATCH_SIZE = 100
        self.NUM_WORKERS = 2

        self.work = dc.Deque(directory=self.WORK_DIR)
        self.pending = dc.Deque(directory=self.PENDING_DIR)

    def search(self, term: str, retmax: int = 1000000000):
        handle = Entrez.esearch(db=self.db, term=term, retmax=retmax)
        results = Entrez.read(handle)
        handle.close()
        ids = results["IdList"]
        return ids

    def recover_pending(self):
        logging.info(f"Recovering {len(self.pending)} pending items")
        while self.pending:
            self.work.append(self.pending.popleft())

    def initialize(self, term: str, retmax: int = 1000000000):
        self.work.clear()
        self.pending.clear()
        if os.path.exists(self.RESULTS_DIR):
            shutil.rmtree(self.RESULTS_DIR)
        os.mkdir(self.RESULTS_DIR)
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

            logging.info(f"Worker {num} processing {len(ids)} items")

            while True:
                try:
                    handle = Entrez.efetch(db=db, id=",".join(ids), rettype=rettype, retmode=retmode)
                    filename = f"worker_{num}_{int(time.time())}.fasta"
                    with open(os.path.join(results_dir, filename), "a") as out_handle:
                        records = SeqIO.parse(handle, "fasta")
                        count = SeqIO.write(records, out_handle, "fasta")
                    handle.close()
                    logging.info(f"Worker {num} wrote {count} items to results")

                    for id in ids:
                        pending.remove(id)

                    break
                except Exception as e:
                    if os.path.exists(os.path.join(results_dir, filename)):
                        os.remove(os.path.join(results_dir, filename))
                    logging.info(f"Worker {num} encountered exception, waiting for retry: {e}")
                    time.sleep(10)

    def download(self, path: str, rettype: str = "fasta", retmode: str = "text"):
        self.recover_pending()
        logging.info(f"Found {len(self.work)} items in work queue")

        time.sleep(5)

        with Pool(processes=self.NUM_WORKERS) as pool:
            args = [(i, self.WORK_DIR, self.PENDING_DIR, self.RESULTS_DIR, self.BATCH_SIZE, self.db, rettype, retmode) for i in range(self.NUM_WORKERS)]
            pool.starmap(self.worker, args)
