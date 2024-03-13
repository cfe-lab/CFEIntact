from dataclasses import dataclass

@dataclass
class BlastRow:
    qseqid: str
    sseqid: str
    sgi: str
    qlen: int
    slen: int
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    pident: float
    nident: float
    sstrand: str

    @staticmethod
    def init(row):
        it = iter(row)
        return BlastRow(
                qseqid=row['qseqid'],
                sseqid=row['sseqid'],
                sgi=row['sgi'],
                qlen=int(row['qlen']),
                slen=int(row['slen']),
                length=int(row['length']),
                qstart=int(row['qstart']),
                qend=int(row['qend']),
                sstart=int(row['sstart']),
                send=int(row['send']),
                evalue=float(row['evalue']),
                bitscore=float(row['bitscore']),
                pident=float(row['pident']),
                nident=float(row['nident']),
                sstrand=row['sstrand'],
        )
