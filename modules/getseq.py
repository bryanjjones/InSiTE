#!/usr/bin/env python
import colorama
import Bio
import Bio.Entrez
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord


# fetch indicated sequence from local 2bit genome. chrom given as number (or X/Y/NT), sense given as '+'/'-'.
def fetch_seq_from_twobit(genome, chrom, sense, loc, lwindow, rwindow, samwindow):
    if sense == '-':
        seq_start = int(loc + samwindow + lwindow)  # range of sequence to return
        seq_stop = int(loc + samwindow - rwindow)
        sequence = genome[f'chr{chrom}'][seq_stop:seq_start]
        record = Bio.SeqRecord.SeqRecord(Bio.Seq.reverse_complement(Bio.Seq.Seq(sequence)),
                                         id=f'chr{chrom}{sense}:{seq_start}-{seq_stop}', description='')
        return record
    elif sense == '+':
        seq_start = int(loc - lwindow)  # range of sequence to return
        seq_stop = int(loc + rwindow)
        sequence = genome[f'chr{chrom}'][seq_start:seq_stop]
        record = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(sequence),
                                         id=f'chr{chrom}{sense}:{seq_start}-{seq_stop}', description='')
        return record


def get_seq(read_item, lwindow, rwindow, samwindow, source="LOCAL", genome=None):
    if source == 'LOCAL':
        record = fetch_seq_from_twobit(genome, read_item.chrom, read_item.sense, read_item.loc, lwindow, rwindow,
                                       samwindow)
        sequence = record.seq
    elif source == 'REMOTE':
        handle = Bio.Entrez.efetch(db="nucleotide",  # fetches the specified sequence as a fasta
                                   id=chromIDS[read_item.chrom],
                                   # id for indicated chromosome from refrence genome assembly
                                   rettype="fasta",  # format
                                   strand=read_item.sensenum,  # which strand
                                   # the behaviour of seq_start and seq_stop might vary depending on sense vs antisense
                                   seq_start=str(read_item.loc - lwindow),  # range of sequence to return
                                   seq_stop=str(read_item.loc + rwindow - 1))
        record = Bio.SeqIO.read(handle, "fasta")  # reads the returned fasta into a sequence object
        handle.close()  # close efetch
        sequence = record.seq  # adds the sequence property as a string
    return record, sequence


def get_seqs(readlist, lwindow, rwindow, samwindow, genome=None, source="LOCAL"):
    line_count = 0
    readlistwseqs = []
    recordlist = []
    if source == 'REMOTE':
        Bio.Entrez.email = "bryan@bmogen.com"  # Always tell NCBI who you are
        Bio.Entrez.api_key = "679368bc3fef47bb440fb743c889befe4e09"  # bryan's personal NCBI key (allows faster
        # requests of their server)
    for read in readlist:
        record, sequence = get_seq(read, lwindow, rwindow, samwindow, genome=genome, source=source)
        read.sequence = sequence
        read.record = record
        readlistwseqs.append(read)
        recordlist.append(record)
        print(f'Sequence ' + colorama.Fore.GREEN + f'{line_count:,} ' + colorama.Style.RESET_ALL +
              'mapped to chrom ' + colorama.Fore.GREEN + f'{read.chrom}' + colorama.Style.RESET_ALL + ' with ' +
              colorama.Fore.GREEN + f'{read.sense}' + colorama.Style.RESET_ALL + ' sense at  ' + colorama.Fore.GREEN +
              f'{read.loc} ' + colorama.Style.RESET_ALL + 'with seq ' + colorama.Fore.BLUE +
              f'{read.sequence[:5]}...{read.sequence[-5:]}' + colorama.Style.RESET_ALL,
              end="\r")
        line_count += 1
    return recordlist, readlistwseqs
