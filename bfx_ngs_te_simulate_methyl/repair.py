import os
from smart_open import open as sopen
import tempfile
import shutil


def stream_read_fastq(fastq_file):
    """Read a fastq file and return an iterable of the sequence ID, the full header line,
    the sequence, and the quality scores. Note that the sequence ID is the header up until the first space,
    while the header is the whole header line.
    """
    # Illumina headers example
    # @NS500726:623:HVFVNBGXF:1:11101:22894:1069 1:N:0:ATTACGCG+NGGCTATA
    # @NS500726:623:HVFVNBGXF:1:11101:22894:1069 2:N:0:ATTACGCG+NGGCTATA
    # BGI headers example
    # @V300052679L1C001R0030000035/1
    # @V300052679L1C001R0030000035/2
    # SRA headers example
    # @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=72
    # @SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=72
    # open file handle stream NOTE: context manager syntax does not play well with generator functions
    fin = sopen(fastq_file, encoding='utf-8')
    while True:
        # read line by line
        header = fin.readline()
        if not header:
            fin.close()
            break
        header = header.strip().lstrip('@')
        seqid = header.split(' ')[0]
        seq = fin.readline().strip()
        _ = fin.readline()
        qualscores = fin.readline().strip()
        yield seqid, header, seq, qualscores


def repair_fastqs_files(read1_fastq, read2_fastq):
    print(f'Repairing fastqs read 1: {read1_fastq} read 2: {read2_fastq}...')
    # move to tempdir
    r1_temp = tempfile.NamedTemporaryFile().name
    shutil.copy(read1_fastq, r1_temp)
    r2_temp = tempfile.NamedTemporaryFile().name
    shutil.copy(read2_fastq, r2_temp)
    # define the output files for paired and unpaired reads
    #read1_paired_fastq = os.path.splitext(read1_fastq.rstrip('.gz'))[0] + '_paired.fastq.gz'
    #read2_paired_fastq = os.path.splitext(read2_fastq.rstrip('.gz'))[0] + '_paired.fastq.gz'
    read1_paired_fastq = read1_fastq
    read2_paired_fastq = read2_fastq
    #read1_unpaired_fastq = os.path.splitext(read1_fastq.rstrip('.gz'))[0] + '_unpaired.fastq.gz'
    #read2_unpaired_fastq = os.path.splitext(read2_fastq.rstrip('.gz'))[0] + '_unpaired.fastq.gz'
    # open files for writing out
    r1_out = sopen(read1_paired_fastq, 'w')
    r2_out = sopen(read2_paired_fastq, 'w')
    #r1_up_out = sopen(read1_unpaired_fastq, 'w')
    #r2_up_out = sopen(read2_unpaired_fastq, 'w')
    # read in the R1 fastq into mem
    r1_seqs = dict()
    for (seqid, header, seq, qual) in stream_read_fastq(r1_temp):
        # standardize seqid by stripping read info
        seqid = seqid.replace('.1', '').replace('/1', '')
        r1_seqs[seqid] = [header, seq, qual]
    # init counters
    paired_count = 0
    unpaired_r1_count = 0
    unpaired_r2_count = 0
    # stream process R1, write out paired reads and unpaired R2 reads
    seen = set()
    for (seqid, header, seq, qual) in stream_read_fastq(r2_temp):
        # standardize seqid by stripping read info
        seqid = seqid.replace('.2', '').replace('/2', '')
        seen.add(seqid)
        if seqid in r1_seqs:
            r1_out.write("@" + r1_seqs[seqid][0] + "\n" + r1_seqs[seqid][1] + "\n+\n" + r1_seqs[seqid][2] + "\n")
            r2_out.write("@" + header + "\n" + seq + "\n+\n" + qual + "\n")
            paired_count += 1
        else:
            #r2_up_out.write("@" + header + "\n" + seq + "\n+\n" + qual + "\n")
            unpaired_r2_count += 1
    # write out unpaired R1 reads
    for seqid in r1_seqs:
        if seqid not in seen:
            #r1_up_out.write("@" + r1_seqs[seqid][0] + "\n" + r1_seqs[seqid][1] + "\n+\n" + r1_seqs[seqid][2] + "\n")
            unpaired_r1_count += 1
    # close files
    r1_out.close()
    r2_out.close()
    #r1_up_out.close()
    #r2_up_out.close()
    # return paths
    print(f'Repaired fastqs read 1: {read1_paired_fastq} read 2: {read2_paired_fastq} total paired count: {paired_count}')
    #print(f'Unpaired fastq read 1: {read1_unpaired_fastq} total unpaired R1 count: {unpaired_r1_count}')
    #print(f'Unpaired fastq read 2: {read2_unpaired_fastq} total unpaired R2 count: {unpaired_r2_count}')
    return read1_paired_fastq, read2_paired_fastq