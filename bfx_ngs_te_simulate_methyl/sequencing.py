#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## MethylFASTQ generates artificial bisulfite data in FASTQ format.
## Copyright (C) 2018, Nicola Licheri (nicola.licheri@gmail.com)
#
## This file is part of MethylFASTQ.
#
## MethylFASTQ is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
#
## MethylFASTQ is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
#
## You should have received a copy of the GNU General Public License
## along with MethylFASTQ.  If not, see <http://www.gnu.org/licenses/>.

from Bio import SeqIO
from io import StringIO
import random, csv
import os, pickle
import multiprocessing as mp
import threading
from timeit import default_timer as timer
import enum
import sys
import math
from scipy.stats import norm, truncnorm

import bfx_ngs_te_simulate_methyl.dna as dna


def format_time(time):
    strtime = "{:03.5f} m".format(time/60) if time >= 60 else "{:03.5f} s".format(time)
    return strtime

class CytosineContext(enum.Enum):
    CG = "CG"
    CHG = "CHG"
    CHH = "CHH"

    def __str__(self):
        return self.value

class Cytosine(object):
    def __init__(self, strand, context=CytosineContext.CHH):
        self.__context = context
        self.__strand = strand
        self.__nmeth = 0
        self.__ncov = 0

    @property
    def context(self):
        return self.__context

    @property
    def strand(self):
        return self.__strand

    @property
    def nmeth(self):
        return self.__nmeth

    @property
    def ncov(self):
        return self.__ncov

    def methylate(self):
        self.__nmeth += 1

    def covered(self):
        self.__ncov += 1


class Stats(object):
    def __init__(self):
        self.__num_reads = 0
        self.__covered_cs = 0
        self.__bp_sequenced = 0

    @property
    def nreads(self):
        return self.__num_reads

    @property
    def ncytosines(self):
        return self.__covered_cs

    @property
    def nbases(self):
        return self.__bp_sequenced

    def increment_reads(self, nr):
        self.__num_reads += nr

    def increment_cytosines(self, nc):
        self.__covered_cs += nc

    def increment_bps(self, nbp):
        self.__bp_sequenced += nbp

    def update(self, stats):
        self.increment_reads(stats.nreads)
        self.increment_cytosines(stats.ncytosines)
        self.increment_bps(stats.nbases)


class ChromosomeSequencer(object):
    """ ChromosomeSequencer legge un cromosoma da un file FASTA e
    scrive su un file FASTQ le read prodotte dal sequenziamento.
    È possibile scegliere il tipo di library (single-end o paired-end),
    se la library è direzionale o non direzionale, la lunghezza dei frammenti
    e delle read """

    def __init__(self, chromosome, params, target_regions=list()):
        """Parserizza il cromosoma individuando gli intervalli da sequenziare e scartando
        i nucleotidi indefiniti (basi N)"""

        chromosome_sequence = str(chromosome.seq).lower()

        self.__chromoId = chromosome.id
        self.__fragment_size_mean = params.fragment_size_mean
        self.__fragment_size_sd = params.fragment_size_sd
        self.__targeted_mode = True if len(target_regions) > 0 else False
        self.__fragments = [(chromosome_sequence[begin:end], begin, end) for (begin, end) in target_regions]
        self.__stats = Stats()

        print("Parsing {} sequence... ".format(chromosome.id), end="", flush=True)
        chr_size = len(chromosome.seq)

        # get salient fragments
        start, i = timer(), 0

        if len(target_regions) == 0:
            # WGBS option
            while i < chr_size:
                while i < chr_size and chromosome_sequence[i] == 'n':
                    i += 1
                begin = i
                while i < chr_size and chromosome_sequence[i] != 'n':
                    i += 1

                t = (chromosome_sequence[begin:i], begin, i)
                self.__fragments.append(t)
            else:
                last = self.__fragments.pop()
                if last[1] < last[2]:
                    self.__fragments.append(last)
        # else:
        #     # Targeted option
        #     expanded_targeted_regions = []
        #     print(' ')
        #     for seq, begin, end in self.__fragments:
        #         print(f'Original start: {begin} end: {end}')
        #         if not len(seq) > self.__fragment_size:
        #             delta = self.__fragment_size - len(seq) + 2
        #             if not (delta % 2) == 0:
        #                 delta = delta + 1
        #             begin = int(begin - (delta/2))
        #             if begin < 0:
        #                 begin = 0
        #             end = int(end + (delta/2))
        #             if end > chr_size:
        #                 end = chr_size
        #             print(f'Expanded start: {begin} end: {end}')
        #             seq = chromosome_sequence[begin:end]
        #         t = (seq, begin, end)
        #         expanded_targeted_regions.append(t)
        #     self.__fragments = expanded_targeted_regions

        tot_time = format_time(timer() - start)

        #sommo dimensioni intervalli
        self.__stats.increment_bps(sum([e-b for _, b, e in self.__fragments]))
        print("{} fragments found: {} bp. Elapsed time {}".format(len(self.__fragments), self.__stats.nbases, tot_time), flush=True)

    def load_balancing(self, num_workers):
        totsize = sum([e-b for _, b, e in self.__fragments])
        average = int(totsize / num_workers)

        fragments = list()
        #sort per visualizzare i frammenti in modo non caotico
        self.__fragments.sort(key=lambda x: x[2]-x[1], reverse=True)

        for i, (sequence, begin, end) in enumerate(self.__fragments, 1):
            size = end-begin
            print("{}. Fragment [{} - {}] of size {} bp".format(i, begin, end, size))

            if size > average: #suddivido il lavoro
                num_pieces = math.ceil(size / average)
                prev, curr = 0, average #offset di inizio e fine

                print(">> {} pieces: ".format(num_pieces), end="")

                for n in range(num_pieces):
                    subseq = (sequence[prev:curr], begin+prev, begin+curr)
                    fragments.append(subseq)

                    #print("[{} - {}] ".format(subseq[1], subseq[2]), end="")

#                    print("\tsubjob {} from {} to {}".format(n+1, begin+prev, begin+curr))
                    prev, curr = curr, curr + average

                    if curr > size:
                        curr = size
                print()
            else:
                fragments.append((sequence, begin, end))

        self.__fragments = sorted(fragments, key=lambda x: x[2]-x[1], reverse=True)

        #for seq, begin, end in self.__fragments:
        #    print(begin, end, end-begin, seq)


    def consumer(self, num_jobs, params, queue):
        """ Thread consumer che legge dalla coda e scrive i dati sull'apposito file"""

        filename = self.__get_output_filename(params)
        single_end = params.seq_mode == "single_end"

        ## apro file di output
        if single_end:
            fastq_file = open("{}.fastq".format(filename), "a")
        else:
            fastq_file1 = open("{}_R1.fastq".format(filename), "a")
            fastq_file2 = open("{}_R2.fastq".format(filename), "a")

        meth_file = open("{}.ch3".format(filename), "a")
        csv_meth = csv.writer(meth_file, delimiter="\t")

        #
        while num_jobs > 0:
            val = queue.get()
            tval = type(val)

            if tval is int:     #segnale di terminazione di un job
                num_jobs -= 1

            elif tval is tuple:
                datatype, data = val
                dsize = len(data)

                if datatype == "fastq_se":
                    for record in data:
                        SeqIO.write(record, fastq_file, "fastq")
                    else:
                        self.__stats.increment_reads(dsize)

                elif datatype == "fastq_pe":
                    for read1, read2 in data:
                        SeqIO.write(read1, fastq_file1, "fastq")
                        SeqIO.write(read2, fastq_file2, "fastq")
                    else:
                        self.__stats.increment_reads(dsize)

                elif datatype == "ch3":
                    for record in data:
                        csv_meth.writerow(record)
                    else:
                        self.__stats.increment_cytosines(dsize)
        else:
            #chiudo file
            meth_file.close()

            if single_end:
                fastq_file.close()
            else:
                fastq_file1.close()
                fastq_file2.close()

    def sequencing(self, params):
        self.load_balancing(params.num_processes)

#        self.__fragments = [{"seq": seq, "from": begin, "to": end, "params": params} \
#                            for seq, begin, end in self.__fragments]

        queue = mp.Queue() #comunicazione tra processi figli e padre
        num_input = len(self.__fragments)   #numero di input da processare
        num_process = params.num_processes  #numero massimo di processi da utilizzare
        single_end = params.seq_mode == "single_end" #variabile di merda, ma vbb

        input_data = [{
                "seq": seq,
                "offset": (begin, end),
                "params": params,
                "queue": queue,
                "process_id": index,
                "targeted_mode": self.__targeted_mode
            }
            for index, (seq, begin, end) in enumerate(self.__fragments)
        ]

        #apro i file di output
        output_filename = self.__get_output_filename(params)

        if single_end:
            fastq_file = open("{}.fastq".format(output_filename), "a")
        else:
            fastq_file1 = open("{}_S1_L001_R1_001.fastq".format(output_filename), "a")
            fastq_file2 = open("{}_S1_L001_R2_001.fastq".format(output_filename), "a")

        meth_file = open("{}.ch3".format(output_filename), "a")
        csv_meth = csv.writer(meth_file, delimiter="\t")

        #inizializzo i processi con i relativi input:
        #parametri: (frammento e offset) + coda + indice per il join
    #    processes = [mp.Process(target=self.create_reads, name="methylPIPPO-{}".format(index), args=(param, params.seq_mode, (index, queue))) \
    #                    for index, param in enumerate(self.__fragments)]
        processes = [mp.Process(target=self.create_reads, args=(data,)) for data in input_data]

        curr = 0        #indice del prossimo processo da startare
        curr_exec = 0   #numero processi attualmente in esecuzione

        #starto primi processi
        while curr_exec < min(num_input, num_process):
            processes[curr].start()
            curr += 1
            curr_exec += 1

        #finchè tutti i processi non sono stati eseguiti e non hanno terminato...
        while not (curr == num_input and curr_exec == 0):
            val = queue.get()
            tval = type(val)

            if tval is int: #segnale di terminazione di un processo
                processes[val].join()

                if curr < num_input:
                    processes[curr].start()
                    curr += 1
                else:
                    curr_exec -= 1
            elif tval is tuple:
                tag, data = val
                dsize = len(data)

                if tag == "fastq_se":
                    for record in data:
                        SeqIO.write(record, fastq_file, "fastq")
                    else:
                        self.__stats.increment_reads(dsize)

                elif tag == "fastq_pe":
                    for read1, read2 in data:
                        SeqIO.write(read1, fastq_file1, "fastq")
                        SeqIO.write(read2, fastq_file2, "fastq")
                    else:
                        self.__stats.increment_reads(dsize)

                elif tag == "ch3":
                    for record in data:
                        csv_meth.writerow(record)
                    else:
                        self.__stats.increment_cytosines(dsize)
        else: #fine while -> chiudo i file
            meth_file.close()

            if single_end:
                fastq_file.close()
            else:
                fastq_file1.close()
                fastq_file2.close()

        return self.__stats, "{}_S1_L001_R1_001.fastq".format(output_filename), "{}_S1_L001_R2_001.fastq".format(output_filename)

    def __sequencing(self, params):
        self.load_balancing(params.num_processes)

        queue = mp.Manager().Queue()

        ###### TERRRIBBBBILEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE!!!!!
        input_data = [{
                "seq": seq,
                "offset": (begin, end),
                "params": params,
                "queue": queue
            }
            for seq, begin, end in self.__fragments
        ]

        num_jobs = len(self.__fragments)

        consumer_thread = threading.Thread(target=self.consumer, args=(num_jobs, params, queue))
        consumer_thread.start()

        with mp.Pool(params.num_processes) as pool:
            pool.map(func=self.create_reads, iterable=input_data, chunksize=1)

        consumer_thread.join()

        return self.__stats



    def __get_output_filename(self, params):
        panelname = os.path.basename(params.regions).rstrip('.bed')
        frag_mean = f'fm{params.fragment_size_mean}'
        frag_sd = f'fsd{params.fragment_size_sd}'
        coverage = f'{params.coverage}X'
        read_len = f'rl{params.read_length}'
        cg = f'cg{params.p_cg}'
        chg = f'chg{params.p_chg}'
        chh = f'chh{params.p_chh}'
        snp = f'snp{params.snp_rate}'
        error = f'err{params.error_rate}'

        filename = f'{panelname}_{frag_mean}_{frag_sd}_{coverage}_{read_len}_{cg}_{chg}_{chh}_{snp}_{error}'

        return "{}/{}".format(params.output_path.rstrip('/'), filename.strip('/'))


    def create_reads(self, input_process):#, seq_mode, queue):
        sequence = input_process["seq"]
        offset_begin, offset_end = input_process["offset"]
        queue = input_process["queue"]
        params = input_process["params"]
        process_id = input_process["process_id"] #new
        targeted_mode = input_process["targeted_mode"]

        pid = os.getpid()

        print("<Process {}>: starting sequencing [{} - {}]: {} bp".format(pid, offset_begin, offset_end, offset_end-offset_begin), flush=True)
        start = timer()

        fs = FragmentSequencer(self.__chromoId, sequence, offset_begin, offset_end, params, targeted_mode=targeted_mode, queue=queue)
        fs.single_end_sequencing() if params.seq_mode == "single_end" else fs.paired_end_sequencing()

        elapsed = format_time(timer() - start)

        print("<Process {}>: sequencing [{} - {}] terminated in {}".format(pid, offset_begin, offset_end, elapsed), flush=True)

        queue.put(process_id)



class FragmentSequencer(object):
    """..."""

    def __init__(self, chr_id, sequence, begin_seq, end_seq, seqparams, targeted_mode, queue):
        #informazioni sulla sequenza
        self.__chromoId = chr_id
        self.__sequence = sequence
#        self.__genome_size = gsize
        self.targeted_mode = targeted_mode
        #posizione del frammento nel genoma
        self.__offset_forward = begin_seq
        self.__begin = begin_seq
        self.__end = end_seq
#        self.__offset_reverse = gsize - end_seq
        #dimensione buffer + parametri vari
        self.__buffer_size = seqparams.buffer_size
        self.__seqparams = seqparams
        self.__p_meth =  {
            CytosineContext.CG: seqparams.p_cg,
            CytosineContext.CHG: seqparams.p_chg,
            CytosineContext.CHH: seqparams.p_chh
        }
        #dati sulle citosine
        self.__cytosines = dict()
        self.__queue = queue

        self.__read_quality = dna.read_quality(seqparams.read_length, seqparams.min_quality, seqparams.max_quality)
        ##############################
        self.__set_snp()
        self.__initialize_cytosines()



    ###################### Sequencing ######################

    def fragmentation(self):
        """ """
        if not self.targeted_mode:
            max_step = 2*int(round(self.__seqparams.fragment_size_mean / self.__seqparams.coverage))
            i = 0

            while (i + self.__seqparams.fragment_size_mean) < len(self.__sequence):
                fragment = self.__sequence[i: i + self.__seqparams.fragment_size_mean]

                if len(fragment) == self.__seqparams.fragment_size_mean:
                    yield dna.fragment(fragment, i, i + self.__seqparams.fragment_size_mean)\
                            .methylate(self.methylate_cytosine)

                i += random.randint(1, max_step) #in media ogni base è coperta da C reads
        else:
            yield dna.fragment(self.__sequence, self.__begin, self.__end).methylate(self.methylate_cytosine)
            #for _ in range(0, self.__seqparams.coverage):
                # isize = 0
                # while isize <= 1:
                #     isize = norm.rvs(loc=self.__seqparams.fragment_size_mean, scale=self.__seqparams.fragment_size_sd, size=1)[0]
                # if isize <= 120:
                #     delta = 120 - isize
                #     left_remove = int(round(delta/2))
                #     new_begin = self.__begin + left_remove
                # else:
                #     delta = isize - 120
                #     left_add = int(round(delta/2))
                #     new_begin = self.__begin - left_add
                # new_end = new_begin + isize
                # fragment = self.__sequence[new_begin:new_end]
                # print(self.__chromoId, new_begin, new_end, isize)
              #  fragment = self.__sequence[self.__begin:self.__end]
              #  yield dna.fragment(fragment, self.__begin, self.__end).methylate(self.methylate_cytosine)


    def single_end_sequencing(self):
        """Produce le read single-end del frammento e le salva su un file fastq temporaneo"""

        fragment_size = self.__seqparams.fragment_size
        read_length = self.__seqparams.read_length
        directional = self.__seqparams.lib_mode == "directional"
        seq_length = len(self.__sequence)
        error_rate = self.__seqparams.error_rate

        reads = list()

        for fragment in self.fragmentation():
            #metilazione?
            bsfs = fragment.forward_strand.bisulfite()
            fastq_fw = bsfs.single_end_sequencing(read_length, self.__read_quality)\
                        .set_sequencing_errors(error_rate)\
                        .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

            bsrs = fragment.reverse_strand.bisulfite()
            fastq_rv = bsrs.single_end_sequencing(read_length, self.__read_quality)\
                        .set_sequencing_errors(error_rate)\
                        .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

            reads.extend([fastq_fw, fastq_rv])

            if not directional:
                fastq_fwrc = bsfs.reverse_complement()\
                            .single_end_sequencing(read_length, self.__read_quality)\
                            .set_sequencing_errors(error_rate)\
                            .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

                fastq_rvrc = bsrs.reverse_complement()\
                            .single_end_sequencing(read_length, self.__read_quality)\
                            .set_sequencing_errors(error_rate)\
                            .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

                reads.extend([fastq_fwrc, fastq_rvrc])


            if len(reads) > self.__buffer_size:
                self.__queue.put(("fastq_se", reads.copy()))
    #            self.__persist_record(reads)
                reads.clear()
        else:
            if len(reads) > 0:
                self.__queue.put(("fastq_se", reads.copy()))
    #            self.__persist_record(reads)
                reads.clear()

        cinfo = self.__format_methylation()
        self.__queue.put(("ch3", cinfo))

#        self.__persist_methylation()


    def paired_end_sequencing(self):
        """Produce le read paired-end del frammento e le salva su un file temporaneo"""

        fragment_size = self.__end - self.__begin
        read_length = self.__seqparams.read_length
        directional = self.__seqparams.lib_mode == "directional"
        seq_length = len(self.__sequence)
        error_rate = self.__seqparams.error_rate

        reads = list()

#        for (fragment, begin, end) in self.__fragmentation():
        for fragment in self.fragmentation():
            bsfs = fragment.forward_strand.bisulfite()
            fastq_fw = bsfs.paired_end_sequencing(read_length, self.__read_quality)\
                        .set_sequencing_errors(error_rate)\
                        .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

            bsrs = fragment.reverse_strand.bisulfite()
            fastq_rv = bsrs.paired_end_sequencing(read_length, self.__read_quality)\
                        .set_sequencing_errors(error_rate)\
                        .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

            reads.extend([fastq_fw, fastq_rv])



            if not directional:
                fastq_fwrc = bsfs.reverse_complement()\
                            .paired_end_sequencing(read_length, self.__read_quality)\
                            .set_sequencing_errors(error_rate)\
                            .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

                fastq_rvrc = bsrs.reverse_complement()\
                            .paired_end_sequencing(read_length, self.__read_quality)\
                            .set_sequencing_errors(error_rate)\
                            .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

                reads.extend([fastq_fwrc, fastq_rvrc])

            if len(reads) > self.__buffer_size:
                self.__queue.put(("fastq_pe", reads.copy()))
        #        self.__persist_record(reads)
                reads.clear()
        else:
            if len(reads) > 0:
                self.__queue.put(("fastq_pe", reads.copy()))
        #        self.__persist_record(reads)
                reads.clear()

        cinfo = self.__format_methylation()
        self.__queue.put(("ch3", cinfo))
    #    self.__persist_methylation()

    ###################### Methylation ######################

    def __initialize_cytosines(self):
        """Parserizza il genoma e indicizza le citosine sui due strand"""

        limit = len(self.__sequence)

        for pos, base in enumerate(self.__sequence):
            #strand +
            if base == 'c':
                context = CytosineContext.CHH
                if pos+1 < limit and self.__sequence[pos+1] == 'g':
                    context = CytosineContext.CG
                elif pos+2 < limit and self.__sequence[pos+2] == 'g':
                    context = CytosineContext.CHG

                self.__cytosines[pos] = Cytosine("+", context)
            #strand -
            elif base == 'g':
                context = CytosineContext.CHH
                if pos-1 >= 0 and self.__sequence[pos-1] == "c":
                    context = CytosineContext.CG
                elif pos-2 >= 0 and self.__sequence[pos-2] == "c":
                    context = CytosineContext.CHG

                self.__cytosines[pos] = Cytosine("-", context)


    def methylate_cytosine(self, base, position):
        state = base.upper()
        if position in self.__cytosines:
            cytosine = self.__cytosines[position]
            cytosine.covered()

            if random.uniform(0, 1) <= self.__p_meth[cytosine.context]:
                state = state.lower()
                cytosine.methylate()

        return state


    ###################### Introduzione mutazioni ######################

    def __set_snp(self):
        """Setta SNP random sulla reference"""

        self.__sequence = "".join([random.sample("actg", 1)[0] if random.uniform(0, 1) < self.__seqparams.snp_rate else base\
                                   for base in self.__sequence])


    def __format_methylation(self):
        """ Returns a list containing all the information about cytosines """
        beta_score = lambda c: 0 if c.ncov == 0 else c.nmeth / c.ncov

        return [(self.__chromoId, abs(self.__offset_forward+position), cytosine.strand, cytosine.context,\
                cytosine.nmeth, cytosine.ncov, beta_score(cytosine)) \
                    for position, cytosine in self.__cytosines.items()]


class TargetedFragmentSequencer(object):
    """..."""

    def __init__(self, chrom, orig_start, orig_end, frag_start, frag_end, isize, sequence, params, minq=10, maxq=40):
        for k, v in params.items():
            params.k = v
        #informazioni sulla sequenza
        self.__chromoId = chrom
        self.__sequence = sequence
        #posizione del frammento nel genoma
        self.__offset_forward = frag_start
        self.__begin = frag_start
        self.__end = frag_end
        self.__frag_size = isize
        #dimensione buffer + parametri vari
        self.__seqparams = params
        self.__p_meth =  {
            CytosineContext.CG: params['cg'],
            CytosineContext.CHG: params['chg'],
            CytosineContext.CHH: params['chh']
        }
        #dati sulle citosine
        self.__cytosines = dict()

        self.__read_quality = dna.read_quality(params.read_length, minq, maxq)
        ##############################
        self.__set_snp()
        self.__initialize_cytosines()


    ###################### Sequencing ######################

    def fragmentation(self):
        """ """
        yield dna.fragment(self.__sequence, self.__begin, self.__end).methylate(self.methylate_cytosine)

    def paired_end_sequencing(self):
        """Produce le read paired-end del frammento e le salva su un file temporaneo"""

        fragment_size = self.__frag_size
        read_length = self.__seqparams.read_length
        directional = True
        seq_length = len(self.__sequence)
        error_rate = self.__seqparams.error_rate

        reads = list()

        for fragment in self.fragmentation():
            bsfs = fragment.forward_strand.bisulfite()
            fastq_fw = bsfs.paired_end_sequencing(read_length, self.__read_quality)\
                        .set_sequencing_errors(error_rate)\
                        .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

            bsrs = fragment.reverse_strand.bisulfite()
            fastq_rv = bsrs.paired_end_sequencing(read_length, self.__read_quality)\
                        .set_sequencing_errors(error_rate)\
                        .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

            reads.extend([fastq_fw, fastq_rv])



            if not directional:
                fastq_fwrc = bsfs.reverse_complement()\
                            .paired_end_sequencing(read_length, self.__read_quality)\
                            .set_sequencing_errors(error_rate)\
                            .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

                fastq_rvrc = bsrs.reverse_complement()\
                            .paired_end_sequencing(read_length, self.__read_quality)\
                            .set_sequencing_errors(error_rate)\
                            .fastqize(self.__chromoId, fragment_size, self.__offset_forward)

                reads.extend([fastq_fwrc, fastq_rvrc])

        cinfo = self.__format_methylation()

        return reads, cinfo

    ###################### Methylation ######################

    def __initialize_cytosines(self):
        """Parserizza il genoma e indicizza le citosine sui due strand"""

        limit = len(self.__sequence)

        for pos, base in enumerate(self.__sequence):
            #strand +
            if base == 'c':
                context = CytosineContext.CHH
                if pos+1 < limit and self.__sequence[pos+1] == 'g':
                    context = CytosineContext.CG
                elif pos+2 < limit and self.__sequence[pos+2] == 'g':
                    context = CytosineContext.CHG

                self.__cytosines[pos] = Cytosine("+", context)
            #strand -
            elif base == 'g':
                context = CytosineContext.CHH
                if pos-1 >= 0 and self.__sequence[pos-1] == "c":
                    context = CytosineContext.CG
                elif pos-2 >= 0 and self.__sequence[pos-2] == "c":
                    context = CytosineContext.CHG

                self.__cytosines[pos] = Cytosine("-", context)


    def methylate_cytosine(self, base, position):
        state = base.upper()
        if position in self.__cytosines:
            cytosine = self.__cytosines[position]
            cytosine.covered()

            if random.uniform(0, 1) <= self.__p_meth[cytosine.context]:
                state = state.lower()
                cytosine.methylate()

        return state


    ###################### Introduzione mutazioni ######################

    def __set_snp(self):
        """Setta SNP random sulla reference"""

        self.__sequence = "".join([random.sample("actg", 1)[0] if random.uniform(0, 1) < self.__seqparams.snp_rate else base\
                                   for base in self.__sequence])


    def __format_methylation(self):
        """ Returns a list containing all the information about cytosines """
        beta_score = lambda c: 0 if c.ncov == 0 else c.nmeth / c.ncov

        return [(self.__chromoId, abs(self.__offset_forward+position), cytosine.strand, cytosine.context,\
                cytosine.nmeth, cytosine.ncov, beta_score(cytosine)) \
                    for position, cytosine in self.__cytosines.items()]