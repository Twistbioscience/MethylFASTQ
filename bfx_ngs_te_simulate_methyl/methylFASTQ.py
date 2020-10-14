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

import argparse, sys
import os, pathlib
import random, string
import csv
from Bio import SeqIO
import bfx_ngs_te_simulate_methyl.sequencing as seq
from pybedtools import BedTool
import multiprocessing as mp
from scipy.stats import norm, truncnorm
import subprocess
import bfx_ngs_te_simulate_methyl.repair as repair

from timeit import default_timer as timer


class MethylFASTQ(object):
    def __init__(self, args):
        self.params = self.__parse_args(args)
        self.regions = self.__load_regions(args.regions)
        self.__stats = seq.Stats()

        self.targeted_mode = True

    def run(self):
        self.targeted_mode = True if self.regions else False
        print(f"Running in targeted sequencing mode: {self.targeted_mode}")

        start = timer()

        self.__run_targeted() if self.targeted_mode else self.__run_wgbs()

        elapsed = timer() - start

        print("Num reads: {}\nNum C: {}\nNum bases: {}\n".format(self.__stats.nreads, self.__stats.ncytosines, self.__stats.nbases))

        print("Elapsed time: {}".format(seq.format_time(elapsed)))


    def __parse_args(self, args):
        #check esistenza directory output
        if not os.path.exists(args.output_path):
            print("Output directory {} does not exist.".format(args.output_path), end=" ")
            os.makedirs(args.output_path)
            print("Now it does.")

        if args.fragment_size_mean == 0:
            args.fragment_size_mean = args.read_length

        if args.num_processes < 1:
            raise ValueError("Number of processes must be greater than zero.")

        return args

    def __run_wgbs(self):
        num_reads, num_c = 0, 0

        for fasta_record in SeqIO.parse(self.params.fasta_file, "fasta"):
            if self.params.chr is None or fasta_record.id in self.params.chr:
                print("Sequencing {}: {} bp".format(fasta_record.id, len(fasta_record)), flush=True)

                stats, r1, r2 = seq.ChromosomeSequencer(fasta_record, self.params).sequencing(self.params)
                self.__stats.update(stats)

        read1_paired_fastq, read2_paired_fastq, read1_unpaired_fastq, read2_unpaired_fastq = repair.repair_fastqs_files(r1, r2)

        for r in [r1, r2]:
            cmd = f'pigz {r}'
            subprocess.check_call(cmd, shell=True)

    def __run_targeted(self):
        for fasta_record in SeqIO.parse(self.params.fasta_file, "fasta"):
            if fasta_record.id in self.regions:
                curr_chr = self.regions[fasta_record.id]
                #print("{} --> {}".format(fasta_record.id, curr_chr))

                stats, r1, r2 = seq.ChromosomeSequencer(fasta_record, params=self.params, target_regions=curr_chr)\
                                   .sequencing(self.params)
                self.__stats.update(stats)

        read1_paired_fastq, read2_paired_fastq = repair.repair_fastqs_files(r1, r2)

        for r in [read1_paired_fastq, read2_paired_fastq]:
            cmd = f'pigz -f {r}'
            subprocess.check_call(cmd, shell=True)

    def __load_regions(self, filename):
        if filename is None or not os.path.exists(filename):
            msg = "ERROR: Regions BED file not supplied or does not exist!"
            raise FileNotFoundError(msg)
        print(f'Loading targeted regions to coverage depth: {self.params.coverage}...')
        target_regions = dict()
        regions = BedTool(filename)
        regions_sorted = regions.sort()
        # regions_sorted_merged = regions_sorted.merge(d=0)
        pd_df = regions_sorted.to_dataframe(names=['chrom', 'start', 'end'])
        probe_count = len(pd_df)
        regions_count = 0
        for index, row in pd_df.iterrows():
            chromosome = row["chrom"]
            for _ in range(self.params.coverage):
                isize = 0
                while isize <= 1:
                    isize = norm.rvs(loc=self.params.fragment_size_mean, scale=self.params.fragment_size_sd, size=1)[0]
                    isize = int(round(isize))
                delta = isize - self.params.probe_length
                left_pad = int(round(delta / 2))
                new_start = int(row["start"]) + left_pad
                new_end = new_start + isize
                #interval = int(row["start"]), int(row["end"])
                interval = new_start, new_end
                if chromosome not in target_regions:
                    target_regions[chromosome] = list()
                target_regions[chromosome].append(interval)
                regions_count += 1
            print(f'Completed probe {index+1} of {probe_count}')
        assert probe_count*self.params.coverage == regions_count
        print(f'Probe count:{probe_count} * Coverage:{self.params.coverage} == Total regions_count:{regions_count}')
        return target_regions


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="methyl-fastq")

    #I/O parameters: input, output
    parser.add_argument("-i", "--in", dest="fasta_file", metavar="fasta-file",
                        action="store", type=str, required=True,
                        help="Path of the FASTA file containing the genome to be sequenced.")
    parser.add_argument("-o", "--out", dest="output_path", metavar="output-file-path",
                        action="store", type=str, required=True,
                        help="Path of output files (.fastq & .cpg)")

    #sequencing mode and library mode
    parser.add_argument("--seq", dest="seq_mode",
                        choices=["single_end", "paired_end"], default="paired_end",
                        help="Type of reads to produce")

    parser.add_argument("--lib", dest="lib_mode",
                        choices=["directional", "non_directional"], default="directional",
                        help="Library type to produce")
    #chromosomes or regions to be sequenced
    parser.add_argument("--chr", dest="chr", metavar="chromosome-id",
                        nargs="*", action="store", type=str,
                        help="List of chromosomes to be sequenced (WGBS)")

    parser.add_argument("--regions", dest="regions", metavar="target-regions",
                        action="store", type=str, default=None,
                        help="Genomic regions to be sequenced (targeted bisulfite sequencing)")
    #coverage
    parser.add_argument("--coverage", dest="coverage", metavar="coverage",
                        action="store", type=int, default=250,
                        help="Depth of coverage")

    #fragment and read length and adapter seqs
    parser.add_argument("--fragment_mean", dest="fragment_size_mean", metavar="fragment-size-mean",
                        action="store", type=int, default=200,
                        help="Mean size of fragments in the fragmenation step")

    parser.add_argument("--fragment_sd", dest="fragment_size_sd", metavar="fragment-size-sd",
                        action="store", type=int, default=70,
                        help="Standard deviation of size of fragments in the fragmenation step")

    parser.add_argument("--read", dest="read_length", metavar="read-length",
                        action="store", type=int, default=76,
                        help="Read length of sequencing step")

    parser.add_argument("--probe", dest="probe_length", metavar="probe-length",
                        action="store", type=int, default=120,
                        help="TE probe length")
    #parallel option
    parser.add_argument("--processes", dest="num_processes", metavar="num-processes",
                        action="store", type=int, default=mp.cpu_count(),
                        help="Number of processes to be used during sequencing step")

    parser.add_argument("--buffer", dest="buffer_size", metavar="buffer-size",
                        action="store", type=int, default=10**5,
                        help="Buffer size of each process during sequencing step")
    #methylation probabilities
    parser.add_argument("--cg", dest="p_cg", metavar="CG-methylation-probability",
                        action="store", type=float, default=1.0,
                        help="Methylation probability in CG context")

    parser.add_argument("--chg", dest="p_chg", metavar="CHG-methylation-probability",
                        action="store", type=float, default=0.00,
                        help="Methylation probability in CHG context")

    parser.add_argument("--chh", dest="p_chh", metavar="CHH-methylation-probability",
                        action="store", type=float, default=0.00,
                        help="Methylation probability in CHH context")

    #sequencing error and snp probabilities
    parser.add_argument("--snp", dest="snp_rate", metavar="snp-probability",
                        action="store", type=float, default=0,
                        help="Probability to set a SNP")

    parser.add_argument("--error", dest="error_rate", metavar="sequencing-error-probability",
                        action="store", type=float, default=0,
                        help="Probability to set a sequencing error")
    #quality of reads
    parser.add_argument("--maxq", dest="max_quality", metavar="max-phred-score",
                        action="store", type=int, default=40,
                        help="Max phred score in the reads")

    parser.add_argument("--minq", dest="min_quality", metavar="min-phred-score",
                        action="store", type=int, default=10,
                        help="Min phred score in the reads (not implemented)") #not implemented yet!!


    args = parser.parse_args()

    MethylFASTQ(args).run()
