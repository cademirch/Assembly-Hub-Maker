import os
import subprocess
import argparse
import urllib.request
import glob
import re
import sys
from pathlib import Path


class Genome:
    def __init__(self, genome_name, path_to_genomes):
        self.genome = genome_name
        self.fa = f"{genome_name}_genomic.fna"
        self.gff = f"{genome_name}_genomic.gff"
        self.rm = f"{genome_name}_rm.out"
        self.path_to_genomes = path_to_genomes
        self.genome_path = f"{self.path_to_genomes}{self.genome}"

    def download(self):
        try:
            os.mkdir(f"{self.genome_path}")
        except FileExistsError:
            print("Error: Genome already exists")
            sys.exit(1)

        link_parts = self.genome.split('_', 2)
        letters = link_parts[0]
        numbers = link_parts[1]
        link = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{letters}/{numbers[0:3]}/{numbers[3:6]}/{numbers[6:9]}/{self.genome}"
        urllib.request.urlretrieve(f"{link}/{self.gff}.gz", f"{self.genome_path}/{self.gff}.gz")
        urllib.request.urlretrieve(f"{link}/{self.fa}.gz", f"{self.genome_path}/{self.fa}.gz")
        urllib.request.urlretrieve(f"{link}/{self.rm}.gz", f"{self.genome_path}/{self.rm}.gz")
        urllib.request.urlretrieve(f"{link}/{self.genome}_assembly_report.txt",
                                   f"{self.genome_path}/assembly_report.txt")
        os.system(f"gunzip {self.genome_path}/*.gz")

    def format(self):
        command = f'faToTwoBit {self.genome_path}/{self.fa} {self.genome_path}/{self.genome}.2bit'
        subprocess.Popen(command, shell=True).wait()
        command = f'twoBitInfo {self.genome_path}/{self.genome}.2bit {self.genome_path}/chrom.sizes'
        subprocess.Popen(command, shell=True).wait()
        command = f'gff3ToGenePred {self.genome_path}/{self.gff} {self.genome_path}/{self.genome}.genePred'
        subprocess.Popen(command, shell=True).wait()
        command = f'genePredToBigGenePred {self.genome_path}/{self.genome}.genePred {self.genome_path}/bigGenePred.txt'
        subprocess.Popen(command, shell=True).wait()
        command = f'bedSort {self.genome_path}/bigGenePred.txt {self.genome_path}/bigGenePred.bed'
        subprocess.Popen(command, shell=True).wait()
        command = f'bedClip {self.genome_path}/bigGenePred.bed {self.genome_path}/chrom.sizes ' \
                  f'{self.genome_path}/clippedBigGenePred.bed'
        subprocess.Popen(command, shell=True).wait()
        command = f'bedToBigBed -type=bed12+8 -tab -as={self.path_to_genomes}bigGenePred.as ' \
                  f'{self.genome_path}/clippedBigGenePred.bed ' \
                  f'{self.genome_path}/chrom.sizes {self.genome_path}/{self.genome}_genePred.bb'
        subprocess.Popen(command, shell=True).wait()

    def make_repeat_files(self):
        repeat_classes = ['SINE', 'LINE', 'LTR', 'DNA', 'Simple_repeat', 'Low_complexity', 'Satellite', 'RNA']
        int_dir = f"{self.genome_path}/intermediate_files"
        os.mkdir(int_dir)
        with open(f"{self.genome_path}/{self.rm}") as f:
            for _ in range(3):
                next(f)

            for line in f:
                line = line.split()

                sw_score, perc_div, perc_del, perc_ins, chrom, begin, end, geno_left, strand = \
                    line[0], float(line[1]), float(line[2]), float(line[3]), line[4], line[5], line[6], line[7], line[8]
                name, family, rep_begin, rep_end, rep_left, = \
                    line[9], line[10], line[11], line[12], line[13]
                geno_left = geno_left.strip('()')
                rep_left = rep_left.strip('()')
                rep_begin = rep_begin.strip('()')
                if "/" in family:
                    repeat_class = family.split("/")[0]
                    repeat_family = family.split("/")[1]

                else:
                    repeat_class = family
                    repeat_family = family

                filename = None

                for repeat_type in repeat_classes:
                    if repeat_type in repeat_class:

                        filename = f'{self.genome_path}/{self.genome}.rmsk.{repeat_type}.tab'

                if filename is None:
                    filename = f'{self.genome_path}/{self.genome}.rmsk.Other.tab'

                if strand is 'C':
                    strand = '-'

                perc_div = int(perc_div * 10)
                perc_del = int(perc_del * 10)
                perc_ins = int(perc_ins * 10)

                bed_line = f"{chrom}\t{begin}\t{end}\t{name}\t0\t{strand}\t{sw_score}\t{perc_div}\t{perc_del}\t" \
                           f"{perc_ins}\t{geno_left}\t{repeat_class}\t{repeat_family}\t{rep_begin}\t{rep_end}" \
                           f"\t{rep_left}\n"
                with open(filename, 'a+') as inp:
                    inp.write(bed_line)

        tab_files = glob.glob(f"{self.genome_path}/*.tab")

        for tab_file in tab_files:
            bed_file = f"{tab_file.rsplit('.', 1)[0]}.bed"
            command = f"sort -k1,1 -k2,2n {tab_file} > {bed_file}"
            subprocess.Popen(command, shell=True).wait()
            bbi_file = f"{tab_file.rsplit('.', 1)[0]}.bb"
            command = f"bedToBigBed -tab -type=bed6+10 -as={self.path_to_genomes}/rmskBed6+10.as " \
                      f"{bed_file} " \
                      f"{self.genome_path}/chrom.sizes " \
                      f"{bbi_file}"
            subprocess.Popen(command, shell=True).wait()
            os.rename(tab_file, f"{int_dir}/{tab_file.rsplit('/', 1)[1]}")
            os.rename(bed_file, f"{int_dir}/{bed_file.rsplit('/', 1)[1]}")

    def make_track_db(self):
        repeat_classes = ['SINE', 'LINE', 'LTR', 'DNA', 'Simple_repeat', 'Low_complexity', 'Satellite', 'RNA', 'Other']
        rmsk_files = glob.glob(f"{self.genome_path}/*.rmsk.*.bb")

        gene_text_block = f"track ncbiGene\n"\
                          f"longLabel ncbiGene - gene predictions delivered with assembly from NCBI\n" \
                          f"shortLabel ncbiGene\n" \
                          f"priority 12\n" \
                          f"visibility pack\n" \
                          f"color 0,80,150\n" \
                          f"altColor 150,80,0\n" \
                          f"colorByStrand 0,80,150 150,80,0\n" \
                          f"bigDataUrl {self.genome}_genePred.bb\n" \
                          f"type bigGenePred\n" \
                          f"group genes\n"

        rmsk_text_block = f"\ntrack RepeatMasker\n" \
                          f"compositeTrack on\n" \
                          f"shortLabel RepeatMasker\n" \
                          f"longLabel Repeating Elements by RepeatMasker\n" \
                          f"group varRep\n" \
                          f"priority 149.1\n" \
                          f"visibility dense\n" \
                          f"type bed 3 .\n" \
                          f"noInherit on\n"

        with open(f"{self.genome_path}/trackDb.txt", "a+") as f:
            f.write(gene_text_block)
            f.write(rmsk_text_block)

            for file in rmsk_files:
                file = f"{file.rsplit('/', maxsplit=1)[1]}"

                for repeat_class in repeat_classes:
                    if repeat_class in file:
                        text_block = f"\n\ttrack RepeatMasker{repeat_class}\n" \
                                     f"\tparent RepeatMasker\n" \
                                     f"\tshortLabel {repeat_class}\n" \
                                     f"\tlongLabel {repeat_class} Repeating Elements by RepeatMasker\n" \
                                     f"\tpriority 1\n" \
                                     f"\tspectrum on\n" \
                                     f"\tmaxWindowToDraw 10000000\n" \
                                     f"\tcolorByStrand 50,50,150 150,50,50\n" \
                                     f"\ttype bigBed 6 +\n" \
                                     f"\tbigDataUrl {file}\n"

                        f.write(text_block)

    def write_to_genome(self):
        command = f"sort -n -k2,2 -r {self.genome_path}//chrom.sizes > {self.genome_path}/sorted.chrom.sizes"
        subprocess.Popen(command, shell=True).wait()
        with open(f"{self.genome_path}/sorted.chrom.sizes") as f:
            line = f.readline().split()
            chrom = line[0]
            length = int(line[1])
            default_pos = f"{chrom}:{length-100000}-{length}"

        with open(f"{self.genome_path}/assembly_report.txt") as l:
            for line in l:
                if line.startswith('#'):

                    if 'Organism' in line:
                        common_name = re.search('\(([^)]+)', line).group(1)
                        line = re.split(r':+|\(+', line)
                        scientific_name = line[1].strip()
                    if 'Date:' in line:
                        line = line.split()
                        date = line[2]

        text = f"\ngenome {self.genome}\n" \
               f"trackDb {self.genome}/trackDb.txt\n" \
               f"groups groups.txt\n" \
               f"description {date} {common_name}\n" \
               f"twoBitPath {self.genome}/{self.genome}.2bit\n" \
               f"organism {scientific_name}\n" \
               f"defaultPos {default_pos}\n"

        with open(f"{self.path_to_genomes}/genomes.txt", 'a+') as h:
            h.write(text)

class CommandLine:

    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description='Creates assembly hub with ncbi gene predictions and repeatMasker tracks. ',
            add_help=True,  # default is True
            prefix_chars='-',
        )

        self.parser.add_argument('-g', '--accession', help='refSeq accession of genome', required=True, type=str,
                                 dest='accession')
        self.parser.add_argument('-p', '--path', help='path where genomes.txt is', required=True, dest='path')
        self.parser.add_argument('-s', '--skip_dl', help='Skip download step', required=False
                                 , dest='skip', action='store_true')

def main():

    myCommandLine = CommandLine()
    args = myCommandLine.parser.parse_args()
    myGenome = Genome(args.accession, args.path)

    if os.getcwd() != args.path:
        os.chdir(args.path)

    if not args.skip:
        myGenome.download()
        myGenome.format()
        myGenome.make_repeat_files()
        myGenome.make_track_db()
        myGenome.write_to_genome()
    else:
        myGenome.format()
        myGenome.make_repeat_files()
        myGenome.make_track_db()
        myGenome.write_to_genome()
if __name__ == "__main__":
    main()
