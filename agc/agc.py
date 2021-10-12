#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    with gzip.open(amplicon_file, "rt") as filin:
        content = ""
        for line in filin:
            if not line.startswith(">"):
                content += line.strip()
            else:
                if len(content) >= minseqlen:
                    yield content
                content = ''
        yield content

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    gen_seq = list(read_fasta(amplicon_file, minseqlen))

    seq_count = [(seq, gen_seq.count(seq)) for seq in gen_seq]
    list_seq = list(get_unique(seq_count))
    list_seq.sort(key = lambda seqlen: seqlen[1])

    for seq in list_seq[::-1]:
        if seq[1] >= mincount:
            yield seq

def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """"""
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size] 
              for i in range(0, len_seq, chunk_size) 
                if i+chunk_size <= len_seq - 1]

def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    kmers = cut_kmer(sequence, kmer_size)

    for kmer in kmers:
        if kmer not in kmer_dict:
            kmer_dict[kmer] = [id_seq]
        else:
            kmer_dict[kmer] += [id_seq]
    return kmer_dict

def search_mates(kmer_dict, sequence, kmer_size):
    kmer_seq = cut_kmer(sequence, kmer_size)
    count_id = []
    for kmer in kmer_seq:
        if kmer in kmer_dict:
            count_id += kmer_dict[kmer]
    parents = Counter(count_id).most_common(2)
    id_parent = [id for id, occ in parents]
    return id_parent

def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    sequence_length = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    nb_seq = len(sequence_length)
    list_otu = []
    mother = []
    for i in range(nb_seq):
        if sequence_length[i] in mother:
            continue
        inter = [sequence_length[i]]            
        for j in range(i+1, nb_seq):
            align = nw.global_align(sequence_length[i][0], sequence_length[j][0], gap_open=-1, gap_extend=-1,
            matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            identity = get_identity(align)
            print(identity)
            if identity > 97: 
                inter.append(sequence_length[j])
        list_otu.append(inter[0])                
        mother += inter
    return list_otu

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    with open(output_file, "w") as filout:
        count = 1
        for seqlen in OTU_list:
            filout.write(f">OTU_{count} occurrence:{seqlen[1]}\n")
            full_seq = fill(seqlen[0])
            filout.write(f"{full_seq}\n")
            count += 1

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """

    # Get arguments
    args = get_arguments()
    # Votre programme ici
    otu = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
    print(otu)
if __name__ == '__main__':
    main()
