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
from nltk.tag import pos_tag
from numpy.lib.shape_base import _put_along_axis_dispatcher
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
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h".format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default=400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default=10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default=100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default=8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    """
    Reads a fasta file and makes a generator containing the sequences having a certain length found in the file
    :param amplicon_file: name or path of the file
    :param minseqlen: the minimum length of the sequences to be used
    :return: a generator of sequences
    """
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
    """
    This function counts the occurences of sequences and return sequences depending on their occurence
    :param amplicon_file: name or path of the file
    :param minseqlen: the minimum length of the sequences to be used
    :param mincount: the minimum number of occurences to be used
    :return: returns a generator with the sequences and their number of occurence
    """
    gen_seq = list(read_fasta(amplicon_file, minseqlen))
    dico_count = {}

    for seq in gen_seq:
        if seq not in dico_count:
            dico_count[seq] = 1
        else:
            dico_count[seq] += 1
    seq_count = dico_count.items()

    list_seq = list(get_unique(seq_count))
    list_seq.sort(key=lambda seqlen: seqlen[1])

    for seq in list_seq[::-1]:
        if seq[1] >= mincount:
            yield seq


def get_unique(ids):
    """
    function used to erase the duplications in a list
    :param ids: a list within which replicates must be erased
    :return: a dereplicated dict
    """
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    """
    return common elements of lists
    :param lst1: a list
    :param lst2: a list
    :return: common elements of lists
    """
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """
    Gets chunks from a sequence
    :param sequence: sequence of interest
    :param chunk_size: size of the chunks
    :return: chunks
    """
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i + chunk_size]
            for i in range(0, len_seq, chunk_size)
            if i + chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """
    Cut sequences into kmers
    """
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i + kmer_size]


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """
    Add kmers and the id of the sequence within which it was found in a dictionnary if the kmer is not in the
    dictionnary already, otherwise it will simply add the id of the sequence to the kmer (key) in the dictionnary
    :param kmer_dict: a dictionnary of kmers
    :param sequence: sequence of interest
    :param id_seq: sequence identifier
    :param kmer_size: size of kmers
    :return: a dictionnary of kmers
    """
    kmers = cut_kmer(sequence, kmer_size)
    for kmer in kmers:
        if kmer not in kmer_dict:
            kmer_dict[kmer] = [id_seq]
        else:
            kmer_dict[kmer] += [id_seq]
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    """
    Searches for parents of a given sequence
    :param kmer_dict: a dictionnary of kmers
    :param sequence: sequence of interest
    :param kmer_size: size of kmers
    :return: the parents identifier
    """
    kmer_seq = cut_kmer(sequence, kmer_size)
    count_id = []
    for kmer in kmer_seq:
        if kmer in kmer_dict:
            count_id += kmer_dict[kmer]
    if len(count_id) > 1:
        parents = Counter(count_id).most_common(2)
        id_parent = [id for id, occ in parents]
        return id_parent
    return count_id


def std(data):
    """
    gets standard deviation
    :param data: data
    :return: standard deviation of the data
    """
    return statistics.stdev(data)


def detect_chimera(perc_identity_matrix):
    """
    Detecting chimera sequences
    :param perc_identity_matrix: matrix of of identity
    :return: a boolean true if chimera, false otherwise
    """
    std_seg = []
    parent1 = False
    parent2 = False

    for seg in perc_identity_matrix:
        std_seg.append(std(seg))
        if seg[0] >= seg[1]:
            parent1 = True
        else:
            parent2 = True
    avg_std = statistics.mean(std_seg)

    if avg_std > 5:
        if parent1 and parent2:
            return True
    return False


def get_identity(alignment_list):
    """Prend en une liste de s??quences align??es au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    Removing chimeras from the dataset
    :param amplicon_file: name or path of the file
    :param minseqlen: the minimum length of the sequences to be used
    :param mincount: the minimum number of occurences to be used
    :param chunk_size: size of the chunks
    :param kmer_size: size of kmers
    :return: a generator of sequences without chimeras
    """
    all_seq = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    nb_seq = len(all_seq)

    kmer_dict = get_unique_kmer({}, all_seq[0][0], 0, kmer_size)
    kmer_dict = get_unique_kmer(kmer_dict, all_seq[1][0], 1, kmer_size)
    yield all_seq[0]
    yield all_seq[1]
    file_match = os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH"))
    for i in range(2, nb_seq):
        pot_parent = []
        perc_identity_matrix = []
        seq_chunk = get_chunks(all_seq[i][0], chunk_size)
        parent1_chunk = []
        parent2_chunk = []

        for chunk in seq_chunk:
            pot_parent += search_mates(kmer_dict, chunk, kmer_size)

        most_parent = Counter(pot_parent).most_common(2)

        parent1_chunk += get_chunks(all_seq[most_parent[0][0]][0],chunk_size)
        parent2_chunk += get_chunks(all_seq[most_parent[1][0]][0],chunk_size)
        for j in range(4):
            align_parent1 = nw.global_align(seq_chunk[j], parent1_chunk[j],
                gap_open=-1, gap_extend=-1, matrix=file_match)
            identity_parent1 = get_identity(align_parent1)

            align_parent2 = nw.global_align(seq_chunk[j], parent2_chunk[j],
                gap_open=-1, gap_extend=-1, matrix=file_match)
            identity_parent2 = get_identity(align_parent2)

            perc_identity_matrix.append([identity_parent1,identity_parent2])

        if not detect_chimera(perc_identity_matrix):
            kmer_dict = get_unique_kmer(kmer_dict, all_seq[i][0], i, kmer_size)
            yield all_seq[i]


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    Finding otus using a greedy algorithm
    :param amplicon_file: name or path of the file
    :param minseqlen: the minimum length of the sequences to be used
    :param mincount: the minimum number of occurences to be used
    :param chunk_size: size of the chunks
    :param kmer_size: size of kmers
    :return: a list of otus
    """
    sequence_length = list(chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size))
    nb_seq = len(sequence_length)
    list_otu = [sequence_length[0]]
    file_match = os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH"))
    for i in range(1,nb_seq):
        non_sim = True
        for j in range(len(list_otu)):
            align = nw.global_align(sequence_length[i][0], list_otu[j][0], gap_open=-1, gap_extend=-1,
            matrix=file_match)
            identity = get_identity(align)
            if identity > 97:
                non_sim = False
                break
        if non_sim:
            list_otu.append(sequence_length[i])

    return list_otu


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i + width] for i in range(0, len(text), width))


def write_OTU(OTU_list, output_file):
    """
    Making a file with all otus
    :param OTU_list: otus list
    :param output_file: name of the output file
    :return: a file
    """
    with open(output_file, "w") as filout:
        count = 1
        for seqlen in OTU_list:
            filout.write(f">OTU_{count} occurrence:{seqlen[1]}\n")
            full_seq = fill(seqlen[0])
            filout.write(f"{full_seq}\n")
            count += 1


# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """

    # Get arguments
    args = get_arguments()
    # Votre programme ici
    seq_otu = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size,
                                          args.kmer_size)
    write_OTU(seq_otu, args.output_file)


if __name__ == '__main__':
    main()
