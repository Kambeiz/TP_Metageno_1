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

__author__ = "Debbah Nagi"
__copyright__ = "Universite de Paris"
__credits__ = ["Debbah Nagi"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Debbah Nagi"
__email__ = "debbah.nagi@gmail.com"
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
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """
    open a gz or fasta file and create a generator from that fasta file 
    according the minimal sequence length chosen
    Parameters
    ----------
    amplicon_file : string
    minseqlen : integer
    
    Yields
    ------
    generator
    """
    if amplicon_file.endswith("gz"): 
        with gzip.open(amplicon_file, "rb") as gz:
            seq = b""
            for line in gz:
                if line.startswith(b">"):
                    if len(seq) >= minseqlen:
                        yield seq.decode('ascii')
                    seq = b""
                else:
                    seq += line.strip()
            yield seq

    elif amplicon_file.endswith("fasta"):
        with open(amplicon_file, "r") as fasta:
            seq = ""
            for line in fasta:
                if line.startswith('>'):
                    if len(seq) >= minseqlen:
                         yield seq
                    seq =""
    
                else:
                    seq = seq + line[:-1]
            if len(seq) >= minseqlen:
                yield seq

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """
    sort sequences according their decreasing count (with minimal thresold 
    being mincount)
    Parameters
    ----------
    amplicon_file : String
    minseqlen : integer
    mincount : integer

    Yields : Generator
    ------
    """

    seq_list = []
    for seq in read_fasta(amplicon_file, minseqlen):
        seq_list.append(seq)

    for o in Counter(seq_list).most_common():
        if o[1] > mincount:
            yield o

def get_chunks(sequence, ck_size):
    """
    Cut a sequence in chunk according a chunk size being at the lowest, 4    

    Parameters
    ----------
    sequence : integer
    ck_size : integer

    Raises
    ------
    ValueError

    Returns
    -------
    list_chunk : list

    """
    
    list_chunk = []
    i=1
    l = len(sequence)
    if l < 4*ck_size:
        raise ValueError("Chunk size should be of 4 at least ")
    for i in range(1, l):
        if i*ck_size < l:
            list_chunk.append(sequence[i*ck_size-ck_size:i*ck_size])
    #while(i*ck_size < l):
        #list_chunk.append(sequence[i*ck_size-ck_size:i*ck_size])
        #i += 1
    return list_chunk

def cut_kmer(sequence, k_mer):
    """
    Take a sequence and cut it according a size given into k_mer.
    Parameters
    ----------
    sequence : Fasta.
        Nucleotidic sequence.
    k_mer : Integer
        Size of the k_mer chosen.

    Returns
    -------
    Iterator.
        Iterator containing all k_mer possible from the sequnece.

    """
    for i in range(0, len(sequence)-k_mer + 1):
        yield sequence[i:i+k_mer]

def get_unique(ids):
    return {}.fromkeys(ids).keys()

def get_unique_kmer(kmer_dict, sequence, seq_id, km_size):
    """
    Create a dictionnary containing a kmer and his position in differents 
    sequences

    Parameters
    ----------
    kmer_dict : dict
    sequence : string
    seq_id : integer
    km_size : integer

    Returns
    -------
    kmer_dict : dict

    """
    for seq in cut_kmer(sequence, km_size):
        if seq not in kmer_dict:
            kmer_dict[seq] = [seq_id]
        elif seq_id not in kmer_dict[seq]:
            kmer_dict[seq].append(seq_id)
    return kmer_dict

def search_mates(kmer_dict, sequence, km_size):    

    return [i[0] for i in Counter([ids for kmer in 
                                   cut_kmer(sequence, km_size) if kmer in
                                   kmer_dict for ids in 
                                   kmer_dict[kmer]]).most_common(8)]

def get_identity(alignement_list):
    """
    Calculate the percent of identity between two sequences in an alignment

    Parameters
    ----------
    alignement_list : list

    Returns
    -------
    float

    """
    nb_base =0
    l = len(alignement_list[0])
    for i in range(0, l):
        if alignement_list[0][i] == alignement_list[1][i]:
            nb_base += 1
    return nb_base/l * 100

def detect_chimera(perc_identity_matrix):
    """
    detect a chimera sequence

    Parameters
    ----------
    perc_identity_matrix : matrix 
        matrix containing percent of identity per segments of sequences

    Returns
    -------
    bool

    """
    std_list = []
    bool_s0 = False
    bool_s1 = False

    for lis in perc_identity_matrix:
        std_list.append(statistics.stdev(lis))
        if lis[0] > lis[1]:
            bool_s0 = True
        if lis[0] < lis[1]:
            bool_s1 = True
                    
    if statistics.mean(std_list) > 5 and bool_s0 and bool_s1:
        return True
    else:
        return False


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    remove chimera sequences

    Parameters
    ----------
    amplicon_file : string
    minseqlen : int
    mincount : int
    chunk_size : int
    kmer_size : int

    Yields
    ------
    sequence : str
    """

    list_nonchim = []
    kmer_dict = {}
    chimera = False

    for i, sequence in enumerate(dereplication_fulllength(amplicon_file, minseqlen, mincount)):
        
        ck_list = get_chunks(sequence[0], chunk_size)
        
        mates_list = []
        for chunk in ck_list:
            mates_list.append(search_mates(kmer_dict, chunk, kmer_size))

        common_list = []
        for i in range(len(mates_list)):
            common_list = common(common_list, mates_list[i])


        if len(common_list) > 1:
            for c in common_list[0:2]:
                ck_ref = get_chunks(list_nonchim[c], chunk_size)
                
                identity_matrix = [[]*4]
                for i in range(len(ck_list)):
                    align = nw.global_align(ck_list[i], ck_ref[i], gap_open=-1, gap_extend=-1, matrix="MATCH")
                    identity_matrix[i].append(get_identity(align))

        if not chimera:
            kmer_dict = get_unique_kmer(kmer_dict, sequence[0], i, kmer_size)
            list_nonchim.append(sequence[0])
            yield sequence



def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))

def fill(text, width=80):
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    return operational taxonomic unit (otu)

    Parameters
    ----------
    amplicon_file : string
    minseqlen : int
    mincount : int
    chunk_size : int
    kmer_size : int

    Returns
    -------
    otu : list

    """
    otu = []

    for seq, number in chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
        otu.append((seq, number))

    return otu



def write_OTU(otu_list, output_file):
    """
    Parameter
    ---------
    otu_list: list
    output_file: path
    """
    with open(output_file, "w") as filout:
        for i, otu in enumerate(otu_list):
            filout.write(">OTU_{} occurrence:{}\n".format(i+1, otu[1]))
            filout.write("{}\n".format(fill(otu[0])))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    otu = abundance_greedy_clustering(args.amplicon_file, args.minseqlen,
                                      args.mincount, args.chunk_size, args.kmer_size)
    write_OTU(otu, args.output_file)

if __name__ == '__main__':
    main()