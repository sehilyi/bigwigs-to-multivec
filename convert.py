# -*- coding: utf-8 -*-

import pyBigWig
import math
from statistics import mean
import os
import os.path as op
import numpy as np
import h5py
import tempfile
import multivec as cmv

chromosome_order = [
    "chr1", "chr2", "chr3", "chr4", "chr5",
    "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20",
    "chr21", "chr22", "chrX", "chrY"
]

def _sort_by_chrom(x):
    """
    Sort a list by the predefined chromosome order

    Parameters
    ----------
    x: tuple (chromosome, size)
    """
    return chromosome_order.index(x[0]) if x[0] in chromosome_order else 999

def bigwig_to_multivec(
    input_files,
    starting_resolution=1
):
    """
    Convert a bigwig file to a multivec file.

    Parameters
    ----------
    input_files: array of strings of file paths
    starting_resolution: int (default 1)
        The starting resolution of the input data
    """

    if len(input_files) is 0:
        print("No enough input files defined.")
        return

    EMPTY_VALUE = np.nan

    with tempfile.TemporaryDirectory() as td:
        print("temporary dir:", td)

        temp_file = op.join(td, "temp.mv5")
        f_out = h5py.File(temp_file, "w")
        
        # Determine chromosome sizes by looking into the first file
        # TODO: 
        bw = pyBigWig.open(input_files[0])
        if not bw.isBigWig(): 
            print("Input file is not in a BigWig format:", input_files[0])
            return
        chromsizes = sorted([(k, v) for k, v in bw.chroms().items()], key=_sort_by_chrom)
        bw.close()

        for cur_input_file in input_files:
            bw = pyBigWig.open(cur_input_file)

            if not bw.isBigWig(): 
                print("One of input file is not in a BigWig format:", cur_input_file)
                return
            
            for chr_and_size in chromsizes:
                (chrom, size) = chr_and_size
                
                f_out.create_dataset(
                    chrom,
                    (
                        math.ceil(size / starting_resolution),
                        len([cur_input_file]) * len([cur_input_file]),
                    ),
                    fillvalue=EMPTY_VALUE,
                    compression="gzip",
                )
                
                data = []
                current_position = 0
                
                for interval in bw.intervals(chrom):
                    (interval_start, interval_end, value) = interval
                    
                    # Fill empty values when not suggested
                    if current_position < interval_start:
                        data += [[EMPTY_VALUE]] * (interval_start - current_position)

                    # Fill values with suggested value
                    data += [[value]] * (interval_end - interval_start)

                    current_position = interval_end
                
                # Fill empty values with the suggested value
                if current_position is not size:
                    data += [[EMPTY_VALUE]] * (size - current_position)
                
                f_out[chrom][0 : len(data)] = np.array(data)
                
                print(chrom, "value loaded...")

            f_out.close()
            tf = temp_file
            f_in = h5py.File(tf, "r")

            # if output_file is None:
            output_file = op.splitext([cur_input_file][0])[0] + ".multires.mv5"
            print("output_file:", output_file)

            # Override the output file if it existts
            if op.exists(output_file):
                os.remove(output_file)
            
            cmv.create_multivec_multires(
                f_in,
                chromsizes=chromsizes,
                agg=lambda x: np.nansum(x.T.reshape((x.shape[1], -1, 2)), axis=2).T,
                starting_resolution=starting_resolution,
                tile_size=1024,
                output_file=output_file
            )

        bw.close()

def main():
    bigwig_to_multivec(["sample_data/6_treat.bw", "sample_data/2_treat.bw", "sample_data/1_treat.bw"])

if __name__ == "__main__":
	main()