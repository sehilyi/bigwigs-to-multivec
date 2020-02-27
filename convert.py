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

def bigwig_to_multivec(
    input_file,
    starting_resolution=1
):
    """
    Convert a bigwig file to a multivec file.

    Parameters
    ----------
    input_file: str
    starting_resolution: int (default 1)
        The starting resolution of the input data
    """

    bw = pyBigWig.open(input_file)

    if not bw.isBigWig(): 
        print("Input file is not in a BigWig format.")

    EMPTY_VALUE = np.nan

    chroms = bw.chroms()

    with tempfile.TemporaryDirectory() as td:
        print("temporary dir:", td)

        temp_file = op.join(td, "temp.mv5")
        f_out = h5py.File(temp_file, "w")
        chromsizes = []

        for chrom in [*chroms]:
            # chrom = "chr20" # TODO: Test with only one chromosome
            f_out.create_dataset(
                chrom,
                (
                    math.ceil(chroms[chrom] / starting_resolution),
                    len([input_file]) * len([input_file]),
                ),
                fillvalue=EMPTY_VALUE,
                compression="gzip",
            )
            
            chromsizes += [(chrom, chroms[chrom])]
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
            
            # Fill empty values with suggested value
            if current_position is not chroms[chrom]:
                data += [[EMPTY_VALUE]] * (chroms[chrom] - current_position)
            
            f_out[chrom][0 : len(data)] = np.array(data)
            
            print(chrom, "value loaded...")
            
            # break   # TODO: Test with only one chromosome

        f_out.close()
        tf = temp_file
        f_in = h5py.File(tf, "r")

        # if output_file is None:
        output_file = op.splitext([input_file][0])[0] + ".multires.mv5"
        print("output_file:", output_file)

        # Override the output file if it existts
        if op.exists(output_file):
            os.remove(output_file)
        
        cmv.create_multivec_multires(
            f_in,
            chromsizes=chromsizes,#[('chr20', 64444167)],#chroms,
            agg=lambda x: np.nansum(x.T.reshape((x.shape[1], -1, 2)), axis=2).T,
            starting_resolution=starting_resolution,
            tile_size=1024,
            output_file=output_file
        )

    bw.close()

def main():
    bigwig_to_multivec("sample_data/6_treat.bw")

if __name__ == "__main__":
	main()