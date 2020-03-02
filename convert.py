# -*- coding: utf-8 -*-

import os
import os.path as op
import pyBigWig
import math
import h5py
import tempfile
from statistics import mean
import numpy as np
import multivec as cmv
import utils

IS_DEBUG = False
test_input_files_1 = ["sample_data/2_treat.bw"]
test_input_files_2 = ["sample_data/1_treat.bw", "sample_data/2_treat.bw"]
test_input_files_3 = [
    "sample_data/6_treat.bw", 
    "sample_data/2_treat.bw", 
    "sample_data/1_treat.bw"
]
test_input_files_58 = ["sample_output/" + str(i + 1) + "_cistrome_db.bw" for i in range(2)]

def bigwig_to_multivec(
    input_files,
    starting_resolution=1   # TODO: Enable accepting non-one resolution
):
    """
    Convert a bigwig file to a multivec file.

    Parameters
    ----------
    input_files: array of file paths
    starting_resolution: int (default 1)
        The starting resolution of the input data
    """

    # TODO: May need to use batch to be able to deal with large data in unstable h5py package.

    ## Handling Errors.
    # None input file.
    if len(input_files) is 0:
        print("No enough input files suggested.")
        return
    print(len(input_files), "files prepared.")

    # None bigwig files.
    for in_file in input_files:
        bw = pyBigWig.open(in_file)
        if not bw.isBigWig(): 
            print("Input files are not in a BigWig format:", in_file)
            bw.close()
            return
        bw.close()

    # Define const variables.
    EMPTY_VALUE = np.nan

    ## Convert
    # Store data in a cache file.
    with tempfile.TemporaryDirectory() as td:
        print("temporary dir:", td)

        temp_file = op.join(td, "temp.mv5")
        f_out = h5py.File(temp_file, "w")
        
        # Store largest sizes of individual chromosomes.
        # Expect to be identical, but make sure.
        chromsizes = []
        for in_file in input_files:
            bw = pyBigWig.open(in_file)
            _chromsizes = bw.chroms() # dictionary
            if len(chromsizes) is 0:
                chromsizes = _chromsizes
            else:
                for (k, v) in _chromsizes.items():
                    if k not in chromsizes:
                        chromsizes[k] = v
                    elif k in chromsizes and chromsizes[k] < v:
                        chromsizes[k] = v
            bw.close()
        
        # Convert dict to a list of tuples to input to multivec function.
        chromsizes = sorted([(k, v) for k, v in chromsizes.items()], key=utils.sort_by_chrom)
        
        # TODO: Remove this line when tested with a single chromosome.
        if IS_DEBUG:
            # chromsizes = [("chr1", 248956422)]
            chromsizes=[("chr9", 138394717)]
            
        # Init a variable & create cache files with name and size
        raw_data = { chrom: np.full((len(input_files),size), EMPTY_VALUE) for (chrom, size) in chromsizes }
        for (chrom, size) in chromsizes:
            f_out.create_dataset(
                chrom,
                (
                    math.ceil(size / starting_resolution),
                    len(input_files),
                ),
                fillvalue=EMPTY_VALUE,
                compression="gzip",
            )
        
        print("Data initialized.")

        for file_index, cur_in_file in enumerate(input_files):
            bw = pyBigWig.open(cur_in_file)

            for (chrom, size) in chromsizes:
                
                # cur_position = 0
                
                for interval in bw.intervals(chrom):
                    (interval_start, interval_end, value) = interval
                    
                    # # Fill empty values when interval skips some positions
                    # if cur_position < interval_start:
                    #     raw_data[chrom][file_index][
                    #         cur_position : interval_start
                    #     ] = [EMPTY_VALUE] * (interval_start - cur_position)

                    # Fill suggested values in the interval
                    raw_data[chrom][file_index][
                        interval_start : interval_end
                    ] = [value] * (interval_end - interval_start)

                    # cur_position = interval_end

                # # Fill empty values when real positions are skipped
                # if cur_position is not size:
                #     raw_data[chrom][file_index][
                #         cur_position : size
                #     ] = [EMPTY_VALUE] * (size - cur_position)
                
                print(chrom, "in", cur_in_file, "processed.")

            bw.close()
            print("File", cur_in_file, "processed.")
        
        # Store raw data to a cache file
        for (chrom, size) in chromsizes:
            # TODO: When using batch, change index range.
            # f_out[chrom][0 : len(raw_data)] = np.array(raw_data)
            f_out[chrom][0 : size] = raw_data[chrom].T
        print("Cache file saved", f_out[chrom][0 : size].shape)
        
        f_out.close()

        tf = temp_file
        f_in = h5py.File(tf, "r")

        # The path of an output file.
        output_file = op.splitext([input_files[0]][0])[0] + ".multires.mv5"
        print("output_file:", output_file)

        # Override the output file if it existts
        if op.exists(output_file):
            os.remove(output_file)
        
        cmv.create_multivec_multires(
            f_in,
            chromsizes=chromsizes,
            agg=lambda x: np.nansum(x.T.reshape((x.shape[1], -1, 2)), axis=2).T, # Default aggregation lamda
            starting_resolution=starting_resolution,
            tile_size=256,
            output_file=output_file
        )

def main():
    # Test
    if False:
        output_file = op.splitext(["sample_data/debug2.bw"][0])[0] + ".multires.mv5"
        
        f = h5py.File('/tmp/blah.h5', 'w')
        f.create_dataset('chr1', (10000, 200), fillvalue=np.nan, compression='gzip')
        f['chr1'][0 : 10000] = np.random.random_sample((10000, 200))
        print("Cache file saved", f['chr1'][0 : 10000].shape)
        cmv.create_multivec_multires(
            f,
            chromsizes=[('chr1', 10000)],
            agg=lambda x: np.nansum(x.T.reshape((x.shape[1], -1, 2)), axis=2).T, # Default aggregation lamda
            starting_resolution=1,
            tile_size=256,
            output_file=output_file
        )
        f.close()
        return
    ######
    
    bigwig_to_multivec(test_input_files_58)

if __name__ == "__main__":
	main()