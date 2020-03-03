# -*- coding: utf-8 -*-

import sys
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
test_input_files_58 = ["sample_output/" + str(i + 1) + "_cistrome_db.bw" for i in range(5)]

def bigwigs_to_multivec(
    input_files,
    output_file=None,
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
    
    utils.set_time()

    ## Handling Errors.
    # None input file.
    if len(input_files) is 0:
        print("No input file suggested.")
        return
    print(len(input_files), "files prepared.")
    
    # Not a bigwig file.
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
        raw_data = { chrom: np.zeros((len(input_files),size)) for (chrom, size) in chromsizes }
        for (chrom, size) in chromsizes:
            f_out.create_dataset(
                chrom,
                (
                    math.ceil(size / starting_resolution),
                    len(input_files),
                ),
                fillvalue=EMPTY_VALUE,
                compression="gzip"
            )
        
        print("Data initialized.", utils.get_time_duration())

        for file_index, cur_in_file in enumerate(input_files):
            bw = pyBigWig.open(cur_in_file)

            for (chrom, size) in chromsizes:
                
                for interval in bw.intervals(chrom):
                    (interval_start, interval_end, value) = interval
                    
                    raw_data[chrom][file_index][
                        interval_start : interval_end
                    ] = np.full((interval_end - interval_start), value)

            bw.close()
            print("File processed:", cur_in_file, utils.get_time_duration())
        
        # Store raw data to a cache file.
        # Use chunk for h5py to enhance performance.
        # https://stackoverflow.com/a/27713489
        chunck_size = 100000000
        for (chrom, size) in chromsizes:
            cur_size = 0
            while cur_size < size:
                next_size = cur_size + chunck_size if cur_size + chunck_size < size else size
                f_out[chrom][cur_size : next_size] = raw_data[chrom].T[cur_size : next_size]
                cur_size = next_size

        print("Cache file saved", f_out[chrom][0 : size].shape, utils.get_time_duration())
        f_out.close()

        tf = temp_file
        f_in = h5py.File(tf, "r")

        # The path of an output file.
        if output_file is None:
            output_file = op.splitext([input_files[0]][0])[0] + ".multires.mv5"
        print("output_file:", output_file, utils.get_time_duration())

        # Override the output file if it existts.
        if op.exists(output_file):
            os.remove(output_file)
        
        cmv.create_multivec_multires(
            f_in,
            chromsizes=chromsizes,
            agg=lambda x: np.nansum(x.T.reshape((x.shape[1], -1, 2)), axis=2).T, # Default aggregation lamda.
            starting_resolution=starting_resolution,
            tile_size=256,
            output_file=output_file
        )
        print("Done Converting.", utils.get_time_duration())

def convert():
    ## Test.
    # Simplest version of conversion.
    if False:
        output_file = op.splitext(["sample_data/debug2.bw"][0])[0] + ".multires.mv5"
        
        f = h5py.File('/tmp/blah.h5', 'w')
        f.create_dataset('chr1', (10000, 200), fillvalue=np.nan, compression='gzip')
        f['chr1'][0 : 10000] = np.random.random_sample((10000, 200))
        print("Cache file saved", f['chr1'][0 : 10000].shape)
        cmv.create_multivec_multires(
            f,
            chromsizes=[('chr1', 10000)],
            agg=lambda x: np.nansum(x.T.reshape((x.shape[1], -1, 2)), axis=2).T, # Default aggregation lamda.
            starting_resolution=1,
            tile_size=256,
            output_file=output_file
        )
        f.close()
        return
    ######
    
    # Check arguments.
    if len(sys.argv) <= 1:
        print("The path of an input file is not suggested.")
        return
    
    # Set output file path.
    output_file = sys.argv[2] if len(sys.argv) >= 3 else None

    # Get input files.
    try:
        input_path = open(sys.argv[1], "r")
    except IOError:
        print("There is no such file:", sys.argv[1])
        return

    input_files = []
    with input_path:
        for line in input_path.readlines():
            input_files += [line[:-1]]  # Remove newlines.

    bigwigs_to_multivec(
        input_files=input_files, 
        output_file=output_file
    )

if __name__ == "__main__":
	convert()