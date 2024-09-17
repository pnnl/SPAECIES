#!/usr/bin/env python

import sys

import numpy as np
import netCDF4 as nc4

filename1 = sys.argv[1]
filename2 = sys.argv[2]

file1 = nc4.Dataset(filename1, 'r')
file2 = nc4.Dataset(filename2, 'r')

files_are_same = True

for dim in file1.dimensions:
    if dim in file2.dimensions:
        size1 = len(file1.dimensions[dim])
        size2 = len(file2.dimensions[dim])
        if size1 != size2:
            files_are_same = False
            print("Dimension {} has mismatched size (file1={}, file2={})."
                  .format(dim, size1, size2))
    else:
        files_are_same = False
        print("Dimension only in first file: {}".format(dim))

for dim in file2.dimensions:
    if dim not in file1.dimensions:
        files_are_same = False
        print("Dimension only in second file: {}".format(dim))

for var in file1.variables:
    if var in file2.variables:
        shape1 = file1[var].shape
        shape2 = file2[var].shape
        if shape1 == shape2:
            if (file1[var][:] != file2[var][:]).any():
                files_are_same = False
                print("Mismatch in variable {} of shape {}.".format(var, shape1))
                diff = file2[var][:] - file1[var][:]
                num_diff = np.count_nonzero(diff)
                first_diff = np.array([np.nonzero(diff)[i][0] for i in range(len(shape1))])
                print("Number of differences: {}".format(num_diff))
                print("Mean difference: {}".format(diff.sum() / num_diff))
                print("Location of first difference: {} ({} vs. {})"
                      .format(first_diff, file1[var][tuple(first_diff)],
                              file2[var][tuple(first_diff)]))
        else:
            files_are_same = False
            print("Variable {} has inconsistent shape (file1={}, file2={})."
                  .format(var, shape1, shape2))
    else:
        files_are_same = False
        print("Variable only in first file: {}".format(var))

for var in file2.variables:
    if var not in file1.variables:
        files_are_same = False
        print("Variable only in second file: {}".format(var))

file1.close()
file2.close()

if files_are_same:
    err_code = 0
else:
    err_code = 1

sys.exit(err_code)
