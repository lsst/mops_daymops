import os
import argparse
import numpy as np
import pandas as pd

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Swap ssmId's that are strings into integers.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputFile", type=str, help="Input file containing the DiaSources.")
    parser.add_argument("-t", "--trackId", type=str, default='ssmId_to_int.dat',
                        help="File tracking the map from the ssmId strings to integers.")
    parser.add_argument("-o", "--outFile", type=str, default=None,
                        help="Output file name for diasources with integer ssmIds. "\
                        "Default value is the input file, with addition of '_newId'.")
    parser.add_argument("-k", "--ssmidKey", type=str, default=None,
                        help="Name of the input file column containing SSMIDs."\
                        "If not specified, will use a list of default key names and look"\
                        "for the correct key.")
    args = parser.parse_args()

    diasources = pd.read_table(args.inputFile, delim_whitespace=True)

    # If ssmidKey is given on command-line then use that.
    # If it's not given, use list of defaults and find key.
    if args.ssmidKey is not None:
        ssmKey = args.ssmidKey
    else:
        ssmKeys = ['objId', 'ssmId', 'ObjID', 'objID', 'ssmID']
        for k in args.ssmidKeys:
            if k in diasources:
                ssmKey = k
    newIds = np.zeros(len(diasources), int)

    if args.outFile is None:
        if len(args.inputFile.split('.')) == 1:
            outFile = args.inputFile + '_newId'
        else:
            outFile = ('.'.join(args.inputFile.split('.')[:-1]) + '_newId'
                       + '.' + args.inputFile.split('.')[-1])
    else:
        outFile = args.outFile

    ssmIds = diasources[ssmKey].unique()

    for i, ssm in enumerate(ssmIds):
        condition = np.where(diasources[ssmKey] == ssm)[0]
        newIds[condition] = i

    # Write new diasource file.
    diasources[ssmKey] = newIds
    diasources.to_csv(outFile, sep=' ', index=False)

    # Write tracking file that maps ssmId to integers.
    with open(args.trackId, 'w') as trackFile:
        for i, ssm in enumerate(ssmIds):
            print >>trackFile, ssm, i
