"""
Splits 'fullerDiaSource' format DiaSources by night and put them into separate files.
This is the first step in preparation for use by findTracklets on a per-night basis.

Also (optionally) writes out a per-obsHist file which holds all dias from a given image.
"""

import sys
import os
import argparse

def getNightNum(mjd, midnight):
    """Determine night number for any MJD."""
    night = int(mjd + 0.5 - midnight)
    return night


if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Split DiaSources in one input file into separate nights.")
    parser.add_argument("inputFile", type=str, help="Input file containing the diaSources.")
    nightly = "nightly"
    parser.add_argument("-n", "--nightlyDir", type=str, default=nightly,
                        help="Output directory containing diaSources split per night. Default is %s"\
                        % nightly)
    parser.add_argument("-o", "--obshistDir", type=str, default=None,
                        help="Output directory containing diaSources split per observation."\
                        " If not specified, diaSources split per observation are not created.")
    lsst_midnight = 0.166
    parser.add_argument("--midnight", type=float, default=lsst_midnight,
                        help="The average MJD value of midnight at the location of the observatory."\
                        " Default value (%.3f) is appropriate for LSST site." % lsst_midnight)
    args = parser.parse_args()

    # Create the output directories if they do not exist.
    if not os.path.isdir(args.nightlyDir):
        os.makedirs(args.nightlyDir)
    if args.obshistDir is not None:
        if not os.path.isdir(args.obshistDir):
            os.makedirs(args.obshistDir)

    # Read and write the input data.
    prev_night = None
    counter = 0
    with open(args.inputFile, 'r') as inFile:
        # Read diasources from input file.
        for line in inFile:
            values = line.split()
            diaId = values[0]
            # Skip a comment line
            try:
                diaId = int(diaId)
            except ValueError:
                # This was probably a comment, which would not translate to int.
                continue
            counter += 1

            # Determine the night number of this particular diasource and write to that file.
            mjd = float(values[5])
            nightNum = getNightNum(mjd, args.midnight)
            print mjd, nightNum

            # Open new output file if needed.
            if nightNum != prev_night:
                try:
                    outfile.close()
                except NameError:
                    # this was just the first night, so outfile doesn't exist yet.
                    pass
                outfile = open(os.path.join(args.nightlyDir, str(nightNum) + ".dias"), "aw")
                prev_night = nightNum
            # Write output line.
            # Since we're writing the whole thing back to disk, we do not need to convert other #'s.
            print>>outfile, line.rstrip()

            # Write to per-obsHist file if desired.
            if args.obshistDir is not None:
                outfile2 = open(os.path.join(args.obshistDir, str(obshistId) + ".dias"), "aw")
                print>>outfile2, line.rstrip()
                outfile2.close()

    print "Read %d lines from input file %s. Wrote to perNight files in %s." \
      % (counter, args.inputFile, args.nightlyDir)
    if args.obshistDir is not None:
        print "Also wrote perObsHist files in %s" % args.obshistDir
