#!/usr/bin/env python
# Time-stamp: <2021-04-20 00:29:38 Tao Liu>

import os
import sys
import re
# ------------------------------------
# Main function
# ------------------------------------
def main():
    if len(sys.argv) < 2:
        sys.stderr.write("need 1+ paras: %s <xxx.peakstat.txt file>xN \n" % sys.argv[0])
        sys.exit(1)

    files = sys.argv[1:]
    print( "sample_name\treplicate\ttotal_reads\ttotal_peaks\tpeaks_in_blacklist\tpeaks_in_promoters\tratio_of_peaks_in_promoters\tpeaks_in_DHSs\tratio_of_peaks_in_DHSs")

    for filename in files:
        total_reads = -1
        total_peaks = -1
        blacklist_peaks = -1
        promoter_peaks = -1
        dhs_peaks = -1

        tmp = os.path.basename(filename).rstrip(".peakstat.txt")
        match = re.match("(.*?)\_r(\d+)",tmp)
        if match:
            (n, r) = match.groups()
        else:
            n = "All"
            r = "-"

        fhd = open( filename, "r" )

        # total reads
        l = fhd.readline()
        total_reads = int(l.split(": ")[1])
        # skip the next line
        fhd.readline()
        # total peaks
        l = fhd.readline()    
        total_peaks = int(l.rstrip())
        # skip the next three lines
        fhd.readline()
        fhd.readline()
        fhd.readline()
        # blacklist peaks
        l = fhd.readline()    
        blacklist_peaks = int(l.rstrip())
        # skip the next line
        fhd.readline()
        fhd.readline()
        fhd.readline()
        # in promoter
        l = fhd.readline()    
        promoter_peaks = int(l.rstrip())
        # skip the next line
        fhd.readline()        
        # in DHS
        l = fhd.readline()    
        dhs_peaks = int(l.rstrip())

        print( f"{n}\t{r}\t{total_reads}\t{total_peaks}\t{blacklist_peaks}\t{promoter_peaks}\t{promoter_peaks/total_peaks:.2f}\t{dhs_peaks}\t{dhs_peaks/total_peaks:2f}")
    
if __name__ == '__main__':
    main()
