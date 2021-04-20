#!/usr/bin/env python
# Time-stamp: <2021-04-20 00:28:56 Tao Liu>

import os
import sys
import re
# ------------------------------------
# Main function
# ------------------------------------
def main():
    if len(sys.argv) < 2:
        sys.stderr.write("need 1+ paras: %s <xxx.stat.txt file>xN \n" % sys.argv[0])
        sys.exit(1)

    files = sys.argv[1:]

    print( "sample_name\treplicate\ttotal_reads\tmapped_reads\tratio_of_total\tduplicated_reads\tratio_of_total\tunique_Q30_non-chrM_reads\tratio_of_total\tunique_Q30_non-chrM_in_promoter\tratio_of_valid_reads\tunique_Q30_non-chrM_in_peaks\tratio_of_valid_reads" )
    
    for filename in files:
        total_reads = -1
        mapped = -1
        duplicated = -1
        unique_non_chrM = -1
        unique_non_chrM_promoter = -1
        unique_non_chrM_peak = -1

        tmp = os.path.basename(filename).rstrip(".stat.txt")
        (n, r) = re.match("(.*?)\_r(\d+)",tmp).groups()
        fhd = open( filename, "r" )
        # skip first line "flagstat:"
        fhd.readline()
        # total
        l = fhd.readline()
        total_reads = int(l.split(" ")[0])
        # skip the next two lines
        fhd.readline()
        fhd.readline()
        # dup
        l = fhd.readline()    
        duplicated = int(l.split(" ")[0])
        # mapped
        l = fhd.readline()    
        mapped = int(l.split(" ")[0])
        # read until 'mapped Q30'
        while unique_non_chrM == -1:
            l = fhd.readline()
            if l.startswith("non chrM reads:"):
                l = fhd.readline()
                unique_non_chrM = int(l.split(" ")[0])
        # skip the next line
        fhd.readline()
        # in promoter
        l = fhd.readline()    
        unique_non_chrM_promoter = int(l.split(" ")[0])
        # skip the next line
        fhd.readline()
        # in peak
        l = fhd.readline()    
        unique_non_chrM_peak = int(l.split(" ")[0])

        print( f"{n}\t{r}\t{total_reads}\t{mapped}\t{mapped/total_reads:.2f}\t{duplicated}\t{duplicated/total_reads:.2f}\t{unique_non_chrM}\t{unique_non_chrM/total_reads:.2f}\t{unique_non_chrM_promoter}\t{unique_non_chrM_promoter/unique_non_chrM:.2f}\t{unique_non_chrM_peak}\t{unique_non_chrM_peak/unique_non_chrM:.2f}" )

if __name__ == '__main__':
    main()
