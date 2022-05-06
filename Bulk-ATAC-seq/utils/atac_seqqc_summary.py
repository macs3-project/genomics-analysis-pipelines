#!/usr/bin/env python
# Time-stamp: <2021-02-25 17:16:42 Tao Liu>

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
        # total
        while total_reads == -1:
            l = fhd.readline()
            if l.find(" in total (") != -1:
                total_reads = int(l.split(" ")[0])
        # dup
        while duplicated == -1:
            l = fhd.readline()
            if l.endswith("duplicates\n"):
                duplicated = int(l.split(" ")[0])
        # map
        while mapped == -1:
            l = fhd.readline()
            if l.find(" mapped (")  != -1:
                mapped = int(l.split(" ")[0])
        # "unique non chrM"
        while unique_non_chrM == -1:
            l = fhd.readline()
            if l.startswith("non chrM reads:"):
                l = fhd.readline()
                unique_non_chrM = int(l.split(" ")[0])
        # "unique non chrM promoter"
        while unique_non_chrM_promoter == -1:
            l = fhd.readline()
            if l.startswith("non chrM reads in promoter:"):
                l = fhd.readline()
                unique_non_chrM_promoter = int(l.split(" ")[0])
        # "unique non chrM peak"
        while unique_non_chrM_peak == -1:
            l = fhd.readline()
            if l.startswith("non chrM reads in peak:"):
                l = fhd.readline()
                unique_non_chrM_peak = int(l.split(" ")[0])

        print( f"{n}\t{r}\t{total_reads}\t{mapped}\t{mapped/total_reads:.2f}\t{duplicated}\t{duplicated/total_reads:.2f}\t{unique_non_chrM}\t{unique_non_chrM/total_reads:.2f}\t{unique_non_chrM_promoter}\t{unique_non_chrM_promoter/unique_non_chrM:.2f}\t{unique_non_chrM_peak}\t{unique_non_chrM_peak/unique_non_chrM:.2f}" )

if __name__ == '__main__':
    main()
