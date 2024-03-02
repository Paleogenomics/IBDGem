#!/usr/bin/bash

printHelp()
{
   # Display help message
   echo "chrarm-stats.sh: Aggregate IBD2/IBD0 and IBD1/IBD0 LLRs over chromosome arms"
   echo "                 and generate plots of the resulting statistics."
   echo
   echo "Usage: ./chrarm-stats.sh -i summary-file-paths.txt -b hg19/hg38 -o out-prefix"
   echo "Options:"
   echo "-i    Tab-delimited file with 2 columns (no header): CHROM (chromosome name with 'chr' prefix)"
   echo "      and SUMMARY_FILE_PATH (path to associated summary file). One line per chromosome/summary file."
   echo "-b    Reference genome build used (hg19/hg38)."
   echo "-o    Output prefix for resulting statistics file and plots."
   echo "-h    Show help message."
   echo
}


# Retrieve options
while getopts ":hi:b:o:" option; do
   case $option in
      h) printHelp
         exit;;
      i) IN=$OPTARG;;
      b) BUILD=$OPTARG;;
      o) OUT=$OPTARG;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

echo "Calculating chromosome-arm LLRs..."
python chrarm-stats.py -s ${IN} -b ${BUILD} -o ${OUT}
echo "Done."
echo "Generating plots..."
python plot-chrarm-stats.py -i ${OUT}.txt -o ${OUT}
echo "Done."
