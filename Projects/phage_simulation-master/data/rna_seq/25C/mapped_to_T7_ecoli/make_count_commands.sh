for file in ../../mapped/mapped_T7_ecoli/*T7W*_sort.bam*; do echo "bedtools multicov -split -bams $file -bed ../../reference/ecoli_T7_ref.gff > ./$(basename ${file%.*})_counts.tsv"; done > count_commands.sh;

