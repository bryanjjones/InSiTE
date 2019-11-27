#!/bin/bash
rm ./summary.csv
./scripts/TcB_ISA.py -q ./ctrl04.fastq -r --append_summary
./scripts/TcB_ISA.py -q ./ctrl05.fastq -r --append_summary
./scripts/TcB_ISA.py -q ./ctrl06.fastq -r --append_summary
./scripts/TcB_ISA.py -q ./ctrl07.fastq -r --barcode 'CTCTCTGATG' --append_summary
./scripts/TcB_ISA.py -q ./ctrl08.fastq -r --barcode 'CTCTCTGATG' --append_summary
./scripts/TcB_ISA.py -q ./ctrl09.fastq -r --barcode 'CTCTCTGATG' --append_summary
./scripts/TcB_ISA.py -q ./ctrl10.fastq -r --barcode 'TATATAGCAC' --append_summary
./scripts/TcB_ISA.py -q ./ctrl11.fastq -r --barcode 'TATATAGCAC' --append_summary
./scripts/TcB_ISA.py -q ./ctrl12.fastq -r --barcode 'TATATAGCAC' --append_summary
./scripts/TcB_ISA.py -q ./01.fastq --append_summary
./scripts/TcB_ISA.py -q ./02.fastq --append_summary
./scripts/TcB_ISA.py -q ./03.fastq --append_summary
./scripts/TcB_ISA.py -q ./04.fastq --append_summary
./scripts/TcB_ISA.py -q ./05.fastq --append_summary
./scripts/TcB_ISA.py -q ./06.fastq --append_summary
./scripts/TcB_ISA.py -q ./07.fastq --barcode 'CTCTCTGATG' --append_summary
./scripts/TcB_ISA.py -q ./08.fastq --barcode 'CTCTCTGATG' --append_summary
./scripts/TcB_ISA.py -q ./09.fastq --barcode 'CTCTCTGATG' --append_summary
./scripts/TcB_ISA.py -q ./10.fastq --barcode 'TATATAGCAC' --append_summary
./scripts/TcB_ISA.py -q ./11.fastq --barcode 'TATATAGCAC' --append_summary
./scripts/TcB_ISA.py -q ./12.fastq --barcode 'TATATAGCAC' --append_summary
./scripts/TcB_ISA.py -q ./13.fastq --barcode 'CTCACGATCT' --append_summary
./scripts/TcB_ISA.py -q ./14.fastq --barcode 'TACTCTGCGT' --append_summary

mv ./summary.csv ./gencode_summary.csv
mv ./*.csv ./results/
mv ./*.bam ./results/
mv ./*.fasta ./results/
mv ./*.log ./results/
mv ./*trimmed.fastq ./results/
mv ./*.sam ./results/
mv ./*.svg ./results/
