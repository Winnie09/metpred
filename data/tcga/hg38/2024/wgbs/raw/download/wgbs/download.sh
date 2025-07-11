#!/bin/bash -l
#SBATCH -A hji7_bigmem
#SBATCH --partition=bigmem
#SBATCH --time=23:58:59
#SBATCH --nodes=1
#SBATCH --ntasks=22
#SBATCH -o %x.o%j

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_BLCA_A13J.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_BLCA_A1AA.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_BLCA_A1AG.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_BLCA_A20V.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_BLCA_A2HQ.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_BLCA_A2LA.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_BLCA_NA20V.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_BRCA_A04X.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_BRCA_A07I.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_BRCA_A0CE.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_BRCA_A0YG.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_BRCA_A15H.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_BRCA_NA0CE.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_COAD_3518.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_COAD_A00R.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_COAD_N3518.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_GBM_0128.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_GBM_1401.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_GBM_1454.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_GBM_1460.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_GBM_1788.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_GBM_3477.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_LUAD_4630.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_LUAD_6148.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_LUAD_6215.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_LUAD_6840.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_LUAD_7156.bed.gz



wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_LUAD_N6148.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_LUSC_1078.bed.gz



wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_LUSC_2600.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_LUSC_2695.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_LUSC_2722.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_LUSC_N2722.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_READ_2689.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_READ_3593.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_READ_N2689.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_STAD_5730.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_STAD_6177.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_STAD_6452.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_STAD_6519.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_STAD_N6452.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_UCEC_A05J.bed.gz


wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_UCEC_A0G2.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_UCEC_A0K6.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_UCEC_A1CI.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_UCEC_A1CK.bed.gz

wget https://zwdzwd.s3.amazonaws.com/trackHubs/TCGA_WGBS/hg38/bed/TCGA_UCEC_NA1CI.bed.gz
