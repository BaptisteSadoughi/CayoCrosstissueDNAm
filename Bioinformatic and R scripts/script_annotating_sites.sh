#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=10GB

# ======= Configuration =======
base_path="${1:-/your/default/path}" # <- define this path

BEDFILES_DIR="$base_path/MASH/bedfiles"
ANNOTATIONS_DIR="$base_path/annotations"
INPUT_BED="$BEDFILES_DIR/all_sites_age_CpG.bed"

# ======= Error checking =======

# Check if input BED file exists
if [[ ! -f "$INPUT_BED" ]]; then
  echo "Error: Input BED file not found: $INPUT_BED"
  echo "Please ensure the file exists or pass a valid BEDFILES_DIR as the first argument."
  exit 1
fi

# Check if BEDFILES_DIR is a directory
if [[ ! -d "$BEDFILES_DIR" ]]; then
  echo "Error: BEDFILES_DIR does not exist: $BEDFILES_DIR"
  exit 1
fi

# Check if ANNOTATIONS_DIR exists
if [[ ! -d "$ANNOTATIONS_DIR" ]]; then
  echo "Error: ANNOTATIONS_DIR does not exist: $ANNOTATIONS_DIR"
  echo "Please ensure the path is correct or pass it as the second argument."
  exit 1
fi

# Check if bedtools is available
if ! command -v bedtools &> /dev/null; then
  echo "Error: 'bedtools' is not installed or not in your PATH."
  exit 1
fi

echo "All paths and dependencies are valid. Running intersections..."

# ======= Genomic annotations =======
bedtools intersect -a "$INPUT_BED" \
-b "$ANNOTATIONS_DIR/rheMac10_refseqGene.bed" -wa > \
"$BEDFILES_DIR/all_sites_age_genebody.bed"

bedtools intersect -a "$INPUT_BED" \
-b "$ANNOTATIONS_DIR/rheMac10_refseqPromoter.bed" -wa > \
"$BEDFILES_DIR/all_sites_age_promoter.bed"

bedtools intersect -a "$INPUT_BED" \
-b "$ANNOTATIONS_DIR/rheMac10_CpGislands_masked.bed" -wa > \
"$BEDFILES_DIR/all_sites_age_CGI.bed"

# ======= ChromHMM annotations =======
TISSUES=(
    "adipose_E063"
    "adrenal_E080"
    "heart_E105"
    "kidney_E086"
    "liver_E066"
    "lung_E096"
    "skeletal_muscle_E108"
    "spleen_E113"
    "thymus_E112"
    "whole_blood_PBMC"
)

for TISSUE in "${TISSUES[@]}"; do
    CHROMHMM_FILE="$ANNOTATIONS_DIR/macaque_chromHMM/${TISSUE}_macaque_liftover.bed"
    OUTPUT_FILE="$BEDFILES_DIR/all_sites_age_${TISSUE}_chromHMM.bed"

    if [[ ! -f "$CHROMHMM_FILE" ]]; then
        echo "Skipping missing annotation file: $CHROMHMM_FILE"
        continue
    fi

    bedtools intersect -a "$INPUT_BED" -b "$CHROMHMM_FILE" -wo > "$OUTPUT_FILE"
    echo "Processed: $TISSUE"
done

echo "Done."
