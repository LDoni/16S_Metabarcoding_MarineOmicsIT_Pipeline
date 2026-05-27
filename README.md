# 16S_Metabarcoding_MarineOmicsIT_Pipeline
---

This repository provides the bioinformatics scripts used to process and analyse the 16S rRNA gene metabarcoding data, developed within the framework of the Italian Omics Observatory Network of Marine Biodiversity (MarineOmicsIT). 
The workflow is designed to support the harmonized analysis of microbial communities collected from long-term marine observatories distributed along the Italian coasts (Portofino Promontory, Gulf of Trieste, Meda Senigallia, and MareChiara – Gulf of Naples), which are also part of the European Long-Term Ecological Research (eLTER) network.


<p align="center">
  <img src="img/Maps-Italian-Osservatory_traspartente-1-922x1024.png" width="400"/>
</p>

---
# 16S rRNA V4–V5 Amplicon Processing Pipeline

## Primers

Target region: **16S rRNA V4–V5**

| Primer | Sequence (5' → 3') |
|---|---|
| 515F-Y | `GTGYCAGCMGCCGCGGTAA` |
| 926R | `CCGYCAATTYMTTTRAGTTT` |



## Primer-based Read Filtering with BBDuk

Before running the downstream pipeline, reads are filtered to retain only pairs containing the expected primers.

Since primer trimming is later handled by DADA2 (`--trim-primers-dada2`), this step is used **only for filtering**, not trimming.

Because `bbduk.sh` works independently on forward and reverse reads, the filtering is performed in two sequential steps:

1. Filter R1 reads for the forward primer (`515F-Y`)
2. Filter R2 reads for the reverse primer (`926R`)



## Recommended Project Structure

```text
project/
├── data/
│   ├── MareChiara_net/
│   ├── MareChiara_water/
│   ├── GulfOfTrieste_run1/
│   ├── GulfOfTrieste_run2/
│   ├── PortofinoPromontory/
│   └── MedaSenigallia/
├── filtered_primers/
├── output_ampwrap/
└── scripts/
```

---

## Primer Filtering Script

```bash
#!/bin/bash

BASE_IN="./data"
BASE_OUT="./filtered_primers"

DIRS=(
"MareChiara_net"
"MareChiara_water"
"GulfOfTrieste_run1"
"GulfOfTrieste_run2"
"PortofinoPromontory"
"MedaSenigallia"
)

# Primer FASTA files
echo -e ">515F\nGTGYCAGCMGCCGCGGTAA" > primer_F.fa
echo -e ">926R\nCCGYCAATTYMTTTRAGTTT" > primer_R.fa

for DIR in "${DIRS[@]}"; do

    INDIR="${BASE_IN}/${DIR}"
    OUTDIR="${BASE_OUT}/${DIR}"

    mkdir -p "$OUTDIR"

    echo "Processing ${DIR}..."

    for R1 in "$INDIR"/*_R1.fastq.gz; do

        R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
        SAMPLE=$(basename "$R1" _R1.fastq.gz)

        echo "  Sample: ${SAMPLE}"

        # STEP 1 — Filter forward primer on R1
        bbduk.sh \
          in1="$R1" in2="$R2" \
          out1="$OUTDIR/tmp_R1.fastq.gz" \
          out2="$OUTDIR/tmp_R2.fastq.gz" \
          ref=primer_F.fa \
          k=19 mink=19 hdist=0 \
          restrictleft=19 rcomp=f ordered=t

        # STEP 2 — Filter reverse primer on R2
        bbduk.sh \
          in1="$OUTDIR/tmp_R1.fastq.gz" \
          in2="$OUTDIR/tmp_R2.fastq.gz" \
          out1="$OUTDIR/${SAMPLE}_R1.fastq.gz" \
          out2="$OUTDIR/${SAMPLE}_R2.fastq.gz" \
          ref=primer_R.fa \
          k=20 mink=20 hdist=0 \
          restrictleft=20 rcomp=f ordered=t

        rm "$OUTDIR/tmp_R1.fastq.gz" "$OUTDIR/tmp_R2.fastq.gz"

    done
done
```

---

## Evaluate Filtering Efficiency

The following script compares raw and filtered read counts using `seqkit`.

```bash
#!/bin/bash

BASE_IN="./data"
BASE_OUT="./filtered_primers"

DIRS=(
"MareChiara_net"
"MareChiara_water"
"GulfOfTrieste_run1"
"GulfOfTrieste_run2"
"PortofinoPromontory"
"MedaSenigallia"
)

echo -e "Dataset\tSample\tRaw_reads\tFiltered_reads\tPercent_kept"

for DIR in "${DIRS[@]}"; do

    INDIR="${BASE_IN}/${DIR}"
    OUTDIR="${BASE_OUT}/${DIR}"

    for R1 in "$INDIR"/*_R1.fastq.gz; do

        SAMPLE=$(basename "$R1" _R1.fastq.gz)

        RAW_R1="$R1"
        FILT_R1="$OUTDIR/${SAMPLE}_R1.fastq.gz"

        if [[ -f "$FILT_R1" ]]; then

            RAW_COUNT=$(seqkit stats -T "$RAW_R1" | awk 'NR==2{print $4}')
            FILT_COUNT=$(seqkit stats -T "$FILT_R1" | awk 'NR==2{print $4}')

            PERC=$(awk -v r="$RAW_COUNT" -v f="$FILT_COUNT" \
                'BEGIN{ if(r>0) printf "%.2f", (f/r)*100; else print 0 }')

            echo -e "${DIR}\t${SAMPLE}\t${RAW_COUNT}\t${FILT_COUNT}\t${PERC}%"

        fi
    done
done
```

---

## Generate General FASTQ Statistics

```bash
#!/bin/bash

DIRS=(
"MareChiara_net"
"MareChiara_water"
"GulfOfTrieste_run1"
"GulfOfTrieste_run2"
"PortofinoPromontory"
"MedaSenigallia"
)

touch seqkit_stats.txt

for DIR in "${DIRS[@]}"; do

    if [ -d "$DIR" ]; then

        echo "Processing directory: $DIR"

        seqkit stats --all "$DIR"/*fastq.gz --tabular >> seqkit_stats.txt

    fi
done
```

---

## Run Amplicon Processing with AmpWrap

```bash
conda run -n ampwrap ampwrap short \
    -i \
    ./filtered_primers/MareChiara_net \
    ./filtered_primers/MareChiara_water \
    ./filtered_primers/GulfOfTrieste_run1 \
    ./filtered_primers/GulfOfTrieste_run2 \
    ./filtered_primers/PortofinoPromontory \
    ./filtered_primers/MedaSenigallia \
    -a GTGYCAGCMGCCGCGGTAA \
    -A CCGYCAATTYMTTTRAGTTT \
    -l 372 \
    -d dada2_silva_genus138 \
    -o ./output_ampwrap \
    --trim-primers-dada2 \
    --bigdata \
    --resource-profile balanced \
    -c 8
```

---
## Notes

Raw sequencing data is available at:  
[NCBI SRA Project PRJNA1406700](https://www.ncbi.nlm.nih.gov/sra/PRJNA1406700)

Results are available at:  
[MarineOmics Results Portal](https://www.marineomics.it)

