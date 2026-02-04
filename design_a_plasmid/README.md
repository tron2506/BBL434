# Plasmid Design & Construction Pipeline

This repository contains a **Python-based bioinformatics pipeline** for the *de novo* identification of **Replication Origins (ORI)** and the **automated assembly of synthetic plasmids**.

---

## Overview

The pipeline follows a **multi-stage workflow** that transforms a raw bacterial genomic sequence (`Input.Fa`) into a **final circularizable plasmid sequence** (`output.txt`).

It integrates:
- Statistical genomic analysis (**GC skew**)
- Sequence motif enrichment (**k-mer analysis**)
- Modular biological part assembly
- Broad Host Range (BHR) plasmid backbone integration

---

## Features

### 1. ORI Identification (*In Silico*)

The Origin of Replication is detected using a **two-tier strategy**:

#### GC Skew Analysis
- Scans the genome using a sliding window
- Computes cumulative skew:  
  \[
  \frac{G - C}{G + C}
  \]
- Identifies the **global minimum**, corresponding to the replication switch point

#### K-mer Enrichment Analysis
- Focuses on the rough ORI candidate region
- Identifies **high-frequency 9-mers** (putative *DnaA boxes*)
- Compares regional vs global k-mer frequencies
- Extracts a **high-confidence 500 bp ORI sequence**

---

### 2. Modular Biological Part Retrieval

The pipeline dynamically assembles plasmid components from multiple sources:

#### Local Databases
- Searches `markers.fa` and `markers.tab`
- Retrieves antibiotic resistance markers  
  *(e.g., Ampicillin, Kanamycin)*

#### NCBI Entrez Integration
- Automatically queries the **NCBI Nucleotide database**
- Downloads missing genetic parts  
  *(e.g., `lacZ_alpha`)*

#### Restriction Site Mapping
- Converts enzyme names from a CSV file
- Maps enzymes (e.g., EcoRI, BamHI) to recognition sequences

---

### 3. Broad Host Range (BHR) Integration

- Incorporates **IncQ-type replicon logic** (i used RSF1010 E.Coli vector from NCBI database)
- Loads complete **Broad Host Range plasmid backbones**
- Example backbone: `AB526842.1`
- Ensures replication across **diverse bacterial species**

---

## File Requirements

Ensure the following files are present in the working directory:

| File | Description |
|-----|-------------|
| `Input.Fa` | Raw bacterial genome (FASTA) for ORI detection |
| `Design.txt` | Comma-separated list of design components |
| `restriction_enzymes.csv` | Restriction enzyme → recognition site mapping |
| `markers.fa` / `markers.tab` | Local antibiotic marker database |
| `BHR_default.fa` | Default Broad Host Range plasmid backbone |

---

## Pipeline Logic Flow

1. **Calculate GC Skew**  
   → Identify the *rough ORI center*

2. **Refine ORI**  
   → Apply 9-mer enrichment to extract a 500 bp ORI

3. **Construct Part 1**  
   → Parse `Design.txt`  
   → Append restriction sites and antibiotic markers

4. **Append Backbone**  
   → Add Broad Host Range plasmid sequence

5. **Finalize Output**  
   → Export full plasmid sequence to `output.txt`

---

## Summary of Results

Based on the current execution:

- **Genomic Origin Identified At:** ~1036 bp  
- **Total Plasmid Length:** 10,729 bp  
- **Restriction Sites Added:**  
  `BamHI`, `HindIII`, `PstI`, `SphI`, `SalI`, `XbaI`, `KpnI`, `SacI`, `SmaI`

---

## Notes

- Designed for **automation-friendly plasmid engineering**
- Suitable for **synthetic biology**, **genome analysis**, and **BHR construct design**
- Easily extendable with additional biological parts or replicon types

---
