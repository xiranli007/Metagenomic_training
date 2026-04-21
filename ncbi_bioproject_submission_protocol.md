# NCBI sequencing data submission protocol

## Important note

- **BioProject** = study-level container
- **BioSample** = metadata for each biological sample
- **SRA** = repository for the actual raw sequencing files

Your sequencing files are deposited to **SRA** and then linked to the relevant **BioProject** and **BioSample** records.

---

## Practical protocol for submitting sequencing data to NCBI

### 1. Prepare the three required components
Before starting the submission, prepare:

- **BioProject information**
  - project title
  - short project description

- **BioSample information**
  - one BioSample for each biologically distinct sample

- **SRA files and run metadata**
  - raw read files
  - library preparation details
  - sequencing platform information


---

### 2. Prepare sequence files in an accepted format
The safest formats for SRA submission are:

- **FASTQ**
- **BAM**

Important notes:

- **FASTA alone is usually not sufficient** because SRA expects per-base quality scores
- For paired-end FASTQ, submit either:
  - forward and reverse reads in **separate files in matching order**, or
  - an **interleaved 8-line FASTQ**
- Do **not** concatenate all forward reads followed by all reverse reads in a single file
- each **FASTQ** file should be **under 100 GB**
- files may be compressed with:
  - **gzip**
  - **bzip2**
- **zip is not accepted**
- studies larger than **5 TB** should be split across multiple submissions

Also avoid placing sensitive information in file names, since file names may become public.

---

### 3. Start the submission in the NCBI Submission Portal

- Go to the [**Submission portal**](https://submit.ncbi.nlm.nih.gov/), login or create an account.
- Go to **My submission** and choose **Biosample** and start a **New submission**
- Fill in the submission steps based on the instruction
- Your submission will receive a temporary **SUB#** identifier. Save that ID for tracking. Now you have you Bioproject created


---

### 4. Manage your data and create BioSamples (You can also create BioSamples on-to-go when submitting SRA, jump to step 5 if you prefer that way)
Create **one BioSample per biologically distinct sample**. Separate BioSamples are usually needed for differences such as:

- organism or metagenome
- time point
- tissue type
- treatment group

For each sample, prepare metadata such as:

- sample name
- organism or metagenome name
- collection date
- geographic location
- tissue or source
- treatment
- replicate information
- to link the BioSample to BioProject created earlier. Fill in the BioProject number you got to the correct column when filling in the data

Choose the correct **BioSample package** for your study based on the package description

---


### 5. Submit SRA and Upload the files

Smilar process as BioProject and BioSample New submission -> fill in the form


NCBI supports several upload options:

- web upload
- **Aspera Connect**
- **FTP / command line**
- cloud import from **Amazon S3**
- cloud import from **Google Cloud**

If you use a preload folder:

1. upload your files into your personal `uploads/...` folder
2. create a submission subfolder
3. wait about **10–15 minutes**
4. return to the wizard and choose **Select Submission Folder**

Example Aspera command pattern:

```bash
id="zdryu2"
ascp -i /path/to/key -QT -l 100m -k1 -d /path/to/folder \
subasp@upload.ncbi.nlm.nih.gov:uploads/your_email_xxxxx
```

Use the exact upload folder name assigned in your NCBI portal.

---

### 6. Final review and submit
On the **Overview** page, verify:

- release date
- BioProject linkage
- BioSample assignments
- run metadata
- file names

Then submit the record.

After submission, NCBI recommends allowing at least **24 hours** for processing before contacting support.

---

### 7. Special caution for human data
If the sequencing data are human and require controlled access, do **not** use the standard public SRA workflow.

Instead, submit through **dbGaP**.

---

## Simple checklist

- Prepare a spreadsheet with all sample metadata
- Confirm the correct BioSample package
- Format sequence files as FASTQ or BAM
- Verify paired-end files are matched correctly
- Keep each FASTQ under 100 GB
- Use gzip or bzip2 if compressing
- Start a new SRA submission
- Request a preload folder if needed
- Link or create BioProject and BioSamples
- Upload files
- Review all metadata carefully
- Submit

---

## Useful links

### Main submission guidance
- SRA submission overview: https://www.ncbi.nlm.nih.gov/sra/docs/submitbio
- SRA submission portal guide: https://www.ncbi.nlm.nih.gov/sra/docs/submitportal/

### BioSample package guidance
- BioSample package list: https://www.ncbi.nlm.nih.gov/biosample/docs/packages/

### Accepted file formats
- SRA accepted submission formats: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats

### File upload instructions
- SRA file upload instructions: https://www.ncbi.nlm.nih.gov/sra/docs/submitfiles/
- NCBI submission portal upload page: https://submit.ncbi.nlm.nih.gov/about/sra/

### Additional overview
- NCBI Insights article on sequence submission: https://ncbiinsights.ncbi.nlm.nih.gov/2023/05/01/sequences-genbank-sra/

---

## Suggested next step
If your dataset is one of the following:

- shotgun metagenomics
- 16S or ITS amplicon sequencing
- isolate whole-genome sequencing
- transcriptome sequencing

the protocol can be further customized to include the exact metadata fields you should prepare before submission.
