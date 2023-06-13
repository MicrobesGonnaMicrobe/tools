# Repository of useful commands

Just some commands, basic to a bit more specialised. Mainly for myself and archival purposes.

These have been collected through years through several helpful resources and some of my own trial and error. Thank you to all the smart folks on the StackOverflow and the like, I am forever grateful.

## Content
- [Repository of useful commands](#repository-of-useful-commands)
  * [Resources to learn from, lists of tools](#resources-to-learn-from-lists-of-tools)
  * [Lists](#lists)
    + [Text file names (with extentions) list of all files](#text-file-names-with-extentions-list-of-all-files)
    + [Text file paths list of all files](#text-file-paths-list-of-all-files)
    + [text file: first name of file, then path of the file - text file list file name and path](#text-file-first-name-of-file-then-path-of-the-file-text-file-list-file-name-and-path)
    + [To obtain filenames and first lines for all files in the directory](#to-obtain-filenames-and-first-lines-for-all-files-in-the-directory)
    + [Comparing two unsorted lists, listing the unique in the second file](#comparing-two-unsorted-lists-listing-the-unique-in-the-second-file)
    + [Extract second column and turn to comma separated list](#extract-second-column-and-turn-to-comma-separated-list)
  * [Subsetting](#subsetting)
    + [Take out only one column](#take-out-only-one-column)
    + [Within the MAGS output folder](#within-the-mags-output-folder)
    + [Copy selected files between folders](#copy-selected-files-between-folders)
  * [Count](#count)
    + [to count how many > (sequences) are there, wc -l counts the lines](#to-count-how-many-sequences-are-there-wc-l-counts-the-lines)
  * [Sequence length of FASTA file](#sequence-length-of-fasta-file)
  * [grep](#grep)
    + [Print each line of a file which contains a match for pattern and keep the headers](#print-each-line-of-a-file-which-contains-a-match-for-pattern-and-keep-the-headers)
    + [Get the header lines of fasta sequence file](#get-the-header-lines-of-fasta-sequence-file)
  * [awk](#awk)
    + [How to merge two files based on column of one file:](#how-to-merge-two-files-based-on-column-of-one-file)
  * [Renaming file fasta header](#renaming-file-fasta-header)
    + [Rename headers/sequence IDs in multi-FASTA file](#rename-headerssequence-ids-in-multi-fasta-file)
  * [Single Line to Extract a Sequence from FASTA](#single-line-to-extract-a-sequence-from-fasta)
    + [Extracting more than one Sequence](#extracting-more-than-one-sequence)
  * [Split one fasta file into several files](#split-one-fasta-file-into-several-files)
  * [sed - replacing text in a file](#sed-replacing-text-in-a-file)
    + [Add filename to fasta headers in a loop: Introduce filename (genome name) in the fasta header in front of other things in the header](#add-filename-to-fasta-headers-in-a-loop-introduce-filename-genome-name-in-the-fasta-header-in-front-of-other-things-in-the-header)
    + [Shorten the fasta header names](#shorten-the-fasta-header-names)
    + [Remove header](#remove-header)
  * [To change tab file to fasta, fasta to tab](#to-change-tab-file-to-fasta-fasta-to-tab)
    + [tab2fasta](#tab2fasta)
    + [fasta2tab](#fasta2tab)
  * [Convert svg to png](#convert-svg-to-png)
  * [How to check if the files in two folders are identical](#how-to-check-if-the-files-in-two-folders-are-identical)
  * [Scripts from Windows to Linux: dos2unix](#scripts-from-windows-to-linux-dos2unix)
- [NCBI MAGs Upload](#ncbi-mags-upload)
  * [Formatting](#formatting)
  * [BioProject](#bioproject)
  * [Metagenomes](#metagenomes)
  * [MAGs](#mags)
    + [BioSample batch (MIMAG table)](#biosample-batch-mimag-table)
    + [Genome](#genome)
  * [SRA](#sra)

## Resources to learn from, lists of tools
* Awesome Bioinformatics - A compilation of Bioinformatic Resources by Daniel E Cook: https://github.com/danielecook/Awesome-Bioinformatics
* Rosalind Project - Learning Bioinformatics Through Problem Solving: https://rosalind.info/problems/list-view/
* Learning Bioinformatics for free at Open Source Society University (OSSU): https://github.com/ossu/bioinformatics
* EMBL-EBI on-demand training: https://www.ebi.ac.uk/training/on-demand
* EliXir e-learning materials: https://tess.elixir-europe.org/elearning_materials
* bio.tools by EliXir: https://bio.tools/
* Biotools - A compilation of Bioinformatic Resources by John Didion: https://github.com/jdidion/biotools

## Lists

### Text file names (with extentions) list of all files
```
ls -1 > MAGs_list.txt
```

### Text file paths list of all files
```
find "$PWD" > DPANN_highmedium_paths.txt 
```

### text file: first name of file, then path of the file - text file list file name and path
```
#!/bin/bash
echo -e "name\tcontigs_db_path" > namepathlist.txt
for bins in *.fa
do
basename=${​​​​bins%.fa}​​​​
id=${​​​​basename}​​​​
echo -e "${id}​​​​\t$PWD/${​​​​id}​​​​.fa"
done >> namepathlist.txt
```

### To obtain filenames and first lines for all files in the directory

```
awk '{print FILENAME" \"" $0"\""; nextfile}' *.fa > fastaheader.txt
```

### Comparing two unsorted lists, listing the unique in the second file
```
grep -Fxv -f MAGs_list.txt List2_GTDB_GCA.txt
```

### Extract second column and turn to comma separated list
```
awk '{print $2}' file.txt | paste -s -d, -
```

## Subsetting

### Take out only one column
Where -f2 extracts the 2, non-zero indexed column, or the second column:
```
cat textfile.tsv | cut -f2 -s > textfile_new.txt
```

This would return the fourth row: 
```
csvtool col 4 csv_file.csv
```

If you want to drop the header row:
```
csvtool col 4 csv_file.csv | sed '1d'
```

### Within the MAGS output folder
```
for i in pwd; do find ~+ -type f -name "*-contigs.fa"; done >> genome_paths
grep "MAG_" genome_paths >> MAGS_for_derep
```

### Copy selected files between folders
```
for file in `cat listoffiles`; do cp "$file" /path/of/destination ; done
```

## Count

### to count how many > (sequences) are there, wc -l counts the lines
```
cat filename.fasta | grep ">" | wc -l
grep -c ">" FILENAME
```

## Sequence length of FASTA file
```
awk '/^>/ {print; next; } { seqlen = length($0); print seqlen}' file.fa
```

## grep

### Print each line of a file which contains a match for pattern and keep the headers
```
for file in Annotation*; do [ "$file" = Annotation_cyc2.txt ] && continue  head -n 1 "$file"  grep -i 'cyc2' "$file"; done > Annotation_cyc2.txt
```

### Get the header lines of fasta sequence file
```
grep ">" FILENAME > HEADERFILE.txt
```

## awk

### How to merge two files based on column of one file:
```
awk -F'[\t|]' 'FNR==NR{a[$1]=$0;next}{$17=a[$3];print}' OFS='\t' BacMet2_PRE.155512.mapping.txt BS4_allcontigs_notbinned_prokka_BacMet_Pred.tsv > BS4_allcontigs_notbinned_prokka_BacMet_Pred_awk.tsv
```

`-F'[\t|]'` : delimiter is tab or pipe
`FNR==NR{a[$1]=$0;next}`: when file number equals awk line number (i.e. only for the first file) create an array a with column 1 as key and the whole line as value, then go to next line
`$17=a[$3]`: as the 17th column of the output, insert the value of the array (from file a) where the key equals column 3 in file b

by replacing `a[$1]=$0` by `a[$1]=$1"\t"$5` the value of the array would be column 1 and 5 separated by tab instead of the whole line.

To understand it, google "awk inner join"

## Renaming file fasta header

### Rename headers/sequence IDs in multi-FASTA file
```
python fa-rename.py --ids new_names.txt FASTA > new.fasta
```
  --ids FILE  specify two column tab-separated file with [oldnames] [newnames]


## Single Line to Extract a Sequence from FASTA

```
awk -v seq="contig_4299" -v RS='>' '$1 == seq {print RS $0}' unbinned_contigs.fasta > contig_4299_pilon.fa
```

### Extracting more than one Sequence
```
seqtk subseq {FASTA/FASTQ/FASTQ_GZ} {LIST_IDS} > {OUTPUT_FILE}
for i in *.ffn; do seqtk subseq $i grep_hydrogenases.txt > ${i%%.ffn}_hya.fa; done
```

## Split one fasta file into several files
```
awk '{if(/^>/){split($1,a,"[_]")}print >> a[2]".fa"}' BS4_checkv_HQvirus.fa
```

## sed - replacing text in a file
```
for i in *.fna; do sed -i 's/_/ontig/g' $i; done
```

### Add filename to fasta headers in a loop: Introduce filename (genome name) in the fasta header in front of other things in the header
```
for f in *.fna; do sed -i "s/^>/>${f%.fna}-/g" "${f}"; done
```

Header >c00001 of a file Faavne_M6_B18.fna is going to become >Faavne_M6_B18-c00001

### Shorten the fasta header names
"replace chromosome names in the format chr_I, chr_II, ... to I, II, ..."
```
sed -e 's/chr_I/I/' -e 's/chr_V/V/' -e 's/chr_X/X/' mySequence.fasta > mySeq.fasta
```
Or even simpler:
```
sed 's/^>chr_/>/' mySequence.fasta > mySeq.fasta
```

### Remove header
```
sed -i '1d' GS19_ROV14_R01_MAG_00004_sfz_information_cut_test.txt

for i in *_cut.txt; do sed -i '1d' $i; done
```

## To change tab file to fasta, fasta to tab

### tab2fasta
```
for genome in *; do python3 /export/work_cgb/Petra/tools/scripts/tab2fasta.py $genome ; done
```

### fasta2tab
```
python fasta2tab.py -i Zeta_16SrRNA_all_tab.fa -o Zeta_16SrRNA_all_tab.tab
```

## Convert svg to png
```
svg2png old.svg new.png
inkscape old.svg -b white --export-png=new.png
```

## How to check if the files in two folders are identical
Use find to list all the files in the directory then calculate the md5 hash for each file and pipe it sorted by filename to a file:
```
find /dir1/ -type f -exec md5sum {} + | sort -k 2 > dir1.txt
```

Do the same procedure to the another directory:
```
find /dir2/ -type f -exec md5sum {} + | sort -k 2 > dir2.txt
```

Then compare the result two files with diff:
```
diff -u dir1.txt dir2.txt
diff -u <(cat ne-shared.txt | cut -f1 -d" ") <(cat shared.txt | cut -f1 -d" ")
```

## Scripts from Windows to Linux: dos2unix
If you write your script in Windows Notes and copy it to server/linux, it will most likely not want to run, since not all characters are the same.
You need to convert your script so linux can read it, with dos2unix:

```
dos2unix ncbi_copy.sh
```

# NCBI MAGs Upload

Submitting MAGs to NCBI.

## Formatting
Ommit any non-ASCII characters (also "" and ~ and °).

## BioProject
"A BioProject is a collection of biological data related to a single initiative, originating from a single organization or from a consortium of coordinating organizations."
- Manually submit a BioProject
- You get a BioProject accession, like for example "PRJNA949439" 

## Metagenomes
- Submit BioSample: Download the template: Packages for metagenome submitters: https://submit.ncbi.nlm.nih.gov/biosample/template/ -> Metagenome or environmental, 1.0
- Or manually open up a "Submit BioSample" for each metagenome
- You will get a metagenome (sample) BioSample SAM accession

## MAGs
### BioSample batch (MIMAG table)
- Submit BioSample: Batch submission (maximum 1000 per submission): Download the template: Packages for metagenome submitters > GSC MIxS packages for genomes, metagenomes, and marker sequences > MIMAG Metagenome-assembled Genome > No environmental package
- Fill in info in the MIMAG table
- You can add the gtdb_taxonomy
- You need to get the appropriate NCBI taxonomy, GTDB taxonomy to NCBI taxonomy might help: https://github.com/nick-youngblut/gtdb_to_taxdump
- Insert the metagenome SAM accession under "derived-from" column into the MIMAG table
- Upload the MIMAG table (batch submission)
- You will get a BioSample MAG SAM accession

Warning messages
- "Warning: Provided taxonomy information was revised according to NCBI Taxonomy database rules. Please contact biosamplehelp@ncbi.nlm.nih.gov if you have any questions. "
- "Warning: Submission processing may be delayed due to necessary curator review. Please check spelling of organism, current information could not be resolved automatically and will require a taxonomy consult. For more information about providing a valid organism, including new species, metagenomes (microbiomes) and metagenome-assembled genomes, see https://www.ncbi.nlm.nih.gov/biosample/docs/organism/."

### Genome
- Info https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/
- Batch submission: Submit Genome > Batch/multiple genomes (maximum 400 per submission / max 300 files when uploading via web browser) > Upload the Genome Info table (GenomeBatch) - https://submit.ncbi.nlm.nih.gov/templates/
- Insert the BioSample MAG SAM accession under "biosample_accession"
- You can check your fasta files for contamination using NCBI Foreign Contamination Screen (FCS): https://github.com/ncbi/fcs

Answer questions on Gaps:
- Appropriate minimum number of Ns in a row (0-10) that represents a gap
- What type of evidence was used to assert linkage across the assembly gaps? Paired-ends.

Answer questions on Genome Info:
- Which type of genome? One or more chromosomes are in multiple pieces and/or some sequences are not assembled into chromosomes.
- Are these submissions exprected to represent the full genomes? Yes.
- Annotate this prokaryotic genome in the NCBI Prokaryotic Annotaton pipeline (PGAP) before its release? Yes. (Running PGAP yourself: https://github.com/ncbi/pgap (webinar: https://www.youtube.com/watch?v=pNn_-_46lpI))

Upload the fasta files: How do you want to provide files for this submission?
- FTP or Aspera Command Line file preload: All files for a submission must be uploaded into a single folder.
- Web browser upload via HTTP or Aspera Connect plugin: Do not use web browser HTTP upload if you are uploading files over 10 GB or more than 300 files. 


- Genome Submission: "Not for complete viral or organellar genomes. Submit those as regular GenBank records by emailing them to GenBank Submissions or using BankIt."

## SRA
"The SRA is a raw data archive, and requires per-base quality scores for all submitted data. SRA accepts binary files such as BAM, SFF, and HDF5 formats and text formats such as FASTQ."

- Submit Sequence Read Archive (SRA)
- Appropriate formats: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/
