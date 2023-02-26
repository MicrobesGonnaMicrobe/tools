# Repository of useful commands

Just some commands, basic to a bit more specilised. Mainly for myself and archival purposes.

These have been collected through years through several helpful resources and some of my own trial and error. Thank you to all the smart folks on the StackOverflow and the like, I am forever grateful.

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