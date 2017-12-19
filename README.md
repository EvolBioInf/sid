# SNPs from diploid genomes

The `sid` program can be used to detect single nucleotide polymorphisms in aligned sequencing reads.
Based on the pileup output of (Samtools)[http://www.htslib.org/], for each genome site the most
likely genotype is computed.

## Installation

To install sid, clone the repository and change into the source directory. To prepare the cloned
repository for compilation, `cd` into the directory and run

```bash
./autogen.sh
./configure
```

The program is then compiled and (optionally) installed to `/usr/local/bin` with

```bash
make
sudo make install
```

## Usage

For preprocessing, a working installation of (`samtools`)[http://www.htslib.org/] is needed. Given
the aligned sequencing data in SAM or BAM format, a pileup needs to be generated:
```
samtools mpileup input.bam > pileup.dat
```
Useful parameters for the `samtools mpileup` program are `-C 50`, which is recommended for
BWA-aligned reads, and `-q 1` which discards reads with ambiguous mapping. The switches `-q` and
`-Q` can be used to set a lower threshold for mapping quality and base quality, respectively.

The pileup is then used as the `sid` input:
```
sid -m local pileup.dat > output.csv
```
The `-m` switch selects the SNP calling method, see `sid -h` for all available ones.

## Output format

The output is given as comma-separated values (CSV).
```
chrom,pos,label,gt,hom_conf,het_conf,conf_type
1,3003229,hom,GG,2.98235e-13,1,p_value
1,3003230,hom,GG,3.67401e-14,1,p_value
1,3003231,het,GA,1,8.01707e-17,p_value
1,3003232,hom,AA,4.36739e-12,1,p_value
```
Each line represents the genotype calling result for one genome site, defined by a chromosome
(`chrom`) and a coordinate (`pos`). The `label` is either `hom` or `het` and indicates a homozygous
or heterozygous site, respecitvely, while the `gt` column lists the called genotype. The next two
columns give the confidence in the call as a P-value.
