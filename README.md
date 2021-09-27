# Reducing Reference Data for Metagenomic Classification
Master's thesis, with the goal of reducing reference data for taxonomic classification.

## Project Updates
- [Sep 27 - k-mer frequency counter](#sep_27)

##  <a name="sep_27"></a> K-mer frequency algorithm finished.
The basic algorithm is done. The steps include
- [x] Reading FASTA files (using klib)
- [x] Syncmer Algorithm 
- [x] Hashtable construction
- [x] Comparing sequences that landed in the same bucket.

#### Results
For the (eukaryotic) reference genome of Drosophila americana (NCBI Taxonomy ID 40366, ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/245/305/GCA_001245305.1_D._americana_W11_strain_genome_assembly) 
- the entire process took 272s (about 5 minutes)
  - 30s FASTA reading time, 
  - 30s Hashing time (this includes finding syncmers), 
  - 208s (about 3 min 30s) to compare sequences that landed in the same bucket.

#### Strengths
- A sequence that is assigned to a certain syncmer is referenced with 3 integers only (sequence number, subsequence position, syncmer position) and there is only one such subsequence per syncmer. This makes the memory overhead relatively small.
- In the last step, comparing sequences that landed in the same bucket, many tricks and shortcuts were used in order to try and minimize the number of operations. In particular checking whether two strings are equal by iterating over the characters is avoided if possible through the use of checksums.

#### Weaknesses
- The algorithm works with type char, even though possibilities only include a,c,g,t. More efficient storage would be bitvectors with two bits per char
- The hastable is still type vector<vector<int>>. This is certainly not optimal. It is questionable whether for this subproject it is helpful to get too complicated. It may be possible to replace it with a hashtable from klib.
- Sequence referencers are all 32 bit integers. 16 bit would certainly be enough, if not even 8 bit.
- The algorithm uses MurmurHash. Is this necessary? Is there something better?
- The algorithm does not allow "infinite" addition of new reference genomes, else it will get a memory overflow. The basic problem is that since the hashtable only stores numbers to find the sequence in the genome, the genome needs to stay in the memory.
