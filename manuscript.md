---
title: Protein k-mer analyses for assembly- and alignment-free sequence analysis
keywords:
- metagenomics
- AAI
- Alignment-free
- Assembly-free
- MinHash
- FracMinHash
- Containment
lang: en-US
date-meta: '2021-12-26'
author-meta:
- N. Tessa Pierce-Ward
- C. Titus Brown
header-includes: |-
  <!--
  Manubot generated metadata rendered from header-includes-template.html.
  Suggest improvements at https://github.com/manubot/manubot/blob/main/manubot/process/header-includes-template.html
  -->
  <meta name="dc.format" content="text/html" />
  <meta name="dc.title" content="Protein k-mer analyses for assembly- and alignment-free sequence analysis" />
  <meta name="citation_title" content="Protein k-mer analyses for assembly- and alignment-free sequence analysis" />
  <meta property="og:title" content="Protein k-mer analyses for assembly- and alignment-free sequence analysis" />
  <meta property="twitter:title" content="Protein k-mer analyses for assembly- and alignment-free sequence analysis" />
  <meta name="dc.date" content="2021-12-26" />
  <meta name="citation_publication_date" content="2021-12-26" />
  <meta name="dc.language" content="en-US" />
  <meta name="citation_language" content="en-US" />
  <meta name="dc.relation.ispartof" content="Manubot" />
  <meta name="dc.publisher" content="Manubot" />
  <meta name="citation_journal_title" content="Manubot" />
  <meta name="citation_technical_report_institution" content="Manubot" />
  <meta name="citation_author" content="N. Tessa Pierce-Ward" />
  <meta name="citation_author_institution" content="Department of Population Health and Reproduction, University of California, Davis" />
  <meta name="citation_author_orcid" content="0000-0002-2942-5331" />
  <meta name="twitter:creator" content="@saltyscientist" />
  <meta name="citation_author" content="C. Titus Brown" />
  <meta name="citation_author_institution" content="Department of Population Health and Reproduction, University of California, Davis" />
  <meta name="citation_author_orcid" content="0000-0001-6001-2677" />
  <meta name="twitter:creator" content="@ctitusbrown" />
  <link rel="canonical" href="https://bluegenes.github.io/2021-paper-protein-kmers/" />
  <meta property="og:url" content="https://bluegenes.github.io/2021-paper-protein-kmers/" />
  <meta property="twitter:url" content="https://bluegenes.github.io/2021-paper-protein-kmers/" />
  <meta name="citation_fulltext_html_url" content="https://bluegenes.github.io/2021-paper-protein-kmers/" />
  <meta name="citation_pdf_url" content="https://bluegenes.github.io/2021-paper-protein-kmers/manuscript.pdf" />
  <link rel="alternate" type="application/pdf" href="https://bluegenes.github.io/2021-paper-protein-kmers/manuscript.pdf" />
  <link rel="alternate" type="text/html" href="https://bluegenes.github.io/2021-paper-protein-kmers/v/51e22610144d9a9c1d74e5a27d0915d6089d309f/" />
  <meta name="manubot_html_url_versioned" content="https://bluegenes.github.io/2021-paper-protein-kmers/v/51e22610144d9a9c1d74e5a27d0915d6089d309f/" />
  <meta name="manubot_pdf_url_versioned" content="https://bluegenes.github.io/2021-paper-protein-kmers/v/51e22610144d9a9c1d74e5a27d0915d6089d309f/manuscript.pdf" />
  <meta property="og:type" content="article" />
  <meta property="twitter:card" content="summary_large_image" />
  <link rel="icon" type="image/png" sizes="192x192" href="https://manubot.org/favicon-192x192.png" />
  <link rel="mask-icon" href="https://manubot.org/safari-pinned-tab.svg" color="#ad1457" />
  <meta name="theme-color" content="#ad1457" />
  <!-- end Manubot generated metadata -->
bibliography:
- content/manual-references.json
manubot-output-bibliography: output/references.json
manubot-output-citekeys: output/citations.tsv
manubot-requests-cache-path: ci/cache/requests-cache
manubot-clear-requests-cache: false
...






<small><em>
This manuscript
([permalink](https://bluegenes.github.io/2021-paper-protein-kmers/v/51e22610144d9a9c1d74e5a27d0915d6089d309f/))
was automatically generated
from [bluegenes/2021-paper-protein-kmers@51e2261](https://github.com/bluegenes/2021-paper-protein-kmers/tree/51e22610144d9a9c1d74e5a27d0915d6089d309f)
on December 26, 2021.
</em></small>

## Authors



+ **N. Tessa Pierce-Ward**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon width=16 height=16}
    [0000-0002-2942-5331](https://orcid.org/0000-0002-2942-5331)
    · ![GitHub icon](images/github.svg){.inline_icon width=16 height=16}
    [bluegenes](https://github.com/bluegenes)
    · ![Twitter icon](images/twitter.svg){.inline_icon width=16 height=16}
    [saltyscientist](https://twitter.com/saltyscientist)<br>
  <small>
     Department of Population Health and Reproduction, University of California, Davis
     · Funded by NSF 1711984, NSF 2018911
  </small>

+ **C. Titus Brown**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon width=16 height=16}
    [0000-0001-6001-2677](https://orcid.org/0000-0001-6001-2677)
    · ![GitHub icon](images/github.svg){.inline_icon width=16 height=16}
    [ctb](https://github.com/ctb)
    · ![Twitter icon](images/twitter.svg){.inline_icon width=16 height=16}
    [ctitusbrown](https://twitter.com/ctitusbrown)<br>
  <small>
     Department of Population Health and Reproduction, University of California, Davis
     · Funded by Moore Foundation GBMF4551
  </small>



## Abstract {.page_break_before}




## Background

As the scale of genomic sequencing continues to grow, alignment-free methods for estimating sequence similarity have become critical for conducting tasks ranging from taxonomic classification to phylogenetic analysis on large-scale datasets [@doi:10.1186/s13059-016-0997-x; @doi:10.1186/gb-2014-15-3-r46].
The majority of alignment-free methods rely upon exact matching of k-mers: subsequences of length k, that can be counted and compared across datasets, with or without use of subsampling methods such as MinHash.
As k-mer based methods rely on exact sequence matches, they can suffer from limited sensitivity when comparing highly polymorphic sequences or classifying organisms from groups that are not well represented in reference databases.

Current best practices methods can still only categorize a fraction of the metagenomic and metatranscriptomic data, especially for understudied and/or diverse habitats (xx% recovery for soil, xx% recovery ocean metagenomes, etc).
Even well-studied environments such as human gut can produce significant uncharacterized metagenome content.
"For example, a reference-based approach failed to map 35% of reads in the iHMP study on inflammatory bowel disease (Supp. Data. of (Franzosa et al., 2019)), omitting them from any further analysis. These reads may belong to unknown microbes, phage or viruses, plasmids, or accessory elements of known microbes, all of which can
play a role in disease.[from RO1]". This phenomenon is not restricted to metagenome samples. Alignment-based estimates can fail at larger evolutionary distances and even rRNA amplicon surveys may underestimate bacterial diversity [@doi:10.1128/AEM.00014-18].

To increase sensitivity of alignment-free methods, modified k-mer approaches have been introduced, including spaced seeds /split k-mers, which accommodate polymorphic sites in highly similar genomes (CITE).
For larger evolutionary distances, protein-based comparisons have long been the gold-standard approach for taxonomic and functional annotation, as protein sequence is more conserved than the underlying DNA sequence [@pubmed:2231712; @doi:10.1038/nmeth.3176].
As microbial and viral genomes are gene-dense, [MinHash-based] alignment-free comparisons of translated protein sequence have been shown to increase sensitivity for taxonomic classification and genome discovery [@doi:10.1038/ncomms11257; @doi:10.1186/s13059-019-1841-x].
Here, we demonstrate the utility of protein k-mer comparisons for phylogenomic reconstruction and taxonomic classification at larger evolutionary distances and across both gene-rich and [gene-sparse] sequences.
We use Scaled Minhash subsampling to facilitate conducting these comparisons at scale [Irber et al., 2021; @https://dib-lab.github.io/2020-paper-sourmash-gather/].

Scaled Minhash is a MinHash variant for selecting and hashing a set of representative k-mers from a sequence dataset [@https://dib-lab.github.io/2020-paper-sourmash-gather/]. Unlike traditional MinHash, Scaled MinHash sketches scale with the size of the dataset, meaning each sketch is comprised of the chosen proportion of k-mers in the input dataset, rather than a chosen number of k-mers.
Downsampling sequencing datasets in this way enables estimation of containment, which has been shown to permit more accurate estimation of genomic distance, particularly for genomes of very different lengths [@doi:10.1016/j.amc.2019.02.018; @doi:10.1186/s13059-019-1875-0].
Streaming containment estimates have been shown to facilitate genome discovery and correlate with Mash Distance, a proxy for Average Nucleotide Identity (ANI) [@doi:10.1186/s13059-019-1841-x; @doi:10.1186/s13059-020-02159-0].

Standardized genomic measures of relatedness such as ANI and its protein counterpart, Average Amino Acid Identity (AAI) have shown lasting utility for genome relatedness and phylogenomic analysis.
Traditional ANI and AAI describe the sequence similarity of all orthologous genes, either in nucleotide or protein space, respectively.
Both been shown to be robust measure of overall pairwise genome relatedness even for highly incomplete datasets, such as those comprised of only ~4% of the genome or 100 genes [@doi:10.1128/AEM.01398-06; @doi:10.1038/ismej.2017.113].
ANI has emerged as the most widely-accepted method for estimating pairwise similarity of microbial genomes and delimiting species boundaries [@doi:10.1073/pnas.0906412106].
Recent research appears to confirm 95% ANI species threshold for prokaryotic species, although there is some debate as to the universality of this threshold [@doi:10.1038/s41467-018-07641-9; @doi:10.1128/mSystems.00731-19; @doi:10.1101/2020.07.27.223511].
AAI thresholds have been proposed for higher taxonomic ranks, <45%, 45-65% and 65-95% for family, genus, and species [@doi:10.1016/j.mib.2007.08.006; @doi:10.1038/ismej.2017.113].
While traditional alignment-based estimation of ANI and AAI are computationally intensive, sketching-based estimates and sketching-facilitated estimates have permitted ANI calculations at the scale of whole-databases [@doi:10.1186/s13059-016-0997-x; @doi:10.1186/s13059-019-1841-x; @doi:10.1038/s41467-018-07641-9].

[Pierce-Ward et al., 2021 (tbd technical paper)] showed that Scaled MinHash containment estimates can be used to approximate both ANI (nucleotide k-mers) and Average Amino Acid Identity (AAI; protein k-mers), while accounting for the non-independence of mutated k-mers [@doi:10.1101/2021.01.15.426881].
Furthermore, Scaled MinHash containment estimates work well for genome pairs of varying lengths and for compositional analysis of metagenome samples.
Taken together, these properties enable robust assembly and alignment-free pairwise relatedness estimation that can be used on sequences separated by a wide range of evolutionary distances.
Here, we demonstrate that the utility of Scaled MinHash protein containment, both used directly and a an approximation of ANI and AAI, for taxonomic classification and phylogenomic reconstruction for species across the tree of life.


#### Notes

- AAI::phylogeny https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1236649/








## Results

K-mer analysis methods enable similarity detection as low as a single shared k-mer between divergent genomes.
As a result, exact matching long nucleotide k-mers can be used for taxonomic classification between closely related genomes, including at the strain, species, and genus level (k-mer lengths 51, 31, and 21, respectively).
At larger evolutionary distances, accumulated nucleotide divergence limits the utility of exact nucleotide k-mer matching.
As protein sequences are more conserved than their corresponding nucleotide sequences, exact matching of protein k-mers can enable similarity analysis across larger evolutionary distances.
To start, we choose protein k-mer sizes that correspond to known informative DNA ksizes, 21 (protein k-mer length: 7) and 31 (closest protein k-mer length: 10).

<!---
NOTE: worth looking at k=17 /51???.
--->

### Protein k-mers enable alignment-free comparisons at increased evolutionary distances

<!---
** gtdb evolpaths (all k-mers)**
As the number of genomes per taxon varies widely across GTDB, comparisons across the entire GTDB database may be impacted by the particular genomes included in the database.
To minimize any database biases, we developed a randomized subset of the GTDB database to assess the utility of protein k-mers across a wide range of evolutionary distances.

[GTDB] This database provides a wide array of genomes for assessing the utility of protein k-mers for bacterial and archaeal similarity estimation and taxonomic classification.

_OR SPECIES VERSION:
This "evolutionary path" consists of eight genomes: an anchor genome, a genome matching anchor taxonomy down to the species level, one matching anchor taxonomy to the genus level, one matching to the family level, and so on.
This creates a gradient of similarity, where comparisons to the anchor genome range from species-level to superkingdom-level._


Path selection using the representative genomes in **GTDB release 95** resulted in 2957 paths comprised of 6690 unique genomes (6543 Bacteria, 237 Archaea).
These paths include genome comparisons across 33 phyla (29 Bacteria, 4 Archaea), covering roughly a quarter of the 129 phyla (111 Bacteria, 18 Archaea) in GTDB release 95.
While paths are limited to taxonomies with at least two GTDB representative genomes for each taxonomic rank, these paths provide a rich resource for comparisons at increasing evolutionary distances. 
--->

The Genome Taxonomy Database (GTDB) provides a genome-based taxonomy for bacterial and archaeal genomes [@doi:10.1038/s41587-020-0501-8].
To assess the utility of protein k-mers for comparisons at an array of evolutionary distances, we selected a subset of these genomes that would allow standardized comparisons across taxonomic ranks.
For each genus with at least two species clusters in GTDB, one representative genome was randomly selected as an "anchor" genome.
Then, one additional genome was selected from the GTDB representative genomes matching the anchor's taxonomy at each higher taxonomic rank.
This "evolutionary path" consists of seven genomes: an anchor genome, a genome matching anchor taxonomy down to the genus level, one matching anchor taxonomy to the family level, one matching to the order level, and so on.
This creates a gradient of similarity, where comparisons to the anchor genome range from genus-level to superkingdom-level.

Path selection using the representative genomes in GTDB rs202 resulted in 4095 paths comprised of 9213 unique genomes (8790 Bacteria, 333 Archaea).
_These paths include genome comparisons across 33 phyla (29 Bacteria, 4 Archaea), covering roughly a quarter of the 129 phyla (111 Bacteria, 18 Archaea) in GTDB release 95._
While paths are limited to taxonomies with at least two GTDB representative genomes for each taxonomic rank, these paths provide a rich resource for comparisons at increasing evolutionary distances. 

**_to do: gtdb rs202 version_**

**WHY no multi-species representatives?**
![**More protein k-mers are shared at genus level** CAPTION](images/pseudomonas_jaccard_vs_containment_prot10.png)

![**Protein k-mers are shared at higher taxonomic ranks** CAPTION](images/anchor-containment.nucl-prot.png)

**genomes in same genus, fraction of k-mers in common at each ksize? ALL THE K-MERS**

### Accurate distance estimation from k-mer containment

Jaccard and Containment of DNA k-mers can be transformed into an estimate of the Average Nucleotide identity between genomes [cite Ondov Mash, Koslicki k-mer paper, koslicki scaled mh paper]. To begin with, we can apply the same equations to protein k-mer comparison statistics to obtain an estimate of Amino Acid Identity.

[Koslicki jaccard k-mer stats paper] showed how to properly transform Jaccard --> ANI assuming a simple mutational model. For protein k-mers, Jaccard/Containment can be transformed into Amino Acid Identity (AAI).

We can apply this same equation to protein k-mers 

**__to do: A. redo with MRCC estimate; use only best graph (or maybe two ksizes)__**
![**Scaled MinHash AAI vs CompareM**
GTDB Evolpaths dataset](images/gtdb95-evolpaths.AAI-concordance.png){#fig:evolpathsAAIvsCompareM height=2in}

** talk with david: max containment? avg conainment?


### **Anchor containment** enables comparisons directly from DNA sequence

In most cases, novel queries are DNA sequences

For protein k-mer comparisons to be useful, DNA queries must be translated into protein sequence, either via assembly and translation or direct 6-frame translation. While assembly-based methods are more accurate, they often only work for a fraction of the available data. In contrast, uninformated 6-frame translation can utilize the entire dataset, but generates a larger set of k-mers, only 1/6th of which are true protein k-mers. Here, the containment index is especially useful: by using only the containment estimate relative to the trusted reference proteomes, we can obtain accurate Amino Acid Identity estimates directly from DNA sequence. We term this "anchor" containment, where the trusted genome is the "anchor" upon which we base the comparison.




** figure: AAI from translated nucleotide --> reference protein**



[reproduce equation? Scaled MinHash Koslicki]
** all k-mers; scaled minhash version **
**[to do: run evolpaths with SCALED!]**



### Robust Taxonomic classification from protein k-mer containment

Anchor containment can also be used to enable robust taxonomic classification from either protein or 6-frame translated DNA queries. 



* talk about full gtdb databases here?
* not sure if need them above so can do the rankinfo assessment?


### Metagenome breakdown using protein k-mers

Anchor containment also enables metagenome breakdown via min-set-cov.

Metagenome sequences can also be translated in 6-frames, and then the anchor containment can be assessed relative to proteomes in a reference database.

*** use a mock metagenome, then a evolutionarily distant metagenome. compare the % of genome recovred with DNA, protein at diff ksizes. Genome grist it, basically.


THIS USES BOTH THE 6-FRAME translation and sourmash tax!!! containment improtant, etc. I think this is the way to go.



_not exactly sure where to put this section_

### GTDB databases

To enable similarity estimation and taxonomic classification.

The most recent GTDB release, `rs202`, encompasses 258,407 genomes from 47,895 species. 
To make analyses at this scale tractable, we downsample k-mers using `sourmash` Scaled MinHash sketches, with a scaling factor of 1000 for nucleotide k-mers (keep ~1/1000 k-mers) and 100 for protein k-mers (keep ~1/100 protein k-mers).
For most genomes, both genomic and protein fastas were available for download from NCBI. In remaining cases (n=36,632), genome fastas were translated into protein sequence via PRODIGAL (CITE) prior to sketching. 

We built `sourmash` databases for GTDB, `rs202` version, both for the species-representatives subset, and for all genomes.
These are available as part of the `Prepared Databases` section of the `sourmash` documentation, and archived on OSF [https://osf.io/t3fqa/] /Zenodo???.

### Different k-mer lengths provide resolution at different taxonomic ranks

`sourmash lca rankinfo`.

For each database and k-mer type, we first assessed the number of k-mers specific to each taxonomic rank: k-mers specific to a geonome were only present in one genome in the entire database, k-mers specific to a species were found in at least two genomes of the same species, etc. K-mers specific to a "superkingdom" were found in genomes spanning at least two phyla. While these characterizations are greatly impacted by the genomes included in the database, the differences observed between k-mer types and sizes suggests that different k-mer sizes may provide resolution at different taxonomic ranks. 

<!---
note 31, 51 --- maybe partially a result of database issues, e.g. not all species have multiple members; sometimes all members are closely related.
--->

For all DNA k-mer sizes, the majority of k-mers are present in only a single species, with only a few k-mers shared across genera.
Long nucleotide k-mers have already been shown to be useful for comparing genomes within the same genus or species.
Only at a dna k-mer size of 21 are a significant fraction of k-mers present in genomes shared across different families or even phyla.
In contrast, all protein k-mer sizes contain a portion of k-mers that are shared across genera and above.
At a protein k-mer size of 7, over 80% of k-mers are present in genomes found in more than one phylum, while at a protein k-size of 10, the number of genome-specific k-mers is more similar to that observed for nucleotide k-mers. 


![**Fraction of k-mers specific to taxonomic rank**
For the GTDB-RS202 database, the majority of nucleotide k-mers are specific to (unique at) a specific genome, species, or genus. Few k-mers are shared across superkingdoms, though these do exist at k=21. K-mers shared at such a high level are indicative of high k-mer homoplasy: the presence of k-mers that are identical by chance rather than evolutionary descent. ](images/gtdb-rs202.lca_f_aggregated_kmers.png){#fig:gtdb-kmers height=2in}

<!---
GTDB rankinfo: 

  - xx% of DNA k-mers (k=21) are shared within-species
  - yy% of protein k-mers are shared within-species
  - zz% of DNA k-mers are shared within-genus ... etc 


Percent/Number of shared k-mers between members of same species/genus?







<!--compare heatmap w/ max containment for subset of gtdb data?-->
<!---
For , e.g. Pseudomonas, XX% of k-mers are 	shared within the chosen/published genomes within species. For all published genomes within the genus, a median of xx% of k-mers are shared between genomes of one species and genomes of the a different species in the same genus.


== median or mean containment at rank?
containent = % of a genome's k-mers that are shared
-- do using ALL of gtdb, BUT, start with just a single set of genomes.. e.g. Pseudomonas? == similar to "shared k-mers" paper [@doi:10.24072/pci.genomics.100001]

![**Protein k-mer containment facilitates genus-level comparisons**
10k pseudomonas genome sequences, median containment at each alphabet](images/pseudomonas_jaccard_vs_containment_prot10.png){#fig:evolpathsContain}

The containment index enables accurate sequence distance estimation between datasets of different sizes, which here provides an added benefit: given a set of trusted proteins (e.g. known proteome), 6-frame translation of DNA query sequences can be used for accurate distance estimation compared with the known/trusted proteome. 

Unlike Jaccard comparisons, which estimate the similarity between sets, containment estimates are relative to each individual set. 
When one set is highly trusted, such as a reference genome or proteome, the containment relative to that set may be most informative.
In these cases, we can consider the trusted genome as an "anchor" upon which we are basing our comparison, and the containment relative to this set as "anchor containment."

This property enables containment comparisons to provide more informative comparisons between sets of different sizes. In cases where one set is more highly trusted to contain accurate k-mers, the containment relative to that set can be more informative. 

This provides some advantages: in cases where one set is more highly trusted to contain accurate k-mers, the containment relative to that set will provide a better estimate than the 

--->















<!---
### Lost Bits


Long dna k-mers ~~ short protein k-mers

while shorter dna k-mers might be shared across more sequence, you increase the risk for result in"shared, non-homologous k-mers" (k-mer homoplasy). A protein k-mer of length `10` coverse 30 base pairs in nuof nucleotide sequence 

### Phylogenetic Reconstruction from k-mer Amino Acid identity
NO, just leave this out

![**K-mer Based Sequence Identity by Lowest Common Taxon**
GTDB Evolpaths dataset](images/anchor-mcANI-AAI.boxen.protnucl.png){#fig:evolpathsANIAAI}
--->


## Discussion

K-mer based estimation of sequence identity has been limited to nucleotide sequences of similar size with high sequence identity (>80%),outside of which MinHash Jaccard is less well correlated with sequence identity [@doi:10.1186/s13059-016-0997-x; @doi:10.1038/s41467-018-07641-9].

By leveraging the Containment Index of Scaled MinHash sketches with both nucleotide and protein k-mers, we can extend accurate k-mer sequence identity to sequences of different sizes and to >50% Amino Acid Identity.


Cricuolo [@doi:10.12688/f1000research.26930.1] (suggests w/ appropriate correction, nucl MinHash Jaccard can be used up to >65% ANI??)

Here, we utilize Scaled MinHash sketches with Containment to overcome size differences between sequences being compared. 

To accurately estimate sequence identity from sequence files of different sizes(genomes, metagenomes, etc), we employ Scaled Minhash sketches, which enables estimation of the Containment Index. 


A number of methods have used discriminatory k-mer analysis for taxonomic classification. However, most rely upon first developing a reference of discriminatory k-mers, e.g. k-mers unique to / diagnostic of a taxonomic group.
Instead, sourmash gather leverages the Containment Index to find the reference match that shares the largest number of k-mers with the query sequence.

At k=21 (dna) and k=7 (protein), many k	-mers are shared across taxonomic groups.
Unlike many k-mer based classifiers, we do not need to explicitly characterize the discriminatory k-mers for each taxonomic group.
The Containment Index uses all matched k-mers between the query and each reference, finding the % of each reference genome present in the query.
Gather then selects the most covered (highest percent contained) reference genome, thus utilizing the combination of shared and discriminatory k-mers to find the most parsimonious match.
After finding the best match, all matched k-mers are removed for the query in order to repeat the analysis to find the next most parsimonious genome match.



While this method is still dependent on a good set of reference genomes, updating the set of references with new data does not require recalculation of discriminatory k=mer sets...

** discussion of k-mer size **

- Scaled Minhash distance estimation is robust to completeness
(unlike standard minhash https://drep.readthedocs.io/en/latest/choosing_parameters.html#importance-of-genome-completeness)





## Conclusions

Containment-based pairwise distance estimation via Scaled Minhash enables accurate assembly-free and alignment-free phylogenomic reconstruction and taxonomic classification across a wide range of evolutionary distances.

## Methods

### Scaled MinHash Sketching with Sourmash

As implemented in sourmash [@https://dib-lab.github.io/2020-paper-sourmash-gather; @doi:10.12688/f1000research.19675.1; @doi:10.21105/joss.00027], Scaled MinHash is a MinHash variant that uses a scaling factor to subsample the unique k-mers in the dataset to the chosen proportion (1/`scaled`).
As k-mers are randomized prior to systematic subsampling, Scaled MinHash sketches are representative subsets that can be used for comparisons, as long as the k-mer size and chosen scaled value remain consistent. 
Unlike traditional MinHash sketches, Scaled MinHash sketches enable similarity estimation with containment, which permits more accurate estimation of genomic distance when genomes or datasets differ in size [@doi:10.1016/j.amc.2019.02.018;@doi:10.1186/s13059-019-1875-0]. 

Sourmash v4.x supports sketching from either nucleotide or protein input sequence.
All genome sequences were sketched with sourmash v4.0 using the `sourmash sketch dna` command, k-mer sizes of 21,31,51, a scaling factor of 1000. 
Sourmash also supports 6-frame translation of nucleotide sequence to amino acid sequence.
To assess the utility of these translated sketches, genome sequences were also sketched with the `sourmash sketch translate` command at protein k-sizes (_kaa-mer sizes?_) of 7-12 and a scaling factor of 100. 
All proteome sequences were sketched with sourmash v4.0 using the `sourmash sketch protein` command at protein k-sizes (_kaa-mer sizes?_) of 7-12 and a scaling factor of 100.
Where higher scaling factors were evaluated, these original sketches were downsampled using the sourmash `downsample` method prior to conducting sequence similarity comparisons.


### Sequence Identity Estimation from Scaled MinHash
_(very DRAFTy)_

Sourmash contains standard implementations of Jaccard Index [@doi:10.1186/s13059-016-0997-x] and Containment Index [@doi:10.1016/j.amc.2019.02.018] set comparisons.

**Estimating Sequence Similarity from Jaccard**
For a comparison between two genomes (genomeA, genomeB), the Jaccard Index represents the k-mers shared between the two genomes (sketch intersection) divided by the k-mers present in both sketches (sketch union).
Thus the Jaccard Index represents the percent of shared k-mers relative to all k-mers across both genomes (intersection/genomeA+genomeB).
MinHash Sketch Jaccard has been shown to correlate well with ANI at high sequence identities (>=90% sequence identity) [@doi:10.1186/s13059-016-0997-x]; (>=80% sequence identity [@doi:10.1038/s41467-018-07641-9].

**Mash Distance from Scaled MinHash Jaccard**

_TBD_



**Estimating Sequence Similarity from Containment**
As the Jaccard Index utilizes the union of all k-mers in a dataset, it is greatly affected by differences in dataset size [@doi:10.1093/bib/bbz083].
The Containment Index instead represents the percent of a genome found in the comparison genome.
Containment is directional: while the number of shared k-mers is fixed for a pairwise comparison, the Containment of each dataset will depend on the unique k-mers found in that particular dataset. Containment for genomeA will be (intersection/genomeA), while Containment for genomeB will be (intersection/genomeB).

Alignment-based ANI represents the sequence similarity of the alignable fraction of two genomes. In this way, ANI only compares the shared sequences, and discounts/ignores all other sequence present in either genome.
Bidirectional containment comparisons use the same numerator (shared k-mers), but may contain different numbers of non-shared k-mers in the denominator.

In cases where both genomes are high-quality and highly complete, we can most closely approximate ANI by using the maximum value between the bidirectional containment values: that is, using the comparison that represents the shared sequence over the genome with the smallest number of non-shared k-mers.

In cases where one genome is more trusted (high quality and highly complete), Containment may be best calculated relative to the trusted genome.
This use case also allows us to estimate sequence identity from larger sequence collections, such as metagenomes.
By definition, metagenomes contain k-mers from many organisms.
We can take advantage of directional Containment by calculating the Containment Index of Reference genomes that share many k-mers with the Metagenome.
We have already shown the utility of Containment for metagenome classification [@https://dib-lab.github.io/2020-paper-sourmash-gather], but now we can report estimated average sequence identity between the matching sequence regions and the reference genome.

**Estimating Sequence Identity from Scaled MinHash**

**_TBD_**

Blanca et al, 2021 [@doi:10.1101/2021.01.15.426881] presented a method to estimate the mutation rate between MinHash sketches while accounting for the non-independence of mutated k-mers. Using [@https://github.com/KoslickiLab/mutation-rate-ci-calculator], we estimate Sequence Identity from Scaled MinHash Containment.

Estimating sequence similarity from Scaled MinHash requires a good estimate of the number of unique k-mers in the sketched sequencing dataset [@https://github.com/dib-lab/sourmash/pull/1270]...







### Scaled MinHash Distance Correlates with Standard Methods

FastANI v1.32 ([@doi:10.1038/s41467-018-07641-9]; run with default parameters)  was used to obtain Average Nucleotide Identity between the anchor genome and each additional genome in its evolutionary path.
FastANI is targeted at ANI values between 80%-100%, so only values in this range are considered "trusted" and used in **assessing the correlation between Scaled MinHash estimates and FastANI._(TBD)_**

CompareM v0.1.2 ([@url:https://github.com/dparks1134/CompareM]; run with `--sensitive` parameter for DIAMOND mapping) was used to obtain Average Amino Acid Identity between the anchor proteome and each additional proteome in its evolutionary path.
CompareM reports the mean and standard deviation of AAI, as well as the fraction of orthologous genes upon which this estimate is based.
Briefly, CompareM calls genes for each genome or proteome using PRODIGAL [@doi:10.1038/nmeth.3176] and conducts reciprocal best-hit mapping via DIAMOND [@doi:10.1186/1471-2105-11-119].
By default, CompareM requires at least 30% percent sequence identity and 70% percent alignment length to identify orthologous genes.
As DIAMOND alignment-based homology identification may correlate less well with BLAST-based homology under 60% sequence identity [@url:https://rodriguez-r.com/blog/aai-blast-vs-diamond/], **we also ran compareM with a percent sequence identity threshold of 60% to obtain a set of high-confidence orthologous genes for AAI estimation. We report correlation between Scaled MinHash AAI estimation and each of these compareM parameter sets in XX _(TBD)_**. _CompareM was also used to obtain AAI values directly from each genome, using PRODIGAL to translate sequences prior to gene calling. These results [were not significantly different from proteome-based AAI estimation??] (Supplemental XX)._


### Taxonomic Classification with Sourmash `Gather` and `Taxonomy`

To take advantage of the increased evolutionary distance comparisons offered by protein k-mers, we apply compositional analysis with sourmash gather [@https://dib-lab.github.io/2020-paper-sourmash-gather] to protein sequences (amino acid input and 6-frame translation from nucleotides).
Sourmash gather is conducted in two parts: 
First (preselection), gather searches the query against all reference genomes, building all genomes with matches into a smaller, in-memory database for use in step 2.
Second (decomposition), gather does iterative best-containment decomposition, where query k-mers are iteratively assigned to the reference genome with best containment match.
In this way, gather reports the minimal list of reference genomes that contain all of the k-mers that matched any reference in the database.

For reference matches with high sequence identity (ANI) to the query, we classify the query sequence as a member of the reference taxonomic group, as in [@https://dib-lab.github.io/2020-paper-sourmash-gather].
**However, when ANI between the query and the top reference match exceeds the taxonomic rank threshold (e.g. species default 95%), we use a least/lowest common ancestor (LCA) approach to report likely taxonomy at a higher taxonomic rank _(TBD)_**.
Briefly, as gather reports non-overlapping genome matches, we can sum the k-mer matches for all genomes with shared taxonomies at the next higher taxonomic rank to report the best query containment at that rank.
As this gather-LCA approach first uniquely assigns k-mers to their best reference genome, it bypasses the impact of increasing database size on taxonomic assignment observed for other LCA-based k-mer classification approaches [@doi:10.1186/s13059-018-1554-6].


### Workflows and Computing Resources

Reproducible workflows associated with this paper are available at XX (gh link + doi for release), with datasets available at OSF (XX). All workflows were executed using snakemake >= 5.26 [@doi:10.12688/f1000research.29032.1)] on the FARM cluster at UC Davis, using practices outlined in [@doi:10.1093/gigascience/giaa140].













## References {.page_break_before}

<!-- Explicitly insert bibliography here -->
<div id="refs"></div>
