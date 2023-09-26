# Single-cell characterization of immune aplastic anemia

Scripts to reproduce figures and analyses related to the manuscript "Single-cell characterization of immune aplastic anemia"

# Installation

Clone this github repository, by e.g., Bash in terminal by

<pre>

git clonehttps://github.com/janihuuh/aa_manu

</pre>

All the other dependencies (R-packages, Python modules) are noted in the code and their installation guides can be found in their own respective guides. 

The typical installation time for all the software on a standard laptop (e.g., Apple M1 16Gb) should in general not exceed 1 hour. 

# Pseudocode for AA manuscript

## Single-cell RNA+TCR&alpha;&beta;-seq analysis; <i>in vivo</i> data

Briefly, for the scRNA+TCR&alpha;&beta;-seq, we profiled the samples in the following way:

1. Preprocessing and quality control (QC) of scRNA+TCR&alpha;&beta;-seq data, including i) AA BM data, ii) healthy BM data, iii) comparison data
2. Merge i) AA BM data with ii) healthy BM data with scVI, annotate with SingleR
3. Identify broader cell lineages (CD4+, CD8+, NK, B, myeloid, HSPC) and delineate
4. To delineated cells, merge iii) comparison data with scVI, cluster with Seurat, annotate manually

For more detailed description, see <b>Methods</b> in the manuscript or the code.
 
### 1. scRNA+TCR&alpha;&beta;-seq preprocessing

Preprocessing and quality control (QC) of scRNA+TCR&alpha;&beta;-seq data, including i) AA BM data, ii) healthy BM data, iii) comparison data

<pre>
<b>for</b> aa_tcr_data, <i><b>do</b></i>  
  <b>read</b> 10X CellRanger output 
  <b>filter</b> based on the quality of cells ## e.g., incomple TCR&alpha;&beta;-seq information, low confidence data  

<b>for</b> public_tcr_data, <i><b>do</b></i>  
  <b>read</b> 10X CellRanger output 
  <b>filter</b> based on the quality of cells ## e.g., incomple TCR&alpha;&beta;-seq information, low confidence data  
  
<b>for</b> aa_rna_data, <i><b>do</b></i>  
  <b>read</b> 10X CellRanger output ## 13 AA BM samples  
  <b>filter</b> based on the quality of cells ## common thresholds for all donors  

<b>for</b> healthy_rna_data, <i><b>do</b></i>  
  <b>read</b> 10X CellRanger output 
  <b>filter</b> based on the quality of cells ## common thresholds for all donors  
 
<b>for</b> public_rna_data, <i><b>do</b></i>  
  <b>read</b> 10X CellRanger output 
  <b>filter</b> based on the quality of cells ## common thresholds for all donors  

</pre>

### 2. Initial scRNA+TCR&alpha;&beta;-seq analysis

Merge i) AA BM data with ii) healthy BM data with scVI, annotate with SingleR. 

These findings can be mainly seen in <b>Figure 1</b>.

<pre>
total_rna_data = <b>merge</b> aa_data with public_data
total_rna_tcr_data = <b>merge</b> total_rna_data with aa_tcr_data
</pre>

<pre>
<b>for</b> total_rna_tcr_data, <i><b>do</b></i> 
  <b>annotate</b> cells with SingleR  
  <b>calculate</b> latent represenation ## with scvi  
  <b>calculate</b> UMAP representation from latents  
  <b>scale</b> data  
</pre>

<pre>
<b>for</b> total_rna_tcr_data, <i><b>do</b></i>  
  <b>calculate</b> DE-genes, pathways between AA and healthy  
  <b>calculate</b> ligand-receptor pairs with CellPhoneDB  
</pre>
	
	
### 3. Delineate scRNA+TCR&alpha;&beta;-seq, add comparison groups

Identify broader cell lineages (CD4+, CD8+, NK, B, myeloid, HSPC) and delineate

<pre>
<b>for</b> total_rna_tcr_data, <i><b>do</b></i>  
	cd4_data = <b>subset</b> CD4 cells
	cd8_data = <b>subset</b> CD8 cells
	nk_data = <b>subset</b> NK cells
	b_data = <b>subset</b> B cells
	myeloid_data = <b>subset</b> Myeloid cells
	hspc_data = <b>subset</b> HSPCs
</pre>

<pre>
delineated_data = list cd4_data, dc8_data, nk_data, b_data, myeloid_data, hspc_data

<b>for</b> delineated_data, <i><b>do</b></i>  
  <b>calculate</b> latent represenation ## with scvi  
  <b>calculate</b> UMAP representation from latents  
  <b>calculate</b> clustering from latents  
  <b>decide</b> optimal clustering ## based on visual inspection  
  <b>scale</b> data  
  
<b>for</b> delineated_data, <i><b>do</b></i>  
  <b>calculate</b> DE-genes, pathways between clusters  
  <b>calculate</b> DE-genes, pathways between AA and healthy  
  <b>calculate</b> ligand-receptor pairs with CellPhoneDB 
</pre>


## Single-cell RNA+TCR&alpha;&beta;-seq analysis; <i>in vitro</i> data

In the manuscript, we performed dynamic co-cultures of immune cells with AA target cells, HSPCs, with multiplexted scRNA+TCR&alpha;&beta;-seq readout. Briefly, 

1. Preprocess samples
2. Demultiplex
3. Merge data with scVI, annotate with SingleR
4. Identify broader cell lineages (CD4+, CD8+, NK, B, myeloid, HSPC) and delineate, and to delineated cells, recluster with Seurat, annotate manually

See the steps 1, 3, and 4 above. For demultiplexing, 

Demultiplexing of hashtag oligonucleotide (HTO) counts to classify cells to samples was performed on the quality-controlled cells.
Additionally, cells with outlier expression of total HTO counts (>10,000 counts) were removed. Centered log-ratio-normalized HTO UMI counts were used with the “HTODemux”-function in Seurat (v. 4.3.0) with a uniquely chosen positive quantile for each sample, ranging from 0.99 to 0.9999. Only cells classified as singlets and cell types corresponding to the experimental design (e.g., only T cells from CD4/CD8 sorted samples) were retained for further analyses.

### 2. Demultiplex
<b>for</b> sample_rna_data, <i><b>do</b></i>  
  <b>remove</b> cells with HTO counts > 10,000
  <b>remove</b> K562 cells based on expression of marker genes
  <b>normalize</b> with HTODemux
  <b>remove</b> cells not corresponding to well ID (e.g., remove B cells from well that contained only T cells)
  <b>pool</b> samples together
</pre>


## TCR&beta;-seq analysis

In the manuscript, we identified AA-associated TCR&beta;-seq motifs. Briefly, it was done as follows:

1. Preprocess TCR&beta;-seq data
2. Take AA discovery cohort samples, split into 10 equal sized train sets, perform GLIPH2 motif serach for each
3. Filter motifs (have to be found 7/10 sets, cannot be found in other diseases, healthy controls)
4. Validate in other sets (validation TCR&beta;-seq, validation scRNA+TCR&alpha;&beta;-seq cohort)

### 1. Preprocess

<pre>
<b>for</b> tcr_data, <i><b>do</b></i>  
	<b>process</b> with VDJtools  
 	<b>subsample</b> to 40k reads ## to avoid bias in sampling depths
</pre>

### 2. GLIPH2

<pre>
<b>for</b> aa_tcr_data, <i><b>do</b></i>  
  <b>split</b> to 10 equal sized sets  
<b>for</b> aa_tcr_data_sets, <i><b>do</b></i>  
  <b>run</b> GLIPH2  
</pre>

<pre>
other_tcr_data = <b>list</b> TCR&beta;-seq from mds, lgl, ra, pnh, healthy
<b>for</b> other_tcr_data, <i><b>do</b></i>
  <b>run</b> GLIPH2
</pre>

### 3. Filter GLIPH2 motifs

<pre>
aa_motifs_filt = <b>filter</b> aa_motifs found in 7/10 sets
aa_motifs_final = <b>filter</b> aa_motifs_filt not found in other_tcr_motifs
</pre>


## To reproduce the results:

### 1) Clone this repository

```
git clone https://github.com/janihuuh/aa_manu
cd path/to/aa_manu/
```

### 2) Obtain the data

* The processed scRNA+TCR&alpha;&beta;-seq data can be received from ArrayExpress (accession number xxx, submission ongoing). 
* The TCR&beta;-seq data can be received from immuneAccess (accession number xxx, submission ongoing).

### 3) Create Seurat-objects for scRNA+TCR&alpha;&beta;-seq data

<pre>
Rscript R/main.R ## init the helper-functions, coloring, etc.
Rscript R/helper/run_preprocessTCRab.R ## preprocess the TCRab-seq
python python/run_scvi_example.py ## obtain the latent embeddings to use in clustering, UMAPs. This is just an example script
Rscript R/rnaseq/run_createSeurat.R ## read the CellRanger output, filter, merge TCR&alpha;&beta;-data, merge scVI latent embeddings, run singleR, DEs, pathways, etc.
</pre>

### 4) Run additional scRNA+TCR&alpha;&beta;-seq analyses

<pre>
Rscript R/rnaseq/run_cellphone.R ## init data for CellPhoneDb
bash bash/run_cellphonedb.sh ## Run CellPhoneDb
</pre>

### 5) Run additional TCR&beta;-seq  analyses 

<pre>
bash bash/run_vdjtools.sh ## preprocess the bulk-TCR&beta;-data, subsample, calcualte diverisities, etc.
bash bash/run_gliph2_expample.sh ## run GLIPH2 collectively and on individual samples; this is just an example script
python python/run_tcrgp.py ## run TCRGP collectively on samples
Rscript R/tcrseq/run_antigen_drive.R ## run antigen-drive based analyses
</pre>
