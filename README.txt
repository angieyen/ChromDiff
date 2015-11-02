## README file created by Angela Yen on March 13, 2015
## Please contact angela@mit.edu with questions or for clarifications.

# Copyright 2015 Angela Yen
# # This file is part of ChromDiff.
# # ChromDiff is free software: you can redistribute it and/or modify
# # it under the terms of the GNU General Public License as published by
# # the Free Software Foundation, either version 3 of the License, or
# # (at your option) any later version.
#
# # ChromDiff is distributed in the hope that it will be useful,
# # but WITHOUT ANY WARRANTY; without even the implied warranty of
# # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# # GNU General Public License for more details.
#
# # You should have received a copy of the GNU General Public License
# # along with ChromDiff.  If not, see <http://www.gnu.org/licenses/>.

###############  Part 0: Setup ####################
This folder "ChromDiff" contains the code and data for running ChromDiff on the examples provided in the manuscript.

Part 0.1: Dependencies for ChromDiff include BedTools, awk, bash, and R. Please make sure these are installed and working before trying to run ChromDiff.

Part 0.2a: Applying ChromDiff to Epigenomics Roadmap data, as described in Yen and Kellis, 2015: 
The only data that was not included in this zipped file is the 15-state ChromHMM state calls. The easiest way to download all these chromatin state calls is 
to download the compressed .tar.gz directory at http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz. 
Then, uncompress the directory and move all the contained files into the ChromDiff subdirectory statecalls/core/, which has already been created within this directory. 
Alternatively, you can download all state call files of the form E***_15_coreMarks_mnemonics.bed.gz from 
http://www.broadinstitute.org/~anshul/projects/roadmap/chromhmmSegmentations/ChmmModels/coreMarks/parallel/set2/final/ into the subdirectory statecalls/core/, which has already been created within this directory. 
Generally, we are using the "Mnemonics bed files" of the "Core 15-state model_ described here: http://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state

Once the state calls have been downloaded an dmoved into the corresponding directory, simply run ./notes_v2.sh to
reproduce the analysis, results, and figures. (NOTE: as described in Yen and
Kellis, 2015, our gene list does not include chrY genes. Therefore, if you
download our gene list in data/gencode_genes_full.txt, this will be a list of all
Version 10 GENCODE protein-coding genes, except for ChrY genes.) All input
files can be found in the data/ folder, which can be used as examples as
described in Part 1 of this README file. 
--> To find results and figures, please see Part 4 of this README file.

Part 0.2b: Applying ChromDiff to other data:
- Part 1 describes preparation and formatting of the data and other input for
  ChromDiff
- Part 2 describes how to specify the input variables and files to ChromDiff.
- Part 3 outlines the various steps of the ChromDiff pipeline
- Part 4 describes where to find results, p-values, and figures.

###########   PART 1: Formatting for other data for ChromDiff input  #############
ChromDiff relies on the following information
1. Epigenome unique IDs and metadata 
2. Gene IDs and locations
3. Chromatin state annotations for the epigenome IDs
4. Expression data in the form of RPKM values for each gene
5. Covariate matrix (based on metadata) and mapping of covariates to metadata
categories

See below for details:
1. Epigenome unique IDs and metadata 
This file should be in the following tab-delimited format:
ID	filepath	color	name	category1	category2	...
id1	filepath/to/id1.bed	#E41A1C	name_of_id1	category1_for_id1	category2_for_id1 ...
id2	filepath/to/id2.bed	#924965	name_of_id2	category1_for_id2	category2_for_id2 ...
.
.
.

That is, the first column is the unique ID, the 3rd column is the hex coding
of the color for the epigenome, and the fourth column is the name
(label that will be used on plots). The second column is the path to the 
chromatin state annotation file (as described in "3. Chromatin state
annotations" below). The remaining columns are any other
categories you will want to compare the data based on. 
The first row must be the "column names" for this matrix, with "ID",
"filepath", and "name" as the first three column names, as shown above. 
NOTE: no periods (.) should be used in the metadata as epigenome IDs, category
names, or category values. That means no periods should be used in 1st column,
or in any columns after the 3rd column. (Periods may and will likely be used
in the filepath column.)
NOTE: every category name should be unique from every other category
name. (i.e. the 5th and 6th column names can not be the same.) 
NOTE: Every option value should be uniquely assigned to one category. (For example, there
can't be a "TRUE" for category "is_female" and also a "TRUE" entry for
category "is_tissue". This can be avoided by using values such as
"TRUEFEMALE" and "TRUETISSUE", or simply "FEMALE" and "TISSUE".)


Example: final_celltype_metadata.txt

2. Gene IDs and locations
The gene file must be modeled in the following tab-delimited format:
#All gencode protein-coding genes (entries with type gene)
#chr	start	end	strand	gencode_gene_id	genesymbols
chr1	34553	36081	-	ENSG00000237613.2	FAM138A
chr1	69090	70008	+	ENSG00000186092.4	OR4F5
chr1	367639	368634	+	ENSG00000235249.1	OR4F29

The 6th column of genesymbols are optional, but are necessary for gene set
enrichment calculations.  
Example: genes/gencode_genes_full.txt

3. Chromatin state annotations 
a) Chromatin state annotations should be in the following tab-delimited bed format: 
chr10	0	119600	15_Quies
chr10	119600	120400	1_TssA
chr10	120400	136200	14_ReprPCWk
chr10	136200	139400	15_Quies
chr10	139400	145200	9_Het
chr10	145200	162800	15_Quies
chr10	162800	165000	9_Het
chr10	165000	176200	15_Quies
chr10	176200	176600	7_Enh

NOTE: the 4th column must be of the form STATENUM_STATEMNEMONIC, and these
must match the numbers/mnemonics given in the chromatin state metadata file
(Part1.3b).

Example: http://www.broadinstitute.org/~anshul/projects/roadmap/segmentations/models/coreMarks/parallel/set2/final/E001_15_coreMarks_mnemonics.bed.gz

b) Chromatin state metadata
A tab-delimited file listing all possible chromatin states with the following information:
STATE NO.       MNEMONIC        DESCRIPTION     COLOR NAME      COLOR CODE
1       TssA    Active TSS      Red     "255,0,0"
2       TssAFlnk        Flanking Active TSS     Orange Red      "255,69,0"
3       TxFlnk  Transcr. at gene 5' and 3'      LimeGreen       "50,205,50"
4       Tx      Strong transcription    Green   "0,128,0"
5       TxWk    Weak transcription      DarkGreen       "0,100,0"
6       EnhG    Genic enhancers GreenYellow     "194,225,5"
7       Enh     Enhancers       Yellow  "255,255,0"
8       ZNF/Rpts        ZNF genes & repeats     Medium Aquamarine
"102,205,170"
9       Het     Heterochromatin PaleTurquoise   "138,145,208"
10      TssBiv  Bivalent/Poised TSS     IndianRed       "205,92,92"
11      BivFlnk Flanking Bivalent TSS/Enh       DarkSalmon      "233,150,122"
12      EnhBiv  Bivalent Enhancer       DarkKhaki       "189,183,107"
13      ReprPC  Repressed PolyComb      Silver  "128,128,128"
14      ReprPCWk        Weak Repressed PolyComb Gainsboro       "192,192,192"
15      Quies   Quiescent/Low   White   "255,255,255"

NOTE: State numbers must go from 1...n, where n is the total number of states.

Example: core_annotation.txt

4. Expression data should be formatted as a matrix in tab-delimited
format, with each row corresponding to a gene and each column corresponding to
a epigenome. The first row should be the labels. An example is below:
gene_id	E000	E003	E004	E005	E006
ENSG00000000003	23.265	43.985	37.413	29.459	21.864
ENSG00000000005	0.872	1.642	6.498	0.000	0.157
ENSG00000000419	55.208	35.259	58.308	48.208	37.477
ENSG00000000457	3.237	2.596	2.345	8.775	2.723
ENSG00000000460	7.299	6.649	7.838	7.324	0.830
ENSG00000000938	0.052	0.211	0.059	0.009	0.012
ENSG00000000971	3.000	0.000	0.003	0.009	41.909
ENSG00000001036	58.371	49.337	25.066	25.679	86.255
ENSG00000001084	7.077	6.966	10.110	4.609	5.870

Example: rnaseq/57epigenomes.RPKM.pc

5.
a) Covariate matrix for factors that will be corrected for. Values can be
factors (which will be automatically converted by R), or numbers for "mixed"
covariates:
UCSD	BI	UCSF-UBC	UW	FEMALE	MALE	UNKNOWN	SOLID_LIQUID	TYPE
E017	1	0	0	0	1	0	0	NONE	CELLLINE
E002	0	1	0	0	1	0	0	NONE	PRIMARYCULTURE
E008	1	0	0	0	1	0	0	NONE	PRIMARYCULTURE
E001	0	1	0	0	1	0	0	NONE	PRIMARYCULTURE
E015	0	1	0	0	1	0	0	NONE	PRIMARYCULTURE

Example: cov.mat.txt

b) Mapping metadata matrix to covariate matrix. This must explain which
metadata columns (i.e. "category1", "category2", ...) correspond to which
covariate columns. Matrix has n rows and m columns, where n is the number of
the metadata columns (categories) that we may correct for, while m is the
number of covariate columns. Each matrix cell is TRUE for elements
matrix[i,j] where metadata column i corresponds to covariate column j, and
FALSE otherwise:
UCSD	BI	UCSF.UBC	UW	FEMALE	MALE	UNKNOWN	SOLID_LIQUID	TYPE
LAB	TRUE	TRUE	TRUE	TRUE	FALSE	FALSE	FALSE	FALSE	FALSE
SOLID_LIQUID	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	TRUE
FALSE
TYPE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	TRUE
SEX	FALSE	FALSE	FALSE	FALSE	TRUE	TRUE	TRUE	FALSE	FALSE


Example: map_vars_covariates.txt

######### PART 2: Specifying input files and variables ###################
ChromDiff requires that input files and variables be set in the file notes_v2.sh. The file currently contains an outline of the pipeline, with inputs currently set to the example files and parameters that generated the sex comparison (Female vs Male) described in Yen and Kellis, 2015. The description of the required inputs is described below, in Part2a-Part 2e. 

These required inputs for ChromDiff can roughly be divided into two parts, as
shown by "Part 1" and "Part 2" in the notes_v2.sh. Part 1, which is referenced
in sections a-b below, has to do with processing the data and generating the
ChromDiff representation. Part 2, which is referenced in sections c-e, has to
do with the parameters for each epigenomic group comparison (such as
Female/Male, Adult/Fetal, etc.).

a. Labels: In notes_v2.sh, you need to specify 2 names/labels. Specify statecalls_label to any name or label that describes your
epigenomic (chromatin state call) information. Set generegions_label to any name or label that
describes your gene regions. These labels are so that you can run the
ChromDiff pipeline on different chromatin state calls or different epigenomic
information, as well as different gene regions. 

Example (as shown in "Part 1" of notes_v2.sh): 
statecalls_label="core"
generegions_label="gencode_v10"

b. File inputs: In notes_v2.sh, you also need to specify the following file
inputs. They must be formatted as described in the corresponding section of
Part 1 above:

	i. metadatafile: epigenomic metadata file as described in Part 1.1 
	ii. genefile: file of gene regions as described in Part 1.2
	iii. state_annotations_file: metadata for chromatin states, as described in
	Part 1.3b. 
	iv. expfile: expression matrix, as specified in Part 1.4.
	v. covariate_mat_file: Covariate matrix specifying the covariates and
	associated epigenomes, as described in Part 1.5a
	vi. map_covariates_file: Matrix specifying relationships between
	covariates (in the covariate matrix) and metadata categories/properties,
	as specified in 1.5b

Example (as shown in "Part 1" of notes_v2.sh):
metadatafile="data/final_celltype_metadata.txt"
genefile="data/gencode_genes_full.txt"
states_info="data/core_annotation.txt"
expfile="data/57epigenomes.RPKM.pc"
covariate_mat_file="data/cov.mat.txt"
map_covariates_file="data/map_vars_covariates.txt"

c. Epigenomic groups: To compare groups of epigenomes to each other, you
need to specify the property and options you want to compare Specifically, the
property is the column name (category) in your epigenome metadata matrix (i.e.
"sex"). The two options are the two (different) values for that column you
want to compare your epigenomes based on (i.e. "Female" vs "Male"). Then,
ChromDiff will compare the group of epigenomes labeled a_option under property to the group of epigenomes labeled b_option under property.

Example (as shown in "Part 2" of notes_v2.sh):
property="sex"
a_option="Female"
b_option="Male"

d. ChromDiff parameters: ChromDiff allows two parameters to be set. 
	i. test_type: The test_type variable specifies the probabilistic test ChromDiff uses to
	calculate p-values, as described in Yen and Kellis, 2015. Options are
	"wilcox", "ttest", and "ftest", for the Mann-Whitney-Wilcoxon test, the
	Student's t-test, and the F-test of the equality of two variances,
	respectively.

	ii. correction: The correction variable specifies the multiple hypothesis
	correction used on the raw p-values, as described in Yen and Kellis, 2015.
	Options are "fdr", "bonferroni", and "BY", for the Benjamini-Hochberg FDR
	correction, the Bonferroni (familywise error rate) correction, and the
	Benjamini-Yekutieli FDR correction, respectively.

Example (as shown in "Part 2" of notes_v2.sh):  
test_type="wilcox"
correction="fdr"

e. Optional setting of dendrogram height cutoff for clustering: heightcutoff
can default to 0, which will produce no clusters, unless there are
features/genes with identical epigenomic information. 

Optionally, after running the analysis the first time (with default height cutoff of 0), 
you can then look at the resulting dendrogram file in 
plots/${curr_label}/${a_option}.${b_option}/perc/${test_type}/${correction}/hclust_domstate.pdf 
to choose a manual cutoff. Then, the height cutoff will generate and analyze gene clusters based on 
the given dendrogram and height cutoff, as shown in clustered plots and
analysis in Yen and Kellis, 2015. 

Example (as shown in "Part 2" of notes_v2.sh):
heightcutoff=0


#### PART 3: Running ChromDiff pipeline ####

After setting up your data for input (as in part 1) and specifying your input
files, variables and parameter (as in part 2), we can proceed to actually
running ChromDiff on your data.

The simple way to run ChromDiff, is simply to run the script notes_v2.sh. 

There, you can see that ChromDiff splits up into 4 steps: 

Step 1: Setting up gene info:
	./process_gene_info.sh $genefile $generegions_label
	
Step 2: Calculate features for all epigenomes, chromatin states, and genes:
	./calculate_all_raw_features.sh $metadatafile $statecalls_label $generegions_label $states_info

	-Note: You can also use ./calculate_raw_features.sh on each epid separately to
parallelize:
	./calculate_raw_features.sh $epid $statecallfile $statecalls_label $generegions_label $states_info

Step 3: Calculate background information, feature names, and check data:
	./featnames_bgvals.sh $statecalls_label $generegions_label $states_info

Step 4: Run ChromDiff on a particular epigenomic group comparison: 
	./perform_analysis.R $property $a_option $b_option $test_type $correction $curr_label $heighcutoff $metadatafile


RESULTS: 
 - figures and plots as presented in Yen and Kellis, 2015 will output in the
   directory: 
   plots/core/${label1}_${label2}/perc/${test}/${correction}/

 - results for comparison of expression can be found in the directory
   exp.pvals/

 - results for msigdb gene set enrichment can be found in msigdb/results/
 
#### PART 4: Finding results generated by ChromDiff #######
There are 4 categories of results that ChromDiff generates, which are listed
below with corresponding directories:
1. P-values of features (all_pvals/)
2. Plots and figures (plots/)
3. Differential expression results (exp.pvals/)
4. Gene set enrichment results (msigdb/)

1. P-values:
a. Corrected p-values for every feature (gene and chromatin state)
for each comparison can be found in the following path:
--> all_pvals/${statecalls_label}_${generegions_label}/perc.${test}.${correction}.${property}.${label1}.${label2}.txt
All ${xxxxx} variables should be substituted for the values used, as
specified in Part 2. 

b. Corrected p-values for significant features that were used in figure generation 
(sampled as necessary) can be found in the following path:
sig_pvals/${statecalls_label}_${generegions_label}/perc.${test}.${correction}.${property}.${label1}.${label2}.txt

The format of the p-value files is two columns, with the left column as the
feature id (GENEID_CHROMATINSTATENUMBER) and the right column as the p-value
after multiple hypothesis correction. 

2. Plots and figures: Plots and figures used in Yen and Kellis, 2015 will
output in the following directory:
plots/${statecalls_label}_${generegions_label}/${label1}.${label2}/perc/${test}/${correction}/

a. Plot for combinations of genes and chromatin states (as shown in Fig. 2a): sig.pvals.pdf  
b. Dendrogram that can be used to choose a height cutoff (as chosen in Table
S1): hclust_domstate.pdf
c. Plot of dominant chromatin state (as shown in Fig. 2c): sig_maj_plot_clust_matched_domstate.pdf
d. Plot of expression (as shown in Fig. 2d): sig_exp_plot_clust_matched_domstate.pdf 
(Supplementary) e. Plot of dominant chromatin state ordered to match gene ordering from
combinations of genes and chromatin states (as shown in Fig. S2b):
sig_maj_plot_clust_matched_combinations.pdf

3. Differential expression: Statistics that show how many genes were
differentially expressed are listed in the file
exp.pvals/${statecalls_label}_${generegions_label}_domstate_percent_all.txt
(Tables are also available in the same directory.)

4. Gene set enrichment: Results for gene set enrichment from MSigDB gene sets
are in the following directory:
msigdb/results/c2_c5/symbols/perc.${test}.${correction}.${property}.${label1}.${label2}_all.txt

The file is in table format with the most significant gene sets listed first.
As shown in the header, the columns represent:
geneset	genesetsize(K)	numoverlaps(k)	k/K	q-value	orig_pval

P-values are generated using the hyper-geometric test, and Storey's q-value
correction is used to generate the q-values shown in the 4th column. 

