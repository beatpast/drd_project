# drd_project
This is my project of drd (dna and rna dynamics) focused on data extracted from Infinium HumanMethylation450 BeadChip experiment. a dataset of DNA methylation at genome-wide CpG levels was processed with R programming language (minfi package) for UNIBO master's degree in Bioinformatics. 
DNA methylation plays an important and dynamic role in regulating gene expression. it is an epigenetic change involving the covalent transfer of a methyl group to the C-5 position of the cytosine ring of DNA 
that could be a biomarker for cancer for example.
This technology has two types of probes: infinium I and infinium II with different design strategies. Every probe contains data in both colors red and green fluorescence intensity even if for example green was not emitted. the output of this data is important for normalization(separate background signal from real signal)
in type I probes the color doesn't give info about the methylation status. AddressA (identify beads with binding unmethyl) and B(=identify beads with binding to methylated version of the genome)
the pipeline follows these steps:
1. Load raw data with minfi and create an object called RGset storing the RGChannelSet object
2. Create the dataframes Red and Green to store the red and green fluorescences respectively
3. Check the Red and Green fluorescences for the address assigned to me (45652402) 
4. Create the object MSet.raw
5. Perform the following quality checks: QCplot, control probes, detection p-values; for each sample check how many probes have a detection p-value higher than the threshold assigned to me(0.01)
6. Calculate raw beta and M values and plot the densities of mean methylation values, dividing the samples in WT (wild type) and MUT (mutation)
7. Normalize the data using the function assigned (preprocessFunnorm) then show the density plots of beta mean values according to the chemistry of the probes, the density plot of beta standard deviation values according to the chemistry of the probes and the boxplot of beta values
8. Perform a PCA on the matrix of normalized beta values generated in step 7, after normalization
9. Using the matrix of normalized beta values generated in step 7, identify differentially methylated probes between group WT and group MUT using the function assigned Mann-Whitney test
10. applying the Benjamini & Hochberg and the Bonferroni corrections
11. Produce a volcano plot and a Manhattan plot of the results of differential methylation analysis
