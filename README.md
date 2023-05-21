# Antibacterial action of amyloidogenic peptide R23I on proteome of T. thermophilus

***Student: Anna Malykhina*** (Institute of Cytology of the Russian Academy of Sciences, Tikhoretsky ave. 4, 194064 St. Petersburg, Russia)

***Supervisor: Oksana Galzitskaya*** (Institute of Protein Research, Russian Academy of Sciences, 142290 Pushchino, Russia)

## Introduction

The massive use of antibiotics has led to the emergence of bacterial resistance, which traditional antibiotics development, high cost and long-time, is failed to overcome. Researches worldwide try to find new ways to deal with this problem, one of which is studying antimicrobial peptides (AMPs). AMPs are an important part of the natural immune defense mechanism of most organisms against pathogens. It has been reported that AMPs have an excellent broad-spectrum antibacterial and low biological toxicity. Of particular interest are AMPs, which, in addition to evident antimicrobial properties, also tend to self-assembly, form supramolecular complexes and induce amyloidogenesis. Despite traditional negative view on amyloidogenic properties as contributors to the neurodegenerative and progressive metabolic diseases, nowadays AMPs with such properties become acknowledged to have an impact on pathogenic bacteria and viruses.

## Methodology

![](https://github.com/anet1223/Project/blob/main/plots/workflow.png)

The amyloidogenic peptide R23I was selected from a number of peptides due to its pronounced antibacterial effect. The aim of the project was to determine the antibacterial mechanism of action of this substance based on changes in the proteomic profile according to mass spectrometry. Bacteria T. thermophilus were treated with peptide R23I in three dosages (20, 50 and 100 μg/ml) in 3 replicas, control was done in 4 replicas. The samples were appropriately prepared and examined by mass spectrometry. The raw data were processed by two programs Peaks and IdentiPy.

Presented data for analysis are in folder [data](https://github.com/anet1223/Project/tree/main/data).

As the result we got information about the proteins in each sample and the number of peptides detected for each protein. According to the theory, the more peptides from a given protein, the higher the content of this protein in the cell. Thus, further analysis was focused on comparing the number of peptides in the control and the three experimental groups.

The comparison of 4 groups was carried out using the **ANOVA** statistical test, which determined whether there were significant changes in the number of peptides in at least one group. According to the p-value proteins were divided into two subsets: proteins with changes for which the p-value was less than or equal to 0.05 and proteins without changes with p-value greater than 0.05.

The next step of the analysis was a pairwise comparison between the groups using the **post-hoc** test with the Tukey correction.

The code in R for analysis could be found in folder [code](https://github.com/anet1223/Project/tree/main/code).

## Results

Before the comparison, the proteins were filtered so that there was information about at least one replica in the control for a particular protein, and one replica in any treatment group. After this selection, we received 278 proteins from the Peaks analysis, and 273 proteins from the IdentiPy.

After ANOVA comparison we obtained proteins with changes (112 proteins for Peaks and 91 for IdentiPy) and proteins without changes (147 for Peaks and 166 for IdentiPy).

Summary tables with statistical information on proteins for both programs are in folder [results](https://github.com/anet1223/Project/tree/main/results).

Since the Peaks program is the leading program for proteomics, further analysis was carried out according to the Peaks data.

Predominantly the changes were in favor of decrease (108 proteins out of 112) which could be noticed on heatmap:

![](https://github.com/anet1223/Project/blob/main/plots/heatmap_proteins_with_difference.png)

On **PCA plots** the two groups of proteins could be visually separated:

<img src="https://github.com/anet1223/Project/blob/main/plots/PCA_proteins_2d.png" width=50% height=50%>

![](https://github.com/anet1223/Project/blob/main/plots/PCA_proteins_3d.gif)

On **PCA of replicas** significant difference between control and treatment groups could be observed. Moreover, despite the small spread between dosages, some clusters still could be noticed.

<img src="https://github.com/anet1223/Project/blob/main/plots/PCA_replicas.png" width=50% height=50%>

<img src="https://github.com/anet1223/Project/blob/main/plots/PCA_biplot.png" width=50% height=50%>

**Dose effect** was also found in pair-wise comparison:

![](https://github.com/anet1223/Project/blob/main/plots/Dose_effect.png)

GO annatation was taken from [here](https://www.ebi.ac.uk/QuickGO/).

<img src="https://github.com/anet1223/Project/blob/main/plots/GO_location.png" width=50% height=50%>  

<img src="https://github.com/anet1223/Project/blob/main/plots/GO_function.png" width=50% height=50%> 

<img src="https://github.com/anet1223/Project/blob/main/plots/GO_processes.png" width=50% height=50%>
 

As the ribosomes were mostly involved in changes, a closer look was taken on ribosomal proteins: 23 ribosomal proteins with difference, 16 ribosomal proteins without changes

<img src="https://github.com/anet1223/Project/blob/main/plots/heatmap_ribosomal_proteins.png" width=50% height=50%>

## Conclusions

-   Expression of about half of the proteins is changes after R23I application (mainly decreases)

-   Expression profile is different in control and treatment groups. Small dose effect was observed.

-   Antibacterial action of R23I is likely through disturbance of translation, ion transport and oxidative stress response.

-   Different software (Peaks and IdentiPy) showed significant disagreement on the same data (common result about 60%)



