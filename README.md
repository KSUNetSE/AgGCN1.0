# AgGCN1.0
Data and code for Anopheles gambiae gene co-expression network, link of the paper: https://www.biorxiv.org/content/10.1101/2022.01.03.474847v1.abstract  

1)	Network related data are stored in the folder Network_Data, containing the following four files:
   
   a)	GAP_GCN1.0.gephi, visualization of the network via gephi (https://gephi.org/). By choosing the attribute in overview/nodes/partition, the viewer can see the communities, core, and the four different node centralities described in the manuscript.
   
   b)	AgGCN1.0_edge_list.csv, edge list of the constructed network, which allows one to make their own analysis.
   
   c)	AgGCN1.0_properties.csv, node properties, including AGAP number, strength, eigenvector, closeness, betweenness, community membership, and presence in cores.
   
   d)	communities.pdf, a visualization (through gephi) of the community distribution of the genes.

2)	Visualization of GO term enrichment data can be found in SupplementalFigures: 
   
   a)	Each figure summarizes the enriched GO terms and KEGG pathways per community and core.  For example, KUANGetal_FigS4_GOplot_community1.pdf shows the GO terms that are enriched in community 1.
   
3)	Raw data and analysis data can be found in SupplementalTables:
   
   a)	KUANGetal_TableS1_Anopheles-gambiae_EXPR-STATS_VB-2019-02_final.xlsx, the raw data used for the network construction, including platform type (microarray, RNA-seq, etc), platform name (EMBL A.gambiae MMC1 20k v1.0, Affymetrix GeneChip® Plasmodium/Anopheles Genome Array), platform identifier (GEO/ArrayExpress), Experiment identifier (ArrayExpress/GEO/bioproject), publications (authors), life stage, tissue, sex, physiology, gene ID, conditions.
   
   b)	KUANGetal_TableS2_GOenrichmentCores.xlsx, a full list of the GO terms that are enriched in the cores.
   
   c)	KUANGetal_TableS3_GOKEGGenrichmentCommunities.xlsx, a full list of GO/KEGG IDs that are enriched in different communities. The spreadsheet also contains a list of enriched genes for the corresponding GO/KEGG ID.
   
   d)	KUANGetal_TableS4_AgGCN1.0_Dataset_conditioninfo.xlsx, the information regarding the conditions we used to construct the network, including the condition number (which can be found in KUANGetal_TableS1_Anopheles-gambiae_EXPR-STATS_VB-2019-02_final.xlsx), platform type, platform name, platform identifier, experiment identifier, number of genes probed, authors, tissue, physiology, condition, etc.

4)	Code
   
   a)	A. PCC_edge_list.py, code to generate edge list.
   
   b)	B. s-core.py, code to find the core of the network.
   
   c)	C. sliding_threshold_fitting.py, find the sliding thresholds. 
