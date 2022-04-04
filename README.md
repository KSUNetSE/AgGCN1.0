# AgGCN1.0
Data and code for Anopheles gambiae gene co-expression network, link of the paper: https://www.biorxiv.org/content/10.1101/2022.01.03.474847v1.abstract  

1. Network related data are stored in the folder Network_Data, the three files in the folder are respectively:

   A. AGAP_GCN1.0.gephi, visualization of the network via gephi (https://gephi.org/). By choosing the atrribute in overview/nodes/partition, we can observe the community distribution, between centrality, eigenvector centrality, etc.
   
   B. AgGCN1.0_edge_list.csv, edge list of the constructed network, which allows one to make their own analysis.
   
   C. AgGCN1.0_properties.csv, properties of the genes, including strength, eigenvector, close, betweenness, community, cores or not.
   
   D. communities.pdf, a visualization (thorugh gephi) of the community distribution of the genes.
   
2. visualization og GO term enrichment data can be found in SupplementalFigures: 

   A. The figures show that the communities are enriched for different GO terms. For example, KUANGetal_FigS4_GOplot_community1.pdf shows the GO terms that are enriched in community 1.
  
3. raw data and analysis data can be found in SupplementalTables:

   A. KUANGetal_TableS1_Anopheles-gambiae_EXPR-STATS_VB-2019-02_final.xlsx, the raw data used for the network construction, including platform type (microarray, RNA-seq, etc), platform name (EMBL A.gambiae MMC1 20k v1.0, Affymetrix GeneChip® Plasmodium/Anopheles Genome Array), platform identifier (GEO/ArrayExpress), Experiment identifier (ArrayExpress/GEO/bioproject), publications (authors), life stage, tissue, sex, physiology, gene ID, conditions.
   
   B. KUANGetal_TableS2_GOenrichmentCores.xlsx, a full list of the GO terms that are enriched in the cores.
   
   C. KUANGetal_TableS3_GOKEGGenrichmentCommunities.xlsx, a full list of GO/KEGG IDs that are enriched in different communities. The spreadsheet also contains a list of enriched genes for the corresponding GO term.
   
   D. KUANGetal_TableS4_AgGCN1.0_Dataset_conditioninfo.xlsx, the information regarding the conditions we used to construct the network, including the condition number (which can be found in KUANGetal_TableS1_Anopheles-gambiae_EXPR-STATS_VB-2019-02_final.xlsx), platform type, platform name, platform identifier, experiment indetifier, number of genes probed, authors, tisse, physiology, condition, etc.
     
