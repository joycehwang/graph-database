# neuroblastoma-graph-database

##Summary
---
This project is focused on building a graph database for the study of pediatric cancer, specifically neuroblastoma. Neuroblastoma is an embryonal tumor of the autonomic nervous system and typically occurs in very young children. It accounts for disproportionate death rate among pediatric cancers and is associated with the highest proportions of regression in all human cancers. Data from pediatric cancers are complex as they try to capture the relationships between genes, pathways, molecular functions, and drug response. Common database structures, such as relationship databases, do not allow for efficient capture of the complicated relationships that exist in the data. Graph databases take advantage of these connections and store data in the form of nodes and relationships which can each contain any number of properties. Graph databases excel at managing highly connected data and complex queries because they prioritize relationships and connections don’t have to be inferred cutting down computation time. We want to capture neuroblastoma cell line data in a Neo4J database and use this infrastructure to discover novel findings. Neo4J is an open source NoSQL graph database that uses Cypher as its open graph query language. The database will be accessed through R, which includes packages that allow for effective visualization and manipulation of the data. The neuroblastoma data contains information on cell lines, genes, pathways, and treatments which will be modeled as nodes while information on treatment response and relevant mutations will be modeled as relationships between the respective nodes. We hope that this method will allow us to better capture the connections present in the data and consequently allow us to query for interesting patterns and relationships.

##Workflow
---
- Learned how to navigate Neo4j and Cypher (Neo4j's query language) along with some graph theory
  - [Online Neo4j tutorial] (https://neo4j.com/graphacademy/online-training/introduction-graph-databases/)
  - Used community version 3.0.1 (Linux/UNIX) and browser interface for quick queries and visualizations
  - On start-up:
    1. Open new terminal with neo4j-commjunity-3.0.1 as directory
    2. Start neo4j server: ./bin/neo4j start
    3. Open shell: bin/neo4j-shell
    4. To stop: ./bin/neo4j stop
  - Place all files to read into the database in the import folder
- Using R as driver
  - [[Video] Neo4j, Graphs R Cool] (https://youtu.be/bdQ90y9Pefo)
  - [[Video] Visualizations with RNeo4j] (https://youtu.be/5u4eT1OgB88)
  - [[Article] Visualize Your Graph with RNeo4j and visNetwork] (https://nicolewhite.github.io/2015/06/18/visualize-your-graph-with-rneo4j-and-visNetwork.html)
  - [[Article] Visualizing Your Graph with RNeo4j] (https://neo4j.com/blog/visualize-graph-with-rneo4j/)
- Obtained treatment information, familiarized myself with the research, and started data modeling
  - [[Paper] Recent advances in neuroblastoma (Review)] (http://www.ncbi.nlm.nih.gov/pubmed/20558371)
  - [[Paper] Cancer genes and the pathways they control] (http://www.ncbi.nlm.nih.gov/pubmed/15286780)
  - [Arrow tool for data modeling] (http://www.apcjones.com/arrows/)
- Formatted and loaded cell line mutations, pathways, expression outliers data into the database (still need to load in copy number data)
  - Mutations(/Amplificatiton/Deletion): raw data in the form of matrix with cell lines and genes as axes and value of 0.5 indicating mutation and value of 0 indicating none -> csv with column of cell line name and column mutated gene name
    - Value of 1 indicates deletion and value of 2 indicated amplification
  - Expression Outliers: raw data in the form of a matrix with cell lines and genes as axes and value of 1 indicating outlier -> csv with column of cell line name and column of outlier gene
  - Pathway: raw data with columns of pathway name, url, list of associated genes -> csv with column of pathway name and column of associated gene
  - [[Data Set] Pathway data] (http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2)
  - [[Article] Data modelling] (https://linkurio.us/the-crunchbase-graph-data-modelling/)
  - [[Article] Importing data into Neo4j] (https://linkurio.us/crunchbase-graph-importing-data-neo4j/)
  - [[Article] Some CSV import tricks in Neo4j] (http://blog.comperiosearch.com/blog/2015/02/04/csv-import-tricks-neo4j/)
- Proof of concept with MEK1/2 inhibitors 
  - [[Paper] Relapsed neuroblastomas show frequent RAS-MAPK pathway mutations] (http://www.ncbi.nlm.nih.gov/pubmed/26121087)
  - [[Paper] MEK1 and MEK2 inhibitors and cancer therapy the long and winding road] (http://www.ncbi.nlm.nih.gov/pubmed/26399658)
  - [[Paper] The clinical development of MEK inhibitors] (http://www.nature.com/nrclinonc/journal/v11/n7/full/nrclinonc.2014.83.html)
  - Look at differential between MEK sensitive and resistant cell lines for gene and pathway counts
  - There should be more connections from sensitive cell lines to targeted genes and less connections from resistant cell lines to target genes -> lower resisitant:sensitive ratio
  - Set IC50 cut off scores as 1000 nM
  - [RACCYCD] (http://software.broadinstitute.org/gsea/msigdb/cards/BIOCARTA_RACCYCD_PATHWAY) pathway yieled significant results
- Proof of concept with CDK4/6 inhibitors
  - [[Paper] Dual CDK4/CDK6 Inhibition Induces Cell Cycle Arrest and Senescence in Neuroblastoma] (http://www.ncbi.nlm.nih.gov/pubmed/24045179)
  - [[Paper] Inhibition of CDK4/6 as a novel therapeutic option for neuroblastoma] (https://cancerci.biomedcentral.com/articles/10.1186/s12935-015-0224-y)
  - Look at differential between CDK sensitive and resistant cell lines for gene and pathway counts
  - Set IC50 cut off scores as 1000 nM
  - [MYCN] (http://www.ncbi.nlm.nih.gov/gene/4613) yieled significant results
- Proof of concept with BET inhibitors
  - [[Paper] Targeting MYCN in neuroblastoma by BET bromodomain inhibition] (http://cancerdiscovery.aacrjournals.org/content/3/3/308.long)
  - [[Abstract] Antitumor activity and sensitivity evaluation of novel BET inhibitors in neuroblastoma] (http://mcr.aacrjournals.org/content/13/10_Supplement/B34)
  - Look at differential between BET sensitive and resistant cell lines for gene and pathway counts
  - Set IC50 cut off scores as 500 nM
