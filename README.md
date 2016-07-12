# neuroblastoma-graph-database

##Summary
---
This project is focused on building a graph database for the study of pediatric cancer, specifically neuroblastoma. Neuroblastoma is an embryonal tumor of the autonomic nervous system and typically occurs in very young children. It accounts for disproportionate death rate among pediatric cancers and is associated with the highest proportions of regression in all human cancers. Data from pediatric cancers are complex as they try to capture the relationships between genes, pathways, molecular functions, and drug response. Common database structures, such as relationship databases, do not allow for efficient capture of the complicated relationships that exist in the data. Graph databases take advantage of these connections and store data in the form of nodes and relationships which can each contain any number of properties. Graph databases excel at managing highly connected data and complex queries because they prioritize relationships and connections don’t have to be inferred cutting down computation time. We want to capture neuroblastoma cell line data in a Neo4J database and use this infrastructure to discover novel findings. Neo4J is an open source NoSQL graph database that uses Cypher as its open graph query language. The database will be accessed through R, which includes packages that allow for effective visualization and manipulation of the data. The neuroblastoma data contains information on cell lines, genes, pathways, and treatments which will be modeled as nodes while information on treatment response and relevant mutations will be modeled as relationships between the respective nodes. We hope that this method will allow us to better capture the connections present in the data and consequently allow us to query for interesting patterns and relationships.

#Workflow
---
- Learned how to navigate Neo4j and Cypher (Neo4j's query language) along with some graph theory
- Obtained treatment data and familiarized myself with the research
- Formatted and loaded cell line mutations, pathways, expression outliers data into the database (still need to load in copy number data)
- Proof of concept with MEK inhibitors 
