# Connectome Library 

##### Contains useful functions for future connectome work

##### Required modules
# Import all necessary libraries
import pandas as pd
import numpy as np
import scipy.stats as stats
import brainconn as con
import bct as bct
import networkx as nx

# A function to generate random graphs and calculate their mean Clustering Coefficient and mean Characteristic Path Length
#### Inputs: number of desired iterations, the number of vertices (nodes) and number of edges (number of connections) of the group averaged graph for the population cohort
#### Outputs: mean clustering coefficient (integer), mean path length (integer)
def random_graph_metrics(iterations, vertices, edges):
	clust_coef = []
	path_len = []
	for i in range(iterations):
		graph = con.reference.makerandCIJ_und(vertices, edges)
		#generate graph
	    #calculate characteristic path length
	    x = con.distance.charpath(graph)
		x = list(x)
	    path_len.append(x[0])
	    #calculate average clustering coefficient
	    y = con.clustering.clustering_coef_bu(graph)
	    y = np.mean(y)
	    clust_coef.append(y)
	    meanCC = np.mean(clust_coef)
	    meanPL = np.mean(path_len)
	return(meanCC, meanPL)

# A function to calculate the average number of edges in a population of graphs
### Inputs: a list of paths for the individual participant connectomes
### Outputs: average number of connections (edges)
def calculate_edges(list_of_paths):
	edges = []
	for i in range(len(List)):
		path = list_of_paths[i]
		data = pd.read_csv(path, header=None, sep=" ")
		x = np.count_nonzero(data.values)
		edges.append(x)
		num = (np.mean(np.array(edges))) / 2
	return(num)

# A function to identify network modules using Louvain community algorithm and consensus approach after 1000 iterations
### Inputs: the path to the average connectome for all participants, desired gamma (as g) and number of desired iterations for the Louvain algorithm
##### G values: small modules (<1) or large modules (>1), default = 1
### Outputs: connectivity matrix and a vector allocating each node to a module
def module_detection(path_averageConnectome, g, iterations):
	vectors = []
	qstats = []
	# Run Louvain algorithm
	for i in range(iterations):
		df = pd.read_csv(path_averageConnectome, header=None) 
		vector, qstat = con.modularity.community_louvain(df.values, gamma=g)
		vectors.append(vector)
		qstats.append(qstat)
	# Construct agreement matrix
	qstats = np.array(qstats)
	vectors = np.array(vectors)
	vectors = np.moveaxis(vectors, 0, -1)
	agr_matrix = con.clustering.agreement(vectors, buffsz=150)
	agr_matrix = con.utils.normalize(agr_matrix)
	# Consensus approach
	cons_matrix = con.clustering.consensus_und(agr_matrix, 0.5, 100)
	return(cons_matrix, vectors)


# A function to calculate the sum strength of connection between modules
### Inputs: graph (in pandas dataframe format), and 2 modules (each a list of the respective nodes allocated to each module)  
### Outputs: sum strength (integer)
def module_connections(graph, module1, module2):
	graph = pd.DataFrame(graph, index=module1)
	graph = graph[module2]
	sumStrength = graph.values.sum()
	return sumStrength

# A function to extract a list of total values from a subselection of a matrix
### Inputs: graph (as a pandas dataframe) and 2 modules
### Outputs: list of connection values from the subselection of the graph
def module_values(graph, module1, module2): 
	graph = pd.DataFrame(graph, index=module1)
	graph = graph[module2]
	values = graph.values
	values = values.flatten()
	return values

# A function to extract average characteristic path length from a subsection of a matrix
### Inputs a graph 
def module_length(graph, module1, module2):
#     module1[:] = [x - 1 for x in module1]
#     module2[:] = [x - 1 for x in module2]
    graph = pd.DataFrame(graph, index=module1)
    graph = graph[module2]
    meanLength = graph.values.sum() / (len(module1, module2))
    return meanLength


# A function to z and tahn transform the values of a connection to atrophy values.
### Inputs: connectivity strength connectome (per individual),  a control average and control std connectome. All connectomes need to be in pandas dataframe format.
### Outputs: connectivity matrix with tahn transformed values inplace of connectivity strength
def tahn_transform(graph, average_connectome, std_connectome):
	data = np.genfromtxt(graph)
	data_z = (average_connectome - data) / std_connectome
	data_transform = np.tahn(data_z)
	return data_transform

# Loading and extracting data from a matlab file containing gene expression data from Allen atlas in Glasser space (extracted using matlab code from Arnatkevic̆i et al. Neuroimage, 2019)
### Inputs: the path to the matlab file
### Outputs: gene list, gene symbol list, 
def annotation(Path):
	annot = sio.loadmat(Path)
	names = annot["probeInformation"] #find probe information on dictionary 
	df = pd.DataFrame(data=names.flatten()) #load as dataframe
	geneList = df["EntrezID"] #find gene symbol list
	geneList= geneList[0].tolist() # transform dataframe column to list
	geneSymbol= df["GeneSymbol"]
	geneSymbol = geneSymbol[0].tolist()
	return geneList, geneSymbol