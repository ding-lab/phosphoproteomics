# Basic setup
library(igraph)
library(RJSONIO)
library(httr)

port.number = 1234
base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")
print(base.url)


# Loading Networks --------------------------------------------------------
# Load list of edges as Data Frame
network.df <- MAPK_trans

# Convert it into igraph object
network <- graph.data.frame(network.df,directed=T)

# Remove duplicate edges & loops
g.tca <- network

# Name it
g.tca$name = "MAPK_trans"

# Convert igraph object into JSON -----------------------------------------

# This function will be published as a part of utility package, but not ready yet.
source('../utility/cytoscape_util.R')

# Convert it into Cytosccape.js JSON
cygraph <- toCytoscape(g.tca)

send2cy(cygraph, "default","circular")


