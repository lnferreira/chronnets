# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

list_packages = c("igraph", "dplyr", "plyr", "zoo")
new.packages = list_packages[!(list_packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) 
    install.packages(new.packages)
for (package in list_packages)
    library(package, character.only = T)

# =================================================================================
# Creates a chronnet from a temporal data set. This method accepts multiple events
#   in the same period of time. Dataset should be a data frame with two columns:
#       cell: The id of the region where an event occured
#       t: time
# =================================================================================
chronnet_create <- function(dataset, self_loops=TRUE, mode="directed",  num_cores=2) {
    time_seq = sort(unique(dataset$t))
    if (length(time_seq) < 2)
        stop("The total time interval in the dataset should be larger than two.")
    max_num_cores = num_cores:1
    num_cores = max(max_num_cores[which((max_num_cores*2) <= length(time_seq))])
    time_list = split(time_seq, ceiling(seq_along(time_seq)/ceiling(length(time_seq)/num_cores)))
    links = mclapply(time_list, function(time_list_part){
        time_combs = rollapply(time_list_part, 2, function(x) x)
        connections = apply(time_combs, 1, function(x){
            cells_before = dataset[dataset$t == x[1], "cell"]
            cells_after = dataset[dataset$t == x[2], "cell"]
            links_part = as.data.frame(expand.grid(cells_before, cells_after))
            links_part$weight = 1
            links_part
        })
        connections = do.call(rbind, connections)
        names(connections) = c("from", "to", "weight")
        connections
    }, mc.cores = num_cores)
    if (length(time_list) > 1){
        for (i in 2:length(time_list)){
            time1 = tail(time_list[[i-1]], 1)
            time2 = head(time_list[[i]], 1)
            cells_before = dataset[dataset$t == time1, "cell"]
            cells_after = dataset[dataset$t == time2, "cell"]
            connections = as.data.frame(expand.grid(cells_before, cells_after))
            connections$weight = 1
            names(connections) = c("from", "to", "weight")
            links[[i=1]] = rbind(links[[i=1]], connections)
        }
    }
    links = do.call(rbind, links)
    links = links %>% group_by(from, to) %>% tally()
    names(links) = mapvalues(names(links), "n", "weight")
    net = graph.empty()
    if (nrow(links)){
        net = graph_from_data_frame(links, directed = TRUE)
        if (!self_loops)
            net = simplify(net, remove.multiple = FALSE)
        if (mode == "undirected")
            net = as.undirected(net, edge.attr.comb = "sum")
    } else {
        warning("Empty graph returned")
    } 
    net
}

# =================================================================================
# Prune network
# =================================================================================
prune_net <- function(net, prune_edge_count, undirected=FALSE){
    if (undirected)
        net = as.undirected(net, edge.attr.comb = "sum")
    net = delete.edges(net, which(E(net)$weight <= prune_edge_count))
    net_list = decompose(net)
    net_list[[which.max(sapply(net_list, vcount))]]
}