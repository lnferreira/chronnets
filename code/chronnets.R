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

library("igraph")

# =================================================================================
# Creates a network from a temporal data set. This version accepts multiple events
#   in the same period of time.
# =================================================================================
chronnet <- function(dataset, self_loops=TRUE, mode="directed",  num_cores=2) {
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

# =================================================================================
# Spatiotemporal clustering
# =================================================================================
spatio_temporal_clustering <- function(data_set, cut_at, min_pnts_change=2){
    net = chronological_net_v2(data_set, self_loops = F, mode="undirected", num_cores = 2)
    comm = fastgreedy.community(net)
    if (!missing(cut_at)) {
        comm = cutat(comm, cut_at)
    } else {
        comm = membership(comm)
    }
    comm = data.frame(cell=as.numeric(V(net)$name), cluster=as.numeric(comm))
    data_set = left_join(data_set, comm, by="cell")
    clustering = data_set$cluster
    max_window_size = 2 + min_pnts_change
    for (window_size in 3:(max_window_size)){
        for (w_index in window_size:length(clustering)) {
            window_indices = (w_index-window_size+1):w_index
            wd = clustering[window_indices]
            middle_wd = wd[c(-1, -length(wd))]
            if (tail(wd,1) == head(wd,1) & any(middle_wd != tail(wd,1)))
                clustering[window_indices[-c(1,window_size)]] = tail(wd,1)
        }
    }
    for (w_index in 3:length(clustering)) {
        window_indices = (w_index-2):w_index
        if (length(unique(clustering[window_indices])) == 3)
            clustering[window_indices[2]] = clustering[window_indices[1]]
    }
    clustering
}

# =================================================================================
# Generate spatiotemporal datasets 
# Grid ids: 
# 
#   1 | 2 | 3 
#   ---------
#   4 | 5 | 6
#   ---------
#   7 | 8 | 9
# =================================================================================
toy_experiment <- function(p_matrix, period=100, one_cell_per_time=FALSE, num_cores=1){
    temporal_grid_states = mclapply(1:period, function(t){
        grid_state = matrix(runif(prod(dim(p_matrix))), nrow = nrow(p_matrix), ncol= ncol(p_matrix), byrow = T) 
        active_cells = which(grid_state < p_matrix, arr.ind = T)
        if (length(active_cells)){
            active_cells = data.frame(t=t, x=active_cells[,1], y=active_cells[,2])
        } else {
            active_cells = data.frame(t=c(), x=c(), y=c())
        }
        active_cells
    }, mc.cores = num_cores)
    if (all(sapply(temporal_grid_states, nrow) == 0)){
        warning("Empty data set generated")
        return(c())
    }
    temporal_grid_states = temporal_grid_states[sapply(temporal_grid_states, nrow) > 0]
    temporal_grid_states = do.call(rbind, temporal_grid_states)
    temporal_grid_states$cell = apply(temporal_grid_states[,-1], 1, function(cell) ((cell[1] - 1) * ncol(p_matrix)) + cell[2])
    if(one_cell_per_time)
        temporal_grid_states$t = 1:nrow(temporal_grid_states)
    rownames(temporal_grid_states) = NULL
    temporal_grid_states
}