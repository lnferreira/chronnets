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

list_packages = c("ggplot2", "viridis")
new.packages = list_packages[!(list_packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) 
    install.packages(new.packages)
for (package in list_packages)
    library(package, character.only = T)

# =================================================================================
# Generate spatiotemporal datasets 
# Return: a data frame with time (t), x, and y in the grid where an event was 
#   generated and the cell ID. Please note that x and y are not the actual 
#   positions of the event, but just the row and column where the event occurred in 
#   the grid system used to represent the modeled area. The actual event position 
#   is not required for the chronnet construction. This information is not 
#   necessary for the experiments using the proposed data set generator. If the 
#   actual position is necessary, check the function toy_experiment_generate_events.
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

# =================================================================================
# Generate events for a specific toy experiment (function toy_experiment)
# Return: Given an artificial data set, it returns the actual (x,y) position for  
#         the events considering that each grid cell has 10 x 10 (area). This 
#         utility function is here used just for plotting purposes. This function 
#         might also be useful to use as input for other spatiotemporal methods.
# =================================================================================
toy_experiment_generate_events <- function(p_matrix, dataset=c(), offset=1){
    min_points = expand.grid(x_min=seq(0,90,10), y_min = rev(seq(0,90, 10)))
    max_points = do.call(rbind, apply(min_points, 1, function(pnt) data.frame(x_max=pnt[1]+10, y_max=pnt[2]+10)))
    rownames(max_points) = NULL
    background = cbind(min_points, max_points)
    background$p = as.vector(t(p_matrix))
    pnts=data.frame(x=c(),y=c())
    if (length(dataset)) {
        pnts = do.call(rbind, lapply(dataset$cell, function(cell){
            border = background[cell, ]
            x = sample((border$x_min+offset):(border$x_max-offset), 1)
            y = sample((border$y_min+offset):(border$y_max-offset), 1)
            data.frame(x=x,y=y)
        }))
    }
    list(pnts=pnts, background=background)
}

# =================================================================================
# Plot the data sets (return of function toy_experiment_generate_events)
# =================================================================================
toy_experiment_plot_matrix <- function(pnts, background, pnts_color_col = "", 
                                       point_size=0.9, legendbar_length=10, show=TRUE) {
    if (pnts_color_col == ""){
        p = ggplot(pnts, aes(x=x,y=y))
    } else {
        p = ggplot(pnts, aes_string(x="x", y="y", color=pnts_color_col))
    }
    p = p +     
        geom_rect(data=background, mapping = aes(x = NULL,y = NULL, 
                                                 xmin=x_min,xmax=x_max, 
                                                 ymin=y_min, ymax=y_max, 
                                                 fill=p, color=NULL)) +
        geom_vline(xintercept = seq(10,90,10), color="white") + 
        geom_hline(yintercept = seq(10,90,10), color="white") + 
        geom_point(size=point_size) +
        scale_x_continuous(breaks = seq(0, 100, 10), limits = c(0,100), expand=c(0,0)) +
        scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0,100), expand=c(0,0)) + 
        theme(
            axis.title = element_blank(), 
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "bottom",
            legend.box.margin=margin(-8,-10,-10,-10), 
            plot.margin = margin(0.1, 0.1, 0.2, 1, "cm"),
            panel.grid.major = element_line(color = "gray80", size = 0.9)) + 
            guides(fill = guide_colorbar(barwidth = legendbar_length, barheight = 1, title="p", frame.colour=c("black"), 
                                     ticks.linewidth = 1.5, nbin=10, raster=T, frame.linewidth=1.2, 
                                     ticks.colour="black", direction="horizontal", title.vjust = 1,
                                     title.theme=element_text(family="Times", face="italic", size=16, colour = "black"))) 
}