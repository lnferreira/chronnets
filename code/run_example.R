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

library(scales)
library(poweRlaw)
library(parallel)
library(magicaxis)
library(autoimage)

source("chronnets.R")
source("dataset_generator.R")

cat("Chronnets\n---------------------\n")
cat("Example of how to construct chronnets using artificial spatiotemporal data sets.\n\n")
devAskNewPage(ask = TRUE)
par(ask=TRUE, mar=c(0,0,0,0))

cat("1.1 Constructing a small data set (10x10) considering a power-law prability distribution...\n")
probs = matrix(scales::rescale(rpldis(n = 10*10, xmin = 1, alpha = 2.3), to=c(0.001,0.05)),10,10)
toy_ds = toy_experiment(p_matrix = probs, period = 500, num_cores = 2)
toy_ds_events = toy_experiment_generate_events(p_matrix = probs, dataset = toy_ds)
p = toy_experiment_plot_matrix(pnts = toy_ds_events$pnts, background = toy_ds_events$background) + 
    scale_fill_viridis(direction = -1, alpha = 0.6, option = "D")
show(p)

cat("1.2 Constructing the chronnet...\n")
chronnet = chronnet_create(dataset = toy_ds, self_loops = FALSE)
plot(chronnet, vertex.size=3, vertex.label=NA, edge.arrow.size=0.4, edge.width=E(chronnet)$weight)

cat("2.1 Constructing now a larger data set (100x100) also with a power-law prability distribution...\n")
probs = matrix(scales::rescale(rpldis(n = 100*100, xmin = 1, alpha = 2.3), to=c(0,0.02)),100,100)
toy_ds = toy_experiment(p_matrix = probs, period = 10000, num_cores = 2)
chronnet = chronnet_create(dataset = toy_ds, self_loops = FALSE)
plot(chronnet, vertex.size=3, vertex.label=NA, edge.arrow.size=0.4, edge.width=E(chronnet)$weight)

reset.par()
par(ask=TRUE)

cat("2.2 Plotting degree distribution...\n")
m_pl = displ$new(degree(chronnet))
m_pl$setXmin(estimate_xmin(m_pl))
plot(m_pl, axes=FALSE, xlab="k", ylab="P(k)", main="Degree distribution")
lines(m_pl, col=2, lwd=2)
magaxis(side=1:2)
box()

par(ask=FALSE)
