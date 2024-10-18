# Proportions plot
#' @export propplot
propplot <- function(
	x, # Vector of event types
	times, # Event times (until censoring or non-censored event)
	xs = seq(0, 10, length.out=1000), # Time points in same unit	
	leftover = "Alive", # Vector of names of variable that ought to be treated as 1-p instead of p
	xlim = range(xs),
	ylim = c(0,1),
	cols = rainbow(length(unique(x))+1),
	...
){
	polys <- polygonprop(x = x, times = times, xs = xs, leftover = leftover)
	
	plot.new()
	plot.window(xlim=xlim, ylim=ylim)

	for(i in 1:length(polys)){
		polygon(x = c(xs, rev(xs)), c(polys[[i]][1,], rev(polys[[i]][2,])), col = cols[i], ...)
	}

	box(); axis(1); axis(2)

	legend <- cols
	names(legend) <- names(polys)

	# Invisibly return the legend and propmat	
	invisible(legend)	
}
