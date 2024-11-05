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
	Ncount = "left",
	...
){
	# Named color vector provided, will use it
	if(!is.null(names(cols))){
		# Compute polygons with (reversed) color ordering, as the direction is from bottom up
		polys <- CompSurv:::polygonprop(x = x, times = times, xs = xs, leftover = leftover, rows = names(cols))
		# Match names
		cols <- cols[match(names(polys), names(cols))]
	}else{
		polys <- CompSurv:::polygonprop(x = x, times = times, xs = xs, leftover = leftover)
	}
	
	# Setup basic base plot window
	plot.new()
	plot.window(xlim=xlim, ylim=ylim)
	
	# Draw proportion polygons
	for(i in 1:length(polys)){
		polygon(x = c(xs, rev(xs)), c(polys[[i]][1,], rev(polys[[i]][2,])), col = cols[i], ...)
	}

	# Surrounding box and axes
	box(); axis(1); axis(2)

	# If non-NA, add a legend for N count
	legend(Ncount, paste0("N=", length(x)))

	# Return invisibly the color legend
	legend <- cols
	#names(legend) <- c(names(polys), leftover)
	names(legend) <- names(cols)

	# Invisibly return the legend and proportion polygons	
	invisible(list(legend = legend, polys = polys))	
}
