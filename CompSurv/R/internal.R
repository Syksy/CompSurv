props <- function(
	x, # Event types
	times, # Event times (until censoring or non-censored event)
	xs = seq(0, 10, length.out=1000),
	col = rainbow(length(unique(x))),
	...
){	
	# Make sure to use x as a factor
	x <- factor(x)
	
	props <- by(cbind(x, times), INDICES=x, FUN=\(z){
			z <- z[order(z[,2], decreasing=FALSE),]
			unlist(lapply(xs, FUN=\(q){
				sum(z[,2] <= q)/nrow(x)
		}))
	})
	names(props) <- levels(x)
	props
}
