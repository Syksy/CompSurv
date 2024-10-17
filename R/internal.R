# Proportions of survival across multiple event types as a function of time
props <- function(
	x, # Vector of event types
	times, # Event times (until censoring or non-censored event)
	xs = seq(0, 10, length.out=1000), # Time points in same unit
	...
){	
	# Make sure to use x as a factor
	x <- factor(x)
	# Aggregate proportion across different factor levels as a function of time
	props <- by(data.frame(x, times), INDICES=x, FUN=\(z){
		z <- z[order(z[,2], decreasing=FALSE),]
		z <- unlist(lapply(xs, FUN=\(q){
			sum(z[,2] <= q)/length(x)
		}))
		z
	})
	# Name by-list elements and return
	names(props) <- levels(x)
	props
}

# Stack proportions and create last class of the remaining proportion as a matrix
stackprop <- function(
	x, # Vector of event types
	times, # Event times (until censoring or non-censored event)
	xs = seq(0, 10, length.out=1000), # Time points in same unit	
	leftover = "Alive", # Vector of names of variables that ought to be treated as 1-p instead of p
	...
){
	# Calculate proportions
	props <- props(x=x, times=times, xs=xs)

	# Make sure to use x as a factor
	x <- factor(x)

	propmat <- matrix(nrow=length(levels(x))+1, ncol=length(xs))
	rownames(propmat) <- c(levels(x), leftover)
	for(r in 1:(nrow(propmat)-1)){
		propmat[r,] <- unlist(props[[r]])
	}
	propmat[nrow(propmat),] <- 1 - apply(propmat[-nrow(propmat),], MARGIN=2, FUN=sum)
	propmat
}
