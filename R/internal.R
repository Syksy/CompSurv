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
	leftover = "Alive", # Vector of names of variable that ought to be treated as 1-p instead of p
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

# Create a list of properly aligned polygons in y in [0,1] and x in [xs]
polygonprop <- function(
	x, # Vector of event types
	times, # Event times (until censoring or non-censored event)
	xs = seq(0, 10, length.out=1000), # Time points in same unit	
	leftover = "Alive", # Vector of names of variable that ought to be treated as 1-p instead of p
	...
){

	# Create stacked proportions matrix
	stackmat <- stackprop(x = x, times = times, xs = xs, leftover = leftover)

	# Create a "zero" line
	stackmat <- rbind(0, stackmat)

	polys <- list()
	for(r in 2:nrow(stackmat)){
		polys[[r-1]] <- rbind(
			# Upper bound
			apply(stackmat[seq(from=1, to=(r-1), by=1),,drop=FALSE], MARGIN=2, FUN=sum),
			# Lower bound
			apply(stackmat[seq(from=1, to=r, by=1),,drop=FALSE], MARGIN=2, FUN=sum)
		)
	}

	# Drop the zero row name and name polygons accordingly
	names(polys) <- rownames(stackmat)[-1]
	# Attach the x-coordinates as an attribute
	attr(polys, "xs") <- xs
	
	polys
}

