# Diagnostic use functions mainly intended for survival models

#' Test for potential collinearity issues in a multivariable Cox model via paired comparisons
#' @export
collin_pairtest <- function(
	# coxph-object
	fit,
	# How many digits to report in p.values
	digits = 4,
	...
){
	if(!inherits(fit, "coxph")){
		stop("Provided fit should inherit class 'coxph' for a Cox proportional hazards model")
	}

	# Data frame, omitting LHS (survival response column)
	df <- model.frame(fit)[,-1,drop=FALSE]
	# Types of variables
	# Aggregate classes of variables to main types
	types <- lapply(df, FUN=\(x){
		if(inherits(x, c("numeric", "integer"))) return("numeric")
		else if(inherits(x, c("factor", "character", "logical"))) return("categorical")
		else return("unknown")
	})

	# Matrix of p-values between testing for collinearity between two variables provided for a Cox model
	pmat <- matrix(NA, nrow=length(types), ncol=length(types))
	# Name rows/cols
	rownames(pmat) <- colnames(pmat) <- names(types)
	
	# First iterator
	for(i in 1:(ncol(df)-1)){
		# Second iterator
		for(j in (i+1):ncol(df)){
			# Values
			val_i <- df[[i]]
			val_j <- df[[j]]
			# Types
			type_i <- types[[i]]
			type_j <- types[[j]]
			
			# Both variables are categorical, use Fisher's Exact Test
			if(type_i == "categorical" && type_j == "categorical"){
				pval <- fisher.test(val_i, val_j)$p.value
			# Both variables are numeric, use Spearman correlation
			}else if(type_i == "numeric" && type_j == "numeric"){
				pval <- cor.test(val_i, val_j, method="spearman", exact=FALSE)$p.value
			# i:th is numeric and j:th is categorical or vice versa, use Kruskal-Wallis
			}else{
				# i:th numeric and j:th categorical
				if(type_i == "numeric"){
					pval <- kruskal.test(val_i ~ as.factor(val_j))$p.value
				# i:th categorical and j:th numeric
				}else{
					pval <- kruskal.test(val_j ~ as.factor(val_i))$p.value
				}
			}
			pmat[j,i] <- pmat[i,j] <- round(pval, digits)
		}
	}	
	# Set unit diagonal and return p-value matrix
	diag(pmat) <- 1
	pmat
}

#' Eigenvalue based collinearity inspection
#' @export
collin_eigen <- function(
	# coxph-object
	fit,
	# Eigenvalue thresholds, arbitrary choice
	threshold = 30,
	# Digits for loadings
	digits = 3,
	...
){
	if(!inherits(fit, "coxph")){
		stop("Provided fit should inherit class 'coxph' for a Cox proportional hazards model")
	}

	X <- model.matrix(fit)
	# Scale to avoid large values dominating
	Xs <- scale(X)
	
	R <- cor(Xs)
	eig <- eigen(R)
	cond_index <- max(eig$values) / eig$values
	# Pick dimensions which have a condition above the threshold
	bad_dims <- which(cond_index > threshold)   
	loadings <- eig$vectors[,bad_dims,drop = FALSE]
	rownames(loadings) <- colnames(X)
	colnames(loadings) <- paste0("dim", bad_dims)
	# Report variable loadings for the problematic dimensions
	round(loadings, digits)
}

#' VIF for Cox model collinearity inspection
#' @export
collin_vif <- function(
	# coxph-object
	fit,
	...
){
	if(!inherits(fit, "coxph")){
		stop("Provided fit should inherit class 'coxph' for a Cox proportional hazards model")
	}

	# TODO
}