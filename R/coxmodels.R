# Functions related to Cox proportional hazards models

#' Fit multivariable Cox ph model for all covariates together with univariable model once per each covariate
#' @export coxambi
coxambi <- function(
	# Left hands side (lhs) is the survival response, right hand side (rhs) are the covariates
	formula, 
	# Data used for fitting
	data, 
	# Parameters passed on to both the coxph call in multivariable and univariable models
	... 
){
	if(!inherits(formula, "formula")){
		stop("First parameter should be of 'formula' class")
	}
	
	# Fit the multivariable model
	coxmulti <- survival::coxph(formula, data = data, ...)
	print(coxmulti)
	
	# Fit each univariable model separately
	# Parse out the rhs
	rhs <- attr(terms(formula, data = data), "term.labels")
	# Parse out the lhs
	lhs <- formula[[2]]
	print(rhs)
	print(lhs)
	
	# lapply across all right hand side term candidates
	coxunis <- lapply(rhs, FUN=\(x){
		# Make the formula
		tmp_formula <- stats::as.formula(paste(deparse(lhs), "~", x))
		survival::coxph(tmp_formula, data = data, ...)
	})
	
	# Return the multivariable fit and univariable fits
	list(
		coxmulti,
		coxunis
	)
}
