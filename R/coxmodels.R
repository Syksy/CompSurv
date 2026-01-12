# Functions related to Cox proportional hazards models

#' Fit multivariable Cox ph model for all covariates together with univariable model once per each covariate
#' @export
coxambi <- function(obj, ...) {
	UseMethod("coxambi")
}

#' Not a suitable class provided as first parameter for coxambi
#' @export
coxambi.default <- function(obj, ...) {
	stop("Unsupported input type for 'coxambi()'. Provide a formula or a survival::coxph object.")
}

#' Wrapper for a formula + data formulation of Cox model
#' @export
coxambi.formula <- function(
	# Left hands side (lhs) is the survival response, right hand side (rhs) are the covariates
	obj, 
	# Data used for fitting
	data, 
	# Parameters passed on to both the coxph call in multivariable and univariable models
	... 
){
	formula <- obj
	if(!inherits(formula, "formula")){
		stop("First parameter should be of 'formula' class")
	}
	if (missing(data)) stop("'data' must be provided.")
	
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

#' Wrapper for a readily fitted Cox model and perform single-variable fits per each fitted covariate
#' @export
coxambi.coxph <- function(
	# Fitted coxph-object
	obj,
	...
){
	if (is.null(obj$model)) {
		stop("Need 'obj$model' for re-fitting. Refit with survival::coxph(..., model = TRUE).")
	}
	combs <- coxambi(obj$formula, data = obj$model)
	combs
}

#' Fetch summary of Cox model coefficients in data.frame class format
#' @export
coxdf <- function(
	# A Cox PH model object
	model
){
	# Extract summary
	summ <- summary(model)
	df <- data.frame(
		term  = rownames(summ$coef),
		coef  = unname(s$coef[,"coef"]),
		hr    = unname(s$coef[,"exp(coef)"]),
		se    = unname(s$coef[,"se(coef)"]),
		z     = unname(s$coef[,"z"]),
		p     = unname(s$coef[,"Pr(>|z|)"]),
		dw95  = unname(s$coef[,"lower .95"]),
		up95  = unname(s$coef[,"upper .95"]),
		stringsAsFactors = FALSE
	)
	df
}
