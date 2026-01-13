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
	coxmulti <- survival::coxph(formula, data = data, model = TRUE, ...)

	# Fit each univariable model separately
	# Parse out the rhs
	rhs <- attr(terms(formula, data = data), "term.labels")
	# Parse out the lhs
	lhs <- formula[[2]]

	# lapply across all right hand side term candidates
	coxunis <- lapply(rhs, FUN=\(x){
		# Make the formula
		tmp_formula <- stats::as.formula(paste(deparse(lhs), "~", x))
		survival::coxph(tmp_formula, data = data, model = TRUE, ...)
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
	combs <- coxambi.formula(obj$formula, data = obj$model)
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
		coef  = unname(summ$coef[,"coef"]),
		hr    = unname(summ$coef[,"exp(coef)"]),
		se    = unname(summ$coef[,"se(coef)"]),
		z     = unname(summ$coef[,"z"]),
		p     = unname(summ$coef[,"Pr(>|z|)"]),
		stringsAsFactors = FALSE
	)
	df
}

#' Tabulate multivariable and single variable Cox fit results next to each other
#' @export
coxtab <- function(
	# Model object; either formula or coxph
	obj,
	# If object is formula, data has to be supplied
	data,
	# Custom formatting; integer values depict preformatted choices
	format,
	# How many digits precision to use; named vector for different formattings
	digits = c("coef" = 3, "p" = 3),
	# Should effective N counts be added
	addN = TRUE,
	...
){
	if(inherits(obj, "coxph")){
		cx <- coxambi.coxph(obj, ...)
	}else if(inherits(obj, "formula")){
		cx <- coxambi.formula(obj, data=data, ...)	
	}else{
		stop(paste0("Invalid class for obj: ", class(obj)))
	}

	# Multivariable fit coef extraction	
	multidf <- coxdf(cx[[1]])
	colnames(multidf) <- paste0("multi_", colnames(multidf))
	# Univariable fits coef extractions
	unidfs <- do.call("rbind", lapply(cx[[2]], FUN=\(x){
		tmp <- coxdf(x)
		colnames(tmp) <- paste0("uni_", colnames(tmp))
		tmp
	}))[,-1]

	# Iterate across model terms and take into account their type (first element is LHS)
	types <- attr(terms(cx[[1]]), "dataClasses")[-1]
	# Iterate and aggregate Ns
	Ns <- lapply(1:length(types), FUN=\(i){
		if(types[i] %in% c("logical")){
			#c(FALSE, TRUE)
			# Should be proper logicals in order FALSE, TRUE
			tab <- table(model.frame(cx[[1]])[,names(types)][,i])
			names(tab) <- paste0(names(types)[i], names(tab))
		}else if(types[i] %in% c("factor", "character")){
			# Should be factor levels in the same order as levels(var)
			tab <- table(model.frame(cx[[1]])[,names(types)][,i])
			names(tab) <- paste0(names(types)[i], names(tab))
		}else if(types[i] %in% c("numeric", "integer")){
			tab <- sum(!is.na(model.frame(cx[[1]])[,names(types)][,i]))
			names(tab) <- names(types)[i]
		}else{
			stop(paste0("Invalid dataClasses instance: ", types[i], " for variable index ", i))
		}
		tab
	})
	# Unlist all the levels of variables, with corresponding N counts
	Ns <- unlist(Ns)

	# Take multivariable fits and univariable side-by-side
	df <- cbind(multidf, unidfs)
	rownames(df) <- df[,1]
	df <- df[,-1]
	
	if(missing(format)){
		df
	# Give as 'rowname ~ "Multi HR [0.95% CI], p", "Uni HR [0.95% CI], p"'
	}else if(format == 1){
		# Multivariable
		multi_hrest <- round(df[,"multi_hr"], digits["coef"])
		multi_low95 <- round(exp(df[,"multi_coef"] - 1.96 * df[,"multi_se"]), digits["coef"])
		multi_top95 <- round(exp(df[,"multi_coef"] + 1.96 * df[,"multi_se"]), digits["coef"])
		multi_p <- ifelse(df[,"multi_p"] < 0.001, "p<0.001", paste0("p=", round(df[,"multi_p"], digits["p"])))
		
		# Univariable(s)
		uni_hrest <- round(df[,"uni_hr"], digits["coef"])
		uni_low95 <- round(exp(df[,"uni_coef"] - 1.96 * df[,"uni_se"]), digits["coef"])
		uni_top95 <- round(exp(df[,"uni_coef"] + 1.96 * df[,"uni_se"]), digits["coef"])
		uni_p <- ifelse(df[,"uni_p"] < 0.001, "p<0.001", paste0("p=", round(df[,"uni_p"], digits["p"])))
		
		df_formatted <- data.frame(
			Multivariable = paste0(multi_hrest, " [", multi_low95, ", ", multi_top95, "], ", multi_p),
			Univariable = paste0(uni_hrest, " [", uni_low95, ", ", uni_top95, "], ", uni_p)
		)
		rownames(df_formatted) <- rownames(df)

		# If user wants to show N counts per level, add those as first column; also add additional reference rows
		if(addN){
			df_formatted <- data.frame(
				# N counts per variable / level
				N = Ns, 
				# Multivariable fit summary
				Multivariable = df_formatted[match(names(Ns), rownames(df_formatted)), "Multivariable"],
				# Aggregated univariable fit summaries
				Univariable = df_formatted[match(names(Ns), rownames(df_formatted)), "Univariable"]
			)
			# Empty cells (NAs) are on the reference rows for logical/factor variables
			df_formatted[is.na(df_formatted)] <- "Reference"
			# Use the extended names
			rownames(df_formatted) <- names(Ns)
		}

		df <- df_formatted
	}else{
		stop(paste0("Invalid 'format', should be an integer: ", format))
	}
		
	df 
}
