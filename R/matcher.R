#' @title quick nearest neighbor search for propensity scores
#'
#' @param ctr_score data frame containing index and propensity scores for control group
#' @param trt_score data frame containing index and propensity scores for treatment group
#'
#' @return
#'
#'
#' @import dplyr, magrittr
closest_index <- function(ctr_score, trt_score){
  n0 <- nrow(ctr_score)
  n1 <- nrow(trt_score)

  ctr_score <- ctr_score[order(ctr_score$score), ]
  trt_score <- trt_score[order(trt_score$score), ]

  results <- data.frame(
    control = rep(0, n1),
    treatment = rep(0, n1)
  )

  index0 <- 1
  index1 <- 1

  while (index0 <= n0 & index1 <= n1) {
    while (index0 + 1 < n0 & ctr_score$score[index0+1] < trt_score$score[index1]) {
      index0 <- index0 + 1
    }

    if (index0 == n0) {
      for (i in index1:n1){
        results[i, ] <- c(ctr_score$index[index0], trt_score$index[i])
      }
      return(results)
    }

    if (ctr_score$score[index0+1] - trt_score$score[index1] >= trt_score$score[index1] - ctr_score$score[index0]) {
      results[index1, ] <- c(ctr_score$index[index0], trt_score$index[index1])
    } else {
      results[index1, ] <- c(ctr_score$index[index0+1], trt_score$index[index1])
      index0 <- index0 + 1
    }
    index1 <- index1 + 1
  }

  return(results)

}

#' @title blazing fast propensity score matching
#'
#' @param data data frame containing data
#' @param treatment_var treatment variable name
#' @param confounding_vars confounding variable names
#' @param treatment_group value corresponding to the treatment group in treatment variable
#'
#' @return list containing two data frame `treatment` and `control`.
#'
#'
#' @import
#' @export
#' @examples
#' data(lalonde)
#' confounding_vars <- c('age', 'educ', 'race', 'nodegree', 'married', 're74', 're75')
#' result <- blazingMatch(data=lalonde, treatment_var='treat', confounding_vars = confounding_vars, treatment_group = 1)
#' summary(result$treatment)
#' summary(result$control)
blazingMatch <- function(data, treatment_var, confounding_vars, treatment_group){
  group_labels <- unique(data[[treatment_var]])
  if ('score' %in% names(data)) {
    stop("data can not have column named as score")
  }
  if (length(group_labels) != 2) {
    stop("currently only support binary treatment variable")
  }
  if (! treatment_group %in% group_labels) {
    stop("treatment_group has to be one of the two groups in treatment variable")
  }

  data[treatment_var] <- as.integer(data[[treatment_var]]==treatment_group)

  formula_string <- paste(treatment_var, paste(confounding_vars, collapse=" + "), sep=" ~ ")
  model <- glm(as.formula(formula_string), data = data, family = "binomial")

  data['score'] <- predict(model)

  trt <- data[data[treatment_var]==1,]
  ctr <- data[data[treatment_var]==0,]

  trt_score <- data.frame(index=c(1:nrow(trt)), score = trt$score)
  ctr_score <- data.frame(index=c(1:nrow(ctr)), score = ctr$score)

  matched_index <- closest_index(ctr_score, trt_score)
  matched_treatment <- trt[matched_index$treatment, ]
  matched_control <- ctr[matched_index$control, ]

  return(list(treatment=matched_treatment, control=matched_control))
}

data("lalonde")

result <- blazingMatch(data=lalonde, treatment_var='treat',
             confounding_vars = c('age', 'educ', 'race', 'nodegree', 'married', 're74', 're75'),
             treatment_group = 1)

summary(result$treatment)
summary(result$control)
