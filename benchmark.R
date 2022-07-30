library(MatchIt)
data(lalonde)

benchmarkClosestIndex <- function(n) {
  trt_score = data.frame(index=c(1:n), score=runif(n))
  ctr_score = data.frame(index=c(1:(5*n)), score=runif(5*n))

  print(system.time(
    blazingMatch:::closestIndexCpp(ctr_score, trt_score)
  ))
  print(system.time(
    closestIndex(ctr_score, trt_score)
  ))
}
benchmarkClosestIndex(100000)

benchmark <- function(n) {
  fake_data <- sapply(lalonde, rep.int, time=n)
  fake_data <- data.frame(fake_data)
  names(fake_data) <- names(lalonde)

  confounding_vars <- c('age', 'educ', 'race', 'nodegree', 'married', 're74', 're75')
  print(system.time(
    result <- blazingMatch(data=fake_data, treatment_var='treat', confounding_vars = confounding_vars, treatment_group = 1)
  ))

  confounding_vars <- c('age', 'educ', 'race', 'nodegree', 'married', 're74', 're75')
  print(system.time(
    result <- blazingMatch(data=fake_data, treatment_var='treat', confounding_vars = confounding_vars, treatment_group = 1, fast=TRUE)
  ))

  print(system.time(
    m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                        married + re74 + re75, data = fake_data)
  ))
}

benchmark(1000)
