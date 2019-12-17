library(intsurv)

## varify results of toy examples
time1 <- c(0.5, 1, 2, 3.2, 4.8, 5)
event1 <- c(0, 1, 1, 1, 0, 0)
risk_score1 <- c(0.1, 0.5, 0.5, - 1, 0, 0.1)

## all one weight
res1 <- cIndex(time1, event1, risk_score1)
expect_equivalent(res1["comparable"], 9.0)
expect_equivalent(res1["tied_risk"], 1.0)
expect_equivalent(res1["concordant"], 6.0)
expect_equivalent(res1["index"], 6.5 / 9)

## equal weights
w0 <- 0.5
w1 <- rep(w0, length(time1))
res2 <- cIndex(time1, event1, risk_score1, w1)
expect_equivalent(res2["comparable"], res1["comparable"] * w0)
expect_equivalent(res2["tied_risk"], res1["tied_risk"] * w0)
expect_equivalent(res2["concordant"], res1["concordant"], 4.0 * w0)
expect_equivalent(res2["index"], 6.5 / 9)

## not equal weights
w2 <- c(0.5, 1, 1, 1, 0.1, 0.5)
res3 <- cIndex(time1, event1, risk_score1, w2)
expect_equivalent(res3["comparable"], 4.8)
expect_equivalent(res3["tied_risk"], 1.0)
expect_equivalent(res3["concordant"], 3.2)
expect_equivalent(res3["index"], 3.7 / 4.8)
