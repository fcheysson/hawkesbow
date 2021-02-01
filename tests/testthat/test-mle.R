test_that("Loglik and its derivatives are correct", {

    # Create exponential Hawkes model
    model = new(Exponential)
    model$param=c(1.0,.5,log(2.0))

    # Create dataset
    events = 1.0:4.0
    end = 5.0

    # Useful outputs
    A = c(0.0, 1.0/2.0, 3.0/4.0,  7.0/8.0,  15.0/16.0)
    C = c(0.0, 1.0/2.0, 5.0/4.0, 17.0/8.0,  49.0/16.0)
    B = c(0.0, 1.0/2.0,     1.0, 11.0/8.0,  26.0/16.0)
    E = c(0.0, 1.0/2.0, 9.0/4.0, 45.0/8.0, 173.0/16.0)
    D = c(0.0, 1.0/2.0, 6.0/4.0, 21.0/8.0,  58.0/16.0)

    # Expected outputs
    denom = 1.0 + .5 * log(2.0) * A[1:4]
    denom2 = denom * denom
    loglik = sum( log(denom) ) - 209.0/32.0
    dloglik = c(
        sum( 1.0 / denom ) - 5.0,
        sum( log(2.0) * A[1:4] / denom ) - 49.0/16.0,
        sum( (.5 * A[1:4] - .5 * log(2.0) * B[1:4]) / denom ) - 26.0/32.0
    )
    ddloglik = matrix(c(
        `11`=- sum( 1 / denom2 ),
        `12`=- sum( log(2.0) * A[1:4] / denom2 ),
        `13`=- sum( (.5 * A[1:4] - .5 * log(2.0) * B[1:4]) / denom2 ),
        `21`=- sum( log(2.0) * A[1:4] / denom2 ),
        `22`=- sum( log(2.0) * log(2.0) * A[1:4] * A[1:4] / denom2 ),
        `23`=+ sum( (A[1:4] - log(2.0) * B[1:4]) / denom2 ) - B[5],
        `31`=- sum( (.5 * A[1:4] - .5 * log(2.0) * B[1:4]) / denom2 ),
        `32`=+ sum( (A[1:4] - log(2.0) * B[1:4]) / denom2 ) - B[5],
        `33`=+ sum( (.5 * (log(2.0) * D[1:4] - 2.0 * B[1:4]) + .25 * (log(2.0) * log(2.0) * (A[1:4] * E[1:4] - C[1:4] * C[1:4]) - A[1:4] * A[1:4])) / denom2 ) + .5 * D[5]
    ), byrow = TRUE, nrow = 3, ncol = 3)

    # expect_equals
    expect_equal(  model$loglik(events, end)    ,   loglik )
    expect_equal( model$dloglik(events, end)[,1],  dloglik )
    expect_equal(model$ddloglik(events, end)    , ddloglik )
    expect_equal(model$loglikngrad(events, end) ,
                 list(objective=model$loglik(events, end),
                      gradient =model$dloglik(events, end)))

})
