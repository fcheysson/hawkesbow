test_that("Methods for spectral densities work", {

    # Create exponential Hawkes model
    model = new(Exponential)
    model$param=c(1,.5,2)

    # For w = 0
    expect_equal(   model$f(0)[1,], 8 )
    expect_equal(  model$df(0)[1,], c(8, 48, 0) )
    expect_equal(model$ddf(0)[,,1], matrix(c( 0,  48, 0,
                                             48, 384, 0,
                                              0,   0, 0),
                                           nrow=3, ncol=3) )

    # For w = 1
    s2 = (2.0 * sin(.5))^2  # sinc²(1.0 / 2.0)
    expect_equal(   model$f(1)[1,], 5*s2)
    expect_equal(  model$df(1)[1,], c(5*s2, 20*s2, 1.5*s2))
    expect_equal(model$ddf(1)[,,1], matrix(c(     0,  20*s2,  1.5*s2,
                                              20*s2, 100*s2,   11*s2,
                                             1.5*s2,  11*s2, -.75*s2),
                                           nrow=3, ncol=3) )

    # # For w = 2
    # s2 = sin(1)^2 # sinc²(2.0 / 2.0)
    # expect_equal(   model$f(2)[1,], 3.2*s2)
    # expect_equal(  model$df(2)[1,], c(3.2*s2, 8.96*s2, .96*s2) )
})
