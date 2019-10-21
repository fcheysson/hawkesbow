test_that("XXXX", {

    # Create discrete sample
    T = runif(1, 50, 200)
    x = sort(runif(100, 0, T))
    length = sample(x = 1:20, size = 1)
    cdata = new(ContinuousData, x, T)
    ddata = cdata$toDiscrete_byLength(length)

    # Correct bins in R
    ti <- seq(0, T, length.out=length+1)
    bin <- sapply(1:length, function(i) {
        sum(x > ti[i] & x < ti[i+1])
    })

    # "expect" methods
    expect_equal( length(ddata$counts), length )
    expect_equal( as.numeric(ddata$counts), bin )
    expect_equal( ddata$binsize, T/length )
})
