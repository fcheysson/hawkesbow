test_that("Methods for spectral densities work", {

    # Create exponential Hawkes model
    model = new(Exponential)
    model$param=c(1,.5,2)

    # Excitation function
    h0 =


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
