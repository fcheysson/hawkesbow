test_that("Exponential integral of imaginary argument", {

    # Can be checked in
    # Nat. Bur. Standards, Appl. Math. Ser. No. 51, Tables of the Exponential
    # Integral for Complex Arguments, U. S. Government Printing Office,
    # Washington 25, D. C, May 15, 1958

    z = seq(0, .9, by = .1)
    reals = c(
        Inf,
        1.727868,
        1.042206,
        .649173,
        .378809,
        .177784,
        .022271,
        -.100515,
        -.198279,
        -.276068
    )
    imags = c(
        - .5 * pi,
        -1.470852,
        -1.371240,
        -1.272292,
        -1.174335,
        -1.077689,
        -.982668,
        -.889574,
        -.798701,
        -.710326
    )
    E1 = sapply(z, E1_imaginary)
    expect_equal(Re(E1), reals, tolerance = 1e-5)
    expect_equal(Im(E1), imags, tolerance = 1e-5)

    z = seq(1, 10, by = .5)
    reals = c(
        -.337404,
        -.470356,
        -.422981,
        -.285871,
        -.119630,
        .032129,
        .140982,
        .193491,
        .190030,
        .142053,
        .068057,
        -.011102,
        -.076695,
        -.115633,
        -.122434,
        -.099431,
        -.055348,
        -.002678,
        .045456
    )
    imags = c(
        -.624713,
        -.246113,
        .034617,
        .207724,
        .277856,
        .262329,
        .187407,
        .083344,
        -.020865,
        -.102072,
        -.146109,
        -.149002,
        -.116200,
        -.060115,
        .003391,
        .058801,
        .094244,
        .103667,
        .087551
    )
    E1 = sapply(z, E1_imaginary)
    expect_equal(Re(E1), reals, tolerance = 1e-5)
    expect_equal(Im(E1), imags, tolerance = 1e-5)

})


test_that("Incomplete gamma function of imaginary argument", {

    # Can be checked using function gamma from R::base, Gradshteyn & Ryzhik, 2007,
    # and Barakat, 1960:
    # inc_gamma_imag(r, b) = i^{-b} \Gamma(b) - r^b \gamma_1(b, ir)
    # where first term from Gradshteyn & Ryzhik, 2007, by contour integration (3.381.7),
    # second term from Barakat, 1960, by Chebyshev polynomials

    b = .6
    r = 1.0
    expect_equal(inc_gamma_imag(r, b),
                 (-1i)^b * gamma(b) - (1.483209 - .580166i) * (r^b),
                 tolerance = 1e-5)

    b = .8
    r = .5
    expect_equal(inc_gamma_imag(r, b),
                 (-1i)^b * gamma(b) - (1.205896 - .272340i) * (r^b),
                 tolerance = 1e-5)

    b = .2
    r = 3.0
    expect_equal(inc_gamma_imag(r, b),
                 (-1i)^b * gamma(b) - (3.613147 - 1.428423i) * (r^b),
                 tolerance = 1e-5)

    b = .4
    r = 5.0
    expect_equal(inc_gamma_imag(r, b),
                 (-1i)^b * gamma(b) - (.750625 - .650291i) * (r^b),
                 tolerance = 1e-5)

})

test_that("PowerLaw::H gives correct results", {

    # theta = 0.4
    xi = seq(0, 3.5, by = .5)
    model = new(PowerLaw)
    model$param[3] = .4
    mu = model$param[2]
    theta = model$param[3]
    a = model$param[4]

    b = 1.0 - theta
    r = xi * a
    gamma1 = c(1.666667 + .000000i,
               1.619153 - .306759i,
               1.483209 - .580166i,
               1.277506 - .792115i,
               1.029608 - .924050i,
               0.771356 - .969553i,
               0.533655 - .934797i,
               0.341618 - .836824i)

    expect_equal(as.vector(model$H(xi)),
                 mu * (1.0 - 1.0i * r * exp(1.0i * r) * ifelse(r == 0, 0, exp(-b * log(r))) * ((-1.0i)^b * gamma(b) - (r^b) * gamma1)),
                 tolerance = 1e-5)

    # theta = 1
    xi = seq(0, 3.5, by = .5)
    model = new(PowerLaw)
    model$param[3] = 1
    mu = model$param[2]
    theta = model$param[3]
    a = model$param[4]

    r = xi * a
    E1 = c(     Inf - .5i * pi ,
            .177784 - 1.077689i,
           -.337404 -  .624713i,
           -.470356 -  .246113i,
           -.422981 +  .034617i,
           -.285871 +  .207724i,
           -.119630 +  .277856i,
            .032129 +  .262329i)

    expect_equal(as.vector(model$H(xi)),
                 mu * (1.0 - ifelse(r == 0, 0, 1.0i * r * exp(1.0i * r) * E1)),
                 tolerance = 1e-5)

    # theta = 1.4
    xi = seq(0, 3.5, by = .5)
    model = new(PowerLaw)
    model$param[3] = 1.4
    mu = model$param[2]
    theta = model$param[3]
    a = model$param[4]

    b = 2.0 - theta
    r = xi * a
    gamma1 = c(1.666667 + .000000i,
               1.619153 - .306759i,
               1.483209 - .580166i,
               1.277506 - .792115i,
               1.029608 - .924050i,
               0.771356 - .969553i,
               0.533655 - .934797i,
               0.341618 - .836824i)

    expect_equal(as.vector(model$H(xi)),
                 mu * (1.0 - 1.0i * r / (theta - 1.0) - (r^2 / (theta - 1.0)) * exp(1.0i * r) * ifelse(r == 0, 0, exp(-b * log(r))) * ((-1.0i)^b * gamma(b) - (r^b) * gamma1)),
                 tolerance = 1e-5)

    # theta = 2
    xi = seq(0, 3.5, by = .5)
    model = new(PowerLaw)
    model$param[3] = 2
    mu = model$param[2]
    theta = model$param[3]
    a = model$param[4]

    r = xi * a
    E1 = c(     Inf - .5i * pi ,
                .177784 - 1.077689i,
                -.337404 -  .624713i,
                -.470356 -  .246113i,
                -.422981 +  .034617i,
                -.285871 +  .207724i,
                -.119630 +  .277856i,
                .032129 +  .262329i)

    expect_equal(as.vector(model$H(xi)),
                 mu * (1.0 + ifelse(r == 0, 0, - 1.0i * r / (theta - 1.0) - (r^2 / (theta - 1.0)) * exp(1.0i * r) * E1)),
                 tolerance = 1e-5)

})
