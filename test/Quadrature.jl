using Quadrature

export _test_all

function _test_all()
    # test sqrt
    (s_ans, l_ans) = test_sqrt()
    @assert((s_ans - 2/3) < 1.0e-6)
    @assert((l_ans - 2/3) < 1.0e-6)
    print("sqrt within tolerance.\n")

    # test a pathological integrand
    (s_ans, l_ans) = test_pathological()
    print("pathological function recursion terminated successfully.\n")
    
    # test 1/(1+x)
    (s_ans, l_ans) = test_1px()
    @assert((s_ans - log(4)) < 1.0e-6)
    @assert((l_ans - log(4)) < 1.0e-6)

    print("1/(1+x) within tolerance.\n")
end

function test_sqrt()
    f = x -> sqrt(x)
    truth = 2/3
    simps_ans = Quadrature.adapt_simpson(f, 0.0, 1.0, 1.0e-7)
    lobat_ans = Quadrature.adapt_lobatto(f, 0.0, 1.0, 1.0e-7)
    simps_ans, lobat_ans
end

function test_pathological()
    f = x -> 1/sqrt(1-x^2)
    truth = pi/2
    simps_ans = Quadrature.adapt_simpson(f, 0.0, 1.0, 1.0e-7)
    lobat_ans = Quadrature.adapt_lobatto(f, 0.0, 1.0, 1.0e-7)
    simps_ans, lobat_ans
end

function test_1px()
    f = x -> 1/(1+x)
    truth = log(4)
    simps_ans = Quadrature.adapt_simpson(f, 1.0, 7.0, 1.0e-7)
    lobat_ans = Quadrature.adapt_lobatto(f, 1.0, 7.0, 1.0e-7)
    simps_ans, lobat_ans
end