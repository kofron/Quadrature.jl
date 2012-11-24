module Quadrature

using Base

export adapt_simpson, adapt_lobatto

# compute quadrature integral estimate via adaptive simpson rule.
function adapt_simpson(f::Function, a::Float64, b::Float64, tol)

    x = [a (a+b)/2 b]
    y = [f(p) for p in x]
    yy = [f(p) for p in (a + [0.9501 0.2311 0.6068 0.4860 0.8913]*(b-a))]

    machine_eps = eps(x[2])
    
    is = (b-a)/8*(sum(y)+sum(yy))
    is = (is == 0 ? b-a : is)*tol/machine_eps
    adapt_simpson1(f,a,b,y[1],y[2],y[3],is)
end

function adapt_simpson1(f,a,b,fa,fm,fb,is)
    m, h = (a+b)/2, (b-a)/4
    x = [a+h, b-h]
    y = [f(p) for p in x]

    i1 = h/1.5 * (fa + 4*fm + fb)
    i2 = h/3 * (fa + 4*(y[1] + y[2]) + 2*fm + fb)
    i1 = (16*i2 - i1)/15

    # check termination criteria and recurse if not met.
    if (is + (i1-i2) == is) || (m <= a) || (b <= m)
        i1
    else
        adapt_simpson1(f,a,m,fa,y[1],fm,is) + adapt_simpson1(f,m,b,fm,y[2],fb,is)
    end 
end

# compute quadrature integral estimate via adaptive lobatto rule.
function adapt_lobatto(f::Function, a::Float64, b::Float64, tol)

    m, h = (a+b)/2, (b-a)/2
    al, bt = sqrt(2/3), 1/sqrt(5)

    x1, x2, x3 = 0.942882415695480, 0.641853342345781, 0.236383199662150

    x = [a,m-x1*h,m-al*h,m-x2*h,m-bt*h,m-x3*h,m,m+x3*h,
         m+bt*h,m+x2*h,m+al*h,m+x1*h,b]
    y = [f(p) for p in x]

    machine_eps = eps(m)

    i2 = (h/6)*(y[1] + y[end] + 5*(y[5]*y[9]))
    i1 = (h/1470)*(77*(y[1]+y[end])+432*(y[3] + y[11]) + 625*(y[5]+y[9])+672*y[7])
    is = 0.0158271919734802*(y[1]+y[end]) 
    is = is + 0.0942738402188500*(y[2]+y[12]) 
    is = is + 0.155071987336585*(y[3]+y[11]) 
    is = is + 0.188821573960182*(y[4]+y[10])
    is = is + 0.199773405226859*(y[5]+y[9])
    is = is + 0.224926465333340*(y[6]+y[8]) 
    is = is + 0.242611071901408*y[7]
    is = h*is

    s = (sign(is) == 0) ? 1 : sign(is)
    ei1 = abs(i1 - is)
    ei2 = abs(i2 - is)
    eR = ei1/ei2
    
    tol = (eR > 0 && eR < 1) ? tol/eR : tol
    is = s*abs(is)*tol/machine_eps

    is = (is == 0) ? b-a : is

    adapt_lobatto1(f, a, b, y[1], y[end], is, al, bt)

end # end adapt_lobatto function

function adapt_lobatto1(f, a, b, fa, fb, is, al, bt)
    h, m = (b-a)/2, (b+a)/2
    
    mll, ml, mr, mrr  = m - al*h, m - bt*h, m + bt*h, m + al*h
    y = [f(p) for p in [mll,ml,m,mr,mrr]]
    
    i2 = (h/6)*(fa + fb + 5*(y[2] + y[4]))
    i1 = (h/1470)*(77*(fa+fb)+432*(y[end] + y[1]) + 625*(y[2]+y[4])+672*y[3])

    if (is + (i1-i2) == is) || (mll <= a) || (b <= mrr)
        i1
    else
        adapt_lobatto1(f,a,mll,fa,y[1],is,al,bt) + adapt_lobatto1(f,mll,ml,y[1],y[2],is,al,bt) + adapt_lobatto1(f,ml,m,y[2],y[3],is,al,bt) + adapt_lobatto1(f,m,mr,y[3],y[4],is,al,bt) + adapt_lobatto1(f,mr,mrr,y[4],y[5],is,al,bt) + adapt_lobatto1(f,mrr,b,y[5],fb,is,al,bt)
    end 

    
end


end # module quadrature