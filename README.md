Quadrature implements the adaptive quadrature integration routines which are
discussed in "Adaptive Quadrature: Revisited" by Gander and Gautschi - 
the adaptive simpson and adaptive lobatto methods.  

To use the routines, simply pass a 1-d integrand as a function handle along with
limits and desired tolerance as in:

integral = adapt_simpson(x->sqrt(x), 0, 1, 1.0e-6)

or

integral = adapt_lobatto(x->sqrt(x), 0, 1, 1.0e-6)

The numerical accuracy and timing performance is about the same as the quad and 
quadgk functions in MATLAB.
