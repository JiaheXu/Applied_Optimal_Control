syms a b c u x_t0 x_t x_tf c1 tf t0 t p_t real;
t0 = 0;
c1 = c * exp(a*(tf-t0)) * x_t0 / ...
    ( exp(-a*tf) + c*b*b*exp(a*tf)*(-exp(-2*a*tf)+exp(-2*a*t0) )/(2*a) );
x_t = exp(a*(t-t0))*x_t0 - c1*b*b*exp(a*t)*(-exp(-2*a*t)+exp(-2*a*t0) )/(2*a);
u1 = b*c1*exp(-a*t);

p_t = 2*a*c / (b*b*c + ( 2*a - b*b*c )*exp( 2*a*(t-tf) ) );
u2 = b * p_t * x_t;
simplify( u1 - u2 )
