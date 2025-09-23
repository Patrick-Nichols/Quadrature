
#include <cstdio>
#include <cmath>
#include <limits>

struct Fun1
{
    double operator()(const double& x) const noexcept
    {
        return exp(-x)*(x+1);
    }
};

struct Fun2
{
    double operator()(const double& x) const noexcept
    {
        return 6./(x*x+1.);
    }
};

struct Fun3
{
    double operator()(const double& x) const noexcept
    {
        if (x == 0) {
            return 0.0;
        }
        return log(x)/sqrt(x);
    }
};

struct Fun4
{
    double operator()(const double& x) const noexcept
    {
        if (x == 0) {
            return 1.e100;
        }
        double vs = (1.0 - exp(-x*2))/x;
        double vp = 4.5/pow(x,4);
        vp *= ( 1 - exp(-0.25*x));
        return vs - vp;
    }
};

template < class Func_t >
double adaptive_integration2(const double& a, const double& b, size_t& ns,Func_t& fun,
    double& rel_err,
    const double& tolerance, int level = -1)
{
    const double eps = 10.0 * std::numeric_limits<double>::epsilon();
    const double max_np = (b-a)/eps;
    size_t np = ns;
    double dx = (b-a)/np;
    if ( np >= max_np) {
        np = max_np;
        level = 0;
        dx = 10. * eps;
    }
    if (level == -1) {
        level = static_cast<size_t>(log2(max_np)); 
    }
    fprintf(stderr,"a = %20.10le b = %20.10le npts = %12lu dx = %20.10le level = %12d\n",a,b,np,dx,level);
    double ends = (fun(b) + fun(a))*0.5;
    // do initial integration
    double sum1 = ends;
    for (size_t j=1;j<np;++j) {
        double x = a + dx * j;
        sum1 += fun(x);
    }
    double int1 = sum1 * dx;
    double sum2 = sum1;
    np = np + np;
    dx *= 0.5;
    for (size_t j=1;j<np;j+=2)
    {
        sum1 += fun((a+dx*j));
    }
    double int2 = sum1 * dx;
    double diff = fabs(int2-int1);
    double den = fabs(int2);
    den = (den > eps) ? den : 1.0;
    rel_err = diff/den;
    if (rel_err < tolerance)
    {
            fprintf(stdout,"converged a = %10.5le b = %10.5le ",a,b); 
            fprintf(stdout,"dx = %15.6le  ns = %12lu err = %20.10le\n",dx,np,rel_err);
            ns = np;
            return int2;
    }
    if (level == 0) {
            fprintf(stdout,"no convergence a = %10.5le b = %10.5le ",a,b); 
            fprintf(stderr,"dx = %15.6le ns = %12lu err= %20.10le\n",dx,np,rel_err);
            return int2;
    }
    double dh = (b-a)*0.25;
    double rel_err_;
    size_t nsum = 0;
    int2 = 0.0;
    for (int k=0;k<4;++k)
    {
        double end1 = a + dh * k;
        double end2 = end1 + dh;
        int2 += adaptive_integration2(end1,end2,ns,fun,rel_err_,tolerance,(level-1));
        rel_err = fmax(rel_err_,rel_err);
        nsum += ns;                 
    }  
    ns = nsum;
    return int2;
}

template < class Func_t >
double adaptive_integration(const double& a, const double& b, size_t& ns, Func_t& fun,
    double& rel_err,
    const double& tolerance)
{
    bool low = false;
    const double eps = std::numeric_limits<double>::epsilon() * 10.0;
    size_t maxlevel = static_cast<size_t>(log2((b-a)/eps/ns));
    double dx = (b-a)/ns;
    if ( dx < eps) {
        dx = eps;
        ns = (b-a)/dx;
        low = true;
    }    
    size_t np = ns;
    double ends = (fun(b) + fun(a))*0.5;
    // do initial integration
    double sum1 = ends;
    for (size_t j=1;j<np;++j) {
        double x = a + dx * j;
        sum1 += fun(x);
    }
    double int1 = sum1 * dx;
    double sum2 = sum1;
    double int2 = 0.0;
    for (size_t m=0;m<maxlevel;++m)
    {
        dx *= 0.5;
        np = np + np;
        for (size_t j=1;j<np;j+=2) 
        {
            double x = a + dx * j;
            sum2 += fun(x);
        }	
        int2 = sum2 * dx;
        double diff = fabs(int2-int1);
        double den = fabs(int2);
        den = (den > eps) ? den : 1.0;
        rel_err = diff/den;
        if (rel_err < tolerance)
        {
            fprintf(stdout,"converged for dx = %15.6le  ns = %12lu\n",dx,np);
            ns = np;
            return int2;
        }
        if ( low ) {
            fprintf(stdout,"too small dx = %15.6le ns = %12lu\n",dx,np);
            ns = np;
            return int2;
        }
        fprintf(stderr,"int2 = %20.10le int1 = %20.10le err = %20.10le npts = %12lu\n",
            int2,int1,rel_err,np);
        sum1 = sum2;
        int1 = int2;            
    }
    fprintf(stderr,"adaptive int did not converge!!!\n");
    ns = np;
    return int2;
}

template < class Func_t >
double trap_integration(const double& a, const double& b, const size_t& ns, Func_t& fun)
{
    double sum{0.};
    double dx = (b-a)/ns;
    sum = (fun(a) + fun(b))*0.5;
    for (size_t j=1;j<ns;++j)
    {
        double x = a + j * dx;
        sum += fun(x);
    }
    return (sum*dx);
}

template < class Func_t >
double simp_integration(const double& a, const double& b, const size_t& ns, Func_t& fun)
{
    size_t np = (ns % 2) ? ((ns/2)+2):ns;
    double sum{0.};
    const double dx = (b-a)/ns;
    double sum1 = 0.0;
    for (size_t j=1;j<np;j+=2)
    {
        double x = a + j * dx;
        sum1 += fun(x);
    }
    double sum2 = 0.0;
    for (size_t j=2;j<np;j+=2)
    {
        double x = a + j * dx;
        sum2 += fun(x);
    }
    return dx * (fun(a) + fun(b) + sum1 * 4. + sum2 * 2.)/3.; 
}

template < class Func_t >
double boole_integration(const double& a, const double& b, const size_t& ns, Func_t& fun)
{
    size_t np = (ns % 4) ? ((ns/4)+4) : ns;
    const double dx = (b-a)/ns;
    double sum1 = 0.0;
    for (size_t j=1;j<np;j+=2)
    {
        double x = a + j * dx;
        sum1 += fun(x);
    }
    sum1 *= 32.;
    double sum2 = 0.0;
    for (size_t j=2;j<np;j+=4)
    {
        double x = a + j * dx;
        sum2 += fun(x);
    }
    sum2 *= 12.0;
    double sum3 = 0.0;
    for (size_t j=4;j<np;j+=4)
    {
        double x = a + j * dx;
        sum3 += fun(x);
    }
    sum3 *= 14.0;
    return dx * 2. * ( 7.*(fun(a) + fun(b)) + sum1 + sum2 + sum3)/45.0; 
}

/** Adaptive Simpson's Rule, Recursive Core */
template < class Func_t >
double adaptiveSimpsonsAux(Func_t& fun, double a, double b, double eps,
                          double whole, double fa, double fb, double fm, int rec) {
    double m   = 0.5*(a + b);
    double h   = 0.5*(b - a);
    double lm = 0.5 *(a + m);
    double rm = 0.5 *(m + b);
    // serious numerical trouble: it won't converge
//    if ((eps/2 == eps) || (a == lm)) { errno = EDOM; return whole; }
    double flm = fun(lm);
    double frm = fun(rm);
    double wt = h/6.0;
    double left = wt * ( fa + 4 * flm + fm);
    double right = wt * ( fm + 4. * frm + fb);
    double delta = left + right - whole;
    if ( rec <= 0 || fabs(delta)<=eps) return (left+right+delta/15.);
    return adaptiveSimpsonsAux(fun,a,m,eps*0.5,left,fa,fm,flm,rec-1) +
           adaptiveSimpsonsAux(fun,m,b,eps*0.5,right,fm,fb,frm,rec-1);
}

/** Adaptive Simpson's Rule Wrapper
 *  (fills in cached function evaluations) */
template < class Func_t >
double adapt_simpsons_integrate(Func_t& fun,
                       double a, double b,      // interval [a,b]
                       double epsilon,         // error tolerance
                       int maxRecDepth) {     // recursion cap
//    errno = 0;
    double h = b - a;
    if (h == 0) return 0;
    double fa = fun(a);
    double fb = fun(b);
    double fm = fun((a + b)/2);
    double S = (h/6.)*(fa + 4.*fm + fb);
    return adaptiveSimpsonsAux(fun, a, b, epsilon, S, fa, fb, fm, maxRecDepth);
}

int main()
{
    double a = 0.0;
    double b = 1.0;
    double err;
    double tol = 1.e-12;
    size_t ns = 16;
    Fun1 f1;
    Fun2 f2;
    Fun3 f3;
    Fun4 f4;
    const size_t ns0 = ns;
        
    double t1 = trap_integration<Fun1>(a,b,ns0,f1);
    double s1 = simp_integration<Fun1>(a,b,ns0,f1);
    double b1 = boole_integration<Fun1>(a,b,ns0,f1);
    ns = ns0;
    double a1 = adaptive_integration<Fun1>(a,b,ns,f1,err,tol);
    double t1b = trap_integration<Fun1>(a,b,ns,f1);
    size_t np = ns;
    ns = ns0;
    double a1_2 = adaptive_integration2<Fun1>(a,b,ns,f1,err,tol);
    double s1a = adapt_simpsons_integrate<Fun1>(f1,a,b,tol,14);

    fprintf(stdout," trap  %12lu %20.10le\n",ns0,t1);
    fprintf(stdout," simp  %12lu %20.10le\n",ns0,s1);
    fprintf(stdout," bool  %12lu %20.10le\n",ns0,b1);
    fprintf(stdout," adapt %12lu %20.10le\n",np,a1);
    fprintf(stdout," adapt %12lu %20.10le\n",ns,a1_2);
    fprintf(stdout," trap  %12lu %20.10le\n",np,t1b);
    fprintf(stdout," asimp              %20.10le\n",s1a);
    fprintf(stdout,"\n");
    
    a = -1.;
    b = 2.;
    double ex_ans = 6. * ( std::atan(2.) - std::atan(-1.));
    ns = ns0;
    double t2 = trap_integration<Fun2>(a,b,ns,f2);
    ns = ns0;
    double s2 = simp_integration<Fun2>(a,b,ns,f2);
    ns = ns0;
    double b2 = boole_integration<Fun2>(a,b,ns,f2);
    ns = ns0;
    double a2 = adaptive_integration<Fun2>(a,b,ns,f2,err,tol);
    double t2b = trap_integration<Fun2>(a,b,ns,f2);
    double s2a = adapt_simpsons_integrate<Fun2>(f2,a,b,tol,14);

    fprintf(stdout," exact = %20.10le\n",ex_ans);
    fprintf(stdout," trap  %12lu %20.10le\n",ns0,t2);
    fprintf(stdout," simp  %12lu %20.10le\n",ns0,s2);
    fprintf(stdout," bool  %12lu %20.10le\n",ns0,b2);
    fprintf(stdout," adapt %12lu %20.10le\n",ns,a2);
    fprintf(stdout," trap  %12lu %20.10le\n",ns,t2b);
    fprintf(stdout," asimp              %20.10le\n",s2a);
    fprintf(stdout,"\n");

    a = 0;
    b = 1.;
    ex_ans = -4.;
    ns = ns0;
    tol = 1.e-4;
    double t3 = trap_integration<Fun3>(a,b,ns,f3);
    double s3 = simp_integration<Fun3>(a,b,ns,f3);
    double b3 = boole_integration<Fun3>(a,b,ns,f3);
    double a3 = adaptive_integration<Fun3>(a,b,ns,f3,err,tol);
    double t3b = trap_integration<Fun3>(a,b,ns,f3);
    np = ns;
    ns = ns0;
    double a3_2 = adaptive_integration2<Fun3>(a,b,ns,f3,err,tol);
    double s3a = adapt_simpsons_integrate<Fun3>(f3,a,b,tol,14);

    fprintf(stdout," exact = %20.10le\n",ex_ans);
    fprintf(stdout," trap  %12lu %20.10le\n",ns0,t3);
    fprintf(stdout," simp  %12lu %20.10le\n",ns0,s3);
    fprintf(stdout," bool  %12lu %20.10le\n",ns0,b3);
    fprintf(stdout," adapt %12lu %20.10le\n",np,a3);
    fprintf(stdout," adap2 %12lu %20.10le\n",ns,a3_2);
    fprintf(stdout," trap  %12lu %20.10le\n",np,t3b);
    fprintf(stdout," asimp              %20.10le\n",s3a);
    fprintf(stdout,"\n");

    a= 1.e-2;
    b= 10.0;
    ns = 512;
    tol = 1.e-12;
    t1 = trap_integration<Fun4>(a,b,ns,f4);
    s1 = simp_integration<Fun4>(a,b,ns,f4);
    b1 = boole_integration<Fun4>(a,b,ns,f4);
    a1 = adaptive_integration<Fun4>(a,b,ns,f4,err,tol);
    t1b = trap_integration<Fun4>(a,b,ns,f4);
    np = ns;
    ns = ns0;
    a1_2 = adaptive_integration2<Fun4>(a,b,ns,f4,err,tol);
    s1a = adapt_simpsons_integrate<Fun4>(f4,a,b,tol,14);

    fprintf(stdout,"scattering prob\n");
    fprintf(stdout," trap  %12lu %20.10le\n",ns0,t1);
    fprintf(stdout," simp  %12lu %20.10le\n",ns0,s1);
    fprintf(stdout," bool  %12lu %20.10le\n",ns0,b1);
    fprintf(stdout," adapt %12lu %20.10le\n",np,a1);
    fprintf(stdout," adap2 %12lu %20.10le\n",ns,a1_2);
    fprintf(stdout," trap  %12lu %20.10le\n",np,t1b);
    fprintf(stdout," asimp              %20.10le\n",s1a);
    fprintf(stdout,"\n");

}
