#include "network_th.h"
#include <mass_th.h>
network_K_N::network_K_N(int Ni, double epsi, double betai, double ai, double di, double Ji)
{
    N = Ni; eps = epsi; beta = betai; a = ai; d = di; J = Ji;
    OMEGA.resize(N, N); OMEGA.zero();
    save.resize(N); save.zero();
    E.resize(N); E.zero();
    BPhi.resize(N); BPhi.zero();
    IPSC.resize(N); IPSC.zero();
    JD.resize(N); JD.zero();
    JX.resize(N); JX.zero();
    cd.resize(N); cd.zero();
    Pinv.resize(N); Pinv.zero();
    h.resize(N); h.zero();
    hr.resize(N); hr.zero();
    r.resize(N); r.zero();
    I.resize(N); I.zero();
    v.resize(N); v.zero();
    epsvar.resize(N); epsvar.zero();

} 

void network_K_N::intialisation_FORCE_method(double Gi, double Qi, double pi, double lambdai)
{
    G = Gi; Q = Qi; p = pi; lambda = lambdai;
    E.random(-1, 1);
    E *= Q;
    save.zero();
    BPhi.zero();
    IPSC.zero();
    h.zero();
    r.zero();
    hr.zero();
    cd.zero();
    JD.zero();
    Pinv.eye();			  
    Pinv *= lambda; 
    v.random(0, .6); // 
    epsvar.random(-.002, .002);
    I.zero();
}
void network_K_N::display_parameters()
{
    cout 
    << "PARAMETRS NETWORK\n"
    << "_________________\n"
    << "N = " << N
    << ", eps = " << eps
    << ", beta = " << beta 
    << ", a = " << a 
    << ", d = " << d
    << ", J = " << J
    << '\n';
}