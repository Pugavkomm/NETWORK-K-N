#ifndef NETWORK_TH
#define NETWORK_TH
#include <mass_th.h>
class network_K_N
{
public:
// Network parameters
    int N;
    double eps, beta, a, d, J, eps_peak;
// FORCE method parameters
    double G, Q, vpeak, p, lambda;
// work time
    unsigned int nt, imin, icrit;
// nt - full time; imin - start learning; icrit - stop learning
// input 
    mass_th f_in;
    unsigned int dim_f_in;
// output
    mass_th f_out;
    unsigned int dim_f_out;
// Network architecture
    mass_th OMEGA;
	mass_th save;	//save memory
	mass_th index;   //special (not work (delete!!!))
	mass_th E;		// Enkoder
	mass_th BPhi;	// Decoder
	mass_th IPSC;	// current 1
	mass_th JD;		// current 2
	mass_th JX;		// full cerrent
	mass_th cd;		// minimisation error
	mass_th Pinv; // RLS - method
	mass_th h;      
	mass_th hr;
	mass_th r;
	// зададим переменные модели (напряжение и ток)
	mass_th I;     // variable Network
	mass_th v;     
	mass_th epsvar;
    network_K_N(int N, double eps, double beta, double a, double d, double J);
    void intialisation_FORCE_method(double G, double Q, double p, double lambda);
    //void get_f_in(mass_th f_in);
    //void get_f_out(mass_th f_out);
    void display_parameters();
};
#endif
