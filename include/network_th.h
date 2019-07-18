#ifndef NETWORK_TH
#define NETWORK_TH
#include <mass_th.h>
class network_K_N
{
public:
// Network parameters
    int N = 0;
    double eps = 0, beta = 0, a = 0, d = 0, J = 0, eps_peak = 0;
// FORCE method parameters
    double G = 0, Q = 0, vpeak = 0, p = 0, lambda = 0, M, M1;
// work time
     int nt = 0, imin = 0, icrit = 0, step_iteration = 300;
// nt - full time; imin - start learning; icrit - stop learning
// input 
    mass_th f_in;
    int dim_f_in = 0;
// output
    mass_th f_out;
    int dim_f_out;
// Network architecture
    mass_th OMEGA;
	mass_th save;	//save memory
	mass_th index;   //special (not work (delete!!!))
	mass_th E;		// Encoder
	mass_th BPhi;	// Decoder output
	mass_th BPhiin; // Decoder input
	mass_th IPSC;	// current 1
	mass_th JD;		// current 2
	mass_th JX;		// full cerrent
	mass_th cd;		// minimization error
	mass_th Pinv; // RLS - method
	mass_th h;      
	mass_th hr;
	mass_th r;
	// зададим переменные модели (напряжение и ток)
	mass_th I;     // variable Network
	mass_th v;     
	mass_th epsvar;
	mass_th save_error; // save error begin f_out and network output 
	string eps_setting;
	network_K_N()
	{}
    network_K_N(int N, double eps, double beta, double a, double d, double J);
	void reinitialisation(int N, double eps, double beta, double a, 
	double d, double J);
    void intialisation_FORCE_method(double G, double Q, double p, double lambda, double M, double M1, string eps_setting, double vpeak);
	void special_omega_rand();
    void get_f_in(mass_th f_in);
    void get_f_out(mass_th f_out);
	void get_time( int,  int,  int, int);
    void display_parameters();
	void FORCE_learning();
	
	void model_neuron_random_eps(mass_th &v, mass_th &I, mass_th &synaptic, int N, double a, double d, double eps,
				  double beta, double J, mass_th &eps_va, double peak, double eps_peak);
	void model_neuron_random_eps_static(mass_th &v, mass_th &I, mass_th &synaptic, int N, double a, double d, double eps,
				  double beta, double J, mass_th &eps_va, double peak);
	void model_neuron_no_random(mass_th &v, mass_th &I, mass_th &synaptic, int N, double a, double d, double eps,
				  double beta, double J, double peak);
	void model_synaps(mass_th &h, mass_th &r, mass_th &hr,
				  mass_th &v, double &vpeak, mass_th &ISPC, mass_th &JD, int N, double M, double M1);
	mass_th ret_error();
};
#endif
