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
    Pinv.resize(N, N); Pinv.zero();
    h.resize(N); h.zero();
    hr.resize(N); hr.zero();
    r.resize(N); r.zero();
    I.resize(N); I.zero();
    v.resize(N); v.zero();
    epsvar.resize(N); epsvar.zero();

} 
void network_K_N::reinitialisation(int Ni, double epsi, double betai, double ai, double di, double Ji)
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
    Pinv.resize(N, N); Pinv.zero();
    h.resize(N); h.zero();
    hr.resize(N); hr.zero();
    r.resize(N); r.zero();
    I.resize(N); I.zero();
    v.resize(N); v.zero();
    epsvar.resize(N); epsvar.zero();
}
void network_K_N::intialisation_FORCE_method(double Gi, double Qi, 
double pi, double lambdai, double Mi, double M1i, string eps_settingi, double vpeaki)
{
    G = Gi; Q = Qi; p = pi; lambda = lambdai; M = Mi; M1 = M1i; eps_setting = eps_settingi;
    vpeak = vpeaki;
    OMEGA.zero();
    special_omega_rand();
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
    v.random(0, .6); // add parameter in future
    epsvar.random(-.002, .002); // add parameter in future
    I.zero();
}
void network_K_N::special_omega_rand()
{
	default_random_engine generator;
	normal_distribution<double> distribution(0.0, 1.0);
	for (int i = 0; i < OMEGA.row; i++)
	{
		for (int j = 0; j < OMEGA.col; j++)
		{
			OMEGA.matrix[i][j] = G * distribution(generator) * ((double)rand() / RAND_MAX < p);

			OMEGA.matrix[i][j] /= (sqrt((double)OMEGA.row) * p);
		}
	}
}
void network_K_N::get_f_out(mass_th f_outi)
{
    f_out = f_outi;
    dim_f_out = f_out.row;
    E.resize(N, dim_f_out);
    BPhi.resize(N, dim_f_out);

}

void network_K_N::get_f_in(mass_th f_ini)
{
    f_in = f_ini;
    dim_f_in = f_in.row;
    BPhiin.resize(N, dim_f_in);
}
void network_K_N::get_time( int nti,  int imini,  int icriti, int step_interationi)
{
    nt = nti;
    imin = imini;
    icrit = icriti;
    step_iteration = step_interationi;
    save_error.resize(dim_f_out, nt);
}
void network_K_N::model_neuron_random_eps(mass_th &v, mass_th &I, mass_th &synaptic, int N, double a, double d, double eps,
				  double beta, double J, mass_th &eps_va, double peak, double eps_peak)
{
	int i;
	double I_save;
#pragma omp parallel shared(v, I, beta, eps, J, a, d, synaptic, N, eps_va, peak) private(i, I_save)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < N; i++)
		{
			I_save = I.matrix[i][0];
			I.matrix[i][0] += ( eps + eps_va.matrix[i][0] * (v.matrix[i][0] > eps_peak)) * (v.matrix[i][0] - J);
			v.matrix[i][0] += v.matrix[i][0] *( (v.matrix[i][0] - a) * (1 - v.matrix[i][0])) -
							 beta * (v.matrix[i][0] > d) - I_save + synaptic.matrix[i][0];
		}
	}
}

void network_K_N::model_neuron_random_eps_static(mass_th &v, mass_th &I, mass_th &synaptic, int N, double a, double d, double eps,
				  double beta, double J, mass_th &eps_va, double peak)
{
	int i;
	double I_save;
#pragma omp parallel shared(v, I, beta, eps, J, a, d, synaptic, N, eps_va, peak) private(i, I_save)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < N; i++)
		{
			I_save = I.matrix[i][0];
			I.matrix[i][0] += (eps + eps_va.matrix[i][0]) * (v.matrix[i][0] - J);
			v.matrix[i][0] += v.matrix[i][0] *( (v.matrix[i][0] - a) * (1 - v.matrix[i][0])) -
							 beta * (v.matrix[i][0] > d) - I_save + synaptic.matrix[i][0];
		}
	}
}

void network_K_N::model_neuron_no_random(mass_th &v, mass_th &I, mass_th &synaptic, int N, double a, double d, double eps,
				  double beta, double J, double peak)
{
	int i;
	double I_save;
#pragma omp parallel shared(v, I, beta, eps, J, a, d, synaptic, N, peak) private(i, I_save)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < N; i++)
		{
			I_save = I.matrix[i][0];
			I.matrix[i][0] += eps * (v.matrix[i][0] - J);
			v.matrix[i][0] += v.matrix[i][0] *( (v.matrix[i][0] - a) * (1 - v.matrix[i][0])) -
							 beta * (v.matrix[i][0] > d) - I_save + synaptic.matrix[i][0];
		}
	}
}
void network_K_N::model_synaps(mass_th &h, mass_th &r, mass_th &hr,
				  mass_th &v, double &vpeak, mass_th &ISPC, mass_th &JD, int N, double M, double M1)
{
#pragma omp parallel shared(h, r, hr, v, vpeak, ISPC, JD,  N, M, M1)
	{
#pragma omp for schedule(dynamic)
		for (int i = 0; i < N; i++)
		{
			r.matrix[i][0] += -r.matrix[i][0] / M + hr.matrix[i][0];
			hr.matrix[i][0] += -hr.matrix[i][0] / M1 + (v.matrix[i][0] > vpeak) / (M * M1);
		}
	}
}
void network_K_N::FORCE_learning()
{
    double start_time;
    mass_th error(dim_f_out, 1);
    mass_th zout(dim_f_out, 1);
    
    for (int step_system = 0; step_system < nt; step_system++)
    {               
                if (dim_f_in == 0)
				    JX = E * zout + IPSC;
                else   
                    JX = E * zout + IPSC + BPhiin * f_in.slice(0, dim_f_in, step_system, step_system + 1);
                if (eps_setting == "random")
					model_neuron_random_eps(v, I, JX, N, a, d, eps, beta, J, epsvar, vpeak, eps_peak);
				else if(eps_setting == "static")
					model_neuron_random_eps_static(v, I, JX, N, a, d, eps, beta, J, epsvar, vpeak);
				else if(eps_setting == "no_random")
					model_neuron_no_random(v, I, JX, N, a, d, eps, beta, J, vpeak);
				
				model_synaps(h, r, hr, v, vpeak, IPSC, JD, N, M, M1);
                
				IPSC = OMEGA * r;
                
				zout = (BPhi--) * r;
				error =  zout - f_out.slice(0, dim_f_out, step_system, step_system + 1);
                #pragma omp parallel
                #pragma omp for
                for (int save_err = 0; save_err < dim_f_out; save_err++)
                    save_error.matrix[save_err][step_system] = error.matrix[save_err][0];
				if (step_system > imin && step_system < icrit)
				{
                    cd = Pinv * r;
                    BPhi -= (cd * error);
                    

					Pinv -= (cd * cd--) / (1.0 + (r-- * cd).ret(0, 0));
				}
				if (step_system % step_iteration == 0)
				{
					start_time = omp_get_wtime() - start_time;
                    for (int c = 0; c < dim_f_out; c++)
                        cout << "error" << c << " = " << f_out.matrix[c][step_system];
					cout << "step = " << step_system << ", осталось: " << (start_time * ((nt - step_system)) / 60. / step_iteration) << " минут" <<'\n';
					start_time = omp_get_wtime();
				}
    }   
}
mass_th network_K_N::ret_error()
{
    return save_error;
}
void network_K_N::display_parameters()
{
    cout 
    << "PARAMETERS NETWORK\n"
    << "_________________\n"
    << "N = " << N
    << ", eps = " << eps
    << ", beta = " << beta 
    << ", a = " << a 
    << ", d = " << d
    << ", J = " << J
    << '\n';
    cout 
    << "PARAMETERS ALGORITHM FORCE\n"
    << "_________________________\n"
    << "G = " << G
    << ", Q = " << Q
    << ", p = " << p
    << ", lambda = " << lambda
    << ", nt = " << nt
    << ", imin = " << imin
    << ", icrit = " << icrit
    << ",vpeak = " << vpeak 
    << '\n';
    cout 
    << "TIME:\n"
    << "full time: " << nt
    << ", start learning: " << imin 
    << ", end learning: " << icrit 
    << '\n'; 
}