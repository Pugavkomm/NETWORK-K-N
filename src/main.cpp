#include <mass_th.h>
#include <network_th.h>
#include <fstream>
using namespace std;


void teaher_f(mass_th &x, double dt, int nt, double hz , 
                            string type, double hz2, int &dim_f_out);
void read_teacher(mass_th &x, int &nt, int &dim_f_out, 
                            string name_file, int &imin, int &icrit);

void read_f_in(mass_th &x, int nt, int &dim_f_in, string name_f_in);

void display_standart_and_write_file(string name_file, string eps_setting,
string name_programm, string name_teacher, int imin, int icrit, int nt, double vpeak,
double eps_peak, double J, double eps, double a, double d, double beta, double G,
double Q, double lambda, double p, double M, double M1, int N);

void display_GQmap_and_write_file(string name_file, string eps_setting,
string name_programm, string name_teacher, int imin, int icrit, int nt, double vpeak,
double eps_peak, double J, double eps, double a, double d, double beta, double Gstart,
double Gstop, double Gstep, double Qstart, double Qstop, double Qstep, double lambda, 
double p, double M, double M1, int N);

int main(int argc, char* argv[])
{
    string name_loader;
    string name_teacher;
    string name_f_in;
    mass_th f_in;
    int dim_f_in = 0;
    cout << argc << endl;
    for (int i = 0; i < argc; i++)
        cout << argv[i] << '\n';
    
    name_loader = argv[1];
    name_teacher = argv[2];
    
    
    cout << "name_loader = " << name_loader << '\n'
    << "name_teacher = " << name_teacher << '\n';
    cout << "Start test programm" << endl;
    string name_programm, type_teacher, eps_setting;
    int N, dim_f_out;
    double vpeak, eps_peak, J, eps, a, d, beta, G, 
    Q, lambda, T, tcrit, tmin, dt, hz1, hz2, M, M1, p;
    ifstream teach(name_teacher);
    teach >> T; teach >> tmin; teach >> tcrit;
    teach.close();
    int imin = round(tmin / dt);
	int icrit = round(tcrit / dt);
	int nt = round(T / dt); // количество итераций
    mass_th teacher;
    read_teacher(teacher, nt, dim_f_out, name_teacher, imin, icrit);
    if (argc == 4)
    {
        name_f_in = argv[3];
        read_f_in(f_in, nt, dim_f_in, name_f_in);
    }
    ifstream in(name_loader);
    getline(in, name_programm);
    getline(in, eps_setting);
    cout << "name_program:__ " << name_programm << endl;
    
    if (name_programm == "standart")
    {

    in >> N; in >> vpeak; in >> eps_peak; 
    in >> J; in >> eps; in >> a; in >> d;
    in >> beta; in >> G; in >> Q; 
    in >> p; in >> M; in >> M1; in >> lambda;
    in.close();
    display_standart_and_write_file("NONE", eps_setting, name_programm, name_teacher,
    imin, icrit, nt, vpeak, eps_peak, J, eps, a, d, beta, G, Q, lambda, p, M, M1, N);
    network_K_N model(N, eps, beta, a, d, J);
   
    model.get_f_out(teacher);
    if (dim_f_in > 0)
        model.get_f_in(f_in);
    model.intialisation_FORCE_method(G, Q, p, lambda, M, M1, eps_setting, vpeak);
    model.get_time(nt, imin, icrit,  300);
    model.display_parameters();
    model.FORCE_learning();
    ofstream out;
    out.open("error");
   
    for (int i = 0; i < nt; i++)
    {
        out << model.ret_error().matrix[0][i] << '\n';
    }
    out.close();
    //outer.display_m();
    }

    if (name_programm == "GQmap")
    {
        double Gstart, Gstop, Gstep, Qstart, Qstop, Qstep;
    
        in >> N; in >> vpeak; in >> eps_peak; 
        in >> J; in >> eps; in >> a; in >> d;
        in >> beta; in >> Gstart; in >> Gstop; in >> Gstep;
        in >> Qstart; in >> Qstop; in >> Qstep; 
        in >> p; in >> M; in >> M1; in >> lambda;
        in.close();
        cout 
        << "type = " << name_teacher << '\n'
        << "eps_setting = " << eps_setting << '\n'
        << "name_program = " << name_programm << '\n'
        << "N = " << N << '\n'
        << "vpeak = " << vpeak << '\n'
        << "eps_peak = " << eps_peak << '\n'
        << "J = " << J << '\n'
        << "eps = " << eps << '\n'
        << "a = " << a << '\n'
        << "d = " << d << '\n'
        << "beta = " << beta << '\n'
        << "Gstart = " << Gstart << '\n'
        << "Gstop = " << Gstop << '\n'
        << "Gstep = " << Gstep << '\n'
        << "Qstart = " << Qstart << '\n'
        << "Qstop = " << Qstop << '\n'
        << "Qstep = " << Qstep << '\n'
        << "lambda = " << lambda << '\n'
        << "p = " << p << '\n'
        << "M = " << M << '\n'
        << "M1 = " << M1 << '\n';
      int Gsteps = round((Gstop - G) / Gstep),
        Qsteps = round((Qstop -Q) / Qstep);
        cout << "Gsteps = " << Gsteps << '\n'
        << "Qsteps = " << Qsteps << '\n'
        << "fullStep = " << (Qsteps * Gsteps) << '\n';
        double L2[dim_f_out];
        ofstream GQmatrix;
        ofstream out;
        display_GQmap_and_write_file("filesaveGQmap.txt", eps_setting, name_programm,
        name_teacher, imin, icrit, nt, vpeak, eps_peak, J, eps, a, d, beta, 
        Gstart, Gstop, Gstep, Qstart, Qstop, Qstep, lambda, p, M, M1, N);
        for (double ig = Gstart; ig < Gstop; ig += Gstep)
        {
            G = ig;
            GQmatrix.open("GQmatrix", ios::app);
            GQmatrix << G << ':' << ' ';
            GQmatrix.close();
            for (double iq = Qstart; iq < Qstop; iq += Qstep)
            {
                for (int i = 0; i < dim_f_out; i++)
                    L2[i] = 0;
                Q = iq;
                cout << "G = " << G << ", Q = " << Q << '\n';
                network_K_N model(N, eps, beta, a, d, J);
                model.get_f_out(teacher);
                if (dim_f_in > 0)
                    model.get_f_in(f_in);
                model.intialisation_FORCE_method(G, Q, p, lambda, M, M1, eps_setting, vpeak);
                model.get_time(nt, imin, icrit,  300);
                model.display_parameters();
                model.FORCE_learning();
                for (int i = icrit; i < nt; i++)
                    for (int j = 0; j < dim_f_out; j++)
                        {
                            L2[j] += model.ret_error().matrix[0][i] * model.ret_error().matrix[0][i];
                        }
                out.open("filesaveGQmap.txt", ios::app);
                for (int i = 0; i < dim_f_out; i++)
                    out << L2[i] << ' ';
                out.close();
                GQmatrix.open("GQmatrix", ios::app);
                GQmatrix << '|' << Q << '|' << ' ';
                GQmatrix.close();
            }
            GQmatrix.open("GQmatrix", ios::app);
            GQmatrix << '\n';
            GQmatrix.close();
            out.open("filesaveGQmap.txt", ios::app);
            out << '\n';
            out.close();
        }
    }

    if (name_programm == "N_size")
    {
        int Nstart, Nstop, Nstep, ex;
        in >> Nstart; in >> Nstop; in >> Nstep;
        in >> ex; in >> vpeak; in >> eps_peak;  in >> J; 
        in >> eps; in >> a; in >> d;
        in >> beta; in >> G; in >> Q;
        in >> p; in >> M; in >> M1; 
        in >> lambda;
        in.close();
        cout 
        << "type = " << name_teacher << '\n'
        << "dt = " << dt << '\n'
        << "T = " << T << '\n'
        << "tmin = " << tmin << '\n'
        << "tcrit  = " << tcrit  << '\n'
        << "hz1 = " << hz1 << '\n'
        << "hz2 = " << hz2 << '\n';
        cout 
        << "eps_setting = " << eps_setting << '\n'
        << "name_program = " << name_programm << '\n'
        << "Nstart = " << Nstart << '\n'
        << "Nstop = " << Nstop << '\n'
        << "Nstep = " << Nstep << '\n'
        << "ex = " << ex << '\n'
        << "vpeak = " << vpeak << '\n'
        << "eps_peak = " << eps_peak << '\n'
        << "J = " << J << '\n'
        << "eps = " << eps << '\n'
        << "a = " << a << '\n'
        << "d = " << d << '\n'
        << "beta = " << beta << '\n'
        << "lambda = " << lambda << '\n'
        << "p = " << p << '\n'
        << "M = " << M << '\n'
        << "M1 = " << M1 << '\n';
        int Nsteps = (Nstop - Nstart) / Nstep;               
        double L2[dim_f_out];
        ofstream GQmatrix;
        ofstream out;
        out.open("filesaveGQmap.txt", ios::app);
        out 
        << "\ntype = " << name_teacher << '\n'
        << "dt = " << dt << '\n'
        << "T = " << T << '\n'
        << "tmin = " << tmin << '\n'
        << "tcrit  = " << tcrit  << '\n'
        << "hz1 = " << hz1 << '\n'
        << "hz2 = " << hz2 << '\n';
        out 
        << "eps_setting = " << eps_setting << '\n'
        << "name_program = " << name_programm << '\n'
        << "N = " << N << '\n'
        << "vpeak = " << vpeak << '\n'
        << "eps_peak = " << eps_peak << '\n'
        << "J = " << J << '\n'
        << "eps = " << eps << '\n'
        << "a = " << a << '\n'
        << "d = " << d << '\n'
        << "beta = " << beta << '\n'
        << "G = " << G << '\n'
        << "Q = " << Q << '\n'
        << "lambda = " << lambda << '\n'
        << "p = " << p << '\n'
        << "M = " << M << '\n'
        << "M1 = " << M1 << '\n';
        out << "___________________________________\n";
        out.close();
        GQmatrix.open("GQmatrix", ios::app);
        GQmatrix 
        << "\ntype = " << name_teacher << '\n'
        << "dt = " << dt << '\n'
        << "T = " << T << '\n'
        << "tmin = " << tmin << '\n'
        << "tcrit  = " << tcrit  << '\n'
        << "hz1 = " << hz1 << '\n'
        << "hz2 = " << hz2 << '\n';
        GQmatrix 
        << "eps_setting = " << eps_setting << '\n'
        << "name_program = " << name_programm << '\n'
        << "N = " << N << '\n'
        << "vpeak = " << vpeak << '\n'
        << "eps_peak = " << eps_peak << '\n'
        << "J = " << J << '\n'
        << "eps = " << eps << '\n'
        << "a = " << a << '\n'
        << "d = " << d << '\n'
        << "beta = " << beta << '\n'
        //<< "Gstart = " << Gstart << '\n'
        //<< "Gstop = " << Gstop << '\n'
        //<< "Gstep = " << Gstep << '\n'
        //<< "Qstart = " << Qstart << '\n'
        //<< "Qstop = " << Qstop << '\n'
        //<< "Qstep = " << Qstep << '\n'
        << "lambda = " << lambda << '\n'
        << "p = " << p << '\n'
        << "M = " << M << '\n'
        << "M1 = " << M1 << '\n';
        GQmatrix << "___________________________________\n";
        GQmatrix << "G1:|Q1|Q2|Q3|Q4|...\nG2:|..|..|..|..|...\nG3:|..|..|..|..|...\n";
        GQmatrix << "...................\n";
        GQmatrix.close();
        //for (double ig = Gstart; ig < Gstop; ig += Gstep)
        //{
        //    G = ig;
        //    GQmatrix.open("GQmatrix", ios::app);
        //    GQmatrix << G << ':' << ' ';
        //    GQmatrix.close();
        //    for (double iq = Qstart; iq < Qstop; iq += Qstep)
        //    {
        //        for (int i = 0; i < dim_f_out; i++)
        //            L2[i] = 0;
        //        Q = iq;
        //        cout << "G = " << G << ", Q = " << Q << '\n';
        //        network_K_N model(N, eps, beta, a, d, J);
        //        model.get_f_out(teacher);
        //        if (dim_f_in > 0)
        //            model.get_f_in(f_in);
        //        model.intialisation_FORCE_method(G, Q, p, lambda, M, M1, eps_setting, vpeak);
        //        model.get_time(nt, imin, icrit,  300);
        //        model.display_parameters();
        //        model.FORCE_learning();
        //        for (int i = icrit; i < nt; i++)
        //            for (int j = 0; j < dim_f_out; j++)
        //                {
        //                    L2[j] += model.ret_error().matrix[0][i] * model.ret_error().matrix[0][i];
        //                }
        //        out.open("filesaveGQmap.txt", ios::app);
        //        for (int i = 0; i < dim_f_out; i++)
        //            out << L2[i] << ' ';
        //        out.close();
        //        GQmatrix.open("GQmatrix", ios::app);
        //        GQmatrix << '|' << Q << '|' << ' ';
        //        GQmatrix.close();
        //    }
        //    GQmatrix.open("GQmatrix", ios::app);
        //    GQmatrix << '\n';
        //    GQmatrix.close();
        //    out.open("filesaveGQmap.txt", ios::app);
        //    out << '\n';
        //    out.close();
        //}
    }

    return 0;
}


void teaher_f(mass_th &x, double dt, int nt, double hz , string type, double hz2, int &dim_f_out)
{
    dim_f_out = 1;
	for (int i = 0; i < nt; i++)
	{
		if (type == "sin")
			x.matrix[dim_f_out - 1][i] = sin(hz * 2 * M_PI * i * dt) ;
		else if(type == "saw")
			x.matrix[dim_f_out - 1][i] = (2/M_PI)*asin(sin(hz * 2 * M_PI * dt * i));
		else if (type == "2sin")
			x.matrix[dim_f_out - 1][i] = sin(hz * 2 * M_PI * i * dt)  * sin(hz2 * 2 * M_PI * i * dt) ;

	}
}

void read_teacher(mass_th &x, int &nt, int &dim_f_out, string name_file, int &imin, int &icrit)
{
    int n, m;
    std::ifstream file(name_file);
    if (!file)
    {
        cout << "File couldn't be open.\n";
    }
    file >> n;
    file >> m;
    file >> imin;
    file >> icrit;
    nt = n;
    dim_f_out = m;
    cout << "nt = " << nt << ", dim_teach = " << m << ", imin = " << imin << "icrit = " << icrit << '\n';
    x.resize(dim_f_out, nt);
    for (int i = 0; i < nt; ++i)
    {
        for (int j = 0; j < dim_f_out; ++j)
        {
            file >> x.matrix[j][i];
        }
    }
    file.close();
}

void read_f_in(mass_th &x, int nt, int &dim_f_in, string name_f_in)
{
    ifstream file(name_f_in);
    file >> dim_f_in;
    x.resize(dim_f_in, nt);
    for (int i = 0; i < nt; ++i)
        for (int j = 0; j < dim_f_in; ++j)
        {
            file >> x.matrix[j][i];
            
        }
    file.close();
}

void display_standart_and_write_file(string name_file, string eps_setting,
string name_programm, string name_teacher, int imin, int icrit, int nt, double vpeak,
double eps_peak, double J, double eps, double a, double d, double beta, double G,
double Q, double lambda, double p, double M, double M1, int N)
{
    cout 
    << "type = " << name_teacher << '\n'
    << "eps_setting = " << eps_setting << '\n'
    << "name_program = " << name_programm << '\n'
    << "N = " << N << '\n'
    << "vpeak = " << vpeak << '\n'
    << "eps_peak = " << eps_peak << '\n'
    << "J = " << J << '\n'
    << "eps = " << eps << '\n'
    << "a = " << a << '\n'
    << "d = " << d << '\n'
    << "beta = " << beta << '\n'
    << "G = " << G << '\n'
    << "Q = " << Q << '\n'
    << "lambda = " << lambda << '\n'
    << "p = " << p << '\n'
    << "M = " << M << '\n'
    << "M1 = " << M1 << '\n';
}

void display_GQmap_and_write_file(string name_file, string eps_setting,
string name_programm, string name_teacher, int imin, int icrit, int nt, double vpeak,
double eps_peak, double J, double eps, double a, double d, double beta, double Gstart,
double Gstop, double Gstep, double Qstart, double Qstop, double Qstep, double lambda, 
double p, double M, double M1, int N)
{
    cout 
        << "type = " << name_teacher << '\n'
        << "eps_setting = " << eps_setting << '\n'
        << "name_program = " << name_programm << '\n'
        << "N = " << N << '\n'
        << "vpeak = " << vpeak << '\n'
        << "eps_peak = " << eps_peak << '\n'
        << "J = " << J << '\n'
        << "eps = " << eps << '\n'
        << "a = " << a << '\n'
        << "d = " << d << '\n'
        << "beta = " << beta << '\n'
        << "Gstart = " << Gstart << '\n'
        << "Gstop = " << Gstop << '\n'
        << "Gstep = " << Gstep << '\n'
        << "Qstart = " << Qstart << '\n'
        << "Qstop = " << Qstop << '\n'
        << "Qstep = " << Qstep << '\n'
        << "lambda = " << lambda << '\n'
        << "p = " << p << '\n'
        << "M = " << M << '\n'
        << "M1 = " << M1 << '\n';
        ofstream file(name_file, ios::app);
        file
        << "\ntype = " << name_teacher << '\n'
        << "eps_setting = " << eps_setting << '\n'
        << "name_program = " << name_programm << '\n'
        << "N = " << N << '\n'
        << "vpeak = " << vpeak << '\n'
        << "eps_peak = " << eps_peak << '\n'
        << "J = " << J << '\n'
        << "eps = " << eps << '\n'
        << "a = " << a << '\n'
        << "d = " << d << '\n'
        << "beta = " << beta << '\n'
        << "Gstart = " << Gstart << '\n'
        << "Gstop = " << Gstop << '\n'
        << "Gstep = " << Gstep << '\n'
        << "Qstart = " << Qstart << '\n'
        << "Qstop = " << Qstop << '\n'
        << "Qstep = " << Qstep << '\n'
        << "lambda = " << lambda << '\n'
        << "p = " << p << '\n'
        << "M = " << M << '\n'
        << "M1 = " << M1 << '\n';
        file << "___________________________________\n";
        file.close();

        file.open("GQmatrix", ios::app);
        file 
        << "\ntype = " << name_teacher << '\n'
        << "eps_setting = " << eps_setting << '\n'
        << "name_program = " << name_programm << '\n'
        << "N = " << N << '\n'
        << "vpeak = " << vpeak << '\n'
        << "eps_peak = " << eps_peak << '\n'
        << "J = " << J << '\n'
        << "eps = " << eps << '\n'
        << "a = " << a << '\n'
        << "d = " << d << '\n'
        << "beta = " << beta << '\n'
        << "Gstart = " << Gstart << '\n'
        << "Gstop = " << Gstop << '\n'
        << "Gstep = " << Gstep << '\n'
        << "Qstart = " << Qstart << '\n'
        << "Qstop = " << Qstop << '\n'
        << "Qstep = " << Qstep << '\n'
        << "lambda = " << lambda << '\n'
        << "p = " << p << '\n'
        << "M = " << M << '\n'
        << "M1 = " << M1 << '\n';
        file << "___________________________________\n";
        file << "G1:|Q1|Q2|Q3|Q4|...\nG2:|..|..|..|..|...\nG3:|..|..|..|..|...\n";
        file << "...................\n";
        file.close();
}