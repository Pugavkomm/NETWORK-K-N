#include <mass_th.h>
#include <network_th.h>
using namespace std;


int main()
{
    //int N = 5;
    //double a, beta, d, eps, J;
	//eps = .0015; beta = .04; d = .5; a = .25; J = .12;
    //network_K_N test(N, eps, beta, a, d, J);
    //test.display_parameters();
    mass_th test(5, "test1");
    mass_th test1(10, 10, "test2");
    test.random(0, 1);
    test1.random(0, 1);
    test.display_m(); test1.display_m();
    test = test1;
    test1.random(0, 1);
    test.display_m(); test1.display_m();



    return 0;
}
