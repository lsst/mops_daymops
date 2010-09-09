#include <vector>
#include <algorithm>


/* jmyers sep 2010 
   
   copy-pasting Kubica's method for calculating quadratic fits,
   replacing dyv_ stuff with C++ style vectors.  use this compare
   performance with GSL.  It's probably quite a bit faster.

   see linkTracklets/obs.c line 651 for the original.

*/


/* Computes the coefficients for the motion equation
   (WITHOUT too much computation... hopefully):

   X = M_0 + M_1 * t + M_2 * 0.5 * t^2

   will default to M_i = 0.0 if i > floor(log2(# points)) + 1
*/

#define NUM_FOR_QUAD 4


std::vector<double> 
mk_fill_obs_moments(std::vector<double> X, std::vector<double> time, 
		    unsigned int num_coefficients) 
{
    double A, B, C, D, E, F, G;
    double a, b, c, t, x;
    double bot, w, sum;
    double dobN, tspread;
    std::vector<double> res;
    int i, Nv;
    unsigned int N;

    /* Count the number of observations and the number of virtual
       observations (i.e. the number of distinct time steps) */
    N    = X.size();
    // Nv   = dyv_count_num_unique(time, 1e-6);
    // jmyers: we will trust that there are no redundant items per time
    Nv = N;
    
    dobN = (double)N;
    tspread = std::max_element( time.begin(), time.end() ) -
        std::min_element( time.begin(), time.end() );
     

    // TODO: should but an exception here  my_assert(N > 0);
    res.resize(num_coefficients, 0.0);
    
    /* There are 4 cases:  
       1) Nv <= 0     -> Return all zeros
       2) Nv == 1     -> Return that point with 0 vel/accel
       3) 1 < Nv < 4  -> Compute POS and VEL
       4) Nv >= 4     -> Compute POS, VELL, and ACCEL
    */

    /* Very quickly filter out the 1 and 2 data points case */
    if(Nv == 1) {
        sum = 0.0;
        bot = 0.0;
        for(i=0;i<N;i++) {
            w = 1.0;
            //sum += (dyv_ref(X,i) * w);
            sum += X.at(i) * w;
            bot += w;
        }
        //dyv_set(res,0,sum/bot);
        res.at(0) = sum/bot;
    } else {
        A = 0.0; B = 0.0; C = 0.0;
        D = 0.0; E = 0.0; F = 0.0;
        G = 0.0;
	  
        for(i=0;i<N;i++) {
            w = 1.0;	       
            t  = time.at(i); //dyv_ref(time,i);
            x  = X.at(i); //dyv_ref(X,i);
            A += ((t*t*t*t)*w);
            B += ((t*t*t)*w);
            C += ((t*t)*w);
            D += ((t)*w);
            E += ((x*t*t)*w);
            F += ((x*t)*w);
            G += ((x)*w);
        }

        /* Default to linear for a few points or */
        /* A very short arc (< 2 hours).         */
        if ((Nv < NUM_FOR_QUAD)||(tspread < 0.1)) {
            bot = D*D - C*dobN;
	       
            a = 0.0;    
            b = (G*D - dobN*F)/bot;
            c = (F*D - C*G)/bot;
        } else {
            bot  = A*D*D - A*dobN*C + dobN*B*B;
            bot += C*C*C - 2.0*C*B*D;
	       
            a  = C*C*G - G*B*D + D*D*E;
            a += B*dobN*F-dobN*E*C - D*F*C;
            a  = 2.0 * (a/bot);
	       
            b  = D*A*G - D*E*C - B*G*C;
            b += B*dobN*E + F*C*C -dobN*A*F;
            b  = (b/bot);
	       
            c  = E*C*C - E*B*D - A*G*C + A*D*F;
            c += B*B*G - C*B*F;
            c  = (c/bot);
        }
	  
        res.at(0) = c; //dyv_set(res,0,c);
        res.at(1) = b; //dyv_set(res,1,b);
        res.at(2) = a; //dyv_set(res,2,a);
    }

    return res;
}
