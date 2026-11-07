// Historical benchmark only; production kernels are in GWLD-R/inst/include/gwld.
#include <iostream>
#include <fstream>
#include <sstream>
#include <armadillo>
#include <unordered_map>
#include <vector>
#include <string>
#include <float.h>
#include <omp.h>
#include "utils_reference.h"
using legacy_reference::allele_freq;
using legacy_reference::table;
using legacy_reference::writearma;

using namespace std;
using namespace arma;

double RMI(arma::Col<int> &g1, arma::Col<int> &g2) {
    // Get unique genotype values in ascending order.
    arma::Col<int> set_g1 = arma::unique(g1); 
    arma::Col<int> set_g2 = arma::unique(g2); 
    // Remove missing genotypes from the category set.
    arma::Col<int> row_names = set_g1.elem(find(set_g1 != -1));
    arma::Col<int> col_names = set_g2.elem(find(set_g2 != -1));
    int len = g1.n_elem;
    // Use floating-point category counts to avoid loss of precision.
    double R = row_names.n_elem, S = col_names.n_elem; 
    std::unordered_map<int, int> row_index;
    for (int i=0; i<R; i++) {
        row_index[row_names[i]] = i;
    }
    std::unordered_map<int, int> col_index;
    for (int j=0; j<S; j++) {
        col_index[col_names[j]] = j;
    }
    // Initialize the contingency table with double-precision zeros.
    arma::Mat<double> m(R, S, fill::zeros);
    for(int k=0; k<len; k++){
        if(g1[k]!= -1 && g2[k]!= -1){
            m(row_index[g1[k]], col_index[g2[k]]) +=1;
        }
    }
    double n = arma::accu(m);// Sum all contingency-table entries.
    arma::Col<double> a = arma::sum(m,1);
    arma::Row<double> b = arma::sum(m,0); // Armadillo dimension 1 sums rows; dimension 0 sums columns.
    double I = (lgamma(n+1) + arma::accu(arma::lgamma(m+1)) - arma::accu(arma::lgamma(a+1)) - arma::accu(arma::lgamma(b+1))) / (n*log(2));
    double w = n / (n + 0.5*R*S);
    arma::Col<double> x = (1-w)/R + w*a/n;
    arma::Row<double> y = (1-w)/S + w*b/n;
    double nu = (S+1)/(S*arma::accu(arma::square(x))) - 1/S;
    double mu = (R+1)/(R*arma::accu(arma::square(y))) - 1/R;
    double logOmega = (R-1)*(S-1)*log(n+0.5*R*S)+
    0.5*(R+nu-2)*arma::accu(arma::log(y))+
    0.5*(S+mu-2)*arma::accu(arma::log(x))+
    0.5*(lgamma(mu*R)+lgamma(nu*S)-R*(lgamma(S)+lgamma(mu))-S*(lgamma(R)+lgamma(nu)));
    double rmi = I-logOmega/(n*log(2));
    return rmi;
}


double MI(arma::Col<int> &g1, arma::Col<int> &g2) {
    // Get unique genotype values in ascending order.
    arma::Col<int> set_g1 = arma::unique(g1); 
    arma::Col<int> set_g2 = arma::unique(g2); 
    // Remove missing genotypes from the category set.
    arma::Col<int> row_names = set_g1.elem(find(set_g1 != -1));
    arma::Col<int> col_names = set_g2.elem(find(set_g2 != -1));
    int len = g1.n_elem;
    // Use floating-point category counts to avoid loss of precision.
    double R = row_names.n_elem, S = col_names.n_elem; 
    std::unordered_map<int, int> row_index;
    for (int i=0; i<R; i++) {
        row_index[row_names[i]] = i;
    }
    std::unordered_map<int, int> col_index;
    for (int j=0; j<S; j++) {
        col_index[col_names[j]] = j;
    }
    // Initialize the contingency table with double-precision zeros.
    arma::Mat<double> m(R, S, fill::zeros);
    for(int k=0; k<len; k++){
        if(g1[k]!= -1 && g2[k]!= -1){
            m(row_index[g1[k]], col_index[g2[k]]) +=1;
        }
    }
    double n = arma::accu(m); // Sum all contingency-table entries.
    arma::Mat<double> freq = m/n;// Convert counts to frequencies.
    arma::Col<double> a = arma::sum(freq,1); // Calculate row sums.
    arma::Row<double> b = arma::sum(freq,0); // Calculate column sums.
    // Calculate entropy in bits, defining 0 * log2(0) as zero.
    arma::Col<double> rowh = a.transform([] (double val) {return(val==0.0 ? 0.0 : -val*log2(val));});
    arma::Row<double> colh = b.transform([] (double val) {return(val==0.0 ? 0.0 : -val*log2(val));});
    arma::Mat<double> tabh = freq.transform([] (double val) {return(val==0.0 ? 0.0 : -val*log2(val));});
    // Calculate entropy at the first locus.
    double H_fg = arma::accu(rowh);
    // Calculate entropy at the second locus.
    double H_sg = arma::accu(colh);
    // Calculate joint entropy.
    double H_fg_sg = arma::accu(tabh);
    // Calculate mutual information.
    double I = H_fg + H_sg - H_fg_sg;
    return I;
}

//______________________________________LD_________________________________________________
double nlogLik(double &pAB, double &pA , double &pB, arma::Mat<double> &m) {
    // Evaluate the negative log-likelihood.
    double y = -((2*m(0,0)+m(0,1)+m(1,0))*log(pAB) +
                 (2*m(0,2)+m(0,1)+m(1,2))*log(pA-pAB) +
                 (2*m(2,0)+m(1,0)+m(2,1))*log(pB-pAB) +
                 (2*m(2,2)+m(2,1)+m(1,2))*log(1-pA-pB+pAB) + 
                m(1,1)*log(pAB*(1-pA-pB+pAB) + (pA-pAB)*(pB-pAB)));
    return y;
}

double fmin(double ax, double bx, double tol, double &pA, double &pB, arma::Mat<double> &matrix) {
	//Author:
	//    Original FORTRAN77 version by Richard Brent.
	//    C++ version by John Burkardt.
	//
	//  Reference:
	//
	//    Richard Brent,
	//    Algorithms for Minimization Without Derivatives,
	//    Dover, 2002,
	//    ISBN: 0-486-41998-3,
	//    LC: QA402.5.B74.
    /*  c is the squared inverse of the golden ratio */
    const double c = (3.0 - sqrt(5.0)) * 0.5;
    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
	/* eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);
    a = ax; b = bx;
    v = a + c * (b - a);
    w = v; x = v;
    d = 0.; e = 0.;/* -Wall */
    fx = nlogLik(x, pA, pB, matrix);
    fv = fx; fw = fx;
    tol3 = tol / 3.;
    /* main loop starts here ----------------------------------- */
    for(;;) {
	xm = (a + b) * .5;
	tol1 = eps * fabs(x) + tol3;
	t2 = tol1 * 2.;
	/* check stopping criterion */
	if (fabs(x - xm) <= t2 - (b - a) * .5) break;
	p = 0.; q = 0.; r = 0.;
	if (fabs(e) > tol1) { /* fit parabola */
	    r = (x - w) * (fx - fv);
	    q = (x - v) * (fx - fw);
	    p = (x - v) * q - (x - w) * r;
	    q = (q - r) * 2.;
	    if (q > 0.) p = -p; else q = -q;
	    r = e; e = d;
	}
	if (fabs(p) >= fabs(q * .5 * r) ||
	    p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */
	    if (x < xm) e = b - x; else e = a - x;
	    d = c * e;
	}
	else { /* a parabolic-interpolation step */
	    d = p / q; u = x + d;
	    /* f must not be evaluated too close to ax or bx */
	    if (u - a < t2 || b - u < t2) {
		d = tol1;
		if (x >= xm) d = -d;
	    }
	}
	/* f must not be evaluated too close to x */
	if (fabs(d) >= tol1)
	    u = x + d;
	else if (d > 0.)
	    u = x + tol1;
	else
	    u = x - tol1;
	fu = nlogLik(u, pA, pB, matrix);
	/* update  a, b, v, w, and x */
	if (fu <= fx) {
	    if (u < x) b = x; else a = x;
	    v = w;    w = x;   x = u;
	    fv = fw; fw = fx; fx = fu;
	} else {
	    if (u < x) a = u; else b = u;
	    if (fu <= fw || w == x) {
            v = w; fv = fw;
            w = u; fw = fu;
	    } else if (fu <= fv || v == x || v == w) {
		    v = u; fv = fu;
	    }
	}
    }
    /* end of main loop */
    return x;
}

double LD(arma::Col<int> &g1, arma::Col<int> &g2, std::string method="r^2") {
    // Calculate allele frequencies at the first locus.
    vector<double> A_freq =  allele_freq(g1);
    double pA = A_freq[0];
    double pa = A_freq[1];
    // Calculate allele frequencies at the second locus.
    vector<double> B_freq =  allele_freq(g2);
    double pB = B_freq[0];
    double pb = B_freq[1];
    // Initialize a 3-by-3 contingency table for the two loci.
    arma::Mat<double> m = table(g1, g2, 3, 3);
    double Dmin = std::max(-pA*pB, -pa*pb);
    double pmin = pA*pB + Dmin;
    double Dmax = std::min(pA*pb, pB*pa);
    double pmax = pA*pB + Dmax;

    double eps = DBL_EPSILON;
    double tol = pow(DBL_EPSILON, 0.25);// Set the optimizer tolerance.
    double pAB = fmin(pmin+eps, pmax+eps, tol, pA, pB, m);
    if(method == "r^2") {
        double corr = pow(pAB - pA*pB, 2)/(pA * pB * pa * pb);// Calculate r-squared.
        return corr;
    } else if(method == "D") {
        double estD = pAB - pA*pB;// Calculate D.
        return estD;
    } else if(method == "D'") {
        double estDp = pAB-pA*pB>0 ? (pAB-pA*pB)/Dmax : (pAB-pA*pB)/Dmin;// Calculate D-prime.
        return estDp;
    } else {
        std::cerr << "unknown method! input one of methods(r^2, D, D')" << endl;
        return 0;
    }
}



double methods(arma::Col<int> &g1, arma::Col<int> &g2, const string &method) {
    double res;
    if(method == "RMI") {
        res = RMI(g1, g2);
    } else if(method == "MI") {
        res = MI(g1, g2);
    } else {
        res = LD(g1, g2, method);
    }
    return res;
}

arma::Mat<double> methods_mat(arma::Mat<int> &geno012, const string &method ="RMI", int thread=1) {
    int n = geno012.n_cols;
    arma::Mat<double> m(n, n);
    m.fill(datum::nan);// Initialize an upper-triangular result with NaN values.
    int i, j;
    double res;
    arma::Col<int> g1, g2;
    omp_set_num_threads(thread);
    #pragma omp parallel for private(i, j, res, g1, g2)
    for (i=0; i< n-1; i++) {
        for (j= i+1; j < n; j++) {
            g1 = geno012.col(i);
            g2 = geno012.col(j);
            res = methods(g1, g2, method);
            m(i,j) = res;
        }
    }     
    return m; 
}

vector<vector<string>> calc_decay(arma::Mat<int> &geno012, vector<vector<string>> &snpinfo , const string &method, int maxdist=300) { 
    int dist, rows = snpinfo.size();
    double res;
    arma::Col<int> g1, g2;
    vector<vector<string>> decay_res;
    vector<string> vec(4);
    for (int i=0; i<rows-1; i++){
        for(int j=i+1; j<rows; j++) {
            dist = std::stoi(snpinfo[j][1]) - std::stoi(snpinfo[i][1]);
            if(snpinfo[j][0]==snpinfo[i][0] && dist <= maxdist*1000) {
                g1 = geno012.col(i);
                g2 = geno012.col(j);
                res = methods(g1, g2, method);
                vec[0] = snpinfo[i][2];//ID
                vec[1] = snpinfo[j][2];//ID
                vec[2] = to_string(dist);
                vec[3] = to_string(res);
                decay_res.push_back(vec);
            }
        }
    }
    return decay_res;
}

vector<vector<string>> circos(arma::Mat<double> &m, vector<vector<string>> &snpinfo , const string &method, double threshold=0.2) { 
    int rows = m.n_rows;
    vector<vector<string>> circos_res;
    vector<string> vec(7);
    for(int i=0; i<rows-1; i++) {
        for(int j=i+1; j<rows; j++) {
            if(snpinfo[j][0]!=snpinfo[i][0] && m(i,j) > threshold) {
                vec[0] = snpinfo[i][0];//CHROM
                vec[1] = snpinfo[i][1];//POS
                vec[2] = snpinfo[i][2];//ID
                vec[3] = snpinfo[j][0];//CHROM
                vec[4] = snpinfo[j][1];//POS
                vec[5] = snpinfo[j][2];//ID
                vec[6] = to_string(m(i,j));
                circos_res.push_back(vec);
            }
        }
    }
    return circos_res;
}

void by_chrom(arma::Mat<double> &m, vector<vector<string>> &chrinfo) {
    int len = chrinfo.size();
    int start, end;
    for (int chr=0; chr<len; chr++) {
        start = std::stoi(chrinfo[chr][1]);
        end = std::stoi(chrinfo[chr][2]);
        arma::Mat<double> sub_mat = m.submat(start,start,end,end);
        string outname= "chr"+ chrinfo[chr][1]+".heatmap";
        writearma(sub_mat, outname, '\t');
    }
}
