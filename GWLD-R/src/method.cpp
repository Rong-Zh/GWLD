// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <float.h>
#include <omp.h>
#include <map>
#include <vector>
#include <string>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double RMIC(arma::Col<int> & g1, arma::Col<int> & g2) {
    if (g1.n_elem != g2.n_elem) {
        throw Rcpp::exception("g1 and g2 have different length");
    }
    //unique去重并按升序排
    arma::Col<int> set_g1 = arma::unique(g1); 
    arma::Col<int> set_g2 = arma::unique(g2); 
    //去掉set中NA值
    arma::Col<int> row_names = set_g1.elem(find(set_g1 != NA_INTEGER));
    arma::Col<int> col_names = set_g2.elem(find(set_g2 != NA_INTEGER));
    int len = g1.n_elem;
    //防止进度丢失
    double R = row_names.n_elem, S = col_names.n_elem; 
    std::unordered_map<int, int> row_index;
    for (int i=0; i<R; i++) {
        row_index[row_names[i]] = i;
    }
    std::unordered_map<int, int> col_index;
    for (int j=0; j<S; j++) {
        col_index[col_names[j]] = j;
    }
    //使用double防止进度丢失，用0.0来填充
    arma::Mat<double> m(R, S, fill::zeros);
    for(int k=0; k<len; k++){
        if(g1[k]!= NA_INTEGER && g2[k]!= NA_INTEGER){
            m(row_index[g1[k]], col_index[g2[k]]) +=1;
        }
    }
    //所有元素和
    double n = accu(m);
    arma::Col<double> a = sum(m,1); 
    arma::Row<double> b = sum(m,0); //1是行, 0是列
    double I = (lgamma(n+1) + accu(lgamma(m+1)) - accu(lgamma(a+1)) - accu(lgamma(b+1))) / (n*log(2));
    double w = n / (n + 0.5*R*S);
    arma::Col<double> x = (1-w)/R + w*a/n;
    arma::Row<double> y = (1-w)/S + w*b/n;
    double nu = (S+1)/(S*accu(square(x))) - 1/S;
    double mu = (R+1)/(R*accu(square(y))) - 1/R;
    double logOmega = (R-1)*(S-1)*log(n+0.5*R*S)+
    0.5*(R+nu-2)*accu(log(y))+
    0.5*(S+mu-2)*accu(log(x))+
    0.5*(lgamma(mu*R)+lgamma(nu*S)-R*(lgamma(S)+lgamma(mu))-S*(lgamma(R)+lgamma(nu)));
    double rmi = I-logOmega/(n*log(2));
    return rmi;
}

// [[Rcpp::export]]
arma::Mat<double> RMIC_Mat(arma::Mat<int> & geno012, int cores=1) {
    //设置核心数
    omp_set_num_threads(cores);
    int n = geno012.n_cols;
    arma::Mat<double> m(n, n);
    m.fill(NA_REAL);//用NA填充
    int i, j;
    double res;
    arma::Col<int> g1, g2;
    #pragma omp parallel for private(i, j, res, g1, g2) //线程私有
    for (i=0; i< n-1; i++) {
        for (j= i+1; j < n; j++) {
            g1 = geno012.col(i);
            g2 = geno012.col(j);
            res = RMIC(g1, g2);
            m(i,j) = res;
        }
    }     
    return m; 
}

// [[Rcpp::export]]
double MIC(arma::Col<int> & g1, arma::Col<int> & g2) {
    if (g1.n_elem != g2.n_elem) {
        throw Rcpp::exception("g1 and g2 have different length");
    }
    //unique去重并按升序排
    arma::Col<int> set_g1 = arma::unique(g1); 
    arma::Col<int> set_g2 = arma::unique(g2); 
    //去掉set中NA值
    arma::Col<int> row_names = set_g1.elem(find(set_g1 != NA_INTEGER));
    arma::Col<int> col_names = set_g2.elem(find(set_g2 != NA_INTEGER));
    int len = g1.n_elem;
    //防止进度丢失
    double R = row_names.n_elem, S = col_names.n_elem; 
    std::unordered_map<int, int> row_index;
    for (int i=0; i<R; i++) {
        row_index[row_names[i]] = i;
    }
    std::unordered_map<int, int> col_index;
    for (int j=0; j<S; j++) {
        col_index[col_names[j]] = j;
    }
    //使用double防止进度丢失，用0.0来填充
    arma::Mat<double> m(R, S, fill::zeros);
    for(int k=0; k<len; k++){
        if(g1[k]!= NA_INTEGER && g2[k]!= NA_INTEGER){
            m(row_index[g1[k]], col_index[g2[k]]) +=1;
        }
    }
    //下面这种方式频繁的调用find函数
    // arma::Mat<double> m(R, S, fill::zeros);//用0来填充
    // for(int i=0; i<len; i++){
    //     if(g1[i]!= NA_INTEGER && g2[i]!= NA_INTEGER){
    //         m(find(r == g1[i]), find(s == g2[i])) +=1;
    //     }
    // }
    double n = accu(m);//所有元素和 
    arma::Mat<double> freq = m/n;//频率table
    arma::Col<double> a = sum(freq,1); //行和
    arma::Row<double> b = sum(freq,0); //列和
    //定义熵的公式,定义0*log2(0)=0,单位为比特.
    arma::Col<double> rowh = a.transform([] (double val) {return(val==0.0 ? 0.0 : -val*log2(val));});
    arma::Row<double> colh = b.transform([] (double val) {return(val==0.0 ? 0.0 : -val*log2(val));});
    arma::Mat<double> tabh = freq.transform([] (double val) {return(val==0.0 ? 0.0 : -val*log2(val));});
    //计算位点g1的信息熵
    double H_fg = accu(rowh);
    //计算位点g2的信息熵
    double H_sg = accu(colh);
    //计算联合熵
    double H_fg_sg = accu(tabh);
    //计算互信息
    double I = H_fg + H_sg - H_fg_sg;
    return I;
}

// [[Rcpp::export]]
arma::Mat<double> MIC_Mat(arma::Mat<int> & geno012, int cores=1) {
    //设置核心数
    omp_set_num_threads(cores);
    int n = geno012.n_cols;
    arma::Mat<double> m(n, n);
    m.fill(NA_REAL);//用NA填充
    int i, j;
    double res;
    arma::Col<int> g1, g2;
    #pragma omp parallel for private(i, j, res, g1, g2)//线程私有
    for (i=0; i< n-1; i++) {
        for (j= i+1; j < n; j++) {
            g1 = geno012.col(i);
            g2 = geno012.col(j);
            res = MIC(g1, g2);
            m(i,j) = res;
        }
    }     
    return m; 
}

double major_allele_freq(arma::Col<int> &g1) {
    std::vector<double> stat(3, 0.0);
    int len = g1.n_elem;//向量的长度
    for(int i=0; i<len; i++) {
        if(g1[i]!= NA_INTEGER) { // NA_INTEGER 表示缺失值
            stat[g1[i]] += 1.0;
        }
    }
    double sum = stat[0] + stat[1] + stat[2];
    double major_freq = (stat[0]+0.5*stat[1])/sum;//major allele
    return major_freq;
}

arma::Mat<double> table(arma::Col<int> &g1, arma::Col<int> &g2, int row, int col) {
    arma::Mat<double> m(row, col, arma::fill::zeros);
    int g1_len = g1.n_elem; //向量的长度
    for(int i=0;i< g1_len; i++){
            if(g1[i]!= NA_INTEGER && g2[i]!= NA_INTEGER){
                m(g1[i], g2[i]) +=1;
        }
    }
    return m;
}

double nlogLik(double &pAB, double &pA , double &pB, arma::Mat<double> &m) {
    //负对数似然函数
    double y =-((2*m(0,0)+m(0,1)+m(1,0))*log(pAB) +
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
	/*  eps is approximately the square root of the relative machine precision. */
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
    /*  main loop starts here ----------------------------------- */
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
	/*  update  a, b, v, w, and x */
	if (fu <= fx) {
	    if (u < x) b = x; else a = x;
	    v = w; w = x;   x = u;
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

// [[Rcpp::export]]
double LDC(arma::Col<int> &g1, arma::Col<int> &g2, std::string method="r^2") {
    if (g1.n_elem != g2.n_elem) {
        throw Rcpp::exception("g1 and g2 have different length");
    }
    //统计第一位点频率
    double pA, pa, pB, pb;
    pA =  major_allele_freq(g1);
    pa = 1.0 -pA;
    //统计第二位点频率
    pB =  major_allele_freq(g2);
    pb = 1.0 -pB;
    //两个位点之间形成3*3的table, 用0来填充
    arma::Mat<double> m = table(g1, g2, 3, 3);
    double Dmin = std::max(-pA*pB, -pa*pb);
    double pmin = pA*pB + Dmin;
    double Dmax = std::min(pA*pb, pB*pa);
    double pmax = pA*pB + Dmax;
    double eps = DBL_EPSILON;
    double tol = pow(DBL_EPSILON, 0.25);//设置精度
    double pAB = fmin(pmin+eps, pmax+eps, tol, pA, pB, m);
    if(method == "r^2") {
        double corr = pow(pAB - pA*pB, 2)/(pA * pB * pa * pb);//计算R^2
        return corr;
    } else if(method == "D") {
        double estD = pAB - pA*pB;//计算D值
        return estD;
    } else if(method == "D'") {
        double estDp = pAB-pA*pB>0 ? (pAB-pA*pB)/Dmax : (pAB-pA*pB)/Dmin;//计算D'值
        return estDp;
    } else {
        throw Rcpp::exception("Unknown method! Input one of methods(D, D', r^2) ");
    }
}

// [[Rcpp::export]]
arma::Mat<double> LDC_Mat(arma::Mat<int> & geno012, std::string method, int cores=1) {
    if(method=="D" || method=="D'" || method=="r^2") {
        //设置核心数
        omp_set_num_threads(cores);
        int n = geno012.n_cols;
        arma::Mat<double> m(n, n);
        m.fill(NA_REAL);//用NA填充
        int i, j;
        double res;
        arma::Col<int> g1, g2;
        #pragma omp parallel for private(i, j, res, g1, g2)
        for (i=0; i< n-1; i++) {
            for (j= i+1; j < n; j++) {
                g1 = geno012.col(i);
                g2 = geno012.col(j);
                res = LDC(g1, g2, method);
                m(i,j) = res;
            }
        }
        return m;
    } else {
      throw Rcpp::exception("Unknown method! Input one of methods(D, D', r^2) ");
    }
}

double methods(arma::Col<int> &g1, arma::Col<int> &g2, std::string method) {
    double res;
    if(method == "RMI") {
        res = RMIC(g1, g2);
    } else if(method == "MI") {
        res = MIC(g1, g2);
    } else {
        res = LDC(g1, g2, method);
    }
    return res;
}

// [[Rcpp::export]]
arma::Mat<double> gwld(arma::Mat<int> &geno012, std::string method="RMI", int thread=1) {
    int n = geno012.n_cols;
    arma::Mat<double> m(n, n);
    m.fill(NA_REAL);//用NA填充,上三角矩阵.
    int i, j;
    double res;
    arma::Col<int> g1, g2;
    omp_set_num_threads(thread);
    #pragma omp parallel for private(i, j, res, g1, g2)
    for (i=0; i<n-1; i++) {
        for (j= i+1; j<n; j++) {
            g1 = geno012.col(i);
            g2 = geno012.col(j);
            res = methods(g1, g2, method);
            m(i,j) = res;
        }
    }     
    return m; 
}

CharacterMatrix cMat2SEXP(std::vector<std::vector<std::string>> &cMat, CharacterVector &setcolnames) {
    int nrows = cMat.size(), ncols=cMat.front().size();
    if(nrows==0) ncols = setcolnames.size();
    CharacterMatrix rMat(nrows, ncols);
    for(int row=0; row<nrows; row++ ){
        for(int col=0; col<ncols; col++) {
            rMat(row, col) = cMat[row][col];
        }
    }
    colnames(rMat) = setcolnames;//设置列名
    return rMat;
}

// [[Rcpp::export]]
CharacterMatrix decay(arma::Mat<double> &m, CharacterMatrix &snpinfo, int maxdist=300) { 
    int dist, rows = snpinfo.nrow(),  m_cols= m.n_cols;
    if(rows != m_cols) {
        throw Rcpp::exception("The columns of m and rows of snpinfo are not equal");
    }
    std::string pos_1, pos_2, snp_1, snp_2;
    std::vector<std::vector<std::string>> decay_res;
    std::vector<std::string> row_v(4);
    for (int i=0; i<rows-1; i++){
        for(int j=i+1; j<rows; j++) {
            pos_1 = Rcpp::as<std::string>(snpinfo(j, 1));//sexp 转换std::string
			pos_2 = Rcpp::as<std::string>(snpinfo(i, 1));
			dist = std::stoi(pos_1) - std::stoi(pos_2);
            if(snpinfo(j, 0)==snpinfo(i, 0) && dist <= maxdist*1000) {
                snp_1 = Rcpp::as<std::string>(snpinfo(i,2));
				snp_2 = Rcpp::as<std::string>(snpinfo(j,2));
                row_v = {snp_1, snp_2, std::to_string(dist), std::to_string(m(i,j))};
                decay_res.push_back(row_v);
            }
        }
    }
    //std::vector<std::vector<std::string>> 转换为CharacterMatrix
    CharacterVector names = {"ID_1", "ID_2", "Dist", "Value"};
    Rcpp::CharacterMatrix final_res = cMat2SEXP(decay_res, names);
    return final_res;
}

// [[Rcpp::export]]
CharacterMatrix calc_decay(arma::Mat<int> &geno012, CharacterMatrix &snpinfo, std::string method="RMI", int thread =1, int maxdist=300) {
    double res;
    int i, j, dist, rows = snpinfo.nrow();
    arma::Col<int> g1, g2;
	std::string pos_1, pos_2, snp_1, snp_2;
    std::vector<std::vector<std::string>> decay_res;
    std::vector<std::string> row_v(4);
    omp_set_num_threads(thread);
    #pragma omp parallel for private(i, j, pos_1, pos_2 ,snp_1, snp_2, dist, g1, g2, res) 
    for (i=0; i<rows-1; i++){
        for(j=i+1; j<rows; j++) {
			pos_1 = Rcpp::as<std::string>(snpinfo(j, 1));
			pos_2 = Rcpp::as<std::string>(snpinfo(i, 1));
			dist = std::stoi(pos_1) - std::stoi(pos_2);
            if(snpinfo(j, 0)==snpinfo(i, 0) && dist <= maxdist*1000) {
				g1 = geno012.col(i);
				g2 = geno012.col(j);
				res = methods(g1, g2, method);
				snp_1 = Rcpp::as<std::string>(snpinfo(i,2));
				snp_2 = Rcpp::as<std::string>(snpinfo(j,2));
				#pragma omp critical
				{ 
					row_v = {snp_1,  snp_2, std::to_string(dist), std::to_string(res)};
					decay_res.push_back(row_v);
				}
            }
        }
    }
    //std::vector<std::vector<std::string>> 转换为CharacterMatrix
    CharacterVector names = {"ID_1", "ID_2", "Dist", "Value"};
    Rcpp::CharacterMatrix final_res = cMat2SEXP(decay_res, names);
    return final_res;
}

// [[Rcpp::export]]
CharacterMatrix calc_circos(arma::Mat<int> &geno012, CharacterMatrix &snpinfo, std::string method="RMI", int thread=1, double threshold=0.2) {
    double res;
    int i, j, rows = snpinfo.nrow();
    arma::Col<int> g1, g2;
	std::string chr_1, chr_2, pos_1, pos_2, snp_1, snp_2;
    std::vector<std::vector<std::string>> circos_res;
    std::vector<std::string> row_v(4);
    omp_set_num_threads(thread);
    #pragma omp parallel for private(i, j, chr_1, chr_2, pos_1, pos_2 ,snp_1, snp_2, g1, g2, res) 
    for (i=0; i<rows-1; i++){
        for(j=i+1; j<rows; j++) {
            if(snpinfo(j, 0)!=snpinfo(i, 0)) {
                g1 = geno012.col(i);
			    g2 = geno012.col(j);
			    res = methods(g1, g2, method);
                if (res > threshold) {
                    chr_1 = Rcpp::as<std::string>(snpinfo(i, 0));
                    chr_2 = Rcpp::as<std::string>(snpinfo(j, 0));
                    pos_1 = Rcpp::as<std::string>(snpinfo(i, 1));
			        pos_2 = Rcpp::as<std::string>(snpinfo(j, 1));
                    snp_1 = Rcpp::as<std::string>(snpinfo(i, 2));
                    snp_2 = Rcpp::as<std::string>(snpinfo(j, 2));
                    #pragma omp critical
                    { 
                        row_v = {chr_1, pos_1, snp_1, chr_2, pos_2, snp_2, std::to_string(res)};
                        circos_res.push_back(row_v);
                    }
                }
            }
        }
    }
    //std::vector<std::vector<std::string>> 转换为CharacterMatrix
    CharacterVector names = {"CHR_1", "POS_1", "ID_1", "CHR_2", "POS_2", "ID_2", "Value"};
    Rcpp::CharacterMatrix final_res = cMat2SEXP(circos_res, names);
    return final_res;
}

// [[Rcpp::export]]
CharacterMatrix circos(arma::Mat<double> &m, CharacterMatrix &snpinfo, double threshold=0.2) {
    int rows = snpinfo.nrow(), m_cols= m.n_cols;
    if(rows != m_cols) {
        throw Rcpp::exception("The columns of m and rows of snpinfo are not equal");
    }
    std::string chr_1, chr_2, pos_1, pos_2, snp_1, snp_2;
    std::vector<std::vector<std::string>> circos_res;
    std::vector<std::string> row_v(7);
    for(int i=0; i<rows-1; i++) {
        for(int j=i+1; j<rows; j++) {
            if(snpinfo(j, 0)!=snpinfo(i, 0) && m(i,j) > threshold) {
                chr_1 = Rcpp::as<std::string>(snpinfo(i, 0));
                chr_2 = Rcpp::as<std::string>(snpinfo(j, 0));
                pos_1 = Rcpp::as<std::string>(snpinfo(i, 1));
                pos_2 = Rcpp::as<std::string>(snpinfo(j, 1));
                snp_1 = Rcpp::as<std::string>(snpinfo(i, 2));
                snp_2 = Rcpp::as<std::string>(snpinfo(j, 2));
                row_v = {chr_1, pos_1, snp_1, chr_2, pos_2, snp_2, std::to_string(m(i,j))};
                circos_res.push_back(row_v);
            }
        }
    }
    //std::vector<std::vector<std::string>> 转换为CharacterMatrix
    CharacterVector names = {"CHR_1", "POS_1", "ID_1", "CHR_2", "POS_2", "ID_2", "Value"};
    Rcpp::CharacterMatrix final_res = cMat2SEXP(circos_res, names);//转换并设置列名
    return final_res;
}