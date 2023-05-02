// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <map>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;
using namespace Rcpp;
using namespace arma;

bool cmp(pair<std::string, double> a, pair<string, double> b) {
     //从大到小排列
    return a.second > b.second;
}

std::vector<pair<std::string, double>> split_str(CharacterVector &v, string sep) {
    map<std::string, double> m;
    std::string str, tmp;
    for(int i=0; i<v.size(); i++) {
        if(v[i] == NA_STRING) continue;
        str = v[i];
        if(sep.size()==0) {
            for(size_t j=0; j<str.size(); j++) {
                std::string temp;
                temp.push_back(str[j]);
                m[temp] +=1;
            }
        } else {
            size_t start = 0;
            while(start < str.size()) {
                size_t index = str.find(sep[0], start);
                for (size_t i = 1; i < sep.size(); ++i) {
                    size_t it = str.find(sep[i],start);
                    if (it < index) index = it;
                }
                if(index == std::string::npos) {
                    tmp = str.substr(start, str.size()-start);
                    m[tmp] += 1;
                    break;
                }
                tmp = str.substr(start, index-start);
                m[tmp] += 1;
                start = index;
                start++;
            }
        }
    }
    std::vector<pair<std::string, double>> sorted(m.begin(), m.end());
    sort(sorted.begin(), sorted.end(), cmp);
    return sorted;
}

//' stats allele frequency
//'  
//' stats allele frequency of character genotypes vector 
//'
//' @param vector a character vector of genotypes
//' @param sep the separator for character genotypes(i.e., "" ,"/", "|")
//' 
//' @return list
//' @export
//' @name stats_freq
//' @rdname stats_freq
//' @examples
//' 
//' g0 <- c('T/A', NA, 'T|T', NA, 'T|A', NA, 'T|T', 'T/A','T/T', 'T/T', 
//'         'T/A', 'A|A', 'T/T', 'T|A', 'T/A', 'T|T', NA, 'T|A', 'T/A', NA)
//' 
//' stats_freq(g0, "/|")

// [[Rcpp::export]]
List stats_freq(CharacterVector &vector, std::string sep) {
    std::vector<pair<std::string, double>> sorted = split_str(vector, sep);
    Rcpp::CharacterVector allele;
    Rcpp::NumericVector counts;// 换成Rcpp::NumericVector
    for(size_t j=0; j<sorted.size(); j++){
        allele.push_back(sorted[j].first);
        counts.push_back(sorted[j].second);
    }
    int NA_counts = vector.size() - 0.5*Rcpp::sum(counts);
    double NA_rate = (double)NA_counts/vector.size();
    Rcpp::NumericVector freq= counts/Rcpp::sum(counts);
    Rcpp::List L = Rcpp::List::create(_["allele"] = allele,
                                      _["allele counts"] = counts,
                                      _["allele freq"] = freq,
                                      _["missing counts"] = NA_counts,
                                      _["missing rate"] = NA_rate);
    return L;
}


string major_allele(CharacterVector &v, std::string sep) {
    std::vector<pair<std::string, double>> sorted = split_str(v, sep);
    string A1 = sorted[0].first;
    return A1;
}


// Character gentype vector recode 0, 1, 2, NA. NA for missing value
arma::Col<int> code_Vec(CharacterVector &vector, std::string sep) {
    int n=vector.size();
    arma::Col<int> geno012_col(n);
    string A1 = major_allele(vector, sep);
    for(int i=0; i<n; i++) {
        int counts=0;
        size_t start = 0;
        if(vector[i] == NA_STRING) {
            geno012_col(i) = NA_INTEGER;//表示缺失
        } else {
            std::string str = Rcpp::as<std::string>(vector[i]);// R语言字符转string
            while (start < str.size()) {
                size_t pos = str.find(A1, start); //有问题
                if(pos == std::string::npos) {
                    break;
                }
                counts +=1;
                start = pos;
                start++;
            }
            if(counts==0) geno012_col(i) = 2;
            if(counts==1) geno012_col(i) = 1;
            if(counts==2) geno012_col(i) = 0;
        }
    }
    return geno012_col;
}

//' Convert character genotype to integer genotype
//' 
//' convert genotypes separated by ("", "/", "|", ...) to integer genotypes(0, 1, 2, NA, NA for missing value)
//'
//' @param vector a character vector of genotypes
//' @param sep the separator for character genotypes(i.e., "" ,"/", "|")
//' 
//' @name recode
//' @rdname recode
//' @export
//' @examples
//' 
//' g0 <- c('T/A', NA, 'T|T', NA, 'T|A', NA, 'T|T', 'T/A','T/T', 'T/T', 
//'         'T/A', 'A|A', 'T/T', 'T|A', 'T/A', 'T|T', NA, 'T|A', 'T/A', NA)
//'         
//' code012_g1 <- recode(g0, "/|")
//' 
//' g1 <- c('T/A', NA, 'T/T', NA, 'T/A', NA, 'T/T', 'T/A','T/T', 'T/T', 
//'         'T/A', 'A/A', 'T/T', 'T/A', 'T/A', 'T/T', NA, 'T/A', 'T/A', NA)  
//' 
//' code012_g1 <- recode(g1, "/")
//' 
//' g2 <- c('01', '01', '00', '01', '00', '01', '01', '01', '01', '00', 
//'         '01', '11', '01', '11', '01', '11','01', '01', '01', '11') 
//'         
//' code012_g2 <- recode(g2, "")

// [[Rcpp::export]]
Rcpp::IntegerVector recode(CharacterVector &vector, std::string sep) {
  arma::Col<int> geno012_col = code_Vec(vector, sep);//列向量转成R语言vector
  return Rcpp::IntegerVector(geno012_col.begin(), geno012_col.end());
}

// [[Rcpp::export]]
arma::Mat<int> code_Mat(CharacterMatrix &Matrix, std::string sep) {
    int rows=Matrix.nrow(), cols=Matrix.ncol();
    arma::Mat<int> geno012(rows, cols);
    for(int i=0; i<cols; i++) {
        CharacterVector v = Matrix.column(i);
        arma::Col<int> geno012_col=code_Vec(v, sep);
        geno012.col(i) = geno012_col;
    }
    return geno012;
}
