// Frozen v1 utility reference for benchmarks.
#include <iostream>
#include <fstream>
#include <sstream>
#include <armadillo>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

namespace legacy_reference {
using namespace std;
using namespace arma;

vector<char> is_any_of(string str) {
    vector<char> res;
    for (auto s : str)
        res.push_back(s);
    return res;
}

void split(vector<string>& result, string str, vector<char> delimiters) {
    auto start = 0;
    while (start < str.size())
    {
        // Split on any of the supplied delimiters.
        auto itRes = str.find(delimiters[0], start);
        for (int i = 1; i < delimiters.size(); ++i)
        {
            auto it = str.find(delimiters[i],start);
            if (it < itRes)
                itRes = it;
        }
        if (itRes == string::npos)
        {
            result.push_back(str.substr(start, str.size() - start));
            break;
        }
        result.push_back(str.substr(start, itRes - start));
        start = itRes;
        ++start;
    }
}

vector<vector<string>> readfile(const string &filename, char sep) {// Use the space character as the delimiter, not a string literal.

    ifstream inFile;
    inFile.open(filename, ios::in);
    if (!inFile.is_open())
	{
		std::cerr << filename <<" 读取文件失败" << endl;
        exit(1);
	}
    // Append each parsed row to the result.
    string str_line;
	vector<vector<string>> str_array;
    while (getline(inFile, str_line)) {
        stringstream ss(str_line);
        string str;
        vector<string> line_array;
        while (getline(ss, str, sep)) {// PED fields are separated by spaces.
            line_array.push_back(str);
        }
        str_array.push_back(line_array);
    }
    inFile.close();

    return str_array;
}

template <typename T> // The output container must be vector<vector<T>>.
void writefile(T &data, string &outfilename, char sep) {
    int rows = data.size();
    int cols = data.front().size();
    ofstream outFile;
    outFile.open(outfilename, ios::out);
    if (outFile.is_open()){
        for(int row=0;row<rows; row++) {
            for(int col=0;col<cols; col++) {
                outFile << data[row][col] << sep;
            }
            outFile <<"\n";// Terminate the row with a newline.
        }
    } else {
        std::cerr << outfilename <<" 创建文件失败" << endl;
		return;
    }
    outFile.close();
}
template void writefile(vector<vector<string>> &data, string &outfilename, char sep);
template void writefile(vector<vector<double>> &data, string &outfilename, char sep);


void writearma(arma::Mat<double> &data, string &outfilename, char sep) {
    int rows = data.n_rows;
    int cols = data.n_cols;
    ofstream outFile;
    outFile.open(outfilename, ios::out);
    if (outFile.is_open()){
        for(int row=0;row<rows; row++) {
            for(int col=0;col<cols; col++) {
                outFile << data(row, col) << sep;
            }
            outFile <<"\n";// Terminate the row with a newline.
        }
    } else {
        std::cerr << outfilename <<" 创建文件失败" << endl;
		return;
    }
    outFile.close();
}

bool cmp(pair<string,int> a, pair<string, int> b) {
    return a.second > b.second; // Sort in descending order.
}

vector<pair<string,int>> counts_allele(vector<string> &allele) {
    map<string, int> M;
    for (int i = 0; i<allele.size(); i++) {
        if (M.find(allele[i]) == M.end()) {
            M[allele[i]] = 1;
        }
        else {
            M[allele[i]]++;
        }
    }
    vector<pair<string,int>> sort_alleles(M.begin(), M.end());// Sort in descending order.
    sort(sort_alleles.begin(), sort_alleles.end(), cmp);// There are at most three allele tokens; "0" denotes a missing allele.
    return sort_alleles;
}

vector<double> allele_freq(arma::Col<int> &g1) {
    vector<double> stat(3, 0.0);
    int len = g1.n_elem;// Get the vector length.
    for(int i=0; i<len; i++) {
        if(g1[i]!= -1) { // -1 denotes a missing genotype.
            stat[g1[i]] += 1.0;
        }
    }
    double sum = stat[0] + stat[1] + stat[2];
    double A1 = (stat[0]+0.5*stat[1])/sum;//major allele
    double A2 = 1.0 - A1;//minor allele
    double NA = 1.0 - sum/len;// Calculate the missing genotype rate.
    stat[0] = A1;
    stat[1] = A2;
    stat[2] = NA;
    return stat;
}

arma::Mat<double> table(arma::Col<int> &g1, arma::Col<int> &g2, int row, int col) {
    arma::Mat<double> m(row, col, arma::fill::zeros);
    int g1_len = g1.n_elem; // Get the vector length.
    for(int i=0;i< g1_len; i++){
            if(g1[i]!= -1 && g2[i]!= -1){
                m(g1[i], g2[i]) +=1;
        }
    }
    return m;
}

vector<vector<double>> stats_freq(arma::Mat<int> &geno012){
    int cols = geno012.n_cols;
    vector<vector<double>> alleles_freq;
    arma::Col<int> g1;
    for(int l=0;l<cols;l++) {
        g1 = geno012.col(l);
        alleles_freq.push_back(allele_freq(g1));
    }
    return alleles_freq;
}

vector<vector<string>> getsnpinfo(vector<vector<string>> &data, vector<int> &col) {
    int rows = data.size();
    int len = col.size();
    vector<vector<string>> snpinfo;
    vector<string> vec(3);//CHROM, POS, ID
    for (int row=0; row<rows; row++) {
        for (int l=0; l<len; l++) {
            vec[l] = data[row][col[l]];
        }
        snpinfo.push_back(vec);
    }
    return(snpinfo);
}

vector<vector<string>> getchrinfo(vector<vector<string>> &snpinfo) {
    vector<vector<string>> chrinfo;  // Store chromosome summaries and SNP counts.
    vector<string> chr(4);           // Chromosome name, start index, end index, and SNP count.
    int index=0, rows = snpinfo.size();
    for(int i=0; i<rows; i++) {
        if(snpinfo[index][0]==snpinfo[i][0] && i!=(rows-1)) {
            continue;
        } else {
            if(i==(rows-1)) i++;
            chr[0] = snpinfo[index][0];
            chr[1] = to_string(index);
            chr[2] = to_string(i - 1);
            chr[3] = to_string(i - index);
            chrinfo.push_back(chr);
            index = i;
        }
    }
    return chrinfo;
}

} // namespace legacy_reference
