#include "gwld/core.hpp"
#include <iostream>
#include <random>
#include <set>

void require(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void near(double a,double b,double tolerance=1e-10) {
    require((std::isnan(a)&&std::isnan(b)) || std::abs(a-b)<=tolerance,"Numeric mismatch");
}
// Independent v1 formula, with variable-size tables and explicit entropy sums.
double reference(const std::vector<int>& a,const std::vector<int>& b,bool reduced) {
    std::set<int> ra,cb;
    for (auto x:a) if (x>=0) ra.insert(x);
    for (auto x:b) if (x>=0) cb.insert(x);
    std::vector<int> r(ra.begin(),ra.end()),s(cb.begin(),cb.end());
    std::vector<double> m(r.size()*s.size()),rows(r.size()),cols(s.size()); double n=0;
    for (std::size_t k=0;k<a.size();++k) if(a[k]>=0 && b[k]>=0) {
        const auto i=std::find(r.begin(),r.end(),a[k])-r.begin(),j=std::find(s.begin(),s.end(),b[k])-s.begin();
        ++m[i*s.size()+j]; ++rows[i]; ++cols[j]; ++n;
    }
    if (!n) return GWLD::nan();
    auto entropy=[n](const std::vector<double>& counts){ double h=0; for(auto v:counts) if(v) h-=v/n*std::log2(v/n); return h; };
    if(!reduced) return entropy(rows)+entropy(cols)-entropy(m);
    double I=std::lgamma(n+1),lx=0,ly=0,sx=0,sy=0;
    const double R=r.size(),S=s.size(),w=n/(n+0.5*R*S);
    for(auto v:m) I+=std::lgamma(v+1);
    for(auto v:rows) { I-=std::lgamma(v+1); double x=(1-w)/R+w*v/n; lx+=std::log(x); sx+=x*x; }
    for(auto v:cols) { I-=std::lgamma(v+1); double y=(1-w)/S+w*v/n; ly+=std::log(y); sy+=y*y; }
    const double nu=(S+1)/(S*sx)-1/S,mu=(R+1)/(R*sy)-1/R;
    const double omega=(R-1)*(S-1)*std::log(n+0.5*R*S)+0.5*(R+nu-2)*ly+0.5*(S+mu-2)*lx
        +0.5*(std::lgamma(mu*R)+std::lgamma(nu*S)-R*(std::lgamma(S)+std::lgamma(mu))-S*(std::lgamma(R)+std::lgamma(nu)));
    return (I-omega)/(n*std::log(2.0));
}
double evaluate(const std::vector<int>& a,const std::vector<int>& b,GWLD::Method method) {
    return GWLD::pair(a.data(),b.data(),a.size(),GWLD::summarize(a.data(),a.size()),GWLD::summarize(b.data(),b.size()),method);
}
int main() {
    try {
        std::mt19937 random(2026);
        for (int test=0;test<500;++test) {
            std::vector<int> a(100),b(100);
            for (int i=0;i<100;++i) { a[i]=static_cast<int>(random()%4)-1; b[i]=static_cast<int>(random()%4)-1; }
            near(evaluate(a,b,GWLD::Method::MI),reference(a,b,false));
            near(evaluate(a,b,GWLD::Method::RMI),reference(a,b,true));
        }
        // Cached integer/half-integer logs must preserve all 1/2/3-category cases,
        // including categories present only outside the complete-pair table.
        std::vector<double> factorials(1025),half_logs(2058);
        for (unsigned i=0;i<factorials.size();++i) factorials[i]=GWLD::log_gamma(i+1.0);
        for (unsigned i=1;i<half_logs.size();++i) half_logs[i]=std::log(0.5*i);
        for (unsigned r=1;r<=3;++r) for (unsigned s=1;s<=3;++s) for (unsigned trial=0;trial<100;++trial) {
            GWLD::Counts table; GWLD::Summary row,column;
            for (unsigned i=0;i<r;++i) row.present[i]=true;
            for (unsigned j=0;j<s;++j) column.present[j]=true;
            for (unsigned i=0;i<r;++i) for (unsigned j=0;j<s;++j) {
                const double value=random()%3 ? random()%100 : 0;
                table.cells[3*i+j]=value; table.rows[i]+=value; table.cols[j]+=value; table.n+=value;
            }
            near(GWLD::information(table,row,column,true),
                 GWLD::information(table,row,column,true,factorials.data(),nullptr,half_logs.data()),1e-12);
        }
        const std::vector<int> a={0,0,1,1,2,2},b={2,2,1,1,0,0},missing(6,-1),mono(6,0);
        near(evaluate(a,a,GWLD::Method::R2),1);
        near(evaluate(a,b,GWLD::Method::R2),1);
        near(evaluate(a,a,GWLD::Method::Dprime),1);
        near(evaluate(a,b,GWLD::Method::Dprime),1);
        near(evaluate(a,a,GWLD::Method::D),0.25);
        require(std::isnan(evaluate(a,missing,GWLD::Method::RMI)),"All-missing RMI must be missing");
        require(std::isnan(evaluate(a,mono,GWLD::Method::R2)),"Monomorphic r2 must be missing");
        near(evaluate(a,mono,GWLD::Method::D),0);
        // Heterozygote-heavy data can have two internal likelihood maxima.
        const int multimodal[3][3]={{0,2,0},{2,262,5},{2,5,2}};
        std::vector<int> mx,my;
        for (int i=0;i<3;++i) for (int j=0;j<3;++j)
            for (int k=0;k<multimodal[i][j];++k) {mx.push_back(i);my.push_back(j);}
        const double mle=evaluate(mx,my,GWLD::Method::D);
        near(mle,-0.23361366206657,1e-8);
        for (auto& x:mx) x=2-x;
        near(evaluate(mx,my,GWLD::Method::D),-mle,1e-10);
        // Check the haplotype MLE against an independent dense likelihood grid.
        // Use complete data here so the feasible marginals are unambiguous.
        for (int trial=0;trial<40;++trial) {
            std::vector<int> x(90),y(90);
            double table[3][3]={},sumx=0,sumy=0;
            for (int k=0;k<90;++k) {
                x[k]=random()%3; y[k]=random()%3; ++table[x[k]][y[k]];
                sumx+=x[k]; sumy+=y[k];
            }
            const double p=1-sumx/180,q=1-sumy/180;
            const double lo=std::max(0.0,p+q-1),hi=std::min(p,q);
            auto likelihood=[&](double h) {
                const double probs[3][3]={
                    {h*h,2*h*(p-h),(p-h)*(p-h)},
                    {2*h*(q-h),2*(h*(1-p-q+h)+(p-h)*(q-h)),2*(p-h)*(1-p-q+h)},
                    {(q-h)*(q-h),2*(q-h)*(1-p-q+h),(1-p-q+h)*(1-p-q+h)}};
                double result=0;
                for (int i=0;i<3;++i) for(int j=0;j<3;++j) if(table[i][j]) {
                    if (probs[i][j]<=0) return -std::numeric_limits<double>::infinity();
                    result+=table[i][j]*std::log(probs[i][j]);
                }
                return result;
            };
            double best=lo,best_score=-std::numeric_limits<double>::infinity();
            for (int k=0;k<=20000;++k) {
                const double h=lo+(hi-lo)*k/20000;
                const double score=likelihood(h);
                if(score>best_score) {best=h; best_score=score;}
            }
            const double D=evaluate(x,y,GWLD::Method::D);
            near(D,best-p*q,1e-4);
            near(evaluate(x,y,GWLD::Method::R2),D*D/(p*q*(1-p)*(1-q)),1e-12);
        }
        // Independent likelihood grids for sparse tables and near-ambiguous modes.
        for (int trial=0;trial<100;++trial) {
            int table[3][3]={};
            std::vector<int> x,y;
            for (int i=0;i<3;++i) for (int j=0;j<3;++j) {
                table[i][j]=i==1 && j==1 ? 20+random()%281 : random()%6;
                for (int k=0;k<table[i][j];++k) {x.push_back(i);y.push_back(j);}
            }
            const double p=GWLD::summarize(x.data(),x.size()).p,q=GWLD::summarize(y.data(),y.size()).p;
            const double lo=std::max(0.0,p+q-1),hi=std::min(p,q);
            const auto likelihood=[&](double h) {
                const double probability[3][3]={{h*h,2*h*(p-h),(p-h)*(p-h)},
                    {2*h*(q-h),2*(h*(1-p-q+h)+(p-h)*(q-h)),2*(p-h)*(1-p-q+h)},
                    {(q-h)*(q-h),2*(q-h)*(1-p-q+h),(1-p-q+h)*(1-p-q+h)}};
                double score=0;
                for (int i=0;i<3;++i) for (int j=0;j<3;++j) if (table[i][j]) {
                    if (probability[i][j]<=0) return -std::numeric_limits<double>::infinity();
                    score+=table[i][j]*std::log(probability[i][j]);
                }
                return score;
            };
            double grid=-std::numeric_limits<double>::infinity();
            for (int k=0;k<=20000;++k) grid=std::max(grid,likelihood(lo+(hi-lo)*k/20000));
            const auto D=evaluate(x,y,GWLD::Method::D);
            require(likelihood(D+p*q)>=grid-1e-8,"Missed a higher likelihood mode");
        }
        std::vector<int> na(6,std::numeric_limits<int>::min());
        require(std::isnan(evaluate(a,na,GWLD::Method::MI)),"R NA must be excluded");
        bool failed=false; try { int bad=3; GWLD::summarize(&bad,1); } catch(const std::invalid_argument&) {failed=true;}
        require(failed,"Invalid genotype accepted");
        std::vector<GWLD::Site> sites={{"1","a",2000},{"2","b",100},{"1","c",1000},{"1","d",5000000000LL},{"1","e",2000}};
        const GWLD::Window window(sites);
        require(window.partners(0,sites[0],1000)==std::vector<std::size_t>({2,4}),"Unsorted window or duplicate position failed");
        require(window.partners(0,sites[0],0)==std::vector<std::size_t>({4}),"Zero-distance window failed");
        require(GWLD::position("5000000000")==5000000000LL,"Large position parsing failed");
        require(GWLD::position("5e9")==5000000000LL,"R scientific position parsing failed");
        for (const auto* value : {"4611686018427387904", "9223372036854775808", "-1", "1.5", "nan"}) {
            bool rejected=false;
            try { GWLD::position(value); } catch (const std::exception&) { rejected=true; }
            require(rejected,"Invalid or unsafe position accepted");
        }
        require(GWLD::parse_method("r2")==GWLD::parse_method("r^2"),"r2 alias failed");
        std::cout<<"Core numerical and window tests passed\n";
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n'; return 1;}
}
