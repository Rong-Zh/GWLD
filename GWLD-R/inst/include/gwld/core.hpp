#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// Shared, allocation-free pair kernels for the CLI and the R package.
namespace GWLD {
enum class Method { RMI, MI, R2, D, Dprime };
inline Method parse_method(const std::string& name) {
    if (name == "RMI") return Method::RMI;
    if (name == "MI") return Method::MI;
    if (name == "r2" || name == "r^2") return Method::R2;
    if (name == "D") return Method::D;
    if (name == "D'") return Method::Dprime;
    throw std::invalid_argument("Unknown method: " + name);
}
inline double nan() { return std::numeric_limits<double>::quiet_NaN(); }
inline double log_gamma(double value) {
    #if defined(__GLIBC__) || defined(__APPLE__) || defined(__FreeBSD__)
    // Avoid lgamma's shared signgam state in parallel RMI evaluations.
    int sign;
    return ::lgamma_r(value,&sign);
    #else
    return std::lgamma(value);
    #endif
}
inline double log_gamma_ratio(double value,unsigned multiplier) {
    if (multiplier==1) return 0;
    return log_gamma(value*multiplier)-multiplier*log_gamma(value);
}
struct Summary {
    bool present[3] = {false, false, false};
    double p = nan();
    double counts[3] = {}, total=0;
};
inline Summary summarize(const int* g, std::size_t n) {
    Summary s;
    std::uint64_t zero=0,one=0,two=0;
    unsigned invalid=0;
    // Independent integer reductions vectorize without floating-point reassociation
    // and avoid indexed counter updates. Missing codes contribute to no category.
    for (std::size_t k=0;k<n;++k) {
        const int code=g[k];
        zero+=code==0; one+=code==1; two+=code==2;
        invalid|=((static_cast<unsigned>(code)+1u)>3u) & (code!=std::numeric_limits<int>::min());
    }
    if (invalid) throw std::invalid_argument("Genotypes must be 0, 1, 2 or missing (-1/NA)");
    s.counts[0]=static_cast<double>(zero); s.counts[1]=static_cast<double>(one); s.counts[2]=static_cast<double>(two);
    s.present[0]=zero>0; s.present[1]=one>0; s.present[2]=two>0;
    s.total=static_cast<double>(zero+one+two);
    if (s.total) s.p=(s.counts[0]+0.5*s.counts[1])/s.total;
    return s;
}
struct Counts {
    double cells[9] = {};
    double rows[3] = {}, cols[3] = {}, n = 0;
    double operator()(int i, int j) const { return cells[3*i+j]; }
};
inline Counts count(const int* a, const int* b, std::size_t n) {
    Counts c;
    for (std::size_t k = 0; k < n; ++k) {
        if (a[k] >= 0 && b[k] >= 0) {
            ++c.cells[3*a[k]+b[k]];
        }
    }
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            c.rows[i] += c(i,j); c.cols[j] += c(i,j); c.n += c(i,j);
        }
    return c;
}
inline double information(const Counts& c, const Summary& a, const Summary& b, bool reduced,
                          const double* log_factorials=nullptr,const double* count_logs=nullptr,
                          const double* half_logs=nullptr) {
    if (!c.n) return nan();
    if (!reduced && count_logs) {
        double value=count_logs[static_cast<std::size_t>(c.n)];
        for (double cell:c.cells) value+=count_logs[static_cast<std::size_t>(cell)];
        for (int i=0;i<3;++i) value-=count_logs[static_cast<std::size_t>(c.rows[i])]+count_logs[static_cast<std::size_t>(c.cols[i])];
        return std::max(0.0,value/c.n);
    }
    if (!reduced) {
        double value = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (c(i,j)) value += c(i,j)/c.n * std::log2(c(i,j)*c.n/(c.rows[i]*c.cols[j]));
        return value;
    }
    // Keep the original RMI convention: category counts use each whole locus,
    // while the contingency table uses complete pairs.
    double R = 0, S = 0;
    for (int i = 0; i < 3; ++i) { R += a.present[i]; S += b.present[i]; }
    const auto factorial=[log_factorials](double count) {
        return log_factorials ? log_factorials[static_cast<std::size_t>(count)] : log_gamma(count+1);
    };
    const auto category_gamma=[log_factorials](double categories) {
        return log_factorials ? log_factorials[static_cast<std::size_t>(categories)-1] : log_gamma(categories);
    };
    double I = factorial(c.n), sx = 0, sy = 0, lx = 0, ly = 0;
    double squared_rows=0,squared_cols=0;
    for (double v : c.cells) I += factorial(v);
    const double denominator=c.n+0.5*R*S;
    const double log_denominator=half_logs ? half_logs[static_cast<std::size_t>(2*c.n+R*S)] : std::log(denominator);
    const double w = c.n/denominator;
    for (int i = 0; i < 3; ++i) {
        if (a.present[i]) {
            I -= factorial(c.rows[i]);
            if (half_logs) {
                const double scaled=2*c.rows[i]+S;
                squared_rows+=scaled*scaled;
                lx+=half_logs[static_cast<std::size_t>(scaled)]-log_denominator;
            } else {
                const double x=(1-w)/R+w*c.rows[i]/c.n;
                sx+=x*x; lx+=std::log(x);
            }
        }
        if (b.present[i]) {
            I -= factorial(c.cols[i]);
            if (half_logs) {
                const double scaled=2*c.cols[i]+R;
                squared_cols+=scaled*scaled;
                ly+=half_logs[static_cast<std::size_t>(scaled)]-log_denominator;
            } else {
                const double y=(1-w)/S+w*c.cols[i]/c.n;
                sy+=y*y; ly+=std::log(y);
            }
        }
    }
    const double denominator_squared=4*denominator*denominator;
    const double nu = half_logs ? ((S+1)*denominator_squared/squared_rows-1)/S : (S+1)/(S*sx)-1/S;
    const double mu = half_logs ? ((R+1)*denominator_squared/squared_cols-1)/R : (R+1)/(R*sy)-1/R;
    const double gamma_term=half_logs ?
        log_gamma_ratio(mu,static_cast<unsigned>(R))+log_gamma_ratio(nu,static_cast<unsigned>(S))
            -R*category_gamma(S)-S*category_gamma(R) :
        log_gamma(mu*R)+log_gamma(nu*S)-R*(category_gamma(S)+log_gamma(mu))-S*(category_gamma(R)+log_gamma(nu));
    const double omega = (R-1)*(S-1)*log_denominator
        +0.5*(R+nu-2)*ly+0.5*(S+mu-2)*lx
        +0.5*gamma_term;
    return (I-omega)/(c.n*std::log(2.0));
}
inline double log_term(double coefficient, double probability) {
    if (coefficient == 0) return 0; // avoid 0 * log(0) at valid boundaries
    return probability > 0 ? coefficient*std::log(probability) : -std::numeric_limits<double>::infinity();
}
inline double nll(double x, double a, double b, const Counts& m) {
    return -(log_term(2*m(0,0)+m(0,1)+m(1,0), x)
        +log_term(2*m(0,2)+m(0,1)+m(1,2), a-x)
        +log_term(2*m(2,0)+m(1,0)+m(2,1), b-x)
        +log_term(2*m(2,2)+m(2,1)+m(1,2), 1-a-b+x)
        +log_term(m(1,1), x*(1-a-b+x)+(a-x)*(b-x)));
}
namespace likelihood_detail {
// The derivative of the fixed-marginal log-likelihood has a numerator of
// degree at most five. Work on u in [0,1] to condition narrow feasible ranges.
struct Polynomial {
    std::array<long double,7> c{};
    int degree=0;
    long double value(long double x) const {
        long double y=c[degree];
        for (int i=degree-1;i>=0;--i) y=y*x+c[i];
        return y;
    }
    Polynomial derivative() const {
        Polynomial result; result.degree=std::max(0,degree-1);
        for (int i=1;i<=degree;++i) result.c[i-1]=i*c[i];
        return result;
    }
    void normalize() {
        long double scale=0;
        for (int i=0;i<=degree;++i) scale=std::max(scale,std::abs(c[i]));
        if (!scale) { degree=0; return; }
        for (int i=0;i<=degree;++i) c[i]/=scale;
        while (degree>0 && std::abs(c[degree])<32*std::numeric_limits<long double>::epsilon()) --degree;
    }
};
inline Polynomial multiply(const Polynomial& a,const Polynomial& b) {
    Polynomial result; result.degree=a.degree+b.degree;
    for (int i=0;i<=a.degree;++i) for (int j=0;j<=b.degree;++j) result.c[i+j]+=a.c[i]*b.c[j];
    return result;
}
struct Roots {
    std::array<long double,5> x{};
    int size=0;
    void add(long double value) {
        if (value<0 || value>1) return;
        const auto tolerance=32*std::numeric_limits<long double>::epsilon();
        for (int i=0;i<size;++i) if (std::abs(value-x[i])<=tolerance) return;
        if (size<5) x[size++]=value;
    }
};
inline Roots roots(Polynomial polynomial) {
    polynomial.normalize();
    Roots result;
    if (!polynomial.degree) return result;
    if (polynomial.degree==1) { result.add(-polynomial.c[0]/polynomial.c[1]); return result; }
    auto critical=roots(polynomial.derivative());
    std::sort(critical.x.begin(),critical.x.begin()+critical.size);
    std::array<long double,7> bounds{};
    int count=0; bounds[count++]=0;
    for (int i=0;i<critical.size;++i) if (critical.x[i]>0 && critical.x[i]<1) bounds[count++]=critical.x[i];
    bounds[count++]=1;
    const auto tolerance=64*std::numeric_limits<long double>::epsilon();
    // A polynomial is monotone between successive roots of its derivative.
    // Isolate each sign-changing root independently, including both modes.
    for (int k=1;k<count;++k) {
        long double left=bounds[k-1],right=bounds[k],fl=polynomial.value(left),fr=polynomial.value(right);
        if (!fl || !fr || (fl<0)==(fr<0)) continue;
        for (int iteration=0;iteration<80;++iteration) {
            const auto mid=(left+right)/2;
            if (mid==left || mid==right) break;
            const auto fm=polynomial.value(mid);
            if (!fm) { left=right=mid; break; }
            if ((fl<0)==(fm<0)) { left=mid; fl=fm; }
            else right=mid;
        }
        result.add((left+right)/2);
    }
    // Add repeated roots after simple roots so roundoff near a flat critical
    // point cannot displace a sign-changing root from the fixed-size buffer.
    for (int k=1;k<count-1;++k)
        if (std::abs(polynomial.value(bounds[k]))<=tolerance) result.add(bounds[k]);
    return result;
}
} // namespace likelihood_detail

inline double minimize(double lo,double hi,double p,double q,const Counts& m) {
    using namespace likelihood_detail;
    if (lo==hi) return lo;
    const double weights[]={2*m(0,0)+m(0,1)+m(1,0),2*m(0,2)+m(0,1)+m(1,2),
                            2*m(2,0)+m(1,0)+m(2,1),2*m(2,2)+m(2,1)+m(1,2),m(1,1)};
    const auto curvature_bound=[](double a,double b,double total) {
        // a/x^2+b/(total-x)^2 >= (a+b+6*min(a,b))/total^2.
        return (a+b+6*std::min(a,b))/(total*total);
    };
    const double negative_bound=std::max(
        curvature_bound(weights[0],weights[1],p)+curvature_bound(weights[2],weights[3],1-p),
        curvature_bound(weights[0],weights[2],q)+curvature_bound(weights[1],weights[3],1-q));
    const double vertex=std::max(lo,std::min(hi,(2*p+2*q-1)/4));
    const double minimum_quad=vertex*(1-p-q+vertex)+(p-vertex)*(q-vertex);
    if (!weights[4] || (minimum_quad>0 && negative_bound>=4*weights[4]/minimum_quad)) {
        // Certified concavity permits a much cheaper, monotone score search.
        // Only ambiguous cases need the complete polynomial root isolation.
        double left=lo,right=hi,x=p*q;
        // Safeguarded Newton steps converge rapidly around ordinary LD values.
        // Difficult boundary cases finish with the original bracketed bisection.
        for (int iteration=0;iteration<80;++iteration) {
            double score=0,curvature=0;
            const double denominators[]={x,p-x,q-x,1-p-q+x};
            for (int i=0;i<4;++i) if (weights[i]) {
                const double term=weights[i]/denominators[i];
                score+=(i==0 || i==3) ? term : -term;
                curvature-=term/denominators[i];
            }
            if (weights[4]) {
                const double quad=x*(1-p-q+x)+(p-x)*(q-x);
                const double ratio=(4*x+1-2*p-2*q)/quad;
                score+=weights[4]*ratio;
                curvature+=weights[4]*(4/quad-ratio*ratio);
            }
            if (score==0) break;
            if (score>0) left=x; else right=x;
            double next=x-score/curvature;
            if (iteration<24 && std::isfinite(next) && next>=left && next<=right &&
                std::abs(next-x)<=2*std::numeric_limits<double>::epsilon()*std::max(std::abs(x),hi-lo)) {
                x=next; break;
            }
            if (iteration>=24 || !std::isfinite(next) || next<=left || next>=right ||
                std::abs(next-x)>(right-left)/2) next=(left+right)/2;
            if (next==x || next==left || next==right) { x=next; break; }
            x=next;
        }
        // A positive exponent on a vanishing haplotype makes the endpoint
        // likelihood zero. Ordinary interior solutions need no log evaluations.
        const bool lower_singular=(lo==0 && weights[0]>0) || (p+q>=1 && weights[3]>0);
        const bool upper_singular=(hi==p && weights[1]>0) || (hi==q && weights[2]>0);
        if (lower_singular && upper_singular) return x;
        double best=x,best_score=nll(x,p,q,m);
        if (!lower_singular) {
            const double score=nll(lo,p,q,m);
            if (score<best_score) { best=lo; best_score=score; }
        }
        if (!upper_singular && nll(hi,p,q,m)<best_score) best=hi;
        return best;
    }
    const long double span=static_cast<long double>(hi)-lo;
    const long double r=1-static_cast<long double>(p)-q;
    std::array<Polynomial,5> factors{};
    factors[0].c={{lo,span}}; factors[0].degree=1;
    factors[1].c={{static_cast<long double>(p)-lo,-span}}; factors[1].degree=1;
    factors[2].c={{static_cast<long double>(q)-lo,-span}}; factors[2].degree=1;
    factors[3].c={{std::max(0.0L,r+lo),span}}; factors[3].degree=1;
    factors[4].c={{lo*(r+lo)+(static_cast<long double>(p)-lo)*(static_cast<long double>(q)-lo),
                   span*(4*lo+1-2*static_cast<long double>(p)-2*static_cast<long double>(q)),2*span*span}};
    factors[4].degree=2;
    for (auto& factor:factors) factor.normalize();
    Polynomial numerator;
    for (int i=0;i<5;++i) if (weights[i]) {
        auto term=factors[i].derivative();
        for (int j=0;j<5;++j) if (j!=i && weights[j]) term=multiply(term,factors[j]);
        numerator.degree=std::max(numerator.degree,term.degree);
        for (int j=0;j<=term.degree;++j) numerator.c[j]+=weights[i]*term.c[j];
    }
    const auto candidates=roots(numerator);
    double best=lo,best_score=nll(lo,p,q,m);
    const auto consider=[&](double x) {
        const double score=nll(x,p,q,m);
        if (score<best_score) { best=x; best_score=score; }
    };
    consider(hi);
    for (int i=0;i<candidates.size;++i) consider(static_cast<double>(lo+span*candidates.x[i]));
    return best;
}
inline double evaluate_counts(const Counts& c,const Summary& sa,const Summary& sb,Method method,
                              const double* log_factorials=nullptr) {
    if (!c.n) return nan();
    if (method==Method::MI || method==Method::RMI)
        return information(c,sa,sb,method==Method::RMI,log_factorials);
    const double p=sa.p, q=sb.p;
    if (!(p>0 && p<1 && q>0 && q<1)) return method==Method::D ? 0 : nan();
    const double lo=std::max(0.0,p+q-1), hi=std::min(p,q);
    double x=minimize(lo,hi,p,q,c);
    const double D=x-p*q;
    if (method==Method::D) return D;
    if (method==Method::R2) return std::min(1.0,D*D/(p*q*(1-p)*(1-q)));
    const double bound=D>0 ? std::min(p*(1-q),q*(1-p)) : std::max(-p*q,-(1-p)*(1-q));
    return std::min(1.0,D/bound);
}
inline double pair(const int* a,const int* b,std::size_t n,
                   const Summary& sa,const Summary& sb,Method method) {
    return evaluate_counts(count(a,b,n),sa,sb,method);
}
struct Site { std::string chr, id; std::int64_t pos; };
class Window {
    typedef std::pair<std::int64_t,std::size_t> Entry;
    std::map<std::string,std::vector<Entry>> groups;
public:
    explicit Window(const std::vector<Site>& sites) {
        for (std::size_t i=0; i<sites.size(); ++i) groups[sites[i].chr].emplace_back(sites[i].pos,i);
        for (auto& group : groups) std::sort(group.second.begin(),group.second.end());
    }
    std::vector<std::size_t> partners(std::size_t i, const Site& site, std::int64_t distance) const {
        const auto& entries=groups.at(site.chr);
        auto first=std::lower_bound(entries.begin(),entries.end(),Entry(site.pos-distance,0));
        auto last=std::upper_bound(entries.begin(),entries.end(),Entry(site.pos+distance,std::numeric_limits<std::size_t>::max()));
        std::vector<std::size_t> result;
        for (; first!=last; ++first) if (first->second>i) result.push_back(first->second);
        std::sort(result.begin(),result.end());
        return result;
    }
};
inline std::int64_t position(const std::string& value) {
    std::int64_t coordinate=0;
    const auto integer=std::from_chars(value.data(),value.data()+value.size(),coordinate);
    if (integer.ec==std::errc() && integer.ptr==value.data()+value.size()) {
        if (coordinate>=0 && coordinate<=std::numeric_limits<std::int64_t>::max()/2) return coordinate;
        throw std::invalid_argument("Invalid site position: "+value);
    }
    // R may format ordinary numeric metadata in scientific notation.
    std::size_t used=0;
    const auto numeric=std::stold(value,&used);
    if (used!=value.size() || !std::isfinite(numeric) || numeric<0 ||
        std::floor(numeric)!=numeric || numeric>=std::ldexp(1.0L,62))
        throw std::invalid_argument("Invalid site position: "+value);
    return static_cast<std::int64_t>(numeric);
}
} // namespace GWLD
