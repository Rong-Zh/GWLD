#pragma once

#include <algorithm>
#include <cmath>
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
struct Summary {
    bool present[3] = {false, false, false};
    double p = nan();
};
inline Summary summarize(const int* g, std::size_t n) {
    Summary s;
    double counts[3] = {};
    for (std::size_t k = 0; k < n; ++k) {
        if (g[k] == -1 || g[k] == std::numeric_limits<int>::min()) continue;
        if (g[k] < 0 || g[k] > 2)
            throw std::invalid_argument("Genotypes must be 0, 1, 2 or missing (-1/NA)");
        s.present[g[k]] = true;
        ++counts[g[k]];
    }
    const double total = counts[0] + counts[1] + counts[2];
    if (total) s.p = (counts[0] + 0.5 * counts[1]) / total;
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
inline double information(const Counts& c, const Summary& a, const Summary& b, bool reduced) {
    if (!c.n) return nan();
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
    double I = std::lgamma(c.n+1), sx = 0, sy = 0, lx = 0, ly = 0;
    for (double v : c.cells) I += std::lgamma(v+1);
    const double w = c.n/(c.n+0.5*R*S);
    for (int i = 0; i < 3; ++i) {
        if (a.present[i]) {
            I -= std::lgamma(c.rows[i]+1);
            const double x = (1-w)/R+w*c.rows[i]/c.n;
            sx += x*x; lx += std::log(x);
        }
        if (b.present[i]) {
            I -= std::lgamma(c.cols[i]+1);
            const double y = (1-w)/S+w*c.cols[i]/c.n;
            sy += y*y; ly += std::log(y);
        }
    }
    const double nu = (S+1)/(S*sx)-1/S, mu = (R+1)/(R*sy)-1/R;
    const double omega = (R-1)*(S-1)*std::log(c.n+0.5*R*S)
        +0.5*(R+nu-2)*ly+0.5*(S+mu-2)*lx
        +0.5*(std::lgamma(mu*R)+std::lgamma(nu*S)
              -R*(std::lgamma(S)+std::lgamma(mu))-S*(std::lgamma(R)+std::lgamma(nu)));
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
// Bounded Brent minimization, retaining the v1 estimator's tolerance.
// Based on the Richard Brent / John Burkardt implementation used in v1.
inline double minimize(double a, double b, double pA, double pB, const Counts& m) {
    const double golden = (3-std::sqrt(5.0))*0.5;
    const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
    const double tol = std::pow(std::numeric_limits<double>::epsilon(), 0.25)/3;
    double x = a+golden*(b-a), w=x, v=x, d=0, e=0;
    double fx=nll(x,pA,pB,m), fw=fx, fv=fx;
    for (int iteration=0; iteration<1000; ++iteration) {
        const double mid=(a+b)*0.5, t=eps*std::abs(x)+tol;
        if (std::abs(x-mid) <= 2*t-(b-a)*0.5) break;
        double p=0, q=0, r=0;
        if (std::abs(e)>t) {
            r=(x-w)*(fx-fv); q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r; q=2*(q-r);
            if (q>0) p=-p; else q=-q;
            r=e; e=d;
        }
        if (std::abs(p)>=std::abs(0.5*q*r) || p<=q*(a-x) || p>=q*(b-x)) {
            e=x<mid ? b-x : a-x; d=golden*e;
        } else {
            d=p/q;
            if (x+d-a<2*t || b-x-d<2*t) d=x>=mid ? -t : t;
        }
        const double u=x+(std::abs(d)>=t ? d : (d>0 ? t : -t));
        const double fu=nll(u,pA,pB,m);
        if (fu<=fx) {
            if (u<x) b=x; else a=x;
            v=w; fv=fw; w=x; fw=fx; x=u; fx=fu;
        } else {
            if (u<x) a=u; else b=u;
            if (fu<=fw || w==x) { v=w; fv=fw; w=u; fw=fu; }
            else if (fu<=fv || v==x || v==w) { v=u; fv=fu; }
        }
    }
    return x;
}
inline double pair(const int* a, const int* b, std::size_t n,
                   const Summary& sa, const Summary& sb, Method method) {
    const Counts c=count(a,b,n);
    if (!c.n) return nan();
    if (method==Method::MI || method==Method::RMI)
        return information(c,sa,sb,method==Method::RMI);
    const double p=sa.p, q=sb.p;
    if (!(p>0 && p<1 && q>0 && q<1)) return method==Method::D ? 0 : nan();
    const double lo=std::max(0.0,p+q-1), hi=std::min(p,q);
    double x=minimize(lo,hi,p,q,c);
    // Compare feasible endpoints explicitly (perfect LD can attain an endpoint).
    if (nll(lo,p,q,c)<nll(x,p,q,c)) x=lo;
    if (nll(hi,p,q,c)<nll(x,p,q,c)) x=hi;
    const double D=x-p*q;
    if (method==Method::D) return D;
    if (method==Method::R2) return std::min(1.0,D*D/(p*q*(1-p)*(1-q)));
    const double bound=D>0 ? std::min(p*(1-q),q*(1-p)) : std::max(-p*q,-(1-p)*(1-q));
    return std::min(1.0,D/bound);
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
    std::size_t used=0;
    const auto pos=std::stold(value,&used);
    if (used!=value.size() || !std::isfinite(pos) || pos<0 || std::floor(pos)!=pos || pos>std::numeric_limits<std::int64_t>::max()/2)
        throw std::invalid_argument("Invalid site position: " + value);
    return static_cast<std::int64_t>(pos);
}
} // namespace GWLD
