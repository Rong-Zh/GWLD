// Compile against two header snapshots to compare packing without changing ISA flags.
#include "gwld/packed.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
struct Matrix {
    std::size_t n_rows,n_cols;
    std::vector<int> values;
    const int* colptr(std::size_t j) const { return values.data()+j*n_rows; }
};
int main(int argc,char** argv) {
    try {
        if (argc!=5) throw std::invalid_argument("Usage: packing_benchmark samples variants missing_percent repeats");
        const auto rows=std::stoul(argv[1]),cols=std::stoul(argv[2]);
        const int missing=std::stoi(argv[3]),repeats=std::stoi(argv[4]);
        if (!rows || cols<2 || missing<0 || missing>100 || repeats<1)
            throw std::invalid_argument("Invalid benchmark dimensions or settings");
        Matrix g{rows,cols,std::vector<int>(rows*cols)};
        std::mt19937 random(2026);
        for (auto& x:g.values) x=random()%100<static_cast<unsigned>(missing) ? -1 : static_cast<int>(random()%3);
        GWLD::PackedGenotypes packed;
        std::vector<double> times;
        for (int repeat=0;repeat<=repeats;++repeat) {
            const auto start=std::chrono::steady_clock::now();
            packed.build(g);
            const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            if (repeat) times.push_back(elapsed);
        }
        std::sort(times.begin(),times.end());
        const auto counts=packed.count(0,1);
        const auto expected=GWLD::count(g.colptr(0),g.colptr(1),rows);
        for (unsigned k=0;k<9;++k) if (counts.cells[k]!=expected.cells[k])
            throw std::runtime_error("Incorrect contingency counts");
        std::cout<<std::setprecision(17)<<times[times.size()/2]<<'\t'<<counts.n<<'\n';
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
