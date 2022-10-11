#include <cstdio>
#include <chrono>
#include <string>
#include "Hits.h"
#include "PageRank.h"
#include "InDegree.h"


int main()
{
    std::vector<std::string> file_names{ "./DataSet/web-NotreDame.txt", "./DataSet/web-Google.txt", "./DataSet/web-BerkStan.txt", "./DataSet/web-Stanford.txt" };
    for (auto f : file_names)
    {
        // name dataset
        std::cout << "Dataset: " << f.substr(10) << std::endl;
        // name algorithm
        std::cout << "-----PAGE RANK-----" << std::endl;
        // start
        auto start = std::chrono::high_resolution_clock::now();
        PageRank pk(f, 0.8);
        pk.compute();
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Time : " << (std::chrono::duration_cast<std::chrono::milliseconds>(end - start)).count() << " milliseconds" << std::endl;
        std::cout << "Number steps to converge: " << pk.n_convergence << std::endl;
        // name algorithm
        std::cout << "------HITS---------" << std::endl;
        // start
        start = std::chrono::high_resolution_clock::now();
        Hits h(f);
        h.compute();
        end = std::chrono::high_resolution_clock::now();
        std::cout << "Time : " << (std::chrono::duration_cast<std::chrono::milliseconds>(end - start)).count() << " milliseconds" << std::endl;
        std::cout << "Number steps to converge: " << h.n_convergence << std::endl;
        std::cout << "--------------------------------------------------------------------------\n";
    }
    return 0;
}