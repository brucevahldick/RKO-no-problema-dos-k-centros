#ifndef _PROBLEM_H
#define _PROBLEM_H

#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <cctype>

struct TProblemData
{
    int n;
    int k;

    std::vector<double> x;
    std::vector<double> y;

    std::vector<std::vector<double>> dist;
};

static inline std::string trim(const std::string &s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}
static inline std::string uppercase(const std::string &s) {
    std::string t = s;
    for (auto &c : t) c = std::toupper((unsigned char)c);
    return t;
}
static inline int find_first_digit(const std::string &s) {
    for (size_t i = 0; i < s.size(); ++i)
        if (std::isdigit((unsigned char)s[i]) || s[i] == '-' ) return (int)i;
    return -1;
}

void ReadData(char name[], TProblemData &data)
{
    std::ifstream in(name);
    if (!in.is_open())
    {
        printf("\nERROR: File (%s) not found!\n", name);
        exit(1);
    }

    std::string line;
    int dimension = -1;
    std::string edge_weight_type = "";

    while (std::getline(in, line))
    {
        std::string s = trim(line);
        if (s.empty()) continue;
        std::string up = uppercase(s);

        if (up.find("NODE_COORD_SECTION") != std::string::npos) break;

        if (up.find("DIMENSION") != std::string::npos)
        {
            int pos = find_first_digit(s);
            if (pos >= 0)
            {
                try { dimension = std::stoi(s.substr(pos)); } catch (...) { dimension = -1; }
                data.n = dimension;
            }
        }

        if (up.find("EDGE_WEIGHT_TYPE") != std::string::npos)
        {
            size_t colon = s.find(':');
            std::string val;
            if (colon != std::string::npos) val = trim(s.substr(colon + 1));
            else {
                std::istringstream iss(s);
                std::string tok;
                iss >> tok; // EDGE_WEIGHT_TYPE
                if (iss >> tok) val = tok;
            }
            edge_weight_type = uppercase(val);
        }
    }

    if (dimension <= 0)
    {
        printf("\nERROR: DIMENSION not found or invalid in file (%s)!\n", name);
        exit(1);
    }

    if (!edge_weight_type.empty() && edge_weight_type.find("EUC_2D") == std::string::npos)
    {
        printf("\nWARNING: EDGE_WEIGHT_TYPE is '%s' (not EUC_2D). Will still compute Euclidean distances.\n", edge_weight_type.c_str());
    }

    data.x.assign(dimension, 0.0);
    data.y.assign(dimension, 0.0);
    int readCount = 0;
    while (readCount < dimension && std::getline(in, line))
    {
        std::istringstream iss(line);
        int idx;
        double cx, cy;
        if (!(iss >> idx >> cx >> cy)) continue;
        if (idx >= 1 && idx <= dimension)
        {
            data.x[idx - 1] = cx;
            data.y[idx - 1] = cy;
            ++readCount;
        }
    }

    in.close();

    data.k = 5;

    data.dist.assign(dimension, std::vector<double>(dimension, 0.0));
    for (int i = 0; i < dimension; ++i)
    {
        for (int j = 0; j < dimension; ++j)
        {
            double dx = data.x[i] - data.x[j];
            double dy = data.y[i] - data.y[j];
            double d = std::sqrt(dx * dx + dy * dy);
            data.dist[i][j] = std::floor(d + 0.5);
        }
    }
}

double Decoder(TSol &s, const TProblemData &data)
{
    std::vector<std::pair<double,int>> rk_idx;
    for (int i = 0; i < data.n; ++i) rk_idx.push_back({s.rk[i], i});
    std::sort(rk_idx.begin(), rk_idx.end());

    int k = data.k;
    if (k > data.n) k = data.n;

    std::vector<int> centers;
    for (int i = 0; i < k; ++i) centers.push_back(rk_idx[i].second);

    double maxDist = 0.0;
    for (int v = 0; v < data.n; ++v)
    {
        double minDist = std::numeric_limits<double>::max();
        for (int c : centers) minDist = std::min(minDist, data.dist[v][c]);
        if (minDist < std::numeric_limits<double>::max()) maxDist = std::max(maxDist, minDist);
    }

    return maxDist;
}

void FreeMemoryProblem(TProblemData &data){
    data.x.clear();
    data.y.clear();
    data.dist.clear();
}

#endif
