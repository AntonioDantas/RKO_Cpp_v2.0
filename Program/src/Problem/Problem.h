#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Data.h"
#include <ranges>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <limits>

//---------------------- DEFINITION OF TYPES OF PROBLEM SPECIFIC --------------------

// Struct with node informations
struct TNode
{
    int id;
    double c; // Cost of security
    double p; // Probability
    double e; // Energy of UAV
};

// Struct to hold all problem data (formerly global variables)
struct TProblemData
{
    int n; // Size of the random-key vector (solution size)

    // Specific problem variables
    int nodeCount;
    int veichesCount;
    int vMin;
    int vMax;
    double f1Max;
    double f2Max;
    double total_probability;

    std::vector<TNode> node;                   // vector of nodes
    std::vector<std::vector<double>> dist;     // matrix with Euclidean distance
    std::vector<std::vector<double>> distNorm; // matrix with distance normalized
};

//------ DEFINITION OF GLOBAL CONSTANTS -----------------------------------------

const double PENALTY_FACTOR = 99999.0;

//----------------------- HELPER FUNCTIONS ---------------------------------------

// Find max value in matrix
double find_max_matrix(const std::vector<std::vector<double>> &m, size_t nRowsCols)
{
    double maxV = -std::numeric_limits<double>::infinity();
    size_t n = nRowsCols;
    for (size_t i = 0; i < n && i < m.size(); ++i)
    {
        for (size_t j = 0; j < n && j < m[i].size(); ++j)
        {
            double val = m[i][j];
            if (!std::isfinite(val))
                continue;
            if (std::fabs(val) <= 0.0)
                continue;
            if (val > maxV)
                maxV = val;
        }
    }
    return maxV;
}

// Normalize matrix in-place
void normalize_matrix(std::vector<std::vector<double>> &m, size_t nRowsCols)
{
    auto maxV = find_max_matrix(m, nRowsCols);

    if (!std::isfinite(maxV) || maxV <= 0.0)
        return;

    size_t n = nRowsCols;
    for (size_t i = 0; i < n && i < m.size(); ++i)
    {
        for (size_t j = 0; j < n && j < m[i].size(); ++j)
        {
            double val = m[i][j];
            if (!std::isfinite(val) || std::fabs(val) <= 0.0)
                continue;
            m[i][j] = val / maxV;
        }
    }
}

// Normalize specific member of NodeType
template <typename NodeType>
void normalize_member(std::vector<NodeType> &nodes, double NodeType::*member, size_t nCount)
{
    double maxV = -std::numeric_limits<double>::infinity();
    size_t n = std::min(nCount, nodes.size());

    for (size_t i = 0; i < n; ++i)
    {
        double val = nodes[i].*member;
        if (!std::isfinite(val))
            continue;
        if (val > maxV)
            maxV = val;
    }

    if (!std::isfinite(maxV) || maxV <= 0.0)
        return;

    for (size_t i = 0; i < n; ++i)
    {
        double val = nodes[i].*member;
        if (!std::isfinite(val))
            continue;
        nodes[i].*member = (val) / maxV;
    }
}

double calculate_route_distance(const std::vector<int> &route, const std::vector<std::vector<double>> &dist)
{
    double total_distance = 0.0;
    if (route.size() < 2)
        return 0.0;

    for (size_t i = 0; i < route.size() - 1; ++i)
    {
        int from_node = route[i];
        int to_node = route[i + 1];
        total_distance += dist[from_node][to_node];
    }
    return total_distance;
}

//----------------------- CORE PROBLEM FUNCTIONS ---------------------------------

/************************************************************************************
 Method: ReadData
 Description: read input data of the problem and populate TProblemData
*************************************************************************************/
void ReadData(char nameTable[], TProblemData &data)
{
    data.nodeCount = 0;
    data.veichesCount = 0;
    char name[400] = "../Instances/";
    strcat(name, nameTable);

    FILE *arq;
    arq = fopen(name, "r");

    if (arq == NULL)
    {
        printf("\nERROR: File (%s) not found!\n", name);
        getchar();
        exit(1);
    }

    // Read instance head
    char temp[100];
    fscanf(arq, "%s %d", temp, &data.nodeCount);
    fscanf(arq, "%s %d", temp, &data.veichesCount);
    fscanf(arq, "%s %d", temp, &data.vMin);
    fscanf(arq, "%s %d", temp, &data.vMax);
    fscanf(arq, "%s %lf", temp, &data.f1Max);
    fscanf(arq, "%s %lf", temp, &data.f2Max);

    printf("\nInstance: %s", nameTable);
    printf("\nNumber of Nodes: %d", data.nodeCount);
    printf("\nNumber of Vehicles: %d", data.veichesCount);
    printf("\nMinimum Visits: %d", data.vMin);
    printf("\nMaximum Visits: %d", data.vMax);
    printf("\nMax value of f1: %lf", data.f1Max);
    printf("\nMin value of f2: %lf\n", data.f2Max);

    fscanf(arq, "%s", temp);

    data.node.clear();
    TNode nodeTemp;

    while (data.node.size() < (size_t)data.nodeCount)
    {
        fscanf(arq, "%d %lf %lf %lf", &nodeTemp.id, &nodeTemp.c, &nodeTemp.p, &nodeTemp.e);
        data.node.push_back(nodeTemp);
    }

    // Read matrix of distances
    data.dist.clear();
    data.distNorm.clear();
    data.dist.resize(data.nodeCount, std::vector<double>(data.nodeCount));
    data.distNorm.resize(data.nodeCount, std::vector<double>(data.nodeCount));

    fscanf(arq, "%s", temp);

    for (int i = 0; i < data.nodeCount; i++)
    {
        for (int j = 0; j < data.nodeCount; j++)
        {
            if (fscanf(arq, "%lf", &data.dist[i][j]) != 1)
            {
                printf("Erro de leitura no elemento [%d][%d]\n", i, j);
            }
            data.distNorm[i][j] = data.dist[i][j];
        }
    }
    fclose(arq);

    data.total_probability = 0.0;
    for (int i = 0; i < data.nodeCount; i++)
    {
        if (data.node[i].e <= 0)
            data.total_probability += data.node[i].p;
    }

    // Normalizations
    normalize_matrix(data.distNorm, data.nodeCount);
    normalize_member(data.node, &TNode::c, data.nodeCount);

    // Set size of the solution vector (RK)
    data.n = data.nodeCount + (data.nodeCount - data.veichesCount) * 2;
}

/************************************************************************************
 Method: Decoder
 Description: mapping the random-key solutions into problem solutions
*************************************************************************************/
double Decoder(TSol &s, const TProblemData &data)
{
    if (s.rk.size() < (size_t)data.n)
    {
        printf("\nERROR: Size of rk vector (%zu) is less than required (%d)!\n", s.rk.size(), data.n);
        return PENALTY_FACTOR;
    }

    for (double val : s.rk)
    {
        if (std::isnan(val))
        {
            printf("RK value: %lf\n", val);
            return PENALTY_FACTOR;
        }
    }

    // Clean the last solution
    s.routes.clear();
    s.routes.reserve(data.veichesCount);
    
    s.ofv = 0.0;
    s.f1 = 0.0;
    s.f2 = 0.0;

    // Local Variables
    int last_vehicle_id = -1;
    std::vector<int> temp_target_pool;
    temp_target_pool.reserve(data.nodeCount);

    std::vector<int> processing_order(data.nodeCount);
    std::vector<int> extra_visits_pool;
    extra_visits_pool.reserve(data.nodeCount * 2);

    // Base on rk order
    std::iota(processing_order.begin(), processing_order.end(), 0);
    std::sort(processing_order.begin(), processing_order.end(),
              [&](int a, int b)
              { return s.rk[a] < s.rk[b]; });

    // Empty routes (Start vehicles at depot/starting nodes)
    for (int i = 0; i < data.nodeCount; ++i)
    {
        if (data.node[i].e > 0)
        {
            s.routes.push_back({i});
        }
    }

    // Group nodes in routes
    int count_check_visit = 0;
    for (int node_id : processing_order)
    {
        if (data.node[node_id].e == 0) // Target node
        {
            double visit_key = s.rk[data.nodeCount + count_check_visit];
            count_check_visit++;

            // Assumes data.vMin/data.vMax are defined in Data.h
            bool is_visited = (data.vMin + floor(visit_key * (data.vMax - data.vMin + 1.0))) > 0;

            if (!is_visited)
                continue;

            temp_target_pool.push_back(node_id);
        }
        else // Vehicle/Depot node
        {
            auto it = std::find_if(s.routes.begin(), s.routes.end(),
                                   [&](const std::vector<int> &route)
                                   { return route[0] == node_id; });
            if (it != s.routes.end())
            {
                it->insert(it->end(), temp_target_pool.begin(), temp_target_pool.end());
            }
            temp_target_pool.clear();
            last_vehicle_id = node_id;
        }
    }

    // Insert pool of targets at the last vehicle
    if (!temp_target_pool.empty() && last_vehicle_id != -1)
    {
        auto it = std::find_if(s.routes.begin(), s.routes.end(),
                               [&](const std::vector<int> &route)
                               { return route[0] == last_vehicle_id; });
        if (it != s.routes.end())
        {
            it->insert(it->end(), temp_target_pool.begin(), temp_target_pool.end());
        }
    }

    // Add extra visits to the pool based on rk values
    int count_extra_visits = 0;
    for (int i = 0; i < data.nodeCount; ++i)
    {
        if (data.node[i].e == 0)
        {
            double revisit_key = s.rk[data.nodeCount + count_check_visit + count_extra_visits];
            count_extra_visits++;

            int extra_visits = floor(revisit_key * (data.vMax - data.vMin));
            for (int k = 0; k < extra_visits; ++k)
            {
                extra_visits_pool.push_back(i);
            }
        }
    }

    // Order extra visits by rk value
    std::sort(extra_visits_pool.begin(), extra_visits_pool.end(),
              [&](int a, int b)
              { return s.rk[a] < s.rk[b]; });

    // Distribute extra visits to vehicles
    if (!extra_visits_pool.empty())
    {
        for (size_t i = 0; i < extra_visits_pool.size(); ++i)
        {
            int target_to_add = extra_visits_pool[i];

            // Find route with minimum distance increase (or total distance)
            auto it = std::min_element(s.routes.begin(), s.routes.end(),
                                       [&data](const std::vector<int> &routeA, const std::vector<int> &routeB)
                                       {
                                           return calculate_route_distance(routeA, data.dist) < calculate_route_distance(routeB, data.dist);
                                       });

            if (it != s.routes.end())
            {
                if (it->back() != target_to_add)
                {
                    it->push_back(target_to_add);
                }
            }
        }
    }

    // Validate and calculate objective functions
    double total_energy_overflow = 0.0;
    double total_distance = 0.0;
    int total_visits = 0;
    std::vector<int> visit_counts(data.nodeCount, 0);

    for (const auto &route : s.routes)
    {
        if (route.size() <= 1)
            continue;

        double current_distance = 0.0;
        double current_distance_norm = 0.0;

        for (size_t i = 1; i < route.size(); ++i)
        {
            int from_node = route[i - 1];
            int to_node = route[i];
            current_distance += data.dist[from_node][to_node];
            current_distance_norm += data.distNorm[from_node][to_node];

            s.f2 += data.node[to_node].p;
            visit_counts[to_node]++;
            total_visits++;
        }
        s.f1 += current_distance_norm;
        total_distance += current_distance;

        // Check energy constraint
        if (current_distance > data.node[route[0]].e)
        {
            total_energy_overflow += (current_distance - data.node[route[0]].e);
        }

        // Check F1 constraint (assumes data.f1Max is defined in Data.h)
        if (data.f2Max == 0 && data.f1Max > 0 && current_distance > data.f1Max)
        {
            total_energy_overflow += (current_distance - data.f1Max);
        }
    }

    // Objective function with penalty if necessary
    // Assumes data.f1Max, data.f2Max are defined in Data.h
    if (data.f1Max > 0 && data.f2Max > 0)
    {
        s.ofv = (data.f1Max * s.f1) - (data.f2Max * s.f2);
    }
    else
    {
        if (data.f1Max == 0)
        {
            s.ofv = s.f1;

            if (s.f2 < data.f2Max * data.total_probability)
                s.ofv += PENALTY_FACTOR * (data.f2Max * data.total_probability - s.f2);
        }
        else
            s.ofv = -1 * s.f2;
    }

    if (total_energy_overflow > 0)
    {
        s.ofv += PENALTY_FACTOR * total_energy_overflow;
    }

    //// Debug print (assumes 'debug', 'print', 'arqSol' are available from Data.h/Global context
    //// or passed implicitly. Keeping structure from original code)
    // if (debug && print)
    //{
    //     printf("\nConstructed Routes: \n");
    //     for (const auto &route : s.routes)
    //     {
    //         printf("Vehicle %d (%f) -> ", data.node[route[0]].id, data.node[route[0]].e);
    //         double current_distance = 0.0;
    //         for (size_t i = 1; i < route.size(); ++i)
    //         {
    //             current_distance += data.dist[route[i - 1]][route[i]];
    //             printf("%d ", data.node[route[i]].id);
    //         }
    //         printf(" = %f \n", current_distance);
    //     }
    //     printf("\nTotal distance: %f \n", total_distance);
    //     if (total_visits > 0)
    //          printf("\nProbabilty AVG: %f \n", s.f2 / total_visits);
    // }

    // if (!debug && print)
    //{
    //     for (const auto &route : s.routes)
    //     {
    //         fprintf(arqSol, "%d ", data.node[route[0]].id);
    //         for (size_t i = 1; i < route.size(); ++i)
    //         {
    //             fprintf(arqSol, "%d ", data.node[route[i]].id);
    //         }
    //     }
    // }

    return s.ofv;
}

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem(TProblemData &data)
{
    data.dist.clear();
    data.distNorm.clear();
    data.node.clear();
}

#endif