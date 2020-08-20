//
// Created by Carlos on 09/08/20
//

#ifndef MS_HEU_H
#define MS_HEU_H

#include "Graph.h"
#include "edmonds_optimum_branching_impl.hpp"
#include "edmonds_optimum_branching.hpp"
#include "Arc.h"

class Heuristic {
  Graph *graph;
  BoostGraph heuristicGraph, edmonds;
  vector<Arc*> branchingEdges;
  property_map<BoostGraph, vertex_index_t>::type indexMap;
  property_map<BoostGraph, edge_weight_t>::type weightMap;
  
public:
  Heuristic(Graph *graph);

  void runEdmonds(vector<vector<vector<double>>> &multipliersRel, vector<vector<double>> &multipliersLeaf);

  void updateEdmonds(vector<vector<vector<double>>> &multipliersRel, vector<vector<double>> &multipliersLeaf);
  
  int initialHeuristic();

  int subgradientHeuristic(vector<vector<vector<double>>> &multipliersRel, vector<vector<double>> &multipliersLeaf);
};

#endif // MS_HEU_H
