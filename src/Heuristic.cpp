//
// Created by Carlos on 08/08/2020
//

#include "../headers/Heuristic.h"

Heuristic::Heuristic(Graph *graph){
  Heuristic::graph = graph;
  int i, j, n = graph->getN();
  
  heuristicGraph = BoostGraph(n);
  edmonds = BoostGraph(n);

  for (i = 0; i < n; i++) {
    for (auto *arc : graph->arcs[i]) {
      j = arc->getD();
      add_edge(i, j, arc->getDelay(), heuristicGraph);
      add_edge(i, j, 0, edmonds);
    }
  }
}

int Heuristic::initialHeuristic() {
    int n = graph->getN();
    property_map<BoostGraph, edge_weight_t>::type weightMapDelay = get(edge_weight, heuristicGraph);
    vector<VertexDescriptor> predecessors = vector<VertexDescriptor>(n);
    vector<int> distance = vector<int>(n);
    vector<bool> notAttended = vector<bool>(n);
    vector<int> jitterDistance = vector<int>(n);

    dijkstra_shortest_paths(heuristicGraph, graph->getRoot(), predecessor_map(
            make_iterator_property_map(predecessors.begin(), get(vertex_index, heuristicGraph))).distance_map(
            make_iterator_property_map(distance.begin(), get(vertex_index, heuristicGraph))));

    double meanValue = 0;
    int countTerm = 0, actual, jitter, count = 0;

    for (auto k : graph->terminals) 
        if (distance[k] > graph->getParamDelay()) {
            countTerm++;
            notAttended[k] = true;
        } else meanValue += distance[k];

    for (auto t : graph->terminals) {
        actual = t, jitter = 0;
        while (actual != graph->getRoot()) {
            jitter += graph->getJitter(predecessors[actual], actual);
            actual = predecessors[actual];
        }
        if (jitter > graph->getParamJitter()) 
            notAttended[t] = true;
        jitterDistance[t] = jitter;
    }

    meanValue = meanValue / double(graph->terminals.size() - countTerm);
    int lessThanAvg = 0, greaThanAvg = 0;
    
    for (auto k : graph->terminals) {
        if (!notAttended[k]) {
            if (distance[k] <= meanValue) lessThanAvg++;
            else greaThanAvg++;
        }
    }

    for (auto k : graph->terminals) {
        if (!notAttended[k]) {
            for (auto l : graph->terminals) {
                if (l != k && !notAttended[l]) {
                    if (distance[k] - distance[l] > graph->getParamVariation()) {
                        if (lessThanAvg < greaThanAvg) {
                            if (distance[k] < distance[l]) {
                                notAttended[k] = true;
                                break;
                            } else notAttended[l] = true;
                        } else {
                            if (distance[k] > distance[l]) {
                                notAttended[k] = true;
                                break;
                            } else notAttended[l] = true;
                        }
                    }
                }
            }
        }
    }

    for (auto t : graph->terminals) 
        if (notAttended[t]) count++;  
    return count;
}

void Heuristic::updateEdmonds(vector<vector<vector<double>>> &multipliersRel, vector<vector<double>> &multipliersLeaf) {
  EdgeDescriptor e;
  bool found;
  int n = graph->getN();
  double edgeCost;
  
  for (auto *arc : graph->arcs[0]) {
    edgeCost = 0.0;
    for (auto k : graph->DuS) {
      edgeCost += multipliersRel[0][arc->getD()][k];
      if (k != arc->getD())
	edgeCost += multipliersLeaf[arc->getD()][k];
    }
    tie(e, found) = edge(0, arc->getD(), edmonds);
    if (found) boost::put(edge_weight_t(), edmonds, e, edgeCost);
  }

  for (int i = 1; i < n; i++) {
    for (auto *arc : graph->arcs[i]) {
      edgeCost = 0.0;
      for (auto k : graph->DuS) edgeCost += multipliersRel[i][arc->getD()][k];

      tie(e, found) = edge(i, arc->getD(), edmonds);
      if (found) boost::put(edge_weight_t(), edmonds, e, edgeCost);
    }
  }
}

void Heuristic::runEdmonds(vector<vector<vector<double>>> &multipliersRel, vector<vector<double>> &multipliersLeaf) {
  updateEdmonds(multipliersRel, multipliersLeaf);
  int n = graph->getN();
  int i, j;
  vector<EdgeDescriptor> branching;
  VertexDescriptor roots[2];
  roots[0] = graph->getRoot();
  roots[1] = graph->getRoot();

  edmonds_optimum_branching<false, true, true>(edmonds,
					       indexMap,
					       weightMap,
					       roots,
					       roots+1,
					       back_inserter(branching));

  for (auto e : branching) {
    i = e.m_source, j = e.m_target;
    branchingEdges.push_back(new Arc(i, j, graph->getDelay(i, j), graph->getJitter(i, j), 0, 0));
  }
}

int Heuristic::subgradientHeuristic(vector<vector<vector<double>>> &multipliersRel, vector<vector<double>> &multipliersLeaf) {
  runEdmonds(multipliersRel, multipliersLeaf);
   
  int obj = 0, i, j, n = graph->getN(), bestCand, root = graph->getRoot();
  vector<int> delayPaths = vector<int>(n), jitterPaths = vector<int>(n),
    predecessors = vector<int>(n), delayPathAux = vector<int>(n),
    jitterPathAux = vector<int>(n), predecessorsAux = vector<int>(n);
  vector<bool> notAttended = vector<bool>(n), notAttendedAux = vector<bool>(n);
  vector<int> changed = vector<int>();

  predecessors[root] = root;
  for (auto arc : branchingEdges) {
    predecessors[arc->getD()] = arc->getO();
  }

  random_shuffle(graph->terminals.begin(), graph->terminals.end());

  int actual, jitter, delay, count = 0;
  for (auto k : graph->DuS) {
    actual = k, jitter = 0, delay = 0;
    while (actual != root) {
      jitter += graph->getJitter(predecessors[actual], actual),
	delay += graph->getDelay(predecessors[actual], actual);
      actual = predecessors[actual];
    }
    delayPaths[k] = delay, jitterPaths[k] = jitter;
  }

  // obj = 0;
  for (auto k : graph->terminals) 
    if (delayPaths[k] > graph->getParamDelay() || jitterPaths[k] > graph->getParamJitter())
      notAttended[k] = true, obj++;

  // Select a terminal to fix as attend
  int diffDelay, diffJitter;
  bool canMove;
  int selected = -1, losts;

  for (auto k : graph->terminals) {
    if (!notAttended[k]) {
      selected = k;
      break;
    }
  }

  if (selected == -1) 
    return graph->terminals.size();

  // get the values of this path
  for (auto k : graph->terminals) {
    if (k != selected) {
      bestCand = -1, delay = graph->getParamDelay()+1, jitter = graph->getParamJitter()+1; 
            
      if (delayPaths[k] < (delayPaths[selected] - graph->getParamVariation()) ||
	  delayPaths[k] > (delayPaths[selected] + graph->getParamVariation())) {
                
	for (auto arc : graph->arcs[k]) {
	  i = arc->getD();
	  // Get the best candidate to move
	  if (i != predecessors[k] && k != predecessors[i]) {
	    if (delayPaths[i] + arc->getDelay() >= (delayPaths[selected] - graph->getParamVariation()) &&
		delayPaths[i] + arc->getDelay() <= (delayPaths[selected] + graph->getParamVariation()) &&
		delayPaths[i] + arc->getDelay() <= graph->getParamDelay() && 
		jitterPaths[i] + arc->getJitter() <= graph->getParamJitter()) {

	      bestCand = i;
	      delay = delayPaths[i] + arc->getDelay();
	      jitter = jitterPaths[i] + arc->getJitter();
	      break;
	    }
	  }
	}

	if (bestCand != -1) {
	  // create the temporary vectors
	  notAttended[k] = false;
	  canMove = true;
	  losts = 0;
	  for (i = 0; i < n; i++) {
	    delayPathAux[i] = delayPaths[i], jitterPathAux[i] = jitterPaths[i], 
	      predecessorsAux[i] = predecessors[i], notAttendedAux[i] = notAttended[i];
	  }

	  // Evaluate the move                    
	  diffDelay = delay - delayPathAux[k], diffJitter = jitter - jitterPathAux[k];
	  predecessorsAux[k] = bestCand;
	  delayPathAux[k] = delay, jitterPathAux[k] = jitter;

	  changed.erase(changed.begin(), changed.end());
	  changed.push_back(k);
	  while (!changed.empty()) {
	    actual = changed.back();
	    changed.pop_back();
	    for (int j : graph->DuS) {
	      if (j != actual && predecessorsAux[j] == actual) {
		delayPathAux[j] += diffDelay, jitterPathAux[j] += diffJitter;
		if (delayPathAux[j] > graph->getParamDelay() || 
		    jitterPathAux[j] > graph->getParamJitter()
		    /*		    delayPathAux[j] < (delayPaths[selected] - graph->getParamVariation()) ||
				    delayPathAux[j] > (delayPaths[selected] + graph->getParamVariation())*/) {
		  notAttendedAux[j] = true;
		  losts++;
		} else notAttendedAux[j] = false;

		if (losts >= 2 || (j == selected && notAttendedAux[j])) {
		  canMove = false;
		  changed.erase(changed.begin(), changed.end());
		  break;
		} else changed.push_back(j);
	      }
	    }
	  }

	  if (canMove) {
	    for (i = 0; i < n; i++) {
	      delayPaths[i] = delayPathAux[i], jitterPaths[i] = jitterPathAux[i], 
		predecessors[i] = predecessorsAux[i], notAttended[i] = notAttendedAux[i];
	    }
	  }
	}  
      } 
    }
  }

  int leqSelec = 0, geqSelec = 0; 
  for (auto k : graph->terminals)
    if (k != selected) {
      if (delayPaths[k] < (delayPaths[selected] - graph->getParamVariation()) ||
	  delayPaths[k] > (delayPaths[selected] + graph->getParamVariation()) ||
	  delayPaths[k] > graph->getParamDelay() || jitterPaths[k] > graph->getParamJitter())
	notAttended[k] = true;
      else notAttended[k] = false;

      if (!notAttended[k]) {
	if (delayPaths[k] <= delayPaths[selected]) leqSelec++;
	else geqSelec++;
      }
    } 
        
  for (auto k : graph->terminals) {
    if (k != selected && !notAttended[k]) {
      for (auto l : graph->terminals) {
	if (k != l && l != selected && !notAttended[l]) {
	  if (delayPaths[k] - delayPaths[l] > graph->getParamVariation()) {
	    if (leqSelec > geqSelec) {
	      if (delayPaths[k] > delayPaths[selected]) {
		notAttended[k] = true;
		break;
	      } else notAttended[l] = true;
	    } else {
	      if (delayPaths[k] < delayPaths[selected]) {
		notAttended[k] = true;
		break;
	      } else notAttended[l] = true;
	    }
	  }
	}
      }
    }
  }

  obj = 0;
  for (auto k : graph->terminals)
    if (notAttended[k]) obj++;
    
  return obj;
}
