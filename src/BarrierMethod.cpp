//
// Created by Carlos 12/08/20
//

#include "../headers/BarrierMethod.h"

BarrierMethod::BarrierMethod(Graph *graph) {
  if (graph != nullptr) {
    this->graph = graph;
    initialize();
  } else exit(EXIT_FAILURE);
}

void BarrierMethod::initialize() {
  int o, d, n = graph->getN(), m = graph->getM();
  try {

    env.set("LogFile", "MS_mip.log");
    env.start();

    f = vector<vector<vector<GRBVar>>>(n, vector<vector<GRBVar>>(n, vector<GRBVar>(n)));
    y = vector<vector<GRBVar>>(n, vector<GRBVar>(n));
    z = vector<GRBVar>(n);

    char name[30];
    for (o = 0; o < n; o++) {
      for (auto *arc : graph->arcs[o]) {
	d = arc->getD();
	sprintf(name, "y_%d_%d", o, d);
	y[o][d] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);

	for (int k: graph->DuS) {
	  sprintf(name, "f_%d_%d_%d", o, d, k);
	  if (!graph->removedF[o][d][k])
	    this->f[o][d][k] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);
	  else this->f[o][d][k] = model.addVar(0.0, 0.0, 0, GRB_CONTINUOUS, name);
	}
      }
    }

    for (auto i : graph->terminals) {
      sprintf(name, "z_%d", i);
      z[i] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);
    }

    model.update();
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
    cout << ex.getErrorCode() << endl;
    exit(EXIT_FAILURE);
  }
}

void BarrierMethod::initModel() {
  cout << "Model Created!" << endl;
  objectiveFunction();
  rootFlow(), flowConservation(), terminalsFlow();
  relXandY(), maxArcs();
  limDelayAndJitter();
  limVariation();
  primeToTerminals();
  nonTerminalsLeafs();
}

void BarrierMethod::objectiveFunction() {
  GRBLinExpr objective;
  for (auto k : graph->terminals) objective += z[k];
  model.setObjective(objective, GRB_MINIMIZE);
  cout << "Objective Function was added successfully!" << endl;
}

void BarrierMethod::rootFlow() {
  int o, d, root = graph->getRoot();
  for (auto k : graph->terminals) {
    GRBLinExpr flowExpr, rootExpr;
    for (o = 0; o < graph->getN(); o++) {
      for (auto *arc : graph->arcs[o]) {
	d = arc->getD();
	if (o == root) flowExpr += f[root][d][k];
	else if (d == root) rootExpr += f[o][root][k];
      }
    }
    model.addConstr((flowExpr - rootExpr) == 1, "root_flow_all_" + to_string(k));
  }
  model.update();
  cout << "Flow on root node" << endl;
}

void BarrierMethod::flowConservation() {
  int o, d, root = graph->getRoot();
  for (auto k : graph->DuS) {
    for (int j = 0; j < graph->getN(); j++) {
      if (j != root && j != k) {
	GRBLinExpr flowIn, flowOut;
	for (o = 0; o < graph->getN(); o++) {
	  for (auto *arc : graph->arcs[o]) {
	    d = arc->getD();
	    if (o == j) flowOut += f[j][d][k];
	    if (d == j) flowIn += f[o][j][k];
	  }
	}
	model.addConstr((flowIn - flowOut) == 0, "flow_conservation_" + to_string(j) + "_" + to_string(k));
      }
    }
  }
  model.update();
  cout << "Flow conservation" << endl;
}

void BarrierMethod::terminalsFlow() {
  int o, d;
  for (auto k : graph->DuS) {
    GRBLinExpr flowIn, flowOut;
    for (o = 0; o < graph->getN(); o++) {
      for (auto *arc : graph->arcs[o]) {
	d = arc->getD();
	if (o == k) flowOut += f[k][d][k];
	if (d == k) flowIn += f[o][k][k];
      }
    }
    model.addConstr((flowOut - flowIn) == -1, "flow_on_terminals_" + to_string(k));
  }
  model.update();
  cout << "Flow on terminals" << endl;
}

void BarrierMethod::relXandY() {
  int o, d;
  for (o = 0; o < graph->getN(); o++) {
    for (auto *arc : graph->arcs[o]) {
      d = arc->getD();
      for (auto k : graph->DuS) {
	model.addConstr(f[o][d][k] <= y[o][d], "f_and_y_relation_" + to_string(o) + "_" + to_string(d) + "_" + to_string(k));
      }
    }
  }
  model.update();
  cout << "f and Y relation" << endl;
}

void BarrierMethod::maxArcs() {
  GRBLinExpr totalArcs;
  for (int o = 0; o < graph->getN(); o++) {
    for (auto *arc : graph->arcs[o]) {
      totalArcs += y[arc->getO()][arc->getD()];
    }
  }
  model.addConstr(totalArcs == (graph->getN() - 1), "maximum_of_arcs");
    
  model.update();
  cout << "maximum of arcs in the tree" << endl;
}

void BarrierMethod::limDelayAndJitter() {
  int o, d, paramDelay, paramJitter;
  paramDelay = graph->getParamDelay(), paramJitter = graph->getParamJitter();
  // Delay
  for (auto k : graph->terminals) {
    GRBLinExpr limDelay, limJitter;
    for (o = 0; o < graph->getN(); o++) {
      for (auto *arc : graph->arcs[o]) {
	d = arc->getD();
	limDelay += arc->getDelay() * f[o][d][k];
	limJitter += arc->getJitter() * f[o][d][k];
      }
    }
    model.addConstr(limDelay <= (paramDelay + (graph->getBigMDelay() - paramDelay) * z[k]), "delay_limit_" + to_string(k));
    model.addConstr(limJitter <= (paramJitter + (graph->getBigMJitter() - paramJitter) * z[k]), "jitter_limit_" + to_string(k));
  }
    
  model.update();
  cout << "Delay and Jitter limits" << endl;
}

void BarrierMethod::limVariation() {
  int o, d, bigMK, bigML;
  for (auto k : graph->terminals) {
    for (auto l : graph->terminals) {
      if (k != l) {
	GRBLinExpr delayVariation;
	for (o = 0; o < graph->getN(); o++) {
	  for (auto *arc : graph->arcs[o]) {
	    d = arc->getD();
	    delayVariation += arc->getDelay() * (f[o][d][k] - f[o][d][l]);
	  }
	}
	bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
	bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
	model.addConstr(delayVariation <= graph->getParamVariation() + bigMK * z[k] + bigML * z[l],
			"limit_of_variation_between_pairs_" + to_string(k) + "_" + to_string(l));
      }
    }
  }
  model.update();
  cout << "Delay variation limits" << endl;
}

void BarrierMethod::primeToTerminals() {
  for (auto k : graph->terminals)
    model.addConstr(z[k] >= f[0][k][k], "prime_to_terminals_" + to_string(k));
  model.update();
  cout << "S' to terminals" << endl;
}

void BarrierMethod::nonTerminalsLeafs() {
  model.addConstr(y[graph->getRoot()][0] == 1);
  for (auto q : graph->DuS) {
    for (auto e : graph->DuS) {
      if (e != q) {
	model.addConstr(f[0][q][e] <= 0, "non_terminals_leafs_" + to_string(q) + "_" + to_string(e));
      }
    }
  }
  model.update();
  cout << "Non terminals are leafs" << endl;
}

void BarrierMethod::solve() {
  try {
    model.set("TimeLimit", "3600.0");
    model.set("OutputFlag", "0");
    model.set(GRB_IntParam_PreDual, 1);
    model.set("Method", "2");
    model.set("Crossover", "0");
    model.update();
    model.write("barrier_model.lp");
    model.optimize();
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
    exit(0);
  }
}

void BarrierMethod::getMultipliersDelay(vector<double> &multipliersDelay) {
  for (auto k : graph->terminals){
    auto constr = model.getConstrByName("delay_limit_" + to_string(k));
    multipliersDelay[k] = abs(constr.get(GRB_DoubleAttr_Pi));
  }
}

void BarrierMethod::getMultipliersJitter(vector<double> &multipliersJitter){
  for (auto k : graph->terminals){
    auto constr = model.getConstrByName("jitter_limit_" + to_string(k));
    multipliersJitter[k] = abs(constr.get(GRB_DoubleAttr_Pi));
  }
}

void BarrierMethod::getMultipliersVariation(vector<vector<double>> &multipliersVar){
  for (auto k : graph->terminals) {
    for (auto l : graph->terminals) {
      if (k != l) {
	auto constr = model.getConstrByName("limit_of_variation_between_pairs_" + to_string(k) + "_" + to_string(l));
	multipliersVar[k][l] = abs(constr.get(GRB_DoubleAttr_Pi));
      }       
    }
  }
}

void BarrierMethod::getMultipliersRelation(vector<vector<vector<double>>> &multipliersRel){
  for (auto k : graph->DuS){
    for (int i = 0; i < graph->getN(); i++) {
      for (auto *arc : graph->arcs[i]) {
	auto constr = model.getConstrByName("f_and_y_relation_" + to_string(i) + "_" + to_string(arc->getD()) + "_" + to_string(k));
	multipliersRel[i][arc->getD()][k] = abs(constr.get(GRB_DoubleAttr_Pi));

      }
    }
  }

}

void BarrierMethod::getMultipliersLeaf(vector<vector<double>> &multipliersLeaf) {
  for (auto q : graph->DuS) {
    for (auto e : graph->DuS) {
      if (e != q) {
	auto constr = model.getConstrByName("non_terminals_leafs_" + to_string(q) + "_" + to_string(e));
	multipliersLeaf[q][e] = abs(constr.get(GRB_DoubleAttr_Pi));
      }       
    }
  }
}
