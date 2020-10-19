//
// Created by carlos on 06/03/19.
//

#include <chrono>
#include "../headers/Model.h"

Model::Model(Graph *graph, int relaxNum, bool heuristics, bool barrierMethod, double lambda, int maxIter, int B, int time) {
  Model::graph = graph;
  Model::lambda = lambda;
  Model::maxIter = maxIter;
  Model::B = B;
  Model::time = time;
  Model::relaxNum = relaxNum;
  Model::heuristics = heuristics;
  Model::barrierM = barrierMethod;

  LB = 0, UB = 0, iter = 0;
  if (heuristics) heuristic = new Heuristic(graph); 
}

void Model::initialize() {
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
	y[o][d] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
	for (int k: graph->DuS) {
	  sprintf(name, "f_%d_%d_%d", o, d, k);
	  if (graph->removedF[o][d][k])
	    this->f[o][d][k] = model.addVar(0.0, 0.0, 0, GRB_BINARY, name); 
	  else
	    this->f[o][d][k] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
	}
      }
    }

    for (auto i : graph->terminals) {
      sprintf(name, "z_%d", i);
      z[i] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
    }
    
    model.update();
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
    cout << ex.getErrorCode() << endl;
    exit(EXIT_FAILURE);
  }
}

void Model::initModelShp() {
  setObjRL1();
  rootFlow(), flowConservation(), terminalsFlow(); 
  maxArcs(), primeToTerminals();
  cout << "Relaxation 01 - Model Created" << endl;
}

void Model::initModel1ResDelayCshp() {
  setObjRL2();
  rootFlow(), flowConservation(), terminalsFlow(); 
  maxArcs(), limDelay(), primeToTerminals();
  cout << "Relaxation 02 - Model Created" << endl;
}

void Model::initModel1ResJitterCshp() {
  setObjRL3();
  rootFlow(), flowConservation(), terminalsFlow(); 
  maxArcs(), limJitter(), primeToTerminals();
  cout << "Relaxation 03 - Model Created" << endl;
}

void Model::initModel2ResCshp() {
  setObjRL4();
  rootFlow(), flowConservation(), terminalsFlow(); 
  maxArcs(), limDelay(), limJitter(), primeToTerminals();
  cout << "Relaxation 04 - Model Created" << endl;
}

void Model::setObjRL1() {
  GRBLinExpr objective;
  int i, j, o, bigMK, bigML, n = graph->getN();

  for (auto k : graph->terminals) {
    for (i = 0; i < n; i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	objective += multipliersDelay[k] * (arc->getDelay() * f[i][j][k] - (graph->getParamDelay() + z[k] * graph->getBigMDelay()));
	objective += multipliersJitter[k] * (arc->getJitter() * f[i][j][k] - (graph->getParamJitter() + z[k] * graph->getBigMJitter()));
      }
    }
  }

  for (auto k : graph->terminals) {
    objective += z[k];
    for (auto l : graph->terminals) {
      if (k != l) {
	bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
	bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);	
	for (i = 0; i < n; i++) {
	  for (auto *arc : graph->arcs[i]) {
	    j = arc->getD();
	    objective += multipliersVar[k][l] * ((arc->getDelay() * (f[i][j][k] - f[i][j][l])) - (graph->getParamVariation() + bigMK * z[k] + bigML * z[l]));
	  }
	}
      }
    }
  }

  for (auto q : graph->DuS) {
    for (i = 0; i < n; i++) 
      for (auto *arc : graph->arcs[i]) objective += multipliersRel[i][arc->getD()][q] * (f[i][arc->getD()][q] - y[i][arc->getD()]);
    for (auto e : graph->DuS) 
      if (q != e) objective += multipliersLeaf[q][e] * (f[0][q][e]);
  }

  model.setObjective(objective, GRB_MINIMIZE);
  model.update();
  cout << "Objective function RL01" << endl;
}

void Model::setObjRL2() {
  GRBLinExpr objective;
  int i, j, o, bigMK, bigML, n = graph->getN();
  double auxVariation = 0, auxVariationMulti = 0;

  for (auto k : graph->terminals) {
    for (i = 0; i < n; i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	objective += multipliersJitter[k] * (arc->getJitter() * f[i][j][k] - (graph->getParamJitter() + z[k] * graph->getBigMJitter()));
      }
    }
  }

  for (auto k : graph->terminals) {
    objective += z[k];
    for (auto l : graph->terminals) {
      if (k != l) {
	bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
	bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);	
	for (i = 0; i < n; i++) {
	  for (auto *arc : graph->arcs[i]) {
	    j = arc->getD();
	    objective += multipliersVar[k][l] * ((arc->getDelay() * (f[i][j][k] - f[i][j][l])) - (graph->getParamVariation() + bigMK * z[k] + bigML * z[l]));
	  }
	}
      }
    }
  }

  for (auto q : graph->DuS) {
    for (i = 0; i < n; i++) 
      for (auto *arc : graph->arcs[i]) objective += multipliersRel[i][arc->getD()][q] * (f[i][arc->getD()][q] - y[i][arc->getD()]);
    for (auto e : graph->DuS) 
      if (q != e) objective += multipliersLeaf[q][e] * (f[0][q][e]);
  }

  model.setObjective(objective, GRB_MINIMIZE);
  model.update();
  cout << "Objective function" << endl;
}

void Model::setObjRL3() {  
  GRBLinExpr objective;
  int j, i, o, bigMK, bigML, n = graph->getN();
  double auxVariation = 0, auxVariationMulti = 0;

  for (auto k : graph->terminals) {
    for (i = 0; i < n; i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	objective += multipliersDelay[k] * (arc->getDelay() * f[i][j][k] - (graph->getParamDelay() + z[k] * graph->getBigMDelay()));
      }
    }
  }

  for (auto k : graph->terminals) {
    objective += z[k];
    for (auto l : graph->terminals) {
      if (k != l) {
	bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
	bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);	
	for (i = 0; i < n; i++) {
	  for (auto *arc : graph->arcs[i]) {
	    j = arc->getD();
	    objective += multipliersVar[k][l] * ((arc->getDelay() * (f[i][j][k] - f[i][j][l])) - (graph->getParamVariation() + bigMK * z[k] + bigML * z[l]));
	  }
	}
      }
    }
  }

  for (auto q : graph->DuS) {
    for (i = 0; i < n; i++) 
      for (auto *arc : graph->arcs[i]) objective += multipliersRel[i][arc->getD()][q] * (f[i][arc->getD()][q] - y[i][arc->getD()]);
    for (auto e : graph->DuS) 
      if (q != e) objective += multipliersLeaf[q][e] * (f[0][q][e]);
  }
  
  model.setObjective(objective, GRB_MINIMIZE);
  model.update();
  cout << "Objective function RL3" << endl;
}

void Model::setObjRL4() {
  GRBLinExpr objective;
  int j, o, i, bigMK, bigML, n = graph->getN();
  double auxVariation = 0, auxVariationMulti = 0;

  for (auto k : graph->terminals) {
    objective += z[k];
  }

  for (auto k : graph->DuS)  {
    for (int i = 0; i < graph->getN(); i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	objective += multipliersRel[i][j][k] * (f[i][j][k] - y[i][j]);
      }
    }
  }

  for (auto q : graph->DuS) {
    for (auto e : graph->DuS) {
      if (q != e) objective += multipliersLeaf[q][e] * f[0][q][e];
    }
  }
  
  /*
  for (auto k : graph->terminals) {
    objective += z[k];
    for (auto l : graph->terminals) {
      if (k != l) {
	bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
	bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);	
	for (i = 0; i < n; i++) {
	  for (auto *arc : graph->arcs[i]) {
	    j = arc->getD();
	    objective += multipliersVar[k][l] * ((arc->getDelay() * (f[i][j][k] - f[i][j][l])) - (graph->getParamVariation() + bigMK * z[k] + bigML * z[l]));
	  }
	}
      }
    }
  }
  
  for (auto q : graph->DuS) {
    for (i = 0; i < n; i++) 
      for (auto *arc : graph->arcs[i]) objective += multipliersRel[i][arc->getD()][q] * (f[i][arc->getD()][q] - y[i][arc->getD()]);
    for (auto e : graph->DuS) 
      if (q != e) objective += multipliersLeaf[q][e] * (f[0][q][e]);
  }
  */
  model.setObjective(objective, GRB_MINIMIZE);
  model.update();
  //cout << "Objective function" << endl;
}

void Model::rootFlow() {
  int o, d, root = graph->getRoot();
  for (auto k : graph->DuS) {
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

void Model::flowConservation() {
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
void Model::terminalsFlow() {
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

void Model::maxArcs() {
  GRBLinExpr totalArcs;
  for (int o = 0; o < graph->getN(); o++)
    for (auto *arc : graph->arcs[o])
      totalArcs += y[arc->getO()][arc->getD()];

  model.addConstr(totalArcs == (graph->getN() - 1), "maximum_of_arcs");

  model.update(); 
  cout << "maximum of arcs in the tree" << endl;
}

void Model::limDelay() {
  int i, j, paramDelay;
  for (auto k : graph->terminals) {
    GRBLinExpr limDelay;
    for (i = 0; i < graph->getN(); i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	limDelay += arc->getDelay() * f[i][j][k];
      }
    }
    paramDelay = graph->getParamDelay();
    model.addConstr(limDelay <= (paramDelay + (graph->getBigMDelay() - paramDelay) * z[k]), "delay_lim_" + to_string(k));
  }

  model.update();
  cout << "Delay limit" << endl;
}

void Model::limJitter() {
  int i, j, paramJitter;
  for (auto k : graph->terminals) {
    GRBLinExpr limJitter;
    for (i = 0; i < graph->getN(); i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	limJitter += arc->getJitter() * f[i][j][k];
      }
    }

    paramJitter = graph->getParamJitter();
    model.addConstr(limJitter <= (paramJitter + (graph->getBigMJitter() - paramJitter) * z[k]), "jitter_limit_" + to_string(k));
  }
  model.update();
  cout << "Jitter limit" << endl;
}

void Model::nonTerminalsLeafs() {
  // model.addConstr(y[graph->getRoot()][0] == 1);
  for (auto q : graph->DuS)
    for (auto e : graph->DuS) 
      if (e != q) 
	model.addConstr(f[0][q][e] == 0, "non_terminals_leafs_" + to_string(q) + "_" + to_string(e));

  model.update();
  cout << "Non terminals are leafs" << endl;
}

void Model::primeToTerminals() {
  model.addConstr(y[graph->getRoot()][0] == 1);
 
  for (auto k : graph->terminals)
    model.addConstr(z[k] >= f[0][k][k], "prime_to_terminals_" + to_string(k));

  model.update(); 
  cout << "S' to terminals" << endl;
}

bool Model::solve() {
  try {
    model.set("TimeLimit", "3600.0");
    model.set("OutputFlag", "0");
    model.update();
    model.write("model.lp");
    model.optimize();
    return true;
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
    return false;   
  }
}

void Model::getGradientDelay(vector<double> &gradientDelay) {
  int i, j, paramDelay;
  for (auto k : graph->terminals) {
    gradientDelay[k] = -graph->getParamDelay();
    if (z[k].get(GRB_DoubleAttr_X) > 0.1)
      gradientDelay[k] = -graph->getBigMDelay();

    for (i = 0; i < graph->getN(); i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	if (f[i][j][k].get(GRB_DoubleAttr_X) > 0.1)
	  gradientDelay[k] += arc->getDelay();
      }
    }
    if (gradientDelay[k] > 0) feasible = false;
  }
}

void Model::getGradientJitter(vector<double> &gradientJitter) {
  int i, j, paramJitter;
  for (auto k : graph->terminals) {
    gradientJitter[k] = -graph->getParamJitter();
    if (z[k].get(GRB_DoubleAttr_X) > 0.1)
      gradientJitter[k] = -graph->getBigMJitter();

    for (i = 0; i < graph->getN(); i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	if (f[i][j][k].get(GRB_DoubleAttr_X) > 0.1)
	  gradientJitter[k] += arc->getJitter();
      }
    }
    if (gradientJitter[k] > 0) feasible = false;
  }
}

void Model::getGradientVar(vector<vector<double>> &gradientVar) {
  int i, j, n = graph->getN(), bigMK, bigML;
  for (int k : graph->terminals) {
    for (int l : graph->terminals) {
      if (k != l) {
	bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
	bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
	gradientVar[k][l] = -graph->getParamVariation();
	
	if (z[k].get(GRB_DoubleAttr_X) > 0.1) gradientVar[k][l] -= bigMK;
	if (z[l].get(GRB_DoubleAttr_X) > 0.1) gradientVar[k][l] -= bigML;

	for (i = 0; i < n; i++) {
	  for (auto *arc : graph->arcs[i]) {
	    j = arc->getD();
	    if (f[i][j][k].get(GRB_DoubleAttr_X) > 0.1) gradientVar[k][l] += arc->getDelay();
	    if (f[i][j][l].get(GRB_DoubleAttr_X) > 0.1) gradientVar[k][l] -= arc->getDelay();
	  }
	}
	if (gradientVar[k][l] > 0) feasible = false;
      }
    }
  }
}

void Model::getGradientRelation(vector<vector<vector<double>>> &gradientRel) {
  int i, j, n = graph->getN();
  for (auto k : graph->DuS) {
    for (i = 0; i < n; i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	gradientRel[i][j][k] = 0;
	if (f[i][j][k].get(GRB_DoubleAttr_X) > 0) gradientRel[i][j][k] += 1;
	if (y[i][j].get(GRB_DoubleAttr_X) > 0) gradientRel[i][j][k] -= 1;
	if (gradientRel[i][j][k] > 0) feasible = false;
      }
    }
  }
}

void Model::getGradientLeaf(vector<vector<double>> &gradientLeaf) {
  int i, j, n = graph->getN();
  for (auto e : graph->DuS) {
    for (auto q : graph->DuS) {
      if (q != e) {
	if (f[0][q][e].get(GRB_DoubleAttr_X) > 0) gradientLeaf[q][e] = 1;
	else gradientLeaf[q][e] = 0;

	if (gradientLeaf[q][e] > 0) feasible = false;
      }
    }
  }
}

double Model::getNormDelaynJitter(vector<double> &gradient) {
  double sum = 0;
  for (auto k : graph->terminals)
    sum += pow(gradient[k], 2);
  return sqrt(sum);
}

double Model::getNormRelation(vector<vector<vector<double>>> &gradient) {
  double sum = 0;
  int i, j, n = graph->getN();
  for (auto k : graph->terminals) {
    for (i = 0; i < n; i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	sum += pow(gradient[i][j][k], 2);
      }
    }
  }
  return sqrt(sum);
}

double Model::getNormVar(vector<vector<double>> &gradient) {
  double sum = 0;
  for (int k : graph->terminals)
    for (int l : graph->terminals)
      if (k != l) sum += pow(gradient[k][l], 2);
  return sqrt(sum);
}

double Model::getNormLeaf(vector<vector<double>> &gradient) {
  double sum = 0;
  for (int e : graph->DuS)
    for (int q : graph->DuS)
      if (q != e) sum += pow(gradient[q][e], 2);
  return sqrt(sum);
}

int Model::getOriginalObjValue() {
  int foValue = 0;
  for (auto t : graph->terminals) 
    if (z[t].get(GRB_DoubleAttr_X) > 0) foValue++; 
  return foValue; 
}

bool Model::isFeasible() {
  if (feasible) return true; 
  feasible = true;
  return false;
}

int Model::lagrangean() {
  int progress = 0, iter = 0, n = graph->getN();
  double thetaDelay, normDelay, thetaJitter, normJitter, thetaRel, normRel, thetaVar, normVar, normLeaf, thetaLeaf, objPpl, originalObj;

  vector<double> gradientDelay = vector<double>(n);
  vector<double> gradientJitter = vector<double>(n);
  vector<vector<double>> gradientVar = vector<vector<double >>(n, vector<double>(n));
  vector<vector<double>> gradientLeaf = vector<vector<double >>(n, vector<double>(n));
  vector<vector<vector<double>>> gradientRel = vector<vector<vector<double >>>(n, vector<vector<double>>(n, vector<double>(n)));

  multipliersDelay = vector<double>(n);
  multipliersJitter = vector<double>(n);
  multipliersVar = vector<vector<double >>(n, vector<double>(n));
  multipliersLeaf = vector<vector<double >>(n, vector<double>(n));
  multipliersRel = vector<vector<vector<double >>>(n, vector<vector<double>>(n, vector<double>(n)));

  if (barrierM) {
    auto start = chrono::steady_clock::now();

    BarrierMethod *bm  = new BarrierMethod(graph);
    bm->initModel();
    bm->solve();    
    if (relaxNum == 1 || relaxNum == 3) bm->getMultipliersDelay(multipliersDelay);
    if (relaxNum <= 2) bm->getMultipliersJitter(multipliersJitter);

    bm->getMultipliersRelation(multipliersRel);
    bm->getMultipliersVariation(multipliersVar);
    bm->getMultipliersLeaf(multipliersLeaf);

    auto end  = chrono::steady_clock::now();
    bmTime = chrono::duration_cast<chrono::seconds>(end - start).count();
  }
  
  auto start = chrono::steady_clock::now();
  auto end  = chrono::steady_clock::now();
  endTime = 0;
  
  if (heuristics) UB = heuristic->initialHeuristic();

  LB = numeric_limits<short>::min();

  firstUB = UB;
  iterBub = 0;
  
  initialize();
  if (relaxNum == 1) initModelShp();
  else if (relaxNum == 2) initModel1ResDelayCshp();
  else if (relaxNum == 3) initModel1ResJitterCshp();
  else initModel2ResCshp();
  
  while(iter < maxIter && endTime < time) {
    if (solve()) {

      if (relaxNum == 1 || relaxNum == 3) getGradientDelay(gradientDelay);
      if (relaxNum <= 2) getGradientJitter(gradientJitter);

      getGradientVar(gradientVar);
      getGradientRelation(gradientRel);
      getGradientLeaf(gradientLeaf);

      objPpl = model.get(GRB_DoubleAttr_ObjVal);
      if (iter == 0) firstLB = objPpl;
      	
      if (objPpl > LB)
	LB = objPpl, progress = 0, iterBlb = iter;
      else {
	progress++;
	if (progress == B) {
	  lambda /= 2;
	  progress = 0;
	}
      }

      originalObj = getOriginalObjValue();
      // cout << "Original Obj: " << originalObj << endl;

      if (isFeasible() && originalObj < UB) {
	UB = originalObj, iterBub = iter;
	if ((UB - LB) / UB <= 0.0001) return UB;
      }

      // Heuristic
      if (heuristics) {
	int heuObj = heuristic->subgradientHeuristic(multipliersRel, multipliersLeaf);
	// cout << "Heuristic Obj: " << heuObj << endl;
	if (heuObj < UB) {
	  UB = heuObj, iterBub = iter;
	  if ((UB - LB) / UB <= 0.0001) return UB;
	}
      }
      
      if (relaxNum <= 2) {
	normJitter = getNormDelaynJitter(gradientJitter);
	if (normJitter == 0) thetaJitter = 0;
	else thetaJitter = lambda * ((UB - objPpl) / pow(normJitter, 2));
      }

      if (relaxNum == 1 || relaxNum == 3) {
	normDelay = getNormDelaynJitter(gradientDelay);
	if (normDelay == 0) thetaDelay = 0;
	else thetaDelay = lambda * ((UB - objPpl) / pow(normDelay, 2));
      }
      
      normVar = getNormVar(gradientVar);
      normRel = getNormRelation(gradientRel);
      normLeaf = getNormLeaf(gradientLeaf);
      
      if (normVar == 0) thetaVar = 0;
      else thetaVar = lambda * ((UB - objPpl) / pow(normVar, 2));

      if (normRel == 0) thetaRel = 0;
      else thetaRel = lambda * ((UB - objPpl) / pow(normRel, 2));

      if (normLeaf == 0) thetaLeaf = 0;
      else thetaLeaf = lambda * ((UB - objPpl) / pow(normLeaf, 2));            

      for (int k : graph->terminals) {
	if (relaxNum == 1 || relaxNum == 3) multipliersDelay[k] = max(0.0, multipliersDelay[k] + (gradientDelay[k] * thetaDelay));

	if (relaxNum <= 2)  multipliersJitter[k] = max(0.0, multipliersJitter[k] + (gradientJitter[k] * thetaJitter));

	for (int l : graph->terminals)
	  if (k != l) multipliersVar[k][l] = max(0.0, multipliersVar[k][l] + (gradientVar[k][l] * thetaVar));
      }
      
      for (auto k : graph->DuS) 
	for (int i = 0; i < n; i++) 
	  for (auto *arc : graph->arcs[i]) 
	    multipliersRel[i][arc->getD()][k] = max(0.0, multipliersRel[i][arc->getD()][k] + gradientRel[i][arc->getD()][k] * thetaRel);

      for (int e : graph->DuS)
	for (int q : graph->DuS)
	  if (e != q) multipliersLeaf[q][e] = max(0.0, multipliersLeaf[q][e] + (gradientLeaf[q][e] * thetaLeaf));

      // cout << "(Feasible) Upper Bound = " << UB << ", (Relaxed) Lower Bound = " << LB << endl;

      if (relaxNum == 1) setObjRL1();
      else if (relaxNum == 2) setObjRL2();
      else if (relaxNum == 3) setObjRL3();
      else setObjRL4();

      iter++;
      end = chrono::steady_clock::now();
      endTime = chrono::duration_cast<chrono::seconds>(end - start).count();
      //      getchar();
    }
  }

  return UB;
}

void Model::showSolution(string instance, int prepTime) {
  ofstream output;
  output.open(instance, ofstream::app);

  output << "Prep. Time: " << prepTime << endl;
  output << "First LB: " << firstLB << "\nLB: " << LB << "\nIter. LB: " << iterBlb << endl;
  output << "First UB: " << firstUB << "\nUB: " << UB << "\nIter. UB: " << iterBub << endl;
    
  if (LB < 0) LB = 0;
  
  output << "gap: " << 100 * (double(UB - ceil(LB)) / double(UB)) << endl;
  output << "BM. Time: " << bmTime << "\nRuntime: " << endTime << endl;
  
  output.close();
}
