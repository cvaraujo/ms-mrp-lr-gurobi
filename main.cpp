#include <sys/stat.h>
#include "headers/Graph.h"
#include "headers/BarrierMethod.h"
#include "headers/Model.h"

int main(int argc, const char *argv[]) {
  if (argc < 4) {
    return 0;
  } else {
    // ./MS relax graph param outp mve sae heur BM lambda maxIter B Time
    // mkdir("results", 0777);
    auto *graph = new Graph(argv[2], argv[3], argv[4]);
    stringstream rn(argv[1]), mve(argv[5]), sae(argv[6]), heu(argv[7]),
      bm(argv[8]), lambdaToken(argv[9]), mi(argv[10]), bToken(argv[11]), timeToken(argv[12]);
    
    int maxIter, b, time, relaxNum;
    double lambda;
    bool prepMve, prepSae, heuristics, barrierM;

    rn >> relaxNum; mve >> prepMve; sae >> prepSae; heu >> heuristics;
    bm >> barrierM; lambdaToken >> lambda; mi >> maxIter; bToken >> b;
    timeToken >> time;

    auto start = chrono::steady_clock::now();
    
    if (prepMve || prepSae) {
      if (prepMve) graph->MVE(argv[4]);
      if (prepSae) graph->SAE(argv[4]);
      graph->finishPreprocessing(argv[4], prepMve, prepSae);
    }
    auto end = chrono::steady_clock::now();
    int prepTime = chrono::duration_cast<chrono::seconds>(end - start).count();
    
    Model *model = new Model(graph, relaxNum, heuristics, barrierM, lambda, maxIter, b, time);
    model->lagrangean();
    model->showSolution(argv[4], prepTime);
  }
  return 0;
}
