#include <system/SFSystem.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  if (argc != 5) {
    cerr << "Usage: optimize-particle-system <input.txt filename> <max # iterations> <first interface> <last interface>\n";
    exit(0);
  }
  SFSystem particle_sys(argv[1]);
  int num_iters = atoi(argv[2]);
  int min_inter = atoi(argv[3]);
  int max_inter = atoi(argv[4]);
  int num_intersections = particle_sys.totalNumberOfSystems();
  
  for (int intersection=0; intersection<num_intersections; intersection++) 
  {
    if (intersection >= min_inter && intersection <= max_inter) 
    {
      unsigned int iter=0;
      bool optimized = false;
      while (iter<num_iters && !optimized) 
      {  
			particle_sys.optimize();
			iter++;
			
			cout << "optimize-particle-system [" << iter << "|" << num_iters << "] iteration." << std::endl;
			
			if (particle_sys.optimized()) 
			{ 
			  particle_sys.resetLambdas();
			  optimized = true;
			}
      }
      cout << "optimize-particle-system finished after [" << iter << "|" << num_iters << "] iterations." << std::endl;
      particle_sys.freezeIntersection();
      particle_sys.writePointFile(intersection);
      //cout << "DEBUG write to file" << endl;
    } 
    else
    {
      particle_sys.freezeIntersection();
    }
  }
  return 0;
}
