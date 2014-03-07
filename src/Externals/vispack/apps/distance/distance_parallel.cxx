#include <volume.h> 
#include <volumefile.h> 
#include <volutils.h>
#include <octree.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <sys/types.h>
#include <sys/wait.h>
//#include <visfield.h>

using namespace std;

pid_t spawn(const char *progname, char *argv[]){
  cout << "Running: ";
  for(int i=0; argv[i]; i++)
    cout << argv[i] << ' ';
  cout <<  endl;

  pid_t child_pid = fork();
  if(child_pid == 0){
    // I am the child.
    execv(progname, argv);
  }
  else{
    // I am parent.
    return child_pid;
  }
}

int main(int argc, char** argv)
{
  VISVolumeFile volume_file;
  VISVolume<float> vol;

  // to restrict compution to a bbox, e.g. for debugging...
  int w, h, d, z_lo, z_hi;

  
  int i, j, k, l;
  if (argc < 5)
  {
    cout << "usage: distance_parallel inputfile material_number epsilon volume_out pts_out num_processes" << endl;
    exit(-1);
  }

  int args = 1;
  string infile(argv[args++]);
  cout << infile << endl;
  
  string material(argv[args++]);
  cout << material << endl;

  string epsilon(argv[args++]);
  cout << epsilon << endl;

  string outname_vol(argv[args++]);

  string outname_pts(argv[args++]);

  string num_procs(argv[args++]);
  int num_proc = atoi(num_procs.c_str());
  cout << "nprocs: " << num_proc << endl;

  ifstream s_infile(infile.c_str());
  float sx, sy, sz;
  int num_files;
  s_infile >>   num_files;
  string inname_vol;
  s_infile >> sx >> sy >> sz >> inname_vol;
  cout << inname_vol << endl;
  s_infile.close();
  vol = VISVolume<float>(volume_file.read(inname_vol.c_str()));
  w = vol.width(); h = vol.height(); d = vol.depth(); 



  int slice_min, slice_max, slices_per_proc = 0;
  if (num_proc > 1) {
    slice_min = slice_max = slices_per_proc = d/(num_proc-1);
  } else { 
    slice_min = slice_max = slices_per_proc = 1;
  }
    

  slice_min = 0;

  cout << w << " " << h << " " << d << endl;

  

  // SILLY ROSS CODE
  //   for (i = 0; i < num_proc; i++)
  //     {
  //       sprintf(outfile_vol, "%s-%d-parallel", outname_vol, i);
  //       sprintf(outfile_pts, "%s-%d-parallel", outname_pts, i);
  //       slice_max = VISmin(slice_min + slices_per_proc + (i < (d%num_proc)) + 1, d);
  //       sprintf(path, "./");
  //       sprintf(process_call, "rm distance_parallel.out.%d", i);
  //       sprintf(process_call, "%sdistance3 %s %d %f %s %s %d %d %d %d %d %d > distance_parallel.out.%d &", 
  // 	      path, infile, material, epsilon, outfile_vol, outfile_pts, 0, w, 0, h, slice_min, slice_max,i);
  //       cout << "About to call process: " << process_call << endl;
  //       system(process_call);
  //       slice_min = slice_max - 2;
  //     }

  // SMART MIRIAH/RONI CODE
  char *my_argv[13];
  for(int i=0; i<12; i++)
    my_argv[i] = new char[80];
      
  for (i = 0; i < num_proc; i++)
    {
      //slice_max = VISmin(slice_min + slices_per_proc + 4, d);
      if ( i != (num_proc-1) )
        slice_max = VISmin(slice_min + slices_per_proc + 2, d);
      else
        slice_max = d;

      sprintf(my_argv[0], "./distance3");
      sprintf(my_argv[1], "%s", infile.c_str());
      sprintf(my_argv[2], "%d", atoi(material.c_str()));
      sprintf(my_argv[3], "%f", epsilon.c_str());
      sprintf(my_argv[4], "%s-%d-parallel", outname_vol.c_str(), i);
      sprintf(my_argv[5], "%s-%d-parallel", outname_pts.c_str(), i);
      sprintf(my_argv[6], "0");
      sprintf(my_argv[7], "%d", w);
      sprintf(my_argv[8], "0");
      sprintf(my_argv[9], "%d", h);
      sprintf(my_argv[10], "%d", slice_min);
      sprintf(my_argv[11], "%d", slice_max);
      my_argv[12] = NULL;

      cout << "Spawning process " << i << endl;
      spawn(my_argv[0], &my_argv[0]);

      slice_min = slice_max - 2;
    }

  for(int i=0; i<12; i++)
    delete my_argv[i];

  for(int i=0; i<num_proc; i++){
    int status;
    wait(&status);
    if(WIFEXITED(status))
      cout << "A child exited = " << i << endl;
  }
}
