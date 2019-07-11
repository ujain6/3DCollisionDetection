#include <iostream>
#include <tuple>
#include <unistd.h>
#include <cstdlib>
#include <fcntl.h>
#include <vector>

#include "../io/csv.h"
using namespace std;


const char * path = "../examples/spherefile/sample_spheres.csv";


int main(){
  io::CSVReader<3> in(path);
  in.read_header(io::ignore_extra_column, "x", "y", "z");
  float x; float y; float z;
  int ID = 0;
  while(in.read_row(x, y, z)){
    if(ID < 10000){
      cout << x << "," << y << "," << z << endl;
    }
    ID++;
  }
}
