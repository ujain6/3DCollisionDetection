/*************************************************
 * collide.cpp
 *  OpenMP program that 
 *  Detect collisions between a triangular mesh
 *  and a uniform list of spheres.
 */

#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <fcntl.h>
#include <omp.h>
#include <list>
#include <tuple>
#include <vector>
#include <unordered_map>
#include "utils/debug.h"
#include "io/csv.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "io/tiny_obj_loader.h"
#include "geometry/hash3d.h"
#include "geometry/geometry.h"
#include "utils/time_omp.h"
#include "utils/time_unix.h"

using namespace std;
using namespace tinyobj;

#pragma mark - Global Variables

/// Number of bins in total
int numBins;
float binSize;
float radius;
float lBound;

vector<tuple<int, Point>> read_sphere_csv(const char * fname){
  vector<tuple<int, Point>> result;
  io::CSVReader<3> in(fname);
  in.read_header(io::ignore_extra_column, "x", "y", "z");
  float x,y,z;

  int sphereID = 0;
  Point maxCenter(0,0,0);
  while(in.read_row(x,y,z)){
    float distFromOrigin =  float(sqrt( (double)(pow(x,2) + pow(y,2) + pow(z,2) )));
    tuple<int, Point> tup = make_tuple(sphereID++, Point(x,y,z));
    result.push_back(tup);
    if(lBound < distFromOrigin){
      lBound = std::max(lBound, distFromOrigin);
      maxCenter = get<1>(tup);
    }
  }
  // We make a cube with diagonal 2*lBound/root3, therefore the edge of cube is
  cout << "Largest center is " << maxCenter << endl;

  // Therefore, the cube should end at lBound + R,
  lBound = (2*lBound + radius);
  cout << "Creating a cube with edge = " << lBound << endl;
  return result;
}

spatialMap myMap;
//From the meshfile, create the face


/// Read Meshes and return the vector of mesh. (each mesh has a vector of <Triangle, tid> pair) 
vector<vector<tuple<int, Triangle>>> read_mesh_shapes(const char* meshfile){
  ObjReader reader;
  reader.ParseFromFile(meshfile);
  const tinyobj::attrib_t attrib = reader.GetAttrib();
  const std::vector<tinyobj::shape_t>& shapes = reader.GetShapes();

  vector<vector<tuple<int, Triangle>>> result;

  //Loop over shapes
  //debug("Read mesh shapes start.");
  for(size_t s = 0; s < shapes.size(); s++){
    size_t offset = 0;
    
    //For each shape, we have a mesh of triangles
    vector<tuple<int, Triangle>> triangle;
    
    //Going over each face
    size_t numOfface = shapes[s].mesh.num_face_vertices.size();
    int triangleID = 0;
    
    //Loop over face of a polygon
    for(size_t f = 0; f < numOfface; f++){
      
      int fv = shapes[s].mesh.num_face_vertices[f];
      
      //For this triangle, get all the 3 points
      float vx, vy, vz;
      tinyobj::index_t idx;
      vector<Point> vertices;
      for(size_t v = 0; v < fv; v++){
        // Access to vertex
        idx = shapes[s].mesh.indices[offset + v];
        vx = (float)(attrib.vertices[3*idx.vertex_index+0]);
        vy = (float)(attrib.vertices[3*idx.vertex_index+1]);
        vz = (float)(attrib.vertices[3*idx.vertex_index+2]);

        //Add this vertex to the vector that holds the triangle
        Point vert(vx,vy,vz);
        vertices.push_back(vert);
        //    debug("|--- Read in Vertex %lu: %f,%f,%f", v, vx, vy,vz);
      }
      //Create a triangle object from this Triangle
      Triangle trig(vertices[0],vertices[1],vertices[2]);
      tuple<int, Triangle> tup = make_tuple(triangleID++, trig);
      triangle.push_back(tup);

      offset += fv;
    }
    //debug("Read mesh shapes finished. Create Triangle in total: %d", triangle.size());
    //After going through each face of a mesh, add the mesh to the final file
    result.push_back(triangle);
  }
  return result;
}


/// Hash the triangle to its corresponding bin
void hashToBinsTriangle(Triangle& trig, int tID){
  
  // Map the coordinates to frame aligned at the bottom of the cube
  trig = trig.shiftIt(lBound/2);

  // Get the maximum bounds of this triangle
  Point maxBound(trig.maxBound());
  Point minBound(trig.minBound());

  // Bin number where maxBound lies
  int p = floor( maxBound.x/ binSize);
  int q = floor( maxBound.y/ binSize);
  int r = floor( maxBound.z/ binSize);

  // Bin number where minBound lies
  int pp = floor( minBound.x/ binSize);
  int qq = floor( minBound.y/ binSize);
  int rr = floor( minBound.z/ binSize);

  // For each object decide which bin it intersects
  for(int k = pp; k <= p; k++){
    for(int l = qq; l <= q; l++){
      for(int m = rr; m <= r; m++){
        myMap[Point(k,l,m)].triangles.push_back(make_tuple(tID, trig));
      }
    }
  }
}

/// Hash a sphere center to its corresponding bin.
void hashToBinsSphere(Point& t, int sID){
  t = t.shiftIt(lBound/2);

  //Bin number where maxBound lies
  int p = floor( (t.x + radius)/ binSize);
  int q = floor( (t.y + radius)/ binSize);
  int r = floor( (t.z + radius)/ binSize);

  //Bin number where minBound lies
  int pp = floor( (t.x - radius)/ binSize);
  int qq = floor( (t.y - radius)/ binSize);
  int rr = floor( (t.z - radius)/ binSize);

  for(int k = pp; k <= p; k++){
    for(int l = qq; l <= q; l++){
      for(int m = rr; m <= r; m++){
        myMap[Point(k,l,m)].spheres.push_back(make_tuple(sID, t));
      }
    }
  }
}

/// Add spheres to the their respective bins
void stage1(vector<tuple<int, Point>> spheres){
  for(size_t i = 0; i < spheres.size(); i++){
    tuple<int, Point> sphere = spheres[i];
    int sID = get<0>(sphere);
    hashToBinsSphere(get<1>(sphere), sID);
  }
}

/// Add Triangles to the their respective bins
void stage2(vector<vector<tuple<int, Triangle>>> shapes){
  //For each mesh
  for(int i = 0; i < shapes.size(); i++){
    //For each triangle of the mesh
    vector<tuple<int, Triangle>> meshOfTriangles = shapes[i];
    cout << "Mesh has " << meshOfTriangles.size() << " faces" << endl;
    for(int j = 0; j < meshOfTriangles.size(); j++){
      //For each triangle of the mesh, hash it to the bins
      tuple<int, Triangle> triangle = meshOfTriangles[j];
      int tID = get<0>(triangle);
      //Call on this Triangle object, TriangleID
      hashToBinsTriangle(get<1>(triangle), tID);
    }
  }
}

/// Calculate collisions between spheres and triangles, write to file
vector<tuple<int,int>> substage4(Point binID, vector<tuple<int, Point>> spheres, vector<tuple<int, Triangle>> triangles, int i){
  vector<tuple<int, Point>>::iterator it;
  vector<tuple<int, Triangle>>::iterator it2;
  vector<tuple<int, int>> coll;
  for (it = spheres.begin() ; it != spheres.end(); ++it){
    for (it2 = triangles.begin() ; it2 != triangles.end(); ++it2){

      const Point p1(get<1>(*it));
      const Triangle t1(get<1>(*it2));
      const num_t r = radius;
      if(Sphere::intersect(p1, r, t1)){
        //Get the midpoint of the vector b/w the centroids of the two bodies
        Point ct = t1.centroid();
        
        //Get the midpoint of the vector
        Point midPoint( (p1.x + ct.x)/2, (p1.y + ct.y)/2, (p1.z + ct.z)/2 );
        
        //Get the bin number of this
        Point collisionBin( floor(midPoint.x/binSize), floor(midPoint.y/binSize), floor(midPoint.z/binSize) );
        int sID = get<0>(*it);
        int tID = get<0>(*it2);
        
        //If midpoint of collision volume lies in my bin, count it
        if(collisionBin == binID){
          //cout << "SID:" << sID << " TID:"<< tID << "( S = Sphere(" << p1 << ", r=" << r << ", T = " << t1 << ")" << endl;
          coll.push_back(make_tuple(sID, tID)); 
        }
      }
    }
  }//For each sphere ends here
  return coll;
}

/// Write the (sid,tid) pairs to the output file.
int writeToFile(float walltime, vector<vector<tuple<int, int>>> collisions, const char* outfile){
  FILE* fp = fopen(outfile, "w+");
  assert(fp != NULL);

  // Write wall time (in float)
  fprintf(fp, "%f\n", walltime);

  // Write all collision pairs
  int count = 0;
  for(int i = 0; i < collisions.size(); i++){
    vector<tuple<int, int>> perBinData = collisions[i];
    for(int j = 0; j < perBinData.size(); j++){
      count++;
      tuple<int, int> tup = perBinData[j];
      fprintf(fp, "%d,%d\n", get<0>(tup), get<1>(tup));
    }
  }
  return count;
}

/// In parallel, detect collision.
void stage3(spatialMap myMap, const char* outfile){
  start_unix_timing(4);
  vector<Point> keys;
  //Create an auxillary data structure to hold the keys
  for (std::pair<Point, hashValue> element : myMap){
    keys.push_back(element.first);
  }

  //Create another array that holds for each bin, the corresponding
  //TID and SID values
  vector<vector<tuple<int, int>>> collisions (keys.size());
  std::cout << "Map has " << myMap.size() << " bins and vector has " << keys.size() << " keys" << endl;
  stop_unix_timing(4);
  std::cout << "Stage 3.1 done (in " << get_time_elapse(4) << " ms)" << endl;

  //Each thread handles a specific bin, applies brute force collison checks
  //Fills the result in a corresponding data structure
  omp_set_num_threads(omp_get_max_threads());
  
  start_omp_timing(1);
  #pragma omp parallel for shared(myMap)
  for(int i = 0; i < keys.size(); i++){
    if(i == 0){ std::cout << "Start Collision Detection with " << omp_get_num_threads() << " threads." << endl; }
    Point bin = keys[i];
    hashValue hv = myMap[bin];
    vector<tuple<int, Point>> spheres = hv.spheres;
    vector<tuple<int, Triangle>> triangles = hv.triangles;
    vector<tuple<int,int>> coll = substage4(bin, spheres, triangles, i);
    //Add this vector to the particular index_t
    collisions[i] = coll;
  }
  stop_omp_timing(1);

  const float walltime = get_time_elapse(1);
  std::cout << "Finished Collision Detection (in " << walltime << " ms)" << endl;

  int count = writeToFile(walltime, collisions, outfile);
  printf("Collision %d in total\n", count);

}

int main(int argc, const char* argv[]) {

  usage<5>(argc, "./collide [meshfile] [spherefile] [radius] [outfile]");

  const char * meshfile   = argv[1]; // meshfile: is the name of the input file that contains the mesh stored in the Wavefront OBJ format
  const char * spherefile = argv[2]; // spherefile: is the name of the CSV input file for all spheres.
  radius                  = atof(argv[3]); // radius: the radius of every sphere
  const char * outfile    = argv[4]; // outfile: the name of the output file
  binSize                 = 2 * radius;

  // Read Meshfile .obj (all the faces from meshfile)
  vector<vector<tuple<int, Triangle>>> result = read_mesh_shapes(meshfile);

  // Read spherefile
  vector<tuple<int, Point>> spheres = read_sphere_csv(spherefile);
  cout << "Number of spheres created = " << spheres.size() << endl;

  // Setup our simulation namespace
  numBins = floor(pow(double(lBound/binSize), 3));

  // Calculate array T, where object i collides with T[i] bins
  start_unix_timing(0);
  stage1(spheres);
  stop_unix_timing(0);
  cout << "Stage 1 done (in " << get_time_elapse(0) << " ms)" << endl;
  

  start_unix_timing(2);
  stage2(result);
  stop_unix_timing(2);
  cout << "Stage 2 done (in " << get_time_elapse(2) << " ms)" << endl;
  
  stage3(myMap, outfile);
  cout << "Stage 3 done " << endl;
  cout << "Program Ended" << endl;


  return 0;
}
