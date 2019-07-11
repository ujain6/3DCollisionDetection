#ifndef _HASH3D_H_
#define _HASH3D_H_

#include <iostream>
#include <ostream>
#include <unordered_map>
#include <unistd.h>
#include <iterator>

#include "geometry.h"
using namespace std;

/* Hash Func: we used the hash function from stackoverflow:
 * https://stackoverflow.com/questions/16792751/hashmap-for-2d3d-coordinates-i-e-vector-of-doubles
 */

struct hashFunc {
  size_t operator()(const Point &k) const{
    size_t h1 = std::hash<double>()(k.x);
    size_t h2 = std::hash<double>()(k.y);
    size_t h3 = std::hash<double>()(k.z);
    return (h1 ^ (h2 << 1)) ^ (h3 << 2);
  }
};

struct equalsFunc{
  bool operator()( const Point & lhs, const Point & rhs ) const{
    return (lhs.x == rhs.x) && (lhs.y == rhs.y) && (lhs.z == rhs.z);
  }
};


//Create an objects
class hashValue{
  public:
    vector<tuple<int, Point>> spheres;
    vector<tuple<int, Triangle>> triangles;

    hashValue(){}
};

typedef std::unordered_map<Point, hashValue, hashFunc, equalsFunc> spatialMap;

#endif
