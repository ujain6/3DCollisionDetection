
#include <iostream>
#include <ostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <list>

#include "../geometry/geometry.h"
#include "../geometry/hash3d.h"
using namespace std;


#pragma mark - Testing Functions
void test_point(){
    // Initialize Point
    Point a(1,2,3);
    cout << a << endl;
    cout << a.x << a.y << a.z << endl;

    // Point Addition
    cout << Point(1,2,3) + Point(4,5,6) << endl;

    // Dot, Cross product
    cout << Point(1,2,3).dot(Point(1,2,3)) << endl;
    cout << Point(1,1,1).dot(Point(-1,-1,2)) << endl;

    cout << Point(1,2,3).cross(Point(1,2,3)) << endl;
    cout << Point(1,2,3).cross(Point(0,0,1)) << endl;

    // Unit vector
    cout << Point(1,2,3).unit() << endl;

    // Mixed vector
    cout << Point(1,2,3).mixed(Point(1,2,3), Point(1,2,3)) << endl;
    cout << Point(1,0,0).mixed(Point(0,0,1), Point(1,1,0)) << endl;
}


void test_segment(){
    // Initialize segment with right reference of points
    cout << Segment(Point(1,2,3), Point(4,5,6)) << endl;

    Point a(0,0,0);
    Point b(0,2,0);
    Segment s(a,b);
    cout << s << endl;

    // Test segment unit vector in right reference objects
    cout << Segment(Point(1,2,3), Point(4,5,6)).unit() << endl;

    // Test projection of point on the segment (ray)
    cout << Segment(Point(1,2,3), Point(2,3,4)).proj(Point(1,2,3)) << endl;
    cout << Segment(Point(1,2,3), Point(2,3,4)).proj(Point(2,3,4)) << endl;
    cout << Segment(Point(1,2,3), Point(2,3,4)).proj(Point(3,4,5)) << endl;


}


void test_triangle(){
    // Init
    cout << Triangle(Point(0,0,0),Point(0,2,0),Point(0,0,2)) << endl;
    cout << Triangle(Point(0,0,0),Point(0,2,0),Point(0,0,2)).a << endl;
    cout << Triangle(Point(0,0,0),Point(0,2,0),Point(0,0,2)).b << endl;
    cout << Triangle(Point(0,0,0),Point(0,2,0),Point(0,0,2)).c << endl;
    
    // Edge
    cout <<Triangle(Point(0,0,0),Point(0,2,0),Point(0,0,2)).edge(0) << endl;
    cout <<Triangle(Point(0,0,0),Point(0,2,0),Point(0,0,2)).edge(1) << endl;
    cout <<Triangle(Point(0,0,0),Point(0,2,0),Point(0,0,2)).edge(2) << endl;

    // Vector
    cout <<Triangle(Point(0,0,0),Point(0,2,0),Point(0,0,2)).vec(0) << endl;
    cout <<Triangle(Point(0,0,0),Point(0,2,0),Point(0,0,2)).vec(1) << endl;
    cout <<Triangle(Point(0,0,0),Point(0,2,0),Point(0,0,2)).vec(2) << endl;

    // Normal Vector
    cout <<Triangle(Point(0,0,0),Point(0,2,0),Point(0,0,2)).normal() << endl;

    // In Triangle Test
    Triangle trig(Point(0,0,0),Point(0,2,0),Point(0,0,2));
    assert(trig.inTriangle(Point(0,0,0)) == false);
    assert(trig.inTriangle(Point(0,2,0)) == false);
    assert(trig.inTriangle(Point(0,0,2)) == false);
    assert(trig.inTriangle(Point(0,0.5,0.5)) == true);

    // dist and directed dist
    cout << trig.directed_vdist(Point(2,2,2)) << endl;

    // nearest distance
    cout << trig.nearest_dist(Point(0,2,2)) << endl;
    assert(approx(trig.nearest_dist(Point(0,2,2)), sqrt(2)));
}


void test_sphere(){
    // Sphere s; // This is forbidden
    Point center(0,0,0);
    num_t radius = 1;
    Triangle trig(Point(0,0,0), Point(0,1,0),Point(0,0,1));
    assert(Sphere::intersect(center, radius, trig) == true);
}

#pragma mark - Main Function

int main(int argc, const char * argv[]) {
    test_point();
    test_segment();
    test_triangle();
    test_sphere();
    return 0;
}
