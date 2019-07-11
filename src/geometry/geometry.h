#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include <iostream>
#include <ostream>
#include <cassert>
#include <vector>
#include <cmath>

using namespace std;
typedef float num_t;
#pragma mark - Math

/// Return if the two value are approximately equal.
template <typename T, typename U> bool approx(T a, U b, const float error = 1e-3) {
    return fabs(a - b) < 1e-3;
}

#pragma mark - Geometry

// TODO: The reference might be a little bit space costly

/// Represent a point or a vector in 3D.
class Point{
public:
    num_t x, y, z;

    Point(){}

    Point(num_t x, num_t y, num_t z){
        this->x = x; this->y = y; this->z = z;
    }

    Point(const Point & p):x(p.x), y(p.y), z(p.z){}

    friend Point operator + (Point lhs, const Point & rhs){
        return Point(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
    }

    friend Point operator - (Point lhs, const Point & rhs){
        return Point(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
    }

    friend Point operator - (Point lhs){
        return Point(-lhs.x, -lhs.y, -lhs.z);
    }

    friend inline bool operator == (Point lhs, const Point & rhs){
        return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
    }

    /// Return cross dot of the Point
    num_t dot(const Point & o) const {
        return x * o.x + y * o.y + z * o.z;
    }

    num_t dot(const Point && o) const {
        return x * o.x + y * o.y + z * o.z;
    }
    
    //Shift a coordinate when origin shifts
    Point shiftIt(float shift) const {
    	return Point(this->x + shift, this->y + shift, this->z + shift);
    }

    /// Return cross product of the Point
    Point cross(const Point & o) const {
        return Point
        (y * o.z - z * o.y,
         z * o.x - x * o.z,
         x * o.y - y * o.x
         );
    }

    Point cross(const Point && o) const {
        return Point
        (y * o.z - z * o.y,
         z * o.x - x * o.z,
         x * o.y - y * o.x
         );
    }

    /// Return mixed product
    // TODO: The Point(c) calls for an constructor
    num_t mixed(const Point && b, const Point && c) const {
        return dot(b.cross(c));
    }

    num_t mixed(const Point & b, const Point & c) const {
        return dot(b.cross(c));
    }

    num_t mixed(const Point && b, const Point & c) const {
        return dot(b.cross(c));
    }

    num_t mixed(const Point & b, const Point && c) const {
        return dot(b.cross(c));
    }

    /// Return norm2 of the vector
    num_t norm2() const {
        return x * x + y * y + z * z;
    }

    /// Return the norm of the vector
    num_t norm() const {
        return sqrt(norm2());
    }

    /// Return the unit vector of the vector
    Point unit() const {
        float dis = x*x + y*y + z*z;
        if(dis == 0){ return *this; }
        dis = sqrt(dis);
        return Point(x/dis, y/dis, z/dis);
    }

    friend std::ostream& operator<< (std::ostream& os, const Point& p);

};

Point operator * (const Point & lhs, num_t rhs){
    return Point(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs);
}

Point operator * (num_t rhs, const Point & lhs){
    return Point(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs);
}


class Segment{
public:
    Point a, b;

    Segment(){}

    Segment(const Point & a, const Point & b): a(a), b(b){}

    Segment(Point && a, Point && b): a(a), b(b){}

    /// Return the unit vector of the segment
    Point unit() const {
        return (this->b - this->a).unit();
    }

    /// Return the projected point of point o on the Segment
    Point proj(const Point && o) const {
        Point n = unit();
        num_t k = n.dot((o - a));
        return a + k * n;
    }

    Point proj(const Point & o) const {
        Point n = unit();
        num_t k = n.dot(Point(o - a));
        return a + k * n;
    }

    /// Return true if the Segment contains the point
    bool contains(const Point & k) const {
        return approx(Point(k - a).cross(k - b).norm(), 0) &&
        (Point(k - a).dot(k - b)) < 0;
    }

    bool contains(const Point && k) const {
        return approx(Point(k - a).cross(k - b).norm(), 0) &&
        (Point(k - a).dot(k - b)) < 0;
    }

    /// Return true if the Segment contains the projection of the point
    bool contains_proj(const Point && k) const {
        return contains(proj(Point(k)));
    }

    bool contains_proj(const Point & k) const {
        return contains(proj(Point(k)));
    }
};


class Triangle{
public:
    Point a,b,c;
    Triangle(){}

    Triangle(const Point & a, const Point & b, const Point & c):a(a), b(b), c(c){}

    Triangle(Point && a, Point && b, Point && c):a(a), b(b), c(c){}

    /// Return the Segment of the edge of the Triangle
    Segment edge(int i) const {
        switch (i) {
            case 0:return Segment(Point(a), Point(b));
            case 1:return Segment(Point(b), Point(c));
            case 2:return Segment(Point(c), Point(a));
            default:assert(false);
        }
        assert(false);
        return Segment();
    }
    Triangle shiftIt(float shift) const{
      //Create a triangle from the shifted coordinates
      return Triangle(this->a.shiftIt(shift), this->b.shiftIt(shift), this->c.shiftIt(shift)); 

    }
    
    Point centroid() const{
      return Point((this->a.x+this->b.x+this->c.x)/3, (this->a.y+this->b.y+this->c.y)/3,(this->a.z+this->b.z+this->c.z)/3);

    } 
   
    //Return the maxBound of the triangle
    Point maxBound() const{
      float max_x,max_y,max_z;

      max_x = max(max(a.x, b.x), c.x);
      max_y = max(max(a.y, b.y), c.y);
      max_z = max(max(a.z, b.z), c.z);

      return Point(max_x,max_y,max_z);
    }

    Point minBound() const{
      float min_x,min_y,min_z;
      min_x = min(min(a.x, b.x), c.x);
      min_y = min(min(a.y, b.y), c.y);
      min_z = min(min(a.z, b.z), c.z);

      return Point(min_x,min_y,min_z);
    }

    /// Return the vector of the edge of the Triangle
    Point vec(int i) const {
        switch (i) {
            case 0:return b - a;
            case 1:return c - b;
            case 2:return a - c;
            default:assert(false);
        }
        assert(false);
        return Point();
    }

    /// Return the normal vector of the Triangle
    Point normal() const {
        return vec(0).cross(vec(1));
    }

    /// Return true if point p is within the triangle (excluding the boundary)
    bool inTriangle(const Point && p) const {
        Point n = normal();
        return  n.mixed(a - p, vec(0)) > 0 &&
        n.mixed(b - p, vec(1)) > 0 &&
        n.mixed(c - p, vec(2)) > 0;
    }

    bool inTriangle(const Point & p) const {
        Point n = normal();
        return  n.mixed(a - p, vec(0)) > 0 &&
        n.mixed(b - p, vec(1)) > 0 &&
        n.mixed(c - p, vec(2)) > 0;
    }

    /// Return the **directinoal verticle distance** between a point p to the triangle plane
    float directed_vdist(const Point && p) const {
        return (p - a).dot(normal().unit());
    }

    float directed_vdist(const Point & p) const {
        return (p - a).dot(normal().unit());
    }

    /// Return the **verticle distance** between a point p to the triangle plane
    float vdist(const Point && p) const {
        return directed_vdist(Point(p));
    }

    float vdist(const Point & p) const {
        return directed_vdist(Point(p));
    }

    /// Return the projected point of point p on the triangle plane
    /// (the point is on the `plane`, not necessarily within the triangle)
    Point proj(const Point && p) const {
        float k = directed_vdist(Point(p));
        const Point n = normal().unit();
        return p - k * n;
    }

    Point proj(const Point & p) const {
        float k = directed_vdist(Point(p));
        const Point n = normal().unit();
        return p - k * n;
    }

    /// Return the nearest distance of a point p to the triangle.
    /// The projected point is actually on the triangle (comapred to `proj`, which could return a point that is on the plane but not in the triangle)
    float nearest_dist(const Point & p) const {
        if (inTriangle(p)) { return 0; }
        float mindist = min(min((p-a).norm(), (p-b).norm()), (p-c).norm());
        for (int i = 0; i < 3; i++) {
            Segment s = edge(i);
            if(s.contains_proj((p))){
                mindist = min(mindist, (s.proj((p)) - p).norm());
            }
        }
        return mindist;
    }
};


class Sphere{
private:
    /// Warning: class Sphere should not actually instantialize the Sphere object.
    /// User should use the static function
    ///     Sphere::intersect(const Point & c, const num_t r, const Triangle & trig)
    /// to judge whether the sphere with (center c, radius r) is intersect with the triangle.
    Sphere(){}
    //    Point c;
    //    num_t r;

public:
    static bool intersect(const Point & c, const num_t r, const Triangle & trig){
        const Point p = trig.proj(c);
        const float perp_norm2 = (p - c).norm2();
        const float perp_norm = sqrt(perp_norm2);
        const float proj_r = sqrt(r * r - perp_norm2);
        
        if(perp_norm > r){
            return false;
        }
        return trig.nearest_dist(p) <= proj_r + 1e-3;
    }
};

#pragma mark - Logging Formatter

std::ostream& operator<< (std::ostream& os, const Point& p){
    os << "Point(" << p.x << "," << p.y << "," << p.z << ")";
    return os;
};


std::ostream& operator<< (std::ostream& os, const Segment& o){
    os << "Segment(" << o.a << "," << o.b << ")";
    return os;
};


std::ostream& operator<< (std::ostream& os, const Triangle& t){
    os << "Triangle(" << t.a << "," << t.b << "," << t.c << ")";
    return os;
};

#pragma mark - Semantically meaningful Data Structures

// TODO: These should all be upper case...

/// Tuple that maps sphereID -> coordinates 
typedef vector<tuple<int, Point>> SphereVector;

/// Tuple that maps triagleID -> coordinates 
typedef tuple<int, Triangle> face_tuple;

/// Vector contains all face_tuples (triangleID -> Coordinates)
typedef vector<face_tuple> mesh;

/// Vector that that holds mesh of each shape
typedef vector<mesh> allShapes;



#endif /* _GEOMETRY_H_ */
