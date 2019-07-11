#include <iostream>
#include <unistd.h>
#include <omp.h>
#include <cstdlib>
#include <fcntl.h>
#include <vector>

#define TINYOBJLOADER_IMPLEMENTATION
#include "io/tiny_obj_loader.h"
#include "io/csv.h"
#include "geometry/geometry.h"
#include "utils/debug.h"
#include "utils/time_unix.h"

using namespace std;
using namespace tinyobj;

/// Paths used for the tests
const char * simple_1_csv_path      = "../examples/spherefile/simple_1.csv";
const char * sample_sphere_csv_path = "../examples/spherefile/sample_spheres.csv";
const char * simple_mesh_obj_path   = "../examples/meshfile/triangles.obj";
const char * sample_mesh_obj_path   = "../examples/meshfile/sample_mesh.obj";
const char * sample_result_path     = "../examples/outfile/sample_result_pairs_only.csv";


vector<tuple<int, Point>> read_sphere_csv(const char * fname){
    vector<tuple<int, Point>> result;
    printf("Read csv file:s %s\n", fname);
    
    io::CSVReader<3> in(fname);
    in.read_header(io::ignore_extra_column, "x", "y", "z");
    
    float x,y,z; 
    int sphereID = 0; 
    float maxSID = -1;
    Point maxCenter(0,0,0);
    while(in.read_row(x, y, z)){
        tuple<int, Point> tup = make_tuple(sphereID++, Point(x,y,z));
        result.push_back(tup);
    }
    return result;
}


void test_read_sphere_csv(){
    vector<tuple<int, Point>> result;
    int id; Point p;

    // Test small triangle file
    start_unix_timing(1);
    result = read_sphere_csv(simple_1_csv_path);
    stop_unix_timing(1);
    for(auto tup: result){
        tie(id, p) = tup;
        cout << id << ", " << p << endl;
    }
    assert(result.size() == 4);
    printf("Finish reading the csv file: %d elements, %.2f (ms)\n", result.size(), get_time_elapse(1));
    
    // Test sample spheres with 1703243 elements
    start_unix_timing(1);
    result = read_sphere_csv(sample_sphere_csv_path);
    stop_unix_timing(1);
    assert(result.size() == 1703244);
    printf("Finish reading the csv file: %d elements, %.2f (ms)\n", result.size(), get_time_elapse(1));
}


vector<vector<tuple<int, Triangle>>> read_mesh_shapes(const char* meshfile){
    ObjReader reader;
    reader.ParseFromFile(meshfile);
    const tinyobj::attrib_t attrib = reader.GetAttrib();
    const std::vector<tinyobj::shape_t> & shapes = reader.GetShapes();

    vector<vector<tuple<int, Triangle>>> result;

    //Loop over shapes
    debug("Read mesh shapes start.");
    for(size_t s = 0; s < shapes.size(); s++){
        
        //Loop over face of a polygon
        size_t offset = 0;
        
        // Global Triangle ID.
        int triangleID = 0;
        
        //For each shape, we have a bunch of (id, Triangle) pairs
        vector<tuple<int, Triangle>> triangle;
        
        //Going over each face
        const size_t numOfface = shapes[s].mesh.num_face_vertices.size();

        for(size_t f = 0; f < numOfface; f++){
            
            int fv = shapes[s].mesh.num_face_vertices[f];
            
            tinyobj::index_t idx;
            
            float vx, vy, vz;

            //For this triangle, get all the 3 points
            vector<Point> vertices;
            for(size_t v = 0; v < fv; v++){
                // Access to vertex
                idx = shapes[s].mesh.indices[offset + v];
                vx = (float)(attrib.vertices[3*idx.vertex_index+0]);
                vy = (float)(attrib.vertices[3*idx.vertex_index+1]);
                vz = (float)(attrib.vertices[3*idx.vertex_index+2]);

                // Add this vertex to the vector that holds the triangle
                Point vert(vx,vy,vz);
                vertices.push_back(vert);
            }
            //Create a triangle object from this Triangle
            Triangle trig(vertices[0],vertices[1],vertices[2]);
            triangle.push_back(make_tuple(triangleID++, trig));

            offset += fv;
        }
        cout << "Number of triangles created " << triangleID << endl;

        //After going through each face of a mesh, add the mesh to the final file
        result.push_back(triangle);
    }

    debug("Read mesh shapes finished.");
    return result;
}

void test_read_mesh_shapes(){

    vector<vector<tuple<int, Triangle>>> result;
    int id; Triangle trig;
    
    // start_unix_timing(1);
    // result = read_mesh_shapes(simple_mesh_obj_path);
    // stop_unix_timing(1);
    // assert(result.size() == 1);
    // assert(result[0].size() == 3);
    // for(auto tup: result[0]){
    //     tie(id, trig) = tup;
    //     cout << id << ", " << trig << endl;
    // }
    // printf("Finish reading the mesh file: %d elements, %.2f (ms)\n", result[0].size(), get_time_elapse(1));

    start_unix_timing(1);
    result = read_mesh_shapes(sample_mesh_obj_path);
    stop_unix_timing(1);
    assert(result.size() == 1);
    assert(result[0].size() == 62976);
    printf("Finish reading the mesh file: %d elements, %.2f (ms)\n", result[0].size(), get_time_elapse(1));
}


vector<tuple<int, int>> read_result_csv(const char * fname){
    vector<tuple<int, int>> result;
    io::CSVReader<2> in(fname);
    
    int sid, tid;
    while(in.read_row(sid, tid)){
        tuple<int, int> tup = make_tuple(sid, tid);
        result.push_back(tup);
    }
    return result;
}



void test_read_result_csv(){
    vector<tuple<int, int>> result;
    start_unix_timing(1);
    result = read_result_csv(sample_result_path);
    stop_unix_timing(1);
    assert(result.size() == 39709);

    int sid, tid;
    for(int i = 0; i < 1000; i++){
        tie(sid, tid) = result[i];
        printf("%d,%d,%d\n", i, sid, tid);
    }
    printf("Finish reading the result file: %d elements, %.2f (ms)\n", result.size(), get_time_elapse(1));
}

/// (Sequential) Match all positive results in the correct answer.
/// If an answer is written on the provided example result,
/// we have to find the match and verify it.
void test_collision_positives(){
    vector<tuple<int, Point>> spheres;
    spheres = read_sphere_csv(sample_sphere_csv_path);

    vector<vector<tuple<int, Triangle>>> meshes;
    meshes = read_mesh_shapes(sample_mesh_obj_path);

    vector<tuple<int, int>> sid_tid_pairs;
    sid_tid_pairs = read_result_csv(sample_result_path);
    
    int sid, tid; const float radius = 3; int false_judges = 0;
    for(int i = 0; i < 39709; i++){
        tie(sid, tid) = sid_tid_pairs[i];
        int v_sid = get<0>(spheres[sid]);
        int v_tid = get<0>(meshes[0][tid]);
        
        assert(v_sid == sid);
        assert(v_tid == tid);

        const Triangle trig = get<1>(meshes[0][tid]);
        const Point p = get<1>(spheres[sid]);
        
        if(Sphere::intersect(p, radius, trig) == false){
            false_judges++;
            printf("%d, nominal=(%d,%d), real=(%d,%d)\n", i, sid, tid, v_sid, v_tid);
            // printf("-c Sphere::intersect(Point(%f,%f,%f), %f, Triangle( Point(%f,%f,%f), Point(%f,%f,%f), Point(%f,%f,%f) ) )\n", 
            //     p.x, p.y, p.z, radius, 
            //     trig.a.x, trig.a.y, trig.a.z,
            //     trig.b.x, trig.b.y, trig.b.z,
            //     trig.c.x, trig.c.y, trig.c.z
            // );
            printf("-c Point p(%f,%f,%f); float radius = %f; Triangle trig(Point(%f,%f,%f), Point(%f,%f,%f), Point(%f,%f,%f)); Sphere::intersect(p, radius, trig);\n", 
                p.x, p.y, p.z, radius, 
                trig.a.x, trig.a.y, trig.a.z,
                trig.b.x, trig.b.y, trig.b.z,
                trig.c.x, trig.c.y, trig.c.z
            );
        }
    }
    printf("False judged in total: %d\n", false_judges);

}

#include <unordered_map>

size_t thash(size_t sid, size_t tid, size_t msize){
    return msize * sid + tid;
}

/// (Sequential) Match all results.
void test_collision_all(){
    vector<tuple<int, Point>> spheres;
    spheres = read_sphere_csv(sample_sphere_csv_path);
    const size_t ssize = spheres.size();
    assert(ssize == 1703244);

    vector<vector<tuple<int, Triangle>>> meshes;
    meshes = read_mesh_shapes(sample_mesh_obj_path);
    const size_t msize = meshes[0].size();
    assert(msize == 62976);

    vector<tuple<int, int>> sid_tid_pairs;
    sid_tid_pairs = read_result_csv(sample_result_path);
    
    // vector<bool> answer;
    bool * answer = (bool *) calloc(ssize * msize, sizeof(bool));
    cout << ssize * msize << endl;
    assert(answer[0] == false);
    assert(answer[rand()] == false);

    int sid, tid; const float radius = 3; int false_pos = 0, false_neg = 0;

    for(int i = 0; i < sid_tid_pairs.size(); i++){
        tie(sid, tid) = sid_tid_pairs[i];
        size_t hval = thash(sid, tid, msize);
        // printf("[%d] %d %d %lld\n", i, sid, tid, hval);
        answer[hval] = true;
        assert(answer[hval] == true);
    }

    omp_set_num_threads(128);
    cout << omp_get_num_threads() << endl;
    start_unix_timing(1);
    for(int j = 0; j < msize; j++){ // mesh
        for(int i = 0; i < ssize; i++){ // sphere
            // if(i == 0){ printf(">%d\t", j);}
            tid = get<0>(meshes[0][j]);
            sid = get<0>(spheres[i]);
            const Triangle & trig = get<1>(meshes[0][j]);
            const Point & p = get<1>(spheres[i]);
            bool flag = Sphere::intersect(p, radius, trig);
            if(answer[thash(sid, tid, msize)] != flag){
                if (flag == false){
                    printf("/*c%s, %d %d*/Point p(%f,%f,%f);\tfloat radius = %f;\tTriangle trig(Point(%f,%f,%f), Point(%f,%f,%f), Point(%f,%f,%f));\n", 
                        flag ? "-" : "+", sid, tid,
                        p.x, p.y, p.z, radius, 
                        trig.a.x, trig.a.y, trig.a.z,
                        trig.b.x, trig.b.y, trig.b.z,
                        trig.c.x, trig.c.y, trig.c.z
                    );
                    false_neg++;

                }else{
                    printf("/*c%s, %d %d*/Point p(%f,%f,%f);\tfloat radius = %f;\tTriangle trig(Point(%f,%f,%f), Point(%f,%f,%f), Point(%f,%f,%f));\n", 
                        flag ? "-" : "+", sid, tid,
                        p.x, p.y, p.z, radius, 
                        trig.a.x, trig.a.y, trig.a.z,
                        trig.b.x, trig.b.y, trig.b.z,
                        trig.c.x, trig.c.y, trig.c.z
                    );
                    false_pos++;
                }   
            }
        }
    }
    stop_unix_timing(1);

    printf("False judged in total: false_pos = %d, flase_neg = %d (in %.2f ms)\n", false_pos, false_neg, get_time_elapse(1));
}

void test_omp(){
    omp_set_num_threads(128);
    cout << omp_get_num_threads() << endl;
    #pragma omp parallel 
    {
        #pragma omp for
        for(int i = 0; i < 10; i++){
            const int tid = omp_get_thread_num();
            printf("[%d] %d\n", i, tid);
            if(i == 0){
                printf("Total Threads: %d\n", omp_get_num_threads());
            }
        }
    }
}

/// (OMP) Match all positive results in the correct answer.
/// If an answer is written on the provided example result,
/// we have to find the match and verify it.
void test_collision_positives_omp(){
    vector<tuple<int, Point>> spheres;
    spheres = read_sphere_csv(sample_sphere_csv_path);

    vector<vector<tuple<int, Triangle>>> meshes;
    meshes = read_mesh_shapes(sample_mesh_obj_path);

    vector<tuple<int, int>> sid_tid_pairs;
    sid_tid_pairs = read_result_csv(sample_result_path);
    
    const float radius = 3; int false_judges = 0;
    #pragma omp for
    for(int i = 0; i < 39709; i++){
        int sid, tid;
        tie(sid, tid) = sid_tid_pairs[i];
        int v_sid = get<0>(spheres[sid]);
        int v_tid = get<0>(meshes[0][tid]);
        
        assert(v_sid == sid);
        assert(v_tid == tid);

        const Triangle trig = get<1>(meshes[0][tid]);
        const Point p = get<1>(spheres[sid]);
        
        if(Sphere::intersect(p, radius, trig) == false){
            false_judges++;
            printf("%d, nominal=(%d,%d), real=(%d,%d)\n", i, sid, tid, v_sid, v_tid);
        }
    }
    printf("False judged in total: %d\n", false_judges);

}

/// (OMP) Match all results.
void test_collision_all_omp(){
    vector<tuple<int, Point>> spheres;
    spheres = read_sphere_csv(sample_sphere_csv_path);
    const size_t ssize = spheres.size();
    assert(ssize == 1703244);

    vector<vector<tuple<int, Triangle>>> meshes;
    meshes = read_mesh_shapes(sample_mesh_obj_path);
    const size_t msize = meshes[0].size();
    assert(msize == 62976);

    vector<tuple<int, int>> sid_tid_pairs;
    sid_tid_pairs = read_result_csv(sample_result_path);
    
    // vector<bool> answer;
    bool * answer = (bool *) calloc(ssize * msize, sizeof(bool));
    cout << ssize * msize << endl;
    assert(answer[0] == false);
    assert(answer[rand()] == false);

    const float radius = 3; int false_pos = 0, false_neg = 0;
    
    
    for(int i = 0; i < sid_tid_pairs.size(); i++){
        int sid, tid;
        tie(sid, tid) = sid_tid_pairs[i];
        size_t hval = thash(sid, tid, msize);
        // printf("[%d] %d %d %lld\n", i, sid, tid, hval);
        answer[hval] = true;
        assert(answer[hval] == true);
    }

    omp_set_num_threads(128);
    cout << omp_get_num_threads() << endl;
    start_unix_timing(1);
    #pragma omp parallel 
    {
        #pragma omp for collapse(2)
        for(int j = 0; j < msize; j++){ // mesh
            for(int i = 0; i < ssize; i++){ // sphere
                if(i == 0){ printf(".");}
                int sid, tid;
                tid = get<0>(meshes[0][j]);
                sid = get<0>(spheres[i]);
                const Triangle & trig = get<1>(meshes[0][j]);
                const Point & p = get<1>(spheres[i]);
                bool flag = Sphere::intersect(p, radius, trig);
                if(answer[thash(sid, tid, msize)] != flag){
                    if (flag == false){
                        printf("/*c%s, %d %d*/Point p(%f,%f,%f);\tfloat radius = %f;\tTriangle trig(Point(%f,%f,%f), Point(%f,%f,%f), Point(%f,%f,%f));\n", 
                            flag ? "-" : "+", sid, tid,
                            p.x, p.y, p.z, radius, 
                            trig.a.x, trig.a.y, trig.a.z,
                            trig.b.x, trig.b.y, trig.b.z,
                            trig.c.x, trig.c.y, trig.c.z
                        );
                        false_neg++;

                    }else{
                        printf("/*c%s, %d %d*/Point p(%f,%f,%f);\tfloat radius = %f;\tTriangle trig(Point(%f,%f,%f), Point(%f,%f,%f), Point(%f,%f,%f));\n", 
                            flag ? "-" : "+", sid, tid,
                            p.x, p.y, p.z, radius, 
                            trig.a.x, trig.a.y, trig.a.z,
                            trig.b.x, trig.b.y, trig.b.z,
                            trig.c.x, trig.c.y, trig.c.z
                        );
                        false_pos++;
                    }   
                }
            }
        }
    }
    stop_unix_timing(1);

    printf("False judged in total: false_pos = %d, flase_neg = %d (in %.2f ms)\n", false_pos, false_neg, get_time_elapse(1));
}


int main(int argc, const char * argv[]){
    test_read_sphere_csv();
    test_read_mesh_shapes();
    test_read_result_csv();
    test_collision_positives();
    // test_collision_all(); // Too time consuming
    test_omp(); 
    test_collision_positives_omp();
    test_collision_all_omp();
}