
// Author(s) : Pierre Alliez

#include <iostream>
#include <fstream>
// #include <CGAL/IO/Polyhedron_iostream.h>

#include <Eigen/Core>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
// #include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/random_perturbation.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/convex_hull_3.h>

#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/property_map.h>
#include <vector>
#include <numeric>



#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>


#include <boost/variant.hpp>



//typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
// typedef K::FT FT;
// typedef K::Point_3 Point;

typedef Kernel::Point_2 Point_2;
typedef Kernel::Vector_2 Vector_2;
typedef Kernel::Point_3 Point;
typedef CGAL::Bbox_3 Bbox;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Facet;
typedef Kernel::Tetrahedron_3 Tetra;
typedef Kernel::Intersect_3 Intersect;

typedef Kernel::FT FT;
typedef CGAL::Surface_mesh<Kernel::Point_3> Triangle_mesh;
typedef Triangle_mesh::Vertex_index vertex_descriptor;
typedef Triangle_mesh::Face_index face_descriptor;
typedef Triangle_mesh::Edge_index edge_descriptor;
typedef Triangle_mesh::Halfedge_index halfedge_descriptor;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
typedef Traits::Barycentric_coordinates Barycentric_coordinates;
typedef boost::graph_traits<Triangle_mesh> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;
// typedef Graph_traits::vertex_descriptor vertex_descriptor;
// typedef Graph_traits::face_descriptor face_descriptor;
// typedef Graph_traits::halfedge_descriptor halfedge_descriptor;
typedef CGAL::Convex_hull_traits_adapter_2<Kernel,
          CGAL::Pointer_property_map<Point_2>::type > Convex_hull_traits_2;

#if 1
typedef Kernel::Segment_3 Segment;
typedef Kernel::Ray_3 Ray;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
// typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
#endif
typedef CGAL::AABB_face_graph_triangle_primitive<Triangle_mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> aabbTraits;
typedef CGAL::AABB_tree<aabbTraits> Tree;
typedef boost::optional< Tree::Intersection_and_primitive_id<Facet>::Type > Facet_intersection;
typedef boost::optional< Tree::Intersection_and_primitive_id<Tetra>::Type > Tetra_intersection;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef std::vector<Point> PolyLine;
typedef std::vector<Point_2> PolyLine2D;

typedef std::pair<edge_descriptor,edge_descriptor> EdgePair;
typedef std::pair<float,float> DistancePair;
typedef std::map<edge_descriptor,DistancePair> PointEdgeMap;
// typedef std::map<EdgePair,DistancePair> EdgePairMap;
typedef std::map<edge_descriptor,PointEdgeMap> EdgeGraph;
typedef std::map<edge_descriptor, edge_descriptor> PredecessorMap;
typedef std::vector<edge_descriptor> EdgePath;

namespace PMP = CGAL::Polygon_mesh_processing;


int main(int argc,char * argv[])
{
    Polyhedron polyhedron;
    Triangle_mesh tmesh;
    if (argc>1) {
        std::ifstream input(argv[1]);
        // input >> polyhedron;
        scan_OFF(input,polyhedron, true);
        input.close();
        // std::ifstream input2(argv[1]);
        // input2 >> tmesh;
        // input2.close();
        printf("Polyhedron mesh: %d faces, %d vertices\n",int(num_faces(polyhedron)),int(num_vertices(polyhedron)));
    } else {
        printf("Usage: %s <mesh file>\n",argv[0]);
        return -1;
    }
    // CGAL::Polygon_mesh_processing::stitch_borders(polyhedron);
    CGAL::Polygon_mesh_processing::triangulate_faces(polyhedron);
    CGAL::copy_face_graph(polyhedron,tmesh);
    // Make sure we don't have really colinear edges... 
    PMP::random_perturbation( tmesh, 1e-3); 
    printf("Triangle mesh: %d faces, %d vertices\n",num_faces(tmesh),num_vertices(tmesh));

#if 0
    auto fnormals = tmesh.add_property_map<face_descriptor, Vector>
        ("f:normals", CGAL::NULL_VECTOR).first;
    auto vnormals = tmesh.add_property_map<vertex_descriptor, Vector>
        ("v:normals", CGAL::NULL_VECTOR).first;
    PMP::compute_normals(tmesh, vnormals, fnormals,
            PMP::parameters::vertex_point_map(tmesh.points()).
            geom_traits(Kernel()));

    // constructs AABB tree and computes internal KD-tree 
    // data structure to accelerate distance queries
    // Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
    std::cout << "Creating AABB tree" << std::endl;
    Tree tree(faces(tmesh).first, faces(tmesh).second, tmesh);
    tree.accelerate_distance_queries();
    std::cout << "\tdone" << std::endl;
#endif

    double resolution = 0.05;
    if (argc>2) {
        resolution = atof(argv[2]);
    }
    Bbox bb;
    for (vertex_descriptor vd : tmesh.vertices()) {
        Point P = tmesh.point(vd);
        bb+=P.bbox();
    }
    Point Pmax(bb.xmax(),bb.ymax(),bb.zmax());
    Point Pmin(bb.xmin(),bb.ymin(),bb.zmin());
    Vector Vdiag = Pmax - Pmin;
    // Vdiag = Vector(Vdiag.x(),Vdiag.y(),0); 
    double l = sqrt(Vdiag.squared_length());
    Vdiag = Vdiag / l;
    bb += (Pmax + (5*resolution)*Vdiag).bbox();
    bb += (Pmin - (5*resolution)*Vdiag).bbox();
    //above is original,
    //here we try to make it bigger to test our comparison
    //bb += (Pmax + (121*resolution)*Vdiag).bbox();
    //bb += (Pmin - (47*resolution)*Vdiag).bbox();
    cv::Size size((bb.xmax()-bb.xmin())/resolution, (bb.ymax()-bb.ymin())/resolution); 



    CGAL::Polygon_mesh_slicer<Triangle_mesh, Kernel> slicer(tmesh);
    Vector vZ(0,0,1);

    std::ofstream of("bbox");
    of << "xmin ymin zmin xmax ymax zmax res" << std::endl;
    of << bb << " " << resolution << std::endl;
    of.close();

    std::cout << "BBox: " << bb << std::endl;
    std::cout << "Resolution: " << resolution << std::endl;

    unsigned int slice_index = 0 ;
    for (double z=bb.zmin();z <= bb.zmax(); z+= resolution) {
        std::cout << "Slice index " << slice_index << ": " << z << std::endl;
        cv::Mat_<uint8_t> slice(size,0);
        Point P((bb.xmax()+bb.xmin())/2,(bb.ymax()+bb.ymin())/2,z);
        Plane pCut(P,vZ);
        std::list<PolyLine> polylines;
        slicer(pCut, std::back_inserter(polylines));

        unsigned int segment_count = 0;
        for (std::list<PolyLine>::const_iterator it=polylines.begin();it!=polylines.end();it++) {
            const PolyLine & pl = *it;
            for (size_t i=1;i<pl.size();i++) {
                // std::cout << i << ": " << pl[i-1] << " -> " << pl[i] << std::endl;
                cv::line(slice,cv::Point((pl[i-1].x()-bb.xmin())/resolution,(pl[i-1].y()-bb.ymin())/resolution),
                        cv::Point((pl[i].x()-bb.xmin())/resolution,(pl[i].y()-bb.ymin())/resolution),255,1);
                segment_count++;
            }
        }
        std::cout << polylines.size() << " polylines; ";
        std::cout << segment_count << " segments" << std::endl;
        char fname[128];
        sprintf(fname,"slice%04d.png",slice_index);
        cv::imwrite(fname,slice);

        slice_index += 1;
    }

    return EXIT_SUCCESS;
}
