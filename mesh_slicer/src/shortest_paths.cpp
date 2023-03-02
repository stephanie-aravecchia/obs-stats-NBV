
// Author(s) : Pierre Alliez

#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
// #include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/random_perturbation.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <boost/variant.hpp>



//typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
// typedef K::FT FT;
// typedef K::Point_3 Point;

typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::FT FT;
typedef CGAL::Surface_mesh<Kernel::Point_3> Triangle_mesh;
typedef Triangle_mesh::Vertex_index vertex_descriptor;
typedef Triangle_mesh::Face_index face_descriptor;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
typedef Traits::Barycentric_coordinates Barycentric_coordinates;
typedef boost::graph_traits<Triangle_mesh> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;
typedef Graph_traits::vertex_descriptor vertex_descriptor;
typedef Graph_traits::face_descriptor face_descriptor;
typedef Graph_traits::halfedge_descriptor halfedge_descriptor;

#if 1
typedef Kernel::Segment_3 Segment;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
// typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
#endif
typedef CGAL::AABB_face_graph_triangle_primitive<Triangle_mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> aabbTraits;
typedef CGAL::AABB_tree<aabbTraits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

namespace PMP = CGAL::Polygon_mesh_processing;

// A model of SurfacemeshShortestPathVisitor storing simplicies
// using boost::variant
struct Sequence_collector
{
    typedef boost::variant< vertex_descriptor,
            std::pair<halfedge_descriptor,double>,
            std::pair<face_descriptor, Barycentric_coordinates> > Simplex;
    std::list< Simplex > sequence;

    void operator()(halfedge_descriptor he, double alpha)
    {

        sequence.push_front( std::make_pair(he, alpha) );
    }

    void operator()(vertex_descriptor v)
    {
        sequence.push_front( v );
    }

    void operator()(face_descriptor f, Barycentric_coordinates alpha)
    {
        sequence.push_front( std::make_pair(f, alpha) );
    }
};

struct PathElement {
    enum {Face, Edge, Vertex} type;
    Point position;
    Vector uz; // x longitudinal, y lateral, z normal
};

struct PathSegment {
    Point start,end;
    Vector ux,uy,uz; // x longitudinal, y lateral, z normal

    void print() {
        std::cout << "Path segment from " << start << " to " << end << std::endl;
        std::cout << "UX: " << ux << std::endl;
        std::cout << "UY: " << uy << std::endl;
        std::cout << "UZ: " << uz << std::endl;
    }
};

// A visitor to print what a variant contains using boost::apply_visitor
struct Print_visitor : public boost::static_visitor<> {
    int i;
    const Triangle_mesh& g;
    const Surface_mesh_shortest_path &sp;
    Triangle_mesh::Property_map<face_descriptor, Vector> fnormals;
    Triangle_mesh::Property_map<vertex_descriptor, Vector> vnormals;
    std::list<PathElement> path;

    Print_visitor(const Triangle_mesh& g, const Surface_mesh_shortest_path & sp) :i(-1), g(g), sp(sp) {
        vnormals = g.property_map<vertex_descriptor,Vector>("v:normals").first;
        fnormals = g.property_map<face_descriptor,Vector>("f:normals").first;
    }

    void operator()(vertex_descriptor v)
    {
        PathElement pe;
        pe.type = PathElement::Face;
        pe.position = sp.point(v);
        pe.uz = vnormals[v];
        path.push_back(pe);
        std::cout << "#" << ++i << " : Vertex : " 
            << get(boost::vertex_index, g)[v] 
            << " -> " << pe.position
            << " / " << pe.uz
            << std::endl;
    }

    void operator()(const std::pair<halfedge_descriptor,double>& h_a)
    {
        PathElement pe;
        pe.type = PathElement::Edge;
        pe.position = sp.point(h_a.first,h_a.second);
        pe.uz = (fnormals[g.face(h_a.first)] + fnormals[g.face(g.opposite(h_a.first))])/2;
        path.push_back(pe);
        std::cout << "#" << ++i << " : Edge : " << get(CGAL::halfedge_index, g)[h_a.first] << " , ("
            << 1.0 - h_a.second << " , "
            << h_a.second << ") -> "
            << pe.position
            << " / " << pe.uz
            << std::endl;
    }

    void operator()(const std::pair<face_descriptor, Barycentric_coordinates>& f_bc)
    {
        PathElement pe;
        pe.type = PathElement::Face;
        pe.position = sp.point(f_bc.first,f_bc.second);
        pe.uz = fnormals[f_bc.first]; 
        path.push_back(pe);
        std::cout << "#" << ++i << " : Face : " << get(CGAL::face_index, g)[f_bc.first] << " , ("
            << f_bc.second[0] << " , "
            << f_bc.second[1] << " , "
            << f_bc.second[2] << ") -> "
            << pe.position
            << " / " << pe.uz
            << std::endl;
    }
};


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
    CGAL::copy_face_graph(polyhedron,tmesh);
    // Make sure we don't have really colinear edges... 
    PMP::random_perturbation( tmesh, 1e-3); 
    printf("Triangle mesh: %d faces, %d vertices\n",num_faces(tmesh),num_vertices(tmesh));

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

    // query point
    Point start(9.5, 0.0, 3.0);
    Point end(9.5, 35.0, 3.0);


    Surface_mesh_shortest_path shortest_paths(tmesh);
    std::cout << "Adding source" << std::endl;
    Surface_mesh_shortest_path::Face_location face_start = shortest_paths.locate<aabbTraits>(start);
    Surface_mesh_shortest_path::Face_location face_end = shortest_paths.locate<aabbTraits>(end);
    shortest_paths.add_source_point(face_start.first,face_start.second);
    std::cout << "Searching path" << std::endl;

    std::vector<Traits::Point_3> points;
    shortest_paths.shortest_path_points_to_source_points(face_end.first, face_end.second, std::back_inserter(points));
    std::cout << points.size() << std::endl;
    for (std::size_t i = 0; i < points.size(); ++i)
        std::cout << points[i] << std::endl;
    std::cout << std::endl;


    // collect the sequence of simplicies crossed by the shortest path
    Sequence_collector sequence_collector;
    shortest_paths.shortest_path_sequence_to_source_points(face_end.first, face_end.second, sequence_collector);

    // print the sequence using the visitor pattern
    Print_visitor print_visitor(tmesh,shortest_paths);
    for (auto it = sequence_collector.sequence.begin(); it != sequence_collector.sequence.end(); it++) {
        boost::apply_visitor(print_visitor, *it);
    }
    printf("Path element list has %d elements\n",int(print_visitor.path.size()));

    std::vector<PathSegment> path;

    // First iteration to add intermediary points on faces
    std::list<PathElement>::iterator it;
    for (it = print_visitor.path.begin(); it != print_visitor.path.end(); it++) {
        std::list<PathElement>::iterator itn=it;
        itn++;
        if (itn == print_visitor.path.end()) {
            break;
        }
        PathSegment ps;
        ps.start = it->position;
        ps.end = itn->position;
        ps.ux = ps.end - ps.start;
        ps.ux /= std::sqrt(ps.ux.squared_length());
        if (it->type == PathElement::Face) {
            ps.uz = it->uz;
        } else if (itn->type == PathElement::Face) {
            ps.uz = itn->uz;
        } else {
            // by construction C must be inside a face
            Point C = CGAL::ORIGIN + ((ps.end-CGAL::ORIGIN)+(ps.start-CGAL::ORIGIN))/2;
            Surface_mesh_shortest_path::Face_location face_center = shortest_paths.locate<aabbTraits>(C);
            ps.uz = fnormals[face_center.first];
        }
        ps.uy = CGAL::cross_product(ps.uz,ps.ux);
        ps.uy /= std::sqrt(ps.uy.squared_length());
        ps.print();
        path.push_back(ps);
    }

    return EXIT_SUCCESS;
}
