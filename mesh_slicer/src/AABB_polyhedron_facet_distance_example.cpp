
// Author(s) : Pierre Alliez

#include <iostream>
#include <fstream>
#include <map>

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
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>


//typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
// typedef K::FT FT;
// typedef K::Point_3 Point;

typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;
typedef CGAL::Surface_mesh<Kernel::Point_3> Triangle_mesh;
// typedef Triangle_mesh::Facet_iterator                   T_Facet_iterator;
// typedef Triangle_mesh::Halfedge_around_facet_circulator T_Halfedge_facet_circulator;
typedef Triangle_mesh::Vertex_index vertex_descriptor;
typedef Triangle_mesh::Face_index face_descriptor;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
typedef boost::graph_traits<Triangle_mesh> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;

#if 1
typedef Kernel::Segment_3 Segment;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Facet_iterator                   P_Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator P_Halfedge_facet_circulator;
// typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
#endif
typedef CGAL::AABB_face_graph_triangle_primitive<Triangle_mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> aabbTraits;
typedef CGAL::AABB_tree<aabbTraits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

namespace PMP = CGAL::Polygon_mesh_processing;



int main(int argc,char * argv[])
{
    Polyhedron polyhedron;
    Triangle_mesh pmesh,tmesh;
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
    CGAL::copy_face_graph(polyhedron,pmesh);
    typedef std::multimap<vertex_descriptor,face_descriptor> VMap;
    VMap vmap;
    std::map<face_descriptor,face_descriptor> fmap;
#if 1
    {
        int i=0;
        std::cout << "TMesh Before: " << std::endl;
        BOOST_FOREACH(face_descriptor f, faces(pmesh)){
            std::cout << i << ": ";
            BOOST_FOREACH(vertex_descriptor vd,vertices_around_face(pmesh.halfedge(f), pmesh)){
                std::cout << vd << " ";
                vmap.insert(VMap::value_type(vd,f));
            }
            std::cout << std::endl;
            i++;
        }
    }
#endif

    CGAL::copy_face_graph(pmesh,tmesh);
    CGAL::Polygon_mesh_processing::triangulate_faces(tmesh);
    
#if 1
    {
        int i = 0;
        std::cout << "TMesh After: " << std::endl;
        BOOST_FOREACH(face_descriptor f, faces(tmesh)){
            std::cout << i << ": ";
            std::map<face_descriptor,unsigned int> fscore;
            BOOST_FOREACH(vertex_descriptor vd,vertices_around_face(tmesh.halfedge(f), tmesh)){
                for (VMap::const_iterator it=vmap.find(vd);(it!=vmap.end())&&(it->first==vd);it++) {
                    std::map<face_descriptor,unsigned int>::iterator sit = fscore.find(it->second);
                    if (sit == fscore.end()) {
                        fscore[it->second] = 1;
                    } else {
                        sit->second += 1;
                    }
                }
                std::cout << vd << " ";
            }
            for (std::map<face_descriptor,unsigned int>::const_iterator sit=fscore.begin();sit!=fscore.end();sit++) {
                if (sit->second == 3) {
                    fmap[f] = sit->first;
                    std::cout << "-> " << sit->first;
                    break;
                }
            }
            std::cout << std::endl;
            i++;
        }
    }
#endif



    //
    // Make sure we don't have really colinear edges... 
    PMP::random_perturbation( tmesh, 1e-3); 
    printf("Triangle mesh: %d faces, %d vertices\n",num_faces(tmesh),num_vertices(tmesh));

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

    // computes squared distance from query
    FT sqd1 = tree.squared_distance(start);
    FT sqd2 = tree.squared_distance(end);
    std::cout << "squared distance: " << sqd1 << " " << sqd2 << std::endl;

    // computes closest point
    Point closest_start = tree.closest_point(start);
    Point closest_end = tree.closest_point(end);
    std::cout << "closest points: " << closest_start << " -> " << closest_end << std::endl;

    // computes closest point and primitive id
#if 0
    Polyhedron::Face_handle f = pp.second; // closest primitive id
    std::cout << "closest triangle: ( "
              << f->halfedge()->vertex()->point() << " , " 
              << f->halfedge()->next()->vertex()->point() << " , "
              << f->halfedge()->next()->next()->vertex()->point()
              << " )" << std::endl;
#endif
    Point_and_primitive_id pp = tree.closest_point_and_primitive(start);
    closest_start = pp.first;
    CGAL::SM_Face_index f_start = pp.second; // closest primitive id
    pp = tree.closest_point_and_primitive(end);
    closest_end = pp.first;
    CGAL::SM_Face_index f_end = pp.second; // closest primitive id
    std::cout << "closest points: " << closest_start << " -> " << closest_end << std::endl;
    std::cout << "Face index " << f_start << " (" << fmap[f_start] << ")" << std::endl;
    BOOST_FOREACH(vertex_descriptor vd,vertices_around_face(tmesh.halfedge(f_start), tmesh)){
        std::cout << vd << ":" << tmesh.point(vd) << std::endl;
    }
    std::cout  << " -> " << f_end << " ("<<fmap[f_end]<<")"<< std::endl;
    BOOST_FOREACH(vertex_descriptor vd,vertices_around_face(tmesh.halfedge(f_end), tmesh)){
        std::cout << vd << ":" << tmesh.point(vd) << std::endl;
    }
    // Traits::Barycentric_coordinates face_location = {{0.25, 0.5, 0.25}};

    Surface_mesh_shortest_path shortest_paths(tmesh);
    std::cout << "Adding source" << std::endl;
    face_iterator face_it = faces(tmesh).first;
    std::advance(face_it,f_start);
    // shortest_paths.add_source_point(*face_it, face_location);
    Surface_mesh_shortest_path::Face_location face_start = shortest_paths.locate<aabbTraits>(start);
    Surface_mesh_shortest_path::Face_location face_end = shortest_paths.locate<aabbTraits>(end);
    shortest_paths.add_source_point(face_start.first,face_start.second);
    std::cout << "Searching path" << std::endl;
    face_it = faces(tmesh).first;
    std::advance(face_it,f_end);
    std::vector<Traits::Point_3> points;
    // shortest_paths.shortest_path_points_to_source_points(*face_it, face_location, std::back_inserter(points));
    shortest_paths.shortest_path_points_to_source_points(face_end.first, face_end.second, std::back_inserter(points));
    std::cout << points.size() << std::endl;
    for (std::size_t i = 0; i < points.size(); ++i) {
        std::cout << points[i] << std::endl;
    }
    std::cout << std::endl;
    return EXIT_SUCCESS;
}
