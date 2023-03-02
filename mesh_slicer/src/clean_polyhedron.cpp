
// Author(s) : Pierre Alliez

#include <iostream>
#include <fstream>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>


typedef CGAL::Simple_cartesian<double> Kernel;
// typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
// typedef K::FT FT;
// typedef K::Point_3 Point;

typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;
typedef CGAL::Surface_mesh<Kernel::Point_3> Triangle_mesh;
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
// typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
#endif
typedef CGAL::AABB_face_graph_triangle_primitive<Triangle_mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> aabbTraits;
typedef CGAL::AABB_tree<aabbTraits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef Polyhedron::Vertex_handle   Vertex_handle;




int main(int argc,char * argv[])
{
    Triangle_mesh tmesh;
    Polyhedron polyhedron, polyhedron2;
    if (argc>2) {
        std::ifstream input(argv[1]);
        input >> polyhedron;
        input.close();
        printf("Polyhedron mesh: %d faces, %d vertices\n",int(num_faces(polyhedron)),int(num_vertices(polyhedron)));
        CGAL::Polygon_mesh_processing::stitch_borders(polyhedron);
        CGAL::Polygon_mesh_processing::triangulate_faces(polyhedron);
        CGAL::copy_face_graph(polyhedron,tmesh);
        printf("Polyhedron mesh: %d faces, %d vertices\n",int(num_faces(tmesh)),int(num_vertices(tmesh)));
        std::ofstream output(argv[2]);
        output  << tmesh;
        output.close();
        output.flush();
        // printf("Trying to read it again...\n");
        // std::ifstream input2(argv[2]);
        // scan_OFF(input2,polyhedron2, true);
        // // input2 >> polyhedron2;
        // input2.close();
        // printf("Polyhedron mesh: %d faces, %d vertices\n",int(num_faces(polyhedron2)),int(num_vertices(polyhedron2)));
    } else {
        printf("Usage: %s in out\n",argv[0]);
        return -1;
    }
    return EXIT_SUCCESS;
}
