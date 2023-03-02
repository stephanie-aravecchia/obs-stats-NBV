
// Author(s) : Pierre Alliez

#include <iostream>
#include <fstream>
// #include <CGAL/IO/Polyhedron_iostream.h>

#include <Eigen/Core>

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




#include <boost/variant.hpp>



//typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
// typedef K::FT FT;
// typedef K::Point_3 Point;

typedef Kernel::Point_2 Point_2;
typedef Kernel::Vector_2 Vector_2;
typedef Kernel::Point_3 Point;
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


// static bool is_valid(face_descriptor f) {
//     return (f != Triangle_mesh::null_face());
// }

void reverse_points(PolyLine& polyline){
    PolyLine temp_line;
    for (auto it = polyline.rbegin(); it != polyline.rend(); it++){
        temp_line.push_back(*it);
    }

    polyline.clear();
    polyline = temp_line;

}

void add_intersect_vals(PolyLine& current, PolyLine& upcoming, Point p_intersect, Point start, Point end){
    if (start != p_intersect){
        current.push_back(start);
    }
    current.push_back(p_intersect);
    upcoming.push_back(p_intersect);

    if (end != p_intersect){
        upcoming.push_back(end);
    }

}


PolyLine combine_polylines(PolyLine& a, PolyLine& b){
    PolyLine result;
    if (b.size() > 1)
        b.pop_back();
    for (auto it = b.cbegin(); it != b.cend(); it++){
        result.push_back(*it);
    }

    for (auto it = a.cbegin(); it != a.cend(); it++){
        result.push_back(*it);
    }
    return result;
}

        
std::vector<PolyLine> plane_cut(const PolyLine& polyline, const Plane & plane){
    std::vector<PolyLine> lines;
    PolyLine path[3];
    unsigned int v_index = 0;

    if (polyline.size() <= 1) {
        lines.push_back(polyline);
        lines.push_back(polyline);
        return lines;
    }
    if (plane.has_on(polyline[0]) && plane.has_on(polyline[1])) {
        // This assumes that polyline is a convex object. 
        lines.push_back(polyline);
        lines.push_back(polyline);
        return lines;
    }
    path[v_index].push_back(polyline[0]);
    for (unsigned int i = 1; i < polyline.size(); i++){
        
        Segment seg(polyline[i-1],polyline[i]);
        if(do_intersect(seg, plane)){
            CGAL::cpp11::result_of<Intersect(Segment, Plane)>::type 
                result = intersection(seg, plane);
            const Point* p = boost::get<Point >(&*result);
            assert(v_index<3);
            add_intersect_vals(path[v_index], path[v_index+1], *p, polyline[i-1], polyline[i]);
            v_index++;
        } else {
            path[v_index].push_back(polyline[i]);
        }

    }
    

    PolyLine combined = combine_polylines(path[0], path[2]);

    // make sure it returns neg-x value size first (left), then positive
    // also make sure points go from front to back of the ship

    if (plane.has_on_positive_side(path[1][0])){
        lines.push_back(path[1]);
        lines.push_back(combined);
    } else {
        lines.push_back(combined);
        lines.push_back(path[1]);
    }

    return lines;

}

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

    FILE * fe = fopen("alledges","w");
    for(edge_descriptor ed : tmesh.edges()){
        Point A = tmesh.point(source(ed,tmesh));
        Point B = tmesh.point(target(ed,tmesh));
        fprintf(fe,"%e %e %e\n%e %e %e\n\n\n",
                A.x(),A.y(),A.z(),
                B.x(),B.y(),B.z());
    }
    fclose(fe);

    double zmin = NAN;
    for (vertex_descriptor vd : tmesh.vertices()) {
        Point P = tmesh.point(vd);
        if (std::isnan(zmin) || (P.z() < zmin)) {
            zmin = P.z();
        }
    }
    printf("Zmin = %f\n",zmin);


    // query point
    Point A(10.0, 0.2, zmin);
    Point B(-10.0, -0.2, zmin);
#if 0
    PolyLine segment;
    segment.push_back(A);
    segment.push_back(B);
    if (!CGAL::Polygon_mesh_processing::do_intersect(tmesh,segment)) {
        std::cout << "Segment is not intersecting the mesh, just do a straight line" << std::endl;
        return EXIT_SUCCESS;
    }
#endif

    Vector vAB = B-A;
    Plane pA(A,vAB);
    // Point_2 A2 = pA.to_2d(A);
    
    std::cout << "Segment: " << A << " , " << B << std::endl;

    CGAL::Polygon_mesh_slicer<Triangle_mesh, Kernel> slicer(tmesh);
    Vector v(0,0,1);
    // Vector u = CGAL::cross_product(vAB,v);
    Plane pCut(A,B,A+v);
    std::list<PolyLine> polylines;
    slicer(pCut, std::back_inserter(polylines));

    PolyLine2D pointlist, ch;
    std::vector<bool> onground;
    for (std::list<PolyLine>::const_iterator it=polylines.begin();it!=polylines.end();it++) {
        const PolyLine & pl = *it;
        for (size_t i=0;i<pl.size();i++) {
            Point pz(pl[i].x(),pl[i].y(),pl[i].z() - zmin);
            Point_2 p2 = pCut.to_2d(pz);
            std::cout << i << " " << p2 << std::endl; 
            pointlist.push_back(p2);
            onground.push_back(pz.z()==0.);
            Point p0(pz.x(),pz.y(),0);
            p2 = pCut.to_2d(p0);
            std::cout << i << "." << p2 << std::endl; 
            pointlist.push_back(p2);
            onground.push_back(true);
        }
    }
    std::vector<std::size_t> indices(pointlist.size()), out;
    std::iota(indices.begin(), indices.end(),0);
    CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
            Convex_hull_traits_2(CGAL::make_property_map(pointlist)));
    PolyLine ch3d,ch3d0;
    std::vector<bool> onground3;
    FILE *fch=fopen("chull","w");
    for (size_t i=0;i<out.size();i++) {
        const Point_2 & pi = pointlist[out[i]];
        ch.push_back(pi);
        Point P = pCut.to_3d(pi);
        ch3d0.push_back(P);
        ch3d.push_back(Point(P.x(),P.y(),P.z() + zmin));
        onground3.push_back(onground[out[i]]);
        fprintf(fch,"%e %e %e %e %d\n",P.x(),P.y(),P.z(),P.z() + zmin,onground[out[i]]?1:0);
    }
    fclose(fch);

    // Search for the first and last point on the ground
    size_t idx = 0;
    enum {Init, SearchingLiftOff, SearchingTouchDown, FoundBoth} state = Init;
    size_t idxLiftOff=0, idxTouchDown=0;
    bool foundLiftOff=false, foundTouchDown=false;
    while (state != FoundBoth) {
        switch (state) {
            case Init:
                if (onground3[idx]) {
                    state = SearchingLiftOff;
                } else {
                    state = SearchingTouchDown;
                }
                break;
            case SearchingLiftOff:
                if (!onground3[idx]) {
                    idxLiftOff = (idx + out.size() -1) % out.size();
                    foundLiftOff = true;
                    state = foundTouchDown?FoundBoth:SearchingTouchDown;
                }
                break;
            case SearchingTouchDown:
                if (onground3[idx]) {
                    idxTouchDown = idx;
                    foundTouchDown = true;
                    state = foundLiftOff?FoundBoth:SearchingLiftOff;
                }
                break;
            case FoundBoth:
                break;
        }
        idx = (idx + 1) % out.size();
    }
    assert(idxLiftOff != idxTouchDown);
    if (onground3[(idxLiftOff+1)%onground3.size()]) {
        std::swap(idxLiftOff,idxTouchDown);
    }
    printf("Idx Lift off %d Touch down %d\n",int(idxLiftOff),int(idxTouchDown));
    std::vector<size_t> path;
    for (size_t i=0;i<out.size();i++) {
        size_t ii = (idxLiftOff + i)%out.size();
        path.push_back(ii);
        if (ii == idxTouchDown) {
            break;
        }
    }
    if ((A - ch3d[path[0]]).squared_length() > (A - ch3d[path[path.size()-1]]).squared_length()) {
        std::reverse(path.begin(),path.end());
    }

    FILE *fp=fopen("path","w");
    fprintf(fp,"%e %e %e -1 -1\n",A.x(),A.y(),A.z());

    for (size_t i=0;i<path.size();i++) {
        const Point & P = ch3d[path[i]];
        fprintf(fp,"%e %e %e %d %d\n",P.x(),P.y(),P.z(),int(path[i]),onground3[path[i]]?1:0);
    }
    fprintf(fp,"%e %e %e -1 -1\n",B.x(),B.y(),B.z());


    fclose(fp);

    return EXIT_SUCCESS;
}
