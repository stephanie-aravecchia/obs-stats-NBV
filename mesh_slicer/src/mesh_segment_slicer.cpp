
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

float dist(const DistancePair & dp) {
    return (dp.first + dp.second)/2;
}

void distanceMap(const PointEdgeMap & start, const EdgeGraph & G, 
        PointEdgeMap & edge_costs, PredecessorMap & predecessor) {
    predecessor.clear();
    edge_costs.clear();
    if ((G.size()==0) || (start.size()==0)) {
        return;
    }
    typedef std::multimap<float, edge_descriptor> Heap;

    Heap heap;
    for (PointEdgeMap::const_iterator it=start.begin();it!=start.end();it++) {
        heap.insert(Heap::value_type(dist(it->second),it->first));
        edge_costs.insert(*it);
    }
    while (!heap.empty()) {
        Heap::iterator ic = heap.begin();
        edge_descriptor current = ic->second;
        //std::cout << "Pop " << current << std::endl;
        EdgeGraph::const_iterator git=G.find(current);
        assert(git!=G.end());
        PointEdgeMap::const_iterator rit = edge_costs.find(current);
        assert(rit != edge_costs.end());
        float current_cost = dist(rit->second);
        for (PointEdgeMap::const_iterator nit=git->second.begin();nit!=git->second.end();nit++) {
            assert(dist(nit->second) > 0);
            float cost_proposal = current_cost + dist(nit->second);
            rit = edge_costs.find(nit->first);
            if ((rit == edge_costs.end()) || (cost_proposal < dist(rit->second))) {
                predecessor[nit->first] = current;
                edge_costs[nit->first] = DistancePair(cost_proposal,cost_proposal);
                //std::cout << "Adding " << current << " -> " << nit->first << std::endl; 
                heap.insert(Heap::value_type(cost_proposal,nit->first));
            }
        }
        heap.erase(ic);
    }
}

void shortestPaths(const PointEdgeMap & start, const EdgeGraph & G, 
        const PointEdgeMap & end, std::vector<EdgePath> & edge_paths, 
        float cost_tolerance, const Triangle_mesh & tmesh) {
    
    edge_paths.clear();
    PointEdgeMap cost_from_start, cost_from_end, combined_cost;
    PredecessorMap pred_from_start, pred_from_end;
    printf("Computing distance maps\n");
    distanceMap(start,G,cost_from_start,pred_from_start);
    printf("Completed start map\n");
    distanceMap(end,G,cost_from_end,pred_from_end);
    printf("Completed end map\n");
    float min_cost = NAN;
    FILE * fe;
    fe = fopen("edgecosts_start","w");
    for (PointEdgeMap::const_iterator sit=cost_from_start.begin();sit!=cost_from_start.end();sit++) {
        Point A1 = tmesh.point(source(sit->first,tmesh));
        Point B1 = tmesh.point(target(sit->first,tmesh));
        fprintf(fe,"%f %f %f %f\n%f %f %f %f\n\n\n",
                A1.x(),A1.y(),A1.z(),dist(sit->second),
                B1.x(),B1.y(),B1.z(),dist(sit->second));
    }
    fe = fopen("edgecosts_end","w");
    for (PointEdgeMap::const_iterator sit=cost_from_end.begin();sit!=cost_from_end.end();sit++) {
        Point A1 = tmesh.point(source(sit->first,tmesh));
        Point B1 = tmesh.point(target(sit->first,tmesh));
        fprintf(fe,"%f %f %f %f\n%f %f %f %f\n\n\n",
                A1.x(),A1.y(),A1.z(),dist(sit->second),
                B1.x(),B1.y(),B1.z(),dist(sit->second));
    }
    fclose(fe);
    fe = fopen("edgecosts","w");
    for (PointEdgeMap::const_iterator sit=cost_from_start.begin();sit!=cost_from_start.end();sit++) {
        PointEdgeMap::const_iterator eit=cost_from_end.find(sit->first);
        if (eit == cost_from_end.end()) {
            // no path from start to end through this edge
            continue;
        }
        Point A1 = tmesh.point(source(sit->first,tmesh));
        Point B1 = tmesh.point(target(sit->first,tmesh));
        float cost = dist(sit->second) + dist(eit->second);
        combined_cost[eit->first] = DistancePair(cost,cost);
        fprintf(fe,"%f %f %f %f %f %f\n%f %f %f %f %f %f\n\n\n",
                A1.x(),A1.y(),A1.z(),dist(sit->second),dist(eit->second),cost, 
                B1.x(),B1.y(),B1.z(),dist(sit->second),dist(eit->second),cost);
        if (std::isnan(min_cost) || (cost < min_cost)) {
            min_cost = cost;
        }
    }
    fclose(fe);
    printf("Minimum Cost: %f\n",min_cost);
    for (PointEdgeMap::const_iterator cit=combined_cost.begin();cit!=combined_cost.end();cit++) {
        if (dist(cit->second) < 0) {
            continue;
        }
        if (dist(cit->second) <= (min_cost + cost_tolerance)) {
            std::list<edge_descriptor> epstart, epend;
            PredecessorMap::const_iterator pit = pred_from_start.find(cit->first);
            while (pit != pred_from_start.end()) {
                combined_cost[pit->second] = DistancePair(-1,-1);
                epstart.push_front(pit->second);
                pit = pred_from_start.find(pit->second);
            }
            pit = pred_from_end.find(cit->first);
            while (pit != pred_from_end.end()) {
                combined_cost[pit->second] = DistancePair(-1,-1);
                epend.push_back(pit->second);
                pit = pred_from_end.find(pit->second);
            }
            EdgePath ep;
            std::copy(epstart.begin(),epstart.end(),std::back_inserter(ep));
            combined_cost[cit->first] = DistancePair(-1,-1);
            ep.push_back(cit->first);
            std::copy(epend.begin(),epend.end(),std::back_inserter(ep));

            edge_paths.push_back(ep);
        }
    }
    printf("Selected %d low cost paths\n",int(edge_paths.size()));
}

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


    // query point
    Point A(10.0, 0.2, 0.0);
    Point B(-10.0, -0.2, 0.0);
    PolyLine segment;
    segment.push_back(A);
    segment.push_back(B);
    if (!CGAL::Polygon_mesh_processing::do_intersect(tmesh,segment)) {
        std::cout << "Segment is not intersecting the mesh, just do a straight line" << std::endl;
        return EXIT_SUCCESS;
    }

    Vector vAB = B-A;
    Plane pA(A,vAB);
    // Point_2 A2 = pA.to_2d(A);
    
    std::cout << "Segment: " << A << " , " << B << std::endl;

#if 0
    FILE * fe2 = fopen("edge2D","w");
    for(edge_descriptor ed : tmesh.edges()){
        std::cout << tmesh.number_of_edges() << std::endl;
        Point a = tmesh.point(source(ed,tmesh));
        Point b = tmesh.point(target(ed,tmesh));
        Vector_2 a2 = pA.to_2d(a) - A2;
        Vector_2 b2 = pA.to_2d(b) - A2;
        double theta_a = atan2(a2.y(),a2.x());
        double theta_b = atan2(b2.y(),b2.x());
        fprintf(fe2,"%e %e %e %e\n%e %e %e %e\n\n\n",a2.x(),a2.y(),cos(theta_a),sin(theta_a),b2.x(),b2.y(),cos(theta_b),sin(theta_b));
        Triangle_mesh lemon(tmesh);
        std::cout << "Edge: " << a << " , " << b << std::endl;
        std::cout << "Cutting plane1: " << A << " , " << B << " , " << a << std::endl;
        std::cout << "Cutting plane2: " << A << " , " << b << " , " << B << std::endl;
        CGAL::Polygon_mesh_processing::clip(lemon,Plane(A,a,B),CGAL::parameters::clip_volume(false));
        CGAL::Polygon_mesh_processing::triangulate_faces(lemon);
        CGAL::Polygon_mesh_processing::clip(lemon,Plane(A,B,b),CGAL::parameters::clip_volume(false));
        CGAL::Polygon_mesh_processing::triangulate_faces(lemon);

        
        std::vector<Point> lemonpoints;
        for(vertex_descriptor vd : lemon.vertices()){
            lemonpoints.push_back(lemon.point(vd));
        }
        // define polyhedron to hold convex hull
        Polyhedron lemonh;
        // compute convex hull of non-collinear points
        CGAL::convex_hull_3(lemonpoints.begin(), lemonpoints.end(), lemonh);

        {
            std::ofstream output("lemon.off");
            output  << lemon;
            output.close();
            output.flush();
        }
        {
            std::ofstream output("lemonh.off");
            output  << lemonh;
            output.close();
            output.flush();
        }
        getchar();
    }
    fclose(fe2);
#endif



    Vector ub1 = pA.base1(), ub2 = pA.base2();
    FILE *fch=fopen("chull","w");
    FILE *fcu=fopen("cull","w");
    FILE *fpl=fopen("slice","w");
    CGAL::Polygon_mesh_slicer<Triangle_mesh, Kernel> slicer(tmesh);
    double theta_range = M_PI;
    double theta_ref = 0;
    for (size_t k=1;k<4;k++) { 
        double best_length = NAN;
        for (size_t ii=0;ii<200;ii++) {
            double theta = theta_ref - theta_range/2. + ii*theta_range/200.;
            Vector u = ub1 * cos(theta) + ub2 * sin(theta);
            Vector v = CGAL::cross_product(vAB,u);
            Plane pCut(A,B,A+u);
            Plane pHalf(A,B,A+v);
            std::list<PolyLine> polylines;
#if 0
            Triangle_mesh lemon(tmesh);
            CGAL::Polygon_mesh_processing::clip(lemon,pHalf,CGAL::parameters::clip_volume(false));
            {
                std::ofstream output("lemon.off");
                output  << lemon;
                output.close();
                output.flush();
            }
            CGAL::Polygon_mesh_slicer<Triangle_mesh, Kernel> slicer(lemon);
#endif
            slicer(pCut, std::back_inserter(polylines));

            PolyLine2D pointlist, ch;
            pointlist.push_back(pCut.to_2d(A));
            pointlist.push_back(pCut.to_2d(B));
            for (std::list<PolyLine>::const_iterator it=polylines.begin();it!=polylines.end();it++) {
                const PolyLine & pl = *it;
                for (size_t i=0;i<pl.size();i++) {
                    Point_2 p2 = pCut.to_2d(pl[i]);
                    // fprintf(fpl,"%e %e %e %e %e %e\n",theta,pl[i].x(),pl[i].y(),pl[i].z(),p2.x(),p2.y());
                    pointlist.push_back(p2);
                }
            }
            // fprintf(fpl,"\n\n");
            // fflush(fpl);
            std::vector<std::size_t> indices(pointlist.size()), out;
            std::iota(indices.begin(), indices.end(),0);
            CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
                    Convex_hull_traits_2(CGAL::make_property_map(pointlist)));
            // CGAL::ch_graham_andrew( pointlist.begin(), pointlist.end(), std::back_inserter(ch));
            int i0 = -1, i1 = -1;
            for (size_t i=0;i<out.size();i++) {
                if (out[i]==0) { i0 = i; }
                if (out[i]==1) { i1 = i; }
                ch.push_back(pointlist[out[i]]);
            }
            // This needs to be managed specially, but for now, we assume we are
            // far enough that the start and end point are garanteed to be on the
            // CH
            assert(i0>=0);
            assert(i1>=0);
            PolyLine ch3d;
            for (size_t i=0;i<ch.size();i++) {
                Point P = pCut.to_3d(ch[i]);
                ch3d.push_back(P);
            }

#if 0
            for (size_t i=0;i<ch3d.size();i++) {
                fprintf(fch,"%e %e %e %e %e %e\n",theta,ch3d[i].x(),ch3d[i].y(),ch3d[i].z(),ch[i].x(),ch[i].y());
            }
            fprintf(fch,"%e %e %e %e %e %e\n",theta,ch3d[0].x(),ch3d[0].y(),ch3d[0].z(),ch[0].x(),ch[0].y());
            fprintf(fch,"\n\n");
            fflush(fch);
#endif

            PolyLine ch0,ch1;
            for (size_t i=0;i<ch3d.size();i++) {
                size_t ii = (i0 + i) % ch3d.size();
                ch0.push_back(ch3d[ii]);
                if (int(ii) == i1) break;
            }
            for (size_t i=0;i<ch3d.size();i++) {
                size_t ii = (i1 + i) % ch3d.size();
                ch1.push_back(ch3d[ii]);
                if (int(ii) == i0) break;
            }

            double l0=0, l1=0;
            for (size_t i=1;i<ch0.size();i++) {
                l0 += sqrt((ch0[i]-ch0[i-1]).squared_length());
            }
            for (size_t i=1;i<ch1.size();i++) {
                l1 += sqrt((ch1[i]-ch1[i-1]).squared_length());
            }
            if (std::isnan(best_length) || (l0 < best_length)) {
                theta_ref = theta;
                best_length = l0;
            }
            if (std::isnan(best_length) || (l1 < best_length)) {
                theta_ref = theta;
                best_length = l1;
            }
            for (size_t i=0;i<ch0.size();i++) {
                fprintf(fcu,"%e %e %e %e %e \n",remainder(theta,2*M_PI),ch0[i].x(),ch0[i].y(),ch0[i].z(),l0);
            }
            fprintf(fcu,"\n\n");
            for (size_t i=0;i<ch1.size();i++) {
                fprintf(fcu,"%e %e %e %e %e \n",remainder(theta+M_PI,2*M_PI),ch1[i].x(),ch1[i].y(),ch1[i].z(),l1);
            }
            fprintf(fcu,"\n\n");
            fflush(fcu);

        }
        theta_range /= 10;
    }
    fclose(fcu);
    fclose(fch);
    fclose(fpl);

    return EXIT_SUCCESS;
}
