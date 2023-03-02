
// Author(s) : Pierre Alliez

#include <iostream>
#include <fstream>

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




#include <boost/variant.hpp>



//typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
// typedef K::FT FT;
// typedef K::Point_3 Point;

typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Facet;
typedef Kernel::Tetrahedron_3 Tetra;
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

static bool is_valid(face_descriptor f) {
    return (f != Triangle_mesh::null_face());
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

    // query point
    Point start(10.0, 0.2, 0.0);
    Point end(-10.0, -0.2, 0.0);
    PolyLine segment;
    segment.push_back(start);
    segment.push_back(end);

    if (CGAL::Polygon_mesh_processing::do_intersect(tmesh,segment)) {
        std::cout << "Segment is intersecting the mesh" << std::endl;
        Ray ray(start,end);
        CGAL::Polygon_mesh_processing::Face_location<Triangle_mesh, FT> fl = CGAL::Polygon_mesh_processing::locate(ray,tmesh);
        std::cout << "Ray casting impactpoints: " << fl.first <<":("
            <<fl.second[0]<<","<<fl.second[1]<<","<<fl.second[2]<<")" << std::endl;
    } else {
        std::cout << "Segment is NOT intersecting the mesh" << std::endl;
    }

    CGAL::Side_of_triangle_mesh<Triangle_mesh, Kernel> inside(tmesh);


    PointEdgeMap startmap, endmap;
    EdgeGraph epm;
    std::set<edge_descriptor> virtual_edges;

    // Instead of the classical for loop one can use
    // the boost macro for a range
    //size_t i=0;
    std::cout << tmesh.number_of_edges() << " Mesh edges:" << std::endl;
    FILE * fe = fopen("alledges","w");
    for(edge_descriptor ed : tmesh.edges()){
        halfedge_descriptor e0 = tmesh.halfedge(ed,0);
        halfedge_descriptor e1 = tmesh.halfedge(ed,1);
        face_descriptor f0 = tmesh.face(e0);
        face_descriptor f1 = tmesh.face(e1);
        Vector n0,n1;
        if (!is_valid(f0) && !is_valid(f1)) {
            // not sure what this means, but this is not a useful edge
            continue;
        } else if (!is_valid(f0)) {
            n1 = fnormals[f1];
            n0 = -n1;
        } else if (!is_valid(f1)) {
            n0 = fnormals[f0];
            n1 = -n0;
        } else {
            n0 = fnormals[f0];
            n1 = fnormals[f1];
        }

        Point A = tmesh.point(source(ed,tmesh));
        Point B = tmesh.point(target(ed,tmesh));
        fprintf(fe,"%e %e %e\n%e %e %e\n\n\n",
                A.x(),A.y(),A.z(),
                B.x(),B.y(),B.z());
        std::cout << "Edge ["<<A<<" , "<<B<<"]"<<std::endl;
        if (CGAL::cross_product(n0,n1).squared_length() < 1e-6) {
            std::cout << "Virtual" << std::endl;
            virtual_edges.insert(ed);
            continue;
        }

        // std::cout << "Normals ["<<n0<<" , "<<n1<<"]"<<std::endl;
        Segment SAB(A,B);
        Point M = A+(B-A)/2;
        // std::cout << "M: " << M << std::endl;

        Point test[2] = {start,end};
        float dp[2] = {-1,-1};
        for (int i=0;i<2;i++) {
            segment[1]=test[i];
            // std::cout << "Testing "<<test[i]<<std::endl;

            Vector vis = test[i] - M;
            // std::cout << "Vis: " << vis << std::endl;
            if (((vis * n0) < 0) && ((vis * n1) < 0)) {
                std::cout << test[i] << " is below the f0 ("<<vis*n0<<") && f1 ("<<vis*n1<<") plane" << std::endl;
                continue;
            }

            Vector u;
            
            u = (test[i]-A);
            double lA = std::sqrt(u.squared_length());
#if 0
            u /= lA;
            segment[0]=A+u*1e-3;
            if (CGAL::Polygon_mesh_processing::do_intersect(tmesh,segment)) {
                // std::cout << "Segment [A,end] is intersecting the mesh" << std::endl;
                continue;
            }
#endif
            u = (test[i]-B);
            double lB = std::sqrt(u.squared_length());
#if 0
            u /= lB;
            segment[0]=B+u*1e-3;
            if (CGAL::Polygon_mesh_processing::do_intersect(tmesh,segment)) {
                //std::cout << "Segment [B,end] is intersecting the mesh" << std::endl;
                continue;
            }
#endif
#if 1
            u = (test[i]-M);
            double lM = std::sqrt(u.squared_length());
            u /= lM;
            segment[0]=M+u*1e-3;
            if (inside(segment[0]) == CGAL::ON_BOUNDED_SIDE) {
                // std::cout << "Segment [M,end] starts inside the mesh" << std::endl;
                continue;
            }
            if (CGAL::Polygon_mesh_processing::do_intersect(tmesh,segment)) {
                // std::cout << "Segment [M,end] is intersecting the mesh" << std::endl;
                continue;
            }
#endif

            dp[i] =  std::max<float>(lA,lB);
        }
#if 0
        std::cout << "Distances : start " << sqrt(CGAL::squared_distance(start,SAB)) << " end " << sqrt(CGAL::squared_distance(end,SAB)) << std::endl;
        if (dp[0] >= 0) {
            startmap[ed] = DistancePair(sqrt(CGAL::squared_distance(start,SAB)),
                    dp[0]);
        }
        if (dp[1] >= 0) {
            endmap[ed] = DistancePair(sqrt(CGAL::squared_distance(end,SAB)),
                    dp[1]);
        }
#else
        std::cout << "Distances : start " << sqrt(CGAL::squared_distance(start,M)) << " end " << sqrt(CGAL::squared_distance(end,M)) << std::endl;
        if (dp[0] >= 0) {
            double d = sqrt(CGAL::squared_distance(start,M));
            startmap[ed] = DistancePair(d,d);
        }
        if (dp[1] >= 0) {
            double d = sqrt(CGAL::squared_distance(start,M));
            endmap[ed] = DistancePair(d,d);
        }
#endif

    }
    fclose(fe);

    std::cout << "Start map size: " << startmap.size() << std::endl;
    std::cout << "End map size: " << endmap.size() << std::endl;


    size_t edge_total = tmesh.number_of_edges()*(tmesh.number_of_edges()-1)/2;
    size_t edge_pairs = 0;
    size_t same_face = 0;
    size_t not_visible = 0;
    size_t intersecting = 0;
    size_t fully_visible = 0;
    float target_percent = 10;
    for(edge_descriptor e1 : tmesh.edges()){
        // if (virtual_edges.find(e1) != virtual_edges.end()) {
        //     continue;
        // }
        halfedge_descriptor e10 = tmesh.halfedge(e1,0);
        halfedge_descriptor e11 = tmesh.halfedge(e1,1);
        face_descriptor f10 = tmesh.face(e10);
        face_descriptor f11 = tmesh.face(e11);
        Vector n10,n11;
        if (!is_valid(f10) && !is_valid(f11)) {
            // not sure what this means, but this is not a useful edge
            continue;
        } else if (!is_valid(f10)) {
            n11 = fnormals[f11];
            n10 = -n11;
        } else if (!is_valid(f11)) {
            n10 = fnormals[f10];
            n11 = -n10;
        } else {
            n10 = fnormals[f10];
            n11 = fnormals[f11];
        }

        Point A1 = tmesh.point(source(e1,tmesh));
        Point B1 = tmesh.point(target(e1,tmesh));
        // Vector vAB1 = B1 - A1;
        Segment SAB1(A1,B1);
        Point M1 = A1+(B1-A1)/2;
        for(edge_descriptor e2 : tmesh.edges()){
            if (e2 == e1) {
                break;
                continue;
            }
            // if (virtual_edges.find(e2) != virtual_edges.end()) {
            //     continue;
            // }
            edge_pairs += 1;
            if (((100.*edge_pairs)/edge_total)>=target_percent) {
                std::cout << target_percent << std::endl;
                target_percent += 10;
            }

            halfedge_descriptor e20 = tmesh.halfedge(e2,0);
            halfedge_descriptor e21 = tmesh.halfedge(e2,1);
            face_descriptor f20 = tmesh.face(e20);
            face_descriptor f21 = tmesh.face(e21);

            bool same_face_edges = false;
            if ((f10==f20)||(f10==f21)||(f11==f20)||(f11==f21)) {
                same_face_edges = true;
                same_face += 1;
            }

            // Vector n20 = fnormals[f20];
            // Vector n21 = fnormals[f21];

            Point A2 = tmesh.point(source(e2,tmesh));
            Point B2 = tmesh.point(target(e2,tmesh));
            Segment SAB2(A2,B2);
            // Vector vAB2 = B2 - A2;
            Point M2 = A2+(B2-A2)/2;

            Vector vA2 = A2 - M1;
            Vector vB2 = B2 - M1;
            if (!same_face_edges) {
                if ((vA2 * n10 < -1e-3) && (vA2 * n11 < -1e-3) && (vB2 * n10 < -1e-3) && (vB2 * n11 < -1e-3)) {
                    not_visible += 1;
                    continue;
                } 

                if ((vA2 * n10 > 0) && (vA2 * n11 > 0) && (vB2 * n10 > 0) && (vB2 * n11 > 0)) {
                    fully_visible += 1;
                }
            }
            Vector vM = M2-M1;
            Point MM = M1 + vM/2;
            double lM = std::sqrt(vM.squared_length());
            Point M1i = MM - vM * (lM-2e-3)/(2*lM);
            Point M2i = MM + vM * (lM-2e-3)/(2*lM);
            if (!same_face_edges) {
                if ((inside(M1i) ==  CGAL::ON_BOUNDED_SIDE) || 
                        (inside(M2i) ==  CGAL::ON_BOUNDED_SIDE)) {
                    continue;
                }
                segment[0] = M1i;
                segment[1] = M2i;
                if (CGAL::Polygon_mesh_processing::do_intersect(tmesh,segment)) {
                    // std::cout << "Segment [M1i,M2i] is intersecting the mesh" << std::endl;
                    continue;
                }
            }

            Segment SM(M1,M2);
            Segment SA1A2(A1,A2);
            Segment SA1B2(A1,B2);
            Segment SB1A2(B1,A2);
            Segment SB1B2(B1,B2);
            Segment Stest[5] = {SM,SA1A2,SA1B2,SB1A2,SB1B2};

            // Calcule de la distance max
            float dmm,dmax=0;
            dmm = sqrt(CGAL::squared_distance(M1,M2));
            // float dmin = sqrt(CGAL::squared_distance(SAB1,SAB2));
            for (int i=0;i<5;i++) {
                dmax = std::max<float>(Stest[i].squared_length(),dmax);
            }
            dmax = sqrt(dmax);
#if 0
            epm[e1][e2] = DistancePair(dmin,dmax);
            epm[e2][e1] = DistancePair(dmin,dmax);
#else
            epm[e1][e2] = DistancePair(dmm,dmm);
            epm[e2][e1] = DistancePair(dmm,dmm);
#endif

#if 0
            double l = std::sqrt(vM.squared_length());
            double ratio = (l-1e-3)/(2*l);
            Segment SM(MM + (M1-MM)*ratio,MM + (M2-MM)*ratio);
            Segment SA1A2(MM + (A1-MM)*ratio,MM + (A2-MM)*ratio);
            Segment SA1B2(MM + (A1-MM)*ratio,MM + (B2-MM)*ratio);
            Segment SB1A2(MM + (B1-MM)*ratio,MM + (A2-MM)*ratio);
            Segment SB1B2(MM + (B1-MM)*ratio,MM + (B2-MM)*ratio);
            Segment Stest[5] = {SM,SA1A2,SA1B2,SB1A2,SB1B2};




            // Facet f1(A1,B1,A2);
            // if (f1.is_degenerate()) {
            //     continue;
            // }

#if 0
            if ((CGAL::cross_product(B1-A1,A2-A1) * (B2-A1)) > 0) {
                tetra.add_face(a1,b1,a2);
                tetra.add_face(a1,a2,b2);
                tetra.add_face(a1,b2,b1);
                tetra.add_face(b1,b2,a2);
            } else {
                tetra.add_face(a1,a2,b1);
                tetra.add_face(a1,b2,a2);
                tetra.add_face(a1,b1,b2);
                tetra.add_face(b1,a2,b2);
            }
#endif


            for (int i=0;i<5;i++) {
                //Facet_intersection intersection = tree.any_intersection(f1);
                if (tree.do_intersect(Stest[i])) {
                    intersecting += 1;
                    break;
                }
            }
#endif
        }
    }
    FILE *fie=fopen("interedges","w");
    for (PointEdgeMap::const_iterator eit=startmap.begin();eit!=startmap.end();eit++) {
        edge_descriptor e2 = eit->first;
        Point A2 = tmesh.point(source(e2,tmesh));
        Point B2 = tmesh.point(target(e2,tmesh));
        Point M2 = A2+(B2-A2)/2;
        const DistancePair & dp = eit->second;
        fprintf(fie,"%.3f %.3f %.3f %.3f %.3f\n%.3f %.3f %.3f %.3f %.3f\n\n\n",
                start.x(),start.y(),start.z(),dp.first,dp.second,
                M2.x(),M2.y(),M2.z(),dp.first,dp.second);
    }
    for (PointEdgeMap::const_iterator eit=endmap.begin();eit!=endmap.end();eit++) {
        edge_descriptor e2 = eit->first;
        Point A2 = tmesh.point(source(e2,tmesh));
        Point B2 = tmesh.point(target(e2,tmesh));
        Point M2 = A2+(B2-A2)/2;
        const DistancePair & dp = eit->second;
        fprintf(fie,"%.3f %.3f %.3f %.3f %.3f\n%.3f %.3f %.3f %.3f %.3f\n\n\n",
                end.x(),end.y(),end.z(),dp.first,dp.second,
                M2.x(),M2.y(),M2.z(),dp.first,dp.second);
    }
    for (EdgeGraph::const_iterator it=epm.begin();it!=epm.end();it++) {
        edge_descriptor e1 = it->first;
        Point A1 = tmesh.point(source(e1,tmesh));
        Point B1 = tmesh.point(target(e1,tmesh));
        Point M1 = A1+(B1-A1)/2;
        for (PointEdgeMap::const_iterator eit=it->second.begin();eit!=it->second.end();eit++) {
            edge_descriptor e2 = eit->first;
            Point A2 = tmesh.point(source(e2,tmesh));
            Point B2 = tmesh.point(target(e2,tmesh));
            Point M2 = A2+(B2-A2)/2;
            const DistancePair & dp = eit->second;
            fprintf(fie,"%.3f %.3f %.3f %.3f %.3f\n%.3f %.3f %.3f %.3f %.3f\n\n\n",
                    M1.x(),M1.y(),M1.z(),dp.first,dp.second,
                    M2.x(),M2.y(),M2.z(),dp.first,dp.second);
        }
    }
    fclose(fie);

    size_t epmsize = 0;
    for (EdgeGraph::const_iterator it=epm.begin();it!=epm.end();it++) {
        epmsize += it->second.size();
    }
    std::cout << "Number of edge pairs: " << edge_pairs << std::endl;
    std::cout << "EPM size: " << epmsize << std::endl;
    std::cout << "Number of same_face edge pairs: " << same_face << " (" << (100*same_face)/edge_pairs << "%)" << std::endl;
    std::cout << "Number of non-visible edge pairs: " << not_visible << " (" << (100*not_visible)/edge_pairs << "%)" << std::endl;
    std::cout << "Number of fully-visible edge pairs: " << fully_visible << " (" << (100*fully_visible)/edge_pairs << "%)" << std::endl;
    std::cout << "Number of intersecting edge pairs: " << intersecting << " (" << (100*intersecting)/edge_pairs << "%)" << std::endl;
    size_t remaining_cases = edge_pairs - not_visible - fully_visible;
    std::cout << "Remaining cases: " << remaining_cases << " (" << (100*remaining_cases)/edge_pairs << "%)" << std::endl;
    std::vector<EdgePath> edge_paths; 
    shortestPaths(startmap, epm, endmap, edge_paths, 1e-2, tmesh);

    std::vector<PolyLine> point_paths;
    for (size_t i=0;i<edge_paths.size();i++) {
        PolyLine pl;
        pl.push_back(start);
        char fname[128]="";
        sprintf(fname,"edges%03d",int(i));
        FILE * fp = fopen(fname,"w");
        for (size_t j=0;j<edge_paths[i].size();j++) {
            edge_descriptor e = edge_paths[i][j];
            Point A1 = tmesh.point(source(e,tmesh));
            Point B1 = tmesh.point(target(e,tmesh));
            Point M1 = A1+(B1-A1)/2;
            fprintf(fp,"%e %e %e\n%e %e %e\n%e %e %e\n\n\n",
                    A1.x(),A1.y(),A1.z(),
                    M1.x(),M1.y(),M1.z(),
                    B1.x(),B1.y(),B1.z());
            pl.push_back(M1);
        }
        fclose(fp);
        pl.push_back(end);
        point_paths.push_back(pl);
    }

    FILE * fp = fopen("paths","w");
    for (size_t i=0;i<point_paths.size();i++) {
        for (size_t j=0;j<point_paths[i].size();j++) {
            const Point & P = point_paths[i][j];
            Point A,B;
            edge_descriptor e = edge_paths[i][j-1];
            if (j==0) {
                A = B = start;
            } else if (j==point_paths[i].size()-1) {
                A = B = end;
            } else {
                A = tmesh.point(source(e,tmesh));
                B = tmesh.point(target(e,tmesh));
            }

            fprintf(fp,"%d %e %e %e %e %e %e %e %e %e\n",int(i),
                    P.x(),P.y(),P.z(),
                    A.x(),A.y(),A.z(),
                    B.x(),B.y(),B.z());
        }
        fprintf(fp,"\n\n");
    }
    fclose(fp);

    return EXIT_SUCCESS;
}
