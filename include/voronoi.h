#ifndef VORONOI_H
#define VORONOI_H

#define CGAL_USE_BASIC_VIEWER

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <tuple>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/draw_voronoi_diagram_2.h>
#include <vector>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Simple_cartesian.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT> AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT,AT,AP> VD;
typedef AT::Site_2 Site_2;
typedef VD::Locate_result Locate_result;	
typedef VD::Face_iterator Face_iterator;
typedef VD::Ccb_halfedge_circulator Ccb_halfedge_circulator;
typedef K::Point_2 Point_2;

namespace bg = boost::geometry;
typedef bg::model::d2::point_xy<double> BoostPoint;
typedef bg::model::polygon<BoostPoint> BoostPolygon;
typedef bg::model::multi_polygon<BoostPolygon> BoostMultiPolygon;

struct Pos {
    double x;
    double y;

    Pos(double x, double y) : x(x), y(y) {}
};

class Voronoi_path {
public:
    int ugv_num;
    std::vector<Pos> waypoints; // Input Data : waypoints
    std::vector<std::vector<Pos>> car_paths; // output Data
    // 생성자
    Voronoi_path(std::vector<Pos> w, int un);

    // 멤버 함수
    double crossProduct(const Pos& p1, const Pos& p2, const Pos& p3);
    double distance(const Pos& p1, const Pos& p2);
    void sortCounterClockwise(std::vector<Pos>& points);
    std::vector<Pos> kmeans(const std::vector<Pos>& points, int k, int maxIterations);
    bool isPointInsidePolygon(const std::vector<Pos>& polygon, const Pos& point);
    std::vector<Pos> fillPolygon(const std::vector<Pos>& polygon, double spacing);
    
};

#endif
