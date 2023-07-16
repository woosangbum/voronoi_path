#include "voronoi.h"
#include <iostream>

Voronoi_path::Voronoi_path(std::vector<Pos> w , int un) : waypoints(w), ugv_num(un) {
    double spacing = 1;  // k-means 생성할 점들의 간격
    int maxIterations = 100;  // 최대 반복 횟수
    std::vector<Pos> filledPoints = fillPolygon(waypoints, spacing);
    std::vector<Pos> centroids = kmeans(filledPoints, ugv_num, maxIterations);
    VD vd;
    std::vector<Site_2> points;
    for (const auto& point : centroids) {
        std::cout << "Centroid: (" << point.x << ", " << point.y << ")" << std::endl;
        points.push_back(Site_2(point.x, point.y));
    }
    points.push_back(Site_2(7000,7000));
    points.push_back(Site_2(-7000,-7000));
    points.push_back(Site_2(7000,-7000));
    points.push_back(Site_2(-7000,7000));
    
    for(const auto& point : centroids){
        points.push_back(Point_2(point.x, point.y));
    }

    for (const auto& point : points) {
        vd.insert(point);
    }

    std::vector<Pos> vertexs;
    // Iterate over Voronoi cells and extract vertices
    for (VD::Vertex_iterator vit = vd.vertices_begin(); vit != vd.vertices_end(); ++vit)
    {
        // Access the coordinates of the vertex
        double x = vit->point().x();
        double y = vit->point().y();

        vertexs.push_back(Pos(x, y));
    }

    for(const auto& vertex : vertexs){
        std::cout << "Vertex: (" << vertex.x << ", " << vertex.y << ")" << std::endl;
    }

    std::vector<std::vector<Pos>> voronoi_polygons;
    for (Face_iterator fit = vd.faces_begin(); fit != vd.faces_end(); ++fit)
    {
        if (fit->is_unbounded())
            continue;

        Ccb_halfedge_circulator circ = fit->ccb();
        std::vector<K::Point_2> face_points;

        do
        {
            if (!circ->has_source())
                break;
            face_points.push_back(circ->source()->point());
        } while (++circ != fit->ccb());

        // 분할된 영역 출력 (예: 좌표 출력)
        std::vector<Pos> voronoi_poly;
        for (const auto& p : face_points)
        {
            voronoi_poly.push_back(Pos(p.x(), p.y()));
        }
        // std::sort(voronoi_poly.begin(), voronoi_poly.end(), comparePoints);
        sortCounterClockwise(voronoi_poly);
        voronoi_polygons.push_back(voronoi_poly);
    }

    for (const auto& voronoi_polygon : voronoi_polygons){
        std::cout << "Region: ";
        for (const auto& p : voronoi_polygon){
            std::cout << "[" << p.x << ", " << p.y << "], ";
        }
        std::cout << std::endl;
    }

    // ===================================================


    BoostPolygon poly1, poly2;
    for(const auto& point : waypoints){
        bg::append(poly1.outer(), BoostPoint(point.x, point.y));
    }
    bg::correct(poly1);


    // 최종 경로 뽑아내는 곳
    for (const auto& voronoi_polygon : voronoi_polygons){
        poly2.clear();
        for (const auto& p : voronoi_polygon){
            bg::append(poly2.outer(), BoostPoint(p.x, p.y));
        }
        bg::correct(poly2);
        BoostMultiPolygon intersection;
        bg::intersection(poly1, poly2, intersection);

        std::cout << "Intersection polygons:" << std::endl;

        std::vector<Pos> car_path;
        for (const auto& polygon : intersection) {
            car_path.clear();
            std::cout << "Vertices: ";
            for (const auto& point : polygon.outer()) {
                std::cout << "[" << bg::get<0>(point) << ", " << bg::get<1>(point) << "], ";
                car_path.push_back(Pos(bg::get<0>(point),bg::get<1>(point)));
            }
            std::cout << std::endl;
        }
        car_paths.push_back(car_path);
    }

    // 벡터의 맨 앞에 있는 요소를 맨 뒤에 추가합니다.
    for (auto& path : car_paths) {
        if (!path.empty()) {
            const Pos& frontPoint = path.front();
            path.push_back(frontPoint);
        }
    }

    CGAL::draw(vd);
}

double Voronoi_path::crossProduct(const Pos& p1, const Pos& p2, const Pos& p3) {
    double x1 = p2.x - p1.x;
    double y1 = p2.y - p1.y;
    double x2 = p3.x - p2.x;
    double y2 = p3.y - p2.y;
    return x1 * y2 - x2 * y1;
}

// 볼록 껍질의 점들을 반시계 방향으로 정렬하는 함수
// 두 점 사이의 거리를 계산하는 함수
double Voronoi_path::distance(const Pos& p1, const Pos& p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return std::sqrt(dx * dx + dy * dy);
}
void Voronoi_path::sortCounterClockwise(std::vector<Pos>& points) {
    // 가장 아래에 있는 점을 찾습니다.
    Pos lowest = points[0];
    for (const auto& point : points) {
        if (point.y < lowest.y || (point.y == lowest.y && point.x < lowest.x)) {
            lowest = point;
        }
    }
    
    // 가장 아래에 있는 점을 기준으로 점들을 정렬합니다.
    auto compareAngles = [this, &lowest](const Pos& p1, const Pos& p2) {
        double angle1 = std::atan2(p1.y - lowest.y, p1.x - lowest.x);
        double angle2 = std::atan2(p2.y - lowest.y, p2.x - lowest.x);
        
        if (angle1 < angle2) {
            return true;
        } else if (angle1 == angle2) {
            double dist1 = distance(lowest, p1);
            double dist2 = distance(lowest, p2);
            return dist1 < dist2;
        } else {
            return false;
        }
    };
    
    std::sort(points.begin(), points.end(), compareAngles);
}


// K-Means 알고리즘을 수행하여 중심점을 찾는 함수
std::vector<Pos> Voronoi_path::kmeans(const std::vector<Pos>& points, int k, int maxIterations) {
    // 랜덤으로 초기 중심점 선택
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, points.size() - 1);
    std::vector<Pos> centroids;
    for (int i = 0; i < k; ++i) {
        centroids.push_back(points[dist(gen)]);
    }

    // 할당 단계와 갱신 단계를 반복하여 수행
    for (int iter = 0; iter < maxIterations; ++iter) {
        // 할당 단계: 각 점을 가장 가까운 중심점에 할당
        std::vector<std::vector<Pos>> clusters(k);
        for (const auto& point : points) {
            int closestCentroid = 0;
            double minDistance = std::numeric_limits<double>::max();
            for (int i = 0; i < k; ++i) {
                double d = distance(point, centroids[i]);
                if (d < minDistance) {
                    closestCentroid = i;
                    minDistance = d;
                }
            }
            clusters[closestCentroid].push_back(point);
        }

        // 갱신 단계: 중심점을 클러스터의 평균 위치로 업데이트
        for (int i = 0; i < k; ++i) {
            double sumX = 0.0;
            double sumY = 0.0;
            for (const auto& point : clusters[i]) {
                sumX += point.x;
                sumY += point.y;
            }
            if (!clusters[i].empty()) {
                centroids[i].x = sumX / clusters[i].size();
                centroids[i].y = sumY / clusters[i].size();
            }
        }
    }

    return centroids;
}

// 다각형 내부에 점이 있는지 확인하는 함수
bool Voronoi_path::isPointInsidePolygon(const std::vector<Pos>& polygon, const Pos& point) {
    int count = 0;
    int n = polygon.size();

    for (int i = 0; i < n; ++i) {
        const Pos& p1 = polygon[i];
        const Pos& p2 = polygon[(i + 1) % n];

        if (((p1.y <= point.y && point.y < p2.y) || (p2.y <= point.y && point.y < p1.y)) &&
            (point.x < (p2.x - p1.x) * (point.y - p1.y) / (p2.y - p1.y) + p1.x)) {
            ++count;
        }
    }

    return count % 2 == 1;
}


// 다각형 내부를 채우는 함수
std::vector<Pos> Voronoi_path::fillPolygon(const std::vector<Pos>& polygon, double spacing) {
    std::vector<Pos> filledPoints;

    // 다각형 경계 좌표 얻기
    std::vector<double> xCoords, yCoords;
    for (const auto& point : polygon) {
        xCoords.push_back(point.x);
        yCoords.push_back(point.y);
    }
    double minX = *std::min_element(xCoords.begin(), xCoords.end());
    double minY = *std::min_element(yCoords.begin(), yCoords.end());
    double maxX = *std::max_element(xCoords.begin(), xCoords.end());
    double maxY = *std::max_element(yCoords.begin(), yCoords.end());

    // 균일한 간격으로 점들 생성하여 다각형 내부에 위치하는지 확인
    for (double y = minY; y <= maxY; y += spacing) {
        for (double x = minX; x <= maxX; x += spacing) {
            Pos point{x, y};
            if (isPointInsidePolygon(polygon, point)) {
                filledPoints.push_back(point);
            }
        }
    }

    return filledPoints;
}