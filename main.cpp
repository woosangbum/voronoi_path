#include "voronoi.h"
#include <iostream>

int main() {
    int ugv_num = 10;  // 클러스터 개수 설정(UGV 개수)
    std::vector<Pos> waypoints = {{269.823, 500}, {450, 312.809}, {438.129, 105.814},
                                    {250.350, 100}, {50, 281.058}, {58.649, 484.560}};
    Voronoi_path vp(waypoints, ugv_num);
    
    for (const auto& path : vp.car_paths) {
        for (const auto& point : path) {
            std::cout << "[" << point.x << ", " << point.y << "], ";
        }
        std::cout << std::endl;
    }
}