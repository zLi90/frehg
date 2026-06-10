#ifndef FREHG2_BC_POLYGON_HPP
#define FREHG2_BC_POLYGON_HPP

#include "core/Grid.hpp"

#include <string>
#include <vector>

namespace frehg2 {

struct Point2 {
    real x = 0.0;
    real y = 0.0;
};

class Polygon {
public:
    Polygon() = default;
    explicit Polygon(std::vector<Point2> vertices);

    static Polygon readFromFile(const std::string& filename);

    const std::vector<Point2>& vertices() const noexcept;
    bool isInside(real x, real y) const;
    bool overlaps(const Polygon& other) const;

    std::vector<index_t> selectSurfaceCells(const Grid& grid) const;
    std::vector<index_t> selectGroundwaterCells(const Grid& grid) const;

private:
    std::vector<Point2> vertices_;

    void validate() const;
};

}  // namespace frehg2

#endif  // FREHG2_BC_POLYGON_HPP
