#include "bc/Polygon.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace frehg2 {

namespace {

constexpr real POINT_TOL = 1.0e-12;

bool pointOnSegment(const Point2& point, const Point2& a, const Point2& b)
{
    const real length_sq = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y);
    if (length_sq <= POINT_TOL) {
        return std::abs(point.x - a.x) <= POINT_TOL && std::abs(point.y - a.y) <= POINT_TOL;
    }
    const real cross = (point.y - a.y) * (b.x - a.x) - (point.x - a.x) * (b.y - a.y);
    if (std::abs(cross) > POINT_TOL) {
        return false;
    }
    const real dot = (point.x - a.x) * (b.x - a.x) + (point.y - a.y) * (b.y - a.y);
    if (dot < -POINT_TOL) {
        return false;
    }
    return dot <= length_sq + POINT_TOL;
}

std::string stripComment(const std::string& line)
{
    const auto comment = line.find('#');
    return comment == std::string::npos ? line : line.substr(0, comment);
}

}  // namespace

Polygon::Polygon(std::vector<Point2> vertices)
    : vertices_(std::move(vertices))
{
    validate();
}

Polygon Polygon::readFromFile(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open polygon file: " + filename);
    }

    std::vector<std::vector<real>> rows;
    std::string line;
    while (std::getline(input, line)) {
        std::istringstream stream(stripComment(line));
        std::vector<real> values;
        real value = 0.0;
        while (stream >> value) {
            values.push_back(value);
        }
        if (!values.empty()) {
            rows.push_back(std::move(values));
        }
    }
    if (rows.empty()) {
        throw std::runtime_error("polygon file contains no vertices: " + filename);
    }

    std::size_t first_vertex = 0;
    if (rows.front().size() == 1) {
        const auto expected = static_cast<std::size_t>(rows.front()[0]);
        if (expected + 1 != rows.size()) {
            throw std::runtime_error("polygon vertex count does not match file rows: " + filename);
        }
        first_vertex = 1;
    }

    std::vector<Point2> vertices;
    for (std::size_t row = first_vertex; row < rows.size(); ++row) {
        if (rows[row].size() < 2) {
            throw std::runtime_error("polygon vertex row must contain x and y: " + filename);
        }
        vertices.push_back(Point2{rows[row][0], rows[row][1]});
    }
    return Polygon(std::move(vertices));
}

const std::vector<Point2>& Polygon::vertices() const noexcept
{
    return vertices_;
}

bool Polygon::isInside(real x, real y) const
{
    const Point2 point{x, y};
    bool inside = false;
    for (std::size_t i = 0, j = vertices_.size() - 1; i < vertices_.size(); j = i++) {
        const Point2& a = vertices_[i];
        const Point2& b = vertices_[j];
        if (pointOnSegment(point, a, b)) {
            return true;
        }
        const bool crosses = ((a.y > y) != (b.y > y)) &&
                             (x < (b.x - a.x) * (y - a.y) / (b.y - a.y) + a.x);
        if (crosses) {
            inside = !inside;
        }
    }
    return inside;
}

bool Polygon::overlaps(const Polygon& other) const
{
    return std::any_of(vertices_.begin(), vertices_.end(), [&other](const Point2& point) {
               return other.isInside(point.x, point.y);
           }) ||
           std::any_of(other.vertices_.begin(), other.vertices_.end(), [this](const Point2& point) {
               return isInside(point.x, point.y);
           });
}

std::vector<index_t> Polygon::selectSurfaceCells(const Grid& grid) const
{
    std::vector<index_t> selected;
    for (int j = 0; j < grid.ny(); ++j) {
        for (int i = 0; i < grid.nx(); ++i) {
            const real x = (static_cast<real>(i) + 0.5) * grid.dx();
            const real y = (static_cast<real>(j) + 0.5) * grid.dy();
            if (isInside(x, y)) {
                selected.push_back(grid.getSurfaceIndex(i, j));
            }
        }
    }
    return selected;
}

std::vector<index_t> Polygon::selectGroundwaterCells(const Grid& grid) const
{
    std::vector<index_t> selected;
    for (int j = 0; j < grid.ny(); ++j) {
        for (int i = 0; i < grid.nx(); ++i) {
            const real x = (static_cast<real>(i) + 0.5) * grid.dx();
            const real y = (static_cast<real>(j) + 0.5) * grid.dy();
            if (!isInside(x, y)) {
                continue;
            }
            for (int k = 0; k < grid.nz(); ++k) {
                selected.push_back(grid.getIndex(i, j, k));
            }
        }
    }
    return selected;
}

void Polygon::validate() const
{
    if (vertices_.size() < 3) {
        throw std::invalid_argument("polygon must contain at least three vertices");
    }
}

}  // namespace frehg2
