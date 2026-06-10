#ifndef FREHG2_CORE_GRID_HPP
#define FREHG2_CORE_GRID_HPP

#include "core/types.hpp"

#include <array>
#include <vector>

namespace frehg2 {

enum class Direction {
    XP,
    XM,
    YP,
    YM,
    ZP,
    ZM
};

struct Index3 {
    int i = 0;
    int j = 0;
    int k = 0;
};

class Grid {
public:
    static constexpr int HALO_WIDTH = 1;

    explicit Grid(GridSpec spec, int mpi_nx = 1, int mpi_ny = 1, int rank = 0);

    const GridSpec& spec() const noexcept;

    int nx() const noexcept;
    int ny() const noexcept;
    int nz() const noexcept;
    real dx() const noexcept;
    real dy() const noexcept;
    real dz(int k = 0) const;

    int mpiNx() const noexcept;
    int mpiNy() const noexcept;
    int rank() const noexcept;
    int rankX() const noexcept;
    int rankY() const noexcept;

    index_t nCell() const noexcept;
    index_t nCellMem() const noexcept;
    index_t nSurfaceCell() const noexcept;
    index_t nSurfaceCellMem() const noexcept;

    index_t getIndex(int i, int j, int k = 0) const;
    index_t getSurfaceIndex(int i, int j) const;
    Index3 getLogicalIndex(index_t index) const;
    std::array<int, 2> getSurfaceLogicalIndex(index_t index) const;

    index_t surfaceGhostIndex(Direction direction, int i, int j) const;
    index_t groundwaterGhostIndex(Direction direction, int i, int j, int k) const;

    bool isInterior(int i, int j, int k = 0) const noexcept;
    bool isActive(int i, int j, int k = 0) const noexcept;
    bool isGlobalBoundary(Direction direction) const noexcept;

    Index3 localToGlobal(int i, int j, int k = 0) const;
    Index3 globalToLocal(int i, int j, int k = 0) const;

    const std::vector<real>& dzLayers() const noexcept;

private:
    GridSpec spec_;
    int mpi_nx_ = 1;
    int mpi_ny_ = 1;
    int rank_ = 0;
    std::vector<real> dz_layers_;

    void validate() const;
};

}  // namespace frehg2

#endif  // FREHG2_CORE_GRID_HPP
