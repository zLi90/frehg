#include "core/Grid.hpp"

#include <stdexcept>

namespace frehg2 {

Grid::Grid(GridSpec spec, int mpi_nx, int mpi_ny, int rank)
    : spec_(spec),
      mpi_nx_(mpi_nx),
      mpi_ny_(mpi_ny),
      rank_(rank)
{
    validate();

    dz_layers_.resize(static_cast<std::size_t>(spec_.nz), spec_.dz);
    real layer_dz = spec_.dz;
    for (auto& value : dz_layers_) {
        value = layer_dz;
        layer_dz *= spec_.dz_multiplier;
    }
}

const GridSpec& Grid::spec() const noexcept
{
    return spec_;
}

int Grid::nx() const noexcept
{
    return spec_.nx;
}

int Grid::ny() const noexcept
{
    return spec_.ny;
}

int Grid::nz() const noexcept
{
    return spec_.nz;
}

real Grid::dx() const noexcept
{
    return spec_.dx;
}

real Grid::dy() const noexcept
{
    return spec_.dy;
}

real Grid::dz(int k) const
{
    if (k < 0 || k >= nz()) {
        throw std::out_of_range("groundwater layer index is out of range");
    }
    return dz_layers_[static_cast<std::size_t>(k)];
}

int Grid::mpiNx() const noexcept
{
    return mpi_nx_;
}

int Grid::mpiNy() const noexcept
{
    return mpi_ny_;
}

int Grid::rank() const noexcept
{
    return rank_;
}

int Grid::rankX() const noexcept
{
    return rank_ % mpi_nx_;
}

int Grid::rankY() const noexcept
{
    return rank_ / mpi_nx_;
}

index_t Grid::nSurfaceCell() const noexcept
{
    return static_cast<index_t>(nx()) * static_cast<index_t>(ny());
}

index_t Grid::nSurfaceCellMem() const noexcept
{
    return static_cast<index_t>(nx() + 2) * static_cast<index_t>(ny() + 2);
}

index_t Grid::nCell() const noexcept
{
    return nSurfaceCell() * static_cast<index_t>(nz());
}

index_t Grid::nCellMem() const noexcept
{
    return nSurfaceCellMem() * static_cast<index_t>(nz() + 2);
}

index_t Grid::getIndex(int i, int j, int k) const
{
    if (!isInterior(i, j, k)) {
        throw std::out_of_range("grid index is outside the interior domain");
    }

    return (static_cast<index_t>(j) * static_cast<index_t>(nx()) + static_cast<index_t>(i)) *
               static_cast<index_t>(nz()) +
           static_cast<index_t>(k);
}

index_t Grid::getSurfaceIndex(int i, int j) const
{
    if (i < 0 || i >= nx() || j < 0 || j >= ny()) {
        throw std::out_of_range("surface index is outside the interior domain");
    }

    return static_cast<index_t>(i) + static_cast<index_t>(j) * static_cast<index_t>(nx());
}

Index3 Grid::getLogicalIndex(index_t index) const
{
    if (index >= nCell()) {
        throw std::out_of_range("flat index is outside the interior domain");
    }

    const auto column = index / static_cast<index_t>(nz());
    return Index3{
        static_cast<int>(column % static_cast<index_t>(nx())),
        static_cast<int>(column / static_cast<index_t>(nx())),
        static_cast<int>(index % static_cast<index_t>(nz())),
    };
}

std::array<int, 2> Grid::getSurfaceLogicalIndex(index_t index) const
{
    if (index >= nSurfaceCell()) {
        throw std::out_of_range("flat surface index is outside the interior domain");
    }

    return {
        static_cast<int>(index % static_cast<index_t>(nx())),
        static_cast<int>(index / static_cast<index_t>(nx())),
    };
}

index_t Grid::surfaceGhostIndex(Direction direction, int i, int j) const
{
    const auto n2ci = nSurfaceCell();
    switch (direction) {
    case Direction::XP:
        if (i != nx() - 1 || j < 0 || j >= ny()) {
            throw std::out_of_range("invalid x+ surface ghost request");
        }
        return n2ci + static_cast<index_t>(2 * nx() + j);
    case Direction::XM:
        if (i != 0 || j < 0 || j >= ny()) {
            throw std::out_of_range("invalid x- surface ghost request");
        }
        return n2ci + static_cast<index_t>(2 * nx() + ny() + j);
    case Direction::YP:
        if (j != ny() - 1 || i < 0 || i >= nx()) {
            throw std::out_of_range("invalid y+ surface ghost request");
        }
        return n2ci + static_cast<index_t>(i);
    case Direction::YM:
        if (j != 0 || i < 0 || i >= nx()) {
            throw std::out_of_range("invalid y- surface ghost request");
        }
        return n2ci + static_cast<index_t>(nx() + i);
    case Direction::ZP:
    case Direction::ZM:
        throw std::invalid_argument("surface grid has no vertical ghost direction");
    }

    throw std::invalid_argument("unknown surface ghost direction");
}

index_t Grid::groundwaterGhostIndex(Direction direction, int i, int j, int k) const
{
    const auto n2ct = nSurfaceCellMem();
    const auto n3ci = nCell();
    const auto top2d = static_cast<index_t>(j) * static_cast<index_t>(nx()) +
                       static_cast<index_t>(i);

    switch (direction) {
    case Direction::XP:
        if (i != nx() - 1 || j < 0 || j >= ny() || k < 0 || k >= nz()) {
            throw std::out_of_range("invalid x+ groundwater ghost request");
        }
        return n3ci + static_cast<index_t>(2 * nx() * nz() + j * nz() + k);
    case Direction::XM:
        if (i != 0 || j < 0 || j >= ny() || k < 0 || k >= nz()) {
            throw std::out_of_range("invalid x- groundwater ghost request");
        }
        return n3ci + static_cast<index_t>(2 * nx() * nz() + ny() * nz() + j * nz() + k);
    case Direction::YP:
        if (j != ny() - 1 || i < 0 || i >= nx() || k < 0 || k >= nz()) {
            throw std::out_of_range("invalid y+ groundwater ghost request");
        }
        return n3ci + static_cast<index_t>(i * nz() + k);
    case Direction::YM:
        if (j != 0 || i < 0 || i >= nx() || k < 0 || k >= nz()) {
            throw std::out_of_range("invalid y- groundwater ghost request");
        }
        return n3ci + static_cast<index_t>(nx() * nz() + i * nz() + k);
    case Direction::ZP:
        if (k != nz() - 1 || i < 0 || i >= nx() || j < 0 || j >= ny()) {
            throw std::out_of_range("invalid z+ groundwater ghost request");
        }
        return n2ct * static_cast<index_t>(nz()) + top2d;
    case Direction::ZM:
        if (k != 0 || i < 0 || i >= nx() || j < 0 || j >= ny()) {
            throw std::out_of_range("invalid z- groundwater ghost request");
        }
        return n2ct * static_cast<index_t>(nz() + 1) + top2d;
    }

    throw std::invalid_argument("unknown groundwater ghost direction");
}

bool Grid::isInterior(int i, int j, int k) const noexcept
{
    return i >= 0 && i < nx() && j >= 0 && j < ny() && k >= 0 && k < nz();
}

bool Grid::isActive(int i, int j, int k) const noexcept
{
    return isInterior(i, j, k);
}

bool Grid::isGlobalBoundary(Direction direction) const noexcept
{
    switch (direction) {
    case Direction::XM:
        return rankX() == 0;
    case Direction::XP:
        return rankX() == mpi_nx_ - 1;
    case Direction::YM:
        return rankY() == 0;
    case Direction::YP:
        return rankY() == mpi_ny_ - 1;
    case Direction::ZP:
    case Direction::ZM:
        return true;
    }

    return true;
}

Index3 Grid::localToGlobal(int i, int j, int k) const
{
    if (!isInterior(i, j, k)) {
        throw std::out_of_range("local index is outside the interior domain");
    }

    return Index3{rankX() * nx() + i, rankY() * ny() + j, k};
}

Index3 Grid::globalToLocal(int i, int j, int k) const
{
    const Index3 local{i - rankX() * nx(), j - rankY() * ny(), k};
    if (!isInterior(local.i, local.j, local.k)) {
        throw std::out_of_range("global index does not belong to this rank");
    }

    return local;
}

const std::vector<real>& Grid::dzLayers() const noexcept
{
    return dz_layers_;
}

void Grid::validate() const
{
    if (spec_.nx <= 0 || spec_.ny <= 0 || spec_.nz <= 0) {
        throw std::invalid_argument("grid dimensions must be positive");
    }
    if (spec_.dx <= 0.0 || spec_.dy <= 0.0 || spec_.dz <= 0.0) {
        throw std::invalid_argument("grid spacing must be positive");
    }
    if (spec_.dz_multiplier <= 0.0) {
        throw std::invalid_argument("vertical spacing multiplier must be positive");
    }
    if (mpi_nx_ <= 0 || mpi_ny_ <= 0) {
        throw std::invalid_argument("MPI decomposition dimensions must be positive");
    }
    if (rank_ < 0 || rank_ >= mpi_nx_ * mpi_ny_) {
        throw std::invalid_argument("rank is outside the MPI decomposition");
    }
}

}  // namespace frehg2
