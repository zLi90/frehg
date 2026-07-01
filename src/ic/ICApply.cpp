#include "ic/ICApply.hpp"

#include <cmath>
#include <stdexcept>

#include "bc/PolygonIndex.hpp"
#include "core/Grid.hpp"
#include "core/MpiComm.hpp"
#include "ic/FormulaParser.hpp"
#include "ic/RasterField.hpp"
#include "re/ReSolver.hpp"
#include "swe/SweSolver.hpp"

namespace frehg2 {

namespace {

RealArr1DHost sliceSurface2D(const RasterField& global, const MpiComm* mc, int gnx, int gny) {
  const int lnx = mc ? mc->localNx() : gnx;
  const int lny = mc ? mc->localNy() : gny;
  RealArr1DHost local("ic2d", static_cast<size_t>(lnx * lny));
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li) {
      const int gi = mc ? mc->i0() + li : li;
      const int gj = mc ? mc->j0() + lj : lj;
      local(static_cast<size_t>(li + lj * lnx)) = global.at(gi, gj, 0);
    }
  return local;
}

std::vector<double> sliceVolume3D(const RasterField& global, const MpiComm* mc, int gnx, int gny,
                                   int gnz) {
  const int lnx = mc ? mc->localNx() : gnx;
  const int lny = mc ? mc->localNy() : gny;
  std::vector<double> local(static_cast<size_t>(lnx * lny * gnz));
  for (int k = 0; k < gnz; ++k)
    for (int lj = 0; lj < lny; ++lj)
      for (int li = 0; li < lnx; ++li) {
        const int gi = mc ? mc->i0() + li : li;
        const int gj = mc ? mc->j0() + lj : lj;
        local[static_cast<size_t>(li + lj * lnx + k * lnx * lny)] = global.at(gi, gj, k);
      }
  return local;
}

RealArr1DHost fillSurfaceConstant(int lnx, int lny, real value) {
  RealArr1DHost f("ic_const", static_cast<size_t>(lnx * lny));
  for (size_t k = 0; k < f.extent(0); ++k) f(k) = value;
  return f;
}

RealArr1DHost fillSurfaceFormula(const ICFieldSpec& spec, const ICApplyContext& ctx, int lnx,
                                 int lny) {
  FormulaParser parser(spec.formula);
  RealArr1DHost f("ic_formula", static_cast<size_t>(lnx * lny));
  FormulaVars vars;
  vars.t = 0.0;
  vars.z = 0.0;
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li) {
      const int gi = ctx.mc ? ctx.mc->i0() + li : li;
      const int gj = ctx.mc ? ctx.mc->j0() + lj : lj;
      vars.x = ctx.x0 + (static_cast<real>(gi) + 0.5) * ctx.dx;
      vars.y = ctx.y0 + (static_cast<real>(gj) + 0.5) * ctx.dy;
      f(static_cast<size_t>(li + lj * lnx)) = parser.eval(vars);
    }
  return f;
}

RealArr1DHost fillSurfacePolygon(const ICFieldSpec& spec, const InitialConditionsConfig& ic,
                                 const PolygonIndex& index, const Grid& grid, int lnx,
                                 int lny) {
  RealArr1DHost f("ic_poly", static_cast<size_t>(lnx * lny));
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li) {
      real val = spec.polygon_default;
      const int surf_idx = grid.getSurfaceIndex(li, lj);
      const int p = index.lookup(surf_idx);
      if (p >= 0 && static_cast<size_t>(p) < ic.regions.size()) {
        val = ic.regions[static_cast<size_t>(p)].value;
        const auto it = spec.polygon_values.find(ic.regions[static_cast<size_t>(p)].polygon.name);
        if (it != spec.polygon_values.end()) val = it->second;
      }
      f(static_cast<size_t>(li + lj * lnx)) = val;
    }
  return f;
}

RealArr1DHost resolveSurface2D(const ICFieldSpec& spec, const InitialConditionsConfig& ic,
                               const ICApplyContext& ctx, const PolygonIndex* region_index,
                               const Grid& surf_grid) {
  const int lnx = ctx.mc ? ctx.mc->localNx() : ctx.gnx;
  const int lny = ctx.mc ? ctx.mc->localNy() : ctx.gny;
  switch (spec.kind) {
    case ICKind::Constant:
      return fillSurfaceConstant(lnx, lny, spec.constant_value);
    case ICKind::Raster: {
      const RasterField global = loadRasterField(*ctx.config, spec.file, spec.format, spec.dataset,
                                                 ctx.gnx, ctx.gny, 1);
      return sliceSurface2D(global, ctx.mc, ctx.gnx, ctx.gny);
    }
    case ICKind::Formula:
      return fillSurfaceFormula(spec, ctx, lnx, lny);
    case ICKind::Polygon:
      if (region_index == nullptr || region_index->empty())
        throw std::runtime_error("ICApply: polygon IC requires initial_conditions.regions");
      return fillSurfacePolygon(spec, ic, *region_index, surf_grid, lnx, lny);
    default:
      throw std::runtime_error("ICApply: unsupported surface IC kind");
  }
}

std::vector<double> resolveVolume3D(const ICFieldSpec& spec, const InitialConditionsConfig& ic,
                                    const ICApplyContext& ctx, const PolygonIndex* region_index,
                                    real column_default) {
  const int lnx = ctx.mc ? ctx.mc->localNx() : ctx.gnx;
  const int lny = ctx.mc ? ctx.mc->localNy() : ctx.gny;
  const int gnz = ctx.gnz;
  switch (spec.kind) {
    case ICKind::Constant: {
      std::vector<double> v(static_cast<size_t>(lnx * lny * gnz), column_default);
      return v;
    }
    case ICKind::Raster: {
      const RasterField global = loadRasterField(*ctx.config, spec.file, spec.format, spec.dataset,
                                                 ctx.gnx, ctx.gny, gnz);
      return sliceVolume3D(global, ctx.mc, ctx.gnx, ctx.gny, gnz);
    }
    case ICKind::Formula: {
      FormulaParser parser(spec.formula);
      std::vector<double> v(static_cast<size_t>(lnx * lny * gnz));
      FormulaVars vars;
      vars.t = 0.0;
      for (int k = 0; k < gnz; ++k)
        for (int lj = 0; lj < lny; ++lj)
          for (int li = 0; li < lnx; ++li) {
            const int gi = ctx.mc ? ctx.mc->i0() + li : li;
            const int gj = ctx.mc ? ctx.mc->j0() + lj : lj;
            vars.x = ctx.x0 + (static_cast<real>(gi) + 0.5) * ctx.dx;
            vars.y = ctx.y0 + (static_cast<real>(gj) + 0.5) * ctx.dy;
            vars.z = ctx.botz - (static_cast<real>(k) + 0.5) * ctx.dz;
            v[static_cast<size_t>(li + lj * lnx + k * lnx * lny)] = parser.eval(vars);
          }
      return v;
    }
    case ICKind::Polygon: {
      if (region_index == nullptr || region_index->empty())
        throw std::runtime_error("ICApply: polygon IC requires initial_conditions.regions");
      std::vector<double> v(static_cast<size_t>(lnx * lny * gnz));
      Grid g(lnx, lny, gnz, ctx.dx, ctx.dy, ctx.dz);
      for (int k = 0; k < gnz; ++k)
        for (int lj = 0; lj < lny; ++lj)
          for (int li = 0; li < lnx; ++li) {
            real val = spec.polygon_default;
            const int p = region_index->lookup(g.getSurfaceIndex(li, lj));
            if (p >= 0 && static_cast<size_t>(p) < ic.regions.size())
              val = ic.regions[static_cast<size_t>(p)].value;
            v[static_cast<size_t>(li + lj * lnx + k * lnx * lny)] = val;
          }
      return v;
    }
    default:
      throw std::runtime_error("ICApply: unsupported volume IC kind");
  }
}

}  // namespace

void buildICRegionIndex(const InitialConditionsConfig& ic, ICApplyContext ctx) {
  if (ic.regions.empty() || ctx.region_index == nullptr) return;
  std::vector<Polygon> polys;
  polys.reserve(ic.regions.size());
  for (const ICRegion& r : ic.regions) polys.push_back(r.polygon);
  const Grid& grid = ctx.mc ? Grid(ctx.mc->localNx(), ctx.mc->localNy(), 1, ctx.dx, ctx.dy, ctx.dz)
                            : Grid(ctx.gnx, ctx.gny, 1, ctx.dx, ctx.dy, ctx.dz);
  ctx.region_index->build(polys, grid, ctx.mc, ctx.x0, ctx.y0);
}

void applyInitialConditions(const InitialConditionsConfig& ic, ICApplyContext ctx, SweSolver* swe,
                            ReSolver* re, RealArr1DHost surface_conc,
                            RealArr1DHost subsurface_conc) {
  if (ic.use_restart) return;

  PolygonIndex local_index;
  if (!ic.regions.empty()) {
    ctx.region_index = &local_index;
    buildICRegionIndex(ic, ctx);
  }
  const PolygonIndex* region_index = ctx.region_index;
  const int lnx = ctx.mc ? ctx.mc->localNx() : ctx.gnx;
  const int lny = ctx.mc ? ctx.mc->localNy() : ctx.gny;
  const Grid surf_grid(lnx, lny, 1, ctx.dx, ctx.dy, ctx.dz);

  if (swe != nullptr) {
    const RealArr1DHost eta =
        resolveSurface2D(ic.surface_eta, ic, ctx, region_index, surf_grid);
    swe->setInitialEtaPhysical(eta);
    if (ic.surface_u.kind != ICKind::Constant || ic.surface_v.kind != ICKind::Constant ||
        ic.surface_u.constant_value != 0.0 || ic.surface_v.constant_value != 0.0) {
      const RealArr1DHost u =
          resolveSurface2D(ic.surface_u, ic, ctx, region_index, surf_grid);
      const RealArr1DHost v =
          resolveSurface2D(ic.surface_v, ic, ctx, region_index, surf_grid);
      swe->setInitialVelocityPhysical(u, v);
    } else {
      swe->setInitialVelocityConstant(0.0, 0.0);
    }
    swe->finalizeInitialState();
  }

  if (re != nullptr) {
    re->initializeGeometry();
    const auto isSpatial = [](ICKind k) {
      return k == ICKind::Raster || k == ICKind::Formula || k == ICKind::Polygon;
    };
    if (isSpatial(ic.groundwater_head.kind)) {
      const std::vector<double> h = resolveVolume3D(ic.groundwater_head, ic, ctx, region_index,
                                                    ic.groundwater_head.constant_value);
      re->setInitialHead(h);
    } else {
      const std::vector<double> wc = resolveVolume3D(ic.groundwater_wc, ic, ctx, region_index,
                                                     ic.groundwater_wc.constant_value);
      re->setInitialWaterContent(wc);
    }
    re->finalizeInitialState();
  }

  // Solute initial concentration (P16-completion). The Orchestrator passes the canonical conc
  // views when solute is enabled; honors the same constant/raster/formula/polygon kinds as the
  // flow ICs. Resolved owned arrays (li + lj*lnx [+ k*lnx*lny]) are scattered into the
  // halo-padded conc storage.
  if (surface_conc.extent(0) > 0) {
    const RealArr1DHost c = resolveSurface2D(ic.solute_surface, ic, ctx, region_index, surf_grid);
    for (int lj = 0; lj < lny; ++lj)
      for (int li = 0; li < lnx; ++li)
        surface_conc(static_cast<size_t>(surf_grid.getSurfaceIndex(li, lj))) =
            c(static_cast<size_t>(li + lj * lnx));
  }
  if (subsurface_conc.extent(0) > 0) {
    const Grid vol_grid(lnx, lny, ctx.gnz, ctx.dx, ctx.dy, ctx.dz);
    const std::vector<double> c = resolveVolume3D(ic.solute_subsurface, ic, ctx, region_index,
                                                  ic.solute_subsurface.constant_value);
    for (int k = 0; k < ctx.gnz; ++k)
      for (int lj = 0; lj < lny; ++lj)
        for (int li = 0; li < lnx; ++li)
          subsurface_conc(static_cast<size_t>(vol_grid.getIndex(li, lj, k))) =
              c[static_cast<size_t>(li + lj * lnx + k * lnx * lny)];
  }
}

}  // namespace frehg2
