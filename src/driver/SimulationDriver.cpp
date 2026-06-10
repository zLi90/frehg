#include "driver/SimulationDriver.hpp"

#include "bc/BoundaryCondition.hpp"
#include "core/InitialCondition.hpp"
#include "core/Monitor.hpp"
#include "coupling/Coupling.hpp"
#include "io/Config.hpp"
#include "io/Hdf5Writer.hpp"
#include "re/ReSolver.hpp"
#include "swe/SweSolver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace frehg2 {

namespace {

using Clock = std::chrono::steady_clock;

double elapsedSeconds(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

struct RuntimeReport {
    std::string simulation_id;
    std::string module;
    double total_seconds = 0.0;
    std::vector<std::pair<std::string, double>> sections;
};

class ProgressReporter {
public:
    ProgressReporter(std::string simulation_id, real tend)
        : simulation_id_(std::move(simulation_id)),
          tend_(tend)
    {
        std::cout << "[" << simulation_id_ << "] progress 0% (t = 0 / " << tend_ << " s)\n";
    }

    void update(real time)
    {
        if (tend_ <= 0.0) {
            return;
        }
        const int percent =
            std::min(100, static_cast<int>(std::floor(100.0 * time / tend_ + 1.0e-9)));
        while (next_percent_ <= 100 && percent >= next_percent_) {
            std::cout << "[" << simulation_id_ << "] progress " << next_percent_
                      << "% (t = " << std::min(time, tend_) << " / " << tend_ << " s)\n";
            next_percent_ += 5;
        }
    }

private:
    std::string simulation_id_;
    real tend_ = 0.0;
    int next_percent_ = 5;
};

void writeRuntimeReport(const std::filesystem::path& output_dir, const RuntimeReport& report)
{
    const auto summary_path = output_dir / "simulation_summary.txt";
    std::ofstream summary(summary_path);
    summary << std::fixed << std::setprecision(6);
    summary << "simulation_id " << report.simulation_id << "\n";
    summary << "module " << report.module << "\n";
    summary << "total_runtime_seconds " << report.total_seconds << "\n";
    for (const auto& [name, seconds] : report.sections) {
        summary << name << "_seconds " << seconds << "\n";
    }

    std::cout << "[" << report.simulation_id << "] simulation complete\n";
    std::cout << "Runtime summary (seconds):\n";
    std::cout << "  total: " << std::fixed << std::setprecision(3) << report.total_seconds << "\n";
    for (const auto& [name, seconds] : report.sections) {
        std::cout << "  " << name << ": " << std::fixed << std::setprecision(3) << seconds
                  << "\n";
    }
    std::cout << "  report: " << summary_path << "\n";
}

std::vector<real> readColumnFile(const std::filesystem::path& path)
{
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("failed to open input file: " + path.string());
    }
    std::vector<real> values;
    real value = 0.0;
    while (input >> value) {
        values.push_back(value);
    }
    return values;
}

real legacySeriesValueAt(const std::vector<real>& values, real time)
{
    if (values.size() < 2) {
        return 0.0;
    }
    for (std::size_t i = 0; i + 3 < values.size(); i += 2) {
        const auto t0 = values[i];
        const auto v0 = values[i + 1];
        const auto t1 = values[i + 2];
        const auto v1 = values[i + 3];
        if (time <= t0) {
            return v0;
        }
        if (time <= t1) {
            const auto fraction = (time - t0) / (t1 - t0);
            return v0 + fraction * (v1 - v0);
        }
    }
    return values.back();
}

std::filesystem::path resolveBenchmarkInput(
    const std::filesystem::path& config_path,
    const std::string& legacy_input_dir,
    const std::string& filename)
{
    const auto from_config = config_path.parent_path() / legacy_input_dir / filename;
    if (std::filesystem::exists(from_config)) {
        return from_config;
    }

    const auto from_legacy = std::filesystem::current_path() / "legacy/benchmarks/b1-sw" /
                             legacy_input_dir / filename;
    if (std::filesystem::exists(from_legacy)) {
        return from_legacy;
    }

    const auto repo_root = config_path.parent_path().parent_path().parent_path();
    const auto from_repo_legacy = repo_root / "legacy/benchmarks/b1-sw" / legacy_input_dir / filename;
    if (std::filesystem::exists(from_repo_legacy)) {
        return from_repo_legacy;
    }

    return from_config;
}

std::filesystem::path resolveOutputPath(
    const Config& config,
    const std::filesystem::path& config_path)
{
    std::filesystem::path output_path = config.get<std::string>("output.filename");
    if (output_path.is_relative()) {
        output_path = config_path.parent_path() / output_path;
    }
    const auto output_dir = output_path.parent_path();
    if (!output_dir.empty()) {
        std::filesystem::create_directories(output_dir);
    }
    return output_path;
}

std::filesystem::path resolveInputPath(
    const std::filesystem::path& config_path,
    const std::string& filename)
{
    std::filesystem::path path(filename);
    if (path.is_relative()) {
        path = config_path.parent_path() / path;
    }
    return path;
}

InitialConditionSource initialConditionSource(const std::string& source)
{
    if (source == "constant") {
        return InitialConditionSource::CONSTANT;
    }
    if (source == "ascii_raster") {
        return InitialConditionSource::ASCII_RASTER;
    }
    if (source == "ascii_raster_3d") {
        return InitialConditionSource::ASCII_RASTER;
    }
    if (source == "hdf5_vector") {
        return InitialConditionSource::HDF5_VECTOR;
    }
    throw std::invalid_argument("unsupported initial condition source: " + source);
}

InitialConditionSpec initialConditionSpec(
    const Config& config,
    const std::filesystem::path& config_path,
    const std::string& path,
    const std::string& variable)
{
    InitialConditionSpec spec;
    spec.variable = variable;
    spec.source = initialConditionSource(config.get<std::string>(path + ".source"));
    spec.value = config.getOr<real>(path + ".value", 0.0);
    if (config.has(path + ".file")) {
        spec.filename = resolveInputPath(config_path, config.get<std::string>(path + ".file")).string();
    }
    spec.dataset = config.getOr<std::string>(path + ".dataset", "");
    return spec;
}

BoundaryConditionType boundaryConditionType(const std::string& value)
{
    if (value == "fixed_water_level") {
        return BoundaryConditionType::FIXED_WATER_LEVEL;
    }
    if (value == "fixed_flow_rate") {
        return BoundaryConditionType::FIXED_FLOW_RATE;
    }
    if (value == "free_outflow") {
        return BoundaryConditionType::FREE_OUTFLOW;
    }
    if (value == "zero_gradient") {
        return BoundaryConditionType::ZERO_GRADIENT;
    }
    if (value == "tidal_water_level") {
        return BoundaryConditionType::TIDAL_WATER_LEVEL;
    }
    if (value == "groundwater_head") {
        return BoundaryConditionType::GROUNDWATER_HEAD;
    }
    if (value == "root_water_uptake") {
        return BoundaryConditionType::ROOT_WATER_UPTAKE;
    }
    if (value == "gravity_drainage") {
        return BoundaryConditionType::GRAVITY_DRAINAGE;
    }
    throw std::invalid_argument("unsupported boundary condition type: " + value);
}

SourceSinkType sourceSinkType(const std::string& value)
{
    if (value == "rainfall") {
        return SourceSinkType::RAINFALL;
    }
    if (value == "evapotranspiration") {
        return SourceSinkType::EVAPOTRANSPIRATION;
    }
    if (value == "inflow") {
        return SourceSinkType::INFLOW;
    }
    if (value == "pumping") {
        return SourceSinkType::PUMPING;
    }
    throw std::invalid_argument("unsupported source/sink type: " + value);
}

PolygonBoundaryCondition loadPolygonBoundaryCondition(
    const Config& config,
    const std::filesystem::path& config_path,
    const std::string& path)
{
    const auto selector_type = config.getOr<std::string>(path + ".selector.type", "polygon");
    if (selector_type != "polygon") {
        throw std::invalid_argument("only polygon boundary selectors are supported in P16");
    }
    const auto file = config.get<std::string>(path + ".selector.file");
    return PolygonBoundaryCondition{
        Polygon::readFromFile(resolveInputPath(config_path, file).string()),
        boundaryConditionType(config.get<std::string>(path + ".type")),
        config.getOr<real>(path + ".value", 0.0),
        config.getOr<real>(path + ".normal_x", 1.0),
        config.getOr<real>(path + ".normal_y", 0.0),
    };
}

std::vector<PolygonBoundaryCondition> loadPolygonBoundaryConditions(
    const Config& config,
    const std::filesystem::path& config_path,
    const std::string& path)
{
    std::vector<PolygonBoundaryCondition> conditions;
    const int count = config.indexedCount(path);
    for (int i = 0; i < count; ++i) {
        conditions.push_back(
            loadPolygonBoundaryCondition(config, config_path, path + "." + std::to_string(i)));
    }
    return conditions;
}

PolygonSourceSink loadPolygonSource(
    const Config& config,
    const std::filesystem::path& config_path,
    const std::string& path)
{
    const auto selector_type = config.getOr<std::string>(path + ".selector.type", "polygon");
    if (selector_type != "polygon") {
        throw std::invalid_argument("only polygon source selectors are supported in P16");
    }
    if (config.getOr<std::string>(path + ".rate.source", "constant") != "constant") {
        throw std::invalid_argument("only constant polygon source rates are runtime-wired in P16");
    }
    const real rate = config.has(path + ".rate.value")
                          ? config.get<real>(path + ".rate.value")
                          : config.getOr<real>(path + ".rate", 0.0);
    const auto file = config.get<std::string>(path + ".selector.file");
    return PolygonSourceSink{
        Polygon::readFromFile(resolveInputPath(config_path, file).string()),
        sourceSinkType(config.get<std::string>(path + ".type")),
        rate,
    };
}

std::vector<PolygonSourceSink> loadPolygonSources(
    const Config& config,
    const std::filesystem::path& config_path,
    const std::string& path)
{
    std::vector<PolygonSourceSink> sources;
    const int count = config.indexedCount(path);
    for (int i = 0; i < count; ++i) {
        const auto source_path = path + "." + std::to_string(i);
        if (config.getOr<std::string>(source_path + ".rate.source", "constant") != "constant") {
            continue;
        }
        sources.push_back(loadPolygonSource(config, config_path, source_path));
    }
    return sources;
}

std::vector<SoilParameters> loadProductionSoilTable(const Config& config)
{
    std::vector<SoilParameters> soils;
    const int count = config.indexedCount("soil.types");
    soils.resize(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
        const std::string path = "soil.types." + std::to_string(i);
        const int id = config.get<int>(path + ".id");
        if (id < 0) {
            throw std::invalid_argument("soil type id must be non-negative");
        }
        if (static_cast<std::size_t>(id) >= soils.size()) {
            soils.resize(static_cast<std::size_t>(id + 1));
        }
        SoilParameters soil;
        soil.soil_a = config.get<real>(path + ".vg.alpha");
        soil.soil_n = config.get<real>(path + ".vg.n");
        soil.wcs = config.get<real>(path + ".vg.theta_s");
        soil.wcr = config.get<real>(path + ".vg.theta_r");
        soil.aev = config.getOr<real>(path + ".vg.aev", config.getOr<real>("groundwater.aev", -0.02));
        soil.ksx = config.get<real>(path + ".conductivity.Ksx");
        soil.ksy = config.get<real>(path + ".conductivity.Ksy");
        soil.ksz = config.get<real>(path + ".conductivity.Ksz");
        soils[static_cast<std::size_t>(id)] = soil;
    }
    return soils;
}

ReParameters parametersForRuntimeCell(
    const ReParameters& parameters,
    const LegacyReState& state,
    int cell)
{
    if (parameters.soil_table.empty()) {
        return parameters;
    }
    if (cell < 0 || cell >= static_cast<int>(state.soil_id.size())) {
        throw std::out_of_range("groundwater runtime IC cell index is out of range");
    }
    const int soil_id = state.soil_id[static_cast<std::size_t>(cell)];
    if (soil_id < 0 || soil_id >= static_cast<int>(parameters.soil_table.size())) {
        throw std::out_of_range("groundwater runtime IC soil id is outside the VG table");
    }
    return ReSolver::parametersForSoil(
        parameters,
        parameters.soil_table[static_cast<std::size_t>(soil_id)]);
}

void syncGroundwaterStateFromHead(const ReParameters& parameters, LegacyReState& state)
{
    for (int cell = 0; cell < static_cast<int>(state.h.size()); ++cell) {
        const auto cell_params = parametersForRuntimeCell(parameters, state, cell);
        state.wc[static_cast<std::size_t>(cell)] =
            ReSolver::waterContentFromHead(state.h[static_cast<std::size_t>(cell)], cell_params);
        state.hn[static_cast<std::size_t>(cell)] = state.h[static_cast<std::size_t>(cell)];
        state.hnm[static_cast<std::size_t>(cell)] = state.h[static_cast<std::size_t>(cell)];
        state.wcn[static_cast<std::size_t>(cell)] = state.wc[static_cast<std::size_t>(cell)];
        state.hwc[static_cast<std::size_t>(cell)] =
            ReSolver::headFromWaterContent(state.wc[static_cast<std::size_t>(cell)], cell_params);
        state.wch[static_cast<std::size_t>(cell)] = state.wc[static_cast<std::size_t>(cell)];
        state.ch[static_cast<std::size_t>(cell)] =
            ReSolver::specificCapacity(state.h[static_cast<std::size_t>(cell)], cell_params);
    }
}

void syncGroundwaterStateFromWaterContent(const ReParameters& parameters, LegacyReState& state)
{
    for (int cell = 0; cell < static_cast<int>(state.wc.size()); ++cell) {
        const auto cell_params = parametersForRuntimeCell(parameters, state, cell);
        state.h[static_cast<std::size_t>(cell)] =
            ReSolver::headFromWaterContent(state.wc[static_cast<std::size_t>(cell)], cell_params);
        state.hn[static_cast<std::size_t>(cell)] = state.h[static_cast<std::size_t>(cell)];
        state.hnm[static_cast<std::size_t>(cell)] = state.h[static_cast<std::size_t>(cell)];
        state.wcn[static_cast<std::size_t>(cell)] = state.wc[static_cast<std::size_t>(cell)];
        state.hwc[static_cast<std::size_t>(cell)] = state.h[static_cast<std::size_t>(cell)];
        state.wch[static_cast<std::size_t>(cell)] =
            ReSolver::waterContentFromHead(state.h[static_cast<std::size_t>(cell)], cell_params);
        state.ch[static_cast<std::size_t>(cell)] =
            ReSolver::specificCapacity(state.h[static_cast<std::size_t>(cell)], cell_params);
    }
}

void updateSurfaceDepthFromEta(LegacySweState& state)
{
    for (std::size_t i = 0; i < state.dept.size(); ++i) {
        state.dept[i] = std::max<real>(0.0, state.eta[i] - state.bottom[i]);
    }
}

void applyProductionSurfaceRuntimeInputs(
    const Config& config,
    const std::filesystem::path& config_path,
    const Grid& grid,
    LegacySweState& state,
    real dt,
    bool apply_initial_conditions)
{
    if (apply_initial_conditions && config.has("initial_conditions.surface.eta.source")) {
        state.eta = InitialCondition::loadSurfaceField(
            grid, initialConditionSpec(config, config_path, "initial_conditions.surface.eta", "eta"));
        updateSurfaceDepthFromEta(state);
    }

    const auto boundary_conditions =
        loadPolygonBoundaryConditions(config, config_path, "boundary_conditions.surface");
    if (!boundary_conditions.empty()) {
        BoundaryCondition::applySurface(grid, boundary_conditions, state.eta, state.bottom);
        updateSurfaceDepthFromEta(state);
    }

    const auto sources = loadPolygonSources(config, config_path, "sources.surface");
    if (!sources.empty()) {
        BoundaryCondition::applySourceSink(grid, sources, state.dept, dt);
        for (std::size_t i = 0; i < state.dept.size(); ++i) {
            state.eta[i] = state.bottom[i] + state.dept[i];
        }
    }
}

void applyProductionGroundwaterRuntimeInputs(
    const Config& config,
    const std::filesystem::path& config_path,
    const Grid& grid,
    const ReParameters& parameters,
    LegacyReState& state)
{
    const auto wc_source =
        config.getOr<std::string>("initial_conditions.groundwater.water_content.source", "");
    const bool wc_loaded = wc_source == "constant" || wc_source == "ascii_raster" ||
                           wc_source == "ascii_raster_3d" || wc_source == "hdf5_vector";
    if (wc_loaded) {
        state.wc = InitialCondition::loadGroundwaterField(
            grid,
            initialConditionSpec(
                config, config_path, "initial_conditions.groundwater.water_content", "water_content"));
    }

    const auto head_source =
        config.getOr<std::string>("initial_conditions.groundwater.hydraulic_head.source", "");
    const bool head_loaded = head_source == "constant" || head_source == "ascii_raster" ||
                             head_source == "ascii_raster_3d" || head_source == "hdf5_vector";
    if (head_loaded) {
        state.h = InitialCondition::loadGroundwaterField(
            grid,
            initialConditionSpec(
                config, config_path, "initial_conditions.groundwater.hydraulic_head", "hydraulic_head"));
    }

    const auto conditions =
        loadPolygonBoundaryConditions(config, config_path, "boundary_conditions.groundwater");
    if (!conditions.empty()) {
        BoundaryCondition::applyGroundwaterHead(grid, conditions, state.h);
    }
    const bool head_bc_applied = std::any_of(
        conditions.begin(),
        conditions.end(),
        [](const PolygonBoundaryCondition& condition) {
            return condition.type == BoundaryConditionType::GROUNDWATER_HEAD;
        });
    if (head_loaded || head_bc_applied) {
        syncGroundwaterStateFromHead(parameters, state);
    } else if (wc_loaded) {
        syncGroundwaterStateFromWaterContent(parameters, state);
    }
}

void applyProductionGroundwaterBoundaryConditions(
    const Config& config,
    const std::filesystem::path& config_path,
    const Grid& grid,
    LegacyReState& state)
{
    const auto conditions =
        loadPolygonBoundaryConditions(config, config_path, "boundary_conditions.groundwater");
    if (!conditions.empty()) {
        BoundaryCondition::applyGroundwaterHead(grid, conditions, state.h);
    }
}

bool hasProductionSurfaceRuntimeInputs(const Config& config)
{
    if (config.getOr<bool>("legacy_compatibility.required_by_current_driver", false) &&
        !config.getOr<bool>("runtime.enable", false)) {
        return false;
    }
    return config.has("initial_conditions.surface.eta.source") ||
           config.indexedCount("boundary_conditions.surface") > 0 ||
           config.indexedCount("sources.surface") > 0 ||
           config.indexedCount("monitoring.polygons") > 0;
}

bool hasProductionGroundwaterRuntimeInputs(const Config& config)
{
    if (config.getOr<bool>("legacy_compatibility.required_by_current_driver", false) &&
        !config.getOr<bool>("runtime.enable", false)) {
        return false;
    }
    return config.has("soil.map.source") ||
           config.indexedCount("soil.types") > 0 ||
           config.has("initial_conditions.groundwater.water_content.source") ||
           config.has("initial_conditions.groundwater.hydraulic_head.source") ||
           config.indexedCount("boundary_conditions.groundwater") > 0 ||
           config.indexedCount("monitoring.polygons") > 0;
}

void writeSurfaceRuntimeSummary(
    const Config& config,
    const std::filesystem::path& config_path,
    Hdf5Writer& writer,
    const std::filesystem::path& output_dir,
    const Grid& grid,
    const std::vector<real>& depth)
{
    const int count = config.indexedCount("monitoring.polygons");
    if (count == 0) {
        return;
    }

    std::ofstream summary(output_dir / "runtime_monitor_summary");
    summary << std::fixed << std::setprecision(12);
    for (int i = 0; i < count; ++i) {
        const std::string path = "monitoring.polygons." + std::to_string(i);
        const auto id = config.getOr<std::string>(path + ".id", "polygon_" + std::to_string(i));
        const auto file = config.get<std::string>(path + ".polygon_file");
        const auto polygon = Polygon::readFromFile(resolveInputPath(config_path, file).string());
        const auto cells = polygon.selectSurfaceCells(grid);
        const auto stats = Monitor::summarizeCells(depth, cells);
        summary << id << " count " << stats.count << " sum " << stats.sum << " mean "
                << stats.mean << " min " << stats.min << " max " << stats.max << "\n";
        writer.writeVector("/monitoring/" + id + "/stats", {
            static_cast<real>(stats.count), stats.sum, stats.mean, stats.min, stats.max});
    }
}

void writeGroundwaterRuntimeSummary(
    const Config& config,
    const std::filesystem::path& config_path,
    Hdf5Writer& writer,
    const std::filesystem::path& output_dir,
    const Grid& grid,
    const std::vector<real>& head,
    const LegacyReState& state)
{
    std::ofstream summary(output_dir / "runtime_monitor_summary");
    summary << std::fixed << std::setprecision(12);

    const int count = config.indexedCount("monitoring.polygons");
    for (int i = 0; i < count; ++i) {
        const std::string path = "monitoring.polygons." + std::to_string(i);
        const auto id = config.getOr<std::string>(path + ".id", "polygon_" + std::to_string(i));
        const auto file = config.get<std::string>(path + ".polygon_file");
        const auto polygon = Polygon::readFromFile(resolveInputPath(config_path, file).string());
        const auto cells = polygon.selectGroundwaterCells(grid);
        const auto stats = Monitor::summarizeCells(head, cells);
        summary << id << " count " << stats.count << " sum " << stats.sum << " mean "
                << stats.mean << " min " << stats.min << " max " << stats.max << "\n";
        writer.writeVector("/monitoring/" + id + "/stats", {
            static_cast<real>(stats.count), stats.sum, stats.mean, stats.min, stats.max});
    }

    if (!state.soil_id.empty()) {
        writer.writeIntVector("/groundwater/soil_id/0", state.soil_id);
    }
}

void writeSurfaceOutput(
    Hdf5Writer& writer,
    const std::filesystem::path& output_dir,
    int time_seconds,
    const std::vector<real>& eta,
    const std::vector<real>& depth_values,
    const std::vector<real>& uu,
    const std::vector<real>& vv,
    real min_depth,
    std::size_t cell_count,
    real surface_offset)
{
    std::vector<real> depth(cell_count, 0.0);
    std::vector<real> eta_out(cell_count, 0.0);
    std::vector<real> uu_out(cell_count, 0.0);
    std::vector<real> vv_out(cell_count, 0.0);
    std::vector<int> inundation(cell_count, 0);
    for (std::size_t i = 0; i < cell_count; ++i) {
        depth[i] = depth_values[i];
        eta_out[i] = eta[i] - surface_offset;
        uu_out[i] = uu[i];
        vv_out[i] = vv[i];
        inundation[i] = depth[i] > min_depth ? 1 : 0;
    }

    const auto suffix = std::to_string(time_seconds);
    writer.writeVector("/surface/water_depth/" + suffix, depth);
    writer.writeVector("/surface/water_surface_elevation/" + suffix, eta_out);
    writer.writeVector("/surface/velocity_x/" + suffix, uu_out);
    writer.writeVector("/surface/velocity_y/" + suffix, vv_out);
    writer.writeIntVector("/surface/inundation/" + suffix, inundation);

    auto writeText = [](const std::string& path, const auto& values) {
        std::ofstream output(path);
        output << std::scientific << std::setprecision(12);
        for (const auto& value : values) {
            output << value << " \n";
        }
    };
    writeText((output_dir / ("depth_" + suffix)).string(), depth);
    writeText((output_dir / ("surf_" + suffix)).string(), eta_out);
    writeText((output_dir / ("uu_" + suffix)).string(), uu_out);
    writeText((output_dir / ("vv_" + suffix)).string(), vv_out);
    writeText((output_dir / ("inun_" + suffix)).string(), inundation);
}

void writeGroundwaterOutput(
    Hdf5Writer& writer,
    const std::filesystem::path& output_dir,
    int time_seconds,
    const LegacyReState& state)
{
    const auto suffix = std::to_string(time_seconds);
    std::vector<real> qx_out = state.qx;
    std::vector<real> qy_out = state.qy;
    std::vector<real> qz_out = state.qz;
    for (auto& value : qx_out) {
        value *= 8.64e7;
    }
    for (auto& value : qy_out) {
        value *= 8.64e7;
    }
    for (auto& value : qz_out) {
        value *= 8.64e7;
    }

    writer.writeVector("/groundwater/hydraulic_head/" + suffix, state.h);
    writer.writeVector("/groundwater/water_content/" + suffix, state.wc);
    writer.writeVector("/groundwater/darcy_flux_x/" + suffix, qx_out);
    writer.writeVector("/groundwater/darcy_flux_y/" + suffix, qy_out);
    writer.writeVector("/groundwater/darcy_flux_z/" + suffix, qz_out);
    writer.writeVector("/groundwater/top_boundary_flow/" + suffix, state.qz_top_surface);

    auto writeText = [](const std::string& path, const auto& values) {
        std::ofstream output(path);
        for (const auto& value : values) {
            output << std::fixed << std::setprecision(6) << value << " \n";
        }
    };
    writeText((output_dir / ("head_" + suffix)).string(), state.h);
    writeText((output_dir / ("moisture_" + suffix)).string(), state.wc);
    writeText((output_dir / ("qx_" + suffix)).string(), qx_out);
    writeText((output_dir / ("qy_" + suffix)).string(), qy_out);
    writeText((output_dir / ("qz_" + suffix)).string(), qz_out);
    writeText((output_dir / ("qtop_" + suffix)).string(), state.qz_top_surface);
}

void writeGroundwaterGeometry(
    Hdf5Writer& writer,
    const std::filesystem::path& output_dir,
    const LegacyReState& state)
{
    writer.writeVector("/groundwater/zcell/0", state.zcntr);
    std::ofstream output(output_dir / "zcell_0");
    for (const auto& value : state.zcntr) {
        output << std::fixed << std::setprecision(6) << value << " \n";
    }
}

bool isCoupledRun(const Config& config)
{
    const bool modules_surface = config.getOr<bool>("modules.surface_water", false);
    const bool modules_groundwater = config.getOr<bool>("modules.groundwater", false);
    const bool coupling_enabled = config.getOr<bool>("coupling.enable", false);
    const auto mode = config.getOr<std::string>("simulation.mode", "");
    return mode == "coupled" || (modules_surface && modules_groundwater && coupling_enabled);
}

}  // namespace

SimulationDriver::SimulationDriver(std::filesystem::path config_path)
    : config_path_(std::move(config_path))
{
}

int SimulationDriver::run() const
{
    const Config config(config_path_.string());
    config.validate();

    if (isCoupledRun(config)) {
        return runCoupledSmoke();
    }
    if (config.get<bool>("surface_water.enable")) {
        return runBenchmarkSurfaceWater();
    }
    if (config.get<bool>("groundwater.enable")) {
        return runBenchmarkGroundwater();
    }
    return 0;
}

int SimulationDriver::runBenchmarkSurfaceWater() const
{
    const auto total_start = Clock::now();
    double setup_seconds = 0.0;
    double solver_seconds = 0.0;
    double runtime_input_seconds = 0.0;
    double output_seconds = 0.0;
    double summary_seconds = 0.0;
    const auto setup_start = Clock::now();

    const Config config(config_path_.string());
    config.validate();
    if (!config.get<bool>("surface_water.enable")) {
        return 0;
    }

    auto grid_spec = config.gridSpec();
    grid_spec.nz = 1;
    const Grid grid(grid_spec);

    SweParameters params;
    params.dt = config.get<real>("time.dt");
    params.gravity = config.get<real>("surface_water.grav");
    params.min_depth = config.get<real>("surface_water.min_depth");
    params.manning = config.get<real>("surface_water.manning");
    params.viscx = config.get<real>("surface_water.viscx");
    params.viscy = config.get<real>("surface_water.viscy");
    params.wtfh = config.get<real>("surface_water.wtfh");
    params.hD = config.get<real>("surface_water.hD");
    for (const auto& condition :
         loadPolygonBoundaryConditions(config, config_path_, "boundary_conditions.surface")) {
        if (condition.type != BoundaryConditionType::FREE_OUTFLOW) {
            continue;
        }
        SweFreeOutflowBoundary boundary;
        boundary.cells = condition.polygon.selectSurfaceCells(grid);
        const real normal_length =
            std::sqrt(condition.normal_x * condition.normal_x + condition.normal_y * condition.normal_y);
        if (normal_length > 0.0) {
            boundary.normal_x = condition.normal_x / normal_length;
            boundary.normal_y = condition.normal_y / normal_length;
        }
        params.free_outflow_boundaries.push_back(std::move(boundary));
    }

    const std::string legacy_input_dir = config.get<std::string>("simulation.legacy_input_dir");
    const auto bath = readColumnFile(resolveBenchmarkInput(config_path_, legacy_input_dir, "bath"));
    const auto rain_series =
        readColumnFile(resolveBenchmarkInput(config_path_, legacy_input_dir, "rain"));

    const auto n = static_cast<std::size_t>(grid.nSurfaceCell());
    if (bath.size() != n) {
        throw std::runtime_error("bath input size does not match surface grid");
    }
    const auto min_bath = *std::min_element(bath.begin(), bath.end());
    const real surface_offset = min_bath < 0.0 ? -min_bath : 0.0;
    std::vector<real> shifted_bath = bath;
    for (auto& value : shifted_bath) {
        value += surface_offset;
    }

    const auto output_path = resolveOutputPath(config, config_path_);
    const auto output_dir = output_path.parent_path();
    Hdf5Writer writer(output_path.string());
    writer.writeMetadata("b1-sw", "phase4");

    const auto tend = config.get<real>("time.Tend");
    const auto dt_out = config.get<real>("time.dt_out");
    SweSolver solver(grid, params);
    auto state = solver.initializeLegacyState(
        shifted_bath, config.get<real>("surface_water.init_eta") + surface_offset);
    const bool use_production_runtime = hasProductionSurfaceRuntimeInputs(config);
    if (use_production_runtime) {
        const auto runtime_start = Clock::now();
        applyProductionSurfaceRuntimeInputs(config, config_path_, grid, state, params.dt, true);
        runtime_input_seconds += elapsedSeconds(runtime_start);
    }
    setup_seconds = elapsedSeconds(setup_start) - runtime_input_seconds;

    const auto initial_output_start = Clock::now();
    writeSurfaceOutput(
        writer, output_dir, 0, state.eta, state.dept, state.uu, state.vv, params.min_depth, n,
        surface_offset);
    output_seconds += elapsedSeconds(initial_output_start);

    ProgressReporter progress(config.get<std::string>("simulation.id"), tend);
    real last_output = 0.0;
    for (real time = 0.0; time < tend - 1.0e-12;) {
        const auto t_current = time + params.dt;
        const auto solver_start = Clock::now();
        solver.advanceLegacyStep(
            state, legacySeriesValueAt(rain_series, t_current),
            config.get<real>("surface_water.q_evap"));
        solver_seconds += elapsedSeconds(solver_start);
        if (use_production_runtime) {
            const auto runtime_start = Clock::now();
            applyProductionSurfaceRuntimeInputs(config, config_path_, grid, state, params.dt, false);
            runtime_input_seconds += elapsedSeconds(runtime_start);
        }

        if (std::abs(t_current - last_output - dt_out) <= params.dt) {
            const auto t_save = std::round(t_current / dt_out) * dt_out;
            last_output = t_save;
            const auto output_start = Clock::now();
            writeSurfaceOutput(
                writer, output_dir, static_cast<int>(std::llround(t_save)), state.eta, state.dept,
                state.uu, state.vv, params.min_depth, n, surface_offset);
            output_seconds += elapsedSeconds(output_start);
        }
        time = t_current;
        progress.update(time);
    }
    if (use_production_runtime) {
        const auto summary_start = Clock::now();
        writeSurfaceRuntimeSummary(config, config_path_, writer, output_dir, grid, state.dept);
        summary_seconds += elapsedSeconds(summary_start);
    }

    RuntimeReport report;
    report.simulation_id = config.get<std::string>("simulation.id");
    report.module = "surface_water";
    report.total_seconds = elapsedSeconds(total_start);
    report.sections = {
        {"setup", setup_seconds},
        {"surface_solver", solver_seconds},
        {"runtime_inputs", runtime_input_seconds},
        {"output", output_seconds},
        {"monitoring_summary", summary_seconds},
    };
    writeRuntimeReport(output_dir, report);

    return 0;
}

int SimulationDriver::runBenchmarkGroundwater() const
{
    const auto total_start = Clock::now();
    double setup_seconds = 0.0;
    double solver_seconds = 0.0;
    double runtime_input_seconds = 0.0;
    double output_seconds = 0.0;
    double summary_seconds = 0.0;
    const auto setup_start = Clock::now();

    const Config config(config_path_.string());
    config.validate();
    if (!config.get<bool>("groundwater.enable")) {
        return 0;
    }

    const auto bath = 0.0;
    const auto bot_z = config.get<real>("domain.botZ");
    auto grid_spec = config.gridSpec();
    grid_spec.nz = ReSolver::legacyNz(
        bath, bot_z, grid_spec.dz, config.get<real>("domain.dz_incre"));
    grid_spec.dz_multiplier = config.get<real>("domain.dz_incre");
    const Grid grid(grid_spec);

    ReParameters params;
    params.dt = config.get<real>("time.dt");
    params.dt_min = config.get<real>("groundwater.dt_min");
    params.dt_max = config.get<real>("groundwater.dt_max");
    params.co_max = config.get<real>("groundwater.Co_max");
    params.ksx = config.get<real>("groundwater.Ksx");
    params.ksy = config.get<real>("groundwater.Ksy");
    params.ksz = config.get<real>("groundwater.Ksz");
    params.ss = config.get<real>("groundwater.Ss");
    params.soil_a = config.get<real>("groundwater.soil_a");
    params.soil_n = config.get<real>("groundwater.soil_n");
    params.wcs = config.get<real>("groundwater.wcs");
    params.wcr = config.get<real>("groundwater.wcr");
    params.init_wc = config.get<real>("groundwater.init_wc");
    params.init_h = config.get<real>("groundwater.init_h");
    params.init_wt_rel = config.get<real>("groundwater.init_wt_rel");
    params.init_wt_abs = config.get<real>("groundwater.init_wt_abs");
    params.qtop = config.get<real>("groundwater.qtop");
    params.qbot = config.get<real>("groundwater.qbot");
    params.htop = config.get<real>("groundwater.htop");
    params.hbot = config.get<real>("groundwater.hbot");
    params.qyp = config.get<real>("groundwater.qyp");
    params.qym = config.get<real>("groundwater.qym");
    params.aev = config.get<real>("groundwater.aev");
    params.use_vg = config.get<bool>("groundwater.use_vg");
    params.use_mvg = config.get<bool>("groundwater.use_mvg");
    params.use_full3d = config.get<bool>("groundwater.use_full3d");
    params.use_corrector = config.get<bool>("groundwater.use_corrector");
    params.dt_adjust = config.get<bool>("groundwater.dt_adjust");
    params.follow_terrain = config.get<bool>("groundwater.follow_terrain");
    const auto bc = config.get<std::vector<int>>("groundwater.bc_type");
    if (bc.size() != params.bc_type.size()) {
        throw std::runtime_error("groundwater.bc_type must contain six entries");
    }
    std::copy(bc.begin(), bc.end(), params.bc_type.begin());
    if (config.indexedCount("soil.types") > 0) {
        params.soil_table = loadProductionSoilTable(config);
    }
    if (config.getOr<std::string>("soil.map.source", "") == "ascii_raster") {
        const auto soil_map_file =
            resolveInputPath(config_path_, config.get<std::string>("soil.map.file"));
        params.soil_id = ReSolver::readSoilMap(grid, soil_map_file.string());
    }
    if (config.getOr<std::string>("soil.map.source", "") == "ascii_raster_3d") {
        const auto soil_map_file =
            resolveInputPath(config_path_, config.get<std::string>("soil.map.file"));
        params.soil_id = ReSolver::readSoilMap(grid, soil_map_file.string());
    }
    if (config.getOr<std::string>("soil.map.source", "") == "constant") {
        const int soil_id = config.getOr<int>("soil.map.value", 0);
        params.soil_id.assign(static_cast<std::size_t>(grid.nCell()), soil_id);
    }
    const auto groundwater_conditions =
        loadPolygonBoundaryConditions(config, config_path_, "boundary_conditions.groundwater");
    if (!groundwater_conditions.empty()) {
        params.qtop_surface.assign(static_cast<std::size_t>(grid.nSurfaceCell()), params.qtop);
        if (BoundaryCondition::applyGroundwaterTopFlux(
                grid,
                groundwater_conditions,
                params.qtop_surface)) {
            params.bc_type[5] = 2;
        }
    }

    const auto output_path = resolveOutputPath(config, config_path_);
    const auto output_dir = output_path.parent_path();
    Hdf5Writer writer(output_path.string());
    writer.writeMetadata("b2-gw", "phase5");

    ReSolver solver(grid, params);
    auto state = solver.initializeLegacyState(bath, bot_z);
    const bool use_production_runtime = hasProductionGroundwaterRuntimeInputs(config);
    if (use_production_runtime) {
        const auto runtime_start = Clock::now();
        applyProductionGroundwaterRuntimeInputs(config, config_path_, grid, params, state);
        runtime_input_seconds += elapsedSeconds(runtime_start);
    }
    setup_seconds = elapsedSeconds(setup_start) - runtime_input_seconds;

    const auto initial_output_start = Clock::now();
    writeGroundwaterGeometry(writer, output_dir, state);
    writeGroundwaterOutput(writer, output_dir, 0, state);
    output_seconds += elapsedSeconds(initial_output_start);

    const auto tend = config.get<real>("time.Tend");
    const auto dt_out = config.get<real>("time.dt_out");
    float last_output = 0.0F;
    std::ofstream timestep_output(output_dir / "timestep");
    ProgressReporter progress(config.get<std::string>("simulation.id"), tend);
    for (float time = 0.0F; time < static_cast<float>(tend);) {
        const auto current_dt = state.dtg;
        const auto t_current = static_cast<float>(time + static_cast<float>(current_dt));
        const auto solver_start = Clock::now();
        solver.advanceLegacyStep(state);
        solver_seconds += elapsedSeconds(solver_start);
        if (use_production_runtime) {
            const auto runtime_start = Clock::now();
            applyProductionGroundwaterBoundaryConditions(config, config_path_, grid, state);
            runtime_input_seconds += elapsedSeconds(runtime_start);
        }
        if (std::abs(t_current - last_output - static_cast<float>(dt_out)) <= state.dtg) {
            const auto t_save = std::round(t_current / dt_out) * dt_out;
            last_output = static_cast<float>(t_save);
            const auto output_start = Clock::now();
            writeGroundwaterOutput(writer, output_dir, static_cast<int>(std::llround(t_save)), state);
            output_seconds += elapsedSeconds(output_start);
        }
        const auto output_start = Clock::now();
        timestep_output << std::fixed << std::setprecision(8) << t_current << " \n";
        output_seconds += elapsedSeconds(output_start);
        time = t_current;
        progress.update(time);
    }
    if (use_production_runtime) {
        const auto summary_start = Clock::now();
        writeGroundwaterRuntimeSummary(config, config_path_, writer, output_dir, grid, state.h, state);
        summary_seconds += elapsedSeconds(summary_start);
    }

    RuntimeReport report;
    report.simulation_id = config.get<std::string>("simulation.id");
    report.module = "groundwater";
    report.total_seconds = elapsedSeconds(total_start);
    report.sections = {
        {"setup", setup_seconds},
        {"groundwater_solver", solver_seconds},
        {"runtime_inputs", runtime_input_seconds},
        {"output", output_seconds},
        {"monitoring_summary", summary_seconds},
    };
    writeRuntimeReport(output_dir, report);

    return 0;
}

int SimulationDriver::runCoupledSmoke() const
{
    const auto total_start = Clock::now();
    double setup_seconds = 0.0;
    double coupling_seconds = 0.0;
    double output_seconds = 0.0;
    const auto setup_start = Clock::now();

    const Config config(config_path_.string());
    config.validate();

    const auto output_path = resolveOutputPath(config, config_path_);
    const auto output_dir = output_path.parent_path();
    Hdf5Writer writer(output_path.string());
    writer.writeMetadata(config.get<std::string>("simulation.id"), "p15-coupled-smoke");

    CouplingParameters parameters;
    parameters.min_depth = config.getOr<real>("coupling.min_depth", 1.0e-8);
    const Coupling coupling(parameters);

    std::vector<AsyncCoupledColumn> columns{
        AsyncCoupledColumn{
            config.getOr<real>("coupled_smoke.surface_depth", 0.01),
            config.getOr<real>("coupled_smoke.surface_eta", 0.01),
            config.getOr<real>("coupled_smoke.groundwater_storage_depth", 0.0),
            config.getOr<real>("coupled_smoke.groundwater_head", -0.02),
            config.getOr<real>("coupled_smoke.dz_top", 0.1),
            config.getOr<real>("coupled_smoke.saturated_conductivity", 1.0e-5),
            config.getOr<real>("coupled_smoke.face_area", 1.0),
        },
    };
    AsyncCouplingClock clock;
    clock.async = config.getOr<std::string>("coupling.mode", "async") == "async";
    setup_seconds = elapsedSeconds(setup_start);

    ProgressReporter progress(config.get<std::string>("simulation.id"), 1.0);
    const auto coupling_start = Clock::now();
    const auto result = coupling.asyncCoupling(
        columns,
        clock,
        config.getOr<real>("coupling.surface_dt", config.get<real>("time.dt")),
        config.getOr<real>("coupling.groundwater_dt", config.get<real>("time.dt")),
        config.getOr<real>("coupled_smoke.rain_rate", 0.0),
        config.getOr<real>("coupled_smoke.evap_rate", 0.0));
    coupling_seconds += elapsedSeconds(coupling_start);
    progress.update(1.0);

    const auto output_start = Clock::now();
    writer.writeVector("/coupling/exchange_rate/0", result.exchange_rate);
    writer.writeVector("/coupling/surface_depth/0", {columns.front().surface_depth});
    writer.writeVector("/coupling/groundwater_storage_depth/0", {columns.front().groundwater_storage_depth});
    writer.writeVector("/coupling/time/0", {result.surface_time, result.groundwater_time});

    std::ofstream summary(output_dir / "coupled_smoke_summary");
    summary << std::fixed << std::setprecision(12);
    summary << "surface_time " << result.surface_time << "\n";
    summary << "groundwater_time " << result.groundwater_time << "\n";
    summary << "groundwater_steps " << result.groundwater_steps << "\n";
    summary << "surface_depth " << columns.front().surface_depth << "\n";
    summary << "groundwater_storage_depth " << columns.front().groundwater_storage_depth << "\n";
    summary << "surface_volume_change " << result.surface_volume_change << "\n";
    summary << "groundwater_volume_change " << result.groundwater_volume_change << "\n";
    output_seconds += elapsedSeconds(output_start);

    RuntimeReport report;
    report.simulation_id = config.get<std::string>("simulation.id");
    report.module = "coupled";
    report.total_seconds = elapsedSeconds(total_start);
    report.sections = {
        {"setup", setup_seconds},
        {"coupling_solver", coupling_seconds},
        {"output", output_seconds},
    };
    writeRuntimeReport(output_dir, report);

    return 0;
}

}  // namespace frehg2
