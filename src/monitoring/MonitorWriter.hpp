#ifndef FREHG2_MONITORING_MONITOR_WRITER_HPP
#define FREHG2_MONITORING_MONITOR_WRITER_HPP

#include <string>

#include "core/Grid.hpp"
#include "core/MpiComm.hpp"
#include "frehg2/core/define.hpp"
#include "monitoring/MonitorSpec.hpp"

namespace frehg2 {

class ReSolver;
class SweSolver;

// CSV monitor writer (P15.3.3). Rank 0 owns the file; each row is fsync'd for restart safety.
class MonitorWriter {
 public:
  MonitorWriter() = default;

  void configure(MonitorBundle bundle, std::string output_dir, std::string run_name, int gnx,
                 int gny, int gnz, double dx, double dy, double dz, double x0, double y0,
                 double botz, const MpiComm* mc);

  // Open (truncate) or resume (append when header matches). Call once after configure().
  void open(bool resume);

  // Sample probes/lines and append one CSV row (no-op when empty or rank != 0). The optional
  // surf_conc / subs_conc halo-padded host concentration fields (P16, owned by the Orchestrator,
  // not the solvers) back the "C"/"conc" probe field; pass nullptr when solute is disabled.
  void writeRow(real time, const SweSolver* swe, const ReSolver* re, const Grid& swe_grid,
                const Grid& gw_grid, const RealArr1DHost* surf_conc = nullptr,
                const RealArr1DHost* subs_conc = nullptr);

  void close();

  bool active() const { return !columns_.empty(); }
  const std::string& path() const { return path_; }

 private:
  struct Column {
    enum class Kind { Probe, Line };
    Kind kind = Kind::Probe;
    size_t probe_index = 0;
    size_t line_index = 0;
    std::string header;
  };

  double sampleProbeField(const ProbeSpec& probe, const std::string& field, const SweSolver* swe,
                          const ReSolver* re, const Grid& swe_grid, const Grid& gw_grid,
                          const RealArr1DHost* surf_conc, const RealArr1DHost* subs_conc) const;

  MonitorBundle bundle_;
  std::string output_dir_;
  std::string run_name_;
  std::string path_;
  std::vector<Column> columns_;
  int gnx_ = 0;
  int gny_ = 0;
  int gnz_ = 0;
  double dx_ = 1.0;
  double dy_ = 1.0;
  double dz_ = 1.0;
  double x0_ = 0.0;
  double y0_ = 0.0;
  double botz_ = 0.0;
  const MpiComm* mc_ = nullptr;
  bool open_ = false;
  bool header_written_ = false;
};

}  // namespace frehg2

#endif  // FREHG2_MONITORING_MONITOR_WRITER_HPP
