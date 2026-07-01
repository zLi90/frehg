#include "frehg2/perf/PerfRecorder.hpp"

#include <algorithm>
#include <iomanip>
#include <ostream>

namespace frehg2 {

void PerfRecorder::addTime(perf::Region region, double seconds) {
  if (seconds <= 0.0) return;
  times_[static_cast<unsigned>(region)] += seconds;
}

double PerfRecorder::timeSeconds(perf::Region region) const {
  return times_[static_cast<unsigned>(region)];
}

double PerfRecorder::totalTimedSeconds() const {
  double sum = 0.0;
  for (unsigned i = 0; i < static_cast<unsigned>(perf::Region::Count); ++i) sum += times_[i];
  return sum;
}

void PerfRecorder::mpiReduce(MPI_Comm comm) {
  int size = 1;
  MPI_Comm_size(comm, &size);
  if (size <= 1) return;
  double g_times[static_cast<unsigned>(perf::Region::Count)] = {};
  MPI_Allreduce(times_, g_times, static_cast<int>(perf::Region::Count), MPI_DOUBLE, MPI_SUM, comm);
  for (unsigned i = 0; i < static_cast<unsigned>(perf::Region::Count); ++i) times_[i] = g_times[i];
  int64_t loc[3] = {counters_.cells_touched, counters_.ksp_iterations, counters_.bytes_staged};
  int64_t glb[3] = {0, 0, 0};
  MPI_Allreduce(loc, glb, 3, MPI_INT64_T, MPI_SUM, comm);
  counters_.cells_touched = glb[0];
  counters_.ksp_iterations = glb[1];
  counters_.bytes_staged = glb[2];
}

void PerfRecorder::writeSummaryLines(std::ostream& out, long long time_steps,
                                     long long gw_substeps) const {
  out << std::setprecision(17);
  for (unsigned i = 0; i < static_cast<unsigned>(perf::Region::Count); ++i) {
    const auto r = static_cast<perf::Region>(i);
    out << "perf_" << perf::regionName(r) << "_seconds " << times_[i] << "\n";
  }
  const long long flow_steps = std::max(time_steps, 1LL);
  out << "perf_sw_assembly_seconds_per_step " << (times_[static_cast<unsigned>(perf::Region::SwAssembly)] / flow_steps) << "\n";
  out << "perf_sw_ksp_seconds_per_step " << (times_[static_cast<unsigned>(perf::Region::SwKsp)] / flow_steps) << "\n";
  out << "perf_sw_update_seconds_per_step " << (times_[static_cast<unsigned>(perf::Region::SwUpdate)] / flow_steps) << "\n";
  if (gw_substeps > 0) {
    out << "perf_gw_assembly_seconds_per_step " << (times_[static_cast<unsigned>(perf::Region::GwAssembly)] / gw_substeps) << "\n";
    out << "perf_gw_ksp_seconds_per_step " << (times_[static_cast<unsigned>(perf::Region::GwKsp)] / gw_substeps) << "\n";
  }
  out << "perf_cells_touched " << counters_.cells_touched << "\n";
  out << "perf_ksp_iterations " << counters_.ksp_iterations << "\n";
  out << "perf_bytes_staged " << counters_.bytes_staged << "\n";
  out << "perf_timed_fraction " << (totalTimedSeconds() > 0.0 ? totalTimedSeconds() : 0.0) << "\n";
}

}  // namespace frehg2
