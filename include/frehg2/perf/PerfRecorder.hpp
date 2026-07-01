// Aggregates named timings and counters for the P21 performance report (Tasks 21.3.1–21.3.2).
#ifndef FREHG2_PERF_PERF_RECORDER_HPP
#define FREHG2_PERF_PERF_RECORDER_HPP

#include <iosfwd>

#include <mpi.h>

#include "frehg2/perf/Counters.hpp"
#include "frehg2/perf/PerfRegions.hpp"

namespace frehg2 {

class PerfRecorder {
 public:
  void addTime(perf::Region region, double seconds);
  perf::Counters& counters() { return counters_; }
  const perf::Counters& counters() const { return counters_; }

  double timeSeconds(perf::Region region) const;
  double totalTimedSeconds() const;

  // Sum timings and counters across MPI ranks (total CPU-seconds / total work units).
  void mpiReduce(MPI_Comm comm);

  // Append perf_* lines to simulation_summary.txt (rank-0 only; caller guards rank).
  void writeSummaryLines(std::ostream& out, long long time_steps, long long gw_substeps) const;

 private:
  double times_[static_cast<unsigned>(perf::Region::Count)] = {};
  perf::Counters counters_;
};

}  // namespace frehg2

#endif  // FREHG2_PERF_PERF_RECORDER_HPP
