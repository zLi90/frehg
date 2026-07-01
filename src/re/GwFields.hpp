// Halo-padded 3D groundwater work arrays (P5).
#ifndef FREHG2_RE_GW_FIELDS_HPP
#define FREHG2_RE_GW_FIELDS_HPP

#include <string>
#include <utility>
#include <vector>

#include "core/Grid.hpp"
#include "frehg2/core/define.hpp"

namespace frehg2 {

struct GwFields {
  explicit GwFields(const Grid& grid) {
    const size_t n = static_cast<size_t>(grid.nCell());
    auto a = [n](const char* name) { return RealArr1DHost(name, n); };
    h = a("h");
    hn = a("hn");
    wc = a("wc");
    wcn = a("wcn");
    wch = a("wch");
    hwc = a("hwc");
    ch = a("ch");
    Kx = a("Kx");
    Ky = a("Ky");
    Kz = a("Kz");
    qx = a("qx");
    qy = a("qy");
    qz = a("qz");
    Gxp = a("Gxp");
    Gxm = a("Gxm");
    Gyp = a("Gyp");
    Gym = a("Gym");
    Gzp = a("Gzp");
    Gzm = a("Gzm");
    Gct = a("Gct");
    Grhs = a("Grhs");
    r_rho = a("r_rho");
    r_rhoxp = a("r_rhoxp");
    r_rhoyp = a("r_rhoyp");
    r_rhozp = a("r_rhozp");
    r_viscxp = a("r_viscxp");
    r_viscyp = a("r_viscyp");
    r_visczp = a("r_visczp");
    room = a("room");
    bot3d = a("bot3d");
    dz3d = a("dz3d");
    actv = a("actv");
  }

  RealArr1DHost h, hn, wc, wcn, wch, hwc, ch;
  RealArr1DHost Kx, Ky, Kz;
  RealArr1DHost qx, qy, qz;
  RealArr1DHost Gxp, Gxm, Gyp, Gym, Gzp, Gzm, Gct, Grhs;
  RealArr1DHost r_rho, r_rhoxp, r_rhoyp, r_rhozp;
  RealArr1DHost r_viscxp, r_viscyp, r_visczp;
  RealArr1DHost room;
  RealArr1DHost bot3d, dz3d, actv;

  // Every field array, paired with a stable name (full-state checkpoint/restart, P7.4).
  std::vector<std::pair<std::string, RealArr1DHost*>> namedViews() {
    return {{"h", &h},         {"hn", &hn},       {"wc", &wc},
            {"wcn", &wcn},     {"wch", &wch},     {"hwc", &hwc},
            {"ch", &ch},       {"Kx", &Kx},       {"Ky", &Ky},
            {"Kz", &Kz},       {"qx", &qx},       {"qy", &qy},
            {"qz", &qz},       {"Gxp", &Gxp},     {"Gxm", &Gxm},
            {"Gyp", &Gyp},     {"Gym", &Gym},     {"Gzp", &Gzp},
            {"Gzm", &Gzm},     {"Gct", &Gct},     {"Grhs", &Grhs},
            {"r_rho", &r_rho}, {"r_rhoxp", &r_rhoxp}, {"r_rhoyp", &r_rhoyp},
            {"r_rhozp", &r_rhozp}, {"r_viscxp", &r_viscxp}, {"r_viscyp", &r_viscyp},
            {"r_visczp", &r_visczp}, {"room", &room}, {"bot3d", &bot3d},
            {"dz3d", &dz3d},   {"actv", &actv}};
  }
};

}  // namespace frehg2

#endif  // FREHG2_RE_GW_FIELDS_HPP
