#include "DeltaRLimiter.h"

class ConstDeltaRLimiter: public matching::DeltaRLimiter {
public:
    ConstDeltaRLimiter(const double thresh,
                       [[maybe_unused]] const double B,
                       [[maybe_unused]] const double C):
        thresh(thresh) {}

    double evaluate(
            [[maybe_unused]] const double pt,
            [[maybe_unused]] const double eta,
            [[maybe_unused]] const double phi) const {
        return thresh;
    }
private:
    const double thresh;
};

class TrackPtDeltaRLimiter : public matching::DeltaRLimiter {
public:
    TrackPtDeltaRLimiter(const double A, const double B, const double C):
        A(A), B(B), C(C) {}

    double evaluate(
            const double pt,
            [[maybe_unused]] const double eta,
            [[maybe_unused]] const double phi) const {

        return std::min(A + B / pt, C);
    }
private:
    const double A, B, C;
};

matching::DeltaRLimiterPtr matching::DeltaRLimiter::get_deltaRlimiter(
        const std::string& mode,
        const double param1,
        const double param2,
        const double param3) {

    if (mode == "Const") {
        return std::make_unique<const ConstDeltaRLimiter>(param1, param2, param3);
    } else if (mode == "TrackPt") {
        return std::make_unique<const TrackPtDeltaRLimiter>(param1, param2, param3);
    } else {
        throw std::invalid_argument("Invalid delta R limiter mode");
    }
}
