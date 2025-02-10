#include "DeltaRLimiter.h"

class ConstDeltaRLimiter: public matching::DeltaRLimiter {
public:
    ConstDeltaRLimiter(const double thresh,
                       [[maybe_unused]] const double B) :
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
    TrackPtDeltaRLimiter(const double A, const double B):
        A(A), B(B) {}

    double evaluate(
            const double pt,
            [[maybe_unused]] const double eta,
            [[maybe_unused]] const double phi) const {

        return A + B / pt;
    }
private:
    const double A, B;
};

matching::DeltaRLimiterPtr matching::DeltaRLimiter::get_deltaRlimiter(
        const std::string& mode,
        const double param1,
        const double param2) {

    if (mode == "Const") {
        return std::make_unique<const ConstDeltaRLimiter>(param1, param2);
    } else if (mode == "TrackPt") {
        return std::make_unique<const TrackPtDeltaRLimiter>(param1, param2);
    } else {
        throw std::invalid_argument("Invalid delta R limiter mode");
    }
}
