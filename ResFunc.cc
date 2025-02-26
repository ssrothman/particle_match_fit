#include "ResFunc.h"
#include <stdexcept>

class ConstRes : public matching::ResFunc {
public:
    ConstRes(const double res,
             [[maybe_unused]] const double B) : res(res) {}

    double evaluate(
            [[maybe_unused]] const double pt, 
            [[maybe_unused]] const double eta, 
            [[maybe_unused]] const double phi, 
            [[maybe_unused]] const int charge) const {
        return res;
    }
private:
    const double res;
};

class ConstFracRes : public matching::ResFunc {
public:
    ConstFracRes(const double res,
                 [[maybe_unused]] const double B) : res(res) {}

    double evaluate(
            const double pt, 
            [[maybe_unused]] const double eta, 
            [[maybe_unused]] const double phi, 
            [[maybe_unused]] const int charge) const {
        return res * pt;
    }
private:
    const double res;
};

class TrackPtRes: public matching::ResFunc {
public:
    TrackPtRes(const double A, const double B) : A(A), B(B) {}

    double evaluate(const double pt, 
                    [[maybe_unused]] const double eta, 
                    [[maybe_unused]] const double phi,
                    [[maybe_unused]] const int charge) const {
        return A + B * pt;
    }
private:
    const double A, B;
};

class TrackAngRes: public matching::ResFunc {
public:
    TrackAngRes(const double A, const double B) : A(A), B(B) {}

    double evaluate(
            const double pt, 
            [[maybe_unused]] const double eta, 
            [[maybe_unused]] const double phi, 
            [[maybe_unused]] const int charge) const {
        return A + B / pt;
    }
private:
    const double A, B;
};

matching::ResFuncPtr matching::ResFunc::get_resfunc(
        const std::string& mode,
        const double param1,
        const double param2) {

    if (mode == "Const") {
        return std::make_unique<const ConstRes>(param1, param2);
    } else if (mode == "ConstFrac") {
        return std::make_unique<const ConstFracRes>(param1, param2);
    } else if (mode == "TrackPt") {
        return std::make_unique<const TrackPtRes>(param1, param2);
    } else if (mode == "TrackAng") {
        return std::make_unique<const TrackAngRes>(param1, param2);
    } else {
        throw std::invalid_argument("Invalid resolution function mode");
    }
}
