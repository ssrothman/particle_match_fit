#include "ChargeFilter.h"

class ChargeMagnitudeFilter : public matching::ChargeFilter {
public:
    ChargeMagnitudeFilter() {}

    bool evaluate(
            const int reco_charge,
            const int gen_charge) const override {
        return std::abs(reco_charge) == std::abs(gen_charge);
    }
};

class ChargeSignFilter : public matching::ChargeFilter{
public:
    ChargeSignFilter() {}

    bool evaluate(
            const int reco_charge,
            const int gen_charge) const override {
        if (reco_charge == 0 && gen_charge == 0) {
            return true;
        } else {
            return reco_charge * gen_charge > 0;
        }
    }
};

class AnyChargeFilter : public matching::ChargeFilter{
public:
    AnyChargeFilter() {}

    bool evaluate(
            [[maybe_unused]] const int reco_charge,
            [[maybe_unused]] const int gen_charge) const override {
        return true;
    }
};

matching::ChargeFilterPtr matching::ChargeFilter::get_charge_filter(
        const std::string& mode) {
    if (mode == "Magnitude") {
        return std::make_unique<const ChargeMagnitudeFilter>();
    } else if (mode == "Sign") {
        return std::make_unique<const ChargeSignFilter>();
    } else if (mode == "Any") {
        return std::make_unique<const AnyChargeFilter>();
    } else {
        throw std::invalid_argument("Invalid charge filter mode");
    }
}
