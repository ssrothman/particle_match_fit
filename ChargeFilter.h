#ifndef SROTHMAN_MATCHING_V2_CHARGEFILTER_H
#define SROTHMAN_MATCHING_V2_CHARGEFILTER_H

#include <string>
#include <memory>

namespace matching{
    class ChargeFilter;
    using ChargeFilterPtr = std::unique_ptr<const ChargeFilter>;

    class ChargeFilter {
    public:
        virtual bool evaluate(
                const int reco_charge,
                const int gen_charge) const = 0;

        static ChargeFilterPtr get_charge_filter(
                const std::string& mode);

        virtual ~ChargeFilter() = default;
    };
};

#endif
