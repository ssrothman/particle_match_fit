#ifndef SROTHMAN_MATCHING_V2_FLAVORFILTER_H
#define SROTHMAN_MATCHING_V2_FLAVORFILTER_H

#include <string>
#include <memory>

namespace matching{
    class FlavorFilter;
    using FlavorFilterPtr = std::unique_ptr<const FlavorFilter>;

    class FlavorFilter {
    public:
        virtual bool evaluate(
                const int gen_charge,
                const int gen_pdgid) const = 0;

        static FlavorFilterPtr get_flavor_filter(
                const std::string& mode);

        virtual ~FlavorFilter() = default;
    };
};

#endif
