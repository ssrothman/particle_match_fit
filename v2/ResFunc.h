#ifndef SROTHMAN_MATCHING_V2_RESFUNC_H
#define SROTHMAN_MATCHING_V2_RESFUNC_H

#include <string>
#include <memory>

namespace matching {
    class ResFunc;
    using ResFuncPtr = std::unique_ptr<const ResFunc>;

    class ResFunc {
    public:
        static ResFuncPtr get_resfunc(const std::string& mode,
                            const double param1,
                            const double param2);

        virtual double evaluate(const double pt,
                                const double eta,
                                const double phi,
                                const int charge) const = 0;

        virtual ~ResFunc() = default;
    };
};

#endif
