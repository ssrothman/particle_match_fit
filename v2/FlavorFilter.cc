#include "FlavorFilter.h"

class AnyFlavorFilter : public matching::FlavorFilter {
public:
    AnyFlavorFilter() {}

    bool evaluate(
            [[maybe_unused]] const int gen_charge,
            [[maybe_unused]] const int gen_pdgid) const override {
        return true;
    }
};

class AnyHadronFilter : public matching::FlavorFilter {
public:
    AnyHadronFilter() {}

    bool evaluate(
            [[maybe_unused]] const int gen_charge,
            const int gen_pdgid) const override {
        return std::abs(gen_pdgid) >= 100;
    }
};

class AnyLeptonFilter : public matching::FlavorFilter {
public:
    AnyLeptonFilter() {}

    bool evaluate(
            [[maybe_unused]] const int gen_charge,
            const int gen_pdgid) const override {
        return gen_pdgid == 11 || gen_pdgid == 13 || gen_pdgid == 15;
    }
};

class ElectronFilter : public matching::FlavorFilter {
public:
    ElectronFilter() {}

    bool evaluate(
            [[maybe_unused]] const int gen_charge,
            const int gen_pdgid) const override {
        return std::abs(gen_pdgid) == 11;
    }
};

class MuonFilter : public matching::FlavorFilter {
public:
    MuonFilter() {}

    bool evaluate(
            [[maybe_unused]] const int gen_charge,
            const int gen_pdgid) const override {
        return std::abs(gen_pdgid) == 13;
    }
};

class ElectronMuonFilter : public matching::FlavorFilter {
public:
    ElectronMuonFilter() {}

    bool evaluate(
            [[maybe_unused]] const int gen_charge,
            const int gen_pdgid) const override {
        return std::abs(gen_pdgid) == 11 || std::abs(gen_pdgid) == 13;
    }
};

class ElectromagneticFilter : public matching::FlavorFilter {
public:
    ElectromagneticFilter() {}

    bool evaluate(
            [[maybe_unused]] const int gen_charge,
            const int gen_pdgid) const override {
        return std::abs(gen_pdgid) == 11 || std::abs(gen_pdgid) == 22 || std::abs(gen_pdgid==111);
    }
};

class AnyChargedFilter : public matching::FlavorFilter {
public:
    AnyChargedFilter() {}

    bool evaluate(
            const int gen_charge,
            [[maybe_unused]] const int gen_pdgid) const override {
        return gen_charge != 0;
    }
};

class AnyNeutralFilter : public matching::FlavorFilter {
public:
    AnyNeutralFilter() {}

    bool evaluate(
            const int gen_charge,
            [[maybe_unused]] const int gen_pdgid) const override {
        return gen_charge == 0;
    }
};

class AnyChargedHadronFilter : public matching::FlavorFilter {
public:
    AnyChargedHadronFilter() {}

    bool evaluate(
            const int gen_charge,
            const int gen_pdgid) const override {
        return gen_charge != 0 && std::abs(gen_pdgid) >= 100;
    }
};

class AnyNeutralHadronFilter : public matching::FlavorFilter {
public:
    AnyNeutralHadronFilter() {}

    bool evaluate(
            const int gen_charge,
            const int gen_pdgid) const override {
        return gen_charge == 0 && std::abs(gen_pdgid) >= 100;
    }
};

matching::FlavorFilterPtr matching::FlavorFilter::get_flavor_filter(
        const std::string& mode) {
    if (mode == "Any") {
        return std::make_unique<const AnyFlavorFilter>();
    } else if (mode == "AnyHadron") {
        return std::make_unique<const AnyHadronFilter>();
    } else if (mode == "AnyLepton") {
        return std::make_unique<const AnyLeptonFilter>();
    } else if (mode == "Electron") {
        return std::make_unique<const ElectronFilter>();
    } else if (mode == "Muon") {
        return std::make_unique<const MuonFilter>();
    } else if (mode == "ElectronMuon") {
        return std::make_unique<const ElectronMuonFilter>();
    } else if (mode == "Electromagnetic") {
        return std::make_unique<const ElectromagneticFilter>();
    } else if (mode == "AnyCharged") {
        return std::make_unique<const AnyChargedFilter>();
    } else if (mode == "AnyNeutral") {
        return std::make_unique<const AnyNeutralFilter>();
    } else if (mode == "AnyChargedHadron") {
        return std::make_unique<const AnyChargedHadronFilter>();
    } else if (mode == "AnyNeutralHadron") {
        return std::make_unique<const AnyNeutralHadronFilter>();
    } else {
        throw std::invalid_argument("Invalid flavor filter mode");
    }
}
