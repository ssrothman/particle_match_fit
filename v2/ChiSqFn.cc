#include "ChiSqFn.h"
#include "SRothman/SimonTools/src/util.h"
#include "SRothman/SimonTools/src/deltaR.h"

matching::ChiSqFn::ChiSqFn(const std::string& ptres_mode,
                           const double ptres_param1,
                           const double ptres_param2,
                           const std::string& angres_mode,
                           const double angres_param1,
                           const double angres_param2,
                           const double opp_charge_penalty,
                           const double no_charge_penalty) :
    ptresfunc(ResFunc::get_resfunc(ptres_mode, 
                               ptres_param1,
                               ptres_param2)),
    angresfunc(ResFunc::get_resfunc(angres_mode, 
                                angres_param1, 
                                angres_param2)),
    opp_charge_penalty(opp_charge_penalty),
    no_charge_penalty(no_charge_penalty) {}

double matching::ChiSqFn::evaluate(const double pt1, 
                                   const double eta1, 
                                   const double phi1, 
                                   const int charge1,
                                   const double pt2, 
                                   const double eta2, 
                                   const double phi2,
                                   const int charge2) const {

    double ptres = ptresfunc->evaluate(pt1, eta1, 
                                       phi1, charge1);
    double angres = angresfunc->evaluate(pt1, eta1, 
                                         phi1, charge1);

    double pt_term = simon::square((pt1 - pt2) / ptres);
    double ang_term = simon::deltaR2(eta1, phi1, eta2, phi2)/simon::square(angres);
    int charge_product = charge1 * charge2;
    double charge_term = 0;
    if (charge_product < 0){
        charge_term = opp_charge_penalty;
    } else if (charge_product == 0){
        charge_term = no_charge_penalty;
    }

    return pt_term + ang_term + charge_term;
}
