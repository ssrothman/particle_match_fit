#include "chisqLossFCN.h"
#include "SRothman/SimonTools/src/isID.h"
#include "matchingUtil.h"

ChisqLossFCN::ChisqLossFCN(const simon::jet& recojet,
                          const simon::jet& genjet,
                          const std::vector<std::pair<unsigned, unsigned>>& locations,
                          const enum spatialLoss type,
                          const std::vector<double>& PUpt0s,
                          const std::vector<double>& PUexps, 
                          const std::vector<double>& PUpenalties):
      recoPT(recojet.ptvec()), 
      recoETA(recojet.etavec()), 
      recoPHI(recojet.phivec()),
      genPT(genjet.ptvec()), 
      genETA(genjet.etavec()), 
      genPHI(genjet.phivec()),
      errPT(recojet.dptvec()), 
      errETA(recojet.detavec()), 
      errPHI(recojet.dphivec()),
      NPReco(recoPT.size()), NPGen(genPT.size()),
      weightedGenETA(genPT.array() * genETA.array()), 
      weightedGenPHI(genPT.array() * genPHI.array()),
      locations(locations),
      type(type),
      PUpt0s(PUpt0s),
      PUexps(PUexps), PUpenalties(PUpenalties),
      ids(){

    for(const auto& part : recojet.particles){
        if(isEM0(part)){
            ids.emplace_back(0);
        } else if(isHAD0(part)){
            ids.emplace_back(1);
        } else if(isHADCH(part)){
            ids.emplace_back(2);
        } else if(isELE(part)){
            ids.emplace_back(3);
        } else if(isMU(part)){
            ids.emplace_back(4);
        } else {
            throw std::runtime_error("Unrecognized particle type in recojet");
        }
    }
}

double ChisqLossFCN::operator()(const std::vector<double>& data) const {
    //printf("top of loss\n");
    //fflush(stdout);
    Eigen::MatrixXd A = fullmat(recoPT.size(), genPT.size(), locations, data);
    //printf("1\n");
    //fflush(stdout);

    double lossPT = 0;
    double lossETA = 0;
    double lossPHI = 0;
    double lossPU = 0;

    Eigen::VectorXd recoPT_pred = A * genPT;
    //printf("2\n");
    //fflush(stdout);

    Eigen::VectorXd PU = (recoPT_pred.array() == 0).cast<double>();
    //printf("3\n");
    //fflush(stdout);

    Eigen::VectorXd lossPTvec = (recoPT_pred - recoPT).array()/ errPT.array();
    lossPTvec.array() *= (1-PU.array());
    lossPT = lossPTvec.dot(lossPTvec);
    //printf("4\n");
    //fflush(stdout);

    for (unsigned i=0; i<recoPT.size(); ++i){
        if(recoPT_pred[i] == 0){
            recoPT_pred[i] = 1;
        }
    }
    //printf("5\n");
    //fflush(stdout);
    if (type == spatialLoss::TYPE1){
        Eigen::VectorXd recoETA_pred = A * weightedGenETA;
        recoETA_pred.array() /= recoPT_pred.array();
        Eigen::VectorXd lossETAvec = (recoETA_pred - recoETA).array()/ errETA.array();
        lossETAvec.array() *= (1-PU.array());
        lossETA = lossETAvec.dot(lossETAvec);
        //printf("6\n");
        //fflush(stdout);

        Eigen::VectorXd recoPHI_pred = A * weightedGenPHI;
        recoPHI_pred.array() /= recoPT_pred.array();
        Eigen::VectorXd lossPHIvec = (recoPHI_pred - recoPHI).array()/ errPHI.array();
        lossPHIvec.array() *= (1-PU.array());
        lossPHI = lossPHIvec.dot(lossPHIvec);
        //printf("7\n");
        //fflush(stdout);
    } else if (type == spatialLoss::TYPE2){
        /*Eigen::MatrixXd diffETA(NPReco, NPGen);
        diffETA.array().colwise() = recoETA.array();
        diffETA.array().rowwise() -= genETA.transpose().array();
        lossETA = diffETA.dot(diffETA);
        //printf("8\n");
        //fflush(stdout);

        Eigen::MatrixXd diffPHI(NPReco, NPGen);
        diffPHI.array().colwise() = recoPHI.array();
        diffPHI.array().rowwise() -= genPHI.transpose().array();
        lossPHI = diffPHI.dot(diffPHI);*/
        throw std::runtime_error("TYPE2 not implemented");
        //printf("9\n");
        //fflush(stdout);
    }

    for(unsigned i=0; i<NPReco; ++i){
        //printf("PU loss for particle %u\n", i);
        //fflush(stdout);
        if(PU[i] == 0) continue;
        //printf("passed PU check\n");
        //fflush(stdout);
        //printf("ids[i] = %u\n", ids[i]);
        //fflush(stdout);
        double pt0 = PUpt0s[ids[i]];
        double exp = PUexps[ids[i]];
        double penalty = PUpenalties[ids[i]];
        //printf("got parameters\n");
        //fflush(stdout);

        lossPU += 2*exp*std::max(std::log(recoPT[i]/pt0), 0.0) + penalty;
    }
    //printf("10\n");
    //fflush(stdout);

    return lossPT + lossETA + lossPHI + lossPU;
};

