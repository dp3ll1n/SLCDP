#include <RcppEigen.h>
#include <random>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
Eigen::MatrixXd generateCloneTraj(int NCell, Eigen::MatrixXd M, Eigen::VectorXd dupRate, Eigen::VectorXd deathRate, Eigen::MatrixXd diffRate, double dt, double tEnd, Eigen::VectorXd E0) {
  int nrow=(tEnd/dt)+1.0;
  double cumDt=dt;
  int rind=0;
  Eigen::MatrixXd DupRate= Eigen::MatrixXd::Zero(1,NCell);
  Eigen::MatrixXd DeathRate= Eigen::MatrixXd::Zero(1,NCell);
  Eigen::MatrixXd X= Eigen::MatrixXd::Zero(1,NCell);
  Eigen::MatrixXd Xt= Eigen::MatrixXd::Zero(1,NCell);
  Eigen::MatrixXd Xout= Eigen::MatrixXd::Zero(nrow,NCell);

  Eigen::MatrixXd thetas= Eigen::MatrixXd::Zero(1,NCell*(NCell+2));
  DupRate.row(0)=dupRate;
  DeathRate.row(0)=deathRate;
  diffRate.resize(1,diffRate.cols()*diffRate.rows());

  thetas<<DupRate,DeathRate,diffRate;
  X.row(0) = E0;
  Xout.row(0) = E0;

  Eigen::MatrixXd rateVec(1,NCell*(NCell+2));

  std::random_device rd1;
  std::mt19937 gen1(rd1());
  std::random_device rd2;
  std::mt19937 gen2(rd2());

  std::vector<double> init(1*(NCell*(NCell+2)));
  double holdTime=0;
  double holdTimes=0;
  while ( holdTimes <= tEnd ) {

    rateVec = X.row(0).replicate(1,NCell+2).array() * thetas.array();
    for (int i = 0; i < NCell; i++){ rateVec(0,(i+NCell))=rateVec(0,(i+NCell))* X(0,i); }

    std::exponential_distribution<double> d1(rateVec.sum());
    holdTime = d1(gen1);
    Eigen::Map<Eigen::MatrixXd>(init.data(), 1,NCell*(NCell+2))=rateVec;
    std::discrete_distribution<> d2(init.begin(),init.end());
    int event=d2(gen2);

    Xt.row(0)=X.row(0)+M.row(event);
    int rowSum=Xt.row(0).sum();
    if(rowSum!=0){
      holdTimes+=holdTime;
      while((holdTimes > cumDt) && (cumDt<=tEnd)) {
        rind++;
        Xout.row(rind) = X.row(0);
        cumDt+=dt;
      }
      X=Xt;
    }
  }

  return(Xout);
}

