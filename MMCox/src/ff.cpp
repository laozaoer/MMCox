#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/

// [[Rcpp::export]]
arma::field<arma::mat> UpdateonceProfileBeta(const int&Indicator,const arma::mat&rules,const arma::mat&beta,const double&theta,const arma::mat&gamma,const arma::mat&Data,
                                             const arma::umat&tL,const arma::umat&tR,const arma::umat&tLR){
  arma::mat WeightMat=WeightFunc(rules,beta,theta,gamma,Data,tL,tR);
  arma::mat index=sort(unique(Data.col(0)));
  int n=index.n_rows;
  int order=rules.n_rows;
  int m=gamma.n_rows;
  arma::mat covariate=Data.submat(0,5,Data.n_rows-1,Data.n_cols-1);
  arma::mat Deletecovariate=covariate.shed_col(Indicator - 1);
  arma::mat betaX=covariate*beta;
  arma::mat ParPart=exp(repmat(betaX,1,order)+repmat(trans(theta*rules.col(0)),Data.n_rows,1));
  
  arma::mat LambdaR=repmat(trans(trans(gamma)*tR),1,order);
  arma::mat LambdaL=repmat(trans(trans(gamma)*tL),1,order);
  arma::mat LambdaLR=repmat(trans(trans(gamma)*tLR),1,order);
  arma::mat expL=LambdaL%ParPart;
  arma::mat expR=LambdaR%ParPart;
  arma::mat expLR=LambdaLR%ParPart;
  
  arma::mat CubegammaR=repmat(trans(pow(gamma,3)),Data.n_rows,1)%trans(tR);
  arma::mat CubegammaLR=repmat(trans(pow(gamma,3)),Data.n_rows,1)%trans(tLR);
  arma::mat InversgammaL=repmat(trans(pow(gamma,-1)),Data.n_rows,1)%trans(tL);
  
  
  arma::mat SharedTerm1=repmat(Data.col(3),1,order)%pow((1+expR),-1);
  arma::mat SharedTerm2=(1-repmat(Data.col(3),1,order)-repmat(Data.col(4),1,order))%pow((1+expLR),-1);
  arma::mat SharedTerm3=(1-repmat(Data.col(3),1,order))%expL;
  
  arma::mat Term1=trans(SharedTerm1%pow((LambdaR),-1));
  arma::mat Term2=trans(SharedTerm2%pow((LambdaLR),-1));
  arma::mat Term3=trans((1-repmat(Data.col(3),1,order))%ParPart);
  // std::cout<<pow((LambdaLR),-1);
  arma::mat Amat(n,m),Bmat(n,m);
  
  arma::mat FirstBigMat=SharedTerm1+SharedTerm2-SharedTerm3;
  
  arma::mat SecondBigMat=FirstBigMat+2*SharedTerm3;
  
  arma::mat FirstDerivMat(n,beta.n_rows);
  arma::mat SecondDeriv=arma::zeros(beta.n_rows,beta.n_rows);
  arma::mat SecondDerivMat=arma::zeros(beta.n_rows,beta.n_rows);
  for(int i=0;i<n;i++){
    arma::umat rowindices=find(Data.col(0)==(i+1));
    arma::mat AclusterI=Term1.cols(rowindices)*CubegammaR.rows(rowindices)+
      Term2.cols(rowindices)*CubegammaLR.rows(rowindices);
    arma::mat BclusterI=Term3.cols(rowindices)*InversgammaL.rows(rowindices);
    
    Amat.row(i)=WeightMat.row(i)*AclusterI;
    Bmat.row(i)=WeightMat.row(i)*BclusterI;
    // std::cout<<Amat;
    arma::mat FirstBigMatij=trans(FirstBigMat.rows(rowindices));
    arma::mat SecondBigMatij=SecondBigMat.rows(rowindices);
    arma::mat FirstDerivI(order,beta.n_rows);
    
    for(int r=0;r<order;r++){
      
      arma::mat Xij=Deletecovariate.rows(rowindices);
      arma::vec onesvec=arma::ones(Xij.n_rows);
      arma::mat Wijr=join_rows(Xij,rules(r,0)*onesvec*theta);
      FirstDerivI.row(r)=FirstBigMatij.row(r)*Wijr;
      arma::mat SecondDerivI=arma::zeros(beta.n_rows,beta.n_rows);
      
      for(int j=0;j<Xij.n_rows;j++){
        SecondDerivI=SecondDerivI+SecondBigMatij(j,r)*(trans(Wijr.row(j))*Wijr.row(j));
      }
      SecondDeriv=SecondDeriv-2*WeightMat(i,r)*SecondDerivI;
      
      
    }
    FirstDerivMat.row(i)=WeightMat.row(i)*FirstDerivI;
    
    
  }
  arma::mat FirstDeriv=trans(sum(FirstDerivMat,0));
  SecondDeriv(beta.n_rows-1,beta.n_rows-1)=SecondDeriv(beta.n_rows-1,beta.n_rows-1)+FirstDeriv(beta.n_rows-1);
  arma::mat updategamma=pow((trans(sum(Amat,0))%pow(trans(sum(Bmat,0)),-1)),0.25);
  
  arma::vec lastpar(beta.n_rows);
  lastpar.subvec(0,beta.n_rows-2)=beta.shed_col(Indicator-1);
  lastpar(beta.n_rows-1)=std::log(theta);
  
  arma::vec updatePar=lastpar-solve(SecondDeriv,FirstDeriv);
  updatePar(beta.n_rows-1)=std::exp(updatePar(beta.n_rows-1));
  arma::mat updateParAgg=undatePar.insert_rows(Indicator-1,beta.row(Indicator-1));
  
  arma::field<arma::mat> result(2);
  result(0)=updategamma;
  result(1)=updatePar;
  return result;
  
}