// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]

arma::umat TmatL(const arma::mat&Inspec,const arma::vec&Timepoints){
    int n=Inspec.n_rows;
    int m=Timepoints.n_elem;
    arma::umat result(m,n);
    for(int i=0;i<n;i++){
        result.col(i)=Timepoints<=Inspec(i,0);
    }
    return result;
}
// [[Rcpp::export]]
arma::umat TmatR(const arma::mat&Inspec,const arma::vec&Timepoints){
    int n=Inspec.n_rows;
    int m=Timepoints.n_elem;
    arma::umat result(m,n);
    for(int i=0;i<n;i++){
        result.col(i)=Timepoints<=Inspec(i,1);
    }
    return result;
}

// [[Rcpp::export]]

arma::umat TmatLR(const arma::mat&Inspec,const arma::vec&Timepoints){
    int n=Inspec.n_rows;
    int m=Timepoints.n_elem;
    arma::umat result(m,n);
    for(int i=0;i<n;i++){
        result.col(i)=Timepoints<=Inspec(i,1)&&Timepoints>Inspec(i,0);
    }
    return result;
}

// [[Rcpp::export]]
arma::mat elementwise_pow(const arma::mat&A,const arma::mat&p){
    int I=A.n_rows;
    int J=A.n_cols;
    arma::vec result=arma::zeros(I,J);
    
    for(int i=0;i<I;i++){
        for(int j=0;j<J;j++){
            result(i,j)=std::pow(A(i,j),p(i,j));
        }
        
    }
    return result;
}

// [[Rcpp::export]]

arma::mat Likelihood_rules(const arma::mat&rules,const arma::vec&beta,const double&theta,const arma::vec&gamma,const arma::mat&Data,
                    const arma::umat&tL,const arma::umat&tR){
    arma::mat index=sort(unique(Data.col(0)));
    int n=index.n_rows;
    int order=rules.n_rows;
    arma::mat covariate=Data.submat(0,5,Data.n_rows-1,Data.n_cols-1);
    arma::mat betaX=covariate*beta;
    arma::mat ParPart=exp(repmat(betaX,1,order)+repmat(trans(theta*rules.col(0)),Data.n_rows,1));
    arma::mat LambdaR=trans(trans(gamma)*tR);
    arma::mat LambdaL=trans(trans(gamma)*tL);
    
    arma::mat expL=exp(-repmat(LambdaL,1,order)%ParPart);
    arma::mat expR=exp(-repmat(LambdaR,1,order)%ParPart);
    arma::mat Likelihood=elementwise_pow(expL,repmat(Data.col(3),1,order))%
        elementwise_pow(expL-expR,1-repmat(Data.col(3),1,order)-repmat(Data.col(4),1,order))%
        elementwise_pow(expR,repmat(Data.col(4),1,order));
    return Likelihood;
}



// [[Rcpp::export]]

arma::mat WeightFunc(const arma::mat&rules,const arma::mat&beta,const double&theta,const arma::mat&gamma,const arma::mat&Data,
                     const arma::umat&tL,const arma::umat&tR){
    arma::vec normvec=arma::normpdf(rules.col(0));
    arma::vec weightvec=rules.col(1)%exp(pow(rules.col(0),2))%normvec;
    int order=rules.n_rows;
    // arma::mat Oriweight=repmat(trans(weightvec),Data.n_rows,1);
    arma::mat index=sort(unique(Data.col(0)));
    int n=index.n_rows;
    arma::mat LikelihoodMat=Likelihood_rules(rules,beta,theta,gamma,Data,tL,tR);
    arma::mat weightMat=arma::zeros(n,order);
    for(int i=1;i<(n+1);i++){
        arma::umat rowindices=find(Data.col(0)==i);
        arma::mat Eachcluster=LikelihoodMat.rows(rowindices);
        arma::mat Rowprod = arma::prod(Eachcluster,0);
        arma::mat Finalweight=Rowprod%trans(weightvec)/accu(Rowprod%trans(weightvec));
        weightMat.row(i-1)=Finalweight;
    }
    return weightMat;
}







// [[Rcpp::export]]
arma::field<arma::mat> Updateonce(const arma::mat&rules,const arma::mat&beta,const double&theta,const arma::mat&gamma,const arma::mat&Data,
                     const arma::umat&tL,const arma::umat&tR,const arma::umat&tLR){
    arma::mat WeightMat=WeightFunc(rules,beta,theta,gamma,Data,tL,tR);
    arma::mat index=sort(unique(Data.col(0)));
    int n=index.n_rows;
    int order=rules.n_rows;
    int m=gamma.n_rows;
    arma::mat covariate=Data.submat(0,5,Data.n_rows-1,Data.n_cols-1);
    arma::mat betaX=covariate*beta;
    arma::mat ParPart=exp(repmat(betaX,1,order)+repmat(trans(theta*rules.col(0)),n,1));
    arma::mat LambdaR=repmat(trans(trans(gamma)*tR),1,order);
    arma::mat LambdaL=repmat(trans(trans(gamma)*tL),1,order);
    arma::mat LambdaLR=repmat(trans(trans(gamma)*tLR),1,order);
    
    arma::mat expL=LambdaL%ParPart;
    arma::mat expR=LambdaR%ParPart;
    arma::mat expLR=LambdaLR%ParPart;
    
    arma::mat CubegammaR=repmat(trans(pow(gamma,3)),n,1)%trans(tR);
    arma::mat CubegammaLR=repmat(trans(pow(gamma,3)),n,1)%trans(tLR);
    arma::mat InversgammaL=repmat(trans(pow(gamma,-1)),n,1)%trans(tL);
        
    arma::mat Term1=trans(repmat(Data.col(3),1,order)%pow((1+expR),-1)%pow((LambdaR),-1));
    arma::mat Term2=trans((1-repmat(Data.col(3),1,order)-repmat(Data.col(4),1,order))%pow((1+expLR),-1)%pow((LambdaLR),-1));
    arma::mat Term3=trans((1-repmat(Data.col(4),1,order))%ParPart);
    arma::mat Amat(n,m),Bmat(n,m);
    
    
    arma::mat FirstBigMat=repmat(Data.col(3),1,order)%pow((1+expR),-1)+
                      (1-repmat(Data.col(3),1,order)-repmat(Data.col(4),1,order))%pow((1+expLR),-1)-
                      (1-repmat(Data.col(3),1,order))%expL;
    
    arma::mat SecondBigMat=repmat(Data.col(3),1,order)%pow((1+expR),-1)+
        (1-repmat(Data.col(3),1,order)-repmat(Data.col(4),1,order))%pow((1+expLR),-1)+
        (1-repmat(Data.col(3),1,order))%expL;
        
        arma::mat FirstDerivMat(n,beta.n_rows+1);
        arma::mat SecondDeriv=arma::zeros(beta.n_rows+1,beta.n_rows+1);
        arma::mat SecondDerivMat=arma::zeros(beta.n_rows+1,beta.n_rows+1);
    for(int i=0;i<n;i++){
        arma::umat rowindices=find(Data.col(0)==(i+1));
        arma::mat AclusterI=Term1.cols(rowindices)*CubegammaR.rows(rowindices)+
                           Term2.cols(rowindices)*CubegammaLR.rows(rowindices);
        arma::mat BclusterI=Term3.cols(rowindices)*InversgammaL.rows(rowindices);
        Amat.row(i)=WeightMat.row(i)*AclusterI;
        Bmat.row(i)=WeightMat.row(i)*BclusterI;
        
        arma::mat FirstBigMatij=FirstBigMat.rows(rowindices);
        arma::mat SecondBigMatij=SecondBigMat.rows(rowindices);
        arma::mat FirstDerivI(order,beta.n_rows+1);
        for(int r=0;r<order;r++){
            
            arma::mat Xij=covariate.rows(rowindices);
            arma::vec onesvec=arma::ones(Xij.n_rows);
            arma::mat Wijr=join_rows(Xij,rules(r,0)*onesvec*theta);
            FirstDerivI.row(r)=FirstBigMatij.row(r)*Wijr;
            arma::mat SecondDerivI=arma::zeros(beta.n_rows+1,beta.n_rows+1);
            
            for(int j=0;j<Xij.n_rows;j++){
                SecondDerivI=SecondDerivI+SecondBigMatij(j,r)*(trans(Wijr.row(j))*Wijr.row(j));
            }
            SecondDeriv=SecondDeriv-2*WeightMat(i,r)*SecondDerivI;
            
            
        }
        FirstDerivMat.row(i)=WeightMat.row(i)*FirstDerivI;
        
        
    }
    arma::mat FirstDeriv=trans(sum(FirstDerivMat,0));
    SecondDeriv(beta.n_rows,beta.n_rows)=SecondDeriv(beta.n_rows,beta.n_rows)+FirstDeriv(beta.n_rows);
    arma::mat updategamma=pow((trans(sum(Amat,0))%pow(trans(sum(Bmat,0)),-1)),1/4);
    arma::vec lastpar(beta.n_rows+1);
    lastpar.subvec(0,beta.n_rows-1)=beta;
    lastpar(beta.n_rows)=std::log(theta);

    arma::vec updatePar=lastpar-solve(SecondDeriv,FirstDeriv);
    updatePar(beta.n_rows)=std::exp(updatePar(beta.n_rows));
    
    
    arma::field<arma::mat> result(2);
    result(0)=updategamma;
    result(1)=updatePar;
    return result;
    
}



// [[Rcpp::export]]
arma::field<arma::mat> MainFunc(const arma::mat&Data,const arma::mat&rules){
    int betadim=Data.n_cols-5;
    arma::mat tktemp=sort(unique(join_rows(Data.col(1),Data.col(2))));
    arma::vec tk=tktemp.col(0);
    tk=tk.subvec(1,tk.n_elem-2);
    arma::mat L=Data.col(1);
    arma::mat R=Data.col(2);
    
    arma::umat tL=TmatL(join_rows(Data.col(1),Data.col(2)),tk);
    arma::umat tR=TmatR(join_rows(Data.col(1),Data.col(2)),tk);
    arma::umat tLR=TmatLR(join_rows(Data.col(1),Data.col(2)),tk);
    arma::mat beta0=arma::ones(betadim,1)*0.5;
    arma::mat gamma=arma::ones(m,0.1);
    double theta=0.5;
    
}

































