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

arma::mat WeightFunc(const arma::mat&rules,const arma::vec&beta,const double&theta,const arma::vec&gamma,const arma::mat&Data,
                     const arma::umat&tL,const arma::umat&tR){
    arma::vec normvec=arma::normpdf(rules.col(0));
    arma::vec weightvec=rules.col(1)%exp(pow(rules.col(0),2))%normvec;
    arma::mat Oriweight=repmat(trans(weightvec),Data.n_rows,1);
    arma::mat index=sort(unique(Data.col(0)));
    int n=index.n_rows;
    
}

// [[Rcpp::export]]
arma::umat testfind(const arma::mat&testmat){
    arma::umat result=find(testmat.col(0)==2);
    return result;
}

// [[Rcpp::export]]
arma::mat findmat(const arma::mat&testmat){
    arma::umat result=find(testmat.col(0)==1);
    arma::mat final=testmat.rows(result);
    arma::mat B = arma::prod(final,0);
    return B;
}















































