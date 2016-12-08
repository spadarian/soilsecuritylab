// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//'Calculate distances
//'
//'@param x matrix
//'@param y matrix
//'@param w matrix
//'@return distance matrix
// [[Rcpp::export]]
arma::mat mahaldist(arma::mat x,arma::mat y, arma::mat w) {

    arma::mat ans(x.n_rows,y.n_rows);

    for (int i = 0; i < y.n_rows; i++) {
         arma::mat Dc = repmat(y.row(i),x.n_rows,1) - x;
         ans.col(i) = real(sum((Dc*w)%conj(Dc),1));
    }
    return ans;
}

//'Run fuzzy clustering algorithm
//'
//'@param alpha temp.
//'@param data a \code{data.frame} with numerical data. Each column corresponds to a dimension of a \code{n}-dimensional clustering domain.
//'@param phi a numerical value representing the degree of fuzziness or overlap of the generated clusters.
//'@param nclass numerical (integer) vector with the number of cluster to generate.
//'@param disttype numerical value representing the type of distance to calculate for the clustering. Possible values are: 1 (Euclidean), 2 (Diagonal), and 3 (default, Mahalanobis).
//'@param maxiter maximum number of iterations.
//'@param toldif temp.
//'@param exp_eg the expected fraction of extragrades. If not provided, it is assumed to be dependant on the number of clusters (\code{1/(1+nclusters)}).
//'@param optim temp.
//'@return \code{FuzzyClusterGroup}
// [[Rcpp::export(name=".fuzzy_extragrades_C")]]
List fuzzy_extragrades_C(double alpha, NumericMatrix data, double phi, int nclass, int disttype, int maxiter, double toldif, double exp_eg, bool optim) {
    // Rcout << "The value is " << ans << std::endl;
    Environment pkg_stats("package:stats");
    Function rnorm_ = pkg_stats["rnorm"];
    Function cov = pkg_stats["cov"];
    Environment pkg_base("package:base");
    Function solve = pkg_base["solve"];
    Function round_ = pkg_base["round"];

    int ndata = data.nrow(), ndim = data.ncol();
    arma::mat dist(ndata,nclass);

    arma::mat U(ndata,nclass);
    for (int i = 0; i < ndata*nclass; i++) {
        U[i] = rnorm(1,1.0/nclass,0.01)[0];
    }

    for (int i = 0; i < ndata; i++) {
        U.row(i) = U.row(i)/sum(U.row(i));
    }

    arma::mat W;

    if (disttype == 1){
        W = arma::mat(ndim, ndim);
        W.eye();
    }
    if (disttype == 2){
        W = arma::mat(ndim, ndim);
        W.eye();
        W = W % as<arma::mat>(cov(data));
    }
    if (disttype == 3){
        W = as<arma::mat>(solve(cov(data)));
    }

    double obj = 0.0;

    arma::vec Ue(ndata);
    Ue = as<arma::vec>(round_(Ue.fill(1.0),15)) - as<arma::vec>(round_(sum(U,1),15));
    arma::mat uphi = pow(U,phi);
    arma::mat uephi = pow(Ue,phi);
    double a1 = (1-alpha)/alpha;

///////// Initialice

    arma::mat c1 = uphi.t() * as<arma::mat>(data);
    arma::mat t1 = sum(uphi);
    t1 = repmat(t1.t(),1,ndim);
    arma::mat centroids = c1/t1;


//////// Calculate distances to centroids

    if (disttype == 1) { // Euclidian

        for (int i = 0; i < centroids.n_rows; i++) {
             dist.col(i) = sqrt(sum(pow(as<arma::mat>(data) - repmat(centroids.row(i),ndata,1),2),1));
        }

    } else { // Diagonal and Mahalanobis
      dist = sqrt(mahaldist(as<arma::mat>(data),centroids,W));
    }

/////// Main loop
    arma::mat ufi;
    int iters;
    for (iters = 0; iters < maxiter; iters++) {

        // Calculate centroids

        ufi = uphi-(arma::mat(ndata,nclass).fill(a1)%pow(dist,-4)%repmat(uephi,1,nclass));
        c1 = ufi.t()*as<arma::mat>(data);
        t1 = sum(ufi);
        t1 = repmat(t1.t(),1,ndim);
        centroids = c1/t1;

        // Calculate distances to centroids
        if (disttype == 1) { // Euclidian

            for (int i = 0; i < centroids.n_rows; i++) {
                dist.col(i) = sqrt(sum(pow(as<arma::mat>(data) - repmat(centroids.row(i),ndata,1),2),1));
            }

        } else { // Diagonal and Mahalanobis
            dist = sqrt(mahaldist(as<arma::mat>(data),centroids,W));
        }

        // Save previous iteration
        arma::mat U_old = U;
        double obj_old = obj;

        // Calculate new membership arma::matrix

        // if (iters == 0) {
        //     Rcout << "dist:\n " << dist << std::endl;
        // }

        arma::mat tmp = pow(dist,-2/(phi-1));
        arma::mat tm2 = pow(dist,-2);
        arma::mat s2 = pow(a1*sum(tm2,1),-1/(phi-1));

        t1 = sum(tmp,1);
        arma::mat t2 = repmat(t1,1,nclass)+repmat(s2,1,nclass);
        U = tmp/t2;
        Ue = as<arma::vec>(round_(Ue.fill(1.0),15)) - as<arma::vec>(round_(sum(U,1),15));
        uphi = pow(U,phi);
        uephi = pow(Ue,phi);

        // Calculate obj function
        arma::mat o1 = pow(dist,2)%uphi;
        arma::mat d2 = pow(dist,-2);
        arma::mat o2 = uephi%sum(d2,1);
        obj = alpha*accu(o1)+(1-alpha)*accu(o2);

        // Check for convergence
        double dif = obj_old-obj;
        arma::mat difU = sqrt(pow(U-U_old,2));
        double Udif = accu(difU);
        if (dif < toldif & Udif < toldif) {
            break;
        }

    }

    ufi = uphi-(arma::mat(ndata,nclass).fill(a1)%pow(dist,-4)%repmat(uephi,1,nclass));
    c1 = ufi.t()*as<arma::mat>(data);
    t1 = sum(ufi);
    t1 = repmat(t1.t(),1,ndim);
    centroids = c1/t1;

    U.resize(U.n_rows,U.n_cols+1);
    U.col(U.n_cols-1) = Ue;

    // colnames(U_end) <- c(1:nclass,'Eg')
    // attributes(U_end)$hard_clust <- colnames(U_end)[apply(U_end,1,which.max)]
    double Ue_mean = mean(Ue);

    arma::vec CI = zeros<arma::vec>(U.n_rows);

    for (int i = 0; i < U.n_rows; i++) {
        arma::vec tmp = sort(U.row(i),"descend").t();
        CI.row(i) = tmp.row(1)/tmp.row(0);
    }

    // attributes(U_end)$CI <- CI

    List resp;
    if (optim){
        resp["Ue_mean-Ue_req"] = Ue_mean-exp_eg;
    } else {
        resp["U"] = U;
        resp["W"] = W;
        resp["centroids"] = centroids;
        resp["Ue_mean-Ue_req"] = Ue_mean-exp_eg;
        resp["iterations"] = iters;
    }

    return resp;

}
