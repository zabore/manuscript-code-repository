################################################################################
# Function to perform a clustering simulation
#
#' @param M number of etiologically distinct disease subtypes
#' @param N sample size
#' @param pi vector of control and subtype prevalences
#' i.e. c(pi0, pi1, pi2, pi3) for controls and 3 subtypes
#' @param P number of risk factors
#' @param mu_m P x M matrix of risk factor means for the cases, by default all
#' risk factors have mean 0 for control subjects
#' @param K_A number of tumor markers that are related to the etiologically 
#' distinct disease subtypes
#' @param lambda_Am K_A x M matrix of means for tumor markers related to the
#' etiologically distinct disease subtypes
#' @param J number of non-etiologically distinct subtypes
#' @param K_B number of tumor markers that are urelated to the etiologically 
#' distinct disease subtypes but have some other structure
#' @param lambda_Bj K_B x J matrix of means for tumor markers unrelated to the
#' etiologically distinct disease subtypes
#' @param K_C number of additional tumor markers that are urelated to the 
#' etiologically distinct disease subtypes
#' @param nstart number of random starts for k-means clustering
#' @param reducedim if TRUE then do dimension reduction using univariable D
################################################################################

devtools::install_github("zabore/riskclustr")
library(riskclustr)
suppressMessages(library(mlogit))
library(gtools)

eh_cluster_sim <- function(M, N, pi, P, mu_m, K_A, lambda_Am, 
                           J, K_B, lambda_Bj, K_C, nstart, reducedim) {
    
    # Set null return objects
    max_d <- NULL
    min_misclass <- NULL
    unique_g <- NULL

    max_d_univd <- NULL
    min_misclass_univd <- NULL
    n_ka_sel <- NULL
    n_kb_sel <- NULL
    n_kc_sel <- NULL
    
    # generate etiologically distinct subtypes indicator matrix
    tAcls <- rep(1:M, N * pi[2:(M + 1)])
    tA <- matrix(0, nrow = length(tAcls), ncol = length(unique(tAcls)))
    tA[cbind(seq_along(tAcls), tAcls)] <- 1
    tA <- rbind(matrix(rep(0, N * pi[1] * M), N * pi[1], M), tA)
    
    # generate risk factor matrix
    x <- matrix(rnorm(N * P), N, P) + tA %*% t(mu_m)
    colnames(x) <- paste0("x", seq(1:ncol(x)))
    
    # matrix of EH-related tumor marker data (N x K_A)
    tmA <- matrix(rnorm(N * K_A), N, K_A) + tA %*% t(lambda_Am)
    
    # generate unrelated subtypes indicator
    tBcls <- as.numeric(cut(runif(N * sum(pi[2:4])), seq(0, 1, by = 1 / J)))
    tB <- matrix(0, nrow = length(tBcls), ncol = length(unique(tBcls)))
    tB[cbind(seq_along(tBcls), tBcls)] <- 1
    tB <- rbind(matrix(rep(0, N * pi[1] * J), N * pi[1], J), tB)
    
    # matrix of EH-unrelated tumor marker data (N x K_B)
    if(K_B == 0) tmB <- NULL else 
        tmB <- matrix(rnorm(N * K_B), N, K_B) + tB %*% t(lambda_Bj)
    
    # matrix of additional noise tumor marker data (N x K_C)
    if(K_C == 0) tmC <- NULL else
        tmC <- matrix(rnorm(N * K_C), N, K_C)
    
    # matrix of all tumor marker data for CASES ONLY with numeric colnames
    y <- tail(cbind(tmA, tmB, tmC), -N * pi[1])
    colnames(y) <- seq(1:ncol(y))
    
    # polytomous model formula
    mform <- mFormula(as.formula(c("cls ~ 1 |", 
                                   paste(colnames(x), collapse = "+"))))
    
    # k-means clustering on the combined true and noise tumor marker data
    # only among the cases
    kres <- sapply(1:nstart, function(l) {
        kk <- kmeans(y, centers = M, algorithm = "MacQueen", iter.max = 30)
        kclust <- kk$cluster
        kg <- kk$betweenss / kk$totss
        return(list(kclust, kg))
    })
    kres_out <- split(kres, 1:2)
    kclust_out <- matrix(unlist(kres_out[[1]]), (N * (1 - pi[1])), nstart)
    kg_out <- unlist(kres_out[[2]])
    tab_kg_out <- data.frame(table(kg_out))
    unique_g <- length(unique(kg_out))
    
    # calculate D for each solution
    d_kres <- sapply(1:nstart, function(l) {
        dat <- data.frame(cbind(cls = c(rep(0, N * pi[1]), kclust_out[, l]), 
                                x), row.names = NULL)
        if(any(table(dat$cls) < 20)) NA else
            dest(mform, "cls", M, dat)
    })
    
    #index of selected sol
    max_res <- which(d_kres == max(d_kres[kg_out >= mean(kg_out)], 
                                   na.rm = TRUE))[1]
    
    # save the results
    max_d <- d_kres[max_res]
    min_misclass <- minmc(table(kclust_out[, max_res], tAcls))$minmisclass
    
    # Dimension reduction (if reducedim = TRUE)
    if(reducedim) {
        
        # Need binary tumor marker data - use median split
        ybin <- apply(y, 2, function(x) ifelse(x < median(x), 1, 2))
        
        # Estimate D for each binary tumor marker
        d_ybin <- sapply(1:ncol(ybin), function(l) {
            dat <- data.frame(cbind(cls = c(rep(0, N * pi[1]), ybin[, l]),
                                    x), row.names = NULL)
            dest(mform, "cls", 2, dat)
        })
        
        # Make an ordered dataframe of the resulting univariate D
        res <- data.frame(d_ybin[order(d_ybin, decreasing = TRUE)])
        rownames(res) <- order(d_ybin, decreasing = TRUE)
        colnames(res) <- "d_ybin"
        
        for(j in seq(5, (K_A + K_B + K_C), 5)){
            
            # Create a matrix with the selected genes
            ysel <- y[, rownames(res)[1:j]]
            
            # Now perform k-means clustering on the reduced tumor marker data
            kres_univd <- sapply(1:nstart, function(l) {
                kk <- kmeans(ysel, centers = M, algorithm = "MacQueen", 
                             iter.max = 30)
                kclust <- kk$cluster
                kg <- kk$betweenss / kk$totss
                return(list(kclust, kg))
            })
            kres_univd_out <- split(kres_univd, 1:2)
            kclust_univd_out <- matrix(unlist(kres_univd_out[[1]]),
                                       nrow = (N * (1 - pi[1])), ncol = nstart)
            kg_univd_out <- unlist(kres_univd_out[[2]])
            
            # calculate D for each solution
            d_kres_univd <- sapply(1:nstart, function(k) {
                dat <- data.frame(cbind(cls = c(rep(0, N * pi[1]), 
                                                kclust_univd_out[, k]), x), 
                                  row.names = NULL)
                if(any(table(dat$cls) < 20)) NA else
                    dest(mform, "cls", M, dat)
            })
            
            #index of selected sol
            max_univd_res <- which(d_kres_univd == max(d_kres_univd[
                kg_univd_out >= mean(kg_univd_out)], na.rm = TRUE))[1]
            
            # save the results
            max_d_univd <- c(max_d_univd, d_kres_univd[max_univd_res])
            min_misclass_univd <- c(min_misclass_univd, minmc(table(
                kclust_univd_out[, max_univd_res], tAcls))$minmisclass)
            n_ka_sel <- c(n_ka_sel, sum(colnames(ysel) %in% 1:K_A))
            n_kb_sel <- c(n_kb_sel, sum(colnames(ysel) %in% (K_A + 1):
                                            (K_A + K_B)))
            n_kc_sel <- c(n_kc_sel, sum(colnames(ysel) %in% (K_B + 1):
                                            (K_A + K_B + K_C)))
        }
    }
    
    return(list(max_d = max_d, 
                min_misclass = min_misclass, 
                unique_g = unique_g, 
                max_d_univd = max_d_univd,
                min_misclass_univd = min_misclass_univd,
                n_ka_sel = n_ka_sel, 
                n_kb_sel = n_kb_sel, 
                n_kc_sel = n_kc_sel))
}