#Group 6
#Aayushi Pandey (2038981/u848144)
#Shreshtha Sharma (2037590/u434770)
#Nasharetty Wiesken (2011820/u876834)
#Carolyn Landbrug (2045788/u901160)


### Title:    Routines for MI-Based Predictions
### Author:   Kyle M. Lang
### Created:  2020-02-06
### Modified: 2020-03-03

source("studentFunctions.R")

###--------------------------------------------------------------------------###

.lambda <- function(x) (x + 1) / (x + 3)

###--------------------------------------------------------------------------###

.fmi <- function(m, b, t) (1 + (1 / m)) * (b / t)

###--------------------------------------------------------------------------###

.miDf <- function(m, b, t, dfCom) {
    fmi   <- .fmi(m, b, t)
    df0   <- (m - 1) * (1 / fmi^2)
    dfObs <- .lambda(dfCom) * dfCom * (1 - fmi)
    
    df0 / (1 + (df0 / dfObs))
}

###--------------------------------------------------------------------------###

splitImps <- function(imps, index) {
    ## Split the imputed dataset by 'index':
    tmp <- lapply(imps, split, f = index)

    ## Create lists of subsets of imputed data:
    out <- list()
    for(i in 1 : length(tmp[[1]]))
        out[[i]] <- lapply(tmp, "[[", x = i)

    ## Name output lists:
    names(out) <- names(tmp[[1]])

    out
}

###--------------------------------------------------------------------------###

predictMi <- function(fits, newData, interval = NA, pooled = TRUE) {
    ## Generate imputation-specific predictions:
    predStuff <- list()
    for(m in 1 : length(fits))
        predStuff[[m]] <-
            predict(fits[[m]], newdata = newData[[m]], se.fit = TRUE)
    
    ## Extract components from prediction objects:
    preds  <- lapply(predStuff, "[[", x = "fit")
    seMat  <- sapply(predStuff, "[[", x = "se.fit")
    scales <- sapply(predStuff, "[[", x = "residual.scale")

    ## Return early when imputation-specific predictions are requested:
    if(!pooled)
        return(preds)
    else
        predMat <- do.call(cbind, preds)
    
    ## Compute pooled predictions:
    preds <- rowMeans(predMat)
    
    ## Compute within-imputation variance:
    if(is.na(interval) | interval == "confidence")
        w <- rowMeans(seMat^2)
    else if(interval == "prediction")
        w <- colMeans(
            apply(seMat^2, 1, function(x, y) x + y, y = scales^2)
        )
    else
        stop(paste0("'",
                    interval,
                    "' is not a valid argument for 'interval'. Should be either 'confidence' or 'prediction'."
                    )
             )

    ## Compute between-imputation variance:
    b <- apply(predMat, 1, var)

    ## Compute total variance:
    t <- w + b + (b / length(fits))
    
    ## Compute df:
    df <- .miDf(length(fits), b, t, predStuff[[1]]$df)

    ## Compute the FMIs:
    out <- cbind(est = preds, fmi = .fmi(length(fits), b, t))

    ## Compute CIs and aggregate output:
    if(!is.na(interval)) {
        moe <- qt(0.975, df) * sqrt(t)
        out <- cbind(out, lwr = preds - moe, upr = preds + moe)
    }
    
    out
}

###--------------------------------------------------------------------------###

mseMi <- function(fits, newData) {
    ## Generate imputation-specific predictions:
    preds <- predictMi(fits = fits, newData = newData, pooled = FALSE)
    
    ## Compute imputation-specific MSEs:
    mse <- c()
    for(m in 1 : length(fits))
        mse[m] <- MSE(y_pred = preds[[m]],
                      y_true = newData[[m]][[dvName(fits[[1]])]]
                      )
    
    ## Return the aggregate MSE:
    mean(mse)
}

