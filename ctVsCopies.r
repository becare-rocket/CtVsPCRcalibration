## Ct Vs. Copies Paper R Code
## Paper Title: Possible Estimation Bias due to High PCR Cycle Threshold Cutoff (Ct)
## version 2021-05-26: For 25 May 2021
## Paper available DOI: 10.13140/RG.2.2.17180.82568
## Code and Paper Author: Michael Halem 
## Code (c) Michael Halem 2021
## This code may be freely used or reused as long as it is attributed to the author.
####################################################################################
## Startup R and source this file, i.e. > source('ctVsCopies.r')
## Then
## > makeEstimatedCutoffs(bUseBig=T) ## Will draw figure 1
## > makeEstimatedCutoffs(bUseBig=F) ## Will draw figure 2
## The bDebug switch will show more detail in various calculations that are not used in the paper
## but that you may be interested in.
## The bReproduceWalker was my attempt to replicate the regressions in Walker as a corroboration.
## The bOSN will display the ORF1ab, the S AB and the N antibody graphs.  This will make it too busy.
## Questions: email michael.halem@becare.net or michael.halem@rocketpwr.com

####################################################################################
## Utilities, constants, and other one-time setup
dateRefPosix = as.POSIXct('2019-12-31 OO:OO', tz='EST')
DateRefDate = as.Date('2019-12-31', tz='EST')
options(width=266, digits=5)
options(stringsAsFactors=F)
printf<-function(...) cat(sprintf(... , sep=''))

###################################################################################

makeEstimatedCutoffs<-function(
                               fName = 'copiesVsCtCal2.csv', ##Change if you move the data elsewhere
                               rho=0, ##correlation of all 3 Ct observations to each other
                               bUseBig = TRUE,
                               bDebug=FALSE,
                               bDisplayOSN=FALSE,
                               bReproduceWalker=FALSE, ##Reproduce Walker et al's regressions (the above doi)
                               nSim=500000,  ##used when using rnorm to create standard errors in X data
                               bRetModel=FALSE
                               )
{
    calTab = read.csv(fName)  ##This csv from CopiesVsCtCalibratioin.xlsx
    ci95sigmas = qnorm(0.975)
    
    nPoints = dim(calTab)[1]
    calTab$Ct.Mean = with(calTab, (Ct.ORF1ab + Ct.N + Ct.S)/3 )
    calTab$SE.ORF1ab=calTab$CI95.ORF1ab/ci95sigmas
    calTab$SE.S=calTab$CI95.S/ci95sigmas
    calTab$SE.N=calTab$CI95.N/ci95sigmas

    if (bReproduceWalker)
    { #walker uses non-linear S01-S08 scales
        calTab$scaleWalker=1:8
        modORF=lm(Ct.ORF1ab~scaleWalker, data=calTab)
        modS=lm(Ct.S~scaleWalker, data=calTab)
        modSW=lm(Ct.S~scaleWalker, data=calTab, weights = 1/SE.S^2)
        modN=lm(Ct.N~scaleWalker, data=calTab)
        print(summary(modS))
        print(summary(modSW))
        with(calTab, plot(Ct.S~scaleWalker, data=calTab))
        for (i in 1:dim(calTab)[1])
        with(calTab[i,], points(y=c(Ct.S+CI95.S, Ct.S-CI95.S),
                                x = rep(scaleWalker,2), col="black", type='l'))
        abline(coef=modS$coefficients)
        r2=1-simVar(x=predict(modS)-calTab$Ct.S, SE.Xi=calTab$SE.S)/simVar(x=calTab$Ct.S,SE.Xi=calTab$SE.S)
        r2Adj = 1 - (1 - r2) * ((nPoints-1)/(nPoints - 2))
        printf("Walker S r2 =%f r2adj=%f\n", r2, r2Adj)
        
        return(list(modSW=modSW, modS=modS, calTab=calTab))     
    }

    
    calTab$SE.Mean = with(calTab, 1/3*sqrt(SE.ORF1ab^2 + SE.S^2 + SE.N^2 + 2*rho*(SE.ORF1ab*SE.S+SE.ORF1ab*SE.N+SE.N*SE.S)) )
    calTab$CI95.Mean = ci95sigmas*calTab$SE.Mean
    mod.ORF1ab = lm(Ct.ORF1ab~log10(copies.ml), data=calTab, weights=1/SE.ORF1ab^2) ##Weighted linear regression    
    mod.S = lm(Ct.S~log10(copies.ml), data=calTab, weights=1/SE.S^2) ##Takes into account wider SE on some data points...
    mod.N = lm(Ct.N~log10(copies.ml), data=calTab, weights=1/SE.N^2)
    mod.Mean = lm(Ct.Mean~log10(copies.ml), data=calTab, weights=1/SE.Mean^2)

    
    uselessLines = c(1:2,4:9,13:15)
    szSummary = character(0)
    if (bDebug)
        {
            szSummary = capture.output(summary(mod.ORF1ab))[-uselessLines]
            szSummary = c(szSummary, capture.output(summary(mod.S))[-uselessLines])
            szSummary = c(szSummary, capture.output(summary(mod.N))[-uselessLines])
        }
    szSummary = c(szSummary, capture.output(summary(mod.Mean))[-uselessLines])
    ##szSummary = szSummary[nchar(szSummary) != 0]
    printf("%s\n", szSummary)
    

    calTabExt = calTab[1:3,]
    calTabExt$copies.ml = c(5,2,1)
    calTabExt[,-1] = NA
    
    calTabExt$Ct.ORF1ab = predict(mod.ORF1ab, newdata=calTabExt)
    calTabExt$Ct.S = predict(mod.S, newdata=calTabExt)
    calTabExt$Ct.N = predict(mod.N, newdata=calTabExt)    
    calTabExt$Ct.Mean = predict(mod.Mean, newdata=calTabExt)

    ##The ordinary model sigmas are wrong because they are biased upwards by 1/SE^2 (so are the rSquared)
    ##An approximate adjustment is to increase sigma from OLS with NO weights by sqrt(n/df * mean(SE^2))
    ##This is because the regression estimate now INCLUDES measurement standard errors (!!!)
    df = nPoints - 2
    ## This is corrected using simVar 500K steps per point estimate

    if (bDebug)
        {
            mod.ORF1ab$sigma = sqrt( (nPoints-1)/df * simVar(x=predict(mod.ORF1ab)-calTab$Ct.ORF1ab, SE.Xi=calTab$SE.ORF1ab, nSim=nSim) ) 
            mod.S$sigma = sqrt( (nPoints-1)/df * simVar(x=predict(mod.S)-calTab$Ct.S, SE.Xi=calTab$SE.S, nSim=nSim) ) 
            mod.N$sigma = sqrt( (nPoints-1)/df * simVar(x=predict(mod.N)-calTab$Ct.N, SE.Xi=calTab$SE.N, nSim=nSim) )
        }
    mod.Mean$sigma = sqrt( (nPoints-1)/df * simVar(x=predict(mod.Mean)-calTab$Ct.Mean, SE.Xi=calTab$SE.Mean, bDebug=bDebug, nSim=nSim) )
    ##Needed for getting the forecast confidence intervals
    mod.Mean$xmean = mean(log10(calTab$copies.ml))
    mod.Mean$xVarP = var(log10(calTab$copies.ml))*(nPoints - 1)/nPoints


    if (bDebug)
    {
        mod.ORF1ab$sigmaApprox = sqrt(1/df*(
            (nPoints-1)*var(predict(lm(Ct.ORF1ab~log10(copies.ml),data=calTab))-calTab$Ct.ORF1ab)+nPoints*mean(calTab$SE.ORF1ab^2)) )
        mod.S$sigmaApprox = sqrt(1/df*(
            (nPoints-1)*var(predict(lm(Ct.S~log10(copies.ml),data=calTab))-calTab$Ct.S)+nPoints*mean(calTab$SE.S^2)) )
        mod.N$sigmaApprox = sqrt(1/df*(
            (nPoints-1)*var(predict(lm(Ct.N~log10(copies.ml),data=calTab))-calTab$Ct.N)+nPoints*mean(calTab$SE.N^2)) )
        mod.Mean$sigmaApprox = sqrt(1/df*(
            (nPoints-1)*var(predict(lm(Ct.Mean~log10(copies.ml),data=calTab))-calTab$Ct.Mean)+nPoints*mean(calTab$SE.Mean^2)) )

        printf("ORF1ab sigma = %f  sigmaApprox=%f\n", mod.ORF1ab$sigma, mod.ORF1ab$sigmaApprox)
        printf("S      sigma = %f  sigmaApprox=%f\n", mod.S$sigma, mod.S$sigmaApprox)
        printf("N      sigma = %f  sigmaApprox=%f\n", mod.N$sigma, mod.N$sigmaApprox)
        printf("Mean   sigma = %f  sigmaApprox=%f\n", mod.Mean$sigma, mod.Mean$sigmaApprox)
    }


    if (bDebug)
    { ## Paper sigma_WLS_Reg
        sigma.noSE = sqrt(1/df * sum( (predict(mod.Mean) - calTab$Ct.Mean)^2) )
        x = numeric(0)
        for (i in 1:nPoints)
            {
            x = c(x, rnorm(n = nSim, mean = predict(mod.Mean)[i] - calTab$Ct.Mean[i], sd = calTab$SE.Mean[i]))
            }
        xmean = mean(x)
        ss = sum((x - xmean)^2)
            
        sigma.SE = sqrt(1/df * 1/nSim * ss)
        sigma.NoMean = sqrt(1/df * 1/nSim * sum(x^2))
        printf( "xmean=%f nPoints=%d df=%d sigma.noSE = %f sigma.SE = %f sigma.NoMean=%f\n", xmean, nPoints, df, sigma.noSE, sigma.SE, sigma.NoMean)
    }


    if (bUseBig)
    {
        calBig = calTab[,c('copies.ml','Ct.Mean')]
        ##add nSim pseudo points with mean and sd of the original.  Note that optionally the original mean can be dropped...
        for (i in 1:nPoints)
            calBig = rbind(calBig, data.frame(copies.ml=calTab$copies.ml[i], Ct.Mean=rnorm(n=nSim, mean=calTab$Ct.Mean[i], sd=calTab$SE.Mean[i])))
        calBig = calBig[-(1:nPoints),]
        modBig = lm(Ct.Mean~log10(copies.ml), data=calBig)
        print(digits=6, summary(modBig))
        modBig$sigma = summary(modBig)$sigma * sqrt(modBig$df.resid/(nSim*nPoints) * (nPoints/df))
        modBig$xmean = mean(log10(calTab$copies.ml))
        modBig$xVarP = var(log10(calTab$copies.ml))*(nPoints - 1)/nPoints
        calTabExt$Ct.Mean = predict(modBig, newdata=calTabExt)

        printf("modBig$xVarP=%f modBig$xmean=%f modBig$sigma = %f on %d points and %d df\n", modBig$xVarP, modBig$xmean, modBig$sigma, nPoints, df)
        modBig$SigmaDoubleChk =  sqrt( sum((predict(modBig) - calBig$Ct.Mean)^2)/dim(calBig)[1] * nPoints/df)
        printf("modBig$SigmaDoubleChk=%f\n", modBig$SigmaDoubleChk)
    }


    
    calTabExt$SE.ORF1ab = mod.ORF1ab$sigma
    calTabExt$SE.S =  mod.S$sigma
    calTabExt$SE.N =  mod.N$sigma
    calTabExt$SE.Mean = with(mod.Mean, sigma * sqrt(1 + 1/nPoints*(1+ (log10(calTabExt$copies.ml) - xmean)^2/xVarP)))
    if (bUseBig) calTabExt$SE.Mean = with(modBig, sigma * sqrt(1 + 1/nPoints*(1+ (log10(calTabExt$copies.ml) - xmean)^2/xVarP)))

   
    if (bDebug)
        {
            calTabExt$CI95.ORF1ab = ci95sigmas * calTabExt$SE.ORF1ab
            calTabExt$CI95.S = ci95sigmas * calTabExt$SE.S
            calTabExt$CI95.N = ci95sigmas * calTabExt$SE.N
        }
    else
        {
            calTabExt$CI95.ORF1ab = calTabExt$CI95.S = calTabExt$CI95.N = rep(NA, dim(calTabExt)[1])
            calTabExt$SE.ORF1ab = calTabExt$SE.S = calTabExt$SE.N = rep(NA, dim(calTabExt)[1])
        }
            
    calTabExt$CI95.Mean = ci95sigmas * calTabExt$SE.Mean 
    
    calTab$comment=""
    calTabExt$comment="Estimated"
    calTabData = calTab
    calTab = rbind(calTab, calTabExt)

    if (bDebug)
        {
            printf("Copies/mL  ORF1ab                  S Gene                  N Gene                  Mean\n")
            for (i in 1:dim(calTab)[1])
                with(calTab[i,],
                     printf("%7d    %4.1f (CI95 %4.1f - %4.1f) %4.1f (CI95 %4.1f - %4.1f) %4.1f (CI95 %4.1f - %4.1f) %4.1f (CI95 %4.1f - %4.1f) %s\n",
                            copies.ml,
                            Ct.ORF1ab, Ct.ORF1ab-CI95.ORF1ab, Ct.ORF1ab+CI95.ORF1ab,
                            Ct.S, Ct.S-CI95.S, Ct.S+CI95.S,
                            Ct.N, Ct.N-CI95.N, Ct.N+CI95.N,
                            Ct.Mean, Ct.Mean-CI95.Mean,Ct.Mean+CI95.Mean,
                            comment))
        }

    offsetMean = 0
    offsetORF = 0.0
    offsetS = 1/2 * offsetORF
    offsetN = -offsetORF
    colORF = 'darkblue'
    colS='darkred'
    colN='darkgreen'

    
    par(mar=c(5,4,4,2))
    
    with(calTab, plot(xlim=c(-.1,6.1), ylim=c(15,40), x=-1, y=-1, col='black', xaxt='n', yaxt='n',  xaxs='i', yaxs='i',
                      main="Mean (of S, N and ORF1ab) PCR Ct vs Direct Copies per mL"
                    , ylab='Ct (Bars: 95% CI)', xlab='Direct Copies per mL'), lwd=1)

    grid(col='black')

    axis(1, at=0:6, cex.axis=0.88, las=1, mgp=c(3.5,0.7,0),labels=sprintf("%.0f", 10^(0:6)))
    extraXs = log10(c(2,5,50,500))
    axis(1, at=extraXs, cex.axis=0.88, las=1, mgp=c(3.5,0.7,0),labels=sprintf("%.0f", 10^extraXs))

    extraYs = c(15,20,25,30,35,40)
    axis(2, at=extraYs, cex.axis=1.0, las=1, mgp=c(3.5,0.7,0),labels=sprintf("%.0f", extraYs))
    extraYs = c(34,36,37,37,38)
    axis(2, at=extraYs, cex.axis=1.0, las=1, mgp=c(3.5,0.7,0),labels=sprintf("%.0f", extraYs))
    

    

    if (bDisplayOSN)
        {
            with(calTabData, points(y=Ct.ORF1ab,x=log10(copies.ml)+offsetORF, col=colORF))
            with(calTabData, points(y=Ct.S, x=log10(copies.ml)+offsetS, col=colS))
            with(calTabData, points(y=Ct.N, x=log10(copies.ml)+offsetN, col=colN))
 
            with(calTabData, points(y=Ct.ORF1ab+CI95.ORF1ab, x=log10(copies.ml)+offsetORF, col=colORF, pch='_'))
            with(calTabData, points(y=Ct.S+CI95.S, x=log10(copies.ml)+offsetS, col=colS, pch='_'))
            with(calTabData, points(y=Ct.N+CI95.N, x=log10(copies.ml)+offsetN, col=colN, pch='_'))
        }
    with(calTab, points(y=Ct.Mean+CI95.Mean, x=log10(copies.ml)+offsetMean, col='black', pch='_', cex=2.0, lwd=1.5))

    if (bDisplayOSN)
        {
            with(calTabData, points(y=Ct.ORF1ab-CI95.ORF1ab, x=log10(copies.ml)+offsetORF, col=colORF, pch='_'))
            with(calTabData, points(y=Ct.S-CI95.S, x=log10(copies.ml)+offsetS, col=colS, pch='_'))
            with(calTabData, points(y=Ct.N-CI95.N, x=log10(copies.ml)+offsetN, col=colN, pch='_'))
        }
    
    with(calTab, points(y=Ct.Mean-CI95.Mean, x=log10(copies.ml)+offsetMean, col='black', pch='_', cex=2.0, lwd=1.5))

    if (!bUseBig)
        {
            abline(coef=mod.Mean$coefficients, col='black', lwd=1.5)
            mod.Mean$x = 0.01*(0:1000) - 1
            mod.Mean$SEforecast = with(mod.Mean, sigma * sqrt(1 + 1/nPoints*(1+ (x - xmean)^2/xVarP)))
            with(mod.Mean, points(y= +ci95sigmas*SEforecast+coefficients[1]+x*coefficients[2], x=x, col='black', type='l', lty='dashed'))
            with(mod.Mean, points(y= -ci95sigmas*SEforecast+coefficients[1]+x*coefficients[2], x=x, col='black', type='l', lty='dashed'))
        }

    if (bUseBig)
        {
        abline(coef=modBig$coefficients, col='black', lwd=1.5)
        modBig$x = 0.01*(0:1000) - 1
        modBig$SEforecast = with(modBig, sigma * sqrt(1 + 1/nPoints*(1+ (x - xmean)^2/xVarP)))
        with(modBig, points(y= +ci95sigmas*SEforecast+coefficients[1]+x*coefficients[2], x=x, col='black', type='l', lty='dashed'))
        with(modBig, points(y= -ci95sigmas*SEforecast+coefficients[1]+x*coefficients[2], x=x, col='black', type='l', lty='dashed'))
        }
    

    offset = 0.005
    for (i in 1:dim(calTabData)[1])
    {
        if (bDisplayOSN)
            {
                with(calTabData[i,], points(y=c(Ct.ORF1ab+CI95.ORF1ab, Ct.ORF1ab-CI95.ORF1ab),
                                            x = rep(log10(copies.ml)+offsetORF,2), col=colORF, type='l'))
                with(calTabData[i,], points(y=c(Ct.S+CI95.S, Ct.S-CI95.S),
                                            x = rep(log10(copies.ml)+offsetS,2), col=colS, type='l'))
                with(calTabData[i,], points(y=c(Ct.N+CI95.N, Ct.N-CI95.N),
                                            x = rep(log10(copies.ml)+offsetN,2), col=colN, type='l'))
            }
        for (i in 1:dim(calTab)[1])
            with(calTabData[i,], points(y=c(Ct.Mean+CI95.Mean, Ct.Mean-CI95.Mean),                               
                                    x = rep(log10(copies.ml)+offsetMean,2), col='black', type='l', lwd=1.5))
    }
        
    with(calTabData, points(y=Ct.Mean, x=log10(copies.ml)+offsetMean, col='black', cex=1.5))

    arrows(x0=log10(5), y0=27, x1=log10(1.2), y1=35, code=2, length=0.1)
    arrows(x0=log10(5), y0=27, x1=log10(2), y1=34, code=2, length=0.1)
    arrows(x0=log10(5), y0=27, x1=log10(5), y1=33, code=2, length=0.1)
    text(x=log10(1), y=26, pos=4,  "Regression Estimates")

    text(x=log10(1000), y=37, pos=4, "Data Points from Walker et al.")
    arrows(x0=log10(10000), y0=36, x1=log10(100), y1=32, code=2, length=0.1)
    arrows(x0=log10(10000), y0=36, x1=log10(900), y1=29, code=2, length=0.1)
    arrows(x0=log10(10000), y0=36, x1=log10(9000), y1=26, code=2, length=0.1)
    arrows(x0=log10(10000), y0=36, x1=log10(90000), y1=23, code=2, length=0.1)
    arrows(x0=log10(10000), y0=36, x1=log10(900000), y1=20, code=2, length=0.1)


    iPlotCalTab = 1:dim(calTabExt)[1]
    for (i in iPlotCalTab)
        with(calTabExt[i,], points(y=c(Ct.Mean+CI95.Mean, Ct.Mean-CI95.Mean),                               
                                x = rep(log10(copies.ml)+offsetMean,2), col='black', type='l', lwd=4))
        
    with( calTabExt, points(y=Ct.Mean, x=log10(copies.ml)+offsetMean, col='black', cex=1.5, pch=19))

    uselessLines = c(1:2,4:9,13:15,18,19)
    useMod = mod.Mean
    if (bUseBig) useMod = modBig
    ##szLegend =  c(szLegend, capture.output(summary(useMod))[-uselessLines])
    if (bUseBig) szModelName = "Normal Monte Carlo Linear" else szModelName = "Weighted Linear"
    szLegend = character(0)
    szLegend = c(szLegend, sprintf("Regression Model: %s", szModelName))
    ##printf("Model %s", format(useMod$call)) -- too wide for graph legend
    szLegend = c(szLegend, sprintf("Pts. x Draws  =  %d x %d = %d (NxM)", nPoints, nSim, nPoints*nSim))
    szLegend = c(szLegend, sprintf("Slope       m = %5.2f Ct/Log10(Direct Copies)", useMod$coefficients[2]))
    szLegend = c(szLegend, sprintf("Intercept   b = %5.2f Ct", useMod$coefficients[1]))
    szLegend = c(szLegend, sprintf("Resid Std Err = %5.2f on %d effective deg. of freedom", useMod$sigma, df))
                       
    szLegend =  c(szLegend, "Estimates:")
    szLegend =  c(szLegend, sprintf("Copies/mL  Ct Mean of ORF1ab/S/N"))
    for (i in 1:dim(calTabExt)[1])
        szLegend =  c(szLegend, with(calTabExt[i,], sprintf("    %d      %4.1f (CI95 %4.1f - %4.1f)",
                                                            copies.ml, Ct.Mean, Ct.Mean-CI95.Mean, Ct.Mean+CI95.Mean)))


    par(family="mono")
    legend(x=log10(1), y=21, cex=0.85, bty='o', bg='white', y.intersp=0.8, szLegend)
    par(family="")
    
    printf("%s\n", szLegend)
                  
    

    if (bRetModel)
        return(list(mod=useMod, calTab=calTab))
}


########################################################################################################################
## simVariance
## Uses Monte Carlo Technique to get variance of a vector each with known standard error

simVar<-function(
                 x, ##vector of items to find variance
                 SE.Xi = 0, ##vector of standard error of EACH x
                 nSim = 500000,
                 bDebug = F
                 
                 )
{
    nPoints = length(x)
    if (length(SE.Xi) != nPoints)
        {
        printf("Standard Error vector different length, assuming constant SE=%f\n", SE.Xi[1])
        SE.Xi = rep(SE.Xi[1], nPoints)
        }

    xsim = numeric(0)

    for (i in 1:nPoints)
        xsim = c(xsim, rnorm(n = nSim, mean = x[i], sd = SE.Xi[i]))
        
    nSim = length(xsim)
    xVar = var(xsim) * (nSim - 1)/nSim * nPoints/(nPoints - 1) 
    xSD = sqrt(xVar)
    xMean = mean(xsim)

    if (bDebug) printf("nPoints=%d nSim=%d xVar=%f  xSD=%f  xMean=%f\n", nPoints, nSim, xVar, xSD, xMean)
        
    return(xVar)

}
