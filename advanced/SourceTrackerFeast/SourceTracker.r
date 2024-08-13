# File: SourceTracker.r
# Author: Dan Knights
# Contact: danknights@gmail.com
# License: GPL
# Copyright: Copyright 2011, Dan Knights
# Version: 0.9.1 (Beta)

# Function "sourcetracker"
#
# Gets prior counts of taxa in source environments for use by SourceTracker.
#
# Params:
# train - a sample x feature observation count matrix
# envs - a factor or character vector indicating the sample environments
# rarefaction_depth - if not NULL, all samples with > rarefaction_depth sequences are rarified;
#            This decreases the influence of high-coverage source samples.
#
# Value: 
# Returns an object of class "sourcetracker". This is a list containing:
# sources - matrix of the total counts of each taxon in each environment
# train - a copy of the input training data
# envs - a copy of the input environment vector
"sourcetracker" <- function(train, envs, rarefaction_depth=1000){
    train <- as.matrix(train)

    # enforce integer data
    if(sum(as.integer(train) != as.numeric(train)) > 0){
        stop('Data must be integral. Consider using "ceiling(datatable)" or ceiling(1000*datatable) to convert floating-point data to integers.')
    }
    envs <- factor(envs)
    train.envs <- sort(unique(levels(envs)))
    
    # rarefy samples above maxdepth if requested
    if(!is.null(rarefaction_depth) && rarefaction_depth > 0) train <- rarefy(train, rarefaction_depth)
    if(!is.null(rarefaction_depth) && rarefaction_depth > 0) cat(sprintf('Rarefying training data at %d\n',rarefaction_depth))
    
    # get source environment counts
    # sources is nenvs X ntaxa
    sources <- t(sapply(split(data.frame(train), envs), colSums)) 
    
    # add an empty row for "Unknown"
    sources <- rbind(sources, rep(0,ncol(train)))
    rownames(sources) <- c(train.envs,"Unknown")
    colnames(sources) <- colnames(train)
    sources <- as.matrix(sources)
    
    ret <- list(sources=sources, train=train, envs=envs)
    class(ret) <- "sourcetracker"
    
    return(invisible(ret))
}


# Function "predict.sourcetracker"
#
# S3-level function to estimate source proportions using a sourcetracker object.
# Returns an object of class "sourcetracker.fit"
#
# Params:
# stobj - output from function "sourcetracker"
# test - test data, a matrix of sample x feature counts
#   if test is NULL, performs leave-one-out predictions of the training
#   samples
# burnin - number of "burn-in" passes for Gibbs sampling
# nrestarts - number of times to restart the Gibbs sampling process
# n.draws.per.restart - number of Gibbs draws to collect per restart
# delay - number of passes between draws (ignored if n.draws.per.restart is 1)
# alpha1 - prior counts of each species in the training environments,
#   Higher values decrease the trust in the training data, and make the 
#   source environment distributions over taxa smoother. By default, this is
#   set to 1e-3, which indicates reasonably high trust in all source environments, even
#   those with few training sequences. This is useful when only a small number 
#   of biological samples are available from a source environment.
#   A more conservative value would be 0.001 or 0.01.
# alpha2 - prior counts of each species in the Unknown environment,
#   Higher values make the Unknown environment smoother and less prone to
#   over-fitting a given training sample. Default is 1e-3.
# beta - prior counts of test sequences in each environment.
#   Higher values cause a smoother distribution over source environments.
#   Default is 1e-2.
# rarefaction_depth - if not NULL, all test samples with > maxdepth sequences are rarified
# verbosity - if > 0, print progress updates while running
# full.results - return full draws from gibbs sampling (all assignments of all taxa)
#
# Value: A list containing:
# draws - an array of dimension (ndraws X nenvironments X nsamples),
#     containing all draws from gibbs sampling
# proportions - the mean proportions over all Gibbs draws for each sample
# proportions_sd - standard deviation of the mean proportions for each sample
# train.envs - the names of the source environments
# samplenames - the names of the test samples
"predict.sourcetracker" <- function(stobj, test=NULL, 
            burnin=100, nrestarts=10, ndraws.per.restart=1, delay=10,
            alpha1=1e-3, alpha2=1e-1, beta=10, rarefaction_depth=1000,
            verbosity=1, full.results=TRUE){  #full.results=FALSE){

    if(!is.null(test)){
        # if test is numeric, cast as a row matrix
        if(class(test)[1] == "numeric" || class(test)[1] == "integer"){
            test <- matrix(test, nrow=1)
        } else {
            test <- as.matrix(test)
        }
        if(sum(as.integer(test) != as.numeric(test)) > 0){
            stop('Data must be integral. Consider using "ceiling(datatable)" or ceiling(1000*datatable) to convert floating-point data to integers.')
        }
        sources <- stobj$sources
        if(verbosity>=1) {
            cat(rep(' ',nrestarts * ndraws.per.restart+1),sep='')
            cat(rep(' ',27),sep='')
            cat(sprintf('%11s',substr(rownames(sources),1,10)),sep='\t'); cat('\n')
        }
        T <- ncol(sources) # number of taxa
        V <- nrow(sources) # number of source envs
        if(is.null(dim(test))) N <- 1
        else N <- nrow(test) # number of sink samples
        samplenames <- rownames(test)
        draws <- run.gibbs(sources, test, V, T, N,
                burnin=burnin, nrestarts=nrestarts, 
                ndraws.per.restart=ndraws.per.restart, delay=delay,
                alpha1=alpha1, alpha2=alpha2, beta=beta, maxdepth=rarefaction_depth,
                verbosity=verbosity, full.results=full.results)
        if(full.results) {
            full.draws <- draws$full.draws
            draws <- draws$draws
        }
    } else {  # leave-one-out    
        samplenames <- rownames(stobj$train)
        envs <- stobj$envs

        T <- ncol(stobj$train) # number of taxa
        V <- nrow(stobj$sources) # number of source envs
        N <- nrow(stobj$train) # number of sink samples
        ndraws <- nrestarts * ndraws.per.restart # total number of draws
        draws <- array(0,dim=c(ndraws, V, N))
        cat(sprintf('ndraws=%d, V=%d, T=%d, N=%d\n',ndraws, V, T, N))
        if(full.results){
            full.draws <- array(0,dim=c(ndraws, V, T, N))
        }
        for(i in (1:N)){
            stobj.i <- sourcetracker(stobj$train[-i,], envs[-i], rarefaction_depth=rarefaction_depth)
            sources <- stobj$sources
            V.i <- nrow(sources) # number of source envs (might be missing one if there's only one sample from this env)
            draws.i <- run.gibbs(sources, stobj$train[i,], V.i, T, 1,
                burnin=burnin, nrestarts=nrestarts, 
                ndraws.per.restart=ndraws.per.restart, delay=delay,
                alpha1=alpha1, alpha2=alpha2, beta=beta, maxdepth=rarefaction_depth,
                verbosity=verbosity, printing.index=i, printing.total=N, full.results=full.results)
            if(full.results){
                full.draws.i <- draws.i$full.draws
                draws.i <- draws.i$draws
            }
            # if(verbosity >= 1) cat(sprintf('%3d of %d: ',i,N))
            # handle case where there are no other samples from this env
            if(sum(envs[-i] == envs[i])==0){
                draws[,-which(rownames(stobj$sources)==envs[i]),i] <- drop(draws.i)
                if(full.results){
                    full.draws[,-which(rownames(stobj$sources)==envs[i]),,i] <- drop(full.draws.i)
                }
            } else {
                draws[,,i] <- drop(draws.i)
                if(full.results){
                    full.draws[,,,i] <- drop(full.draws.i)
                }
            }
        }
    }

    proportions <- matrix(nrow=N, ncol=V)
    proportions_sd <- matrix(nrow=N, ncol=V)
    for(i in 1:N){
        proportions[i,] <- apply(matrix(draws[,,i], ncol=V),2,mean)
        proportions_sd[i,] <- apply(matrix(draws[,,i], ncol=V),2,sd)
    }
    rownames(proportions) <- samplenames
    colnames(proportions) <- rownames(stobj$sources)
    rownames(proportions_sd) <- samplenames
    colnames(proportions_sd) <- rownames(stobj$sources)
    
    res <- list(draws=draws, proportions=proportions,
                proportions_sd=proportions_sd,
                train.envs=rownames(sources), samplenames=samplenames)
    dimnames(full.draws) <- list(
        sprintf('draw%05d',1:dim(full.draws)[1]),
        rownames(sources),
        colnames(sources),
        samplenames
    )
    if(full.results) res$full.results <- full.draws
    class(res) <- "sourcetracker.fit"
    return(invisible(res))
}


# Function "plot.sourcetracker.fit"
#
# S3-level function to plot the SourceTracker output.
#
# Params:
# stresult - output from function "predict.sourcetracker"
# labels - Labels for samples; if NULL, uses rownames of data table
# type - One of 'pie', 'bar', or 'dist' (distribution plots). Default is 'pie'.
# gridsize - number of samples to plot per row; if NULL, will be estimated
# ... - Additional graphical parameters
"plot.sourcetracker.fit" <- function(stresult, labels=NULL, 
        type=c('pie','bar','dist')[1], gridsize=NULL, env.colors=NULL, 
        titlesize=NULL, indices=NULL, include.legend=FALSE, ...){
    if(is.null(env.colors)){
        env.colors <- std.env.colors
        # always set 'Unknown' to grey
        env.colors[stresult$train.envs=='Unknown'] <- std.env.colors[length(std.env.colors)]
    }
    
    if(is.null(indices)) indices <- 1:dim(stresult$draws)[3]
    N <- length(indices)
    V <- dim(stresult$draws)[2]

    if(include.legend) N <- N + 1

    if(!is.null(gridsize) && gridsize**2 < N)
        stop(sprintf('Please choose a gridsize of at least %d.',ceiling(sqrt(N))))
        
    if(is.null(labels)) labels <- stresult[['samplenames']]
    if(is.null(gridsize)) gridsize <- ceiling(sqrt(N))
    if(is.null(titlesize)){
        if(gridsize > 1){
            titlesize <- .3/log10(gridsize)
        } else {
            titlesize=1
        }
    } 

    ngridrows <- ceiling(N / gridsize)
    par(mfrow=c(ngridrows,gridsize))
    par(oma=c(1,1,1,1), mar=c(0,0,titlesize,0))

    # legend will occupy one full plot in the upper left
    if(include.legend){
        plot(0,0,xlim=c(0,1), ylim=c(0,1), type='n', axes=FALSE)
        leg.cex <- 0.1
        leg <- legend('topleft',stresult$train.envs, fill=env.colors, bg='white', cex=leg.cex, plot=FALSE)
        maxdim <- max(leg$rect$w, leg$rect$h)

        # resize legend to be just big enough to fill the plot (80%), or to have maximum text size 2
        while(maxdim <.8 && leg.cex <= 2){
            leg.cex <- leg.cex + 0.01
            leg <- legend('topleft',stresult$train.envs, fill=env.colors, bg='white', cex=leg.cex, plot=FALSE)
            maxdim <- max(leg$rect$w, leg$rect$h)
        }
        leg <- legend('topleft',stresult$train.envs, fill=env.colors, bg='white', cex=leg.cex, border=NA)
    }

    if(type=='pie') plot.sourcetracker.pie(stresult, labels, gridsize, env.colors, titlesize, indices=indices, ...)
    if(type=='bar') plot.sourcetracker.bar(stresult, labels, gridsize, env.colors, titlesize, indices=indices, ...)
    if(type=='dist') plot.sourcetracker.dist(stresult, labels, gridsize, env.colors, titlesize, indices=indices, ...)
}


######### Internal functions below ####################


# Internal SourceTracker function to run Gibbs sampling
# total.n is used to supply the total number of samples in leave-one-out
# predictions for printing status updates
# full.results returns the source env. x taxon counts for every draw
"run.gibbs" <- function(sources, test, V, T, N,
        burnin=100, nrestarts=10, ndraws.per.restart=10, delay=10,
        alpha1=1e-3, alpha2=.1, beta=10, maxdepth=NULL,
        verbosity=1, printing.index=NULL, printing.total=NULL,
        full.results=TRUE){
        #full.results=FALSE){

    if(is.null(printing.total)) printing.total <- N
    
    train.envs <- rownames(sources)
    ndraws <- nrestarts * ndraws.per.restart # total number of draws
    npasses <- burnin + (ndraws.per.restart-1) * delay + 1 # passes per restart

    # draws will hold all draws (ndraws x V x N)
    draws <- array(dim=c(ndraws, V, N))
    if(full.results){
        full.draws <- array(dim=c(ndraws, V, T, N))
    }
    
    # rarefy samples above maxdepth if requested
    if(!is.null(maxdepth))    test <- rarefy(test, maxdepth)

    # sink samples must have integer counts
    if(is.null(dim(test))) test <- matrix(test,ncol=T)
    test <- round(test)
    
    # store original prior counts for "Unknown"
    unknown.prior <- sources[V,]
    # sources[-V,] <- sweep(sources[-V,],1,alpha1 * rowSums(sources[-V,]),'+')  # add relative alpha prior counts
    sources[-V,] <- sources[-V,] + alpha1 # add absolute alpha prior counts
    
    # for each sink sample
    for(i in 1:N){
        sink <- test[i,]
        D <- sum(sink) # sink sample depth
        
        # precalculate denominator for Pr(env in sample)
        p_v_denominator = max(1,(D-1) + V*beta)
        
        # get taxon index for each sequence
        tax.cumsum <- cumsum(sink)
        tax.ix <- sapply(1:D,function(x) min(which(x<=tax.cumsum)))
        
        drawcount <- 1 # keeps running count of draws for this sample
        # for each restart
        for(j in 1:nrestarts){
            if(verbosity>=1) cat('.')
            options(warn=-1)
            z <- sample(V,D,replace=TRUE) # random env assignments
            options(warn=0)

            sources[V,] <- unknown.prior # prior counts of taxa in Unknown
            sources[V,] <- sources[V,] + alpha2 * D # add relative alpha prior counts
            
            # tally counts in envs
            # count all assignments to the "other" environment
            # other environments don't get incremented because they are fixed from training data
            envcounts <- rep(beta, V)
            for(ix in 1:D){
                if(z[ix] == V)    sources[V,tax.ix[ix]] <- sources[V,tax.ix[ix]] + 1
                envcounts[z[ix]] <- envcounts[z[ix]] + 1
            }

            for(rep in 1:npasses){
                rand.ix <- sample(D) # random order for traversing sequence
                # temporary: not random
                # rand.ix <- 1:D
                
                cnt <- 0
                for(ix in rand.ix){
                    taxon <- tax.ix[ix]
                     # remove this sequence from all counts
                    envcounts[z[ix]] <- envcounts[z[ix]] - 1
                    if(z[ix] == V)    sources[V,taxon] <- sources[V,taxon] - 1

                    # get relative PDF over env assignments
                    p_tv <- sources[,taxon] / rowSums(sources) # Pr(taxon | env)
                    p_v <- envcounts/p_v_denominator# Pr(env in sample)

                    # re-sample this sequence's env assignment
                    z[ix] <- sample(1:V, prob=p_tv * p_v, size=1)

                    # replace this sequence in all counts
                    envcounts[z[ix]] <- envcounts[z[ix]] + 1

                    # if this sequence is assigned to "other", increase count
                    if(z[ix] == V)    sources[V,taxon] <- sources[V,taxon] + 1
                }
                
                # take sample
                if(rep > burnin && (((rep-burnin) %% delay)==1 || delay<=1)){
                        
                    # save current mixing proportions
                    draws[drawcount,,i] <- round((envcounts - beta) / D,7)
                    draws[drawcount,,i] <- draws[drawcount,,i] / sum(draws[drawcount,,i])
                    # save full taxon-source assignments if requested
                    if(full.results){
                        # for each environment, save taxon counts
                        for(j in 1:V){
                            full.draws[drawcount,j,,i] <- sapply(1:T,function(x) sum(tax.ix[z==j]==x))
                        }
                    }
                    drawcount <- drawcount + 1
                }
            }
        }

        if(verbosity>=1){
            if(is.null(printing.index)){
                cat(sprintf('%4d of %4d, depth=%5d: ', i, printing.total, D))
            } else {
                cat(sprintf('%4d of %4d, depth=%5d: ', printing.index, printing.total, D))
            }
            props <- colMeans(matrix(draws[,,i],ncol=V))
            prop_devs <- apply(matrix(draws[,,i], ncol=V), 2, sd)
            cat(' ')
            cat(sprintf('%.2f (%.2f)', props, prop_devs),sep='\t')
            cat('\n')
        }
    }
    if(full.results){
        return(list(draws=draws, full.draws=full.draws))
    } else {
        return(draws=draws)
    }
}



# tries all values of alpha1 and alpha2 for best RMSE
# if individual.samples, tries to predict mixtures of single samples
# instead of mixtures of the environment means
# ntrials is the number of simulated samples per fit.
# nrepeats is the number of times to repeat the entire experiment at each alpha value
# verbosity > 1 means inner loop will print
"tune.st" <- function(otus, envs, individual.samples=TRUE, ntrials=25,
            rarefaction_depth=1000, alpha1=10**(-3), 
            alpha2=10**(-3:0), beta=10, verbosity=0, ...){
    results <- list()
    alphas <- expand.grid(rev(alpha1), alpha2)
    colnames(alphas) <- c('alpha1','alpha2')
    rmse <- numeric(nrow(alphas))
    rmse.sem <- numeric(nrow(alphas))
    for(i in 1:nrow(alphas)){
        cat(sprintf('Loop %d of %d, alpha1=%f, alpha2=%f ',i,nrow(alphas),alphas[i,1], alphas[i,2]))
        if(verbosity > 2) cat('\n')
        results[[i]] <- eval.fit(otus, envs, individual.samples=individual.samples,
                                ntrials=ntrials, rarefaction_depth=rarefaction_depth,
                                alpha1=alphas[i,1], alpha2=alphas[i,2], beta=beta, verbosity=verbosity-1, ...)
        rmse[i] <- results[[i]]$rmse
        rmse.sem[i] <- results[[i]]$rmse.sem
        if(verbosity > 0) cat(sprintf('RMSE = %.3f +/- %.3f\n',rmse[i], rmse.sem[i]))
    }
    # choose alpha as most conservative value of alpha2 (smallest)
    # then most conservative value of alpha1 (largest)
    # that gives pseudo-r2 within 1 sem of the max.
    # best.ix <- min(which(rmse <= min(rmse + rmse.sem)))
    # Alternative: simply choose the lowest rmse
    best.ix <- which.min(rmse)
    best.rmse <- rmse[best.ix]
    best.alpha1 <- alphas[best.ix,1]
    best.alpha2 <- alphas[best.ix,2]
    return(list(alphas=alphas, rmse=rmse, rmse.sem=rmse.sem, best.rmse=best.rmse, 
                best.alpha1=best.alpha1, best.alpha2=best.alpha2, 
                results=results))
}

# train SourceTracker object on training data
# returns RMSE
# ... are additional params to pass to sourcetracker predict object
"eval.fit" <- function(otus, envs, individual.samples=TRUE,
            ntrials=25, rarefaction_depth=1000, verbosity=1, ...){
    train.envs <- sort(unique(envs))
    V <- length(train.envs)
    env.sizes <- table(envs)
    
    # make sure each pair of envs gets picked
    # build up all pairs of samples, each column is a pair
    # each source env gets to be first and second sample once
    # pairs <- expand.grid(1:V,1:V)
    # pairs <- pairs[pairs[,1]!=pairs[,2],]
    # make nreps pairs randomly
    pairs <- NULL
    for(i in 1:ntrials){
        pairs <- rbind(pairs, sample(V,size=2))
    }
    
    mixtures <- runif(ntrials)
    y <- matrix(0,nrow=ntrials, ncol=V+1)
    yhat <- matrix(0,nrow=ntrials, ncol=V+1)
    yhat.sd <- matrix(0,nrow=ntrials, ncol=V+1)
    colnames(y) <- c(as.character(train.envs),'Unknown')
    colnames(yhat) <- c(as.character(train.envs),'Unknown')
    newsamples <- NULL
    allenvs <- NULL
    for(i in 1:ntrials){
        env1 <- pairs[i,1]
        env2 <- pairs[i,2]
        allenvs <- rbind(allenvs, c(env1, env2))
        if(verbosity > 1){
           cat(sprintf('%d of %d: %.2f*%s + %.2f*%s: \n',i,ntrials,mixtures[i], train.envs[env1],1-mixtures[i], train.envs[env2]))
        } else if(verbosity > 0){
            cat('.')
        }

        # all indices of each environment
        env1.ix.all <- which(envs == train.envs[env1])
        env2.ix.all <- which(envs == train.envs[env2])
        
        if(individual.samples){
            # get one sample from each env
            # cast as list so that sample doesn't misinterpret a length-1 vector
            env1.ix <- sample(as.list(env1.ix.all),size=1)[[1]]
            env2.ix <- sample(as.list(env2.ix.all),size=1)[[1]]
            
            # train sourcetracker, hold out entire second env. and first env. sample
            # note: don't hold out first sample if that env has only one sample
            if(length(env1.ix.all) == 1){
                st <- sourcetracker(otus[-env2.ix.all,], envs[-env2.ix.all])
            } else {
                st <- sourcetracker(otus[-c(env1.ix,env2.ix.all),], envs[-c(env1.ix,env2.ix.all)])
            }
            
            # make fake sample, weighted mixture of two source samples
            s1 <- otus[env1.ix,]
            s2 <- otus[env2.ix,]
            
        } else {
            # train sourcetracker, hold out entire second env.
            st <- sourcetracker(otus[-env2.ix.all,], envs[-env2.ix.all])
            
            # make fake sample as mixture of _environment_ means
            s1 <- colSums(rarefy(otus[env1.ix.all,], maxdepth=rarefaction_depth))
            s2 <- colSums(rarefy(otus[env2.ix.all,], maxdepth=rarefaction_depth))
        }
        
        newsample <- mixtures[i] * s1/sum(s1) + (1-mixtures[i]) * s2/sum(s2)
        newsample <- round(100000 * newsample)
        newsample <- matrix(newsample, nrow=1)
        newsample <- rarefy(newsample,maxdepth=ceiling(sum(s1+s2)/2))
        newsamples <- rbind(newsamples, newsample)
        y[i,env1] <- mixtures[i]
        y[i,V+1] <- 1-mixtures[i]
        
        # test on fake sample
        results <- predict(st, newsample, rarefaction_depth=rarefaction_depth, verbosity=verbosity-1, ...)
        for(j in 1:ncol(results$proportions)){
            whichenv <- which(colnames(yhat) == colnames(results$proportions)[j])
            yhat[i,whichenv] <- results$proportions[,j]
            yhat.sd[i,whichenv] <- results$proportions_sd[,j]
        }
    }

    # calculate RMSE
    se <- as.numeric((y[,-V] - yhat[,-V])**2)
    mse <- mean(se)
    se.sem <- sd(se)/sqrt(length(se))
    rmse <- mse**.5
    rmse.sem <- se.sem**.5
    
    return(list(y=y,yhat=yhat,yhat.sd=yhat.sd,newsamples=newsamples, 
            env.pairs=allenvs, train.envs=train.envs, rmse=rmse, rmse.sem=rmse.sem))
}


# Internal SourceTracker function to perform rarefaction analysis
"rarefy" <- function(x,maxdepth){
    if(is.null(maxdepth)) return(x)
    
    if(!is.element(class(x)[1], c('matrix', 'data.frame','array')))
        x <- matrix(x,nrow=1)
    nr <- nrow(x)
    nc <- ncol(x)

    for(i in 1:nrow(x)){
        if(sum(x[i,]) > maxdepth){
            options(warn=-1)
            s <- sample(nc, size=maxdepth, prob=x[i,], replace=T)
            options(warn=0)
            x[i,] <- hist(s,breaks=seq(.5,nc+.5,1), plot=FALSE)$counts
        }
    }
    return(x)
}

# Internal SourceTracker function to plot pie charts
"plot.sourcetracker.pie" <- function(stresult, labels, 
        gridsize, env.colors, titlesize, indices, ...){
    V <- length(stresult$train.envs)
    for(i in indices){
        props <- apply(matrix(stresult$draws[,,i], ncol=V),2,mean)
        pie(props, labels=NA, col=env.colors, main=labels[i],cex.main=titlesize, ...)
    }
}

# Internal SourceTracker function to plot bar plots
"plot.sourcetracker.bar" <- function(stresult, labels    , 
        gridsize, env.colors, titlesize, indices, ...){
    V <- length(stresult$train.envs)
    # add extra space at top for title
    new.margins <- par('mar')
    new.margins[3] <- max(.5, new.margins[3] * 1.5)
    new.margins[2] <- .25 * titlesize
    par(mar=new.margins)
    for(i in indices){
        props <- apply(matrix(stresult$draws[,,i], ncol=V),2,mean)
        prop_devs <- apply(matrix(stresult$draws[,,i], ncol=V),2,sd)
        centers <- barplot(props, col=env.colors, main=labels[i],
                cex.main=titlesize, axes=FALSE, axisnames=FALSE,
                ylim=c(0,1.5), ...)
        sourcetracker.error.bars(centers, props, prop_devs)
        for(j in 1:4)  axis(j, at=c(-100,100),labels=FALSE)
    }
}

# Internal SourceTracker function to plot distribution plots
"plot.sourcetracker.dist" <- function(stresult, labels, 
        gridsize, env.colors, titlesize, indices, sortmethod=c('divergence', 'multilevel')[1], ...){
    # stop conditions
    if(dim(stresult$draws)[1] < 2)
        stop('Distribution plots require more than one draw.')
    V <- length(stresult$train.envs)
    for(i in indices){
        x <- stresult$draws[,,i]
        rownames(x) <- 1:nrow(x)
        if(sortmethod=='multilevel'){
            # sort by size of column
            sortby.ix <- sort(colMeans(x),index=T, dec=T)$ix
            x <- matrix(x, ncol=V)
            ix <- sortmatrix(x[,sortby.ix])
            x <- x[ix,]
        } else {
            d <- jsdmatrix(x)
            if(sum(d) > 0){
                # ensure that all mixtures were not identical (e.g. 100% unknown)
                ix <- cmdscale(jsdmatrix(x),k=1)
                ix <- sort(ix,index=T)$ix
                x <- x[ix,]
            }
        }
        centers <- barplot(t(x), beside=FALSE, col=env.colors,
                space=-1/ncol(x), border=NA, axes=FALSE, axisnames=FALSE,
                main=labels[i],cex.main=titlesize,
                ylim=c(-.05,1.05), ...)
        bounds <- c(0, min(centers)-.5, 1, max(centers)+.5)
        
        lines(c(bounds[2], bounds[4]), c(bounds[1], bounds[1]), lty=1, lwd=1)
        lines(c(bounds[2], bounds[4]), c(bounds[3], bounds[3]), lty=1, lwd=1)
        lines(c(bounds[2], bounds[2]), c(bounds[1], bounds[3]), lty=1, lwd=1)
        lines(c(bounds[4], bounds[4]), c(bounds[1], bounds[3]), lty=1, lwd=1)
    }
}

# Internal SourceTracker function to plot error bars on a bar chart
"sourcetracker.error.bars" <- function(x,centers,spread,...){
    width = .01
    
    xlim <- range(x)
    barw <- diff(xlim) * width
    
    upper <- centers + spread
    lower <- centers - spread
    keepix <- which(spread > 0)
    x <- x[keepix]
    upper <- upper[keepix]
    lower <- lower[keepix]
    
    segments(x, upper, x, lower, lwd=1.5, ...)
    segments(x - barw, upper, x + barw, upper, lwd=1.5, ...)
    segments(x - barw, lower, x + barw, lower, lwd=1.5, ...)
}

# plot type 1: predicted fit verses actual proportion of that site (when source is included)
# plot type 2: predicted fit versus actual proportion of every other site, when other source is hidden
"plot.eval" <- function(res, plot.all=FALSE,plot.type=c(1,2)[1], filename=NULL){
    
    V <- length(res$train.envs)
    if(plot.type==1){
        pch <- 16
        mycolors <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#A65628', '#F781BF', '#999999', '#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#D9D9D9', '#BC80BD', '#CCEBC5')
        mycolors <- sprintf('%sBB',mycolors)
        if(!is.null(filename)) pdf(filename,width=6,height=6)
        plot(0,0,xlim=0:1,ylim=0:1,type='n',xlab='True proportion', ylab='Estimated proportion')
        abline(0,1)
        for(i in 1:(V+1)){
            # sourcetracker.error.bars(res$y[,i], res$yhat[,i], res$yhat.sd[,i],col='#33333399')
            points(res$y[,i],res$yhat[,i],col=mycolors[i], pch=pch)
        }
        legend('topleft',legend=c(as.character(res$train.envs),'Unknown'),col=mycolors,pch=pch, cex=.75)
        if(!is.null(filename)) dev.off()
    } else {
        if(!is.null(filename)) pdf(filename,width=10,height=10)
        par(mfrow=c(V,V),mar=rep(.5,4))
        par(mfrow=c(V,V),mar=c(4.5,4.5,.5,.5))
        for(i in 1:V){
            for(j in 1:V){
                # keep only those points where j is the hidden source and i
                # is not present
                keepix <- which(res$env.pairs[,1] != i & res$env.pairs[,2] == j)
                # don't plot same-same plots
                if(i==j) {
                    plot(0,0,xlim=c(-1,1), ylim=c(-1,1),type='n',xaxt='none',yaxt='none',xlab='', ylab='')
                    text(0,0,res$train.envs[i], cex=2)
                } else if(length(keepix) < 3){
                    plot(0,0,type='n',xaxt='none', yaxt='none',xlab='', ylab='')
                    text(0,0,'Too few samples')
                } else if (i!=j){
                    a <- res$y[keepix,V+1] # actual hidden proportion
                    b <- res$yhat[keepix,i]; # recovered proportion for this other source
                    rmse <- mean(b**2)**.5
                    plot.color <- 'black'
                    # plot red if RMSE > .1
                    if(rmse > .1) plot.color='red'
                    plot(a,b,type='p',xaxt='none',yaxt='none',
                         xlim=0:1,ylim=0:1,col=plot.color,
                         xlab='', ylab='')
                    # draw axis labels close to the axes
                    cex = .65
                    # mtext(sprintf('Prop. %s',res$train.envs[j]), side=1, line=1, cex=cex)
                    # mtext(sprintf('Prop. %s when %s absent',res$train.envs[i],res$train.envs[j]), side=2, line=1, cex=cex)
                    mtext(sprintf('Prop. Env 1',res$train.envs[j]), side=1, line=1, cex=cex)
                    mtext(sprintf('Prop. Env 2 when 1 absent',res$train.envs[i],res$train.envs[j]), side=2, line=1, cex=cex)
                    cex=1
                    text(0.05,.90,sprintf('RMSE = %.3f',rmse),adj=c(0,0),cex=cex,col=plot.color)
                }
            }
        }
        if(!is.null(filename)) dev.off()
    }
}

"save.mapping.file" <- function(results, map, filename='map.txt',
        include.contamination.predictions=FALSE){
    # also append results to a new version of the mapping file in output directory
    # add columns filled with NA to the mapping file, add values where available
    for(j in 1:ncol(results$proportions)){
        colname <- sprintf('Proportion_%s',colnames(results$proportions)[j])
        colname_sd <- sprintf('Proportion_SD_%s',colnames(results$proportions)[j])
        map[,colname] <- rep(NA,nrow(map))
        map[rownames(results$proportions), colname] <- results$proportions[,j]
        map[rownames(results$proportions), colname_sd] <- results$proportions_sd[,j]
    }
    
    # add this environment and its s.d. to mapping file
    thisenv <- rep(NA,nrow(results$proportions))
    thisenv_sd <- rep(NA,nrow(results$proportions))
    for(i in 1:nrow(results$proportions)){
        sampleid <- rownames(results$proportions)[i]
        whichenv <- which(results$train.envs==map[sampleid,'Env'])
        if(length(whichenv) > 0){
            thisenv[i] <-  results$proportions[i,whichenv]
            thisenv_sd[i] <- results$proportions_sd[i,whichenv]
        }
    }
    map[,'Proportion_This_Env'] <- rep(NA, nrow(map))
    map[,'Proportion_SD_This_Env'] <- rep(NA, nrow(map))
    map[rownames(results$proportions), 'Proportion_This_Env'] <- thisenv
    map[rownames(results$proportions), 'Proportion_SD_This_Env'] <- thisenv_sd
    
    if(include.contamination.predictions){
        for(threshold in seq(0.05, 0.95, .05)){
            header <- sprintf('contaminated_at_%.02f',threshold)
            map[,header] <- rep(NA, nrow(map))
            map[rownames(results$proportions), header] <- as.character(thisenv < threshold)
        }
    }
        
    # write new mapping file
    sink(filename)
    cat('#SampleID\t')
    write.table(map,sep='\t',quote=F)
    sink(NULL)
}

# sorts a matrix by first column, breaks ties by the 2nd, 3rd, etc. columns
# returns row indices
"sortmatrix" <- function(x){
    # sort by last column, then 2nd-to-last, etc.
    ix <- 1:nrow(x)
    for(j in ncol(x):1){
        ixj <- sort(x[ix,j], index=T)$ix
        ix <- ix[ixj]
    }
    return(ix)
}

"jsdmatrix" <- function(x){
    d <- matrix(0,nrow=nrow(x),ncol=nrow(x))
    for(i in 1:(nrow(x)-1)){
        for(j in (i+1):nrow(x)){
            d[i,j] <- jsd(x[i,], x[j,])
            d[j,i] <- d[i,j]
        }
    }
    return(d)
}

"jsd" <- function(p,q){
    m <- (p + q)/2
    return((kld(p,m) + kld(q,m))/2)
}

"kld" <- function(p,q){
    nonzero <- p>0 & q>0
    return(sum(p[nonzero] * log2(p[nonzero]/q[nonzero])))    
}

# global definition of standard env colors
std.env.colors <- c(
'#885588',
'#CC6666',
'#47697E',
'#5B7444',
'#79BEDB',
'#663333',
'#3F52A2',
'#128244',
'#e93E4A',
'#1DBDBC',
'#A43995',
'#FFCC33',
'#B1BDCD',
'#A3C586',
'#6B78B4',
'#266A2E',
'#FCF1D1',
'#660F57',
'#272B20',
'#003366',
'#656565'
)

