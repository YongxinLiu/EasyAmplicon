# Loads the trait table and 16s-normalized otu table.
# Calculates the variance of trait abundance for thresholds 0 to 1.
# Sets the default variance to the max variance for each trait (if no user
#   threshold was set).
# Does single cell predictions for all otus and traits in the input tables
# Required: trait table, 16S normalized otu table
# Optional: user-defined threshold threshold_set, use coefficient of variance 
#   (default is variance)
# Returns: prediction array, final predictions (both are written to tables)

"single.cell.predictions" <- function(trait_table_fp, 
                                      otu_table_fp,
                                      test_trait, 
                                      threshold_set,
                                      use_cov,
                                      clr_trans,
                                      thresholds=c(seq(0, 0.01, 0.001), 
                                                   seq(0.02, 0.1, 0.01), 
                                                   seq(0.2, 1.0, 0.1))){
  
  #load the trait table and 16S normalized otu table
  #trait_table is otu x trait
  trait_table <- as.matrix(read.table(trait_table_fp, 
                                      sep='\t', 
                                      head=TRUE, 
                                      row=1, 
                                      check=F,
                                      quote=""))
  
  if(!is.null(threshold_set)){
    if(! 0 < threshold_set && threshold_set <= 1){
      stop("Error:Threshold must be between 0 and 1")
    }
  }
  #define traits to test (either all or defined)
  if(is.null(test_trait)){
    traits <- colnames(trait_table)
  } else {
    test_trait <- strsplit(test_trait, ",")[[1]]
    trait_missing <- which(! test_trait %in% colnames(trait_table))
    if(length(trait_missing) > 0){
      print("Error: The following specified traits do not exist:")
      print(test_trait[trait_missing])
      print("These are the available traits:")
      print(colnames(trait_table))
      stop("Error is specified above")
    } else {
      traits <- test_trait
    }
  }
  ###############################Process colnames, give map for kegg traits
  
  #subset trait table if specific trait is wanted
  trait_keep <- intersect(colnames(trait_table),traits)
  trait_table <- trait_table[,trait_keep, drop=F]
  
  #otu_table is otu x samples
  #otu_table <- as.matrix(read.table(otu_table_fp, 
  #sep='\t', 
  #head=T, 
  #row=1, 
  #check=F, 
  #comment='', 
  #skip=1))
  otu_table <- otu_table_fp

  if(is.null(clr_trans)){
    #convert to samples as rows
    otu_table <- t(otu_table)
  
    #relative abundance conversion of otu table
    otu_table <- sweep(otu_table, 1, rowSums(otu_table), FUN='/')
  } else {
    #Convert any 0 to 0.65 to allow for CLR transform
    #Ref: Palarea-Albaladejo J, et al. 2014. JOURNAL OF CHEMOMETRICS. A bootstrap estimation scheme for chemical compositional data with nondetects. 28;7:585–599.
    #otu_table[otu_table == 0] <- 0.65
    
    #Centered log-ratio transform for compositions
    #Ref: Gloor GB, et al. 2016. ANNALS OF EPIDEMIOLOGY. It's all relative: analyzing microbiome data as compositions. 26;5:322-329.
    #Multiplicative:
    #Mart´ın-Fern´andez JA, et al. 2003. MATHEMATICAL GEOLOGY. Dealing With Zeros and Missing Values in Compositional Data Sets Using Nonparametric Imputation. 35;3:253-278. 
    
    #convert to samples as rows
    otu_table <- t(otu_table)
    #otu_table <- clr(otu_table)   # Centered log-ratio transform for compositions 
    eps = 0.5
    otu_table = otu_table*(1 - rowSums(otu_table==0)*eps/rowSums(otu_table))
    otu_table[otu_table==0]=eps
    otu_table = sweep(otu_table,1,rowSums(otu_table),'/')
    ls = log(otu_table)
    otu_table = ls - rowMeans(ls)
    otu_table = otu_table[!is.nan(rowSums(otu_table)),]
  }

  #define which otus are in both tables
  otus_keep <- intersect(colnames(otu_table),rownames(trait_table))

  #make sure otus_keep > 0
  if(length(otus_keep)==0){
    stop("Error: no OTU overlap between OTU table and trait table.")
  }
  
  print(paste(length(otus_keep), "OTUs from the input table matched the", nrow(trait_table), "available database OTUs"))

  #subset pathway table and otu table to include only otus within the 
  # otu table
  otu_table <- otu_table[,otus_keep,drop=F]
  trait_table <- trait_table[otus_keep,,drop=F]
  
  #if threshold is user-defined, set it here
  if(!is.null(threshold_set)){
    thresholds <- threshold_set
  }
  
  nsamples <- nrow(otu_table)
  notus <- nrow(trait_table)
  ntraits <- ncol(trait_table)
  nthresholds <- length(thresholds)
  
  #prediction is samples x traits x thresholds
  prediction <- array(0, dim=c(nsamples, ntraits, nthresholds))
 
  for(i in 1:length(thresholds)){
    th <- thresholds[i]
    trait_table_gt <- trait_table > th
    trait_table_eq <- apply(trait_table, c(1,2), function(xx) 
      isTRUE(all.equal(xx,th))) # Floating pt issue
    trait_table_th <- trait_table_gt | trait_table_eq
    prediction[,,i] <- otu_table %*% trait_table_th
  }
  #name the rows (samples), columns (traits) and layers (thresholds)
  # of the prediction array
  dimnames(prediction)[[1]] <- rownames(otu_table)
  dimnames(prediction)[[2]] <- colnames(trait_table)
  dimnames(prediction)[[3]] <- thresholds

  #make variance table showing the variance in trait abundance
  #variances is trait x threshold

  if(is.null(use_cov)){
    variances <- apply(prediction, c(2,3), var)
  } else {
    variances <- apply(prediction, c(2,3), function(xx) sd(xx)/mean(xx))
  }
  rownames(variances) <- colnames(trait_table)
  colnames(variances) <- thresholds
  
  #populate final prediction table
  #final_prediction is samples x traits
  final_prediction <- matrix(0, nsamples, ntraits)
  
  if(is.null(threshold_set)){
    #use the max variance as the threshold for each trait
    which.threshold <- apply(variances, 1, which.max)
    which.threshold[which.threshold == 1] <- 2
    #create a matrix of traits and thresholds used
    trait_thresholds <- matrix(0, ntraits, 1)
    rownames(trait_thresholds) <- rownames(variances)
    trait_thresholds[,1] <- as.numeric(thresholds[which.threshold])
    colnames(trait_thresholds) <- "Threshold"
  } else {
    trait_thresholds <- matrix(0, ntraits, 1)
    rownames(trait_thresholds) <- rownames(variances)
    trait_thresholds[,1] <- as.numeric(threshold_set)
    colnames(trait_thresholds) <- "Threshold"
  }
  
  #create a matrix of otus by traits (boolean) using the thresholds used for
  # each trait
  otus_contributing <- matrix(0, notus, ntraits)
  rownames(otus_contributing) <- rownames(trait_table)
  colnames(otus_contributing) <- colnames(trait_table)
  otus_gt <- otus_contributing
  otus_eq <- otus_contributing
  
  #otus_contributing will be otus x traits
  for(i in 1:ntraits){
    for(x in 1:notus){
      otus_gt[x,i] <- trait_table[x,i] > trait_thresholds[i,1]
      otus_eq[x,i] <- isTRUE(all.equal(trait_table[x,i], 
                                       trait_thresholds[i,1]))
    }
  }
  otus_contributing <- otus_gt | otus_eq
  
  #fill in the final prediction table
  #this is samples x traits
  for(i in 1:ntraits){
    if(is.null(threshold_set)){
      final_prediction[,i] <- prediction[,i,which.threshold[i]]
    } else {
      final_prediction[,i] <- prediction[,i,1]
    }
  }
  
  rownames(final_prediction) <- rownames(otu_table)
  colnames(final_prediction) <- colnames(trait_table)
  
  #write the prediction file and variance file to text outputs
  var_name <- paste(output, "thresholds/variances.txt", sep='/')
  write.table(variances, var_name, sep='\t', quote=F, 
              col.names=NA)
  prediction_name <- paste(output, "predicted_phenotypes/predictions.txt", 
                           sep='/')
  write.table(final_prediction, prediction_name, sep='\t', quote=F, 
              col.names=NA)
  
  thresholds_name <- paste(output, "thresholds/thresholds_used.txt", sep='/')
  write.table(trait_thresholds, thresholds_name, sep='\t', quote=F, 
              col.names=NA)
  
  return(list(
    predictions=prediction,
    final_predictions=final_prediction,
    thresholds_used=trait_thresholds,
    otus_contributing=otus_contributing,
    otu_table_subset=otu_table
  ))
}

