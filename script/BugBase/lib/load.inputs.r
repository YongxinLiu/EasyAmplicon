# Loads the user-inputs and makes sure they are valid
# Required: otu table, mapping file,
# Optional: map column, groups
# Returns: loaded otu_table, mapping file

"load.inputs" <- function(otu_table, map=NULL, map_column=NULL, groups=NULL){
  
  #check for an otu table
  if(is.null(otu_table)){
    stop("\nError: No otu table specified.\n")
  }
  
  #determine file type of otu table (json .biom or .txt)
  otu_ext <- tail(strsplit(otu_table, ".", fixed=T)[[1]], n=1)
  
  #otu_table is otus (rows) x samples (columns)
  if(!file.exists(otu_table)){
    cat("\nError: OTU table specified does not exist. Check file path.\n")
    stop("\nError listed above\n")
  } else {
      if(otu_ext == "txt"){
      line1 <-readLines(otu_table,n=1)
      if(line1=="# Constructed from biom file") {
        otu_table <- as.matrix(read.table(otu_table, sep='\t', head=T, row=1, 
                                      check=F, comment='', skip=1))
      } else {
        otu_table <-as.matrix(read.table(otu_table, sep='\t', head=T, row=1, 
                                      check=F, comment=''))
      }
    } else {
      if(otu_ext == "biom"){
        otu_table <- as.matrix(biom_data(read_biom(otu_table)))
      } else {
        stop("\nError: otu table must be either .txt or .biom (json)\n")
      }
    }
  }

  if(is.null(map)){
    cat("\nEither \'predict only\' or \'plot all\' specified. All samples will be predicted.\n")

    #Drop any samples with no OTUs
    otu_table <- otu_table[,colSums(otu_table) > 1,drop=F]

    new_otu <- otu_table
    new_map <- NULL

  } else {
    #map is samples x metadata
    if(!file.exists(map)){
      cat("\nError: Map specified does not exist. Check file path.\n")
      stop("\nError listed above\n")
    } else {
      map <- read.table(map, sep='\t', head=T, row=1 ,check=F, comment='')
    }
    
    if(is.null(map_column)){
      stop("\nError: No map column specified. To run BugBase without a mapping file use '-a'\n")
    }
    
    if(! map_column %in% colnames(map)){
      cat("\nError: Map column specified does not exist\n")
      cat("\nThe map columns available are:\n")
      cat(colnames(map),"\n\n")
      stop("\nError listed above\n")
    }
    
    #print the number of samples matching between the otu table and map
    # cat("\nThe number of samples matching between the otu table and mapping file:")
    # cat("\n",length(intersect(rownames(map),colnames(otu_table)))[1],"\n")
    
    #define treatment groups
    if(is.null(groups)){
      if(isTRUE(opts$continuous)){
        groups <- as.numeric(map[,map_column])
      } else {
        groups <- unique(map[,map_column])
        groups <- lapply(groups, as.character)
        if(length(groups) <= 1){
          stop("\nError: a minimum of two groups must be tested\n")
        }
      }
    } else {
      #use user-defined groups
      groups <- strsplit(groups, ",")[[1]]
      #remove any duplicates
      groups <- unique(groups)
      if(length(groups) <= 1){
        stop("\nError: a minimum of two groups must be tested\n")
      }
      #define the groups that exist in that mapping column
      groups_avail <- unique(map[,map_column])
      #check is groups specified exist
      for(g in groups){
        if(! g %in% groups_avail){
          cat("\nError:\n")
          cat(g," ")
          cat("is not a group listed in the mapping file.\n")
          cat("These are the groups available:\n")
          cat(groups_avail,"\n\n")
          stop("\nError listed above\n")
        }
      }
    #factor groups so they appear in user-listed order
    map[,map_column] <- factor(map[,map_column],groups)  
    }

    #Drop any samples with no OTUs
    otu_table <- otu_table[,colSums(otu_table) > 1]

    #Change those less than 1/1 millionth of read depth to 0
    otu_table[otu_table < sum(colSums(otu_table))/1000000] <- 0

    #Change singletons to 0 (needed for low depth OTU tables)
    otu_table[otu_table < 2] <- 0

    ##Filter the OTU table to keep OTUs in at least 2 sample
    otu_table <- otu_table[rowSums(otu_table > 0) > 2,]
    #otu_table <- otu_table[rowSums(otu_table > 0) > (0.05*ncol(otu_table)),]

    #get indices of which rows to keep
    ix.keep <- map[,map_column] %in% groups
    #keep only subset of samples belonging to requested groups
    new_map <- droplevels(as.data.frame(map[ix.keep,,drop=F]))
    
    #keep only samples that intersect between the map and otu table
    intersect_btwn <- intersect(rownames(map),colnames(otu_table))
    new_map <- map[intersect_btwn,,drop=F]
    new_otu <- droplevels(as.data.frame(otu_table[,intersect_btwn]))

    #print the groups to be tested
    # cat("\nThe number of samples in each group are:")
    # print(table(new_map[,map_column]))
    # map <- new_map
    # otu_table <- new_otu
  }
  rownames(new_otu) <- gsub("\ ", "_", rownames(new_otu))

  return(list(
    otu_table=new_otu,
    map=new_map,
    groups=groups,
    map_column=map_column))
}

