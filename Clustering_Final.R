##############################################################################################################
# Code developed by Eric Roden for research done by Erin Jewett in the Hegeman Lab at the University of
# Minnesota - Twin Cities. Please properly site our work if you use this script in your own analysis. 

# Also feel free to contact the developer at roden026@umn.edu with suggestions for code improvement or if any
# bugs are found.

# Please note this code has been written for readability for non-technical users, improved efficiency may 
# be accomplished with edits to the code if speed and storage are or high importance to you. 
##############################################################################################################

##############################################################################################################
# Loading Library - May need to install missing package first but will only need to do that once.
##############################################################################################################

# Install if needed
# install.packages(plyr)
# install.packages("RColorBrewer")

# Load library
library(plyr)
library(RColorBrewer)

##############################################################################################################
# User Input - Check that all the values are correct and then run this secton.
# Note: filenames should be unique within the folder. Non-unique names will lead
# to files being overwritten or the script throwing errors. 
##############################################################################################################

# Location for plots and tables to be output by this script
# remember if using windows you need to change all \ to / in the file path
outputPathCluster <- "/home/erin/Desktop/Thesis Data/No_Suc_Ghostweed_MZML/No_Suc/Normal/Res_No_4Rep/Clusters/"

## Note: All filenames are output filenames
# Set filename
unlabeledFilename <- "unlabeled_table2.csv"

# Set filename
labeledFilename <- "labeled_table2.csv"

# Set filename
allDataFilename <- "allData_table2.csv"

# k finder plot filename
kfinderFilename <- "k_finder2.pdf"

# clustered line plot filename
clusterFilename <- "clusters.pdf"

# (optional) Dendrogram filename for saving graphical result of hierarchical clustering
dendroFilename <- "dendrogram2.pdf"
##

# point distance calculation. See dist() function documentation for other options
pointDistance <- 'euclidean'

# cluster distance calculation methods for kmeans() and hclust() functions.
# see their documentation for other options. Defaults are shown here. 
kmeansAlgorithm <- "Hartigan-Wong"
hclustMethod <- "ward.D"

# Pick which data set to analyze. 'labeled' or 'unlabeled' or 'all'
typeOfAnalysis <- 'unlabeled'

# Pick 'h' or 'k' to use either hierarchical clustering or k-means clustering
typeOfClustering <- 'k'

# cluster by regression values ('slope') or abundance ratios ('ratio'). 
# used in aggregateData function and is 'ratio' by default. 
clusterBy <- 'slope'

# Set to TRUE to show aggregate lines on the plot for each cluster or set to FALSE to ignore. 
# If TRUE, aggregate lines will show up as a thick, dashed, black line. 
showAggregateLine <- FALSE

##############################################################################################################
# Load Needed Functions - This creates all the needed functions. Run everything in this section.
# Note: If looking at functions, some make reference to a variable unlabeledAmino and this actually
# refers to any unlabeled compound but analysis was originally run using only aminos. 
##############################################################################################################

# takes output from generate_output (other script) and generates aggregated data split into a list of  
# labeled [1], unlabeled [2], and all compounds separated [3]
aggregateData <- function(turnoverData, valueToAggregateOn = 'ratio'){
  
  # Sums based on time, set, RT (scan), and compound
  aggregated <- aggregate(turnoverData$Count,
                          list(turnoverData$unlabeledAmino, turnoverData$RT, turnoverData$time, turnoverData$set),
                          sum, na.action(na.pass))
  
  # add names back in for easy referencing
  colnames(aggregated) <- c("unlabeledAmino", "RT", "time", "set", "total.count")
  
  # make aggregated a dataframe
  aggregated <- as.data.frame(aggregated)
  
  # joins aggregated data onto turnoverData (join from plyr library)
  joinedTable <- join(turnoverData, aggregated, by = c("unlabeledAmino", "RT", "time", "set"))
  
  # Removed rows where zeros were in the total count row. 
  joinedNoZeros <- joinedTable[joinedTable$total.count != 0, ]
  
  # remove bigger table to reduce memory usage as this object can get quite large depending on the dataset
  rm(aggregated)
  
  # Split data into labeled and unlabeled tables and joinedNoZeros has all data
  unlabeled <- joinedNoZeros[joinedNoZeros$name == joinedNoZeros$unlabeledAmino, ]
  labeled <- joinedNoZeros[joinedNoZeros$name != joinedNoZeros$unlabeledAmino, ]
  allData <- joinedNoZeros
  
  # if analysis is focusing on ratio of count/total count
  if(valueToAggregateOn == 'ratio'){
    
    # generate ratios (count/total.count)
    unlabeled$ratio <- unlabeled$Count/unlabeled$total.count
    labeled$ratio <- labeled$Count/labeled$total.count
    allData$ratio <- allData$Count/allData$total.count
    
    # Aggregate unlabeled based on unlabeled compound and time 
    aggUnlabeled <- aggregate(unlabeled$ratio, by = list(unlabeled$unlabeledAmino, unlabeled$time),
                              FUN = mean)
    colnames(aggUnlabeled) <- c("name", "time", "ratio")
    
    # Aggregate labeled based on unlabeled compound and time (this combines all labeled versions of a compound)
    aggLabeled <- aggregate(labeled$ratio, by = list(labeled$unlabeledAmino, labeled$time),
                            FUN = mean)
    colnames(aggLabeled) <- c("name", "time", "ratio")
    
    # Aggregate all based on compound name and time 
    aggAll <- aggregate(allData$ratio, by = list(allData$name, allData$time),
                        FUN = mean)
    colnames(aggAll) <- c("name", "time", "ratio")
    
    # creates list with data that is to be returned. 
    returnData <- list(aggLabeled, aggUnlabeled, aggAll)
    
    # name the list
    names(returnData) <- c('labeled', 'unlabeled', 'all')
    
  }
  
  # if analysis is focusing on the slope of regression analyses
  else if (valueToAggregateOn == 'slope'){
    
    # list out data for analysis
    dataList <- list(labeled, unlabeled, allData)
    
    # analysis for labeled [1] and unlabed [2] is the same and aggregates by unlabeled compound name (unlabeledAmino)
    for (d in 1:2){
      
      # assign data to current data variable
      currentData <- dataList[[d]]
      
      # set up slope column in data frame
      currentData$slope <- 0.0
      
      # for each unlabeled compound (unlabeledAmino)
      for (n in unique(currentData$unlabeledAmino)) {
        
        # take subset of data where the name is equal to the unlabeled compound
        nameSub <- currentData[currentData$unlabeledAmino == n, ]
        
        # for each time point in the subset
        for (t in unique(nameSub$time)) {
          
          # make a new subset of just those times
          timeSub <- nameSub[nameSub$time == t, ]
          
          # for each set in that subset
          for (s in unique(timeSub$set)) {
            
            # make a subset of just that set 
            setSub <- timeSub[timeSub$set == s, ]
            
            # run a regression on that set's data
            regression <- lm(Count ~ total.count, data = setSub)
            
            # store the slope value and assign it to every record that it applies to (correct unlabeled compound, time, and set). 
            slope <- regression$coefficients[[2]]
            currentData$slope[currentData$unlabeledAmino == n & currentData$time == t & currentData$set == s] <- slope 
            
          }
          
        }
        
      }
      
      # store data with slopes in dataList
      dataList[[d]] <- data.frame(currentData)
      
    }
    
    # Address all compound dataset.
    currentData <- dataList[[3]]
    currentData$slope <- 0.0
    
    # for each compound name
    for (n in unique(currentData$name)) {  
      
      nameSub <- currentData[currentData$name == n, ]
      
      # for each time point
      for (t in unique(nameSub$time)) {
        
        timeSub <- nameSub[nameSub$time == t, ]
        
        # for each set
        for (s in unique(timeSub$set)) {
          
          setSub <- timeSub[timeSub$set == s, ]
          
          # run the regression 
          regression <- lm(Count ~ total.count, data = setSub)
          
          # store the slope and assign it to the appropriate records
          slope <- regression$coefficients[[2]]
          currentData$slope[currentData$name == n & currentData$time == t & currentData$set == s] <- slope
          
        }
        
      }
      
    }
    
    # store new data frame in dataList
    dataList[[3]] <- data.frame(currentData)
    
    # rename elements of the dataList for easier access
    names(dataList) <- c('labeled', 'unlabeled', 'all')
    
    # Aggregate unlabeled based on name and time 
    aggUnlabeled <- aggregate(dataList$unlabeled$slope, by = list(dataList$unlabeled$unlabeledAmino, dataList$unlabeled$time),
                              FUN = mean)
    colnames(aggUnlabeled) <- c("name", "time", "slope")
    
    # Aggregate labeled based on name and time 
    aggLabeled <- aggregate(dataList$labeled$slope, by = list(dataList$labeled$unlabeledAmino, dataList$labeled$time),
                            FUN = mean)
    colnames(aggLabeled) <- c("name", "time", "slope")
    
    # Aggregate all based on name and time 
    aggAll <- aggregate(dataList$all$slope, by = list(dataList$all$name, dataList$all$time),
                        FUN = mean)
    colnames(aggAll) <- c("name", "time", "slope")
    
    # creates list with data that is to be returned. 
    returnData <- list(aggLabeled, aggUnlabeled, aggAll)
    
    # name the list
    names(returnData) <- c('labeled', 'unlabeled', 'all')
    
  }
  
  # input different from 'ratio' or 'slope' 
  else {
    
    # assign returnData to something to show that the input was invalid but that can still be returned
    returnData <- NaN
    
    # Tell user that the input was invalid.
    print('Please enter a valid valueToAggregateOn value. (either "ratio" or "slope")')
    
  }
  
  # return all of the data
  returnData
  
}

# takes input data transforms it so each row contains all time points for a compound
# Note, timeseries need to have the same number of points
makeLineData <- function(inputData, analysis = 'labeled'){
  
  # select the correct subset of the aggregated data
  selectedData <- inputData[[analysis]]
  
  # set up data holder for what will eventually be each plotted line
  lineData <- data.frame()
  
  # loop through each unique name in the selected data
  for (name in unique(selectedData$name)){
    
    # new data row for new line
    newLine <- data.frame(name, stringsAsFactors = FALSE)
    
    # subsets data for each compound
    compoundData <- selectedData[selectedData$name == name, ]
    
    # loop over each time and append the desired value to the new line
    for (value in compoundData[[3]]){
      
      newLine <- cbind(newLine, value, stringsAsFactors = FALSE)
      
    }
    
    # bind line to data holder
    lineData <- rbind(lineData, newLine, stringsAsFactors = FALSE)
  }
  
  # return line data to eventually be used for plotting
  lineData
  
}

# generates plot of bss and wss for given data to allow user to determine best k (k is the number of clusters)
compoundKFinder <- function(lineData, filename){
  
  # Set up lists of values for within and between sum of squares
  # wss measures cohesion within each cluster, bss measures separation of clusters
  wss <- c(rep(0, nrow(lineData)-1))
  bss <- c(rep(0, nrow(lineData)-1))
  
  # loop over every option for k
  for(i in 1:(nrow(lineData)-1)){
    
    # calculate k-means 
    kmeansTemp <- kmeans(lineData[2:ncol(lineData)], centers = i)
    
    # store results
    wss[i] <- kmeansTemp$tot.withinss
    bss[i] <- kmeansTemp$betweenss
    
  }
  
  # set up file so that the plot saves
  pdf(file = filename)
  
  # plot wss as a blue line
  plot(1:(nrow(lineData)-1), wss, type = "b", xlab = "Number of Clusters",
       ylab = "Distance", col="blue")
  
  # to overlay the bss plot
  par(new=T)
  
  # plot bss as a red line
  plot(1:(nrow(lineData)-1), bss, type = "b", axes = F, xlab = "",
       ylab = "", col = "red")
  
  # stop overlaying plots
  par(new=F)
  
  # Add in a legend
  legend("right",legend = c("BSS", "WSS"), title = "Separation and Cohesion", xpd = TRUE, horiz = FALSE, 
         col = c("red", "blue"), lty = 1)
  
  # add a title
  title(main = "BSS and WSS for All Possible K Values")
  
  # stop adding to saved plot
  dev.off()
  
  ## Show Plot - saved plot will not display to users so plot is regenerated in display 
  # plot wss as a blue line
  plot(1:(nrow(lineData)-1), wss, type = "b", xlab = "Number of Clusters",
       ylab = "Distance", col="blue")
  
  # overlay the bss plot
  par(new=T)
  
  # plot bss as a red line
  plot(1:(nrow(lineData)-1), bss, type = "b", axes = F, xlab = "",
       ylab = "", col = "red")
  
  # Add a legend
  legend("right",legend = c("BSS", "WSS"), title = "Separation and Cohesion", xpd = TRUE, horiz = FALSE, 
         col = c("red", "blue"), lty = 1)
  
  # Add a title
  title(main = "BSS and WSS for All Possible K Values")
  
  # stop overlaying plots
  par(new=F)
  
}

# shows hierarchical clustering and returns groupings if k provided
compoundHCluster <- function(lineData, kvalue = 0, clusterMethod = 'ward.D', distance = 'euclidean',
                             filename = NaN){
  
  # initialize groups
  groupings <- NaN
  
  # calculate distances between each line based on user supplied distance method. 'euclidean' is default. 
  # ((un)labeled compounds or all individual compounds depending on type of analysis)
  d <- dist(lineData[, 2:ncol(lineData)], method = distance)
  
  # cluster using user supplied cluster method. 'ward.D' is default if no method provided. 
  h <- hclust(d, method = clusterMethod)
  
  # if filename provided, save plot
  if(!is.nan(filename)){
    
    # create file to save plot to
    pdf(file = filename)
    
    # plot the dendrogram
    plot(h, labels = lineData$name, xlab = "Compounds")
    
    # if a kvalue has been provided and is greater than 0
    if(kvalue > 0) {
      
      # draws rectangles around each cluster
      rect.hclust(h, k=kvalue, border = "blue")
      
    }
    
    dev.off()
    
  }
  
  # show plot
  plot(h, labels = lineData$name, xlab = "Compounds")
  
  # if a kvalue has been provided and is greater than 0
  if(kvalue > 0) {
    
    # draws rectangles around each cluster
    rect.hclust(h, k=kvalue, border = "blue")
    
    # Stores the groupings to line data
    groupings <- cutree(h, k=kvalue)
    lineData$groups <- groupings
    
  }
  
  # returns the line data with groupings added if a k value was provided. 
  lineData
  
}

# generates the groups based on k means clustering and returns the lineData with cluster numbers appended and a separate dataframe of the centers. 
# providing a kvalue is NOT optional as it was with compoundHCluster
compoundKMeansCluster <- function(lineData, kvalue, clusteringAlgorithm) {
  
  # run k means clustering and store results. 
  partition <- kmeans(lineData[2:ncol(lineData)], centers = kvalue, algorithm = clusteringAlgorithm, nstart = 5)
  
  # assign groups
  lineData$groups <- partition$cluster
  
  # assign centers
  centers <- partition$centers
  
  # return the results
  list(lineData, centers)
  
}

# Takes aggregated data with clusters provided (from either compoundHCluster or compoundKMeansCluster)
# and generates plots for each group. Uses k-means clusters by default
# ** limited to < 74 compounds due to limitations on diverging colors for plotting. 
plotMatrix <- function(aggregatedClusteredData, typeOfClustering = 'k', filename, valueToAggregateOn = 'ratio',
                       plotAggregateLine = TRUE){
  
  # get k-means or hierarchical clustering data. k-means used by default
  if(typeOfClustering == 'h') {aggregatedClusterData <- aggregatedClusteredData[[2]]}
  else {aggregatedClusterData <- aggregatedClusteredData[[1]]}
  
  # get ranges for the axes 
  xrange <- range(aggregatedClusterData$time)
  yrange <- range(aggregatedClusterData[[valueToAggregateOn]])
  
  # get the number of groups that were assigned to the provided data
  numGroups <- length(unique(aggregatedClusterData$groups))
  
  # save plot
  pdf(file = filename)
  
  # set plotting area to be the correct size and number sub plots
  par(mfrow=c(ceiling(numGroups/2), 2))
  par(mar=c(5.1, 4.8, 4.1, 6.2))
  # note the mar sets the margins for the plotting area and this portion of the code is imperfect. With more than 6 plots the area 
  # and plots may look odd but can be corrected in a simple image editor such as inkscape. The legends also have a similar issue with the same solution.
  
  # get color spectrum for plotting
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # initialize color index variable
  c <- 1
  
  # for each group 
  for(i in 1:numGroups){
    
    # store group subset
    groupSub <- aggregatedClusterData[aggregatedClusterData$groups == i,]
    
    # set up plotting area
    plot(xrange, yrange, type = "n", xlab="Time (hrs)", 
         ylab = valueToAggregateOn, xpd = TRUE)
    
    # add a title
    title(paste("Group", i))
    
    # get number of lines to draw for that group
    numLines <- length(unique(groupSub$name))
    
    # set up ploting variables
    names <- unique(groupSub$name)
    subColors <- colors[c:(c+numLines)]
    
    # increment color index so no duplicates are used
    c <- numLines + c +1
    
    # add lines 
    for (j in 1:numLines) {
      
      # get compound name
      new_name <- names[j]
      
      # get compound data
      compound <- subset(groupSub, name == new_name)
      
      # plot compound with time as the x value and the desired value (ratio or slope) as the y
      lines(compound$time, compound[[valueToAggregateOn]], type = "b", lwd = 1.5, col=subColors[j])
      
    }
    
    # check if aggregate plot is desired
    if(plotAggregateLine) {
      
      # Plot aggregate line using LOESS approach (local regressions)
      loess_fit <- loess(groupSub[[valueToAggregateOn]] ~ groupSub$time)
      lines(groupSub$time, predict(loess_fit), col = "black", lty='dashed', lwd=3) # thick black line. Dashed. 
      
    }
    
    # plotting parameters for the legend
    cex_value <- 1.2/ceiling(numGroups/2)
    inset_value <- -.6 + (.16 * ceiling(numGroups/2))
    
    # add a legend, this should save as being off the plot to the right but still a bit buggy. Fix with inkscape or similar. 
    legend("right",legend = names, title = "Compounds", xpd = TRUE, horiz = FALSE, 
           col = subColors, lty = 1, inset = c(inset_value,0), cex = cex_value)
    
  }
  
  # turn plotting device off to finalize saved pdf file
  dev.off()
  
  ## redo plot to show plot to user now that a file has been saved
  # get ranges for the axes 
  xrange <- range(aggregatedClusterData$time)
  yrange <- range(aggregatedClusterData[[valueToAggregateOn]])
  
  # get the number of groups that will need to be plotted. 
  numGroups <- length(unique(aggregatedClusterData$groups))
  
  # set plotting area to be the correct size and number sub plots
  par(mfrow=c(ceiling(numGroups/2), 2))
  par(mar=c(5.1, 4.8, 4.1, 6.2))
  # note the mar sets the margins for the plotting area and this portion of the code is imperfect. With more than 6 plots the area 
  # and plots may look odd but can be corrected in a simple image editor such as inkscape. The legends also have a similar issue with the same solution.
  
  # color spectrum for plotting
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # initialize color index
  c <- 1
  
  # for each group in the data
  for(i in 1:numGroups){
    
    # generate group subset
    groupSub <- aggregatedClusterData[aggregatedClusterData$groups == i,]
    
    # set up plotting area (subplot)
    plot(xrange, yrange, type = "n", xlab="Time (hrs)", 
         ylab = valueToAggregateOn, xpd = TRUE)
    
    # add a title to the subplot
    title(paste("Group", i))
    
    # get number of lines to draw
    numLines <- length(unique(groupSub$name))
    
    # set up ploting variables
    names <- unique(groupSub$name)
    subColors <- colors[c:(c+numLines)]
    
    # increment color index for next pass through the loop
    c <- numLines + c +1
    
    # add lines 
    for (j in 1:numLines) {
      
      # get name
      new_name <- names[j]
      
      # additional subset for plotting each individual line
      compound <- subset(groupSub, name == new_name)
      
      # Plot lines
      lines(compound$time, compound[[valueToAggregateOn]], type = "b", lwd = 1.5, col=subColors[j])
      
    }
    
    # check if aggregate plot is desired
    if(plotAggregateLine) {
      
      # Plot aggregate line using LOESS approach (local regressions)
      loess_fit <- loess(groupSub[[valueToAggregateOn]] ~ groupSub$time)
      
      # line will be large, dashed black line. 
      lines(groupSub$time, predict(loess_fit), col = "black", lty='dashed', lwd=3) 
      
    }
    
    # plotting parameters
    cex_value <- 1.2/ceiling(numGroups/2)
    inset_value <- -.6 + (.16 * ceiling(numGroups/2))
    
    # add a legend, this should save as being off the plot to the right but still a bit buggy. Fix with inkscape or similar.  
    legend("right",legend = names, title = "Compounds", xpd = TRUE, horiz = FALSE, 
           col = subColors, lty = 1, inset = c(inset_value,0), cex = cex_value)
    
  }
}

# Takes arguments for the output filenames and also takes output generated from 
# the generate_output() function in the turnover script and the type of value to aggregate on. 'ratio' by default 
generateCSVs <- function(uFilename, lFilename, aFilename, valueToAggregateOn = 'ratio', turnoverData) {
  
  # Sums based on time, set, RT (scan), and compound
  aggregated <- aggregate(turnoverData$Count,
                          list(turnoverData$unlabeledAmino, turnoverData$RT, turnoverData$time, turnoverData$set),
                          sum, na.action(na.pass))
  
  # add names back in for easy referencing
  colnames(aggregated) <- c("unlabeledAmino", "RT", "time", "set", "total.count")
  
  # make aggregated a dataframe
  aggregated <- as.data.frame(aggregated)
  
  # joins aggregated data onto turnoverData then creates ratios
  joinedTable <- join(turnoverData, aggregated, by = c("unlabeledAmino", "RT", "time", "set"))
  
  # Removed rows where zeros were in the total count row. 
  joinedNoZeros <- joinedTable[joinedTable$total.count != 0, ]
  
  # remove bigger table to reduce memory usage
  rm(aggregated)
  
  # Split data into labeled, unlabeled, and allData tables
  unlabeled <- joinedNoZeros[joinedNoZeros$name == joinedNoZeros$unlabeledAmino, ]
  labeled <- joinedNoZeros[joinedNoZeros$name != joinedNoZeros$unlabeledAmino, ]
  allData <- joinedNoZeros
  
  # if desired analysis is to look at the ratio (vs slope)
  if(valueToAggregateOn == 'ratio'){
    
    # generate ratios
    unlabeled$ratio <- unlabeled$Count/unlabeled$total.count
    labeled$ratio <- labeled$Count/labeled$total.count
    allData$ratio <- allData$Count/allData$total.count
    
    # Aggregate unlabeled compounds based on unlabeled name and time 
    aggUnlabeled <- aggregate(unlabeled$ratio, by = list(unlabeled$unlabeledAmino, unlabeled$time, unlabeled$set),
                              FUN = mean)
    colnames(aggUnlabeled) <- c("name", "time", "set", "ratio")
    
    # Aggregate labeled compounds based on unlabeled name and time 
    aggLabeled <- aggregate(labeled$ratio, by = list(labeled$unlabeledAmino, labeled$time, labeled$set),
                            FUN = mean)
    colnames(aggLabeled) <- c("name", "time", "set", "ratio")
    
    # Aggregate all compounds based on specific compound name and time 
    aggAll <- aggregate(allData$ratio, by = list(allData$name, allData$time, allData$set),
                        FUN = mean)
    colnames(aggAll) <- c("name", "time", "set", "ratio")
    
    # creates list with data that is to be returned. 
    returnData <- list(aggLabeled, aggUnlabeled, aggAll)
    
    # name the list
    names(returnData) <- c('labeled', 'unlabeled', 'all')
    
  }
  
  # if analysis is focusing on the slope of regression analyses
  else if (valueToAggregateOn == 'slope'){
    
    # list out data for analysis
    dataList <- list(labeled, unlabeled, allData)
    
    # analysis for labeled and unlabed is the same
    for (d in 1:2){
      
      # assign data to current data variable
      currentData <- dataList[[d]]
      
      # set up slope column in data frame
      currentData$slope <- 0.0
      
      # for each unlabeled compound
      for (n in unique(currentData$unlabeledAmino)) {
        
        # take subset of data where the name is equal to the unlabeled compound
        nameSub <- currentData[currentData$unlabeledAmino == n, ]
        
        # for each time point
        for (t in unique(nameSub$time)) {
          
          # make a new subset of just those times
          timeSub <- nameSub[nameSub$time == t, ]
          
          # for each set
          for (s in unique(timeSub$set)) {
            
            # make a subset of just that set 
            setSub <- timeSub[timeSub$set == s, ]
            
            # run a regression on that sets data
            regression <- lm(Count ~ total.count, data = setSub)
            
            # store the slope value and assign it to every record that it applies to. 
            slope <- regression$coefficients[[2]]
            currentData$slope[currentData$unlabeledAmino == n & currentData$time == t & currentData$set == s] <- slope 
            
          }
          
        }
        
      }
      
      # store data with slopes in dataList
      dataList[[d]] <- data.frame(currentData)
      
    }
    
    # Hard coded which is bad.... but seems like the easiest way to do this.
    currentData <- dataList[[3]]
    currentData$slope <- 0.0
    
    # for each compound name
    for (n in unique(currentData$name)) {  
      
      nameSub <- currentData[currentData$name == n, ]
      
      # for each time point
      for (t in unique(nameSub$time)) {
        
        timeSub <- nameSub[nameSub$time == t, ]
        
        # for each set
        for (s in unique(timeSub$set)) {
          
          setSub <- timeSub[timeSub$set == s, ]
          
          # run the regression 
          regression <- lm(Count ~ total.count, data = setSub)
          
          # store the slope and assign it to the appropriate records
          slope <- regression$coefficients[[2]]
          currentData$slope[currentData$name == n & currentData$time == t & currentData$set == s] <- slope
          
        }
        
      }
      
    }
    
    # store new data frame in dataList
    dataList[[3]] <- data.frame(currentData)
    
    # rename elements of the dataList for easier access
    names(dataList) <- c('labeled', 'unlabeled', 'all')
    
    # Aggregate unlabeled based on name and time 
    aggUnlabeled <- aggregate(dataList$unlabeled$slope, by = list(dataList$unlabeled$unlabeledAmino, dataList$unlabeled$time,
                                                                  dataList$unlabeled$set), FUN = mean)
    colnames(aggUnlabeled) <- c("name", "time", "set", "slope")
    
    # Aggregate labeled based on name and time 
    aggLabeled <- aggregate(dataList$labeled$slope, by = list(dataList$labeled$unlabeledAmino, dataList$labeled$time,
                                                              dataList$labeled$set), FUN = mean)
    colnames(aggLabeled) <- c("name", "time", "set", "slope")
    
    # Aggregate all based on name and time 
    aggAll <- aggregate(dataList$all$slope, by = list(dataList$all$name, dataList$all$time, dataList$all$set),
                        FUN = mean)
    colnames(aggAll) <- c("name", "time", "set", "slope")
    
    # creates list with data that is to be returned. 
    returnData <- list(aggLabeled, aggUnlabeled, aggAll)
    
    # name the list
    names(returnData) <- c('labeled', 'unlabeled', 'all')
    
  }
  
  # sort output for easier use. Order by name
  returnData$labeled <- returnData$labeled[order(returnData$labeled$name), ]
  returnData$unlabeled <- returnData$unlabeled[order(returnData$unlabeled$name), ]
  returnData$all <- returnData$all[order(returnData$all$name), ]
  
  # write the CSVs
  write.csv(returnData['labeled'], file = lFilename, row.names = FALSE)
  write.csv(returnData['unlabeled'], file = uFilename, row.names = FALSE)
  write.csv(returnData['all'], file = aFilename, row.names = FALSE)
  
  # return data
  returnData
  
}

##############################################################################################################
# Workflow - Need to step through the sections this sequentially. 
##############################################################################################################

# set directory to output location
setwd(outputPathCluster)

# aggregate the output from the generate_output function in the metabolite_turnover script. Requires 
# the user to provide a value for aggre which should either be 'ratio' or 'slope'. 
aggregatedDataList <- aggregateData(out, clusterBy)

# transform data to make it usable for cluster analysis
labeledLines <- makeLineData(aggregatedDataList, typeOfAnalysis)

# store names and then remove incomplete data. 
writeLines(paste("removing", labeledLines[!complete.cases(labeledLines), 'name']))
incompletes <- labeledLines[!complete.cases(labeledLines), 'name']
labeledLines <- labeledLines[complete.cases(labeledLines), ]

# plots wss and bss for all possible k values (number of clusters)
# graph is interpreted by  evealuating where both lines reach an asymptote, this is the ideal k value. 
compoundKFinder(labeledLines, kfinderFilename)

##  ********************    THIS NEEDS TO BE ENTERED HERE AFTER DETERMINING K VALUE    ******************** ##

k <- 8

## ******************************************************************************************************** ##


# using k from compound k finder, plot hierarchical clustering
# and return the line data with groups appended. 
hClusters <- compoundHCluster(lineData = labeledLines, kvalue = k, clusterMethod =  hclustMethod, 
                              filename = dendroFilename, distance = pointDistance)

# run k-means clustering and get resulting groups. Output stored as a list with the linedata in the first
# slot and the centers (averages of each cluster) in the second slot
results <- compoundKMeansCluster(labeledLines, kvalue = k, clusteringAlgorithm = kmeansAlgorithm)

# store clusters and centers separately
kClusters <- as.data.frame(results[1]) 
centers <- as.data.frame(results[2])  

# add groups to time series (line format) data
groupedAggregatedData <- list(join(aggregatedDataList[[typeOfAnalysis]], kClusters[c("name", "groups")]),
                              join(aggregatedDataList[[typeOfAnalysis]], hClusters[c("name", "groups")]))

# remove incompletes
groupedAggregatedData[[1]] <- groupedAggregatedData[[1]][complete.cases(groupedAggregatedData[[1]]), ]
groupedAggregatedData[[2]] <- groupedAggregatedData[[2]][complete.cases(groupedAggregatedData[[2]]), ]

# generate plots by groupings (only showing for k-means results)
plotMatrix(groupedAggregatedData, typeOfClustering = typeOfClustering, filename = clusterFilename, valueToAggregateOn = clusterBy)

# reset plotting area
dev.off()

# (optional) Generate .csv files with raw data for both labeled and unlabeled datasets
generateCSVs(uFilename = unlabeledFilename, lFilename = labeledFilename,aFilename = allDataFilename,
             valueToAggregateOn = clusterBy, turnoverData = out)

#output overall results tables
write.csv(results[[1]], file='Average_Lines.csv')
write.csv(results[[2]], file='Consensus_Lines.csv')
