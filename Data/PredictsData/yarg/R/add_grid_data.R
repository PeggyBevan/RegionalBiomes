AddGridData<-function(gridData,dataFrame,columnName,silent=FALSE){

  if (class(gridData)!="SpatialGridDataFrame"){
    stop("Error: gridData needs to be a SpatialGridDataFrame object")
  }
  if (class(dataFrame)!="data.frame"){
    stop("Error: dataFrame needs to be a data frame")
  }
  if (class(columnName)!="character"){
    stop("Error: columnName needs to be a character string")
  }
  
  # Get the cell sizes of the grid
  cell.size.y<-slot(slot(gridData,"grid"),"cellsize")[1]
  cell.size.x<-slot(slot(gridData,"grid"),"cellsize")[2]
  
  # Get the longitude and latitude of the top left of the grid 
  long.top.left<-bbox(gridData)[1,1]
  lat.top.left<-bbox(gridData)[2,2]
  
  # Get the numbers of cells in each dimension of the grid
  ncols<-slot(slot(gridData,"grid"),"cells.dim")[1]
  nrows<-slot(slot(gridData,"grid"),"cells.dim")[2]
  
  # Make a matrix containing the grid data
  gridData.matrix<-matrix(gridData$band1,nrow=nrows,ncol=ncols,byrow=T)
  
  # Add a new column to the data frame for the grid data
  dataFrame$newCol<-NA
  # Loop over rows in the data frame
  if(!silent) cat("Getting grid data\n")
  for (r in 1:dim(dataFrame)[1]){
    if (!silent) cat(paste("\rGetting grid data (",columnName,"): row ",r," of ",dim(dataFrame)[1],sep=""))
    
    # Get the coordinates of the cell
    long<-dataFrame$Longitude[r]
    lat<-dataFrame$Latitude[r]
    
    # Find the appropriate cell in the grid
    # Create a running longitude value, initially set as the longitude of the top-left corner
    long.search<-long.top.left
    # Create a variable for finding the longitudinal cell index
    long.cell<-0
    while(long.search<long)
    {
      long.search<-long.search+cell.size.x
      long.cell<-long.cell+1
    }
    
    # Create a running latitude value, initially set as the latitude of the top-left corner
    lat.search<-lat.top.left
    # Create a variable for finding the latitudinal cell index
    lat.cell<-0
    while(lat.search>lat)
    {
      lat.search<-lat.search-cell.size.x
      lat.cell<-lat.cell+1
    }
    
    # Add data from the appropriate cell in the grid to the data frame
    dataFrame$newCol[r]<-gridData.matrix[lat.cell,long.cell]
    
    
  }
  
  if (!silent) cat("\n")
  
  # Assign the specified name to the new column
  names(dataFrame)[length(names(dataFrame))]<-columnName
  return(dataFrame)
}

AddGridDataYear<-function(gridDataDir,dataFrame,columnName,silent=FALSE){
  
  if (class(dataFrame)!="data.frame"){
    stop("Error: dataFrame needs to be a data frame")
  }
  if (class(columnName)!="character"){
    stop("Error: columnName needs to be a character string")
  }
  
  # Find out which years are required
  dataFrame$EndYear<-as.integer(format(dataFrame$Sample_end_latest, "%Y"))
  years<-unique(dataFrame$EndYear)
  
  yr.data<-lapply(as.list(years),function(x) readGDAL(paste(gridDataDir,"/",x,".asc",sep=""),
                                                      silent=TRUE))
  names(yr.data)<-years
  
  # TODO - add check that grids are identical for all years
  
  # Get the cell sizes of the grid
  cell.size.y<-slot(slot(yr.data[[1]],"grid"),"cellsize")[1]
  cell.size.x<-slot(slot(yr.data[[1]],"grid"),"cellsize")[2]
  
  # Get the longitude and latitude of the top left of the grid 
  long.top.left<-bbox(yr.data[[1]])[1,1]
  lat.top.left<-bbox(yr.data[[1]])[2,2]
  
  # Get the numbers of cells in each dimension of the grid
  ncols<-slot(slot(yr.data[[1]],"grid"),"cells.dim")[1]
  nrows<-slot(slot(yr.data[[1]],"grid"),"cells.dim")[2]
  
  # Make a matrix containing the grid data
  yr.data.mat<-lapply(yr.data,function(x) return(matrix(x$band1,nrow=nrows,ncol=ncols,byrow=T)))
  
  # Add a new column to the data frame for the grid data
  dataFrame$newCol<-NA
  # Loop over rows in the data frame
  if(!silent) cat("Getting grid data\n")
  for (r in 1:dim(dataFrame)[1]){
    if (!silent) cat(paste("\rGetting grid data (",columnName,"): row ",r," of ",dim(dataFrame)[1],sep=""))
    
    # Get the coordinates of the cell
    long<-dataFrame$Longitude[r]
    lat<-dataFrame$Latitude[r]
    
    # Find the appropriate cell in the grid
    # Create a running longitude value, initially set as the longitude of the top-left corner
    long.search<-long.top.left
    # Create a variable for finding the longitudinal cell index
    long.cell<-0
    while(long.search<long)
    {
      long.search<-long.search+cell.size.x
      long.cell<-long.cell+1
    }
    
    # Create a running latitude value, initially set as the latitude of the top-left corner
    lat.search<-lat.top.left
    # Create a variable for finding the latitudinal cell index
    lat.cell<-0
    while(lat.search>lat)
    {
      lat.search<-lat.search-cell.size.x
      lat.cell<-lat.cell+1
    }
    
    # Add data from the appropriate cell in the grid to the data frame
    dataFrame$newCol[r]<-yr.data.mat[[paste(dataFrame$EndYear[r])]][lat.cell,long.cell]
    
    
  }
  
  cat("\n")
  
  # Assign the specified name to the new column
  names(dataFrame)[length(names(dataFrame))]<-columnName
  return(dataFrame)
}
