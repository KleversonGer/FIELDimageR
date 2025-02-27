#' fieldCountWithRealArea
#' 
#' @title Calculating number of objects per plot
#' 
#' @description The mask from function \code{\link{fieldMask}} is used to identify the number of objects per plot..
#' 
#' @param mosaic object mask of class stack from the function \code{\link{fieldMask}}.
#' @param fieldShape plot shape file.
#' @param value referent value to vegetation pixels in the mask. If "HUE" was used on 
#'   \code{\link{fieldMask}} the value=0.
#' @param minSize used to set the minimum size percentage of plant canopy  (to remove weeds and more).
#' @param n.core number of cores to use for multicore processing (Parallel). 
#' @param pch point symbol, please check \code{help("points")}.
#' @param cex character (or symbol) expansion: a numerical vector, please check \code{help("points")}.
#' @param col color code or name, please check \code{help("points")}.
#' @param na.rm logical. Should missing values (including NaN) be removed?.
#' 
#' @importFrom raster projection extract xyFromCell
#' @importFrom EBImage distmap watershed
#'
#' 
#' @return A list with five element
#' \itemize{
#'   \item \code{fieldCount} is the number of objects per plot represented in DataFrame.
#'   \item \code{fieldShape} is the new shapeFile with stand count.
#'   \item \code{mosaic} is the Watershed layer.
#'   \item \code{objectSel} is the objects area per plot.
#'   \item \code{objectReject} is the objects position in the image.
#' }
#' 
#'
#' @export

fieldCountWithRealArea <- function(mosaic, fieldShape, value = 0, minSize = 0.01, n.core = NULL, pch = 16, 
                       cex = 0.7, col = "red", na.rm = FALSE) {
  if(!is.na(projection(fieldShape))&is.na(projection(mosaic))){
    if(projection(fieldShape)!=projection(mosaic)){stop("fieldShape and mosaic must have the same projection CRS, strongly suggested to use fieldRotate() for both files.")}}
  mosaic <- stack(mosaic)
  num.band<-length(mosaic@layers)
  if(num.band>1){stop("Only mask with values of 1 and 0 can be processed, use the mask output from fieldMask()")}
  if(!value%in%c(1,0)){stop("Values in the mask must be 1 or 0 to represent the objects, use the mask output from fieldMask()")}
  if(!all(c(raster::minValue(mosaic),raster::maxValue(mosaic))%in%c(1,0))){stop("Values in the mask must be 1 or 0 to represent the objects, use the mask output from fieldMask()")}
  mosaic <- crop(x = mosaic, y = fieldShape)
  print("Identifying objects... ")
  par(mfrow=c(1,1))
  raster::plot(mosaic, col=grey(1:100/100), axes=FALSE, box=FALSE, legend=FALSE)
  sp::plot(fieldShape, add=T)
  names(mosaic)<-"mask"
  if(na.rm){
    mosaic$mask[is.na(mosaic$mask)] <- c(0,1)[c(0,1)!=value] 
  }
  mask <- raster::as.matrix(mosaic$mask) == value
  dd <- distmap(mask)
  mosaic$watershed <- watershed(dd)
  if(is.null(n.core)){extM <- extract(x = mosaic$watershed, y = fieldShape)}
  if (!is.null(n.core)){
    if(n.core>detectCores()){stop(paste(" 'n.core' must be less than ",detectCores(),sep = ""))}
    cl <- parallel::makeCluster(n.core, output = "", setup_strategy = "sequential")
    registerDoParallel(cl)
    extM <- foreach(i = 1:length(fieldShape), .packages = c("raster")) %dopar% {
      single <- fieldShape[i, ]
      CropPlot <- crop(x = mosaic$watershed, y = single)
      extract(x = CropPlot, y = single)
    }
    names(extM) <- 1:length(fieldShape)
    parallel::stopCluster(cl)
  }
  objects<-lapply(extM, function(x){table(x)})
  cent <- lapply(objects, function(x){as.numeric(names(x))[-1]})
  objectsPosition<- lapply(cent, function(x){
    if(length(x)==0){return(NULL)}
    pos<-NULL
    for(i in 1:length(x)){pos<-rbind(pos,colMeans(xyFromCell(mosaic$watershed, which(mosaic$watershed[]==x[i]))))}
    if(abs(max(pos[,1])-min(pos[,1]))>=abs(max(pos[,2])-min(pos[,2]))){ord<-order(pos[,1])}
    if(abs(max(pos[,1])-min(pos[,1]))<abs(max(pos[,2])-min(pos[,2]))){ord<-order(pos[,2])}
    pos<-pos[ord,]
    return(list(seqName=ord,Position=pos))})
  objectSel<-list()
  objectReject<-list()
  for(j in 1:length(cent)){
    x1<-objects[[j]]
    y<-objectsPosition[[j]]
    x<-NULL
    PS<-NULL
    PR<-NULL
    if(!is.null(y)){
      for (i in 2:length(x1)) {x<-c(x,round(x1[i],3))}
      x<-x[y$seqName]
      
      if(dim(as.matrix(y$Position))[2]==1){
        PS <- data.frame(objectArea = x[x >= minSize], x = y$Position[1][x >= minSize], y = y$Position[2][x >= minSize])
        PR <- data.frame(objectArea = x[x < minSize], x = y$Position[1][x < minSize], y = y$Position[2][x < minSize])
      }
      if(dim(as.matrix(y$Position))[2]!=1){
        PS <- data.frame(objectArea = x[x >= minSize], 
                         x = y$Position[x >= minSize, 1], 
                         y = y$Position[x >= minSize, 2])
        PR <- data.frame(objectArea = x[x < minSize], 
                         x = y$Position[x < minSize, 1], 
                         y = y$Position[x < minSize, 2])
      }}
    rownames(PS)<-NULL
    rownames(PR)<-NULL
    objectSel[[j]]<-PS
    objectReject[[j]]<-PR
  }
  if(length(objectSel)!=length(cent)){
    objectSel[[length(cent)+1]] <- NA
    objectReject[[length(cent)+1]] <- NA
    objectSel[[length(cent)+1]] <- NULL
    objectReject[[length(cent)+1]] <- NULL
  }
  field <- unlist(lapply(objectSel, function(x){length(x$objectArea)}))
  fieldShape@data$fieldCount <- field
  print(paste("Number of objects: ", sum(field), sep = ""))
  graphics::points(do.call(rbind,objectSel)[,c(2,3)], pch=pch, cex=cex, col=col)
  sp::plot(fieldShape, add=T)
  Out <- list(fieldCount=field, 
              fieldShape=fieldShape, 
              mosaic=mosaic$watershed, 
              objectSel=objectSel, 
              objectReject=objectReject)
  return(Out)
}
# 