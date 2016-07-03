#' Return sequence numbers of records that are selceted as validation.
#' 
#' @param total Total number of records.
#' @param group fold number for cross validation.
#' @param seq which fold.
#' @return A vector consists of number of records.
validition_range <- function(total, group, seq)
{
  size = round(total/group)
  end = seq * size
  start = (seq - 1) * size + 1
  if (seq == group)
    return(start:total)
  else
    return(start:end)
}


#' compute statistics based on values of prediction and observation.
#' 
#' @param prediction denote values of prediction.
#' @param observation denote ture values
#' @return A list consist of statistics.
measures <- function(prediction, observation)
{
  difference = prediction - observation
  RMSE = sqrt(sum((difference * difference)/length(difference)))
  ME = sum(difference)/length(difference)
  MAE = sum(abs(difference))/length(difference)
  result = list()
  result$RMSE = RMSE
  result$MAE = MAE
  result$ME = ME
  return (result)
}


#' k-fold cross validation for inverst distance weighting.
#' 
#' @param target value to be estimated.
#' @param data A data.frame.
#' @param neighbors the number of data points used to estimate target
#' @param powers A value of power
#' @param crs A Coordinates reference system
#' @param loc A vector indicate column names representing spatial information
#' @param k indicates the number of folds
#' @return A list consist of average value of RMSE, MAE and ME.
idw.kfold <- function(target, data, neighbors, powers, 
                      crs, loc, k)
{
  
  accumulation = list("RMSE" = 0, "MAE" = 0, "ME" = 0)
  relation = as.formula(paste(target, "~", "1"))
  for (seq in 1:k)
  {
    data$tag = 0
    selection = validition_range(dim(data)[1], k, seq)
    
    data[selection,]$tag = 1
    train.data = data[data$tag == 0,]
    validation = data[data$tag == 1,]
    coordinates(train.data)<-loc
    coordinates(validation)<-loc
    proj4string(train.data)<-crs
    proj4string(validation)<-crs
    prediction = idw(relation, train.data, validation, idp = powers,
                     nmax = neighbors)
    measure = measures(prediction$var1.pred, as.data.frame(validation)[,target])
    len = length(accumulation)
    for (i in 1:len)
    {
      accumulation[[i]] = accumulation[[i]] + measure[[i]]
    }
  }
  for (i in 1:len)
  {
    accumulation[[i]] = accumulation[[i]]/k
  }
  
  return(accumulation)
  
}


#' k-fold cross validation for ordinary kriging.
#' 
#' @param target value to be estimated.
#' @param data A data.frame.
#' @param neighbors the number of data points used to estimate target
#' @param vgm variogram model
#' @param crs A Coordinates reference system
#' @param loc A vector indicate column names representing spatial information
#' @param k indicates the number of folds
#' @return A list consist of average value of RMSE, MAE and ME.
krige.kfold <- function(target, data, neighbors, vgm,
                        crs, loc, k)
{
  
  accumulation = list("RMSE" = 0, "MAE" = 0, "ME" = 0)
  relation = as.formula(paste(target, "~", "1"))
  for (seq in 1:k)
  {
    data$tag = 0
    selection = validition_range(dim(data)[1], k, seq)
    
    data[selection,]$tag = 1
    train.data = data[data$tag == 0,]
    validation = data[data$tag == 1,]
    coordinates(train.data)<-loc
    coordinates(validation)<-loc
    proj4string(train.data)<-crs
    proj4string(validation)<-crs
    prediction = krige(relation, train.data, validation,
                     nmax = neighbors, model = vgm)
    measure = measures(prediction$var1.pred, as.data.frame(validation)[,target])
    len = length(accumulation)
    for (i in 1:len)
    {
      accumulation[[i]] = accumulation[[i]] + measure[[i]]
    }
    
    
  }
  for (i in 1:len)
  {
    accumulation[[i]] = accumulation[[i]]/k
  }
  
  return(accumulation)
  
}


#' k-fold cross validation for trend surface analysis.
#' 
#' @param target value to be estimated.
#' @param data A data.frame.
#' @param degree of trend surface
#' @param crs A Coordinates reference system
#' @param loc A vector indicate column names representing spatial information
#' @param k indicates the number of folds
#' @return A list consist of average value of RMSE, MAE and ME.
surf.kfold <- function(target, data, degree,
                       crs, loc, k)
{
  
  accumulation = list("RMSE" = 0, "MAE" = 0, "ME" = 0)
  relation = as.formula(paste(target, "~", "1"))
  for (seq in 1:k)
  {
    data$tag = 0
    selection = validition_range(dim(data)[1], k, seq)
    
    data[selection,]$tag = 1
    train.data = data[data$tag == 0,]
    validation = data[data$tag == 1,]
    coordinates(train.data)<-loc
    coordinates(validation)<-loc
    proj4string(train.data)<-crs
    proj4string(validation)<-crs
    prediction = krige(relation, train.data, validation,
                       degree = degree)
    measure = measures(prediction$var1.pred, as.data.frame(validation)[,target])
    len = length(accumulation)
    for (i in 1:len)
    {
      accumulation[[i]] = accumulation[[i]] + measure[[i]]
    }
    
    
  }
  for (i in 1:len)
  {
    accumulation[[i]] = accumulation[[i]]/k
  }
  
  return(accumulation)
  
}


average_measures <- function(measures)
{
  count = length(measures)
  
  accumulation = list("RMSE" = 0, "MAE" = 0, "ME" = 0)
  for (i in 1:count)
  {
    accumulation$RMSE = accumulation$RMSE + measures[[i]]$RMSE
    accumulation$MAE = accumulation$MAE + measures[[i]]$MAE
    accumulation$ME = accumulation$ME + measures[[i]]$ME
    
  }
  
  accumulation$RMSE = accumulation$RMSE/count
  accumulation$MAE = accumulation$MAE/count
  accumulation$ME = accumulation$ME/count
  
  return (accumulation)
}


#' select parameters using k-fold cross validation for 
#' inverse distance weighting.
#' 
#' @param target value to be estimated.
#' @param data A data.frame.
#' @param neighbors A vector whose elementes are the number of data points
#'  used to estimate target
#' @param powers A vector whose elementes are power value of power functions
#' @param crs A Coordinates reference system
#' @param loc A vector indicate column names representing spatial information
#' @param k indicates the number of folds
#' @return A list whose elements are list including optimal values of 
#' neighbor and power coresponding to constraints RMSE, MAE and ME.

idw.para <- function(target, data, neighbors, powers,
                     crs, loc, k)
{
  
  first = TRUE
  optimal = list()
  paraset = list()

  
  for (neighbor in neighbors)
  {
    for (power in powers)
    {
      error = idw.kfold(target,data, neighbor, power,
                        crs, loc, k)
      if (first)
      {
        first = FALSE
        optimal = error
        paraset$RMSE = list("neighbor" = neighbor, "power" = power)
        paraset$MAE = list("neighbor" = neighbor, "power" = power)
        paraset$ME = list("neighbor" = neighbor, "power" = power)
      }
      else
      {

        if (error$RMSE<optimal$RMSE)
        {
          optimal$RMSE = error$RMSE
          paraset$RMSE$neighbor = neighbor
          paraset$RMSE$power = power
        }
        if (error$MAE<optimal$MAE)
        {
          optimal$MAE = error$MAE
          paraset$MAE$neighbor = neighbor
          paraset$MAE$power = power
        }
        if (abs(error$ME)<abs(optimal$ME))
        {
          optimal$ME = error$ME
          paraset$ME$neighbor = neighbor
          paraset$ME$power = power
        }

      }
    }
  }
  return(paraset)

}

#' select parameters using k-fold cross validation for kriging
#' 
#' @param target value to be estimated.
#' @param data A data.frame.
#' @param neighbors A vector whose elementes are the number of data points
#'  used to estimate target
#' @param vgm Semi-variogram for spatial distribution
#' @param crs A Coordinates reference system
#' @param loc A vector indicate column names representing spatial information
#' @param k indicates the number of folds
#' @return A list whose elements are list including optimal values of 
#' neighbor coresponding to constraints RMSE, MAE and ME.

krige.para <- function(target, data, neighbors, vgm,
                     crs, loc, k)
{
  
  first = TRUE
  optimal = list()
  paraset = list()
  
  
  for (neighbor in neighbors)
  {
      error = krige.kfold(target,data, neighbor, vgm,
                        crs, loc, k)
      if (first)
      {
        first = FALSE
        optimal = error
        paraset$RMSE = list("neighbor" = neighbor)
        paraset$MAE = list("neighbor" = neighbor)
        paraset$ME = list("neighbor" = neighbor)
      }
      else
      {
        
        if (error$RMSE<optimal$RMSE)
        {
          optimal$RMSE = error$RMSE
          paraset$RMSE$neighbor = neighbor
        }
        if (error$MAE<optimal$MAE)
        {
          optimal$MAE = error$MAE
          paraset$MAE$neighbor = neighbor
        }
        if (abs(error$ME)<abs(optimal$ME))
        {
          optimal$ME = error$ME
          paraset$ME$neighbor = neighbor
        }

    }
  }
  return(paraset)
  
}


#' select parameters using k-fold cross validation for
#' trend surface analysis
#' 
#' @param target value to be estimated.
#' @param data A data.frame.
#' @param degrees A vector whose elementes are dependent degree of
#'  spatial location.
#' @param crs A Coordinates reference system
#' @param loc A vector indicate column names representing spatial information
#' @param k indicates the number of folds
#' @return A list whose elements are list including optimal values of 
#' neighbor coresponding to constraints RMSE, MAE and ME.

surf.para <- function(target, data, degrees,
                       crs, loc, k)
{
  
  first = TRUE
  optimal = list()
  paraset = list()
  
  
  for (degree in degrees)
  {
      error = surf.kfold(target,data, degree,
                          crs, loc, k)
      if (first)
      {
        first = FALSE
        optimal = error
        paraset$RMSE = list("degree" = degree)
        paraset$MAE = list("degree" = degree)
        paraset$ME = list("degree" = degree)
      }
      else
      {
        
        if (error$RMSE<optimal$RMSE)
        {
          optimal$RMSE = error$RMSE
          paraset$RMSE$degree = degree
        }
        if (error$MAE<optimal$MAE)
        {
          optimal$MAE = error$MAE
          paraset$MAE$degree = degree
        }
        if (abs(error$ME)<abs(optimal$ME))
        {
          optimal$ME = error$ME
          paraset$ME$degree = degree
        }
        
      }
  }
  return(paraset)
  
}

#' Convert list of parameters into dataframe of parameters.
#' Constraints: Elements of the list are identical
#' @param para: list of parameters.
#' @return A dataframe consists parameters.
para_dataframe <- function(para)
{
  len <- length(para)
  for (i in 1:len)
  {
    if (i == 1)
      para_dataframe <- as.data.frame(para[[i]])
    else
      para_dataframe <- rbind(para_dataframe, as.data.frame(para[[i]]))
  }
  rownames(para_dataframe) <- names(para)
  return (para_dataframe)

}

#' Convert list of parameters into dataframe of parameters.
#' Constraints: Elements of the list are identical
#' @param para: list of parameters.
#' @return A dataframe consists parameters.

sampledist <- function(data, locname, bin, method = "euclidean", labels = TRUE)
{
  location = as.matrix(data[,locname])
  distance = dist(location, method = method)
  hist(distance, breaks = bin, labels = labels)
  return(distance)
}