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

#' generate a list containg experimental variogram and fitted variogram.
#' @param target: string defining the response vector.
#' @param data: spatialdataframe including target
#' @param model: string represents model of variogram.
#' 
#' @return A list contains experimental variogram and fitted variogram.
#' 
generate_vgm <- function(target, data, name, crs, loc, ratio)
{
  relation <- as.formula(paste(target, '~', '1'))
  coordinates(data) <- loc
  proj4string(data) <- crs
  vertexes <- bbox(data)
  diagonal <- max(dist(t(as.matrix(vertexes))))
  vgm <- variogram(relation, data, cutoff = diagonal * ratio)
  fitted_vgm <- fit.variogram(vgm, vgm(name))
  model <- list("vgm" = vgm, "fitted_vgm" = fitted_vgm)
  return (model)
}

#' generate a list containg experimental variogram and fitted variogram.
#' @param source_list: composed of dataframes of measures.
#' @param rownames: dataframe row names of elements of list

#' 
#' @return A list contains dataframes composed of measures of various interpolation.
#' 
dataframe_list_convertion <- function(source_list, rownames)
{
  result <- list()
  if (check_list(source_list))
  {
    len <- length(source_list)
    rowcounts <- dim(source_list[[1]])[1]
    for (i in 1:rowcounts)
    {
      temp_dataframe <- data.frame()
      for (j in 1:len)
      {
        cur_row <- source_list[[j]][i,]
        # if (j > 1)
        #   names(cur_row) <- substring(rownames(cur_row), 1,
        #                               nchar(rownames(cur_row)))
        temp_dataframe <- rbind(temp_dataframe, cur_row)
      }
      #row_names <- cons_names(temp_dataframe, rownames)
      row_names <- rownames
      rownames(temp_dataframe) <- row_names
      result[[i]] <- temp_dataframe
      
    }
    return (result)
    
  }
  else
    return (NULL)

}

#' Check if the dimensions of elements of source_list are identical and if source_list
#' include element.
#' @param source_list: composed of dataframes of measures.
#' 
#' @return Boolean value.
#'
check_list <- function(source_list)
{
  len <- length(source_list)
  if (len <= 0)
    return (FALSE)
  else
  {
    dim_dataframe <- dim(source_list[[1]])
    for (i in 2:len)
    {
      dim_i = dim(source_list[[i]])
      if (sum(dim_dataframe == dim_i) != 2)
        return (FALSE)
    }
    return (TRUE)
  }
}

#' construct rownames of a dataframe 
#' @param dataframe: a dataframe of measures.
#' @param names: prefix of names of dataframe.
#' 
#' @return array of names.
#'
cons_names <- function(dataframe, names)
{
  len <- length(names)
  if (len != dim(dataframe)[1])
    return (NULL)
  row_names <- rownames(dataframe)
  for (i in 2:len)
  {
    row_names[i] <- substring(row_names[i], 1, nchar(row_names[i]) - 1)
  }
  return (paste(names, row_names))
  
}

#' barplot for dataframe elements of a list 
#' @param dataframe_list: a list of dataframes.
#'
measure_barplot <- function(dataframe_list, titles)
{
  element_counts <- length(dataframe_list)
  rows <- dim(dataframe_list[[1]])[1]
  for (i in 1:element_counts)
  {
    measures_matrix <- as.matrix(dataframe_list[[i]])
    peak <- max(measures_matrix)
    bplt <- barplot(measures_matrix, col = 1:rows, 
                    main = titles[[i]], beside = TRUE)
    legend("topright", legend = rownames(measures_matrix),
           fill = 1:rows, ncol = 2, cex = 0.75)
    ypos <- apply(measures_matrix, 1:2, function(x) if (x>0) x + 0.1*peak else x - 0.1*peak)
    text(x= bplt, y= ypos,
         labels=as.character(round(as.matrix(dataframe_list[[i]]), 2))
         , xpd=TRUE, cex = 0.75)  
  }
}

sampledist <- function(data, locname, bin, method = "euclidean", labels = TRUE)
{
  location = as.matrix(data[,locname])
  distance = dist(location, method = method)
  hist(distance, breaks = bin, labels = labels)
  return(distance)
}

