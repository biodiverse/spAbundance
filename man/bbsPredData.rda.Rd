\name{bbsPredData}

\docType{data}

\alias{bbsPredData}

\title{Covariates and coordinates for prediction of relative warbler abundance in Pennsylvania, USA}

\description{
Bioclimatic and land cover variables extracted at a 12km resolution across 
the state of Pennsylvania, USA for use in predicting relative abundance of six 
warbler species across the state. Land cover data come from USGS EROS, while 
climate data come from PRISM.
} 

\usage{
data(bbsPredData)
}

\references{
Daly, C., Halbleib, M., Smith, J. I., Gibson, W. P., Doggett, M. K., Taylor, G. H., Curtis, J., and
Pasteris, P. P. (2008). Physiographically sensitive mapping of climatological temperature and
precipitation across the conterminous united states. International Journal of Climatology: a
Journal of the Royal Meteorological Society, 28(15):2031–2064

Sohl, T., Reker, R., Bouchard, M., Sayler, K., Dornbierer, J., Wika, S., ... & Friesz, A. (2016). Modeled historical land use and land cover for the conterminous United States. Journal of Land Use Science, 11(4), 476-499.
}

\format{
  \code{bbsPredData} is a data frame with seven columns: 

  \code{bio2}: bioclim variable 2.
  
  \code{bio8}: bioclim variable 8.

  \code{bio18}: bioclim variable 18.

  \code{forest}: proportion of forest cover within a 5km radius.

  \code{devel}: proportion of developed land cover within a 5km radius.

  \code{x}: the x coordinate of the point. The proj4strig is "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs".

  \code{y}: the y coordinate of the point. The proj4string is "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs".
}

\keyword{datasets}


