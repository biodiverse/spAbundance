\name{neonPredData}

\docType{data}

\alias{neonPredData}

\title{Land cover covariates and coordinates at a 1ha resolution across Disney Wilderness Preserve}

\description{
Land cover covariates (forest cover and grassland cover) extracted at 
a 1km resolution across the Disney Wildnerss Preserve for use in predicting density 
across the park. Land cover data come from USGS EROS.
} 

\usage{
data(neonPredData)
}

\source{USGS Earth Resources Observation and Science Ceneter https://www.usgs.gov/centers/eros}

\format{
  \code{neonPredData} is a data frame with four columns: 

  \code{forest}: proportion of forest cover within 1km radius.

  \code{grass}: proportion of grassland cover within 1km radius.

  \code{Easting}: the x coordinate of the point. The EPSG is 32617 (UTM Zone 17N).

  \code{Northing}: the y coordinate of the point. The EPSG is 32617 (UTM Zone 17N).
}

\keyword{datasets}


