\name{bbsData}

\docType{data}

\alias{bbsData}

\title{Count data for six warbler species in Pennsylvania, USA}

\description{
Count data on 6 warblers in Pennsylvania, USA in 2018. Data come from the North American
Breeding Bird Survey. Data indicate the total numer of individuals for each of 6 species
counted at 50 stops along a 40km route (95 routes in the data set). The six species included 
in the data set are: (1) American Redstart (AMRE); (2) Blackburnian Warbler (BLBW); (3)
Black-throated Blue Warbler; (4) Black-throated Green Warbler; (5) Hooded Warbler; and
(6) Magnolia Warbler. Covariate data include three bioclimatic variables derived from 
PRISM. Two landcover variables (forest cover and developed land cover) come from USGS 
Earth Resources Observation and Science (EROS) Center.
}

\usage{
data(bbsData)
}

\source{U.S. Geological Survey. Downloaded from https://www.sciencebase.gov/catalog/item/52b1dfa8e4b0d9b325230cd9 on August 25, 2023.}

\references{
Daly, C., Halbleib, M., Smith, J. I., Gibson, W. P., Doggett, M. K., Taylor, G. H., Curtis, J., and
Pasteris, P. P. (2008). Physiographically sensitive mapping of climatological temperature and
precipitation across the conterminous united states. International Journal of Climatology: a
Journal of the Royal Meteorological Society, 28(15):2031–2064

Ziolkowski Jr., D.J., Lutmerding, M., English, W.B., Aponte, V.I., and Hudson, M-A.R., 2023, North American Breeding Bird Survey Dataset 1966 - 2022: U.S. Geological Survey data release, https://doi.org/10.5066/P9GS9K64.

Sohl, T., Reker, R., Bouchard, M., Sayler, K., Dornbierer, J., Wika, S., ... & Friesz, A. (2016). Modeled historical land use and land cover for the conterminous United States. Journal of Land Use Science, 11(4), 476-499.
}

\format{
  \code{bbsData} is a list with three elements: 

   \code{y}: a two-dimensional matrix of count data with rows corresponding
     to species (6) and columns corresponding to sites (95). 

   \code{covs}: a data frame with 95 rows and 8 columns consisting of 
     covariates for use in modeling relative abundance.

   \code{coords}: a numeric matrix with 95 rows and two columns containing the 
     site coordinates. The proj4string is "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
}

\keyword{datasets}

