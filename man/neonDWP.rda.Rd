\name{neonDWP}

\docType{data}

\alias{neonDWP}

\title{Distance sampling data of 16 bird species  
  observed in the Disney Wilderness Preserve in 2018 in Florida, USA}

\description{
Distance sampling data of 16 bird species in 2018 in the Disney Wilderness Preserve 
in Florida, USA. These data were collected as part of the National Ecological 
Observatory Network (NEON). Data were collected at 90 sites where observers recorded
the number of all bird species observed during a six minute, unlimited radius point count survey once
during the breeding season. Distance of each individual bird to the observer was 
recorded using a laser rangefinder. For modeling, we binned the distance measurements into 
4 distance bins, and only used observations within 250m. The 16 species included in the data set are as follows: 
(1) EATO (Eastern Towhee); (2) EAME (Eastern Meadowlark); (3) AMCR (American Crow);
(4) BACS (Bachman's Sparrow); (5) CARW (Carolina Wren); (6) COGD (Common Ground Dove); 
(7) CONI (Common Nighthawk); (8) COYE (Common Yellowthroat); (9) EABL (Eastern Bluebird);
(10) GCFL (Great-crested Flycatcher); (11) MODO (Mourning Dover); (12) NOCA (Northern Cardinal);
(13) NOMO (Northern Mockingbird); (14) RBWO (Red-bellied Woodpecker); (15) RHWO (Red-headed Woodpecker); 
(16) WEVI (White-eyed Vireo).
}

\usage{
data(neonDWP)
}

\source{NEON (National Ecological Observatory Network). Breeding landbird point 
	counts (DP1.10003.001), RELEASE-2023. https://doi.org/10.48443/00pg-vm19. 
	Dataset accessed from https://data.neonscience.org on May 25, 2023}

\references{
  Barnett, D. T., Duffy, P. A., Schimel, D. S., Krauss, R. E., Irvine, K. M., 
  Davis, F. W.,Gross, J. E., Azuaje, E. I., Thorpe, A. S., Gudex-Cross, D., et al. (2019). 
  The terrestrial organism and biogeochemistry spatial sampling design for 
  the national ecological observatory network. Ecosphere, 10(2):e02540.
}

\format{
  \code{neonDWP} is a list with five elements: 

   \code{y}: a three-dimensional array of distance sampling data with 
     dimensions of species (16), sites (90) and distance bin (4). 

   \code{covs}: a data frame with 90 rows and four columns consisting of 
     covariates for use in modeling abundance and/or detection.

   \code{dist.breaks}: a vector of five values indicating the break points 
     of the four distance bands in the data.

   \code{offset}: an offset used to scale the 250m radius point count surveys 
     to ha, such that resulting estimates will be the number of individuals 
     per ha.

   \code{coords}: a numeric matrix with 373 rows and two columns containing the 
     site coordinates (Easting and Northing) in UTM Zone 17N. The EPSG is 32617. 
}

\keyword{datasets}

