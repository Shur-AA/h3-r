% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{h3_bilinear_interpolation}
\alias{h3_bilinear_interpolation}
\title{Gets stars raster object, extracts values and interpolates them into
the net of H3 hexagons' centers which covers the raster
with bilinear interpolator from stars library}
\usage{
h3_bilinear_interpolation(raster_obj, hex_tab)
}
\arguments{
\item{hex_tab}{}
}
\value{
df with hex's centers coords and assigned values (z, x, y)
}
\description{
Gets stars raster object, extracts values and interpolates them into
the net of H3 hexagons' centers which covers the raster
with bilinear interpolator from stars library
}
