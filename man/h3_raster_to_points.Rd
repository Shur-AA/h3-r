% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{h3_raster_to_points}
\alias{h3_raster_to_points}
\title{Gets centers of raster cells (stars object) and converts them into the
df with their coords and values}
\usage{
h3_raster_to_points(raster_obj)
}
\arguments{
\item{raster_obj}{}
}
\value{
rast_pnt_grid with columns: ('lon_p', 'lat_p', 'value')
}
\description{
Gets centers of raster cells (stars object) and converts them into the
df with their coords and values
}
