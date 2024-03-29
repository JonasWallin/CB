% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{build_mesh_divergence}
\alias{build_mesh_divergence}
\title{Mesh construction for example on GP regression}
\usage{
build_mesh_divergence(
  grid_,
  expand,
  n,
  const.by = 3,
  soft_dirchlet = TRUE,
  soft_derivatieve = TRUE,
  alpha
)
}
\arguments{
\item{grid_}{(2 x 1) start end location of grid}

\item{n}{(int) number of mesh points}

\item{const.by}{- (int) often to evalute gradient}

\item{soft_dirchlet}{(bool) use to remove boundary effects}

\item{expan}{(double) expand the boundary to remove boundary effect}

\item{soft_derivative}{(bool) use soft constraint}
}
\value{
input_data (list) for likelhiood or mean construction
}
\description{
Function used to construct a mesh for the constrained GP regression
}
