% Generated by roxygen2 (4.0.1): do not edit by hand
\name{eye.foodweb.random}
\alias{eye.foodweb.random}
\title{Builds a food web according to the Erdos-Renyi random graph model}
\usage{
eye.foodweb.random(size = 100, connectance = 0.25)
}
\arguments{
\item{size}{The number of species (integer)}

\item{connectance}{The desired connectance (between 0 and 1)}
}
\value{
A connected (one piece) food web of connectance close to the desired value.
The food web is a list containing:
\describe{
\item{$links}{The trophic links from resource to consumer (L x 2 matrix)}
\item{$size}{The number of species}
\item{$L}{The number of links}
\item{$connectance}{The connectance L / (size choose 2)}
\item{$model}{A code for the model used to generate the food web}
\item{$cycles}{Boolean, whether cycles are present}
}
All the self-loops and double arrows (a eats b, b eats a) are removed. If the food web contains
no cycles, then it is sorted that all the coefficients would be in the upper-triangular part
of the corresponding adjacency matrix.

If the resulting network has connectance too different from the desired one, or it is not connected,
the function will call itself recursively. If too many recursions are performed, it will fail.
}
\description{
Any two species are connected with probability connectance
}
\examples{
FW <- eye.foodweb.random() # size = 100, connectance = 0.25
FW <- eye.foodweb.random(30, 0.2)
FW <- eye.foodweb.random(size = 250, connectance = 0.314)
}

