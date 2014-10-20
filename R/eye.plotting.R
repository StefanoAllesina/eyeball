eye.ellipse.df <- function(centerx, radiusx, radiusy){
  thetas <- seq(pi / 2.0, 0.0, length = 1000)
  xbase <- radiusx * cos(thetas)
  ybase <- radiusy * sin(thetas)
  x <- c(xbase, rev(xbase), -xbase, rev(-xbase))
  y <- c(ybase, rev(-ybase), -ybase, rev(ybase))
  return(data.frame(Real = x + centerx, Imaginary = y))
}


eye.plot.approximations <- function(Approx, only.bulk = FALSE){
  ## Check if the eigenvalues of M have been computed already
  if (is.null(Approx$M.eigenvalues)){
    Approx$M.eigenvalues <- eigen(Approx$M, only.values = TRUE, symmetric = FALSE)$values
  }
  ## Data frame for the eigenvalues of M
  ev <- data.frame(Real = Re(Approx$M.eigenvalues), Imaginary = Im(Approx$M.eigenvalues))
  pl <- ggplot(ev, aes(Real, Imaginary)) + geom_point(alpha = 0.7) + theme_bw()
  ## Data frame for approximations
  ev.polygon <- data.frame()
  ev.extra <- data.frame()
  ## Circle and extra eigenvalue for May's stability criterion
  st <- Approx$May.stats
  radius <- sqrt((st$S - 1) * st$V)
  center <- st$d - st$E
  tmp <- eye.ellipse.df(center, radius, radius)
  tmp$Type <- "May"
  ev.polygon <- rbind(ev.polygon, tmp)
  ev.extra <- rbind(ev.extra, data.frame(Real = (st$S-1) * st$E + st$d, Imaginary = 0, Type = "May"))
  ## Ellipse and extra eigenvalue for Tang et al. criteria
  st <- Approx$TangEtAl.stats
  radius.h <- sqrt((st$S - 1) * st$V * (1 + st$rho))
  radius.v <- sqrt((st$S - 1) * st$V * (1 - st$rho))
  center <- st$d - st$E
  tmp <- eye.ellipse.df(center, radius.h, radius.v)
  tmp$Type <- "Tang et al."
  ev.polygon <- rbind(ev.polygon, tmp)
  ev.extra <- rbind(ev.extra, data.frame(Real = (st$S-1) * st$E + st$d, Imaginary = 0, Type = "Tang et al."))
  ## For cutting the graph
  bulkRe <- range(ev.polygon$Real)
  bulkIm <- range(ev.polygon$Imaginary)
  ## Circle and Ellipse for eyeball approximation
  st <- Approx$eyeball.stats
  tmp <- eye.ellipse.df(st$center.A + st$d, st$radius.A, st$radius.A)
  tmp$Type <- "eyeball, A"
  ev.polygon <- rbind(ev.polygon, tmp)
  tmp <- eye.ellipse.df(st$radius.A + st$center.A + st$d, st$radius.B.h, st$radius.B.v)
  tmp$Type <- "eyeball, shifted B"
  ev.polygon <- rbind(ev.polygon, tmp)
  ## Now draw the eigenvalues and the approximations
  pl <- pl + geom_polygon(data = ev.polygon, aes(Real, Imaginary, colour = Type), fill = NA) +
    geom_point(data = ev.extra, aes(Real, Imaginary, colour = Type), shape = 1, size = 2, show_guide = FALSE)
  ## and cut the plot
  if (only.bulk == FALSE) {
    pl <- pl + coord_cartesian(xlim = 1.05 * range(c(bulkRe, Re(Approx$M.eigenvalues))),
                              ylim = 1.05 * range(c(bulkIm, Im(Approx$M.eigenvalues))))
  } else {
    ## Cut according to the polygons
    pl <- pl + coord_cartesian(xlim = 1.05 * bulkRe,
                               ylim = 1.05 * bulkIm)
  }
  ## Legend stuff
  pl <- pl + labs(colour = "Approximation") + theme(legend.position = "bottom")
  ## colour stuff
  pl <- pl + scale_colour_brewer(palette = "Set1")
  print(pl)
  return(pl)
}
