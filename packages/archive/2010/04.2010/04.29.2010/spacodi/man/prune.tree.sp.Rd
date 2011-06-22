\name{prune.tree.sp}

\alias{prune.tree.sp}

\title{pruning a tree to match datasets used in SPACoDi}

\description{
 \code{prune.tree.sp} is a utility for automating the process of pruning (especially large) trees for use in \pkg{spacodi}. Although seldom required, users
 also have the option to randomly resolve a tree with multifurcations (\code{resolve=TRUE}. Note that nearly all \pkg{spacodi} functions require that
 trees have present (at least) all species that are contained within the \code{spacodi.object}.}

\usage{prune.tree.sp(phy, spacodi.object, resolve=FALSE)}

\arguments{
 \item{phy}{a phylogenetic tree with (at least some) tips matching row names of the \code{spacodi.object}}
 \item{spacodi.object}{a dataset (community data or traits, for instance), where row names match tips in \code{phy}}
 \item{resolve}{a logical indicator determining whether multifurcations are resolved with zero-length branches; this is not required for many of the functions in the remainder of this package}
}

\value{A pruned phylogeny}

\author{Jonathan Eastman}

\seealso{\code{\link[ape]{multi2di}} for random tree-resolution}
