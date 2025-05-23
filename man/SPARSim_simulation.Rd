% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPARSim.R
\name{SPARSim_simulation}
\alias{SPARSim_simulation}
\title{Simulate ST Data Using SPARSim}
\usage{
SPARSim_simulation(
  estimated_result,
  de_prop = 0.2,
  fc_group = 2,
  nBatches = 1,
  prop_batch = 1,
  verbose = FALSE,
  seed,
  return_format,
  ...
)
}
\arguments{
\item{estimated_result}{The SCST object after estimating parameters using [SPARSim::SPARSim_estimate_parameter_from_data()]}

\item{de_prop}{The proportion of DEGs over all genes between different simulated cell groups. SPARSim can only simulated cell groups if users provide group labels of reference cells to estimate parameters. Otherwise, the cell groups can not be simulated in this step and `de_prop` parameter is useless. Default is 0.2.}

\item{fc_group}{The fold change of the generated DEGs. This parameter works when the reference cell groups are provided in the estimation step. Default is 2.}

\item{nBatches}{The number of cell batches to be simulated. Default is 1.}

\item{prop_batch}{The proportion of cells in each batch. Default is 1.}

\item{verbose}{Whether to return messages or not}

\item{seed}{Random seed}

\item{return_format}{The format of returned simulation data. Choices: list, SingleCellExperiment and Seurat.}

\item{...}{Other parameters represented in SPARSim, see [SPARSim::srtsim_count()]}
}
\description{
Simulate ST Data Using SPARSim
}
\references{
Baruzzo G, Patuzzi I, Di Camillo B. SPARSim single cell: a count data simulator for scRNA-seq data. Bioinformatics, 2020, 36(5): 1468-1475. <https://doi.org/10.1093/bioinformatics/btz752>
}
