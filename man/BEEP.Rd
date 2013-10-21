\name{BEEP}
\alias{BEEP}
\title{BEEP: Bimodal Expression Exclusive with Pharmacology}
\description{Detects bimodality and non-normality in all genes across the dataset.  Integrates second matrix with discreet resistance/ sensitivity calls for any given compound for samples present in the gene expression matrix.  Indicates if there is an enrichment of resistance in either the high or low expression modes of any of the genes
}
\usage{
BEEP(
	data = data, 
	sensData = sensData, 
	confC = confC
	)
}
\arguments{
  \item{data}{This should be a log2 gene expression matrix with genes as rownames and samples as column names.  Suitable for gene expression data from any platform - NGS datasets should be RPKM or RSEM values.
}
  \item{sensData}{This should be a 2 column matrix where column 1 contains the sample names from the gene expression matrix and column 2 gives a set of discreet values describing the sensitivity of the sample to compound.  A number of these values need to be called 'Resistant' to run the tool
}
  \item{confC}{The confidence with which you would like to call bimodality in the gene expression dataset.  Either 'high' 'med' or 'low' can be input.
}
}
\details{
The lower confidence calls will dramatically affect the number of gene pairs that the tool produces and increase the false positive rate.  The tool will take approximately 10 minutes to run a 5,000 row and 200 column input gene expression matrix using a 'medium,' confidence interval.
}
\value{
A list of 2 matrices.
%%  \item{highEnriched}{Contains genes where resistance is enriched in the high expression mode of the gene}
%%  \item{lowEnriched}{Contains genes where resistance is enriched in the low expression mode of the gene}
}
\author{
Mark Wappett
}

\examples{
data(BEEP_input_data)
data(BEEP_sens_data)
outputBEEP <- BEEP(data=BEEP_input_data, sensData=BEEP_sens_data, confC='high')
}
\keyword{ Bimodality }
\keyword{ Pharmacology }
