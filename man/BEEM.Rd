\name{BEEM}
\alias{BEEM}
\title{BEEM: Bimodal Expression Exclusive with Mutation}
\description{
Takes a gene expression matrix and mutation matrix as input.  The mutation matrix samples (columns) must mirror or significantly overlap with the gene expression matrix.  The data in the mutation matrix must be a discreet 'WT' or 'MUT' call based on the status of each gene with each sample. Detects bimodality and non-normality in all genes across the dataset. Detects mutations of genes enriched in either the high or low gene expression modes.
}
\usage{
BEEM(
	data = data, 
	mutData = mutData, 
	confC = confC, 	
	minMut = minMut
	)
}
\arguments{
  \item{data}{This should be a log2 gene expression matrix with genes as rownames and samples as column names.  Suitable for gene expression data from any platform - NGS datasets should be RPKM or RSEM values.
}
  \item{mutData}{This should be a matrix with genes rownames and samples as column names.  All cells should be made up of a discreet 'WT' or 'MUT' call.  There should be good overlap (by sample) with the gene expression matrix
}
  \item{confC}{The confidence with which you would like to call bimodality in the gene expression dataset.  Either 'high' 'med' or 'low' can be input.
}
  \item{minMut}{The minimum number of mutations you for a gene would consider for analysis.
}
}
\details{
The lower confidence calls will dramatically affect the number of gene pairs that the tool produces and increase the false positive rate.  The tool will take approximately 10 minutes to run a 5,000 row and 200 column input gene expression matrix using a 'medium,' confidence interval.
}
\value{
A matrix containing 10 columns. Column 1 contains the bimodal genes from the expression data (gene 1) and column 2 contains the mutated candidate synthetic lethal gene pair (gene 2).  Columns 3 and 4 contain the number of mutations of gene 2 in the low and high expression modes of gene 1.  Column 5 contains the fishers p value that evaluates enrichment of mutation in either the high or low mode (indicated by column 10).  Columns 6 and 7 contain the percentage of samples in the low and high expression modes of gene 1 that are mutated for gene 2.  Columns 8 and 9 contain information on the overall size (in terms of sample) of the low and high expression modes of gene 1. 
}

\author{
Mark Wappett
}

\examples{
data(BEEM_EXPRS)
data(BEEM_MUT)
outBEEM <- BEEM(data=BEEM_EXPRS, mutData=BEEM_MUT, confC='high', minMut=10)
}

\keyword{ Mutation }
\keyword{ SyntheticLethal }
