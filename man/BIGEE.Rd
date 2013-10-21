\name{BIGEE}
\alias{BIGEE}
\title{BIGEE: Bimodal Gene Expression Exclusivity.
}
\description{
Part of the Synthetic Lethality detection in Genomics toolkit.  Detects bimodality and non-normality in all genes across the dataset.  Compares all pairwise combinations of bimodal genes and searches for mutually exclusive low expression as evidence of potential synthetic lethality.  Scores gene-pairs based on the presence of mutual exclusive bimodality and the distribution of signal intensity across the rest of the dataset.
}
\usage{
BIGEE(
	data = data, 
	confC = confC
	)
}
\arguments{
  \item{data}{This should be a log2 gene expression matrix with genes as rownames and samples as column names.  Suitable for gene expression data from any platform - NGS datasets should be RPKM or RSEM values.
}
  \item{confC}{The confidence with which you would like to call bimodality in the gene expression dataset.  Either 'high' 'med' or 'low' can be input.
}
}
\details{
The lower confidence calls will dramatically affect the number of gene pairs that the tool produces and increase the false positive rate.  The tool will take approximately 10 minutes to run a 5,000 row and 200 column input matrix using a 'medium,' confidence interval.
}
\value{
A matrix containing three columns. Columns 1 and 2 are the gene symbols that make up the candidate synthetic lethal gene pairs.  Column 3 is the score calculated the tool to rank the statistical significance of the gene pairs.
}
\author{
Mark Wappett
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(BIGEE_example_input)
outBIGEE <- BIGEE(data=BIGEE_example_input, confC='high')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SyntheticLethal }
\keyword{ Bimodality }
