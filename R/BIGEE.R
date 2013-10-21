# Bimodal Gene Expression Exclusivity Build V3.0
#
# Date of V1.0: November - December 2012
# Date of V2.0: April 2013
# Date of V3.0: May - June 2013	
# Author: M Wappett
# Decription: Input a matrix of continuous data with genes as rows and samples as columns.  Algorithm may take some considerable time to run (1 hour for 20,500 genes across 800 samples)

# Load relevant libraries
BIGEE <- function(data=data, confC=confC)
{
# Define confidence criteria
for(i in 1:1)
{
if(confC == "high")
{
	pI <- 0.5
	dTA <- 3.5
	bI <- 1
	numALL <- 10
}
else if(confC == "med")
{
	pI <- 0.5
	dTA <- 2.5
	bI <- 0.9
	numALL <- 5
}
else if(confC == "low")
{
	pI <- 0.55
	dTA <- 2
	bI <- 0.6
	numALL <- 1.5
}
}
print(pI)
print(dTA)
print(bI)
print(numALL)
# Deal with any negative values (< 1 min)
data2 <- apply(data, 1, function(x) { 
			w1 <-  which(x[1:dim(data)[2]] < 0)
			len1 <- length(w1)
			vals1 <- runif(len1)
			x[w1] <- vals1
			x
})
data2 <- t(data2)

# Work out bimodality indexes of genes, and then rank based on this index (< 10 mins)
bimodalIndex <- function(dataset, verbose = TRUE) 
{
	# From Coombes,KR (2012) ClassDiscovery: Classes and methods for "class discovery" with microarrays or proteomics. 
    bim <- matrix(NA, nrow = nrow(dataset), ncol = 6)
    if (verbose) 
        cat("1 ")
    for (i in 1:nrow(dataset)) {
        if (verbose && 0 == i%%100) 
            cat(".")
        if (verbose && 0 == i%%1000) 
            cat(paste("\n", 1 + i/1000, " ", sep = ""))
        x <- as.vector(as.matrix(dataset[i, ]))
        if (any(is.na(x))) 
            next
        mc <- Mclust(x, G = 2, modelNames = "E")
        sigma <- sqrt(mc$parameters$variance$sigmasq)
        delta <- abs(diff(mc$parameters$mean))/sigma
        pi <- mc$parameters$pro[1]
        bi <- delta * sqrt(pi * (1 - pi))
        bim[i, ] <- c(mc$parameters$mean, sigma = sigma, delta = delta, 
            pi = pi, bim = bi)
    }
    if (verbose) 
        cat("\n")
    dimnames(bim) <- list(rownames(dataset), c("mu1", "mu2", 
        "sigma", "delta", "pi", "BI"))
    bim <- as.data.frame(bim)
    bim
}
biIndex <- bimodalIndex(data2)
subBiIndex <- subset(biIndex, biIndex[,6] > bI & biIndex[,5] < pI & biIndex[,4] > dTA)
subBiIndex2 <- subBiIndex[order(-subBiIndex[,6]),]

# Set up bimodality in gene expression (BIG) function and run
indexing<-function(x.sort,test.bnd,alpha)
{
n<-length(x.sort);
low<-x.sort[1:test.bnd];
hgh<-x.sort[(test.bnd+1):n];
n1<-length(low);
n2<-length(hgh);

qr.low<-quantile(low,seq(0,1,0.01));
qr.hgh<-quantile(hgh,seq(0,1,0.01));

v1<-abs(qr.low[[75]][1]-qr.low[[25]][1])/1.34896;# statistics of low cluster
v2<-abs(qr.hgh[[75]][1]-qr.hgh[[25]][1])/1.34896;# statistics of high cluster

gap.bet.2.cls<-abs(max(low)-min(hgh));
dst.bet.2.cls<-abs(qr.low[[75]][1]-qr.hgh[[25]][1]);
het.std<-sqrt((n1*v1^2+n2*v2^2)*(1/n1+1/n2)/n);
rtv<-(alpha*gap.bet.2.cls+(1-alpha)*dst.bet.2.cls)/het.std;

return(rtv);
}

big<-function(x,alpha,top)
#get index for one gene
#x: expression vector of a gene
#alpha: weight on gap between two clusters
#top: stop rule for scanning
{
x.sort<-sort(x);
x.sort.or<-order(diff(x.sort),decreasing=T);#order of gaps

max.index<-0;#default index
my.boundary<-0;#default boundary

prediction<-c();
for(i in 1:top)
	{
	new.ind<-indexing(x.sort,x.sort.or[i],alpha);
	new.bnd<-mean(c(x.sort[x.sort.or[i]],x.sort[x.sort.or[i]+1]));
	prediction<-rbind(prediction,c(new.ind,new.bnd));
	}

sel<-which(prediction[,1]==max(prediction[,1]))[1];#find the highest BIG index

return(list(prediction,c(prediction[sel,1],prediction[sel,2])));
}

Besag.Clifford<-function(BI,H,B)#Beasg-Clifford algorithm for p value
#BI: bimodal index list
#H: rejection number (100: reasonable)
#B: random sampling times (1e5: reasonable)
{
p<-c();
for(i in 1:length(BI))
	{
	bs<-c();
	h<-0;
	while(h<H & length(bs)<B)
		{
		bs<-c(bs,BI[sample(1:length(BI),1)]);
		h<-sum(1*(bs>=BI[i]));
		}
	p<-c(p,h/length(bs));
	}
return(p);
}

BIG<-function(dat)
#the subroutine to use BIG algorithm
#dat is a list
#X: expression matrix (logarithm or normalisation pre-processed)
#H: rejection number used by Beasg.Clifford algorithm
#B: random sampling times used by Besag.Clifford algorithm
#alpha: weight on gap between two clusters
#top: stop rule when scanning a gap vector
{
if(length(dat$X)==0)
	{
	print("no expression matrix!");
	return();
} else	{
	X<-dat$X;
	H<-100;
	if(length(dat$H)>0)
		H<-dat$H;
	B<-10000;
	if(length(dat$B)>0)
		B<-dat$B;
	alpha<-0.75;
	if(length(dat$alpha)>0)
		alpha<-dat$alpha;
	top<-10;
	if(length(dat$top)>0)
		top<-dat$top;
	BI<-c();
	all<-c();
	for(i in 1:nrow(X))
		{
		buf<-big(X[i,],alpha,top);
		all<-c(all,buf);
		BI<-rbind(BI,buf[[2]]);
		}
	p<-Besag.Clifford(BI[,1],H,B);
	model<-list(BI=BI[,1],BND=BI[,2],p=p);
	return(model);
	}
}

data<-function(fn)
{
X<-as.matrix(read.csv(fn,header=FALSE));
return(X);
}

# Run the big method (will take ~2 mins per 1,000 samples)
big.model<-BIG(list(X=data2))

# Calculate the mid-points for all genes identified as being bi-modal
midpoint <- big.model[[2]]
names(midpoint) <- rownames(data2)
w1 <- rownames(subBiIndex)
midpointBI <- midpoint[w1]

data3 <- data2[w1,]
rownames(data3) <- w1

med1 <- apply(data3, 1, function(x) mean(as.numeric(x[1:dim(data3)[2]])))
w1 <- which(med1 > 2)
data3 <- data3[w1,]
w1 <- rownames(data3)
midpointBI <- midpointBI[w1]

# Set data up for pure synthetic lethal detection
list3 <- as.list(as.data.frame(t(data3)))
outList <- list(data.frame=(0))
ind1 <- 0
names1 <- names(list3)

# Run pure synthetic lethal detection
subBiIndex3 <- subBiIndex[rownames(data3),]
rangeTab <- rbind(c(0.55, 0.5), c(0.5, 0), c(0.45, 0.5), c(0.4, 1), c(0.35, 1.5), c(0.30, 2), c(0.25, 2.5), c(0.20, 3), c(0.15, 3.5), c(0.10, 4), c(0.05, 4.5), c(0.00, 5))
for(i in 1:(length(list3)))
{
	range1 <- range(list3[[i]])
	range1_2 <- ((range1[2] - range1[1])/100)* (rangeTab[which.min(abs(rangeTab[,1] -subBiIndex3[i,5])),2])
	midpoint1 <- midpointBI[i]
	for(j in 1:(length(list3)))
	{
		range2 <- range(list3[[j]])
		range2_2 <- ((range2[2] - range2[1])/100)* (rangeTab[which.min(abs(rangeTab[,1] -subBiIndex3[j,5])),2])
		midpoint2 <- midpointBI[j]
		midpoint3 <- midpoint1 - range1_2
		midpoint4 <- midpoint2 - range2_2
		num1 <- which(list3[[i]] < midpoint3)
		pair2 <- list3[[j]][num1]
		num2 <- which(pair2 < midpoint4)
		num3 <- which(list3[[j]] < midpoint4)
		num4 <- (sum(length(num1), length(num3))/(length(list3[[1]])*2)) * 100
		len1 <- c(length(num1), length(num3))
		perLL <- length(num2)/length(list3[[i]])*100
		sd1 <- sd(len1)
		rank1 <- (num4 / sd1)
		if(perLL < 2 & num4 > numALL & num4 < 95)
		{
			ind1 <- ind1 + 1
			BBMM <- c(names1[i], names1[j])
			outList[[ind1]] <- BBMM
		}
		else
		{
			# Do Nothing
		}		
	}
	print(i)
}
x2 <- do.call("rbind", outList)
x3 <- as.data.frame(x2[,1:2], stringsAsFactors=FALSE)
x4 <- x3[x3[,1] <= x3[,2],]
colnames(x4) <- c("gene", "gene2")

biIndex$gene <- rownames(biIndex)
colnames(biIndex)[7] <- "gene"
x5 <- merge(x4, biIndex, by.x="gene", by.y="gene")
colnames(biIndex)[7] <- "gene2"
x6 <- merge(x5, biIndex, by.x="gene2", by.y="gene2")
x7 <- x6[,1:2]
score1 <- ((x6[,8] + x6[,14]) * (x6[,7] + x6[,13]))* (x6[,6] + x6[,12])
x7$score <- score1
x7 <- x7[order(-x7[,3]),]
return(x7)
}
