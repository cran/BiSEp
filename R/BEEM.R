# Mutation and Gene Expression Automated Detected Algorithm Build V2.0
#
# Date of V1.0: February - March 2013
# Date of V2.0: May - June 2013	
# Author: M Wappett
# Decription: Input a matrix of continuous data with genes as rows and samples as columns.  Input a matrix of discreet mutation data ("WT" or "MUT" calls) with genes as rows and samples as columns.  Algorithm may take some considerable time to run (1 hour for 20,500 genes across 800 samples).

BEEM <- function(data=data, mutData=mutData, confC=confC, minMut=minMut)
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

# Deal with gene expression matrix negative values (< 1 min)
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

midpoint <- big.model[[2]]
names(midpoint) <- rownames(data2)
w1 <- rownames(subBiIndex)
midpointBI <- midpoint[w1]
w1 <- which(midpointBI > 4)
midpointBI <- midpointBI[w1]
w1 <- names(midpointBI)

data3 <- data2[w1,]
rownames(data3) <- w1

# Perform analysis of mutation data, and then integrate with the gene expression analysis outcomes
c1 <- apply(mutData, 1, function(x) length( x [x == "MUT"] ))
mutData$mutCount <- c1
mutData2 <- subset(mutData, mutData$mutCount >= minMut)


# Evaluate the significance of this
w1 <- which(colnames(mutData2) %in% colnames(data3))
mutData3 <- mutData2[,w1]
w1 <- which(colnames(data3) %in% colnames(mutData2))
data4 <- data3[,w1]

# Order mutation data columns the same as expression data columns.
mutData3 <- mutData3[,colnames(data4)]

list1 <- as.list(as.data.frame(t(data4)))
names1 <- colnames(data4)
names2 <- colnames(mutData3)
outList <- list(data.frame(0))
ind1 <- 0
for(i in 1:length(list1))
{
	data5 <- list1[[i]][order(list1[[i]])]
	names1_1 <- names1[order(list1[[i]])]
	mutData4 <- mutData3[,names1_1]
	mpNum <- midpointBI[[i]] - (5*((max(data5) - min(data5))/100))
	sampleMidIndex <- which.min(abs(data5 - mpNum))
	lower <- data5[1:sampleMidIndex]
	upper <- data5[(sampleMidIndex+1): length(data5)]
	contTab <- matrix(ncol=2, nrow=2)
	valVec <- numeric(0)
	list2 <- as.list(as.data.frame(t(mutData4)))
	for(j in 1:length(list2))
	{
		lowerM <- list2[[j]][1:sampleMidIndex]
		upperM <- list2[[j]][(sampleMidIndex+1):length(list2[[j]])]
		count1 <- length(which(lowerM == "MUT"))
		count2 <- length(which(upperM == "MUT"))
		count3 <- count1 + count2
		contTab[1,1] <- count1
		contTab[1,2] <- count2
		contTab[2,1] <- length(lower) - count1
		contTab[2,2] <- length(upper) - count2
		per1 <- (count1/length(lower))*100
		per2 <- (count2/length(upper))*100
		per3 <- per1 - per2
		pV <- fisher.test(contTab)
		ind1 <- ind1 + 1
		if(per3 > 0)
		{
			outList[[ind1]] <- c(names(list1)[i], rownames(mutData4)[j], count1, count2, pV[[1]], per1, per2, length(lowerM), length(upperM), "LowEnriched")
		}
		else
		{
			outList[[ind1]] <- c(names(list1)[i], rownames(mutData4)[j], count1, count2, pV[[1]], per1, per2, length(lowerM), length(upperM), "HighEnriched")
		}
	}
	print(i)
}
genePairs2 <- do.call("rbind", outList)
genePairs2 <- as.data.frame(genePairs2, stringsAsFactors=FALSE)
genePairs2 <- genePairs2[order(genePairs2[,5]),]
genePairs3 <- subset(genePairs2, genePairs2[,10] == "HighEnriched" & genePairs2[,5] < 0.05)
genePairs3 <- as.data.frame(genePairs3, stringsAsFactors=FALSE)
genePairs3[,5] <- as.numeric(genePairs3[,5])
genePairs4 <- genePairs3[order(genePairs3[,5]),]
colnames(genePairs4) <- c("Gene1", "Gene2", "LowerExpressionMutationCount", "HighExpressionMutationCount", "Fishers P Value", "Percentage of lower samples mutated", "Percenage of high samples mutated", "Size of low expression population", "Size of high expression population", "Enrichment Status")
return(genePairs4)
}
