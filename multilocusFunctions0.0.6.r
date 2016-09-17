#==============================================================================#
# This function will place any list of files with any pattern into a vector
# in an object. The way it should be used for the functions on this page
# is to search for .dat files. 
#==============================================================================#
getfile <- function(multFile=TRUE, pattern=NULL){
# the default option is to grab all of the files in a directory matching a 
# specified pattern. If no pattern is set, all files will be listed.
	if (multFile==TRUE){
		# this sets the path variable that the user can use to set the path
		# to the files with setwd(x$path), where x is the datastructure 
		# this function dumped into.
		path <- sub("(^.+?/).+?\\.[a-z]{3}$", "\\1", file.path(file.choose()))		
		if (!is.null(pattern)){
			pat <- pattern
			x <- list.files(path, pattern=pat)
		}
		else {
			x <- list.files(path)
		}
	}
	else {
		# if the user chooses to analyze only one file, a pattern is not needed
		csv <- file.choose()
		path <- sub("(^.+?/).+?\\.[a-z]{3}$", "\\1", file.path(csv))
		csv <- sub("^.+?/(.+?\\.[a-z]{3}$)", "\\1", csv)
		x <- csv
	}
	filepath <- list(files=x, path=path)
	return(filepath)
}
#==============================================================================#
# The calculation of the Index of Association and standardized Index of 
# Association.
#==============================================================================#
stdIa <- function(pop, subpops=FALSE, sample=0){	
	# Testing for vaild filetypes and saving the filename.	
	x <- .file.type(pop)
	# converting the population into genind class if it already isn't.	
	pop <- x$POP	
	x <- x$X	
	# Separating out all of the subpopulations and creating an indicator for
	# multiple subpopulations.
	if (subpops==TRUE){
		pops <- .pop.divide(pop)$POP
		MPI <- .pop.divide(pop)$MPI
	}
	else {
		MPI <- NULL
	}
	# calculating indecies for total population	
	Iout <- as.data.frame(.Ia.Rd(pop, x))
	# calculating incecies for subpopulations
	if (!is.null(MPI)){
		for (r in pops){
			Iout <- rbind(Iout, .Ia.Rd(r, x))
		}
	}	
	if (sample != 0){
		Iout <- rbind(Iout, .sampling(pop, x, sample))
		if (!is.null(MPI)){
			for (r in pops){
				Iout <- rbind(Iout, .sampling(r, x, sample))
			}
		}		
	}
	return(Iout)
}
#==============================================================================#
# This function allows the user to analyze a list of files automatically
#==============================================================================#
stdIa.all <- function(x, subpops=FALSE, sample=0) {
	Iout <- NULL
	for(a in x){
		Iout <- rbind(Iout, stdIa(a, subpops=subpops, sample=sample))
	}
	Iout <- extract.info(Iout)
	return(Iout)
}
#==============================================================================#
# This function will extract all relevant information from the files
#==============================================================================#
extract.info <- function(x) {
	if (length(grep("^clone.+?dat$", x$File)) != 0){
		# Rate of clonal reproduction
		x$Clone <- as.numeric(sub("^clone.(\\d{3}.\\d{2}).+?dat$","\\1", x$File))
		# Rate of Sexual reproduction
		x$Sex.Rate <- (100-x$Clone)/100
	}
	if (length(grep(".+?rep.+?dat$", x$File)) != 0){
		# Replicate indicators
		x$Replicate <- sub(".+?rep.(\\d{2}).+?dat$", "\\1", x$File)
	}
	if (length(grep(".+?pop.+?dat$", x$File)) != 0){
		# Population size indicators
		x$Pop.Size <- sub(".+?pop.(\\d+?).+?dat$","\\1", x$File)
	}
	if (length(grep(".+?sam.+?dat$", x$File)) != 0){
		# Sample size
		x$Samp.Size <- sub(".+?sam.+?(\\d{2,3}).+?dat$", "\\1", x$File)
	}
	return(x)
}
#==============================================================================#
# .file.type will make sure that the data structure entering into the stdIa
# function is a genind object
#==============================================================================#
.file.type <- function(pop){
  if (!is.genind(pop)){
    x <- pop
    pop <- import2genind(pop)
  }
  else { 
    x <- as.character(pop@call)[2]		
  }
  return(list(X=x, POP=pop))
}
#==============================================================================#
# .pop.divide will attempt to deal with separate populations
#==============================================================================#
.pop.divide <- function(x) {
	if (!is.null(x@pop)) {
		pop <- seppop(x)
		mult.pop.ind <- 1
	}
	else {
		mult.pop.ind <- NULL
	}
	return(list(MPI=mult.pop.ind, POP=pop))
}
#==============================================================================#
# .pairwise.differences will calculate three vectors that will be used for the
# calculation of the Index of Association and standardized Index of Association
# Later.
# pop = genind object 
# count = running count of the number of loci
# numIsolates = number of isolates in the sample
# temp.d.vector = temporary vector to store the differences
# d.vector = a vector of the sum of the differences at each locus. The length
# 			 of this vector will be the same as the number of loci.
# d2.vector = the same as d.vector, except it's the sum of the squares
# D.vector = a vector of the the pairwise distances over all loci. The length
#			 of this vector will be the same as n(n-1)/2, where n is number of
# 			isolates.
#==============================================================================#
.pairwise.differences <- function(pop,numIsolates){	
	j <- 0
	temp.d.vector <- NULL
	d.vector <- NULL
	d2.vector <- NULL
	D.vector <- NULL
	for (loc in names(pop@loc.nall)){
		# loc is the list of all loci in the population
		n <- j + 1
		# n is the starting point for the range of alleles in each locus
		j <- pop@loc.nall[loc] + j
		# j is the final point for the range of alleles in each locus
		for(m in 1:numIsolates){
			# Loop for the rows of the distance matrix
			for(i in 1:numIsolates){
				z <- 0
				# Setting the contraint for the pairwise comparisons by the
				# equation (n(n-1))/2 where n = numIsolates 
				if(m < i && i <= numIsolates){
					z <- sum(abs(pop@tab[m,n:j] - pop@tab[i,n:j]))
				# The value of z (0, 1, or 2) is pushed into the temp vector
					temp.d.vector <- append(temp.d.vector, z)
				}
			}
		}
		d.vector <- append(d.vector, sum(temp.d.vector))
		d2.vector <- append(d2.vector, sum(temp.d.vector^2))
		# the values are added onto D.vector for pairwise comparison of all
		# of the isolates over all loci as opposed to each locus.
		if (is.null(D.vector)){
			D.vector <- temp.d.vector
		} 
		else {
			D.vector <- D.vector + temp.d.vector
		}
	temp.d.vector <- NULL
	}
	vectors <- list(d.vector=d.vector, d2.vector=d2.vector, D.vector=D.vector)
	return(vectors)
}
#==============================================================================#
# To calculate rbarD, the pairwise variances for each locus needs to be
# caluclated. 
#==============================================================================#
.pairwise.variances <- function(vard.vector, V){	
	# Here the roots of the products of the variances are being produced and
	# the sum of those values is taken.
	vardpair.vector <- NULL
	for (b in 1:length(V$d.vector)){
		for (d in 1:length(V$d.vector)){
			if(b < d && d <= length(V$d.vector)){
				vardpair <- sqrt(vard.vector[b]*vard.vector[d])
				vardpair.vector <- append(vardpair.vector, vardpair)
			}
		}
	}
	return(vardpair.vector)
}
#==============================================================================#
# The actual calculation of Ia and rbarD. This allows for multiple populations
# to be calculated.
#==============================================================================#
.Ia.Rd <- function(pop, x, sample=0){
	if (!exists("Ia.vector")){
		Ia.vector <- NULL
		rbarD.vector <- NULL
		file.vector <- NULL
		population.vector <- NULL
		sample.vector <- NULL
	}
	vard.vector <- NULL
	numAlleles <- length(pop@loc.names)
	numIsolates <- length(pop@ind.names)
	# Creating this number is necessary because it is how the variance is
	# calculated.
	np <- (numIsolates*(numIsolates-1))/2
	V <- .pairwise.differences(pop,numIsolates)
	# Now to begin the calculations. First, set the variance of D	
	varD <- ((sum(V$D.vector^2)-((sum(V$D.vector))^2)/np))/np
	# Next is to create a vector containing all of the variances of d (there
	# will be one for each locus)
	vard.vector <- ((V$d2.vector-((V$d.vector^2)/np))/np)
	vardpair.vector <- .pairwise.variances(vard.vector, V)
	# The sum of the variances necessary for the calculation of Ia is calculated
	sigVarj <- sum(vard.vector)
	rm(vard.vector)
	# Finally, the Index of Association and the standardized Index of associati-
	# on are calculated.
	Ia <- (varD/sigVarj)-1
	rbarD <- (varD - sigVarj)/(2*sum(vardpair.vector))
	# Prints to screen as loop progresses
	print(paste("File Name:", x, sep=" "))
	print(paste("Index of Association:", Ia, sep=" "))
	print(paste("Standardized Index of Association (rbarD):", rbarD, sep=" "))
	if (length(pop@pop.names) == 1){
		print(paste("Population:", pop@pop.names, sep=" "))
		population.vector <- append(population.vector, pop@pop.names)
	}
	else {
		print(paste("Population: Total"))
		population.vector <- append(population.vector, "Total")
	}
	if (sample != 0){
		print(paste("Sub Sample:", sample, sep=" "))
		sample.vector <- append(sample.vector, "Yes")
	}
	else {
		print(paste("Sub Sample: Original"))
		sample.vector <- append(sample.vector, "Original")
	}
	# Saves the values of Ia, rbarD, and the filename into datastructures
	# that will be listed in a dataframe
	file.vector <- append(file.vector, x)
	Ia.vector <- append(Ia.vector, Ia)
	rbarD.vector <- append(rbarD.vector, rbarD)
	return(as.data.frame(list(Ia=Ia.vector, rbarD=rbarD.vector, 
			File=file.vector, Population=population.vector, Sub.Sample=sample.vector)))
}
################################################################################
#==============================================================================#
# Testing version of stdIa calculation
#==============================================================================#
################################################################################
#==============================================================================#
# .sampling will reshuffle the alleles per individual, per locus. 
# Future plans will incorporate sampling that will shuffle based on populations
# as some of the populations will not have certain alleles. 
#==============================================================================#
.sampling <- function(pop, filename, iterations){
	# creating a copy of the object	
	pop.new <- pop
	sample.data <- NULL
	for (c in 1:iterations){
		j <- 0
		# making a loop for each locus
		for (loc in names(pop@loc.nall)){
			# loc is the list of all loci in the population
			n <- j + 1
			# n is the starting point for the range of alleles in each locus
			j <- pop@loc.nall[loc] + j
			# j is the final point for the range of alleles in each locus
			plus <- n:j
			# a vector for all possible alleles in the locus
			minus <- NULL
			# a vector for the alleles to remove from sampling
			for (i in n:j){
				if (sum(pop@tab[,i])==0){
					# if the sum of all the individuals at a specified allele
					# is zero, then it is not present in the population, and 
					# should be removed as a potential resampling target.
					minus <- append(minus, -i)
				}
				else {
					minus <- append(minus, 0)
				}
			} 
			for(m in 1:length(pop@tab[,1])){
				# eg. plus= 1,2,3,4; minus= -1,0,-3,0; plus+minus=0,2,0,4
				# this means that only the 2 and 4 alleles will be subsampled
				pop.new@tab[m,select=c(plus+minus)] <- sample(pop@tab[m,select=c(plus+minus)])
			}	
		}
		# calculating Ia and rbarD for the resampled data
		sample.data <- rbind(sample.data, as.data.frame(.Ia.Rd(pop.new, filename, sample=c)))
	}	
	return(sample.data)
}
histpdf <- function(data){
    for (i in 1:length(summary(data$Population))-1){
        popn <- subset(data, Population=="i" & Sub.Sample=="Yes")
        #pdf(file=sprintf("%s_%d.pdf",data$File[1],i))
        hist(popn$rbarD, breaks=c(-data$rbarD[1],data$rbard[1]))
        abline(v=data$rbarD[i+1])
        #dev.off()
    }
}
