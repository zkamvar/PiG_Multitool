#--------------------------------------------------------------------------#
# This function will create a list of csv files that contain the isolates in 
# rows and loci in columns. This currently only works for unlabeled files,
# but what it will do is allow the user to specify an entire directory to 
# analyze. 
#
# It will change in the future, but it works for my purposes for now.
#--------------------------------------------------------------------------#
getfile <- function(multFile=NULL, pattern=NULL){
# the default option is to grab all of the files in a directory matching a 
# specified pattern. If no pattern is set, all files will be listed.
	if (is.null(multFile) || multFile == "yes"){
		# this sets the path variable that the user can use to set the path
		# to the files with setwd(x$path), where x is the datastructure 
		# this function dumped into.
		path <- sub("(^.+?/).+?\\.[a-z]{3}", "\\1", file.path(file.choose()))		
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
		path <- sub("(^.+?/).+?\\.[a-z]{3}", "\\1", file.path(csv))
		csv <- sub("^.+?/(.+?\\.[a-z]{3})", "\\1", csv)
		x <- csv
	}
	filepath <- list(files=x, path=path)
	return(filepath)
}

#--------------------------------------------------------------------------#
# Here is the function that will produce the descrete distance matrix for 
# diploid organisms.
# pop.matrix is a matrix of rows of individuals with columns of Loci.
# Currently the loci are every two columns in the matrix.
# Specifically, this compares two individuals at one locus. A loop for each 
# type of distance matrix will need to be constructued separately
# 
# i is the row for each new individual in the pairwise comparison
# j is the row for the reference individual
# m is the column for the first allele in the locus
# n is the column for the second allele in the locus
# z is the measure of distance. It should be zero before the loop
#
# Again, this is for diploids. position of the allele at the locus does not
# matter as the inheritance is unknown, so this algorithm takes that into
# account when doing the calculations.
#--------------------------------------------------------------------------#	
difference.test <- function(pop.matrix, i, j, m, n, z){
	# This if loop is analyzing allele m of individual j against that of
	# individual i. If they are equal, then that means the distance is equal
	# to zero at that allele.
	if(pop.matrix[j,m] == pop.matrix[i,m]){
		# given that the distance between allele m in individuals j and i
		# are equal, there is only one more comparison to do. If allele n
		# of individual j is equal to that of individual i, then the
		# distance is equal to 0. If they are not equal, z at the locus is 1
		if(pop.matrix[j,n] != pop.matrix[i,n]){
			z <- z+1
		}
	}
	# If allele m of individual j is not equal to that of individual i, then
	# a test to see if allele m is equal to allele n in individuals j and i,
	# respectively.
	else if(pop.matrix[j,m] == pop.matrix[i,n]){
		# given that the distance between alleles m and n in individuals j
		# and i respectively are equal, the only comparison left is to 
		# compare the distances between alleles m and n in individuals i and
		# j respectively. If they are not equal, z at the locus is 1
		if(pop.matrix[j,n] != pop.matrix[i,m]){
			z <- z+1
		}
	}
	# If neither allele in individual i is equal to allele m in individual j
	# z is equal to one and a tests if the n allele of individual j is equal
	# to either allele in individual i.
	else{
		z <- z+1
		# Testing if the n allele in individual j is not equal to either allele
		# in individual i. If neither is equal, z at the locus increases by 1.
		if(pop.matrix[j,n] != pop.matrix[i,m] && pop.matrix[j,n] != pop.matrix[i,n]){
			z <- z+1
		}	
	}
	# The final value of z is returned. At this point it can take on the
	# values of 0, 1, or 2 if doing a pairwise comparison at a single locus.
	# For pairwise comparisons over multiple loci, z can take on any value
	# between 0 and 2*M where M is the number of loci sampled.
	return(z)
}

stdIa <- function(x){
	#--------------------------------------------------------------------------#
	# This is the beginning of the analysis. It will take the list of files 
	# provided by the user and pull the files out of the directory that has 
	# been set by the user. 
	#--------------------------------------------------------------------------#
	for(a in 1:length(x)){
		pop <- read.table(file(description= x[a], open="r")) 
		pop.matrix <- as.matrix(pop)
		numAlleles <- ncol(pop.matrix)
		numIsolates <- nrow(pop.matrix)
		np <- (numIsolates*(numIsolates-1))/2
		# Creation of the datastructures
		if (!exists("Ia.vector")){
			Ia.vector <- NULL
			rbarD.vector <- NULL
			file.vector <- NULL
		}
		D.matrix <- matrix(0, nrow=numIsolates, ncol=numIsolates)
		d.matrix <- D.matrix
		d.vector <- NULL
		d2.vector <- NULL
		vard.vector <- NULL
		vardpair.vector <- NULL
		to.remove <- NULL
	#--------------------------------------------------------------------------#
	# This is the loop for analyzing the pairwise distances at each locus,
	# also known as d. These values will be placed into two vectors representing
	# the sum of d for each locus and the sum of d^2 for each locus. 
	# 
	# This will also calculate D, the pairwise comparison of all isolates over
	# all loci. 
	#--------------------------------------------------------------------------#
		# Initiating the loop over the columns of the population matrix
		for(m in seq(1, (numAlleles-1), 2)){
			n <- m+1
			# Loop for the columns of the distance matrix
			for(j in 1:numIsolates){
				# Loop for the rows of the distance matrix
				for(i in 1:numIsolates){
					z <- 0
					# Setting the contraint for the pairwise comparisons by the
					# equation (n(n-1))/2 where n = numIsolates 
					# If this loop did not exist, the construction of the
					# distance matrix would take twice as long.
					if(j < i && i <= numIsolates){
						z <- difference.test(pop.matrix, i, j, m, n, z)
					}
					# The value of z (0, 1, or 2) is pushed into the matrix
					d.matrix[i,j] <- d.matrix[i,j]+z
					# This value is added to the matrix for pairwise comparison
					# of all the isolates over all loci as opposed to each locus
					D.matrix[i,j] <- D.matrix[i,j]+d.matrix[i,j]
				}
			}
			# placing the sum of the resulting matrix into a vector for later
			# use
			d.vector <- append(d.vector, sum(d.matrix))
			# placing the sum of the squares of the resulting matrix into a 
			# vector for later use
			d2.vector <- append(d2.vector, sum(d.matrix^2))
			# zeroing out the matrix
			d.matrix[1:numIsolates,1:numIsolates] <- 0
		}
		# removing the matrix from the namespace as it is no longer needed.
		rm(d.matrix)
	#--------------------------------------------------------------------------#
	# Now to begin the calculations. First, set the variance of D
	#--------------------------------------------------------------------------#
		varD <- ((sum(D.matrix^2)-((sum(D.matrix))^2)/np))/np
	#--------------------------------------------------------------------------#
	# Next is to create a vector containing all of the variances of d (there
	# will be one for each locus)
	#--------------------------------------------------------------------------#
		vard.vector <- ((d2.vector-((d.vector^2)/np))/np)
	#--------------------------------------------------------------------------#
	# Here the roots of the products of the variances are being produced and
	# the sum of those values is taken.
	#--------------------------------------------------------------------------#
		for (b in 1:length(d.vector)){
			for (d in 1:length(d.vector)){
				# As pairwise multiplication is required, a pairwise constriant
				# loop must be set up. 
				if(b < d && d <= length(d.vector)){
					vardpair <- sqrt(vard.vector[b]*vard.vector[d])
					vardpair.vector <- append(vardpair.vector, vardpair)
				}
			}
		}
	#--------------------------------------------------------------------------#
	# The sum of the variances necessary for the calculation of Ia is calculated
	#--------------------------------------------------------------------------#
		sigVarj <- sum(vard.vector)
		rm(vard.vector)
	#--------------------------------------------------------------------------#
	# Finally, the Index of Association and the standardized Index of associati-
	# on are calculated.
	#--------------------------------------------------------------------------#
		Ia <- (varD/sigVarj)-1
		rbarD <- (varD - sigVarj)/(2*sum(vardpair.vector))
		# Prints to screen as loop progresses
		print(paste("File Name:", x[a], sep=" "))
		print(paste("Index of Association:", Ia, sep=" "))
		print(paste("Standardized Index of Association (rbarD):", rbarD, sep=" "))
		# Saves the values of Ia, rbarD, and the filename into datastructures
		# that will be listed in a dataframe
		file.vector <- append(file.vector, x[a])
		Ia.vector <- append(Ia.vector, Ia)
		rbarD.vector <- append(rbarD.vector, rbarD)
	}
	# Creating the data frame for output.
	Iout <- list(Ia=Ia.vector, rbarD=rbarD.vector, File=file.vector)
	return(as.data.frame(Iout))
}