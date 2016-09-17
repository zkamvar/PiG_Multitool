#==============================================================================#
# This function will place any list of files with any pattern into a vector
# in an object. The way it should be used for the functions on this page
# is to search for .dat files. 
#==============================================================================#
getfile <- function(multFile=NULL, pattern=NULL){
# the default option is to grab all of the files in a directory matching a 
# specified pattern. If no pattern is set, all files will be listed.
	if (is.null(multFile) || toupper(multFile) == "YES"){
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
		csv <- sub("^.+?/(.+?\\.[a-z]{3})$", "\\1", csv)
		x <- csv
	}
	filepath <- list(files=x, path=path)
	return(filepath)
}
#==============================================================================#
# The calculation of the Index of Association and standardized Index of 
# Association.
#==============================================================================#
stdIa <- function(x){
#------------------------------------------------------------------------------#
# Testing for valid filetypes
#------------------------------------------------------------------------------#
	if (toupper(.readExt(x)) == "DAT") {
		pop <- read.fstat(x, missing = 0, quiet=T) 
	}
	else if (toupper(.readExt(x)) %in% c("STR", "STRU")) {
		pop <- read.structure(x, missing = 0, ask=FALSE, quiet=TRUE)
	}
	else if (toupper(.readExt(x)) == "GTX") {
		pop <- read.genetix(x, missing = 0, quiet=T) 
	}
	else if (toupper(.readExt(x)) == "GEN") {
		pop <- read.genepop(x, missing = 0, quiet=T) 
	}
	else {
		stop("File ext .dat, .str, .gtx, or .gen expected")
	}		
	numAlleles <- length(pop@loc.names)
	numIsolates <- length(pop@ind.names)
	# Creating this number is necessary because it is how the variance is
	# calculated.
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
	count <- 0
#--------------------------------------------------------------------------#
# This is the loop for analyzing the pairwise distances at each locus,
# also known as d. These values will be placed into two vectors representing
# the sum of d for each locus and the sum of d^2 for each locus. 
# 
# This will also calculate D, the pairwise comparison of all isolates over
# all loci. 
#--------------------------------------------------------------------------#
	for (loc in names(pop@loc.nall)){
		n <- count + 1
		count <- pop@loc.nall[loc] + count
		m <- count
		#print(paste(n, m, sep=":"))
		for(j in 1:numIsolates){
			# Loop for the rows of the distance matrix
			for(i in 1:numIsolates){
				z <- 0
				# Setting the contraint for the pairwise comparisons by the
				# equation (n(n-1))/2 where n = numIsolates 
				if(j < i && i <= numIsolates){
					z <- sum(abs(pop@tab[j,m:n]  
								-pop@tab[i,m:n]))
				}
				# The value of z (0, 1, or 2) is pushed into the matrix
				d.matrix[i,j] <- d.matrix[i,j]+z
				# This value is added to the matrix for pairwise comparison
				# of all the isolates over all loci as opposed to each locus
				D.matrix[i,j] <- D.matrix[i,j]+d.matrix[i,j]
			}
		}
		d.vector <- append(d.vector, sum(d.matrix))
		d2.vector <- append(d2.vector, sum(d.matrix^2))
		d.matrix[1:numIsolates,1:numIsolates] <- 0
	}
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
	print(paste("File Name:", x, sep=" "))
	print(paste("Index of Association:", Ia, sep=" "))
	print(paste("Standardized Index of Association (rbarD):", rbarD, sep=" "))
	# Saves the values of Ia, rbarD, and the filename into datastructures
	# that will be listed in a dataframe
	file.vector <- append(file.vector, x)
	Ia.vector <- append(Ia.vector, Ia)
	rbarD.vector <- append(rbarD.vector, rbarD)
	Iout <- list(Ia=Ia.vector, rbarD=rbarD.vector, File=file.vector)
	return(as.data.frame(Iout))
}

#==============================================================================#
# This function allows the user to analyze many separate files automatically
#==============================================================================#
stdIa.all <- function(x) {
	Iout <- NULL
	for(a in x){
		Iout <- rbind(Iout, stdIa(a))
	}
	return(Iout)
}

extract.info <- function(x) {
	x$clone <- as.factor(sub("^clone.(\\d{3}.\\d{2}).+?$)","\\1", x$File))
	x$rep <- as.numeric(sub("^clone.+?rep.(\\d{2}).+?$)","\\1", x$File))
	x$pop.size <- as.numeric(sub("^clone.+?pop.(\\d+?).+?$)","\\1", x$File))
	x$sex.rate <- (100-x$clone)/100
	return(x)
}