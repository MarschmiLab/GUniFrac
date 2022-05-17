
dICC <- function (dist.mat, strata) {
	n <- length(strata)
	strata <- factor(strata)
	m <- nlevels(strata)
	ns <- table(strata)
	weights <- strata
	levels(weights) <- ns
	weights <- as.numeric(as.character(weights))
	
	dist.mat <- as.matrix(dist.mat^2)
	
	with.mask <- dist.mat[, ]
	with.mask[, ] <- 0
	for (ID in levels(strata)) {
		ind <- strata == ID
		with.mask[ind, ind] <- 1
	}
	btw.mask <- !with.mask
	btw.mask[upper.tri(btw.mask)] <- with.mask[upper.tri(with.mask)] <- 0
	diag(btw.mask) <- diag(with.mask) <- 0
	
	SSW <- sum(dist.mat * with.mask / weights)
	SST <- (sum(dist.mat * btw.mask) + sum(dist.mat * with.mask)) / n
	SSB <- SST - SSW
	
	R2 <- SSB / SST
	
	SE <- SSW / (n - m)
	SB <- (SST - (n-1) * SE) / (n - sum(ns^2) / n)
	
	ICC <- SB / (SB + SE)
	
	return(list(ICC = ICC))
}



# For N by M design only
dICC.SE.asympt <- function (dist.mat, strata) {
	
	strata <- factor(strata)
	N <- nlevels(strata)
	if (length(strata) %% N != 0) {
		stop("Asmptotic SE calculation only supports equal numbers of technical replicates! Consider using bootstrap SE instead!\n")
	}

	M <- length(strata) / N
	
	M.exp <- outer(1:M, 1:M, paste)
	M.exp <- M.exp[lower.tri(M.exp)]
	M.exp <- strsplit(M.exp, ' ')
	
	# first order the mat
	strata <- as.numeric(strata)
	obj <- sort(strata, index.return = TRUE)
	mat <- dist.mat[obj$ix, obj$ix]
	strata <- obj$x
	index <- rep(1:M, N)
	
	mat <- as.matrix(mat^2)
	
	# Calculate Tn and Gn
	TN <- 0
	for (i in 1:length(M.exp)) {
		index1 <- as.numeric(M.exp[[i]])
		TN <- TN + sum(mat[cbind(which(index == index1[1]), which(index == index1[2]))])
	}
	TN <- TN * 2 / N / M / (M - 1)
	
	GN <- 0
	for (i in 1:M) {
		for (j in 1:M) {
			mat1 <- mat[which(index == i), which(index == j)]
			diag(mat1) <- 0
			GN <- GN + sum(mat1)
		}
	}
	GN <- GN / N / (N - 1) / M^2
	
	# Calculate variance component s11
	s11 <- 0
	for (i in 1:length(M.exp)) {
		index1 <- as.numeric(M.exp[[i]])
		for (j in 1:length(M.exp)) {
			index2 <- as.numeric(M.exp[[j]])
			s11 <- s11 + mean(mat[cbind(which(index == index1[1]), which(index == index1[2]))] * 
							mat[cbind(which(index == index2[1]), which(index == index2[2]))]) - TN^2
		}
	}
	s11 <- 4 * s11 / M^2 / (M - 1)^2 
	
	# Calculate variance component s22
	s22 <- 0
	for (s in 1:M) {
		for (u in 1:M) {
			t1 <- 0
			t2 <- 0
			for (v in 1:M) {
				mat1 <- (mat[which(index == s), which(index == v)])
				mat2 <- (mat[which(index == u), which(index == v)]) 
				t1 <- t1 + (rowSums(mat1) - diag(mat1))
				t2 <- t2 + (rowSums(mat2) - diag(mat2))			
			}
			t1 <- t1 / (N - 1) / M
			t2 <- t2 / (N - 1) / M
			
			s22 <-  s22 +  mean(t1 * t2) - GN^2
		}
	}
	s22 <- 4  / M^2 * s22
	
	# Calculate variance component s12
	s12 <- 0
	for (i in 1:length(M.exp)) {
		index1 <- as.numeric(M.exp[[i]])
		vec <- mat[cbind(which(index == index1[1]), which(index == index1[2]))]
		for (s in 1:M) {
			t1 <- 0
			for (u in 1:M) {
				mat1 <- (mat[which(index == s), which(index == u)])
				t1 <- t1 + (rowSums(mat1) - diag(mat1))
			}
			t1 <- t1 / (N - 1) / M
			s12 <- s12 + mean(vec * t1) - TN * GN
		}
	}
	s12 <- 4 * s12 / M^2 / (M - 1)
	
	# http://www.stat.cmu.edu/~hseltman/files/ratio.pdf
	SE <- sqrt(TN^2 / GN^4 * s22 + s11 / GN^2 - 2 * TN / GN^3 * s12) / sqrt(N)
#	return(list(s11 = s11, s22 = s22, s12 = s12, SE = SE,
#					TN = TN, GN = GN, ICC = 1 - TN / GN))
	return(list(ICC = 1 - TN / GN, SE = SE))
}




dICC.SE.bt <- function (dist.mat, strata, B = 199) {
	
	paste2 <- function(x, collapse) {
		if (length(x) == 1) {
			return(paste(x))
		} else {
			return(paste(x, collapse=collapse))
		}
	}
	
# subject-based ID
	boot.ind <- function (ind, subject) {
		subject <- factor(subject)
		temp <- tapply(ind, subject, paste2, collapse = ' ')
		temp <- temp[sample(levels(subject), replace = TRUE)]
		strata <- paste0(names(temp), 1:length(temp))
		temp <- strsplit(temp, ' ')
		return(list(ind=as.numeric(unlist(temp)), strata=factor(rep(strata, sapply(temp, length)))))
	}
	
	ind <- 1:nrow(dist.mat)
	res <- sapply(1:B, function (i) {
				temp <- boot.ind(ind, strata)
				ind0 <- temp$ind
				strata0 <- temp$strata
				dICC(dist.mat[ind0, ind0], strata0)$ICC	
			})
	return(list(ICC = dICC(dist.mat, strata)$ICC, SE = sd(res)))
}

#
#rdirichlet1 <- function (n = 1, alpha) {
#	p <- length(alpha)
#	Gam <- matrix(rgamma(n * p, shape = alpha), ncol = p, byrow = TRUE)
#	Gam / rowSums(Gam)
#}
#
#rdirichlet2 <- function(Pi, t2) {
#	
#	n <- nrow(Pi)
#	p <- ncol(Pi)
#	alpha <- t(Pi * t2)
#	Gam <- matrix(rgamma(n * p, shape = alpha), ncol = p, byrow = TRUE)
#	Gam / rowSums(Gam)
#	
#}
