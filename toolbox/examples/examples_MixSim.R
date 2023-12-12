library(MixSim)

### Test main functions

set.seed(1234, kind='Mersenne-Twister' , normal.kind = 'Inversion')
repeat {
   Q <- MixSim(BarOmega = 0.01, MaxOmega = 0.05, K = 5, p = 7)
   if (Q$fail == 0) break
}

overlap(Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
pdplot(Pi = Q$Pi, Mu = Q$Mu, S = Q$S)

### Example 1 of Section 3.1

set.seed(1234, kind='Mersenne-Twister' , normal.kind = 'Inversion')
(ex.1 <- MixSim(BarOmega = 0.05, MaxOmega = 0.15, K = 4, p = 5))
summary(ex.1)

### Example 2 of Section 3.1

set.seed(1234, kind='Mersenne-Twister' , normal.kind = 'Inversion')
(ex.2 <- MixSim(MaxOmega = 0.1, K = 3, p = 2, sph = TRUE, PiLow = 0.1))
summary(ex.2)

### Example 3 of Section 3.1

set.seed(1234, kind='Mersenne-Twister' , normal.kind = 'Inversion')
(ex.3 <- MixSim(BarOmega = 0.05, K = 2, p = 4, sph = TRUE, hom = TRUE,
                int = c(0, 10), eps = 1e-10))
summary(ex.3)

### Example for Section 3.2

data("iris", package = "datasets")
p <- ncol(iris) - 1
id <- as.integer(iris[, 5])
K <- max(id)

# estimate mixture parameters
Pi <- prop.table(tabulate(id))
Mu <- t(sapply(1:K, function(k){ colMeans(iris[id == k, -5]) }))
S <- sapply(1:K, function(k){ var(iris[id == k, -5]) })
dim(S) <- c(p, p, K)

overlap(Pi = Pi, Mu = Mu, S = S)

### Example 1 of Section 3.3 plot (a)

set.seed(1235, kind='Mersenne-Twister' , normal.kind = 'Inversion')
Q <- MixSim(MaxOmega = 0.20, BarOmega = 0.05, K = 5, p = 2)
A <- simdataset(n = 500, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
colors <- c("red", "green", "blue", "brown", "magenta")

par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(A$X, col = colors[A$id], pch = 19, cex = 0.8, 
     xlab = "", ylab = "", axes = FALSE)
box()

### Example 1 of Section 3.3 plot (b)

set.seed(1238, kind='Mersenne-Twister' , normal.kind = 'Inversion')
Q <- MixSim(MaxOmega = 0.20, BarOmega = 0.05, K = 5, p = 2)
A <- simdataset(n = 500, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
colors <- c("red", "green", "blue", "brown", "magenta")

par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(A$X, col = colors[A$id], pch = 19, cex = 0.8, 
     xlab = "", ylab = "", axes = FALSE)
box()

### Example 2 of Section 3.3 plot (c)

set.seed(1238, kind='Mersenne-Twister' , normal.kind = 'Inversion')
Q <- MixSim(MaxOmega = 0.1, K = 3, p = 2, int = c(0.2, 1))
A <- simdataset(n = 300, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, lambda = c(0.1, 10))
colors <- c("red", "green", "blue")

par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(A$X, col = colors[A$id], pch = 19, cex = 0.8, 
     xlab = "", ylab = "", axes = FALSE)
box()


### Example 2 of Section 3.3 plot (d)

set.seed(1238, kind='Mersenne-Twister' , normal.kind = 'Inversion')

Q <- MixSim(MaxOmega = 0.1, K = 3, p = 2, int = c(0.2, 1))
A <- simdataset(n = 300, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, lambda = c(10, 10))
colors <- c("red", "green", "blue")

par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(A$X, col = colors[A$id], pch = 19, cex = 0.8, 
     xlab = "", ylab = "", axes = FALSE)
box()

### Example 3 of Section 3.3 plot (a)

set.seed(1234, kind='Mersenne-Twister' , normal.kind = 'Inversion')

Q <- MixSim(BarOmega = 0.01, K = 4, p = 2)
A <- simdataset(n = 500, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, n.out = 10)
colors <- c("red", "green", "blue", "brown", "magenta")

par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(A$X, col = colors[A$id+1], pch = 19, cex = 0.8, 
     xlab = "", ylab = "", axes = FALSE)
box()


### Example 3 of Section 3.3 plot (b)

set.seed(1237, kind='Mersenne-Twister' , normal.kind = 'Inversion')

Q <- MixSim(BarOmega = 0.01, K = 4, p = 2)
A <- simdataset(n = 500, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, n.out = 10)
colors <- c("red", "green", "blue", "brown", "magenta")

par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(A$X, col = colors[A$id+1], pch = 19, cex = 0.8, 
     xlab = "", ylab = "", axes = FALSE)
box()


### Example 4 of Section 3.3 plot (c)

set.seed(1235, kind='Mersenne-Twister' , normal.kind = 'Inversion')

Q <- MixSim(MaxOmega = 0.1, K = 4, p = 1)
A <- simdataset(n = 300, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, n.noise = 1)
colors <- c("red", "green", "blue", "brown")

par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(A$X, col = colors[A$id], pch = 19, cex = 0.8, 
     xlab = "", ylab = "", axes = FALSE)
box()

### Example 4 of Section 3.3 plot (d)

set.seed(1236, kind='Mersenne-Twister' , normal.kind = 'Inversion')

Q <- MixSim(MaxOmega = 0.1, K = 4, p = 1)
A <- simdataset(n = 300, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, n.noise = 1)
colors <- c("red", "green", "blue", "brown")

par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(A$X, col = colors[A$id], pch = 19, cex = 0.8, 
     xlab = "", ylab = "", axes = FALSE)
box()


### Example of Section 3.4 plot (a)
data("iris", package = "datasets")
p <- ncol(iris) - 1
id <- as.integer(iris[, 5])
K <- max(id)

# estimate mixture parameters
Pi <- prop.table(tabulate(id))
Mu <- t(sapply(1:K, function(k){ colMeans(iris[id == k, -5]) }))
S <- sapply(1:K, function(k){ var(iris[id == k, -5]) })
dim(S) <- c(p, p, K)

pdplot(Pi = Pi, Mu = Mu, S = S)


### Example of Section 3.4 plot (b)

set.seed(1234, kind='Mersenne-Twister' , normal.kind = 'Inversion')
Q <- MixSim(BarOmega = 0.001, K = 6, p = 4)
pdplot(Pi = Q$Pi, Mu = Q$Mu, S = Q$S)


### Example of Section 3.4 plot (c)

set.seed(1232, kind='Mersenne-Twister' , normal.kind = 'Inversion')
Q <- MixSim(BarOmega = 0.05, K = 6, p = 4)
pdplot(Pi = Q$Pi, Mu = Q$Mu, S = Q$S)

###### SONO ARRIVATO QUI'


### Example of Section 4

library("mclust")
library("cluster")

p <- 4
K <- 6
n <- 200
n.files <- 100
omega <- c(0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.01, 0.005, 0.001)

set.seed(1234)

AR.Mclust <- rep(NA, n.files)
AR.Kmeans <- rep(NA, n.files)
AR.PAM <- rep(NA, n.files)
AR.Ward <- rep(NA, n.files)
	
P.Mclust <- rep(NA, n.files)
P.Kmeans <- rep(NA, n.files)
P.PAM <- rep(NA, n.files)
P.Ward <- rep(NA, n.files)
	
VI.Mclust <- rep(NA, n.files)
VI.Kmeans <- rep(NA, n.files)
VI.PAM <- rep(NA, n.files)
VI.Ward <- rep(NA, n.files)

for (bo in omega){
	
	for (i in 1:n.files){

		Q <- MixSim(BarOmega = bo, K = K, p = p)
		A <- simdataset(n = n, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)

		# Use try(...) to catch Mclust failure instances
		id.Mclust <- try(Mclust(A$X, G = K, model = "VVV")$class)
		# Provide error ids if a failure occurs
		if (inherits(id.Mclust, what = "try-error") || (is.null(id.Mclust)))
			id.Mclust <- rep(1, n)

		id.Kmeans <- kmeans(A$X, K, iter.max = 1000, nstart = 10)$cluster
		id.PAM <- pam(A$X, K)$clustering
		d <- dist(A$X, method = "euclidean")
		id.Ward <- cutree(hclust(d = d, method = "ward"), k = K)

		AR.Mclust[i] <- RandIndex(A$id, id.Mclust)$AR
		AR.Kmeans[i] <- RandIndex(A$id, id.Kmeans)$AR
		AR.PAM[i] <- RandIndex(A$id, id.PAM)$AR
		AR.Ward[i] <- RandIndex(A$id, id.Ward)$AR

		P.Mclust[i] <- ClassProp(A$id, id.Mclust)
		P.Kmeans[i] <- ClassProp(A$id, id.Kmeans)
		P.PAM[i] <- ClassProp(A$id, id.PAM)
		P.Ward[i] <- ClassProp(A$id, id.Ward)

		VI.Mclust[i] <- VarInf(A$id, id.Mclust)
		VI.Kmeans[i] <- VarInf(A$id, id.Kmeans)
		VI.PAM[i] <- VarInf(A$id, id.PAM)
		VI.Ward[i] <- VarInf(A$id, id.Ward)

	}

	cat("Omega =", bo, "\n\n")

	cat("Mclust: AR mean =", mean(AR.Mclust), " sd =", sd(AR.Mclust), "\n")
	cat("Mclust: P mean =", mean(P.Mclust), " sd =", sd(P.Mclust), "\n")
	cat("Mclust: VI mean =", mean(VI.Mclust), " sd =", sd(VI.Mclust), "\n\n")

	cat("Kmeans: AR mean =", mean(AR.Kmeans), " sd =", sd(AR.Kmeans), "\n")
	cat("Kmeans: P mean =", mean(P.Kmeans), " sd =", sd(P.Kmeans), "\n")
	cat("Kmeans: VI mean =", mean(VI.Kmeans), " sd =", sd(VI.Kmeans), "\n\n")

	cat("PAM: AR mean =", mean(AR.PAM), " sd =", sd(AR.PAM), "\n")
	cat("PAM: P mean =", mean(P.PAM), " sd =", sd(P.PAM), "\n")
	cat("PAM: VI mean =", mean(VI.PAM), " sd =", sd(VI.PAM), "\n\n")

	cat("Ward: AR mean =", mean(AR.Ward), " sd =", sd(AR.Ward), "\n")
	cat("Ward: P mean =", mean(P.Ward), " sd =", sd(P.Ward), "\n")
	cat("Ward: VI mean =", mean(VI.Ward), " sd =", sd(VI.Ward), "\n\n")

}



