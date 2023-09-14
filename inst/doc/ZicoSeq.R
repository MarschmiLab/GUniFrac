## -----------------------------------------------------------------------------
# install.packages("GUniFrac")

## ---- message=FALSE-----------------------------------------------------------
library(GUniFrac)

## ---- results = FALSE---------------------------------------------------------
data(throat.otu.tab)
data(throat.meta)
comm <- t(throat.otu.tab)
meta.dat <- throat.meta
meta.dat


## ---- results = FALSE---------------------------------------------------------
ZicoSeq.obj <- ZicoSeq(meta.dat = meta.dat, feature.dat = comm, 
                    grp.name = 'SmokingStatus', adj.name = 'Sex', feature.dat.type = "count",
                    # Filter to remove rare taxa
                    prev.filter = 0.2, mean.abund.filter = 0,  
                    max.abund.filter = 0.002, min.prop = 0, 
                    # Winsorization to replace outliers
                    is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                    # Posterior sampling 
                    is.post.sample = TRUE, post.sample.no = 25, 
                    # Use the square-root transformation
                    link.func = list(function (x) x^0.5), stats.combine.func = max,
                    # Permutation-based multiple testing correction
                    perm.no = 99,  strata = NULL, 
                    # Reference-based multiple stage normalization
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                    # Family-wise error rate control
                    is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)

## ---- fig.retina = 4, fig.width= 8, fig.height=8, results=FALSE, message=FALSE, warning=FALSE----
ZicoSeq.plot(ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = 0.1, text.size = 10,
             out.dir = NULL, width = 10, height = 6)

## ---- results = FALSE---------------------------------------------------------
comm.p <- t(t(comm) / colSums(comm))
ZicoSeq.obj.p <- ZicoSeq(meta.dat = meta.dat, feature.dat = comm.p, 
                    grp.name = 'SmokingStatus', adj.name = 'Sex', feature.dat.type = "proportion",
                    # Filter to remove rare taxa
                    prev.filter = 0.2, mean.abund.filter = 0,  max.abund.filter = 0.002, min.prop = 0, 
                    # Winsorization to replace outliers
                    is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                    # Posterior sampling will be automatically disabled
                    is.post.sample = FALSE, post.sample.no = 25, 
                    # Use the square-root transformation
                    link.func = list(function (x) x^0.5, function (x) x^0.25), stats.combine.func = max,
                    # Permutation-based multiple testing correction
                    perm.no = 99,  strata = NULL, 
                    # Reference-based multiple stage normalization
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                    # Family-wise error rate control
                    is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)
suppressWarnings(ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj.p, pvalue.type = 'p.adj.fdr', 
             cutoff = 0.1, text.size = 10, out.dir = NULL, width = 10, height = 6))                    


## ---- fig.retina = 4, fig.width= 14, fig.height=8, results=FALSE, message=FALSE, warning=FALSE----
comm.o <- comm[rowMeans(comm != 0) >= 0.2, ] + 1
comm.o <- log(t(t(comm.o) / colSums(comm.o)))
ZicoSeq.obj.o <- ZicoSeq(meta.dat = meta.dat, feature.dat = comm.o, 
		grp.name = 'SmokingStatus', adj.name = 'Sex', feature.dat.type = "other",
		# Filter will not be applied
		prev.filter = 0, mean.abund.filter = 0,  max.abund.filter = 0, min.prop = 0, 
		# Winsorization the top end 
		is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
		# Posterior sampling will be automatically disabled
		is.post.sample = FALSE, post.sample.no = 25, 
		# Identity function is used
		link.func = list(function (x) x), stats.combine.func = max,
		# Permutation-based multiple testing correction
		perm.no = 99,  strata = NULL, 
		# Reference-based multiple-stage normalization will not be performed
		ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
		# Family-wise error rate control
		is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)
    ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj.o,  pvalue.type = 'p.adj.fdr', 
             cutoff = 0.1, text.size = 10, out.dir = NULL, width = 10, height = 6)    

## -----------------------------------------------------------------------------
sessionInfo()

