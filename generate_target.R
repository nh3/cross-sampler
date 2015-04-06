#!/usr/bin/env Rscript
library(getopt);

opt.spec <- matrix(ncol=5, byrow=T, data=c(
    'input',  'i', 1, 'character', 'input intervals with GC and DP info',
    'ref',    'r', 1, 'character', 'reference intervals with GC and DP info',
    'output', 'o', 1, 'character', 'output intervals with DP info',
    'method', 'm', 1, 'character', '"parametric" or "empirical"',
    'scale',  's', 1, 'integer',   'scale of auto-correlation for the parametric method, integer from 1(weak) to 10(strong), default 5',
    'help',   'h', 0, 'logical',   'print this message'
));
script.name <- sub('^.*/', '', strsplit(commandArgs(FALSE)[4],"=")[[1]][2]);
usage <- getopt(spec=opt.spec, command=script.name, usage=T);
opt <- getopt(spec=opt.spec);

if (!is.null(opt$help)) {
    cat(usage);
    q(status=0);
}

if (is.null(opt$input)) {
    cat('missing --input\n');
    q(status=1);
}

if (is.null(opt$ref)) {
    cat('missing --ref\n');
    q(status=1);
}

if (is.null(opt$scale)) {
    opt$scale <- 5;
} else {
    if (opt$scale < 1) opt$scale <- 1;
    if (opt$scale > 10) opt$scale <- 10;
}

if (is.null(opt$output)) opt$output <- '';

if (is.null(opt$method)) opt$method <- 'empirical';
choice <- pmatch(opt$method, c('parametric', 'empirical'));
if (is.na(choice)) {
    cat('unsupported --method\n');
    q(status=1);
}

fix.negative <- function(x) {x[x<0] <- 0; x;};

preprocess <- function(dat) {
    if ('dp' %in% colnames(dat)) {
        dat$dp <- fix.negative(dat$dp);
        dat$dp <- log10(dat$dp+1);
    }
    dat$gc <- fix.negative(dat$gc);
    return(dat);
}

make.gc.breaks <- function(gc) {
    gc.breaks = unique(quantile(gc, seq(0,1,0.01)));
    if (gc.breaks[1]>0) gc.breaks[1] = 0;
    if (gc.breaks[length(gc.breaks)]<100) gc.breaks[length(gc.breaks)] = 100;
    gc.breaks;
}

parametric.resampling <- function(src, tgt, acf.scale=1) {
    N.src <- nrow(src);
    N.tgt <- nrow(tgt);
    src.n.bins <- length(unique(src$gc.bins));
    tgt.n.bins <- length(unique(tgt$gc.bins));

    calc.mean <- function(x) mean(x[x>0]);
    calc.sd <- function(x) sd(x[x>0]);
    calc.p <- function(x) sum(x==0)/length(x);
    gc.stats <- data.frame(vapply(list(mu=calc.mean,sigma=calc.sd,p0=calc.p), function(f) tapply(src$dp, src$gc.bins, f), numeric(src.n.bins)));

    K <- tapply(seq(N.src-1), paste(src$gc.bins[-N.src],src$gc.bins[-1]), function(x) x);
    acf.by.gc <- vapply(seq(src.n.bins), function(i) {
        vapply(seq(src.n.bins), function(j) {
            k <- K[[paste(i,j)]];
            if (is.null(k)) {
                0;
            } else {
                cor(src$dp[k], src$dp[k+1]);
            }
        }, numeric(1));
    }, numeric(src.n.bins));
    coef <- vapply(seq(src.n.bins), function(i) {
        vapply(seq(src.n.bins), function(j) {
            rho <- acf.by.gc[i,j]^acf.scale;
            c(gc.stats$sigma[i]/gc.stats$sigma[j]*rho, gc.stats$mu[i]-gc.stats$sigma[i]/gc.stats$sigma[j]*rho*gc.stats$mu[j], sqrt(1-rho^2)*gc.stats$sigma[i])
        }, numeric(3))
    }, array(0,dim=c(3,src.n.bins)));

    dp <- rep(0, N.tgt);
    R <- runif(N.tgt);
    prev.chr <- 0;
    for (i in 1:N.tgt) {
        chr <- tgt$chr[i];
        k <- tgt$gc.bins[i];
        if (R[i]<gc.stats$p0[k]) {
            dp[i] <- 0;
        } else {
            if (chr != prev.chr || dp[i-1]==0) {
                dp[i] <- max(0, rnorm(1,mean=gc.stats$mu[k],sd=gc.stats$sigma[k]));
            } else {
                k0 = tgt$gc.bins[i-1];
                dp[i] <- max(0, rnorm(1, mean=coef[1,k0,k]*dp[i-1]+coef[2,k0,k], sd=coef[3,k0,k]));
            }
        }
        prev.chr <- chr;
    }
    return(dp);
}

empirical.resampling <- function(src, tgt) {
    N.src <- nrow(src);
    N.tgt <- nrow(tgt);

    dp.by.gc <- tapply(src$dp, src$gc.bins, function(x) x);
    dp = vapply(seq(N.tgt), function(i) sample(dp.by.gc[[tgt$gc.bins[i]]], size=1), numeric(1));

    return(dp);
}

# read reference
src = read.tsv(opt$ref);
colnames(src) <- c('chr','start','end','gc','dp');
src <- preprocess(src);
orig.src <- src;
src.k.gc0 <- src$gc == 0;
src <- src[!src.k.gc0,];
gc.breaks <- make.gc.breaks(src$gc);
src$gc.bins <- findInterval(src$gc, gc.breaks);

# read input
tgt = read.tsv(opt$input);
colnames(tgt) <- c('chr','start','end','gc','dp');
tgt <- preprocess(tgt);
orig.tgt <- tgt;
tgt.k.gc0 <- tgt$gc == 0;
tgt <- tgt[!tgt.k.gc0,]
tgt$gc.bins <- findInterval(tgt$gc, gc.breaks);

# re-sampling
dp <- list(parametric.resampling, empirical.resampling)[[choice]](src, tgt, 1/opt$scale);

# write output
output <- data.frame(chr=orig.tgt$chr,start=orig.tgt$start,end=orig.tgt$end,dp=round(10^orig.tgt$dp)-1,tgt.dp=integer(nrow(orig.tgt)));
output$tgt.dp[!tgt.k.gc0] <- round(10^dp)-1;
write.tsv(output, opt$output, col.names=F, row.names=F);
