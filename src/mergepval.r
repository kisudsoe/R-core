# These codes are statistics for p-value meta-analysis

fishersMethod = function(x) {
    pchisq(-2*sum(log(x)), df=2*length(x), lower=F)
}
psumunif = function(x,n) {
    1/factorial(n)*sum(
        sapply(0:n, function(k) (-1)^k*choose(n,k)*ifelse(x>k,x-k,0)^(n))
    )
}
sumPvalsMethod = function(x) {
    n = length(x)
    if(n<10) psumunif(sum(x),n)
    else pnorm(sum(x),n/2,sqrt(n/12),lower=T)
}
binomMethod = function(x) {
    pbinom(sum(x<.05),size=length(x),prob=.05,lower=F)
}