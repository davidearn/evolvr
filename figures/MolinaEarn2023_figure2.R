library(evolvr)

## MODELS
mtanh <- create_snowdrift_game(nPlayers=5, nGroups=2)
mq <- create_snowdrift_game("quadratic")
mgerf <-create_snowdrift_game(type="gerf", tauTurn=15, L=10, mmf=1.5, c1=1,
c2=0, nben=1, nGroups=20, nPlayers=2, nRepetitions = 10 )

msig <- mgerf

library(tikzDevice)
library(RColorBrewer)
if (!interactive()) tikz("figure2.tex", standAlone=TRUE, width=5, height=5)
## margins:  [default is c(b,l,t,r) = c(5, 4, 4, 2) + 0.1]
par(mar = c(4,4,1,1))
logscale <- ""  # "" or "x"

## AS A FUNCTION OF Npop

##NpopVec <- 1:1024
NpopVec <- if (logscale == "x") c( 1:2^5, floor(2^seq(5.1,10,by=0.1)) ) else 1:32
nPlayersVec <- 2^(1:4)

## sigmoid
x <- NpopVec
xlimmin <- if (logscale == "x") 2 else 0
yinfmax <- max(XstarInf_fun(msig, nPlayers=nPlayersVec))
plot(NA,NA, xlim=c(xlimmin,1.02*max(x)), ylim=c(0,1.02*yinfmax),
     xaxs="i", yaxs="i", log=logscale,
     type="n", bty="L", las=1,
     xlab="Population size $N$",
     ylab=if (dev_is_tikz()) "$X^*_\\infty$ and $X^*_N$" else "XstarInf and XstarN")
if (logscale == "x") {
    abline(h=1:4, v=c(5,10,20,50,100,200,500,1000), col="lightgrey", lty="dotted")
} else {
    my_grid()
}

cols <- brewer.pal(n=length(nPlayersVec)+1, name="Blues")
cols <- rev(cols[-1]) # faintest is too faint
yinf <- c(XstarInf_fun(msig, nPlayers=nPlayersVec), 0) # 0 needed for i+1 below
for (i in seq_along(nPlayersVec)) {
    nPlayers <- nPlayersVec[i]
    y <- XstarN_fun(msig, nPlayers=nPlayers, Npop=x)
    abline(h=yinf[i], lwd=3, col=cols[i])
    NpopCrit <- critical_Npop(msig, nPlayers=nPlayers)
    yyy <- c(yinf[i+1] + (1/3)*(yinf[i]-yinf[i+1]), yinf[i])
    lines(c(NpopCrit,NpopCrit), yyy, col=cols[i])
    coop <- ifelse(x > NpopCrit, TRUE, FALSE)
    points(x[coop], y[coop], type="o", pch=21, bg=cols[i], xpd=NA)
}
pu <- par("usr")
pux <- if (par("xlog")) 10^pu[1:2] else pu[1:2]
puy <- pu[3:4]
xleg <- pux[2]
yleg <- puy[1] + 0.9*(puy[2]-puy[1])
legend(x=xleg, y=yleg, xjust=1,
       legend=nPlayersVec,
       ##title=if (dev_is_tikz()) "$n_{{}_{\\rm players}}$" else "nPlayers",
       title=if (dev_is_tikz()) "$n$" else "nPlayers",
       bty="n", lty=1, lwd=3, col=cols)
if (dev_is_tikz()) {
    adjfac <- if (logscale == "x") 0.325 else 0.063
} else {
    adjfac <- if (logscale == "x") 0.25 else 0.045
}
xleg <- pux[2] - adjfac*(pux[2]-pux[1])
legend(x=xleg, y=yleg, xjust=1,
       legend=rep("",length(nPlayersVec)), title="",
       bty="n", pch=21, pt.bg=cols)

if (!interactive()) {
    dev.off()
    system("pdflatex figure2.tex > figure2.texlog")
}
