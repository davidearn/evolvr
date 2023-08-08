library(evolvr)

## MODELS
## original version of figure was made with tanh:
mtanh <- create_snowdrift_game(nPlayers=2, nGroups=2, K=0.0042)
## for gerf, we need L tauTurn, mmf instead: 
convert_David_to_L_tauturn_mmf(A=1, L=1000, K=0.0042)
## which yields:
L <- 1000
tauTurn <- 7 # 6.907755
mmf <- 0.05 # 0.05105
## which we use to set up the gerf model:
mgerf <-create_snowdrift_game(type="gerf", tauTurn=tauTurn, L=L, mmf=mmf, c1=1,
c2=0, nben=1, nGroups=2, nPlayers=2, nRepetitions = 10 )
## choose which model to use and call it msig:
msig <- mgerf

library(tikzDevice)
library(RColorBrewer)
if (!interactive()) tikz("figure3.tex", standAlone=TRUE, width=5, height=5)
par(mfrow=c(2,3))

## margins:  [default is c(b,l,t,r) = c(5, 4, 4, 2) + 0.1]
par(mar = c(4,4,1,1))
logscale <- "x"  # "" or "x"

## AS A FUNCTION OF Npop

##NpopVec <- 1:1024
##NpopVecOrig <- if (logscale == "x") c( 1:2^5, floor(2^seq(5.1,10,by=0.1)) ) else 1:32
NpopVecOrig <- seq(1,1024,by=1)
nGroupsVec <- c(4, 12, 16, 18, 20, 21)

start_plot <- function() {
    ## sigmoid
    x <- NpopVecOrig
    xlimmin <- if (logscale == "x") 2 else 0
    ## yinfmax <- max(XstarInf_fun(msig, nPlayers=x/max(nGroupsVec)))
    yinfmax <- if (msig$type == "tanh") 7.5 else 235
    xlim <- c(xlimmin,600) #1.02*max(x))
    ylim <- c(0.01,1.02*yinfmax)
    plot(NA,NA, xlim=xlim, ylim=ylim,
         xaxs="i", yaxs="i", log=logscale,
         type="n", bty="L", las=1,
         xlab="Population size $N$",
         ylab=if (dev_is_tikz()) "$X^*_\\infty$ and $X^*_N$" else "XstarInf and XstarN"
         )
    if (logscale == "x") {
        abline(h=seq(50,200,by=50), #1:7,
               v=c(5,10,20,50,100,200,500,1000),
           col="lightgrey", lty="dotted")
    } else {
        my_grid()
    }
}
add_legend <- function(msig) {
    p <- msig$parms
    if (msig$type == "tanh") {
        legend("top",bty="n",legend=paste(c("$K$","$L$","$A$"),"=",c(p$K,p$L,p$A)))
    } else if (msig$type == "gerf") {
        legend.text <-
            if (dev_is_tikz()) {
                paste(c("$m$","$L$","$\\tau_{\\rm turn}$"),"=",c(p$mmf,p$L,p$tauTurn))
            } else {
                paste(c("m","L","tauTurn"),"=",c(p$mmf,p$L,p$tauTurn))
            }
        legend("right",bty="n",legend=legend.text)
    }
}

if (msig$type == "tanh") {
    cols <- brewer.pal(n=length(nGroupsVec)+1, name="Greens")
    cols <- rev(cols[-1]) # faintest is too faint
} else { # assuming gerf
    ##cols <- brewer.pal(n=length(nGroupsVec), name="Set1")
    ##cols <- c("grey", "red", "darkgreen", "darkgrey")
    cols <- rep("black",length(nGroupsVec))
}

for (i in seq_along(nGroupsVec)) {
    start_plot()
    if (i==1) add_legend( msig )
    nGroups <- nGroupsVec[i]
    NpopCrit <- critical_Npop_from_nGroups(msig, nGroups=nGroups)
    ##abline(v=NpopCrit, col="yellow")
    ##abline(v=NpopCrit, col=transparent_colour(cols[i],alpha=1/4), lwd=5)
    legend.text <-
        if (dev_is_tikz()) {
            ##c(paste0("$n_{{}_{\\rm groups}} =", nGroups, "$"),
            c(paste0("$G =", nGroups, "$"),
              paste0("$N_{{}_{\\rm crit}} = ",
                     if (NpopCrit==Inf) "\\infty" else signif(NpopCrit,3),"$"))
        } else {
            c(paste0("nGroups = ", nGroups),
              paste0("NpopCrit = ", signif(NpopCrit,3)))
        }
    legend("topright", legend=legend.text, bty="n")
    yinfNpopCrit <- XstarInf_fun(msig, nPlayers=NpopCrit/nGroups)
    lines(x=c(NpopCrit,NpopCrit), y=c(0,yinfNpopCrit), lwd=5,
          ##col=transparent_colour(cols[i],alpha=1/4))
          col=transparent_colour("red",alpha=1/2))
    ##x <- ifelse(NpopVecOrig >= 2*nGroups, NpopVecOrig, NA)
    x <- ifelse(NpopVecOrig >= nGroups+1, NpopVecOrig, NA)
    ##yinf <- c(XstarInf_fun(msig, nPlayers=NpopVec/nGroupsVec), 0) # 0 needed for i+1 below
    yinf <- XstarInf_fun(msig, nPlayers=x/nGroups)
    lines(x, yinf, lwd=5, col=cols[i])
    y <- XstarN_fun(msig, nPlayers=x/nGroups, Npop=x)
    ##yyy <- c(yinf[i+1] + (1/3)*(yinf[i]-yinf[i+1]), yinf[i])
    ##lines(c(NpopCrit,NpopCrit), yyy, col=cols[i])
    ##abline(v=x[which.min(yinf)], col=cols[i])
    coop <- ifelse(x < NpopCrit, TRUE, FALSE)
    mycol <- brewer.pal(n=3, name="Blues")[3]
    lines(x[coop], y[coop])
    points(x[coop], y[coop], cex=0.7)
    points(x[coop], y[coop], pch=16, cex=0.5, col=mycol)
    ##points(x[coop], y[coop], pch=21, xpd=FALSE, cex=0.7,
           ####bg=transparent_colour(cols[i],alpha=0.5))
           ##bg=transparent_colour(cols[i],alpha=0.5),
           ##col=mycol)
    ##points(x, y, type="o", pch=21, bg=cols[i], xpd=NA)
    message("i=", i, ", nGroups=", nGroups, ", mmf=", msig$parms$mmf, ", NpopCrit=", NpopCrit)
}

if (FALSE) {

pu <- par("usr")
pux <- if (par("xlog")) 10^pu[1:2] else pu[1:2]
puy <- pu[3:4]
xleg <- pux[2]
yleg <- puy[1] + 0.9*(puy[2]-puy[1])
legend(x=xleg, y=yleg, xjust=1,
       legend=nGroupsVec,
       title=if (dev_is_tikz()) "$n_{{}_{\\rm groups}}$" else "nGroups",
       bty="n", lty=1, lwd=3, col=cols)
if (dev_is_tikz()) {
    adjfac <- if (logscale == "x") 0.325 else 0.063
} else {
    adjfac <- if (logscale == "x") 0.25 else 0.045
}
xleg <- pux[2] - adjfac*(pux[2]-pux[1])
legend(x=xleg, y=yleg, xjust=1,
       legend=rep("",length(nGroupsVec)), title="",
       bty="n", pch=21, pt.bg=cols)

} # end if (FALSE)

if (!interactive()) {
    dev.off()
    system("pdflatex figure3.tex > figure3.texlog")
}

## CHECK RELATIONSHIP BTW ESS AND ESSN AS nGroups --> Inf

f <- function(Npop) {
    return(XstarInf_fun(msig, nPlayers = Npop/nGroups) /
           XstarN_fun(msig, nPlayers = Npop/nGroups, Npop=Npop) )
}

par(mfrow=c(1,1))
Npop <- 2^(1:32)
plot(Npop, f(Npop), type="o", log="x", 
     xlab="Population size Npop",
     ylab="XstarInf / XstarN",
     xaxt="n")
message("nGroups = ", nGroups)
sfsmisc::eaxis(side=1)
abline(h=c(0:4,5*(1:7)),col="lightgrey",lty="dotted")
nGroupsVec <- c(20,20.99,21.01,22,23,32)
colvec <- c("blue","darkgrey","darkgreen","red","orange","cyan")
for (i in seq_along(nGroupsVec)) {
  nGroups <- nGroupsVec[i]
  points(Npop, f(Npop), type="o", col=colvec[i])  
}
