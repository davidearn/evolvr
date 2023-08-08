library(evolvr)

## original default parameters using type=tanh:
nsg <- create_snowdrift_game()
## parameters used for selmut simulations:
##nsg <-create_snowdrift_game(type="gerf", tauTurn=50, L=.5, mmf=0.1, c1=1, c2=0, nben=10, nGroups=10, nPlayers=15, nRepetitions = 10 )
## gerf parameters that yield a nice sigmoid:
nsg <-create_snowdrift_game(type="gerf", tauTurn=15, L=10, mmf=1.5, c1=1,
c2=0, nben=1, nGroups=20, nPlayers=2, nRepetitions = 10 )

summary(nsg)

myxmax <- 30
asp <- NA

if (!interactive()) {
    library(tikzDevice)
    tikz("figure1.tex", standAlone=TRUE, width=3, height=8)
}
##par(mfrow=c(3,1))
layout( matrix(c(1,2,3,3), 4, 1, byrow=TRUE) )
## margins:  [default is c(b,l,t,r) = c(5, 4, 4, 2) + 0.1]
par(mar=c(4.5,3,1,1))
lwd <- 5
plot(nsg, what=c("cost","benefit","fitness"),
     cex.axis=1.5, cex.lab=1.5, yaxs="i", asp=asp,
     residentStrategy=0, xmax=myxmax, xmax.group=myxmax,
     ylim.fitness=c(-10.5,9.5),
     show.main=FALSE,
     show.annotation=FALSE,
     show.leftlab=TRUE,
     col.resident="darkgreen",
     label.resident="",
     draw.special.strategies = TRUE,
     draw.XstarInf=FALSE,
     draw.XstarN=FALSE,
     ##leftlab.fitness=paste0("Fitness\n","{\\tiny$n_{{\\rm players}} = ", nsg$parms$nPlayers, "$}"),
     xlab.total="Total group contribution $\\tau$",
     xlab.focal="Focal individual's contribution $x$"
)
##legend("bottomleft", bty="n", cex=1.25, legend=c("Residents defect","",""))
mycol <- transparent_colour("darkgreen",alpha=1/2)
faintcol <- transparent_colour("darkgreen",alpha=1/4)
plot(nsg, what=c("fitness"), add=TRUE,
     ##cex.axis=1.5, cex.lab=1.5, yaxs="i", asp=asp,
     ##ylim=c(-4,9), xpd.fitness=FALSE,
     residentStrategy=5,
     xmax=myxmax,
     show.main=FALSE,
     show.annotation=FALSE,
     show.leftlab=FALSE,
     draw.special.strategies = TRUE,
     draw.XstarInf=FALSE,
     draw.XstarN=FALSE,
     draw.grid=FALSE,
     col.fitness=mycol,
     col.resident=mycol,
     ##lwd.resident=lwd/2,
     label.resident="", ##"\\qquad residents play $<$ ESS$_N$",
     cex.label.resident=0.6,
     xlab.focal="Focal individual's contribution $x$"
)
plot(nsg, what=c("fitness"), add=TRUE, xpd.fitness=FALSE,
     residentStrategy=nsg$parms$XstarN, # play ESSN
     xmax=myxmax,
     show.main=FALSE,
     show.annotation=FALSE,
     show.leftlab=FALSE,
     draw.special.strategies = TRUE,
     draw.XstarInf=FALSE,
     draw.XstarN=FALSE,
     draw.grid=FALSE,
     col.fitness=faintcol,
     col.resident=faintcol, #"lightgrey",
     ##lwd.resident=lwd/2,
     label.resident="", ##"residents play ESS$_N$",
     cex.label.resident=0.6,
     ##xlab.focal="Focal individual's contribution $x$"
)
if (FALSE) {
legend("bottomleft", bty="n", cex=1.25,
       legend=c(##"Residents cooperate",
                ##paste0("Resident strategy $=X^*_{\\infty}=\\frac{\\displaystyle\\tau_{\\rm max}}{\\displaystyle n_{\\rm players}}=", signif(nsg$parms$XstarInf,3), "$"),
                paste0("$\\!\\!\\!\\!\\!\\!n_{{\\rm players}} = ", nsg$parms$nPlayers, "$")))
}
legend("topright", bty="n", cex=1,
       legend=c("", "Residents play 9.6 (ESS$_N$)",
                "Residents play 5 ($<$ ESS$_N$)",
                "Residents play 0 (defect)"),
       col=c(NA, faintcol,mycol,"darkgreen"),
       lty="solid", lwd=lwd,
       ##title=paste0("\\LARGE$n_{{\\rm players}} = ", nsg$parms$nPlayers, "$")
       )
##legend.text <- paste0("{$n_{{\\rm players}} = ", nsg$parms$nPlayers, "$}")
legend.text <- paste0("{$n = ", nsg$parms$nPlayers, "$}")
pu <- par("usr")
xleg <- pu[1] + 0.03*(pu[2]-pu[1])
yleg <- pu[3] + 0.93*(pu[4]-pu[3])
legend(x=xleg,y=yleg, bty="n", cex=1.5,
       legend=legend.text)

if (!interactive()) {
    dev.off()
    system("pdflatex figure1 > figure1.texlog")
}

