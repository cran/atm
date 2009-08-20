\name{term}
\alias{clean.explanatory}
\alias{clean.fitted}
\alias{clean.residuals}
\alias{clean.response}
\alias{clean.weights}
\alias{force.residuals}
\alias{gf.term}
\alias{oospfitted.term}
\alias{oosppredict.term}
\alias{oospprint.summary.term}
\alias{oospprint.term}
\alias{oospprint.term.design}
\alias{oospprint.term.estimate}
\alias{oospresiduals.term}
\alias{oospsummary.term}
\alias{oospweights.term}
\alias{pstars}
\alias{seq.term}
\alias{simple.gf}
\alias{term}
\alias{term.confband}
\alias{term.design}
\alias{termenv}
\alias{term.estimate}
\alias{term.label}
\alias{term.predband}
\alias{term.pump}
\alias{term.reweight}
\alias{term.validate}
\title{Term Objects}
\description{Create term objects.}
\usage{
clean.explanatory (t)
clean.fitted (t)
clean.residuals (t, y)
clean.response (t, y)
clean.weights (t)
force.residuals (t, y)
\method{gf}{term}(t, y, method=c ("rsq", "mar"), \dots)
\method{oospfitted}{term}(t, \dots)
\method{oosppredict}{term}(t, y, u=t$x, type=c ("argument", "value"), cl=0.95, \dots)
\method{oospprint}{summary.term}(s, \dots)
\method{oospprint}{term}(t, \dots)
\method{oospprint}{term.estimate}(e, \dots)
\method{oospresiduals}{term}(t, y, \dots)
\method{oospsummary}{term}(t, \dots)
\method{oospweights}{term}(t, \dots)
pstars (ps)
\method{seq}{term}(t, n=200, \dots)
simple.gf (yb, yh, r, method=c ("rsq", "mar") )
term (d, E, e, s, x, z=NULL, validate=TRUE, clean=TRUE,
	nr=length (x), nc=1, nv, valid=NULL, w=NULL)
term.confband (t, y, u=seq (t), type=standard, cl=0.95, \dots)
term.design ()
termenv ()
term.estimate ()
term.label (t)
term.predband (t, y, \dots)
term.pump (t, y, sc=NULL, \dots)
term.reweight (t, w)
term.validate (t, vz)
}
\arguments{
\item{d}{}
\item{E}{}
\item{e}{}
\item{s}{}
\item{x}{}
\item{z}{}
\item{validate}{}
\item{clean}{}
\item{nr}{}
\item{nc}{}
\item{nv}{}
\item{valid}{}
\item{w}{}
\item{t}{}
\item{y}{}
\item{u}{}
\item{sc}{}
\item{vz}{}
\item{type}{}
\item{cl}{}
\item{n}{}
\item{method}{}
\item{yb}{}
\item{yh}{}
\item{r}{}
\item{ps}{}
\item{\dots}{}
}




