\name{smooth}
\alias{atm.hs}
\alias{atm.ss}
\alias{evaluate.smooth}
\alias{fit.smooth}
\alias{hardwf}
\alias{oospprint.smooth.estimate}
\alias{oospprint.summary.smooth}
\alias{oospsummary.smooth}
\alias{reset.smooth}
\alias{smooth}
\alias{smooth.design}
\alias{smooth.estimate}
\alias{stdwf}
\title{Smooth Terms}
\description{Create smooth terms.}
\usage{
atm.hs (x, y, ns, sx, deg, wf, smoothness, w)
atm.ss (x, y, ns, sx, deg, wf, smoothness, w)
\method{evaluate}{smooth}(t, u, \dots)
\method{fit}{smooth}(t, y, pump, \dots)
hardwf (a, b, wf, m, w, x)
\method{oospprint}{smooth.estimate}(e, \dots)
\method{oospprint}{summary.smooth} (s, \dots)
\method{oospsummary}{smooth}(t, y, \dots)
\method{reset}{smooth}(t, \dots)
smooth (x, y=NULL, ns=50, refit=2, deg=2, E, wf, smoothness=NULL, s, w=NULL)
smooth.design (ns, refit, deg, wf, smoothness, x)
smooth.estimate (ns, sy=rep (0, ns) )
stdwf (m, w, x)
}
\arguments{
\item{x}{}
\item{y}{}
\item{ns}{}
\item{sx}{}
\item{deg}{}
\item{wf}{}
\item{smoothness}{}
\item{w}{}
\item{t}{}
\item{u}{}
\item{sy}{}
\item{pump}{}
\item{a}{}
\item{b}{}
\item{m}{}
\item{e}{}
\item{s}{}
\item{refit}{}
\item{E}{}
\item{\dots}{}
}



