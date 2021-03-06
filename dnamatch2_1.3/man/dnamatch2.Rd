% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dnamatch2.R
\name{dnamatch2}
\alias{dnamatch2}
\title{dnamatch2}
\usage{
dnamatch2(fn, freqfile, reffold = NULL, sameCID = FALSE,
  betweensamples = TRUE, Thist = Inf, threshMAC = 0.75, threshLR = c(10,
  100), threshHeight = 200, threshStutt = 0.1, threshMaj = 0.6,
  minLocStain = 3, minLocMaj = 3, pC = 0.05, lambda = 0.01,
  kit = "ESX17", N = 5000, searchtime = Sys.time(), SIDvec = NULL,
  TAvec = NULL, CIDvec = NULL, timediff = NULL, TAptrn = NULL,
  SIDpat = NULL, BIDpat = NULL, printHistPlots = FALSE)
}
\arguments{
\item{fn}{A folder with stain files (possible a vector). Full directory must be given.}

\item{freqfile}{A file containing population frequencies for alleles. Full directory must be given.}

\item{reffold}{A folder with stored personal references. Default is no references. Full directory must be given.}

\item{sameCID}{Boolean whether matches within same case ID number should be allowed.}

\item{betweensamples}{A boolean of whether between samples are searched. Default is TRUE.}

\item{Thist}{Number of days back in time for a TA-file to be imported (stains).}

\item{threshMAC}{Threshold for a allele match. A ratio [0,1]}

\item{threshLR}{Threshold for a match. Can be a vector (qual,quan).}

\item{threshHeight}{Acceptable peak height (rfu) in stains.}

\item{threshStutt}{Acceptable stutter ratio in stains (relative peak heights) .}

\item{threshMaj}{If second largest allele has ratio (relative to the largest allele) above this threshold, then second allele is part of the major profile (used for extracting major from mixture). If the relative peak height between second and third largest allele has ratio greater than this threshold, no major is assigned.}

\item{minLocStain}{Number of minimum loci in stain profiles in order to be evaluated.}

\item{minLocMaj}{Number of minimum loci which are required to be considered for a major component (extracted from each stain profile).}

\item{pC}{Assumed drop-in rate per marker (parameter in the LR models).}

\item{lambda}{Assumed hyperparameter lambda for peak height drop-in model (parameter in quantiative model).}

\item{kit}{The shortname of the kit used for the samples. This allows for taking into account degradation. Use getKit function in euroformix R-package.}

\item{N}{Database size used for estimating population frequencies. Used to assign new allele frequenceis. Default is 5000.}

\item{searchtime}{An object with format as returned from Sys.time(). Default is Sys.time().}

\item{SIDvec}{A vector with Sample-ID numbers which are considered in the search (i.e. a search filter).}

\item{TAvec}{A vector with TA-files, which are considered in the search (i.e. a search filter).}

\item{CIDvec}{A vector with Case-ID numbers which are considered in the search (i.e. a search filter).}

\item{timediff}{Timedifference (in days) allowed between matching reference and target. Default is NULL (not used).}

\item{TAptrn}{Filename structure of stain files to evaluate. For instance TAptrn="TA-".}

\item{SIDpat}{Pattern of name which is sample ID (SID_BID_CID)=("-S0001_BES00001-14_2014234231"). Here SIDpat="-S". Can also be a vector.}

\item{printHistPlots}{Boolean of showing plots of the scores (number of matching alleles or LRs) for each comparisons.}
}
\description{
dnamatch2 is a function for doing large scale DNA database search between trace samples (also mixture profiles) and between trace samples and reference samples.
}
\details{
dnamatch2 automatically imports DNA profiles from genemapper-format and feeds it into a structure for doing effective comparison matching. Before the comparison matches are carried out, every trace samples are optionally filtered: Alleles with peak heights below a speficied threshold (threshHeight) or alleles in stutter position having peak heights below a specified threshold (threshStutt) are removed. The comparison match algorithm first counts the number of alleles of the reference profiles which are included in each of the trace samples. Then a candidate list of matching situations is created, whom satisfy being over a given treshold (threshMAC). Here all candidate matches having the same CID (case ID) can be optionally be removed. If wanted, also candidate matches which have a timedifference (based on last edited file dates) outside a specified time difference (timediff) can be removed. The second part of the comparison match algorithm first estimates the number of contributors of the remaining trace samples in the candidate list using the likelihood function of the quatliative model (likEvid in forensim), based on the AIC. After, the Likelihood Ratio for all comparisons (between all reference and trace samples) are calculated based on the same qualitative model, and an updated candidate list is created by considering all comparisons with LR greater than a specified threshold (threshLR[1]). Last, the Likelihood Ratio for all remaining comparisons in the candidate list are calculated based on the same quantitative model (likEvidMLE in euroformix). The final match candidates are the comparisons which has a LR greater than a specified threshold (threshLR[2]). Finally, detailed results of these match candidates are automatically printed to files.

To search between trace samples, a major component of each trace sample is extracted (based on relative peak heights): The allele with largest peak height belongs to the major. If the second largest allele has ratio (relative to the largest allele) above a specified threshold (threshMaj), then second allele is part of the major profile as well. If the relative peak height between second and third largest allele has ratio greater than this threshold, no profile for the major is assigned.

It is optional whether between trace samples matches should be considered (boolean betweensamples). If this is turned off, the speed is drastically increased (because of many less comparisons).

Encoding idea: References coded with primenumbers, stains alleles/heights coded with original values but collapsed(/) strings per marker. 

Note 1: Trace samples must have this format: "SID_BID_CID", seperated by underscores. The trace name must be unique and consist of SID=Sample ID, BID=Unique ID for all traces, CID=Case ID. Hence the IDs must not contain underscores themselves. 
Note 2: The names of the files containing trace samples must start with pattern specified by variable TAptrn. The variable SIDpat is required to recognize specific sample type (i.e. it can be used as a filter).  
Note 3: Marker names of Amelogen must contain the name "AMEL" (not case-sensitive).

Timestamp="YY-MM-DD-HH-MM-SS" is used as identification key of when the dnamatch2 function was run.
Match results are stored as a table in matchfile.csv with the timestamp. Matches which are earlier found in matchfile.csv will not be stored again. The column "Checked" can be used for comments.
More details about a given dnamatch2 run are stored in the folder "session" with corresponding timestamp.
}
\references{
dnamatch2: An open source software to carry out large scale database searches of mixtures using qualitative and quantitative models.
}
\author{
Oyvind Bleka <oyvble.at.hotmail.com>
}
