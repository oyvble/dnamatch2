#' @title dnamatch2
#' @author Oyvind Bleka <oyvble.at.hotmail.com>
#' @description dnamatch2 is a function for doing large scale DNA database search between trace samples (also mixture profiles) and between trace samples and reference samples. 
#' @details dnamatch2 automatically imports DNA profiles from genemapper-format and feeds it into a structure for doing effective comparison matching. Before the comparison matches are carried out, every trace samples are optionally filtered: Alleles with peak heights below a speficied threshold (threshHeight) or alleles in stutter position having peak heights below a specified threshold (threshStutt) are removed. The comparison match algorithm first counts the number of alleles of the reference profiles which are included in each of the trace samples. Then a candidate list of matching situations is created, whom satisfy being over a given treshold (threshMAC). Here all candidate matches having the same CID (case ID) can be optionally be removed. If wanted, also candidate matches which have a timedifference (based on last edited file dates) outside a specified time difference (timediff) can be removed. The second part of the comparison match algorithm first estimates the number of contributors of the remaining trace samples in the candidate list using the likelihood function of the quatliative model (likEvid in forensim), based on the AIC. After, the Likelihood Ratio for all comparisons (between all reference and trace samples) are calculated based on the same qualitative model, and an updated candidate list is created by considering all comparisons with LR greater than a specified threshold (threshLR[1]). Last, the Likelihood Ratio for all remaining comparisons in the candidate list are calculated based on the same quantitative model (likEvidMLE in euroformix). The final match candidates are the comparisons which has a LR greater than a specified threshold (threshLR[2]). Finally, detailed results of these match candidates are automatically printed to files.
#' 
#' To search between trace samples, a major component of each trace sample is extracted (based on relative peak heights): The allele with largest peak height belongs to the major. If the second largest allele has ratio (relative to the largest allele) above a specified threshold (threshMaj), then second allele is part of the major profile as well. If the relative peak height between second and third largest allele has ratio greater than this threshold, no profile for the major is assigned.
#' 
#' It is optional whether between trace samples matches should be considered (boolean betweensamples). If this is turned off, the speed is drastically increased (because of many less comparisons).
#' 
#' Encoding idea: References coded with primenumbers, stains alleles/heights coded with original values but collapsed(/) strings per marker. 
#' 
#' Note 1: Trace samples must have this format: "SID_BID_CID", seperated by underscores. The trace name must be unique and consist of SID=Sample ID, BID=Unique ID for all traces, CID=Case ID. Hence the IDs must not contain underscores themselves. 
#' Note 2: The names of the files containing trace samples must start with pattern specified by variable TAptrn. The variable SIDpat is required to recognize specific sample type (i.e. it can be used as a filter).  
#' Note 3: Marker names of Amelogen must contain the name "AMEL" (not case-sensitive).
#' 
#' Timestamp="YY-MM-DD-HH-MM-SS" is used as identification key of when the dnamatch2 function was run.
#' Match results are stored as a table in matchfile.csv with the timestamp. Matches which are earlier found in matchfile.csv will not be stored again. The column "Checked" can be used for comments.
#' More details about a given dnamatch2 run are stored in the folder "session" with corresponding timestamp.
#'
#' @param fn A folder with stain files (possible a vector). Full directory must be given. 
#' @param freqfile A file containing population frequencies for alleles. Full directory must be given.
#' @param reffold A folder with stored personal references. Default is no references. Full directory must be given.
#' @param sameCID Boolean whether matches within same case ID number should be allowed.
#' @param betweensamples A boolean of whether between samples are searched. Default is TRUE.
#' @param Thist Number of days back in time for a TA-file to be imported (stains). 
#' @param threshMAC Threshold for a allele match. A ratio [0,1]
#' @param threshLR Threshold for a match. Can be a vector (qual,quan).
#' @param threshHeight Acceptable peak height (rfu) in stains.
#' @param threshStutt Acceptable stutter ratio in stains (relative peak heights) .
#' @param threshMaj If second largest allele has ratio (relative to the largest allele) above this threshold, then second allele is part of the major profile (used for extracting major from mixture). If the relative peak height between second and third largest allele has ratio greater than this threshold, no major is assigned.
#' @param minLocStain Number of minimum loci in stain profiles in order to be evaluated.
#' @param minLocMaj Number of minimum loci which are required to be considered for a major component (extracted from each stain profile).
#' @param pC Assumed drop-in rate per marker (parameter in the LR models).
#' @param lambda Assumed hyperparameter lambda for peak height drop-in model (parameter in quantiative model).
#' @param kit The shortname of the kit used for the samples. This allows for taking into account degradation. Use getKit function in euroformix R-package.
#' @param N Database size used for estimating population frequencies. Used to assign new allele frequenceis. Default is 5000.
#' @param searchtime An object with format as returned from Sys.time(). Default is Sys.time().
#' @param SIDvec A vector with Sample-ID numbers which are considered in the search (i.e. a search filter).
#' @param TAvec A vector with TA-files, which are considered in the search (i.e. a search filter).
#' @param CIDvec A vector with Case-ID numbers which are considered in the search (i.e. a search filter).
#' @param timediff Timedifference (in days) allowed between matching reference and target. Default is NULL (not used).
#' @param TAptrn Filename structure of stain files to evaluate. For instance TAptrn="TA-".
#' @param SIDpat Pattern of name which is sample ID (SID_BID_CID)=("-S0001_BES00001-14_2014234231"). Here SIDpat="-S". Can also be a vector.
#' @param printHistPlots Boolean of showing plots of the scores (number of matching alleles or LRs) for each comparisons. 
#' @references dnamatch2: An open source software to carry out large scale database searches of mixtures using qualitative and quantitative models.
#' @export

#library(roxygen2);roxygenize("dnamatch2")

dnamatch2 <- function (fn, freqfile, reffold = NULL, sameCID = FALSE ,betweensamples=TRUE, 
    Thist = Inf, threshMAC=0.75, threshLR = c(10,100), threshHeight = 200,threshStutt = 0.1, threshMaj = 0.6, 
    minLocStain=3,minLocMaj=3, pC = 0.05, lambda=0.01,kit="ESX17",
    N =5000, searchtime=Sys.time(),SIDvec=NULL,TAvec=NULL,CIDvec=NULL,timediff=NULL,TAptrn=NULL, SIDpat=NULL, BIDpat=NULL,printHistPlots=FALSE) 
{

    library(forensim) #requires at least v4.3
    library(euroformix) #requires at least v0.6.4
    newf0 <- 5/(2 * N) #allele-frequence for new alleles

	###############
	#HELPFUNCTIONS#
	###############
    tableReader = function(filename) { #Robust function for reading tables:
        tab <- read.table(filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        tryCatch({
            if (ncol(tab) == 1) tab <- read.table(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE)
        }, error = function(e) e)
        tryCatch({
            if (ncol(tab) == 1) tab <- read.table(filename, header = TRUE, sep = ";", stringsAsFactors = FALSE)
        }, error = function(e) e)
        if (ncol(tab) == 1)  tab <- read.table(filename, header = TRUE, sep = ";", stringsAsFactors = FALSE)
        return(tab)
    }
    strsplit2 <- function(x, spl) {
        if (nchar(x) == 0)  return("")
        txt <- x
        for (j in 1:length(spl)) {
            txt <- unlist(strsplit(txt, split = spl[j]))
        }
        return(txt)
    }
    makematrix <- function(x) { #not necessary if you use X[1,,drop=F]
        if (is.null(dim(x))) x <- t(x)
        return(x)
    }
    prim = as.integer(c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 
        37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 
        101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 
        157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 
        223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 
        277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 
        349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 
        419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 
        479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 
        563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 
        619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 
        691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 
        769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 
        853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 
        929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 
        1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 
        1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 
        1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 
        1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 
        1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 
        1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 
        1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 
        1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 
        1523, 1531, 1543, 1549)) #max number of alleles (per marker) is 244

    #1) Configurate Kit:
    tab = tableReader(freqfile)
    locs <- toupper(colnames(tab[, -1]))
    popFreq <- popFreqP <- list()
    for (loc in locs) {
        freqs <- tab[, which(loc == locs) + 1]
        popFreq[[loc]] <- popFreqP[[loc]] <- tab[!is.na(freqs), which(loc == locs) + 1]
        names(popFreq[[loc]]) <- tab[!is.na(freqs), 1]
        names(popFreqP[[loc]]) <- prim[1:length(popFreq[[loc]])]
    }
    names(popFreq) <- names(popFreqP) <- locs
	
    DBcolN <- c("ID","SID","CID","TID","Time") #column name of DBrefN and DBstainN: (id,sampleID,"CIDnr","TA-file","Time")
    #2) IMPORT References:
    DBref <- matrix(nrow = 0, ncol = length(locs)) #encoded matrix
    DBrefN <- numeric()  #corresponding names

    Rfiles <- NULL
    if (!is.null(reffold)) Rfiles <- list.files(reffold)
    for (Rfile in Rfiles) {
#Rfile=Rfiles[1]
        X = tableReader(paste0(reffold, "/", Rfile))
        cn = colnames(X)
        sind = grep("ID", toupper(cn)) #sample col-ind
        if (length(sind) == 0)  sind = grep("SAMPLE", toupper(cn))
        if (length(sind) > 1)   stop("More than one SID was find in the file")
        sn = unique(as.character(X[, sind]))
        if (length(sn) == 1) {  #only one sample was found
            sn <- strsplit(Rfile, "\\.")[[1]][1] #get samplename of filename
            X[, sind] <- sn #update sample name within file!
        }
        lind = grep("marker", tolower(cn)) #locus col-ind
        X[, lind] <- toupper(X[, lind]) #make all loci upper-case
        ln = unique(X[, lind]) #locus names
        if (length(grep("AM", ln)) > 0) ln <- ln[-grep("AM", ln)]  #remove AMEL
        if (!all(ln %in% locs)) {
            print(paste0("Reference-File '", Rfile, "' had another kit-format than ESX17!"))
            next
        }
        Aind = grep("allele", tolower(cn))[1:2] #allele col-ind. Only consider the two first
        ishom <- X[, Aind[2]] == "" | is.na(X[, Aind[2]])
        X[ishom, Aind[2]] <- X[ishom, Aind[1]] #add homozygote to both

	  #Introduce fast encoding:
        smallREF <- matrix(1, ncol = length(locs), nrow = length(sn))
        for (loc in locs) {
            locind <- which(X[, lind] == loc) #get index of loci
            tmp <- popFreq[[loc]]
            tmpP <- popFreqP[[loc]]
            newA1 <- X[locind, Aind[1]][!X[locind, Aind[1]] %in% names(popFreq[[loc]])] #get new Alleles
            newA2 <- X[locind, Aind[2]][!X[locind, Aind[2]] %in% names(popFreq[[loc]])] #get new Alleles
            newA <- unique(c(newA1, newA2))
            if (length(newA) > 0) {
                tmp <- popFreq[[loc]]
                tmpP <- popFreqP[[loc]]
                popFreqP[[loc]] <- popFreq[[loc]] <- c(popFreq[[loc]], rep(newf0, length(newA))) #insert new freqs.
                names(popFreq[[loc]]) <- c(names(tmp), newA) #insert new a-names
                names(popFreqP[[loc]]) <- c(names(tmpP), prim[(length(tmp) + 1):length(popFreq[[loc]])]) #insert new p-names
            }
            Pname <- as.integer(names(popFreqP[[loc]]))  #get primenumbers
            Aname <- names(popFreq[[loc]])  #get allelenames
            for (an in Aname) { #for each alleles in population
                for (ai in Aind) {
                  aind <- X[locind, ai] == an #get rows which has corresponding allele
                  srow <- which(sn %in% X[locind, sind]) #get belonging rows in smallREF
                  smallREF[srow[aind], which(locs == loc)] <- smallREF[srow[aind], which(locs == loc)] * Pname[which(Aname == an)]
                }
            }
        }
        smallREF[smallREF == 1] = NA #insert missing loci as NA
        DBref <- rbind(DBref, smallREF)
        DBrefN <- rbind(DBrefN, cbind(sn,"0","0","0","0") )
    } #end for each file
    if (length(DBrefN) == 0) {
       print("NOTE: No reference samples were imported!")
    } else {
       colnames(DBrefN) <- DBcolN 
       #Check for duplicated R-ID:
       if (length(unique(DBrefN[,DBcolN=="ID"])) != nrow(DBrefN))  stop("Some imported references had same R-ID")
    }
    colnames(DBref) <- locs

    #2) IMPORT Stains:
    DBstainA <- DBstainH <- DBstainN <- numeric() #matrices to store all stains, DBstainN stores SID,BID,KID,TAnames
    for(ff in fn) {
#ff <- fn[1]
	 TAfiles <- list.files(ff, pattern = TAptrn) #Uses TA-pattern in files.
	for (TAfile in TAfiles) {	
#TAfile = TAfiles[1]
		if(!is.null(TAvec) && !any(grepl(TAfile,TAvec))) next #skip TAfile if not specified    
		fname <- paste0(ff, "/", TAfile)
		X = tableReader(fname)
		cn = colnames(X)
		sind = grep("sample", tolower(cn)) #sample col-ind
		if (length(sind) > 1) sind <- sind[grep("name", tolower(cn)[sind])] #if multiple sind
		sn = unique(X[, sind]) #sample names
		if (length(sn) == 0)   next  #this file did not contain valid samples!
		lind = grep("marker", tolower(cn))  #locus col-ind
		X[, lind] <- toupper(X[, lind]) #make all loci upper-case
		ln = unique(X[, lind]) #locus names:
		ln <- ln[-grep("AMEL", ln)] #remove AMEL
		if (!all(ln %in% locs)) { #check that marker names are same as in the freq-file
			print(paste0("Mixture-File '", TAfile, "' had another kit-format!"))
			next
		}
		Aind = grep("allele", tolower(cn)) #allele col-ind
		Hind = grep("height", tolower(cn)) #height col-ind
		if (length(Aind) == 0 | length(Hind) == 0)  next  #if no allele or peak height info
			
		#checking stain-id pattern:
            sidind <- tind <- 1
            if(!is.null(SIDpat)) {
	 		tmpS <- NULL
			for (pat in SIDpat) {
				tmp <- sapply(strsplit(sn, pat), function(x) x[2])
				if (all(!is.na(tmp))) {
					tmpS <- tmp
					break
				}
			}
			if (is.null(tmpS)) next
			sidind <- which(SIDpat==pat) #get index of SID pattern
		}

		#get time-info from file:
		Tfile <- file.info(fname)$mtime #time when file was last modified!
		Fdate <- format(Tfile)
		Tdiff <- difftime(searchtime, Tfile, units = "days")
		Tdiff <- as.integer(strsplit(format(Tdiff, digits = 1)," ")[[1]][1])
		if(length(Thist)==length(SIDpat)) tind <- sidind  #Thist index to use
		if (Tdiff > Thist[tind]) next #don't import if file was modified more than Thist days ago	

		SID <- sapply(strsplit(sn, "_"), function(x) x[1]) #extract SID
		BID <- sapply(strsplit(sn, "_"), function(x) x[2]) #extract BID
		KID <- sapply(strsplit(sn, "_"), function(x) x[3]) #extract CID
		tmp <- strsplit2(TAfile, c(TAptrn, "\\."))[1] #TA-id
		TID <- paste0(TAptrn, tmp) #TA-id
		for (ss in 1:length(sn)) { #for each combination of SID and BID (unique stains)
			if(!is.null(SIDvec) && !any(grepl(SID[ss],SIDvec) )) next #skip sample if not specified
			if(!is.null(CIDvec) && !any(grepl(KID[ss],CIDvec) )) next #skip sample if not specified
			subX <- X[X[, sind] == sn[ss], ] #get submatrix
			SSvec <- rep(NA, length(locs))  #place to store allele-info into MAJOR/single source-matrix
	
			Avec <- Hvec <- rep(NA, length(locs))   #place to store allele/height-info into evidence-matrix
			for (ii in 1:nrow(subX)) { #for each marker-row
				loc <- toupper(subX[ii, lind])
				locind <- grep(loc, locs) #find correct locus
				if (length(locind) == 0)  next
				Ainfo <- as.character(subX[ii, Aind]) #extract allele info
				Hinfo <- as.numeric(subX[ii, Hind]) #extract height info
				usedInd <- !(is.na(Hinfo) | Ainfo == "" | Ainfo == "NA" | Ainfo == "OL")
				Ainfo <- Ainfo[usedInd]
				Hinfo <- Hinfo[usedInd]

				#(1) Peak height above threshold
				keep <- Hinfo >= threshHeight #required minumum peak height
				Ainfo <- Ainfo[keep]
				Hinfo <- Hinfo[keep]
				if (length(Ainfo) == 0) next
				  
				#(2) Stutter-filter
				AllsNum <- as.numeric(Ainfo) #convert to numbers
				stuttindL <- which(as.character(AllsNum) %in% as.character(AllsNum - 1)) #stutter-alleles one low bp
				stuttindH <- which(as.character(AllsNum - 1) %in% as.character(AllsNum)) #stutter-alleles one high bp
				stuttR <- Hinfo[stuttindL]/Hinfo[stuttindH] #stutter-ratio is comparing observed peak heights
				remove <- stuttindL[stuttR < threshStutt]  #alleles to remove
				if (length(remove) > 0) {
					Ainfo <- Ainfo[-remove]
					Hinfo <- Hinfo[-remove]
				}
				
				#Insert new alleles
				Anew <- unique(Ainfo[!Ainfo %in% names(popFreq[[loc]])])
				if (length(Anew) > 0) {
					tmp <- popFreq[[loc]]
					tmpP <- popFreqP[[loc]]
					popFreqP[[loc]] <- popFreq[[loc]] <- c(popFreq[[loc]], rep(newf0, length(Anew)))
					names(popFreq[[loc]]) <- c(names(tmp), Anew)
					names(popFreqP[[loc]]) <- c(names(tmpP), prim[(length(tmp) + 1):length(popFreq[[loc]])])
				}
				Pname <- as.integer(names(popFreqP[[loc]])) #remember to update Pname (used for SS-vec)
				 
                        #Collapse alleles and heights and save in vector 
				Avec[locind] <- paste0(Ainfo,collapse="/")
				Hvec[locind] <- paste0(Hinfo,collapse="/")
	
				#Extract major profiles: Put into an own matrix:
				if(betweensamples) {
					ord <- order(Hinfo, decreasing = TRUE)
					Ainfo <- Ainfo[ord]
					Hinfo <- Hinfo[ord]
					SS <- Ainfo[1] #largest allele always included
					nA <- length(Ainfo)
					if (nA > 1 && Hinfo[2]/Hinfo[1] > threshMaj) { #if more than 1 allele and ratios of two largest suffice
						SS <- c(SS, Ainfo[2]) #add second allele
						if (nA > 2 && Hinfo[3]/Hinfo[2] > threshMaj) { #if more than 2 major alleles and also 3. largest suffice
							SS <- NA #couldn't deside anything for marker
						}
					}
					if (!any(is.na(SS))) {
                                tmp <- prod(Pname[names(popFreq[[loc]]) %in% SS])
                                if(length(SS)==1) tmp <- tmp*tmp #homozygous multiplied twice
					  SSvec[locind] <- tmp #store SS-info
                 			}
				} #end if between samples
 			} #end for each markers

                  if( sum(!is.na(Avec)) < minLocStain ) { #if too few accepted markers (below a given threshold)
				print(paste0("Sample ",sn[ss]," was not included: Num.mark=",sum(!is.na(Avec)))) #print number of markers
				next
			}
                  stainN <- cbind(sn[ss],SID[ss],KID[ss], TID, Fdate) #stain names (include date for file also)
			DBstainN <- rbind(DBstainN, stainN ) 
			DBstainA <- rbind(DBstainA, Avec ) 
			DBstainH <- rbind(DBstainH, Hvec ) 

			if(betweensamples && sum(!is.na(SSvec))>=minLocMaj ) {
				DBref <- rbind(DBref, SSvec ) #add single source to DBref
				DBrefN <- rbind(DBrefN , stainN) #add full sample name (as in v1.0)
			}
		}  #end for each sample
		if (which(TAfiles == TAfile)%%10 == 0)  print(paste0(round(which(TAfiles == TAfile)/length(TAfiles) * 100), "% import complete"))
	 } #end for each TAfiles
    } #end for each fn path
    if (length(DBstainN) == 0)  stop("No evidence was imported! Program stops.")
    colnames(DBstainN) <- DBcolN  
    colnames(DBstainA) <- colnames(DBstainH) <- locs

    #Check for duplicates of tracenamesnames:
    dupind <- duplicated(DBstainN[, DBcolN=="ID"])
    if (sum(dupind) > 0) {
        print(DBstainN[dupind,])
        warning(paste0("It was ", sum(dupind), " duplicated stain samples. Program continuous..."))
        DBstainN <- DBstainN[!dupind, ,drop=F]
        DBstainA <- DBstainA[!dupind, ,drop=F]
        DBstainH <- DBstainH[!dupind, ,drop=F]
   }

	
############################################################################################################################################
######################################COMPARISON#########################################################################
#############################################################################################################################################

######################################################################
#Part 1: Compare number of alleles in references included in evidence#
######################################################################
#Create full MAC matrix: all combinations of rows in DBstainN and rows in DBrefN:
#DBrefLocs <- rowSums(!is.na(DBref))
nS <- nrow(DBstainN)
nR <- nrow(DBrefN)
nL <- ncol(DBstainA) #number of loci
print(paste0("Calculating MAC for all ",nS*nR," comparisons: All refs against all stains"))
bigMAC <- rep(0,nS*nR) #keep MAC in a vector (sample1-ref1,sample1-ref2,...,sample2-ref1 etc.)
#all(colnames(DBstainA)==colnames(DBref))

systime <- system.time( {
for(ss in 1:nS) { #for each sample: Limited in how the samples are looking
#ss=1
 bigInd <-  nR*(ss-1) + 1:nR  #index in bigMAC matrix
 macS <- nLocs <- rep(0,nR) #make vector for all references 
 for(ll in 1:ncol(DBstainA)) { #for each locus: Vectorizing as much as possible!
#ll=2
  if(is.na(DBstainA[ss,ll]))  next #skip if no data
  sttmp <- unlist(strsplit(DBstainA[ss,ll],"/")) #get alleles
  sttmpP <- as.integer(names(popFreqP[[ll]][match( sttmp , names(popFreq[[ll]]) )])) #get prime numbs
  isna <- is.na(DBref[,ll]) #get which refs are non-zero
  numWithin <- rep(0,sum(!isna)) #number of unique alleles of reference that are in stain
  for(aa in sttmpP) numWithin <- numWithin + as.integer(DBref[!isna,ll]%%aa==0) #for each alleles in stain
  isHom <- sqrt(DBref[!isna,ll])
  isHom <- abs(round(isHom) - isHom)<1e-6 #which are homozygous
  macS[!isna] <- macS[!isna] + numWithin #add number of matching alleles
  macS[!isna][isHom & numWithin==1] <- macS[!isna][isHom & numWithin==1] + 1 #add homozygous twice IF it was matching with one
  nLocs[!isna] <- nLocs[!isna] + 1  #add locus  
 } #end for each locus
 bigMAC[bigInd] <- macS/(2*nLocs)  #divide number of matching-alleles in refernce by maximum possible
#hist(macS)
 #insert MAC to long vector
 if (ss%%100 == 0)  print(paste0(round(ss/nS* 100), "% MAC calculation complete")) 
} #end number of samples
})[3]
print(paste0("Calculating MAC scores took ",ceiling(systime), " seconds"))
#End compare with MAC
if(printHistPlots && length(bigMAC)>0) hist(bigMAC)
#max(bigMAC)

###########
#FILTER 1:# 
###########
keepInd <- which(bigMAC>=threshMAC)
#bigMAC[keepInd]
#(1:length(bigMAC)-1)%%nR + 1
#floor( (1:length(bigMAC)-1)/nR ) + 1

Ctab <- cbind(  floor((keepInd-1)/nR) + 1, (keepInd-1)%%nR + 1 , bigMAC[keepInd]) #convert back indices
colnames(Ctab) <- c("stain","ref","score") #note the order: stainID first
#nR*(Ctab[,1]-1) + Ctab[,2] #check
#cbind(DBstainN[Ctab[,1],4],DBrefN[Ctab[,2]]) #
#Ctab <- Ctab[rowSums(is.na(DBref[Ctab[,2],]))>0,] #consider these only


#FILTER 1b: remove for other reasons:
#(1) Remove because it was the same stain:
sameSID <- DBrefN[ Ctab[,2],DBcolN=="SID"]==DBstainN[ Ctab[,1],DBcolN=="SID"]  #check for mathces which have same sampleID
Ctab <- Ctab[!sameSID,,drop=F] #remove those who are the same
    
#(2) Remove because it was the same CID
if (!sameCID) {
 sameCID2 <- DBrefN[ Ctab[,2],DBcolN=="CID"]==DBstainN[ Ctab[,1],DBcolN=="CID"]  #check for mathces which have same sampleID
 Ctab <- Ctab[!sameCID2,,drop=F] #remove those who are the same
}

#(3) Remove because it was outside time difference
if (betweensamples && !is.null(timediff)) {
  isREF <-  DBrefN[Ctab[,2],DBcolN=="Time"]=="0" #get samples which are ref
  if(!all(isREF)) {
   diff <- abs(difftime(DBrefN[Ctab[!isREF,2],DBcolN=="Time"], DBstainN[ Ctab[!isREF,1],DBcolN=="Time"], units = "days"))
   isREF[!isREF] <- diff<=timediff #if ref-profile or inside time
   Ctab <- Ctab[isREF,,drop=F] #keep only ref-profiles or inside time
  }
}
nK <- nrow(Ctab)
print(paste0("Number of comparisons satisfying (after filters) threshMAC=",threshMAC,": ", nK))
if(nK==0) stop("There were no more comparisons to do. Program stops!")

#####################################################
#Part 2: Calculate LRqual for relevant combinations:#
#####################################################
unStain <- unique(Ctab[,1]) #unique stain profiles: 
nU <- length(unStain) #number of unique
print(paste0("Estimating num. contr. for ", nU," stains"))
print(paste0("Calculating LRqual for ", nK," combinations"))

#2.1) Estimate number of contr. using qual model 
logLik <- function(pDvec,joint=TRUE) { 
          val <- numeric()
          for(loc in names(data$popFreq)) { #for each locus
           Ei <- NULL #get evidence
           for(s in 1:length(data$samples)) { #fix samples
            if(s>1) Ei <- c(Ei,0) #seperate with 0  
            adata <- data$samples[[s]][[loc]]$adata
            if(length(adata)==0) adata=0 #is empty
            Ei <- c(Ei,adata)
           } 
           Ti <- unlist(data$refData[[loc]])
           if( length(Ti)==0 ) Ti=NULL
           val <- c(val , log( likEvid( Ei,T=Ti,V=NULL,x=data$nU[loc==names(data$popFreq)],theta=0, prDHet=pDvec, prDHom=pDvec^2, prC=pC, freq=data$popFreq[[loc]]) ) )
         } #end for each markers
         if(joint) val <- sum(val)
         return(val)
}
negloglik <- function(pD) {
	return( -logLik( pDvec=rep(1/(1+exp(-pD)),K) )) #K is outer var.
}

Kqual <- log10LRqual <- rep(0,nK) #for storing num. contr and logLd=log P(E|Hd)
systime <- system.time( {
for(ss in unStain) { #for each unique stain we estimate number of contr and calculate the LR for all the references
# ss = unStain[1]
 dat <- DBstainA[ss,]
 locUseS <- !is.na(dat) | TRUE #loci to consider for sample: Assumes kit to popFreq. Important to get correct image of whole sample
 locs <- names(popFreq)[locUseS] #update locs to use in all calculations
 #if(any(is.na(dat))) stop("sasd")
 samples <- list(lapply(dat[locUseS], function(x) { #create sample data list
  a0 <- unlist(strsplit(x,"/"))
  if(all(is.na(a0))) a0 <- numeric()
  list(adata=a0)
 })) #strsplit(dat[locUseS],"/")
 data <- Qassignate(samples, popFreq[locUseS],incS=FALSE)
 data$samples <- samples
 K <- ceiling(max(sapply(dat,function(x) length( unlist(strsplit(x,"/")) )))/2) #lowest number of contr
 
 done <- FALSE
 hdval <- matrix(,ncol=4,nrow=0) #columns: #contr,,dropprob,loglik_max,loglik_max-K
 while(!done) { #model selection for LRmix: Find MLE under hd
      if(K>3) break #maximum found (at most 3 contributors)
      data$nU <- rep(K,length(locs)) #equal for all loci
      foo <- nlm(Vectorize(negloglik),log(0.1/(1-0.1)))
      pdhat <- 1/(1+exp(-foo$est)) #maximum likelihood dropout estimate 
      hdval <- rbind(hdval,c(K,pdhat,-foo$min, -foo$min - K))
      #print(hdval)
      nval <- nrow(hdval)
      if(nval>1 && as.numeric(hdval[nval,4])<as.numeric(hdval[nval-1,4]) ) {
       done <- TRUE
      } else {
       K <- K + 1  
      }
  }
  indAIC <- which.max(hdval[,4]) 
  Khat <-  hdval[indAIC,1] #estimated contr.
  pDhat <-  hdval[indAIC,2] #estimated dropout
  loghd <- hdval[indAIC,3] #Problem: Not always global since ref may miss some markers.
  Kqual[ ss == Ctab[,1] ] <- Khat #estimated number used in quanLR directly (speeds up)
  K <- Khat #NB: K is necessary variable in negloglik function

  #Calc. Hp to get LR:
  whatR <- which(ss == Ctab[,1]) #what ref ind to consider
  subRef <- DBref[ Ctab[whatR,2], ,drop=F] #get references

  for(r in 1:length(whatR) ) { #for each  reference to compare with selected 
   data$nU <- rep(Khat - 1,length(locs)) #reduce number of unknown with one under Hp! Consider as vector, values may change

   #fix refData:
   subRef2 <- subRef[r,]
   refD <- list()
   for(loc in locs ) {
      pri <- subRef2[loc==names(popFreq)] #get prime number
      if(is.na(pri)) {
        data$nU[which(loc==locs)] <- data$nU[which(loc==locs)] + 1 #add an extra unknown if missing
#        next #skip if empty 
	   a0 <- numeric()
      } else {
       an <- names(popFreq[[loc]])
       anP <- as.integer( names(popFreqP[[loc]]) )
       a0 <- an[pri%%anP==0]
      }
      if(length(a0)==1) a0 <- rep(a0,2)
	refD[[loc]] <- list(ref1=a0)     
   }
   data2 <- Qassignate(samples, popFreq[locs],refD,incS=FALSE,incR=FALSE)
   data$refData <- data2$refData 
   data$popFreq <- data2$popFreq #needed to update relevant loci
   foo <- nlm(Vectorize(negloglik),log(0.1/(1-0.1)))    #Optimize under Hp:
   log10LRqual[whatR[r]] <- (-foo$min - loghd)/log(10) #get calculated LR
 } #end for each references given unique stain
 ii <- which(ss==unStain)  
 if (ii%%20 == 0)  print(paste0(round(ii/nU* 100), "% LR qual calculation complete")) 
} #end for each stains
#Kqual 
})[3]
print(paste0("Calculating LR (qual) took ",ceiling(systime), " seconds"))

if(printHistPlots && length(log10LRqual)>0) hist(log10LRqual)

###########
#FILTER 2:# 
###########
keepInd2 <- which(log10LRqual>=log10(threshLR[1]))
Ctab2 <- cbind(Ctab[keepInd2,,drop=F], 10^log10LRqual[keepInd2]) #include LRqual (on original scale)
nK2 <- nrow(Ctab2) #update nK (number of combinations)
print(paste0("Number of comparisons satisfying threshLRqual>",threshLR[1],": ", nK2))
if(nK2==0) stop("There were no more comparisons to do. Program stops!")

#####################################################
#Part 3: Calculate LRquan for relevant combinations:#
#####################################################
#Using estimated contr. in Kqual for quan LR:
#ALSO IT REQUIRES ALL MARKERS FOR REFERENCES: ADVANCED NOT IMPLEMENTED!
log10LRquan <- rep(0,nK2)
locs <- colnames(DBstainA) #loci to consider
nDone <- 2 #number of optimizations
unStain <- unique(Ctab2[,1]) #unique stain profiles
nU <- length(unStain) #number of unique stains

print(paste0("Calculating LRquan for ", nK2," combinations (",nU," unique samples)"))
#Kqual[Ctab[,1]%in%unStain]
#nrow(unique(DBstainA[Ctab2[,1],]))

#Notice: Empty markers important because of information about allele dropouts.
systime <- system.time( {
for(ss in unStain) { #for each unique stain we estimate number of contr.
# ss = unStain[1]
 datA <- DBstainA[ss,]
 datH <- DBstainH[ss,] #get peak heights

 samples <- list()
 for(loc in locs) {
  ind <- which(locs==loc)
  a0 <- unlist(strsplit(datA[ind],"/"))
  h0 <- as.numeric(unlist(strsplit(datH[ind],"/")))
  if(is.na(datA[ind])) a0 <- h0 <- numeric()
  samples[[loc]] <- list(adata=a0,hdata=h0)
 }
 samples <- list(samples)
 data <- Qassignate(samples, popFreq[locUseS],incS=FALSE) #don't include stutters
 nC <- Kqual[ ss == Ctab[,1] ][1] #number of contr as for qual model
 
 fitHd <- contLikMLE(nC,samples,data$popFreq,xi=0,prC=pC,lambda=lambda,nDone=nDone,threshT=threshHeight,kit=kit,verbose=FALSE)
 #fitHd$fit$thetahat2
 loghd <- fitHd$fit$loglik #used under all combination for this stain.

  #Calc. Hp to get LR:
  whatR <- which(ss == Ctab2[,1]) #what ref ind to consider
  subRef <- DBref[ Ctab2[whatR,2], ,drop=F] #get references
  for(r in 1:length(whatR) ) {
#r=1
   #fix refData:
   subRef2 <- subRef[r,]
   refD <- list()
   for(loc in locs) {
      pri <- subRef2[loc==locs] #get prime number
      if(is.na(pri)) {
	   a0 <- numeric() #if empty
      } else {
       an <- names(popFreq[[loc]])
       anP <- as.integer( names(popFreqP[[loc]]) )
       a0 <- an[pri%%anP==0]
       if(length(a0)==1) a0 <- rep(a0,2)
      }
	refD[[loc]] <- list(ref1=a0)     
   }
   #locUseR <- names(popFreq)%in%names(refD) #loci in ref.
   data2 <- Qassignate(samples, popFreq[locs],refD,incS=FALSE,incR=FALSE) #popFreq must be given with correct order?
   data$refData <- data2$refData 
   data$popFreq <- data2$popFreq #needed to update relevant loci

   #Optimize under Hp:
   fitHp <- contLikMLE(nC,samples,data$popFreq,data$refData,condOrder=1,xi=0,prC=pC,lambda=lambda,nDone=nDone,threshT=threshHeight,kit=kit,verbose=FALSE)

   #Optimized Values under Hd:
   log10LRquan[whatR[r]] <- (fitHp$fit$loglik  - loghd)/log(10) 
 } #end for each references given unique stain
 ii <- which(ss==unStain)  
 if (ii%%5 == 0)  print(paste0(round(ii/nU* 100), "% LR quan calculation complete")) 
} #end for each stains
})[3]
print(paste0("Calculating LR (quan) took ",ceiling(systime), " seconds"))
#hist(log10LRquan)
if(printHistPlots && length(log10LRquan)>0) hist(log10LRquan)

###########
#FILTER 3:# 
###########
threshQuan <- threshLR[1]
if(length(threshLR)>1) threshQuan <- threshLR[2]
keepInd3 <- which(log10LRquan>=log10(threshQuan))
Ctab3 <- cbind(Ctab2[keepInd3,,drop=F], 10^log10LRquan[keepInd3]) #include LRquan (on original scale)
nK3 <- nrow(Ctab3)
print(paste0("Number of comparisons satisfying threshLRquan>",threshQuan,": ", nK3))
if(nK3==0) stop("There were no more comparisons to do. Program stops!")


###################################################################################################################################
######################################LR CALCULATIONS DONE#########################################################################
###################################################################################################################################

  Ctab0 <- Ctab3 #Final tab to use
  Ctab0 <- Ctab0[order(Ctab0[,ncol(Ctab0)],decreasing=T),,drop=F] #order results

	nLstain <- rowSums(!is.na(DBstainA[Ctab0[,1],,drop=F]) & !is.na(DBref[Ctab0[,2],,drop=F])) #get number of non-zero loci in both stain and reference
	stamp <- format(searchtime, "%y-%m-%d-%H-%M-%S")  #timestamp: Year-month-day-hour-minute-second
	cn <- c("refTA","tarTA","refID","tarID","nLocs","MAC","LRqual","LRquan","Searchtime") #colnames
 
    #Create result matrix:
	resMat <- cbind( DBrefN[Ctab0[,2], DBcolN=="TID"], DBstainN[Ctab0[,1], DBcolN=="TID"], DBrefN[Ctab0[,2], DBcolN=="ID"], DBstainN[Ctab0[,1], DBcolN=="ID"],nLstain,Ctab0[,-(1:2),drop=F],stamp)
	colnames(resMat) <- cn

    #FILTER AWAY SYMMETRIC MATCHES 
    allSID <- cbind(resMat[,cn=="refID"],resMat[,cn=="tarID"] )
    indrev <- as.character(allSID[, 1]) > as.character(allSID[, 2])
    allSID[indrev, ] <- allSID[indrev, 2:1]  #sort columns
    isdup <- duplicated(allSID) #find (redudance+symmetrical) matches already found before
    resMat <- resMat[!isdup,,drop=F] #get only unique matches!
 
    #5) STORE MATCHES
    sfolder <- "sessions"
    if (length(grep(sfolder, list.dirs(full.names = FALSE, recursive = FALSE))) == 0)   dir.create(sfolder) #create folder if not existing
    #save.image(file = paste0(sfolder,"/matching_", stamp, ".Rdata")) #dont save its too large!
    write.table( resMat[ resMat[,cn=="refTA"]!="0",,drop=F] , file = paste0(sfolder,"/stainmatches_", stamp, ".csv"), row.names = FALSE, sep = ";") #store matches not personell:
    write.table( resMat[ resMat[,cn=="refTA"]=="0",,drop=F] , file = paste0(sfolder,"/personelmatches_", stamp, ".csv"), row.names = FALSE, sep = ";") #store matches only personell:

 
	matchingprofiles = function(resMat) {
        if (length(resMat) == 0) return()
        cn <- colnames(resMat)
        refind <- which("refID" == cn)
        tarind <- which("tarID" == cn)
        refindTA <- which("refTA" == cn)
        tarindTA <- which("tarTA" == cn)

        fulltab <- numeric()
        for (ss in 1:nrow(resMat)) {
            refID <- resMat[ss,refind]
            refA <- DBref[DBrefN[,DBcolN=="ID"]==refID,] #get reference data
            tarID <- resMat[ss,tarind]
            tarA <- DBstainA[DBstainN[,DBcolN=="ID"]==tarID,] #get reference data
            tarH <- DBstainH[DBstainN[,DBcolN=="ID"]==tarID,] #get reference data


            nL <- paste0("nLocs=",resMat[ss,cn=="nLocs"])
            mac <- paste0("MAC=", round( as.numeric(resMat[ss,cn=="MAC"]),2))
            lr1 <- paste0("LRqual=",signif(as.numeric(resMat[ss,cn=="LRqual"]),digits=5))
            lr2 <- paste0("LRquan=",signif(as.numeric(resMat[ss,cn=="LRquan"]),digits=5))
            score = paste0(c(mac,lr1,lr2),collapse=" - ")
            rowtab <- rbind(c(resMat[ss, refind],resMat[ss, tarind]),c(nL,score),c("----REFERENCE----","---------TARGET------"))
            for (loc in locs) {
                Pname <- as.integer(names(popFreqP[[loc]]))
		    Ei <- paste0(tarA[loc==locs]," - ",tarH[loc==locs])
                RPi <- as.integer(refA[loc == locs])
                Ri <- ifelse(!is.na(RPi),paste0(names(popFreq[[loc]])[RPi%%Pname ==  0], collapse = "/"),NA)
                if(is.na(Ri)) next #don't show if ref is missing
                if(!grepl("/",Ri)) Ri <- paste0(Ri,"/",Ri)
                rowtab <- rbind(rowtab, c(paste0(loc,": ",Ri), Ei)) #v1.4 updated
            }
            fulltab <- rbind(fulltab, rep(paste0("----------", ss, "----------"), 2))
            fulltab <- rbind(fulltab, rowtab)
            fulltab <- rbind(fulltab, rep(paste0("------------------------------"), 2))
        }
        return(fulltab)
    }
	#store details of all matches:
    write.table(matchingprofiles(resMat), file = paste0(sfolder,"/matchinfo_", stamp, ".csv"), col.names = FALSE,row.names = FALSE, sep = ";")

    #Number of matches:
    print(paste0("Number of matches=", nrow(resMat)))
	
    #THIS BLOCK: Loads matches in matchfile. Makes sure that new matches don't duplicate with old (not even the symmetry)
    matchfile <- paste0("matchfile.csv")
    if (!file.exists(matchfile)) { #create file if not found
        x <- matrix(, ncol = length(cn), nrow = 0)
        colnames(x) <- cn
        write.table(x, file = matchfile, row.names = FALSE, col.names = TRUE, sep = ";")
    }
    if(length(resMat) > 0) {
        matchlist2 <- read.table(file = matchfile, sep = ";", header = TRUE, stringsAsFactors = FALSE) 
        cn2 <- colnames(matchlist2)
        SID <- cbind(resMat[,cn=="refID"],resMat[,cn=="tarID"] )
        SID2 <- cbind(matchlist2[,cn2=="refID"],matchlist2[,cn2=="tarID"] )

        allSID <- rbind(SID2, SID) #combine SID from file and SID from new matches
        allSID <- as.matrix(allSID) #vectors in matrix is already converted to factors, necessary to make string agaain!
        indrev <- as.character(allSID[, 1]) > as.character(allSID[,2])
        allSID[indrev, ] <- allSID[indrev, 2:1] #sort columns
        isdup <- duplicated(allSID) #find (redudance+symmetrical) matches already found before
        resMat2 <- resMat[!isdup[(nrow(SID2) + 1):length(isdup)], ,drop=F] #get only NEW matches!
        print(paste0("Number of new matches=", nrow(resMat2)))

        if (length(resMat2) > 0) {
            write.table(resMat2, file = matchfile, row.names = FALSE, col.names = FALSE, sep = ";", append = TRUE)
            write.table(matchingprofiles(resMat), file = paste0(sfolder,"/matchinfoNEW_", stamp, ".csv"), row.names = FALSE, col.names = FALSE,  sep = ";")
        }
    }
} #end dnamatch
