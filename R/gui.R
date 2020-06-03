#' @title gui
#' @author Oyvind Bleka <oyvble.at.hotmail.com>
#' @description A simple GUI for using dnamatch2
#' @details The user can select settings in a GUI for running dnamatch2
#' @param envirfile A file to a saved environment of a project
#' @export

#envirfile=NULL
#source("K:\\Sensitivt\\Klinikk04\\RESP_sikker\\Beregningsverktøy\\dnamatch2\\dnamatch2_2.0.0\\gui.R")

gui= function(envirfile=NULL) {
 require(dnamatch2)
 #size of main window
 mwH <- 800
 mwW <- 1000

 #Required for GUI:
 require(gWidgetstcltk) #requires only gWidgets.

 #type of gwidgets-kit
 options(guiToolkit="tcltk")

 #version:
 version =  packageVersion("dnamatch2") #follows same version as package number

 #software name:
 softname <- paste0("dnamatch2 v",version)

 #Spacing between widgets
 spc <- 5
 emptyName = "none" #Text Indicate that nothing is selected
 longspace = "                                                                                 "
 shortspace = "                          "
 timestamp = "%y-%m-%d-%H-%M-%S" #Format of time stamp (YEAR-MONTH-DAY-HOUR-MONTH-SECOND)

 ###############
 #HELPFUNCTIONS#
 ###############

 #helpfunction to get environment and file name for different data types
 getEnvirFileNames = function(type="EVID") {
   if(type=="EVID") {
    envirvar = "evidfolds"
    fname = evidFile 
   } else if(type=="REF") {
    envirvar = "reffolds"
    fname = refFile 
   } else if(type=="SIDsel") {
    envirvar = "SIDsel"
    fname = SIDselFile 
   } else if(type=="BIDsel") {
    envirvar = "BIDsel"
    fname = BIDselFile 
   } else if(type=="CIDsel") {
    envirvar = "CIDsel"
    fname = CIDselFile 
   } else if(type=="SIDpat") {
    envirvar = "SIDpat"
    fname = SIDpatFile 
   } else if(type=="BIDpat") {
    envirvar = "BIDpat"
    fname = BIDpatFile 
   } else {
    envirvar = ""
    fname = ""
   }
   return(c(envirvar,fname ))
 }

 #helpfunction to get folder/IDs vector from list (taken from list envir list only)
  getFolds = function(type="EVID") { #get vector of folders from list
   envirvar = getEnvirFileNames(type)[1]
   return( unlist(get(envirvar,envir=mmTK)) ) #return vector of folds
  }

  #helpfunction to set folders/IDs to environment and file (both envir list and file is updated with new info)
  addFold = function(foldadd,type="EVID") {
   tmp =  getEnvirFileNames(type)
   envirvar = tmp[1]
   fname = tmp[2]
 
   #Update environment variable:
   X = get(envirvar,envir=mmTK) #get list
   X[[length(X)+1]] = foldadd #add file to envir list
   assign(envirvar,X,envir=mmTK) #Store setup values

   #Store to file (update file):
   X = unlist(X)
   write(X,file=fname )    #save to file in installation folder 
   return(X) #return vector of folderes
  }

  saveFolds = function(folds,type="EVID") {
   tmp =  getEnvirFileNames(type)
   envirvar = tmp[1]
   fname = tmp[2]

   #Store to file (update file):
   write(folds,file=fname )    #save to file in installation folder 

   #Update environment variable:
   assign(envirvar,as.list(folds),envir=mmTK) #Store setup values
  }

  #Helpfunction to store/load values in optList to file
  saveSetup = function(opt) { 
   write(unlist(opt),file=setupFile)
  }
  openSetup = function() { 
   #vars = names(optL)
   vars = c("freqfile","workdir","IDsep","sameCID","betweensamples","Thist","timediff","searchtime","threshMAC","threshLRqual","threshLRquan","threshHeight","threshStutt","threshMaj","minLocStain","minLocMaj","pC","lambda","kit","minFreq","printHistPlots","writeScores","maxKqual","maxKquan","matchfile","sessionfold","printGraph")
   vart = c("s","s","s","b","b","d","d","s","d","d","d","d","d","d","i","i","d","d","s","d","b","b","i","i","s","s","b") #variable types (s=string,b=boolean,d=double,i=integer)
   dat = scan(file=setupFile,what=character(),quiet=TRUE,sep="\n")
   opt = list() #init list
   for(i in 1:length(vars)) { #
      x = dat[i] #string is standard
      if(vart[i]=="b") {
        x = as.logical(x)
      } else if(vart[i]=="d") {
	   x = as.numeric(x)
      } else if(vart[i]=="i") {
	   x = as.integer(x)
      }
      opt[[vars[i]]] = x  #insert correct type of variable for each list element
   }
   return(opt)
  }

 #helpfunction which checks that at value is in interval of [0,1]
 checkProb = function(x,what) {
  if(x < 0 || x>1) {
   gWidgets::gmessage(message=paste0(what," must be specified in interval [0,1] "),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
  return(x)
 }
 checkPositive = function(x,what,strict=FALSE) {
  if(x < 0 ) {
   gWidgets::gmessage(message=paste0(what," cannot be a negative number"),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
  if(strict && x==0) {
   gWidgets::gmessage(message=paste0(what," cannot be zero"),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
  return(x)
 }
 checkPosInteger = function(x,what,strict=TRUE) {
  if(x < 1 || round(x)!=x ) {
   gWidgets::gmessage(message=paste0(what," must be a positive integer"),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
  if(strict && x==0) {
   gWidgets::gmessage(message=paste0(what," cannot be zero"),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
  return(x)
 }
 val = function(wid) { #helpfunction to get value from widget
  tmp = as.numeric(gWidgets::svalue(wid))
  if(is.na(tmp )) {
   gWidgets::gmessage(paste0("Text found where number was expected:\n",gWidgets::svalue(wid)))
   return()
  }
  return(tmp )
 }

 errorMessage = function(msg) {  #Helpfunction to throw error message to user + stop running 
  gWidgets::gmessage(msg,title="Error",icon="error")
  stop(msg)
 }

 NullIfEmpty = function(x) {
   if(length(x)==0) {
    return(NULL)
   } else {
    return(unlist(x)) #input argument in dnamatch2 must be vectors instead of lists
   }
 }

############################################################
#FUNCTION TO RUN ANALYSIS (running with specified settings)#
############################################################
  runAnalysis = function(h,...) {
    #Step 1: Read search setup from envir vars
    #Step 2: Run dnamatch2 with setup
    #Step 3: Future update: Make availble Post-processing

	evidList = get("evidfolds",envir=mmTK) #get evid search folders
	refList = NullIfEmpty(get("reffolds",envir=mmTK))
	SIDselList = NullIfEmpty(get("SIDsel",envir=mmTK))
	BIDselList = NullIfEmpty(get("BIDsel",envir=mmTK))
	CIDselList = NullIfEmpty(get("CIDsel",envir=mmTK))
	SIDpatList = NullIfEmpty(get("SIDpat",envir=mmTK))
	BIDpatList = NullIfEmpty(get("BIDpat",,envir=mmTK))
	opt = get("setup",envir=mmTK)  #receive settings from envir (preassigned or from file)

     #prechecks:
	if(opt$freqfile==emptyName)  errorMessage("The freqency file has not been specified.\nPlease select one!") 
     kitUse=opt$kit #get kit
     if(kitUse==emptyName) kitUse = NULL #Set back to NULL if none specified

  status <- dnamatch2(
	evidfold=evidList, 
     freqfile=opt$freqfile, 
	reffold=refList, 
	sameCID=opt$sameCID,
	betweensamples=opt$betweensamples, 
	Thist=opt$Thist, 
	threshMAC=opt$threshMAC, 
	threshLR=c(opt$threshLRqual,opt$threshLRquan), 
	threshHeight=opt$threshHeight,
	threshStutt=opt$threshStutt, 
	threshMaj=opt$threshMaj, 
	minLocStain=opt$minLocStain,
	minLocMaj=opt$minLocMaj, 
	pC=opt$pC, 
	lambda=opt$lambda,
	kit=kitUse,
	minFreq=opt$minFreq, 
	searchtime=as.POSIXct(opt$searchtime,format=timestamp), # Sys.time(). Searchtime specified by timestamp!
	SIDvec=SIDselList,
	BIDvec=BIDselList,
	CIDvec=CIDselList,
	timediff=opt$timediff,
	IDsep = opt$IDsep,
	BIDptrn=BIDpatList,
	SIDptrn=SIDpatList,
	printHistPlots=opt$printHistPlots,
	writeScores=opt$writeScores,
	maxK=c(opt$maxKqual,opt$maxKquan),
	matchfile=opt$matchfile,
	sessionfold=opt$sessionfold
   )
  
  if(!status) gWidgets::gmessage("Search completed:\nNo match candidates were found!")

   #Step 3: Show graph of match candidates
   if(status && require(igraph) && opt$printGraph) { 
	tab <- read.table(file=opt$matchfile,header=TRUE,sep=";",stringsAsFactors=FALSE)
	fra <- tab$refID
	til <- tab$tarID 
	rem <- duplicated(cbind(fra,til)) #indices to remove
	rels <- data.frame(from=fra[!rem],to=til[!rem],weight=sqrt(log10(as.numeric(tab$LRquan[!rem]))))
	gg <- graph.data.frame(rels,directed=FALSE)

     dev.new() #avoid overriding existing plots
     plot(gg,edge.width=E(gg)$weight,vertex.color="white",vertex.size=10,vertex.label.color="black",vertex.label.cex=0.8)
     op <- par(no.readonly = TRUE)
     par(op)
   }
 } #end run function

 #####################
 #create environment #
 #####################
 pgkPath <- path.package("dnamatch2", quiet = FALSE) # Get package path.
 .sep <- .Platform$file.sep # Platform dependent path separator. 

 #STORING INFO IN BOTH Files and Environment (easy user access)
 #File variable - envir variable:
 #evidFile (configEvid) - evidfolds
 #refFile (configRef) - reffolds
 #SIDselFile (configSIDsel) - SIDsel  #Filter on specific SampleID names. Prefix included.  
 #BIDselFile (configBIDsel) - BIDsel  #Filter on specific Batch/TA-files. Prefix included. 
 #CIDselFile (configCIDsel) - CIDsel  #Filter on specific Batch/TA-files. Prefix included. 
 #SIDpatFile (configSIDpat) - SIDpat  #Pattern before SampleID nr (prefix). Used to recognize types. 
 #BIDpatFile (configBIDpat) - BIDpat  #Pattern before BatchID nr (prefix). Used to recognize types. 
 #setupFile (config) - setup


 evidFile <- paste(pgkPath,"configEvid",sep=.sep) #Setting file for evidence folders
 refFile <- paste(pgkPath,"configRef",sep=.sep) #Setting file for reference folders
 SIDselFile <- paste(pgkPath,"configSIDsel",sep=.sep) #Setting file for evidence folders
 BIDselFile <- paste(pgkPath,"configBIDsel",sep=.sep) #Setting file for reference folders
 CIDselFile <- paste(pgkPath,"configCIDsel",sep=.sep) #Setting file for reference folders
 SIDpatFile <- paste(pgkPath,"configSIDpat",sep=.sep) #Setting file for evidence folders
 BIDpatFile <- paste(pgkPath,"configBIDpat",sep=.sep) #Setting file for reference folders
 setupFile <- paste(pgkPath,"config",sep=.sep) #Create a file with all settings (not lists)

 if(is.null(envirfile)) {
  mmTK = new.env( parent = emptyenv() ) #create new environment object (globalenv?)

  scanFileToList = function(file) { #helpfunction to convert data from file to list
    as.list( scan(file=file,what=character(),quiet=TRUE,sep="\n") )  #convert to list 
  }
  
  evidList <- refList <- SIDselList <- BIDselList <- CIDselList <- SIDpatList <- BIDpatList <- list() #This is default
  if(file.exists(evidFile)) evidList <- scanFileToList(evidFile)
  if(file.exists(refFile)) refList <- scanFileToList(refFile)
  if(file.exists(SIDselFile)) SIDselList <- scanFileToList(SIDselFile)
  if(file.exists(BIDselFile)) BIDselList <- scanFileToList(BIDselFile)
  if(file.exists(CIDselFile)) CIDselList <- scanFileToList(CIDselFile)
  if(file.exists(SIDpatFile)) SIDpatList <- scanFileToList(SIDpatFile)
  if(file.exists(BIDpatFile)) BIDpatList <- scanFileToList(BIDpatFile)

  #Default set (empty) of folds:
  assign("evidfolds",evidList,envir=mmTK)
  assign("reffolds",refList,envir=mmTK)
  assign("SIDsel",SIDselList,envir=mmTK)
  assign("BIDsel",BIDselList,envir=mmTK)
  assign("CIDsel",CIDselList,envir=mmTK)
  assign("SIDpat",SIDpatList,envir=mmTK)
  assign("BIDpat",BIDpatList,envir=mmTK)

  #Default Settings set if setupfile not found:
  if(file.exists(setupFile)) {
   opt <- openSetup() #get setup list
  } else { #SETUP VALUES:
   opt = list() #list of options
   opt$freqfile = emptyName #name of freq file (should be full path?). Use file selector 
   opt$workdir = getwd() #default is work directory

   #Patterns for recognizing SID,BID,CID: 
   opt$IDsep= "_" #sepator for SID,BID,CID (Default)

   #Search type:
   opt$sameCID = FALSE #Search within same Cases (CIDs)? 
   opt$betweensamples = TRUE #Search between samples? 

   #Time windows: 
   opt$Thist = Inf #The number of days back in search time (using date of files)
   opt$timediff = Inf #required number of timedifference between match candidates (Can also be NULL by default)
   opt$searchtime = format(Sys.time(),timestamp) #Set current search time as when opening GUI (a specific string format), can be modified

   #Search thresholds:
   opt$threshMAC = 0.75
   opt$threshLRqual = 10
   opt$threshLRquan = 100

   #Prefiltering settings
   opt$threshHeight = 50 #Detection threshold specified
   opt$threshStutt = 0.1 #10% stutter rate is threshold
   opt$threshMaj = 0.6 #P.H. balancy to be included as major allele

   #Data quality settings:
   opt$minLocStain = 3
   opt$minLocMaj = 3

   #Model settings (dropin/kit):
   opt$pC = 0.05
   opt$lambda = 0.01
   opt$kit =  emptyName #Requires kitname (can be NULL also), use euroformix::getKit()
 
   #Small option settings:
   opt$minFreq = 0.001 #minimum frequency used
   opt$printHistPlots = TRUE
   opt$writeScores = TRUE
   opt$maxKqual=4 #max number of contributors under QUAL
   opt$maxKquan=3 #max number of contributors under QUAN
   opt$matchfile="matchfile.csv" #file name for storing results
   opt$sessionfold ="session" #folder name for storing results

   #GUI option settings:
   opt$printGraph = TRUE #Should a graph tree of matches in the matchfile be shown?

  } #end if file not found
  assign("setup",opt,envir=mmTK) #Store setup values 
 } else {
  load(envirfile) #loading environment
 }
 #optL = get("setup",envir=mmTK)  #receive settings from envir (preassigned or from file)

###################################################################
###########################GUI#####################################
###################################################################

 #Menu bar file-lists:
 f_setwd = function(h,...) {
  dirsel = gWidgets::gfile(text="Select folder",type="selectdir")
  if(!is.na(dirfile)) {
   setwd(dirfile)
   opt = get("setup",envir=mmTK) #get
   opt$workdir = dirsel
   assign("setup",opt,envir=mmTK) #set
  }
 }
 f_openproj = function(h,...) {
  projfile = gWidgets::gfile(text="Open settings",type="open")
  if(!is.na(projfile)) {
   gWidgets::dispose(mainwin)
   dnamatch2::gui(projfile) #send environment into program
  }
 }
 f_saveproj = function(h,...) {
  projfile = gWidgets::gfile(text="Save settings",type="save")
  if(!is.na(projfile)) {
   save(mmTK,file=projfile) #save environment
   print(paste("Settings saved in ",projfile,sep=""))
  }
 }

 #helpfunction for adding folder/ID when clicking button
  f_addFolder = function(h,...) {
    dirfile = gWidgets::gfile(text="Select folder",type="selectdir")
    if(!is.na(dirfile)) {
      fv = addFold(foldadd=dirfile,h$action) #Add folder to environment and file, h=list(action="EVID")
      if(h$action=="EVID") {
        tab2a[1,2][] = fv #update combolist
        gWidgets::enabled(tab2a[2,2]) = TRUE
      }
      if(h$action=="REF") {
        tab2b[1,2][] = fv #update combolist
        gWidgets::enabled(tab2b[2,2]) = TRUE
      }
    }
  }
  f_addID = function(h,...) {
    if(grepl("sel",h$action)) txt = "Insert ID to add for searching (include prefix)"
    if(grepl("pat",h$action)) txt = "Insert prefix pattern for recognizion"
    idADD = gWidgets::ginput(message=txt)
    if( is.na(idADD) || idADD==FALSE || idADD=="") return() #no input given

    fv = addFold(foldadd=idADD,h$action) #Add folder to environment and file, h=list(action="EVID")
    if(h$action=="SIDsel") {
      tab2c[1,2][] = fv #update combolist
      gWidgets::enabled(tab2c[2,2]) = TRUE
    }
    if(h$action=="BIDsel") {
      tab2d[1,2][] = fv #update combolist
      gWidgets::enabled(tab2d[2,2]) = TRUE
    }
    if(h$action=="CIDsel") {
      tab2f[1,2][] = fv #update combolist
      gWidgets::enabled(tab2f[2,2]) = TRUE
    }
    if(h$action=="SIDpat") {
      tab4b[1,2][] = fv #update combolist
      gWidgets::enabled(tab4b[2,2]) = TRUE
    }
    if(h$action=="BIDpat") {
      tab4c[1,2][] = fv #update combolist
      gWidgets::enabled(tab4c[2,2]) = TRUE
    }

  }

  #helpfunction for deleting marked folder when clicking button
  f_delFolder = function(h,...) {
      if(h$action=="EVID") {
        sel = gWidgets::svalue(tab2a[1,2])
        vals = tab2a[1,2][]
        folds = setdiff(vals,sel)
        tab2a[1,2][] = folds  #folds #update combolist
        if(length(folds)==0) {
          tab2a[1,2][] = longspace #folds #update combolist
          gWidgets::svalue(tab2a[1,2]) = longspace 
          gWidgets::enabled(tab2a[2,2]) = FALSE
        } else {
          gWidgets::svalue(tab2a[1,2]) = folds[1]
        }
      }
      if(h$action=="REF") {
        sel = gWidgets::svalue(tab2b[1,2])
        vals = tab2b[1,2][]
        folds = setdiff(vals,sel)
        tab2b[1,2][] = folds #update combolist
        if(length(folds)==0) {
          tab2b[1,2][] = longspace #folds #update combolist
          gWidgets::svalue(tab2b[1,2]) = longspace 
          gWidgets::enabled(tab2b[2,2]) = FALSE
        } else {
          gWidgets::svalue(tab2b[1,2]) = folds[1]
        }
      }
      saveFolds(folds,h$action) #get  h=list(action="EVID")
  } 

 f_delID = function(h,...) {
      if(h$action=="SIDsel") {
        sel = gWidgets::svalue(tab2c[1,2])
        vals = tab2c[1,2][]
        folds = setdiff(vals,sel)
        tab2c[1,2][] = folds  #folds #update combolist
        if(length(folds)==0) {
          tab2c[1,2][] = longspace #folds #update combolist
          gWidgets::svalue(tab2c[1,2]) = longspace 
          gWidgets::enabled(tab2c[2,2]) = FALSE
        } else {
          gWidgets::svalue(tab2c[1,2]) = folds[1]
        }
      }
      if(h$action=="BIDsel") {
        sel = gWidgets::svalue(tab2d[1,2])
        vals = tab2d[1,2][]
        folds = setdiff(vals,sel)
        tab2d[1,2][] = folds #update combolist
        if(length(folds)==0) {
          tab2d[1,2][] = longspace #folds #update combolist
          gWidgets::svalue(tab2d[1,2]) = longspace 
          gWidgets::enabled(tab2d[2,2]) = FALSE
        } else {
          gWidgets::svalue(tab2d[1,2]) = folds[1]
        }
      }
      if(h$action=="CIDsel") {
        sel = gWidgets::svalue(tab2f[1,2])
        vals = tab2f[1,2][]
        folds = setdiff(vals,sel)
        tab2f[1,2][] = folds #update combolist
        if(length(folds)==0) {
          tab2f[1,2][] = longspace #folds #update combolist
          gWidgets::svalue(tab2f[1,2]) = longspace 
          gWidgets::enabled(tab2f[2,2]) = FALSE
        } else {
          gWidgets::svalue(tab2f[1,2]) = folds[1]
        }
      }
      if(h$action=="SIDpat") {
        sel = gWidgets::svalue(tab4b[1,2])
        vals = tab4b[1,2][]
        folds = setdiff(vals,sel)
        tab4b[1,2][] = folds  #folds #update combolist
        if(length(folds)==0) {
          tab4b[1,2][] = longspace #folds #update combolist
          gWidgets::svalue(tab4b[1,2]) = longspace 
          gWidgets::enabled(tab4b[2,2]) = FALSE
        } else {
          gWidgets::svalue(tab4b[1,2]) = folds[1]
        }
      }
      if(h$action=="BIDpat") {
        sel = gWidgets::svalue(tab4c[1,2])
        vals = tab4c[1,2][]
        folds = setdiff(vals,sel)
        tab4c[1,2][] = folds #update combolist
        if(length(folds)==0) {
          tab4c[1,2][] = longspace #folds #update combolist
          gWidgets::svalue(tab4c[1,2]) = longspace 
          gWidgets::enabled(tab4c[2,2]) = FALSE
        } else {
          gWidgets::svalue(tab4c[1,2]) = folds[1]
        }
      }

      saveFolds(folds,h$action) #get  h=list(action="EVID")
  } 



##################################################################################################
########### Program starts #######################################################################
##################################################################################################

 sysdat = Sys.info() #get user info: Will be added to log
 #sysdat[['user']] #USER INFO

 #change working directory to the one stored in mmTK-environment
 wd=get("setup",envir=mmTK)$workdir #assign working directory to mmTK-environment
 if(!is.null(wd)) {
   tryCatch( { setwd(wd) }, error=function(e) print("Folder not found. Using existing") )
 }
 
 #Main window:
 mainwin <- gWidgets::gwindow(softname, visible=FALSE, width=mwW,height=mwH)
 nb = gWidgets::gnotebook(container=mainwin)
 tabanalyse = gWidgets::ggroup(expand=TRUE,spacing=spc,container=nb,label="Analyse") #tab1: (select project and file storage)
 tabdata = gWidgets::ggroup(expand=TRUE,spacing=spc,container=nb,label="Data setup") #tab2: (select what data to search)
 tabsearch = gWidgets::ggroup(expand=TRUE,spacing=spc,container=nb,label="Search setup") #tab3: (select prefilter thresholds and model settings)
 tabpattern = gWidgets::ggroup(expand=TRUE,spacing=spc,container=nb,label="Pattern setup") #tab4: (User specified pattern setup)

 gWidgets::svalue(nb) <- 1 #initial start in first tab


#####################################################
###############Tab 1: Analysis:######################
#####################################################

  tab1 <- gWidgets::glayout(spacing=spc,container=tabanalyse ) 

  tab1a = gWidgets::glayout(spacing=spc,container=(tab1[1,1] <-gWidgets::gframe("Save/Load settings",container=tab1))) 
  tab1a[1,1] <- gWidgets::gbutton("Save settings to file",container=tab1a,handler=f_saveproj)
  tab1a[2,1] <- gWidgets::gbutton("Load settings from file",container=tab1a,handler=f_openproj)
  tab1a[3,1] <- gWidgets::gbutton("RESTART",container=tab1a,handler=
	function(h,...) {
  	 gWidgets::dispose(mainwin)
	 dnamatch2::gui(envirfile) #restart GUI with current environment (useful if projects are considered)
     }
  )

  tab1b = gWidgets::glayout(spacing=spc,container=(tab1[2,1] <-gWidgets::gframe("Directories",container=tab1))) 
  tab1b[1,1] <- gWidgets::gbutton("Selected working directory:",container=tab1b,handler=
	function(h,...) { 
      fsel = gWidgets::gfile(text="Select directory",type="selectdir")
      if(!is.na(fsel)) {
	  setwd(fsel) #actually set work directory
       opt = get("setup",envir=mmTK) #get option vals
       opt$workdir = fsel
       assign("setup",opt,envir=mmTK) #set to envir again
       saveSetup(opt) #Save to file
       gWidgets::svalue(tab1b[1,2]) = fsel
      }
  })
  tab1b[1,2] <- gWidgets::glabel(get("setup",envir=mmTK)$workdir,container=tab1b)

  tab1b[2,1] <- gWidgets::gbutton("Selected name for matchfile:",container=tab1b,handler=
	function(h,...) { 
      opt = get("setup",envir=mmTK) #get option vals
      fn = gWidgets::ginput(message="Select name",text=opt$matchfile)
      if( is.na(fn) || fn==FALSE || fn=="") return() #no input given
      if(!grepl("\\.",fn)) fn = paste0(fn,".csv") #automatically add extension
      opt$matchfile = fn
      assign("setup",opt,envir=mmTK) #set to envir again
      saveSetup(opt) #Save to file
      gWidgets::svalue(tab1b[2,2]) = fn
  })
  tab1b[2,2] <- gWidgets::glabel(get("setup",envir=mmTK)$matchfile,container=tab1b)

  tab1b[3,1] <- gWidgets::gbutton("Selected name for session folder:",container=tab1b,handler=
	function(h,...) { 
      fn = gWidgets::ginput(message="Select name",text=optL$sessionfold)
      if( is.na(fn) || fn==FALSE || fn=="") return() #no input given
      dir.create(fn, showWarnings = FALSE) #create folder
      optL = get("setup",envir=mmTK) #get option vals
      optL$sessionfold = fn
      assign("setup",optL,envir=mmTK) #set to envir again
      saveSetup(optL) #Save to file
      gWidgets::svalue(tab1b[3,2]) = fn
  })
  tab1b[3,2] <- gWidgets::glabel(get("setup",envir=mmTK)$sessionfold,container=tab1b)

  tab1c = gWidgets::glayout(spacing=spc,container=(tab1[3,1] <-gWidgets::gframe("Analyse",container=tab1))) 
  tab1b[1,1] <- gWidgets::glabel("Remember to check the SETUP before searching:",container=tab1c)
  tab1b[2,1] <- gWidgets::gbutton("Perform search",container=tab1c, handler= runAnalysis )


#########################################################
###############Tab 2: Setup for format ##################
#########################################################

  tab2 <- gWidgets::glayout(spacing=spc,container=tabdata ) 

  tab2a = gWidgets::glayout(spacing=spc,container=(tab2[2,1] <-gWidgets::gframe("Selecting evidence folders",container=tab2))) 
  tab2a[1,1] <- gWidgets::glabel("Selected folders:",container=tab2a)
  tab2a[2,1] <- gWidgets::gbutton("Add a folder",container=tab2a,handler=f_addFolder,action="EVID")
  tab2a[2,2] <- gWidgets::gbutton("Remove marked folder",container=tab2a,handler=f_delFolder,action="EVID")
  folds = getFolds("EVID") 
  if(is.null(folds) || length(folds)==0)  folds = longspace #numeric()
  tab2a[1,2] <- gWidgets::gcombobox(items=folds,container=tab2a)
  if(folds[1] == longspace) gWidgets::enabled(tab2a[2,2]) = FALSE
 
  tab2b = gWidgets::glayout(spacing=spc,container=(tab2[3,1] <-gWidgets::gframe("Selecting reference folders",container=tab2))) 
  tab2b[1,1] <- gWidgets::glabel("Selected folders:",container=tab2b)
  tab2b[2,1] <- gWidgets::gbutton("Add a folder",container=tab2b,handler=f_addFolder,action="REF")
  tab2b[2,2] <- gWidgets::gbutton("Remove marked folder",container=tab2b,handler=f_delFolder,action="REF")
  folds = getFolds("REF") 
  if(is.null(folds) || length(folds)==0)  folds = longspace #numeric()
  tab2b[1,2] <- gWidgets::gcombobox(items=folds,container=tab2b)
  if(folds[1] == longspace) gWidgets::enabled(tab2b[2,2]) = FALSE


  tab2c = gWidgets::glayout(spacing=spc,container=(tab2[4,1] <-gWidgets::gframe("Selecting specific SampleIDs (SIDs)",container=tab2))) 
  tab2c[1,1] <- gWidgets::glabel("Selected SIDs:",container=tab2c)
  tab2c[2,1] <- gWidgets::gbutton("Add an ID",container=tab2c,handler=f_addID,action="SIDsel")
  tab2c[2,2] <- gWidgets::gbutton("Remove marked ID",container=tab2c,handler=f_delID,action="SIDsel")
  folds = getFolds("SIDsel") 
  if(is.null(folds) || length(folds)==0)  folds = shortspace #numeric()
  tab2c[1,2] <- gWidgets::gcombobox(items=folds,container=tab2c)
  if(folds[1] == shortspace) gWidgets::enabled(tab2c[2,2]) = FALSE

  tab2d = gWidgets::glayout(spacing=spc,container=(tab2[5,1] <-gWidgets::gframe("Selecting specific BatchIDs (BIDs)",container=tab2))) 
  tab2d[1,1] <- gWidgets::glabel("Selected BIDs:",container=tab2d)
  tab2d[2,1] <- gWidgets::gbutton("Add an ID",container=tab2d,handler=f_addID,action="BIDsel")
  tab2d[2,2] <- gWidgets::gbutton("Remove marked ID",container=tab2d,handler=f_delID,action="BIDsel")
  folds = getFolds("BIDsel") 
  if(is.null(folds) || length(folds)==0)  folds = shortspace #numeric()
  tab2d[1,2] <- gWidgets::gcombobox(items=folds,container=tab2d)
  if(folds[1] == shortspace) gWidgets::enabled(tab2d[2,2]) = FALSE

  tab2f = gWidgets::glayout(spacing=spc,container=(tab2[6,1] <-gWidgets::gframe("Selecting specific CaseIDs (CIDs)",container=tab2))) 
  tab2f[1,1] <- gWidgets::glabel("Selected CIDs:",container=tab2f)
  tab2f[2,1] <- gWidgets::gbutton("Add an ID",container=tab2f,handler=f_addID,action="CIDsel")
  tab2f[2,2] <- gWidgets::gbutton("Remove marked ID",container=tab2f,handler=f_delID,action="CIDsel")
  folds = getFolds("CIDsel") 
  if(is.null(folds) || length(folds)==0)  folds = shortspace #numeric()
  tab2f[1,2] <- gWidgets::gcombobox(items=folds,container=tab2f)
  if(folds[1] == shortspace) gWidgets::enabled(tab2f[2,2]) = FALSE


  tab2e = gWidgets::glayout(spacing=spc,container=(tab2[1,1] <-gWidgets::gframe("Population frequencies",container=tab2))) 
  tab2e[1,1] <- gWidgets::gbutton("Selected frequency file:",container=tab2e,handler = 
	function(h,...) { 
      fsel = gWidgets::gfile(text="Select frequency file",type="open")
      if(!is.na(fsel)) {
       opt = get("setup",envir=mmTK) #get
       opt$freqfile = fsel
       assign("setup",opt,envir=mmTK) #set to envir
       saveSetup(opt) #Save to file
       gWidgets::svalue(tab2e[1,2]) = fsel
      }
  })
  tab2e[1,2] <- gWidgets::glabel(get("setup",envir=mmTK)$freqfile,container=tab2e)


####################################################
###############Tab 3: Search setup #################
####################################################

  txtbool = c("NO","YES")
  tab3 <- gWidgets::glayout(spacing=spc,container=tabsearch ) 

  tab3a = gWidgets::glayout(spacing=spc,container=(tab3[1,2] <-gWidgets::gframe("Search options",container=tab3))) 
  tab3a[1,1] <- gWidgets::glabel("Search within same cases (CID):",container=tab3a)
  tab3a[1,2] <- gWidgets::gradio(items=txtbool,container=tab3a,horizontal = TRUE,selected=sum(get("setup",envir=mmTK)$sameCID)+1 )
  tab3a[2,1] <- gWidgets::glabel("Search between stains:",container=tab3a)
  tab3a[2,2] <- gWidgets::gradio(items=txtbool,container=tab3a,horizontal = TRUE,selected=sum(get("setup",envir=mmTK)$betweensamples)+1)

  tab3b = gWidgets::glayout(spacing=spc,container=(tab3[2,2] <-gWidgets::gframe("Time windows",container=tab3))) 
  tab3b[1,1] <- gWidgets::glabel("Number of days back (days):",container=tab3b)
  tab3b[1,2] <- gWidgets::gedit(get("setup",envir=mmTK)$Thist,container=tab3b)

  tab3b[2,1] <- gWidgets::glabel("Time difference between matches (days):",container=tab3b)
  tab3b[2,2] <- gWidgets::gedit(get("setup",envir=mmTK)$timediff,container=tab3b)

  tab3b[3,1] <- gWidgets::glabel("Search time (YY-MM-DD-HH-MM-SS)",container=tab3b)
  tab3b[3,2] <- gWidgets::gedit(get("setup",envir=mmTK)$searchtime,container=tab3b)
  tab3b[4,2] <- gWidgets::gbutton("Update time stamp (current time)",container=tab3b,handler=function(h,...) {
    gWidgets::svalue(tab3b[3,2]) = format(Sys.time(),format=timestamp) #set time stamp to now
  })
  
  tab3c = gWidgets::glayout(spacing=spc,container=(tab3[1,1] <-gWidgets::gframe("Score thresholds",container=tab3))) 
  tab3c[1,1] <- gWidgets::glabel("Matching allele counting (MAC):",container=tab3c)
  tab3c[1,2] <- gWidgets::gedit(get("setup",envir=mmTK)$threshMAC,container=tab3c)

  tab3c[2,1] <- gWidgets::glabel("Qualitative LR:",container=tab3c)
  tab3c[2,2] <- gWidgets::gedit(get("setup",envir=mmTK)$threshLRqual,container=tab3c)

  tab3c[3,1] <- gWidgets::glabel("Quantitative LR:",container=tab3c)
  tab3c[3,2] <- gWidgets::gedit(get("setup",envir=mmTK)$threshLRquan,container=tab3c)


  tab3d = gWidgets::glayout(spacing=spc,container=(tab3[2,1] <-gWidgets::gframe("Model setup",container=tab3))) 
  tab3d[1,1] <- gWidgets::glabel("Set kit:",container=tab3d)
  kits <- c( emptyName,euroformix::getKit()) #include empty kit as a possibility (this is default)
  tab3d[1,2] <- gWidgets::gcombobox(kits,container=tab3d,selected=which(get("setup",envir=mmTK)$kit==kits))
  tab3d[2,1] <- gWidgets::glabel("Drop-in prob=",container=tab3d)
  tab3d[2,2] <- gWidgets::gedit(get("setup",envir=mmTK)$pC,container=tab3d)
  tab3d[3,1] <- gWidgets::glabel("Lambda param=",container=tab3d)
  tab3d[3,2] <- gWidgets::gedit(get("setup",envir=mmTK)$lambda,container=tab3d)
  tab3d[4,1] <- gWidgets::glabel("Min Freq=",container=tab3d)
  tab3d[4,2] <- gWidgets::gedit(get("setup",envir=mmTK)$minFreq,container=tab3d)

  tab3e = gWidgets::glayout(spacing=spc,container=(tab3[3,1] <-gWidgets::gframe("Prefilter thresholds",container=tab3))) 
  tab3e[1,1] <- gWidgets::glabel("Analytical threshold (AT)",container=tab3e)
  tab3e[1,2] <- gWidgets::gedit(get("setup",envir=mmTK)$threshHeight,container=tab3e)
  tab3e[2,1] <- gWidgets::glabel("Stutter rate threshold",container=tab3e)
  tab3e[2,2] <- gWidgets::gedit(get("setup",envir=mmTK)$threshStutt,container=tab3e)
  tab3e[3,1] <- gWidgets::glabel("Major extraction rate threshold",container=tab3e)
  tab3e[3,2] <- gWidgets::gedit(get("setup",envir=mmTK)$threshMaj,container=tab3e)
  tab3e[4,1] <- gWidgets::glabel("Minimum loci requirement (Evid)",container=tab3e)
  tab3e[4,2] <- gWidgets::gedit(get("setup",envir=mmTK)$minLocStain,container=tab3e)
  tab3e[5,1] <- gWidgets::glabel("Minimum loci requirement (Maj)",container=tab3e)
  tab3e[5,2] <- gWidgets::gedit(get("setup",envir=mmTK)$minLocMaj,container=tab3e)

  tab3f = gWidgets::glayout(spacing=spc,container=(tab3[3,2] <-gWidgets::gframe("Other options",container=tab3))) 
  tab3f[1,1] <- gWidgets::glabel("Plot score histogram in search",container=tab3f)
  tab3f[1,2] <- gWidgets::gradio(txtbool,horizontal=TRUE,container=tab3f,selected=sum(get("setup",envir=mmTK)$printHistPlots)+1 )
  tab3f[2,1] <- gWidgets::glabel("Write detailed score info to file",container=tab3f)
  tab3f[2,2] <- gWidgets::gradio(txtbool,horizontal=TRUE,container=tab3f,selected=sum(get("setup",envir=mmTK)$writeScores)+1 )
  tab3f[3,1] <- gWidgets::glabel("Print graph of matches:",container=tab3f)
  tab3f[3,2] <- gWidgets::gradio(items=txtbool,container=tab3f,horizontal = TRUE,selected=sum(get("setup",envir=mmTK)$printGraph)+1)

  tab3f[4,1] <- gWidgets::glabel("Maximum number of contributors:",container=tab3f)
  tab3f[5,1] <- gWidgets::glabel("Qualitative LR:",container=tab3f)
  tab3f[5,2] <- gWidgets::gedit(get("setup",envir=mmTK)$maxKqual,container=tab3f)
  tab3f[6,1] <- gWidgets::glabel("Quantitative LR:",container=tab3f)
  tab3f[6,2] <- gWidgets::gedit(get("setup",envir=mmTK)$maxKquan,container=tab3f)

  tab3[4,1] = gWidgets::gbutton("Save settings",container=tab3,handler = 
	function(h,...) { 
	#CHECKS: checkProb[0,1];checkPositive[(0;checkPosInteger[(1;  (val,what,strict)
	#names(optL)
     opt = get("setup",envir=mmTK) #get uptodate settings
 
	#Values should be zero or positive:
	opt$Thist = checkPositive(val(tab3b[1,2]),"Search history setting(days)",strict=FALSE)
	opt$timediff = checkPositive(val(tab3b[2,2]),"Time-difference setting (days)",strict=FALSE)
	opt$threshLRqual = checkPositive(val(tab3c[2,2]),"Qualitative likelihood threshold",strict=FALSE)
	opt$threshLRquan = checkPositive(val(tab3c[3,2]),"Quantitative likelihood threshold",strict=FALSE)

	#Values should be positive:
	opt$lambda = checkPositive(val(tab3d[3,2]),"Lambda parameter",strict=TRUE)
	opt$threshHeight = checkPositive(val(tab3e[1,2]),"Detection threshold (AT)",strict=TRUE)

	#Values should be between zero and one:
	opt$threshMAC = checkProb(val(tab3c[1,2]),"MAC threshold")
	opt$pC =  checkProb(val(tab3d[2,2]),"Drop-in probability")
	opt$minFreq = checkProb(val(tab3d[4,2]),"Minimum frequency")
	opt$threshStutt = checkProb(val(tab3e[2,2]),"Stutter rate threshold")
	opt$threshMaj = checkProb(val(tab3e[3,2]),"Major extraction rate threshold")

	#values should be integers (positive or zero) 
	opt$minLocStain = checkPosInteger(val(tab3e[4,2]),"Minimum loci (EVID) threshold",strict=FALSE)
	opt$minLocMaj = checkPosInteger(val(tab3e[5,2]),"Minimum loci (MAJ) threshold",strict=FALSE)
	opt$maxKqual = checkPosInteger(val(tab3f[5,2]),"Maximum number of contributors (QUAL)",strict=FALSE)
	opt$maxKquan = checkPosInteger(val(tab3f[6,2]),"Maximum number of contributors (QUAN)",strict=FALSE)

	#NO/YES choice
	opt$sameCID = gWidgets::svalue(tab3a[1,2])==txtbool[2]
	opt$betweensamples = gWidgets::svalue(tab3a[2,2])==txtbool[2]
	opt$printHistPlots = gWidgets::svalue(tab3f[1,2])==txtbool[2]
	opt$writeScores = gWidgets::svalue(tab3f[2,2])==txtbool[2]
	opt$printGraph = gWidgets::svalue(tab3f[3,2])==txtbool[2]

	#Text format:
	opt$kit = gWidgets::svalue(tab3d[1,2])
     tmp = gWidgets::svalue(tab3b[3,2]) #get search time
     tmp2 = as.POSIXct(tmp,format=timestamp) #convert to time format
     if(is.na(tmp) || is.na(tmp2)) errorMessage("The search time was not correctly specified.\nPlease check format!") 
     opt$searchtime = tmp #keep timestamp only if successful

 	#Store settings:
     assign("setup",opt,envir=mmTK) #set to envir
   	saveSetup(opt)
	gWidgets::gmessage("Settings were saved successfully!",title="Message") 
  }) #Done saving settings

#########################################################
###############Tab 4: Setup for format ##################
#########################################################
  tab4 <- gWidgets::glayout(spacing=spc,container=tabpattern) 

  tab4a = gWidgets::glayout(spacing=spc,container=(tab4[1,1] <-gWidgets::gframe("Patterns (prefix) of IDs",container=tab4))) 
  tab4a[1,1] <- gWidgets::gbutton("Set pattern for Separating IDs:",container=tab4a,handler=
    function(h,...) {
     pat = gWidgets::ginput(message="Set separator pattern for ID",text=get("setup",envir=mmTK)$IDsep)
     if( is.na(pat) || pat==FALSE || pat=="") return() #no input given
     opt = get("setup",envir=mmTK) #get
     opt$IDsep= pat 
     assign("setup",opt,envir=mmTK) #set to envir
     saveSetup(opt) #Save to file
     gWidgets::svalue(tab4a[1,2]) = pat
    })
  tab4a[1,2] <- gWidgets::glabel(get("setup",envir=mmTK)$IDsep,container=tab4a)

  tab4b = gWidgets::glayout(spacing=spc,container=(tab4[2,1] <-gWidgets::gframe("Set pattern for SampleIDs (SIDs):",container=tab4))) 
  tab4b[1,1] <- gWidgets::glabel("Required pattern(s):",container=tab4b)
  tab4b[2,1] <- gWidgets::gbutton("Add a pattern",container=tab4b,handler=f_addID,action="SIDpat")
  tab4b[2,2] <- gWidgets::gbutton("Remove a pattern",container=tab4b,handler=f_delID,action="SIDpat")
  folds = getFolds("SIDpat") 
  if(is.null(folds) || length(folds)==0)  folds = shortspace #numeric()
  tab4b[1,2] <- gWidgets::gcombobox(items=folds,container=tab4b)
  if(folds[1] == shortspace) gWidgets::enabled(tab4b[2,2]) = FALSE

  tab4c = gWidgets::glayout(spacing=spc,container=(tab4[3,1] <-gWidgets::gframe("Set pattern for Batch files (BIDs):",container=tab4))) 
  tab4c[1,1] <- gWidgets::glabel("Required pattern(s):",container=tab4c)
  tab4c[2,1] <- gWidgets::gbutton("Add a pattern",container=tab4c,handler=f_addID,action="BIDpat")
  tab4c[2,2] <- gWidgets::gbutton("Remove a pattern",container=tab4c,handler=f_delID,action="BIDpat")
  folds = getFolds("BIDpat") 
  if(is.null(folds) || length(folds)==0)  folds = shortspace #numeric()
  tab4c[1,2] <- gWidgets::gcombobox(items=folds,container=tab4c)
  if(folds[1] == shortspace) gWidgets::enabled(tab4c[2,2]) = FALSE

  
  gWidgets::visible(mainwin) = TRUE

} #End function

