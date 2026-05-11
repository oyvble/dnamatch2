#This function is written since the encoding in  gWidgets2::gfile is fixed to UTF-8 which doesn't handle special letters
mygfile <- function(text,type,filter=list(),initf=NULL) { #Possible bug: text ignored when type="selectdir"
  file <- gWidgets2::gfile(text=text,type=type,filter=filter,initial.filename=initf)
  Encoding(file ) <- options()$encoding #Set default to local encoder: Makes it possible to use Language-text files in own language
  return(file)
}

#helpfunction to convert data from file to list
scanFileToList = function(file) { 
  as.list( scan(file=file,what=character(),quiet=TRUE,sep="\n") )  #convert to list 
}

#helpfunction which checks that at value is in interval of [0,1]
checkProb = function(x,what) {
  if(x < 0 || x>1) {
    gWidgets2::gmessage(paste0(what," must be specified in interval [0,1] "),title="Wrong input",icon="error")
    stop("Wrong user-input")
  }
  return(x)
}
checkPositive = function(x,what,strict=FALSE) {
  if(x < 0 ) {
    gWidgets2::gmessage(paste0(what," cannot be a negative number"),title="Wrong input",icon="error")
    stop("Wrong user-input")
  }
  if(strict && x==0) {
    gWidgets2::gmessage(paste0(what," cannot be zero"),title="Wrong input",icon="error")
    stop("Wrong user-input")
  }
  return(x)
}
checkPosInteger = function(x,what,strict=TRUE) {
  if(x < 1 || round(x)!=x ) {
    gWidgets2::gmessage(paste0(what," must be a positive integer"),title="Wrong input",icon="error")
    stop("Wrong user-input")
  }
  if(strict && x==0) {
    gWidgets2::gmessage(paste0(what," cannot be zero"),title="Wrong input",icon="error")
    stop("Wrong user-input")
  }
  return(x)
}
val = function(wid) { #helpfunction to get value from widget
  tmp = as.numeric(gWidgets2::svalue(wid))
  if(is.na(tmp )) {
    gWidgets2::gmessage(paste0("Text found where number was expected:\n",gWidgets2::svalue(wid)))
    return()
  }
  return(tmp )
}

errorMessage = function(msg) {  #Helpfunction to throw error message to user + stop running 
  gWidgets2::gmessage(msg,title="Error",icon="error")
  stop(msg)
}

NullIfEmpty = function(x) {
  if(length(x)==0) {
    return(NULL)
  } else {
    return(unlist(x)) #input argument in dnamatch2 must be vectors instead of lists
  }
}
