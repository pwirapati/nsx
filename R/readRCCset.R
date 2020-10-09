readRCCset <-
function(
  path,                    # path(s) to be recursively searched for RCC files
  ext = "\\.RCC(\\.gz)?$", # regexp for recognizing RCC files
  trim = c(ext, paste0(path,"/")),      # regexps for trimming filenames -> row.names
  drop = TRUE)             # if there is only one RCCset, drop the list
  {
  
  #: find RCC files recursively, read and parse into intermediate
  #: format (list of vectors for the blocks)
  ssdata <-
  sapply(
    list.files(path, ext, full.names=TRUE, recursive=TRUE),
    function(fn) read1File(fn),
    simplify=FALSE)
  if( length(ssdata) == 0 ) stop("No RCC is processed")

  #: trim the filenames to create shorter row.names
  for(pat in trim)
    if (is.character(pat))
      names(ssdata) <- gsub( pat, "", names(ssdata))
  
  GeneRLF <- sapply(ssdata, function(u) u$header["GeneRLF"])

  RCCsetlist <-
  sapply( 
    sort(unique(GeneRLF)),
    function(u) collate1GeneRLF(ssdata[GeneRLF == u]),
    simplify=FALSE)
  
  return(
    if (drop == TRUE && length(RCCsetlist) == 1 )
      RCCsetlist[[1]]
    else
      RCCsetlist)
  }

#| read and parse single RCC files. Produce intermediate data to
#| obtain some header information first.
#|
read1File <-
function(filename)
  {
  if(length(filename)==0) return(NULL)

  ### read all lines
  fi <- gzfile(filename)
  open(fi)
  li <- readLines(fi)
  close(fi)
  #? TODO: check for IO error as well as wrong format

  #: c(i[1],i[2])  divide lines into header, data, the rest
  i <- grep("^<.?Code_Summary>$", li) 
  
  #: extract headers
  v <- strsplit(
         grep("^[^<]", li[1:i[1]], value=TRUE),
         ",",
         fixed=TRUE )
  header <- sapply(v, function(u) u[2])
  names(header) <- sapply(v, function(u) u[1])
  names(header)[grep("^ID$", names(header))] <- c("SampleID","LaneID")
    # This resolves uniqueness of key strings

  #: extract count data for a given reporter class
  getDataClass <-
    function( pattern )
    {
    v <- simplify2array(
           strsplit(
             grep(pattern, li[i[1]:i[2]], value=TRUE),
             ",",
             fixed=TRUE ))
    data <- as.numeric(v[4, ])
    names(data) <- gsub(" \\(\\+\\+\\+.*$", "", v[2, ])
      # also filter '(+++ See Message below ...)'

    return(data)
    }

  list(
    filename=filename,
    header=header,
    endog=getDataClass("^Endogenous"),
    hk=getDataClass("^Housekeeping"),
    pos=getDataClass("^Positive"),
    neg=getDataClass("^Negative"))
  }


#| collate single-sample data from the same codeset (as identified by
#| GeneRLF string) into a dataframe or matrix.
#|
collate1GeneRLF <-
function(ssdata)
  {
  collateField <- function(f) t(sapply(ssdata,function(u) u[[f]]))

  #: turn into matrices with samples as rows
  x <- data.frame(
    filename   = I(c(collateField("filename"))),
    header = I(as.data.frame(
                    collateField("header"),
                    stringsAsFactors=FALSE)),
    endog     = I(collateField("endog")),
    pos       = I(collateField("pos")),
    neg       = I(collateField("neg")),
    hk        = I(collateField("hk"))
  )

  for(k in c("LaneID", "StagePosition", "FovCount", "FovCounted"))
    x$header[[k]] <- as.integer(x$header[[k]])
  x$header$BindingDensity <- as.numeric(x$header$BindingDensity)
  x$header$Date <- as.Date(x$header$Date,"%Y%m%d")
  
  class(x) <- c("RCCset","data.frame")
  return(x)
  }

