#!/usr/bin/Rscript


cmd.help <- function() {
    cat("\n")
    cat("## Mandatory parameters:\n")
    cat("  -F   The factor name given by the -f option in macs2-idr.sh \n")
    cat("  -D   The macs2-idr output directory as given by the -o option in macs2-idr.sh \n")
    cat("\n")
    cat("\n")
}


# Parse command line arguments and extract them into an associate array.
# Check if the required arguments are all satisfied.

parseArgs <- function(args, manditories) {
    if(length(args) %% 2 == 1 || length(args) == 0) {
        cat('Unpaired argument and value.\n')
        return(NULL)
    }
    n.i <- seq(1, length(args), by=2)
    v.i <- seq(2, length(args), by=2)
    args.name <- args[n.i]
    args.value <- args[v.i]

    # Check if required argument values are supplied.
    miss_tag <- F
    man.bool <- manditories %in% args.name
    if(!all(man.bool)){
        cat(paste('Missing argument: ', paste(manditories[!man.bool],
                                              collapse=','), '.', sep='')
           )
        miss_tag <- T
    }
    if(miss_tag){
        res <- NULL
    }else{
        res <- args.value
        names(res) <- args.name
    }
    res
}


# Read command arguments.
args <- commandArgs(T)
if(length(args) < 1) {
    cmd.help()
    stop("No arguments provided.\n")
}


# Parse arguments.
args.tbl <- parseArgs(args, c('-I', '-O'))
if(is.null(args.tbl)){
    cmd.help()
    stop('Error in parsing command line arguments. Stop.\n')
}
factor <- args.tbl['-F']  # iname: input zip file name
baseDir <- args.tbl['-D']  # oname: output file root name

peakDir = paste(baseDir, factor, sep="/")



# Dependencies
if(!suppressMessages(require(GenomicRanges, warn.conflicts=F))) {
    source("http://bioconductor.org/biocLite.R",echo=FALSE,verbose=FALSE)
    biocLite("GenomicRanges",suppressUpdates=TRUE,verbose=FALSE)
    if(!suppressMessages(require(GenomicRanges, warn.conflicts=F))) {
        stop('Loading package GenomicRanges failed!')
    }
}

#import IDR bed files to GRanges
#The extra columns in the bed files are not supported by rtracklayer
importPeaks = function(file) {
  raw = read.table(file)
  peaks = GRanges(seqnames = raw$V1, ranges = IRanges(start=raw$V2, end=raw$V3), score=raw$V5)
  return(peaks)
}

#converts the colour(s) to the MCRI version of the colour, if available, otherwise returns the input.
#al sets the opaqueness.
mcri = function(col=0, al=1) {
  if ( col[1] == 0 ) {
    cat('Use: mcri(\'colour\'), returning an official MCRI colour.\nAvailable MCRI colours are:\n\ndarkblue\nblue\nlightblue\nazure\ngreen\norange\nviolet\ncyan\nred\nmagenta (aka rose).\n\nReturning default blue.\n')
    return(mcri('blue'))
  }
  if ( length(col) > 1 ) return(sapply(col, function(c) mcri(c, al)))
  if ( is.numeric(col) ) {
    if ( col == 1 ) col = 'blue'
    else if ( col == 2 ) col = 'orange'
    else if ( col == 3 ) col = 'green'
    else if ( col == 4 ) col = 'magenta'
    else if ( col == 5 ) col = 'cyan'
    else if ( col == 6 ) col = 'red'
    else if ( col == 7 ) col = 'violet'
    else if ( col == 8 ) col = 'darkblue'
    else if ( col == 9 ) col = 'azure'
    else if ( col == 10 ) col = 'lightblue'
    else col = 'black'
  }
  ret = 0
  if ( col == 'darkblue') ret = rgb(9/255, 47/255, 94/255, al)
  if ( col == 'blue') ret = rgb(0, 83/255, 161/255, al)
  if ( col == 'lightblue') ret = rgb(0, 165/255, 210/255, al)
  if ( col == 'azure') ret = rgb(0, 173/255, 239/255, al)
  if ( col == 'green') ret = rgb(141/255, 198/255, 63/255, al)
  if ( col == 'orange') ret = rgb(244/255, 121/255, 32/255, al)
  if ( col == 'violet') ret = rgb(122/255, 82/255, 199/255, al)
  if ( col == 'cyan') ret = rgb(0/255, 183/255, 198/255, al)
  if ( col == 'red') ret = rgb(192/255, 80/255, 77/255, al)
  if ( col == 'magenta' | col == 'rose') ret = rgb(236/255, 0/255, 140/255, al)
  if ( ret == 0 ) ret = do.call(rgb, as.list(c(col2rgb(col)/255, al)))
  return(ret)
}


compare_IDR_peaks = function(factor, peakDir) {

	peaks1pc = importPeaks( paste(peakDir,  "/", factor, ".MACS2_IDR-0.01_peaks.bed", sep="") )
	peaks2pc = importPeaks( paste(peakDir,  "/", factor, ".MACS2_IDR-0.02_peaks.bed", sep="") )
	peaks3pc = importPeaks( paste(peakDir,  "/", factor, ".MACS2_IDR-0.03_peaks.bed", sep="") )
	peaks4pc = importPeaks( paste(peakDir,  "/", factor, ".MACS2_IDR-0.04_peaks.bed", sep="") )
	peaks5pc = importPeaks( paste(peakDir,  "/", factor, ".MACS2_IDR-0.05_peaks.bed", sep="") )
	peaks6pc = importPeaks( paste(peakDir,  "/", factor, ".MACS2_IDR-0.06_peaks.bed", sep="") )
	peaks7pc = importPeaks( paste(peakDir,  "/", factor, ".MACS2_IDR-0.07_peaks.bed", sep="") )
	peaks8pc = importPeaks( paste(peakDir,  "/", factor, ".MACS2_IDR-0.08_peaks.bed", sep="") )
	peaks9pc = importPeaks( paste(peakDir,  "/", factor, ".MACS2_IDR-0.09_peaks.bed", sep="") )
	peaks10pc = importPeaks( paste(peakDir,  "/", factor, ".MACS2_IDR-0.1_peaks.bed", sep="") )

	pdf( file=paste0(peakDir, "/", factor, "_IDR_peakCalling_comparison.pdf", sep="") )

	col=c( mcri('red',1), mcri('red',0.9), mcri('red',0.8), mcri('red',0.7), mcri('red',0.6), mcri('red',0.5), mcri('red',0.4), mcri('red',0.3), mcri('red',0.2), mcri('red',0.1) )
	names=c( '1', '2', '3', '4', '5', '6', '7', '8', '9', '10' )
	xlab='Irreproducible Discovery Rate (%)'

	boxplot ( log10(peaks1pc$score), log10(peaks2pc$score), log10(peaks3pc$score), log10(peaks4pc$score), log10(peaks5pc$score),
	          log10(peaks6pc$score), log10(peaks7pc$score), log10(peaks8pc$score), log10(peaks9pc$score), log10(peaks10pc$score),
	          col=col,
	          names=names,
	          xlab=xlab, ylab=expression('log'[10]*'(peak score)'), main=paste0('Peak scores of ', factor, ' peaks'), pch=16)

	bp = barplot ( c( length(peaks1pc), length(peaks2pc), length(peaks3pc), length(peaks4pc), length(peaks5pc),
	             length(peaks6pc), length(peaks7pc), length(peaks8pc), length(peaks9pc), length(peaks10pc) ),
	          col=col,
	          names=names,
	          xlab=xlab, ylab='Number of peaks', main=paste0('Number of peaks called for ', factor), pch=16, )
	text(bp, 0, c( length(peaks1pc), length(peaks2pc), length(peaks3pc), length(peaks4pc), length(peaks5pc),
	             length(peaks6pc), length(peaks7pc), length(peaks8pc), length(peaks9pc), length(peaks10pc) ),cex=1, pos=3, offset=1.5,  srt=90)

	dev.off()
}

compare_IDR_peaks(factor, peakDir)
