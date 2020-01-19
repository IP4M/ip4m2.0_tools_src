

#@author Y. Guitton
getTIC <- function(file,rtcor=NULL) {
    object <- xcmsRaw(file)
    cbind(if (is.null(rtcor)) object@scantime else rtcor, rawEIC(object,mzrange=range(object@env$mz))$intensity)
}

##
##  overlay TIC from all files in current folder or from xcmsSet, create pdf
##
#@author Y. Guitton
getTICs <- function(xcmsSet=NULL,files=NULL, pdfname="TICs.pdf",rt=c("raw","corrected")) {
    cat("Creating TIC pdf...\n")

    if (is.null(xcmsSet)) {
        filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
        filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
        if (is.null(files))
            files <- getwd()
        info <- file.info(files)
        listed <- list.files(files[info$isdir], pattern = filepattern, recursive = TRUE, full.names = TRUE)
        files <- c(files[!info$isdir], listed)
    } else {
        files <- filepaths(xcmsSet)
    }

    phenoDataClass<-as.vector(levels(xcmsSet@phenoData[,1])) #sometime phenoData have more than 1 column use first as class
    classnames<-vector("list",length(phenoDataClass))
    for (i in 1:length(phenoDataClass)){
        classnames[[i]]<-which( xcmsSet@phenoData[,1]==phenoDataClass[i])
    }

    N <- length(files)
    TIC <- vector("list",N)

    for (i in 1:N) {
        if (!is.null(xcmsSet) && rt == "corrected")
            rtcor <- xcmsSet@rt$corrected[[i]] else
        rtcor <- NULL
        TIC[[i]] <- getTIC(files[i],rtcor=rtcor)
    }

    pdf(pdfname,w=16,h=10)
    cols <- rainbow(N)
    lty = 1:N
    pch = 1:N
    #search for max x and max y in TICs
    xlim = range(sapply(TIC, function(x) range(x[,1])))
    ylim = range(sapply(TIC, function(x) range(x[,2])))
    ylim = c(-ylim[2], ylim[2])


    ##plot start
    if (length(phenoDataClass)>2){
        for (k in 1:(length(phenoDataClass)-1)){
            for (l in (k+1):length(phenoDataClass)){
                #print(paste(phenoDataClass[k],"vs",phenoDataClass[l],sep=" "))
                plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",phenoDataClass[k]," vs ",phenoDataClass[l], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
                colvect<-NULL
                for (j in 1:length(classnames[[k]])) {
                    tic <- TIC[[classnames[[k]][j]]]
                    # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
                    points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[k]][j]])
                }
                for (j in 1:length(classnames[[l]])) {
                    # i=class2names[j]
                    tic <- TIC[[classnames[[l]][j]]]
                    points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[l]][j]])
                }
                legend("topright",paste(basename(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)
            }
        }
    }#end if length >2
    if (length(phenoDataClass)==2){
        k=1
        l=2

        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",phenoDataClass[k],"vs",phenoDataClass[l], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
        colvect<-NULL
        for (j in 1:length(classnames[[k]])) {
            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }
        for (j in 1:length(classnames[[l]])) {
            # i=class2names[j]
            tic <- TIC[[classnames[[l]][j]]]
            points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[l]][j]])
        }
        legend("topright",paste(basename(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==2

    #case where only one class
    if (length(phenoDataClass)==1){
        k=1
        ylim = range(sapply(TIC, function(x) range(x[,2])))

        #plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",phenoDataClass[k], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = "Raw Total Ion Chromatograms", xlab = "Retention Time (min)", ylab = "TIC")
        colvect<-NULL
        for (j in 1:length(classnames[[k]])) {
            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }

        legend("topright",paste(basename(files[c(classnames[[k]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==1

    dev.off() #pdf(pdfname,w=16,h=10)

    invisible(TIC)
}


##
##  overlay TIC from all files in current folder or from xcmsSet, create pdf
##
#@author Y. Guitton
getTICs2 <- function(xcmsSet=NULL,files=NULL, pdfname="TICs.pdf",rt=c("raw","corrected")) {
    cat("Creating TIC pdf...\n")

    if (is.null(xcmsSet)) {
        filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
        filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
        if (is.null(files))
            files <- getwd()
        info <- file.info(files)
        listed <- list.files(files[info$isdir], pattern = filepattern, recursive = TRUE, full.names = TRUE)
        files <- c(files[!info$isdir], listed)
    } else {
        files <- filepaths(xcmsSet)
    }

    phenoDataClass<-as.vector(levels(xcmsSet@phenoData[,1])) #sometime phenoData have more than 1 column use first as class
    classnames<-vector("list",length(phenoDataClass))
    for (i in 1:length(phenoDataClass)){
        classnames[[i]]<-which( xcmsSet@phenoData[,1]==phenoDataClass[i])
    }

    N <- length(files)
    TIC <- vector("list",N)

    for (i in 1:N) {
        if (!is.null(xcmsSet) && rt == "corrected")
            rtcor <- xcmsSet@rt$corrected[[i]] else
        rtcor <- NULL
        TIC[[i]] <- getTIC(files[i],rtcor=rtcor)
    }

    pdf(pdfname,w=16,h=10)
    cols <- rainbow(N)
    lty = 1:N
    pch = 1:N
    #search for max x and max y in TICs
    xlim = range(sapply(TIC, function(x) range(x[,1])))
    ylim = range(sapply(TIC, function(x) range(x[,2])))
    ylim = c(-ylim[2], ylim[2])


    ##plot start
    if (length(phenoDataClass)>2){
        for (k in 1:(length(phenoDataClass)-1)){
            for (l in (k+1):length(phenoDataClass)){
                #print(paste(phenoDataClass[k],"vs",phenoDataClass[l],sep=" "))
                plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",phenoDataClass[k]," vs ",phenoDataClass[l], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
                colvect<-NULL
                for (j in 1:length(classnames[[k]])) {
                    tic <- TIC[[classnames[[k]][j]]]
                    # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
                    points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[k]][j]])
                }
                for (j in 1:length(classnames[[l]])) {
                    # i=class2names[j]
                    tic <- TIC[[classnames[[l]][j]]]
                    points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[l]][j]])
                }
                legend("topright",paste(basename(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)
            }
        }
    }#end if length >2
    if (length(phenoDataClass)==2){
        k=1
        l=2

        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",phenoDataClass[k],"vs",phenoDataClass[l], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
        colvect<-NULL
        for (j in 1:length(classnames[[k]])) {
            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }
        for (j in 1:length(classnames[[l]])) {
            # i=class2names[j]
            tic <- TIC[[classnames[[l]][j]]]
            points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[l]][j]])
        }
        legend("topright",paste(basename(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==2

    #case where only one class
    if (length(phenoDataClass)==1){
        k=1
        ylim = range(sapply(TIC, function(x) range(x[,2])))

        #plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",phenoDataClass[k], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = "RT Corrected Total Ion Chromatograms", xlab = "Retention Time (min)", ylab = "TIC")
        colvect<-NULL
        for (j in 1:length(classnames[[k]])) {
            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }

        legend("topright",paste(basename(files[c(classnames[[k]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==1

    dev.off() #pdf(pdfname,w=16,h=10)

    invisible(TIC)
}



#@author Y. Guitton
getBPC <- function(file,rtcor=NULL, ...) {
    object <- xcmsRaw(file)
    sel <- profRange(object, ...)
    cbind(if (is.null(rtcor)) object@scantime[sel$scanidx] else rtcor ,xcms:::colMax(object@env$profile[sel$massidx,sel$scanidx,drop=FALSE]))
    #plotChrom(xcmsRaw(file), base=T)
}

#@author Y. Guitton
getBPCs <- function (xcmsSet=NULL, pdfname="BPCs.pdf",rt=c("raw","corrected"), scanrange=NULL) {
    cat("Creating BIC pdf...\n")

    if (is.null(xcmsSet)) {
        cat("Enter an xcmsSet \n")
        stop()
    } else {
        files <- filepaths(xcmsSet)
    }

    phenoDataClass<-as.vector(levels(xcmsSet@phenoData[,1])) #sometime phenoData have more than 1 column use first as class

    classnames<-vector("list",length(phenoDataClass))
    for (i in 1:length(phenoDataClass)){
        classnames[[i]]<-which( xcmsSet@phenoData[,1]==phenoDataClass[i])
    }

    N <- dim(phenoData(xcmsSet))[1]

    TIC <- vector("list",N)


    for (j in 1:N) {

        TIC[[j]] <- getBPC(files[j])
        #good for raw
        # seems strange for corrected
        #errors if scanrange used in xcmsSetgeneration
        if (!is.null(xcmsSet) && rt == "corrected")
            rtcor <- xcmsSet@rt$corrected[[j]]
        else
            rtcor <- NULL

        TIC[[j]] <- getBPC(files[j],rtcor=rtcor)
        # TIC[[j]][,1]<-rtcor
    }



    pdf(pdfname,w=16,h=10)
    cols <- rainbow(N)
    lty = 1:N
    pch = 1:N
    #search for max x and max y in BPCs
    xlim = range(sapply(TIC, function(x) range(x[,1])))
    ylim = range(sapply(TIC, function(x) range(x[,2])))
    ylim = c(-ylim[2], ylim[2])


    ##plot start

    if (length(phenoDataClass)>2){
        for (k in 1:(length(phenoDataClass)-1)){
            for (l in (k+1):length(phenoDataClass)){
                #print(paste(phenoDataClass[k],"vs",phenoDataClass[l],sep=" "))
                plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",phenoDataClass[k]," vs ",phenoDataClass[l], sep=""), xlab = "Retention Time (min)", ylab = "BPC")
                colvect<-NULL
                for (j in 1:length(classnames[[k]])) {
                    tic <- TIC[[classnames[[k]][j]]]
                    # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
                    points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[k]][j]])
                }
                for (j in 1:length(classnames[[l]])) {
                    # i=class2names[j]
                    tic <- TIC[[classnames[[l]][j]]]
                    points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[l]][j]])
                }
                legend("topright",paste(basename(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)
            }
        }
    }#end if length >2

    if (length(phenoDataClass)==2){
        k=1
        l=2
        colvect<-NULL
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",phenoDataClass[k],"vs",phenoDataClass[l], sep=""), xlab = "Retention Time (min)", ylab = "BPC")

        for (j in 1:length(classnames[[k]])) {

            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }
        for (j in 1:length(classnames[[l]])) {
            # i=class2names[j]
            tic <- TIC[[classnames[[l]][j]]]
            points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[l]][j]])
        }
        legend("topright",paste(basename(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==2

    #case where only one class
    if (length(phenoDataClass)==1){
        k=1
        ylim = range(sapply(TIC, function(x) range(x[,2])))
        colvect<-NULL
        #plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",phenoDataClass[k], sep=""), xlab = "Retention Time (min)", ylab = "BPC")
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = "Raw Base Peak Chromatograms", xlab = "Retention Time (min)", ylab = "BPC")
        for (j in 1:length(classnames[[k]])) {
            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }

        legend("topright",paste(basename(files[c(classnames[[k]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==1

    dev.off() #pdf(pdfname,w=16,h=10)

    invisible(TIC)
}




#@author Y. Guitton
getBPCs2 <- function (xcmsSet=NULL, pdfname="BPCs.pdf",rt=c("raw","corrected"), scanrange=NULL) {
    cat("Creating BIC pdf...\n")

    if (is.null(xcmsSet)) {
        cat("Enter an xcmsSet \n")
        stop()
    } else {
        files <- filepaths(xcmsSet)
    }

    phenoDataClass<-as.vector(levels(xcmsSet@phenoData[,1])) #sometime phenoData have more than 1 column use first as class

    classnames<-vector("list",length(phenoDataClass))
    for (i in 1:length(phenoDataClass)){
        classnames[[i]]<-which( xcmsSet@phenoData[,1]==phenoDataClass[i])
    }

    N <- dim(phenoData(xcmsSet))[1]

    TIC <- vector("list",N)


    for (j in 1:N) {

        TIC[[j]] <- getBPC(files[j])
        #good for raw
        # seems strange for corrected
        #errors if scanrange used in xcmsSetgeneration
        if (!is.null(xcmsSet) && rt == "corrected")
            rtcor <- xcmsSet@rt$corrected[[j]]
        else
            rtcor <- NULL

        TIC[[j]] <- getBPC(files[j],rtcor=rtcor)
        # TIC[[j]][,1]<-rtcor
    }



    pdf(pdfname,w=16,h=10)
    cols <- rainbow(N)
    lty = 1:N
    pch = 1:N
    #search for max x and max y in BPCs
    xlim = range(sapply(TIC, function(x) range(x[,1])))
    ylim = range(sapply(TIC, function(x) range(x[,2])))
    ylim = c(-ylim[2], ylim[2])


    ##plot start

    if (length(phenoDataClass)>2){
        for (k in 1:(length(phenoDataClass)-1)){
            for (l in (k+1):length(phenoDataClass)){
                #print(paste(phenoDataClass[k],"vs",phenoDataClass[l],sep=" "))
                plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",phenoDataClass[k]," vs ",phenoDataClass[l], sep=""), xlab = "Retention Time (min)", ylab = "BPC")
                colvect<-NULL
                for (j in 1:length(classnames[[k]])) {
                    tic <- TIC[[classnames[[k]][j]]]
                    # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
                    points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[k]][j]])
                }
                for (j in 1:length(classnames[[l]])) {
                    # i=class2names[j]
                    tic <- TIC[[classnames[[l]][j]]]
                    points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[l]][j]])
                }
                legend("topright",paste(basename(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)
            }
        }
    }#end if length >2

    if (length(phenoDataClass)==2){
        k=1
        l=2
        colvect<-NULL
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",phenoDataClass[k],"vs",phenoDataClass[l], sep=""), xlab = "Retention Time (min)", ylab = "BPC")

        for (j in 1:length(classnames[[k]])) {

            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }
        for (j in 1:length(classnames[[l]])) {
            # i=class2names[j]
            tic <- TIC[[classnames[[l]][j]]]
            points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[l]][j]])
        }
        legend("topright",paste(basename(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==2

    #case where only one class
    if (length(phenoDataClass)==1){
        k=1
        ylim = range(sapply(TIC, function(x) range(x[,2])))
        colvect<-NULL
        #plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",phenoDataClass[k], sep=""), xlab = "Retention Time (min)", ylab = "BPC")
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = "RT Corrected Base Peak Chromatograms", xlab = "Retention Time (min)", ylab = "BPC")
        for (j in 1:length(classnames[[k]])) {
            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }

        legend("topright",paste(basename(files[c(classnames[[k]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==1

    dev.off() #pdf(pdfname,w=16,h=10)

    invisible(TIC)
}


getDeviation <- function (xcmsSet=NULL, pdfname="RT_Deviation_Plot.pdf", scanrange=NULL) {
    cat("Creating RT Devition pdf...\n")

    if (is.null(xcmsSet)) {
        cat("Enter an xcmsSet \n")
        stop()
    } else {
        files <- filepaths(xcmsSet)
    }

    phenoDataClass<-as.vector(levels(xcmsSet@phenoData[,1])) #sometime phenoData have more than 1 column use first as class

    classnames<-vector("list",length(phenoDataClass))
    for (i in 1:length(phenoDataClass)){
        classnames[[i]]<-which( xcmsSet@phenoData[,1]==phenoDataClass[i])
    }

    N <- length(files)
    deviation <- vector("list",N)

    pdf(pdfname,w=16,h=10)
    cols <- rainbow(N)
    lty = 1:N
    pch = 1:N

    
    for (i in 1:N) {
        rtcor <- xcmsSet@rt$corrected[[i]]
	rt <- xcmsSet@rt$raw[[i]]
        dev = rt - rtcor
	deviation[[i]] = cbind( rtcor, dev  )
    }

    xlim = range(sapply(deviation, function(x) range(x[,1])))
    ylim = range(sapply(deviation, function(x) range(x[,2])))


    plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, xlab = "Corrected Retention Time (min)", ylab = "RT-RTcor (s)")
    k=1
    colvect<-NULL
    for (j in 1:length(classnames[[k]])) {
	dev =  deviation[[classnames[[k]][j]]]
	points(dev[,1]/60, dev[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
	colvect<-append(colvect,cols[classnames[[k]][j]])
    }
    legend("topright",paste(basename(files[c(classnames[[k]])])), col = colvect, lty = lty, pch = pch)
    dev.off()
}


#@author Y. Guitton
mygetEIC <- function(file,rtcor=NULL, mzrange=NULL) {
    object <- xcmsRaw(file)
    cbind(if (is.null(rtcor)) object@scantime else rtcor, rawEIC(object,mzrange=mzrange)$intensity)
}

getEICs <- function (xcmsSet=NULL, pt=NULL, scanrange=NULL) {
    if (is.null(xcmsSet)) {
        cat("Enter an xcmsSet \n")
        stop()
    } 

    N <- dim(phenoData(xcmsSet))[1]
    cols <- rainbow(N)
    lty = 1:N
    pch = 1:N

    mzrange = cbind( pt$mzmin, pt$mzmax )
    rtrange = cbind( pt$rt-60, pt$rt+60 )
    xcmsEIC = getEIC(xcmsSet, mzrange = mzrange, rtrange = rtrange, rt = "corrected")
    eic = xcmsEIC@eic
   
    M <- dim(xcmsEIC@mzrange)[1]
    sampnames = sampnames(xcmsEIC)

    for( i in 1:M){

        xlim = xcmsEIC@rtrange[i,];

	pdfname = paste("EICs/",  "LC", i , ".pdf", sep="");
	pdf(pdfname,w=16,h=10)
	
	ymin = min(eic[[1]][[i]][,2])
	ymax = max(eic[[1]][[i]][,2])
	for(j in 2:length(sampnames)){
		ymin = min( c( ymin, min(eic[[j]][[i]][,2]) ) )
		ymax = max( c( ymax, max(eic[[j]][[i]][,2]) ) )
	}
	
	ylim = c( ymin, ymax );

	plot(0, 0, type="n", xlim = xlim, ylim = ylim, main = paste( "LC", i , " extracted ion chromatogram for m/z values: ", round(pt$mzmin[i], 3), "-", round(pt$mzmax[i], 3) ,"\nrtmin:", round(pt$rtmin[i], 1), "  rt:", round(pt$rt[i], 1),   "  rtmax:", round(pt$rtmax[i], 1), sep=""), xlab = "Retention Time (s)", ylab = "EIC")

	colvect<-NULL
        for ( j in 1:length(sampnames) ) {

            points( eic[[j]][[i]][,1], eic[[j]][[i]][,2], col = cols[j], pch = pch[j], type="l")
            colvect<-append(colvect,cols[j])
        }
	abline( v = pt$rt[i], col = "black", lwd = "2", lty=2)
	abline( v = pt$rtmin[i], col = "black", lwd = "2" ,lty=2)
	abline( v = pt$rtmax[i], col = "black", lwd = "2", lty=2)
	legend("topright",paste( sampnames), col = colvect, lty = lty, pch = pch)
	dev.off()
    }
}


getEICs_v2 <- function (xcmsSet=NULL, pt=NULL, sortorder=sortorder ) {
	groupidx = 1:nrow(xcmsSet@groups)
	xeic.raw <- getEIC(xcmsSet, rt = "raw", groupidx= groupidx)
	xeic.corrected <- getEIC(xcmsSet, rt = "corrected", groupidx= groupidx)
	sampnames = names(xeic.raw@eic)
	N = length(sampnames)
	cols <- rainbow(N)
	lty = 1:N
	pch = 1:N

	raw_eic = xeic.raw@eic
	cor_eic = xeic.corrected@eic

	M <- dim(xeic.raw@mzrange)[1]

	for ( i in groupidx){
		i_new = which(sortorder==i)
		pdfname = paste("EICs/",  "LC", i_new , ".pdf", sep="");
		pdf(pdfname,w=16,h=12)
		layout(matrix(1:2,2,1))

		xlim = xeic.raw@rtrange[i,];
		ymin = min(raw_eic[[1]][[i]][,2])
		ymax = max(raw_eic[[1]][[i]][,2])
		for(j in 2:length(sampnames)){
			ymin = min( c( ymin, min(raw_eic[[j]][[i]][,2]) ) )
			ymax = max( c( ymax, max(raw_eic[[j]][[i]][,2]) ) )
		}
		ylim = c( ymin, ymax );

		plot(0, 0, type="n", xlim = xlim, ylim = ylim, main = paste( "LC", i_new , " extracted ion chromatogram for m/z values: ", round(xeic.raw@mzrange[i,1], 3), "-", round(xeic.raw@mzrange[i,2], 3), sep=""), xlab = "Retention Time (s)", ylab = "EIC")

		colvect<-NULL
		for ( j in 1:length(sampnames) ) {
		    points( raw_eic[[j]][[i]][,1], raw_eic[[j]][[i]][,2], col = cols[j], pch = pch[j], type="l")
		    colvect<-append(colvect,cols[j])
		}
		legend("topright",paste( sampnames), col = colvect, lty = lty, pch = pch)


		xlim = xeic.corrected@rtrange[i,];
		ymin = min(cor_eic[[1]][[i]][,2])
		ymax = max(cor_eic[[1]][[i]][,2])
		for(j in 2:length(sampnames)){
			ymin = min( c( ymin, min(cor_eic[[j]][[i]][,2]) ) )
			ymax = max( c( ymax, max(cor_eic[[j]][[i]][,2]) ) )
		}
		ylim = c( ymin, ymax );

		plot(0, 0, type="n", xlim = xlim, ylim = ylim, main = paste( "LC", i_new , " extracted ion chromatogram for m/z values: ", round(xeic.corrected@mzrange[i,1], 3), "-", round(xeic.corrected@mzrange[i,2], 3), sep=""), xlab = " Corrected Retention Time (s)", ylab = "EIC")

		colvect<-NULL
		for ( j in 1:length(sampnames) ) {
		    points( cor_eic[[j]][[i]][,1], cor_eic[[j]][[i]][,2], col = cols[j], pch = pch[j], type="l")
		    colvect<-append(colvect,cols[j])
		}
		abline( v = pt$rt[i], col = "black", lwd = "2", lty=2)
		legend("topright",paste( sampnames), col = colvect, lty = lty, pch = pch)

		dev.off()
	}
}

