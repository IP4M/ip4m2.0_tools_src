
##ADDITIONS FROM Y. Guitton
getBPC <- function(file,rtcor=NULL, ...) {
    object <- xcmsRaw(file)
    sel <- profRange(object, ...)
    cbind(if (is.null(rtcor)) object@scantime[sel$scanidx] else rtcor ,xcms:::colMax(object@env$profile[sel$massidx,sel$scanidx,drop=FALSE]))
}

getBPC2s <- function (files, xset = NULL, pdfname="BPCs.pdf", rt = c("raw","corrected"), scanrange=NULL) {
    require(xcms)

    #create sampleMetadata, get sampleMetadata and class
    if(!is.null(xset)) {
    	#When files come from XCMS3 directly before metaMS
    	sampleMetadata <- xset@phenoData
    } else {
    	#When files come from a zip directory with raw files before metaMS
	    sampleMetadata<-xcms:::phenoDataFromPaths(files)
    }
    class<-unique(sampleMetadata[,"class"]) #create phenoData like table
    classnames<-vector("list",length(class))
    for (i in 1:length(class)){
        classnames[[i]]<-which( sampleMetadata[,"class"]==class[i])
    }

    N <- dim(sampleMetadata)[1]
    BPC <- vector("list",N)

    for (j in 1:N) {
        BPC[[j]] <- getBPC(files[j])
        #good for raw 
        # seems strange for corrected
        #errors if scanrange used in xcmsSetgeneration
        if (!is.null(xcmsSet) && rt == "corrected"){
            rtcor <- xcmsSet@rt$corrected[[j]] 
        }else{
            rtcor <- NULL
        }
        BPC[[j]] <- getBPC(files[j],rtcor=rtcor)
    }

    pdf(pdfname,w=16,h=10)
    cols <- rainbow(N)
    lty = 1:N
    pch = 1:N
    #search for max x and max y in BPCs

    xlim = range(sapply(BPC, function(x) range(x[,1])))
    ylim = range(sapply(BPC, function(x) range(x[,2])))

    ylim = c(-ylim[2], ylim[2])

    ##plot start
    if (length(class)>2){
        for (k in 1:(length(class)-1)){
            for (l in (k+1):length(class)){
                cat(paste(class[k],"vs",class[l],sep=" ","\n")) 
                plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",class[k]," vs ",class[l], sep=""), xlab = "Retention Time (min)", ylab = "BPC")
                colvect<-NULL
                for (j in 1:length(classnames[[k]])) {
                    bpc <- BPC[[classnames[[k]][j]]]
                    # points(bpc[,1]/60, bpc[,2], col = cols[i], pch = pch[i], type="l")
                    points(bpc[,1]/60, bpc[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[k]][j]])
                }
                for (j in 1:length(classnames[[l]])) {
                    # i=class2names[j]
                    bpc <- BPC[[classnames[[l]][j]]]
                    points(bpc[,1]/60, -bpc[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[l]][j]])
                }
                legend("topright",paste(gsub("(^.+)\\..*$","\\1",basename(files[c(classnames[[k]],classnames[[l]])]))), col = colvect, lty = lty, pch = pch)
            }
        }
    }#end if length >2

    if (length(class)==2){
        k=1
		l=2
        colvect<-NULL
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",class[k],"vs",class[l], sep=""), xlab = "Retention Time (min)", ylab = "BPC")

        for (j in 1:length(classnames[[k]])) {
            bpc <- BPC[[classnames[[k]][j]]]
            # points(bpc[,1]/60, bpc[,2], col = cols[i], pch = pch[i], type="l")
            points(bpc[,1]/60, bpc[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }
        for (j in 1:length(classnames[[l]])) {
            # i=class2names[j]
            bpc <- BPC[[classnames[[l]][j]]]
            points(bpc[,1]/60, -bpc[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[l]][j]])
        }
        legend("topright",paste(gsub("(^.+)\\..*$","\\1",basename(files[c(classnames[[k]],classnames[[l]])]))), col = colvect, lty = lty, pch = pch)
    }#end length ==2
    
    if (length(class)==1){
        k=1

		ylim = range(sapply(BPC, function(x) range(x[,2])))

        colvect<-NULL
        #plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",class[k], sep=""), xlab = "Retention Time (min)", ylab = "BPC")
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = "Base Peak Chromatograms", xlab = "Retention Time (min)", ylab = "BPC")
        for (j in 1:length(classnames[[k]])) {
            bpc <- BPC[[classnames[[k]][j]]]
            # points(bpc[,1]/60, bpc[,2], col = cols[i], pch = pch[i], type="l")
            points(bpc[,1]/60, bpc[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }
        legend("topright",paste(gsub("(^.+)\\..*$","\\1",basename(files[c(classnames[[k]])]))), col = colvect, lty = lty, pch = pch)
    }#end length ==1
    dev.off()
}

getTIC <- function(file,rtcor=NULL) {
    object <- xcmsRaw(file)
    cbind(if (is.null(rtcor)) object@scantime else rtcor, rawEIC(object,mzrange=range(object@env$mz))$intensity)
}

##  overlay TIC from all files in current folder or from xcmsSet, create pdf
getTIC2s <- function(files, xset=NULL, pdfname="TICs.pdf", rt=c("raw","corrected")) {
    require(xcms)

    #create sampleMetadata, get sampleMetadata and class
    if(!is.null(xset)){
		#When files come from XCMS3 before metaMS treatment
    	sampleMetadata<-xset@phenoData
    } else {
		#When files come from a zip directory with raw files before metaMS
    	sampleMetadata<-xcms:::phenoDataFromPaths(files)
    }
    class<-as.vector(levels(sampleMetadata[,"class"])) #create phenoData like table
    classnames<-vector("list",length(class))
    for (i in 1:length(class)){
        classnames[[i]]<-which( sampleMetadata[,"class"]==class[i])
    }
        
    N <- dim(sampleMetadata)[1]
    TIC <- vector("list",N)

    for (i in 1:N) {
        if (!is.null(xcmsSet) && rt == "corrected")
            rtcor <- xcmsSet@rt$corrected[[i]]
        else
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
    if (length(class)>2){
        for (k in 1:(length(class)-1)){
            for (l in (k+1):length(class)){
                cat(paste(class[k],"vs",class[l],"\n",sep=" ")) 
                plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",class[k]," vs ",class[l], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
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
                legend("topright",paste(gsub("(^.+)\\..*$","\\1",basename(files[c(classnames[[k]],classnames[[l]])]))), col = colvect, lty = lty, pch = pch)
            }
        }
    }#end if length >2

    if (length(class)==2){
        k=1
        l=2
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",class[k],"vs",class[l], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
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
        legend("topright",paste(gsub("(^.+)\\..*$","\\1",basename(files[c(classnames[[k]],classnames[[l]])]))), col = colvect, lty = lty, pch = pch)
    }#end length ==2

    if (length(class)==1){
        k=1
        ylim = range(sapply(TIC, function(x) range(x[,2])))
        #plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",class[k], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = "Total Ion Chromatograms", xlab = "Retention Time (min)", ylab = "TIC")
        colvect<-NULL
        for (j in 1:length(classnames[[k]])) {
            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }
        legend("topright",paste(gsub("(^.+)\\..*$","\\1",basename(files[c(classnames[[k]])]))), col = colvect, lty = lty, pch = pch)
    }#end length ==1
    dev.off()
}



#Update by J.Saint-Vanne
##addition for quality control of peak picking
#metaMS EIC and pspectra plotting option
#version 20190520
#only for Galaxy 
plotUnknowns<-function(resGC, unkn="", DB=NULL, fileFrom=NULL){

    ##Annotation table each value is a pcgrp associated to the unknown 
    ##NOTE pcgrp index are different between xcmsSet and resGC due to filtering steps in metaMS
    ##R. Wehrens give me some clues on that and we found a correction
    
	#correction of annotation matrix due to pcgrp removal by quality check in runGCresult
	#matrix of correspondance between an@pspectra and filtered pspectra from runGC
    #Select only pspectra which correpond to them select in metaMS
    # col1 = filtered spectra from runGC and col2 = an@spectra
	allPCGRPs <-lapply(1:length(resGC$xset),
				function(i) {
					an <- resGC$xset[[i]]
					huhn <- an@pspectra[which(sapply(an@pspectra, length) >=
					metaSetting(resGC$settings,
					"DBconstruction.minfeat"))]
					matCORR<-cbind(1:length(huhn), match(huhn, an@pspectra))
				})
    #Build a new annotation list with sampnames and pseudospectra number from xset
    helpannotation <- list()
    for(j in 1:length(resGC$xset)){
        helpannotation[[j]] <- resGC$annotation[[j]][1:2]
        pspvector <- vector()
        for(i in 1: nrow(helpannotation[[j]])){
            #Find corresponding pspec
            psplink <- allPCGRPs[[j]][match(helpannotation[[j]][i,1],allPCGRPs[[j]]),2]
            pspvector <- c(pspvector,psplink)
            #Change the annotation column into sampname column
            if(helpannotation[[j]][i,2] < 0){
                #It's an unknown
                new_name <- paste("Unknown",abs(as.integer(helpannotation[[j]][i,2])))
                helpannotation[[j]][i,2] <- new_name
            }else{
                #It has been found in local database
                for(k in 1:length(DB)){
                    if(helpannotation[[j]][i,2] == k){
                        helpannotation[[j]][i,2] <- DB[[k]]$Name
                        break
                    }
                }
            }
        }
        helpannotation[[j]] <- cbind(helpannotation[[j]],pspvector)
        names(helpannotation)[j] <- names(resGC$annotation[j])
    }
    peaktable <- resGC$PeakTable
		
	par (mar=c(5, 4, 4, 2) + 0.1)
	#For each unknown
	for (l in 1:length(unkn)){
		#recordPlot
		perpage=3 #if change change layout also!
		dev.new(width=21/2.54, height=29.7/2.54, file=paste("EICs/MC",unkn[l],".pdf", sep="")) #A4 pdf
		# par(mfrow=c(perpage,2))
		layout(matrix(c(1,1,2,3,4,4,5,6,7,7,8,9), 6, 2, byrow = TRUE), widths=rep(c(1,1),perpage), heights=rep(c(1,5),perpage))
		# layout.show(6)
		oma.saved <- par("oma")
		par(oma = rep.int(0, 4))
		par(oma = oma.saved)
		o.par <- par(mar = rep.int(0, 4))
		on.exit(par(o.par))
		#For each sample
		for (c in 1:length(resGC$xset)) {	
			#get sample name
			sampname<-basename(resGC$xset[[c]]@xcmsSet@filepaths)
			#remove .cdf, .mzXML filepattern
			filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
									"[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
			filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
			sampname<-gsub(filepattern, "",sampname)					 
			title1<-paste("MC", unkn[l], "in",sampname, sep = " ")
			an<-resGC$xset[[c]]
    		if(fileFrom == "zipfile") {
				an@xcmsSet@filepaths <- paste0("./",an@xcmsSet@phenoData[,"class"],"/",basename(an@xcmsSet@filepaths))
				}#else {
        			#print(an@xcmsSet@filepaths)
					#an@xcmsSet@filepaths <- paste0("./",basename(an@xcmsSet@filepaths))
				#}
			#Find the good annotation for this sample
            for(a in 1:length(helpannotation)){    
                #if(gsub(filepattern, "", names(helpannotation)[a]) == paste0("./",sampname)){
		if(gsub(filepattern, "", names(helpannotation)[a]) == paste0("",sampname)){
                    #Find the unkn or the matched std in this sample
                    findunkn <- FALSE
                    for(r in 1:nrow(helpannotation[[a]])){
                        if(helpannotation[[a]][r,"annotation"] == peaktable[unkn[l],1]){
                            findunkn <- TRUE
                            pcgrp <- helpannotation[[a]][r,"pspvector"]
							par (mar=c(0, 0, 0, 0) + 0.1)
							#Write title
							plot.new()
							box()
							text(0.5, 0.5, title1, cex=2)						
							par (mar=c(3, 2.5, 3, 1.5) + 0.1)
							#Window for EIC
							plotEICs(an, pspec=pcgrp, maxlabel=0)
							#Window for pseudospectra
							plotPsSpectrum(an, pspec=pcgrp, maxlabel=3)
						}
					}
					if(!findunkn) {
						par (mar=c(0, 0, 0, 0) + 0.1)
                   		#Write title
						plot.new()
						box()
						text(0.5, 0.5, title1, cex=2)
                   		#Window for EIC
                   		plot.new()
                   		box()
						text(0.5, 0.5, "NOT FOUND", cex=2)
						#Window for pseudospectra
						plot.new()
						box()
						text(0.5, 0.5, "NOT FOUND", cex=2)
					}
					break
				}
    		}
		}
		graphics.off()
	}#end  for unkn[l]
}#end function
