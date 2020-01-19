# Title     : TODO
# Objective : TODO
# Created by: Administrator
# Created on: 2019/8/26
library(tidyverse)
library(data.table)
library(igraph)
library(xlsx)
library(ggrepel)
library(Rgraphviz)
library(Cairo)

.on.public.web <- F

.set.mSet <- function(mSetObj = NA) {
  if (.on.public.web) {
    mSet <<- mSetObj;
    return(1);
  }
  return(mSetObj);
}

.get.mSet <- function(mSetObj = NA) {
  if (.on.public.web) {
    return(mSet)
  }else {
    return(mSetObj);
  }
}

Setup.MapData <- function(mSetObj = NA, nm.mat)
  {
  mSetObj <- .get.mSet(mSetObj)
  mSetObj$dataSet$cmpd <- nm.mat$kegg
  mSetObj$dataSet$nm.mat <- nm.mat
  return(.set.mSet(mSetObj))
}

.read.metaboanalyst.lib <- function(filenm, dataDir) {
  lib.path <- paste(dataDir, "/", filenm, sep = "");
  return(readRDS(lib.path));
}

MetaboliteMappingExact <- function(mSetObj = NA, q.type)
  {
  mSetObj <- .get.mSet(mSetObj)
  qvec <- mSetObj$dataSet$cmpd
  hit.inx <- vector(mode = "numeric", length = length(qvec))
  names(hit.inx) <- qvec
  gc()
  # mSetObj$name.map$query.vec <- qvec
  # mSetObj$name.map$hit.inx <- hit.inx
  # mSetObj$name.map$hit.values <- match.values
  # mSetObj$name.map$match.state <- match.state
  print(mSetObj$name.map)
  return(.set.mSet(mSetObj))
}

CrossReferencing <- function(mSetObj = NA, q.type, hmdb = T, pubchem = T, chebi = F,
                             kegg = T, metlin = F)
  {
  mSetObj <- .get.mSet(mSetObj)
  mSetObj$return.cols <- c(hmdb, pubchem, chebi, kegg, metlin)
  if (!exists("name.map", where = mSetObj)) {
    mSetObj$name.map <- list()
  }
  mSetObj$dataSet$q.type <- q.type
  if (.on.public.web) {
    .set.mSet(mSetObj)
    MetaboliteMappingExact(mSetObj, q.type)
    mSetObj <- .get.mSet(mSetObj)
  }
  else {
    mSetObj <- MetaboliteMappingExact(mSetObj, q.type)
  }
  todo.inx <- which(is.na(mSetObj$name.map$hit.inx))
  if (length(todo.inx) > 15) {
    mSetObj$msgSet$nmcheck.msg <- c(2, "There are >15 compounds without matches. You can either proceed or if necessary, update these compound IDs and upload again.")
  }
  else {
    mSetObj$msgSet$nmcheck.msg <- c(1, "Name matching OK, please inspect (and manual correct) the results then proceed.")
  }
  return(.set.mSet(mSetObj))
}

InitDataObjects <- function(data.type, anal.type, paired = FALSE)
  {
  dataSet <- list()
  dataSet$type <- data.type
  dataSet$design.type <- "regular"
  dataSet$cls.type <- "disc"
  dataSet$format <- "rowu"
  dataSet$paired <- paired
  analSet <- list()
  analSet$type <- anal.type
  mSetObj <- list()
  mSetObj$dataSet <- dataSet
  mSetObj$analSet <- analSet
  mSetObj$imgSet <- list()
  mSetObj$msgSet <- list()
  mSetObj$msgSet$msg.vec <- vector(mode = "character")
  mSetObj$cmdSet <- vector(mode = "character")
  msg.vec <<- ""
  module.count <<- 0
  smpdbpw.count <<- 0
  data.org <<- NULL
  if (.on.public.web) {
    lib.path <<- "../../data/"
    library(BiocParallel)
    register(SerialParam())
  }
  else {
    lib.path <<- "https://www.metaboanalyst.ca/resources/data/"
  }
  peakFormat <<- "mpt"
  mdata.all <<- list()
  mdata.siggenes <<- vector("list")
  meta.selected <<- TRUE
  anal.type <<- anal.type
  Cairo::CairoFonts(regular = "Arial:style=Regular", bold = "Arial:style=Bold",
                    italic = "Arial:style=Italic", bolditalic = "Arial:style=Bold Italic",
                    symbol = "Symbol")
  print("MetaboAnalyst R objects initialized ...")
  return(.set.mSet(mSetObj))
}

.load.metaboanalyst.lib <- function(libtype, libname, dataDir) {
  destfile <- paste(libname, ".rda", sep = "");
  destfile <- paste(dataDir, "/", libtype, "/", libname, ".rda", sep = "");
  load(destfile, .GlobalEnv);
}

SetKEGG.PathLib <- function(mSetObj = NA, kegg.rda, dataDir)
  {
  mSetObj <- .get.mSet(mSetObj)
  mSetObj$msgSet$lib.msg <- paste("Your selected pathway library code is \\textbf{",
                                  kegg.rda, "}(KEGG organisms abbreviation).")
  kegglib <- .load.metaboanalyst.lib("kegg", kegg.rda, dataDir = dataDir)
  mSetObj$pathwaylibtype <- "KEGG"
  if (.on.public.web) {
    .set.mSet(mSetObj)
    return(1)
  }
  return(.set.mSet(mSetObj))
}

SetMetabolomeFilter <- function(mSetObj = NA, TorF)
  {
  mSetObj <- .get.mSet(mSetObj)
  mSetObj$dataSet$use.metabo.filter <- TorF
  return(.set.mSet(mSetObj))
}

GetFinalNameMap <- function(mSetObj = NA)
  {
  mSetObj <- .get.mSet(mSetObj)
  hit.inx <- mSetObj$name.map$hit.inx
  hit.values <- mSetObj$name.map$hit.values
  qvec <- mSetObj$dataSet$cmpd
  nm.mat <-  mSetObj$dataSet$nm.mat
  return(as.data.frame(nm.mat))
}

GetORA.pathNames <- function(mSetObj = NA)
  {
  mSetObj <- .get.mSet(mSetObj)
  hit.inx <- match(rownames(mSetObj$analSet$ora.mat), metpa$path.ids)
  return(names(metpa$path.ids)[hit.inx])
}

SetupKEGGLinks <- function(smpdb.ids)
  {
  kegg.vec <- metpa$path.keggs[match(smpdb.ids, names(metpa$mset.list))]
  GetKEGGLinks(kegg.vec)
}

GetKEGGLinks <- function(kegg.vec)
  {
  lk.len <- length(kegg.vec)
  all.lks <- vector(mode = "character", length = lk.len)
  for (i in 1:lk.len) {
    lks <- strsplit(kegg.vec[i], "; ")[[1]]
    if (!is.na(lks[1])) {
      all.lks[i] <- paste("http://www.genome.jp/kegg-bin/show_pathway?",
                          lks, sep = "")
    }
  }
  return(all.lks)
}

SetupSMPDBLinks <- function(kegg.ids)
  {
  smpdb.vec <- names(metpa$path.smps)[match(kegg.ids, metpa$path.smps)]
  GetSMPDBLinks(smpdb.vec)
}

GetSMPDBLinks <- function(smpdb.vec)
  {
  lk.len <- length(smpdb.vec)
  all.lks <- vector(mode = "character", length = lk.len)
  for (i in 1:lk.len) {
    lks <- strsplit(smpdb.vec[i], "; ")[[1]]
    if (!is.na(lks[1])) {
      all.lks[i] <- paste("http://www.smpdb.ca/view/",
                          lks, sep = "", collapse = " ")
    }
  }
  return(all.lks)
}

AddErrMsg <- function(msg)
  {
  msg.vec <<- ""
  msg.vec <<- c(msg.vec, msg)
  print(msg)
}

GetFisherPvalue <- function(numSigMembers, numSigAll, numMembers, numAllMembers)
  {
  z <- cbind(numSigMembers, numSigAll - numSigMembers, numMembers -
    numSigMembers, numAllMembers - numMembers - numSigAll +
               numSigMembers)
  z <- lapply(split(z, 1:nrow(z)), matrix, ncol = 2)
  z <- lapply(z, fisher.test, alternative = "greater")
  p.values <- as.numeric(unlist(lapply(z, "[[", "p.value"),
                                use.names = FALSE))
  return(p.values)
}

CalculateOraScore <- function(mSetObj = NA, nodeImp, method, dataDir, fcData, hasSmp)
  {
  mSetObj <- .get.mSet(mSetObj)
  nm.map <- GetFinalNameMap(mSetObj)
  nmMap <- nm.map %>%
    as.data.frame()
  if (mSetObj$pathwaylibtype == "KEGG") {
    valid.inx <- !(is.na(nm.map$kegg) | duplicated(nm.map$kegg))
    ora.vec <- nm.map$kegg[valid.inx]
  }
  else if (mSetObj$pathwaylibtype == "SMPDB") {
    valid.inx <- !(is.na(nm.map$hmdbid) | duplicated(nm.map$hmdbid))
    ora.vec <- nm.map$hmdbid[valid.inx]
  }
  q.size <- length(ora.vec)
  if (is.na(ora.vec) || q.size == 0) {
    if (mSetObj$pathwaylibtype == "KEGG") {
      AddErrMsg("No valid KEGG compounds found!")
    }
    else if (mSetObj$pathwaylibtype == "SMPDB") {
      AddErrMsg("No valid SMPDB compounds found!")
    }
    return(0)
  }
  current.mset <- metpa$mset.list
  uniq.count <- metpa$uniq.count
  if (mSetObj$dataSet$use.metabo.filter && !is.null(mSetObj$dataSet$metabo.filter.kegg)) {
    current.mset <- lapply(current.mset, function(x) {
      x[x %in% mSetObj$dataSet$metabo.filter.kegg]
    })
    mSetObj$analSet$ora.filtered.mset <- current.mset
    uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)))
  }
  # print(current.mset)
  print(ora.vec)
  hits <- lapply(current.mset, function(x) {
    x[x %in% ora.vec]
  })
  hit.num <- unlist(lapply(hits, function(x) {
    length(x)
  }), use.names = FALSE)
  set.size <- length(current.mset)
  set.num <- unlist(lapply(current.mset, length), use.names = FALSE)
  # res.mat <- matrix(0, nrow = set.size, ncol = 8)
  res.mat <- tibble(name = names(current.mset))
  # print(names(current.mset))
  # rownames(res.mat) <- names(current.mset)
  # colnames(res.mat) <- c("Total_In_Pathway", "Expected", "Hits", "Raw p",
  # "-log(p)", "Holm Adjusted P", "FDR", "Impact")
  # colnames(res.mat) <- c("name","Total_In_Pathway", "Expected", "Hits", "Raw p",
  # "-log(p)", "Holm Adjusted P", "FDR", "Impact")
  if (nodeImp == "rbc") {
    imp.list <- metpa$rbc
    mSetObj$msgSet$topo.msg <- "Your selected node importance measure for topological analysis is \\textbf{relative betweenness centrality}."
  }
  else if(nodeImp == "dgr") {
    imp.list <- metpa$dgr
    mSetObj$msgSet$topo.msg <- "Your selected node importance measure for topological analysis is \\textbf{out degree centrality}."
  }
  else if(nodeImp == "tdc") {
    imp.list <- metpa$degree_all
    mSetObj$msgSet$topo.msg <- "Your selected node importance measure for topological analysis is \\textbf{total degree centrality}."
  }
  else if(nodeImp == "occ") {
    imp.list <- metpa$outcloseness
    mSetObj$msgSet$topo.msg <- "Your selected node importance measure for topological analysis is \\textbf{out closeness centrality}."
  }
  else if(nodeImp == "icc") {
    imp.list <- metpa$incloseness
    mSetObj$msgSet$topo.msg <- "Your selected node importance measure for topological analysis is \\textbf{in closeness centrality}."
  }
  else if(nodeImp == "tcc") {
    imp.list <- metpa$closeness_all
    mSetObj$msgSet$topo.msg <- "Your selected node importance measure for topological analysis is \\textbf{total closeness centrality}."
  }
  else if(nodeImp == "ec") {
    imp.list <- metpa$evcent
    mSetObj$msgSet$topo.msg <- "Your selected node importance measure for topological analysis is \\textbf{evcent centrality}."
  }
  # res.mat[, 1] <- set.num
  res.mat <- res.mat %>%
    mutate(`Total_In_Pathway` = set.num)
  # res.mat[, 2] <- q.size * (set.num / uniq.count)
  res.mat <- res.mat %>%
    mutate(`Expected` = q.size * (set.num / uniq.count))
  # res.mat[, 3] <- hit.num
  res.mat <- res.mat %>%
    mutate(Hits = hit.num)
  if (method == "fisher") {
    # res.mat[, 4] <- GetFisherPvalue(hit.num, q.size, set.num,
    #                                 uniq.count)
    rawP <- GetFisherPvalue(hit.num, q.size, set.num, uniq.count)
    mSetObj$msgSet$rich.msg <- "The selected over-representation analysis method is \\textbf{Fishers' exact test}."
  }
  else {
    # res.mat[, 4] <- phyper(hit.num - 1, set.num, uniq.count -
    #   set.num, q.size, lower.tail = F)
    rawP <- phyper(hit.num - 1, set.num, uniq.count - set.num, q.size, lower.tail = F)
    mSetObj$msgSet$rich.msg <- "The selected over-representation analysis method is \\textbf{Hypergeometric test}."
  }
  res.mat <- res.mat %>%
    mutate(`Raw p` = rawP)
  # res.mat[, 5] <- -log(res.mat[, 4])
  res.mat <- res.mat %>%
    mutate(`-log(p)` = -log(`Raw p`))
  # res.mat[, 6] <- p.adjust(res.mat[, 4], "holm")
  res.mat <- res.mat %>%
    mutate(`Holm Adjusted P` = p.adjust(`Raw p`, "holm"))
  # res.mat[, 7] <- p.adjust(res.mat[, 4], "fdr")
  res.mat <- res.mat %>%
    mutate(`FDR` = p.adjust(`Raw p`, "fdr"))
  # res.mat[, 8] <- mapply(function(x, y) {
  #   sum(x[y])
  # }, imp.list, hits)
  impact <- mapply(function(x, y) {
    sum(x[y])
  }, imp.list, hits)
  res.mat <- res.mat %>%
    mutate(Impact = impact)
  # res.mat <- res.mat[hit.num > 0, , drop = FALSE]
  res.mat <- res.mat %>%
    filter(Hits > 0)
  # res.mat <- res.mat[!is.na(res.mat[, 8]), , drop = FALSE]
  res.mat <- res.mat %>%
    filter(!is.na(Impact))
  res.mat <- res.mat %>%
    arrange(`Raw p`, Impact)

  print(res.mat)
  # if (nrow(res.mat) > 1) {
  #   ord.inx <- order(res.mat[, 4], res.mat[, 8])
  #   res.mat <- res.mat[ord.inx,]
  # }
  res.mat <- res.mat %>%
    distinct(name, .keep_all = T) %>%
    column_to_rownames("name") %>%
    as.matrix()
  # print(head(res.mat))

  mSetObj$analSet$ora.mat <- signif(res.mat, 5)
  mSetObj$analSet$ora.hits <- hits
  mSetObj$analSet$node.imp <- nodeImp
  .set.mSet(mSetObj)
  save.mat <- mSetObj$analSet$ora.mat
  rownames(save.mat) <- GetORA.pathNames(mSetObj)
  if (mSetObj$pathwaylibtype == "KEGG") {
    hit.inx <- match(rownames(mSetObj$analSet$ora.mat), metpa$path.ids)
    extraKegg <- c()
    redNames <- c()
    setName <- names(metpa$path.ids)[hit.inx]
    for (msetNm in setName) {
      pathid <- metpa$path.ids[msetNm]
      mset <- metpa$mset.list[[pathid]]
      hits <- mSetObj$analSet$ora.hits
      red.inx <- which(mset %in% hits[[pathid]])
      cName = mset[red.inx]
      nms <- names(mset);
      eachRowkeggID <- c()
      for (value in nms[red.inx]) {
        keggId = cName[value]
        eachRowkeggID <- c(eachRowkeggID, keggId)
      }
      eachCol <- fcData %>%
        filter(KEGG %in% eachRowkeggID) %>%
        .$col
      eachRowkegg <- paste("/", eachRowkeggID, "%09", eachCol, collapse = "", sep = "")
      extraKegg <- c(extraKegg, eachRowkegg)
      redName <- paste(sort(nms[red.inx]), collapse = " || ")
      redNames <- c(redNames, redName)
    }
    save.mat <- cbind(save.mat, redNames)
    kegg.vec <- rownames(mSetObj$analSet$ora.mat)
    keggId <- kegg.vec
    keggLink <- paste("http://www.genome.jp/kegg-bin/show_pathway?", kegg.vec, extraKegg, sep = "")

    saveMat <- cbind(save.mat, keggLink)
    save.mat <- saveMat %>%
      as.data.frame(stringsAsFactors = F) %>%
      rownames_to_column("Metabolite")

    if (hasSmp) {
      save.mat <- save.mat %>%
        mutate(keggId = keggId) %>%
        mutate_at(vars("keggId"), function(x) {
          smpdb.vec <- names(metpa$path.smps)[match(x, metpa$path.smps)]
          print(smpdb.vec)
          link <- smpdb.vec %>%
            map_chr(function(smpId) {
              lks <- strsplit(smpId, "; ")[[1]]
              if (!is.na(lks[1])) {
                paste("http://www.smpdb.ca/view/",
                      lks, sep = "", collapse = " ")
              }else ""
            })
          link
        })
    }
  }
  else if (mSetObj$pathwaylibtype == "SMPDB") {
    hit.inx <- match(rownames(mSetObj$analSet$ora.mat), metpa$path.ids)
    extraKegg <- c()
    redNames <- c()
    setName <- names(metpa$path.ids)[hit.inx]

    hmdbIds <- setName %>%
      map_chr(function(x) {
        pathid <- metpa$path.ids[x]
        mset <- metpa$mset.list[[pathid]]
        red.inx <- which(mset %in% hits[[pathid]])
        ids <- mset[red.inx]
        str_c(ids, collapse = ",")
      })
    hmdbDf <- tibble(Name = setName, hmdbId = hmdbIds)
    write_tsv(hmdbDf, "red_hmdb_id.txt")

    for (msetNm in setName) {
      pathid <- metpa$path.ids[msetNm]
      mset <- metpa$mset.list[[pathid]]
      hits <- mSetObj$analSet$ora.hits
      red.inx <- which(mset %in% hits[[pathid]])
      nms <- names(mset)
      eachRowkeggID <- c()
      for (value in nms[red.inx]) {
        keggId = nmMap %>%
          filter(query == value) %>%
          head(1) %>%
          .$kegg %>%
          as.character()
        eachRowkeggID <- c(eachRowkeggID, keggId) %>%
          discard(is.na)
      }
      eachCol <- fcData %>%
        filter(KEGG %in% eachRowkeggID) %>%
        .$col
      eachRowkegg <- if (length(eachRowkeggID) > 0) {
        eachRowkegg <- paste("/", eachRowkeggID, "%09", eachCol, collapse = "", sep = "")
      } else ""
      extraKegg <- c(extraKegg, eachRowkegg)
      redName <- paste(sort(nms[red.inx]), collapse = "\n")
      redNames <- c(redNames, redName)
    }
    save.mat <- cbind(save.mat, redNames)
    smpdb.ids <- rownames(mSetObj$analSet$ora.mat)
    kegg.vec <- metpa$path.keggs[match(smpdb.ids, names(metpa$mset.list))]
    keggLink <- kegg.vec %>%
      map_chr(function(x) {
        unlist(strsplit(x, ";"))[1]
      }) %>%
      imap_chr(function(v, i) {
        extraKegg <- extraKegg[i]
        if (is.na(v)) {
          ""
        }else {
          paste0("http://www.genome.jp/kegg-bin/show_pathway?", v, extraKegg)
        }
      })
    smpdId <- smpdb.ids
    saveMat <- cbind(save.mat, keggLink, smpdId)
    save.mat <- saveMat %>%
      as.data.frame(stringsAsFactors = F) %>%
      rownames_to_column("Metabolite") %>%
      mutate_at(vars("smpdId"), function(x) {
        GetSMPDBLinks(x)
      })
  }

  colNames <- if (hasSmp) {
    c("Metabolite", "Total_In_Pathway", "Expected", "Hits", "Raw P", "-log(p)", "Holm P", "FDR", "Impact",
      "Enriched_Compounds", "KEGGLink", "SMPLink")
  }else {
    c("Metabolite", "Total_In_Pathway", "Expected", "Hits", "Raw P", "-log(p)", "Holm P", "FDR", "Impact",
      "Enriched_Compounds", "KEGGLink")
  }

  save.mat <- save.mat %>%
    set_colnames(colNames) %>%
    column_to_rownames("Metabolite")

  write.csv(save.mat, "Pathway_Result.csv")
  write.table(save.mat, "Pathway_Result.txt", quote = FALSE, sep = "\t")

  if (.on.public.web) {
    return(1)
  }
  else {
    return(.set.mSet(mSetObj))
  }
}

xytrans2 <- function(xy, par) {
  cbind(par[1] + ((par[2] - par[1]) * xy[, 1]),
        par[3] + ((par[4] - par[3]) * xy[, 2]))
}

xytrans <- function(xy, par) {
  cbind((xy[, 1] - par[1]) / (par[2] - par[1]),
        (xy[, 2] - par[3]) / (par[4] - par[3]))
}

plt2fig <- function(xy, dev = dev.cur()) {
  olddev <- dev.cur()
  dev.set(dev)
  plt <- par("plt")
  dev.set(olddev)
  xytrans2(xy, plt)
}

fig2dev <- function(xy, dev = dev.cur()) {
  olddev <- dev.cur()
  dev.set(dev)
  fig <- par("fig")
  dev.set(olddev)
  xytrans2(xy, fig)
}

usr2plt <- function(xy, dev = dev.cur()) {
  olddev <- dev.cur()
  dev.set(dev)
  usr <- par("usr")
  dev.set(olddev)
  xytrans(xy, usr)
}

usr2dev <- function(xy, dev = dev.cur()) {
  fig2dev(plt2fig(usr2plt(xy, dev), dev), dev)
}

CalculateCircleInfo <- function(x, y, r, width, height, lbls) {
  jscode <- paste("leftImgWidth = ", width, "\n", "leftImgHeight = ", height, sep = "");
  dot.len <- length(r);
  for (i in 1:dot.len) {
    xy <- cbind(c(x[i], r[i], 0), c(y[i], 0, 0));
    xy <- usr2dev(xy, dev.cur());
    xyrc <- cbind(ceiling(xy[, 1] * width), ceiling((1 - xy[, 2]) * height));
    radius <- abs(xyrc[2, 1] - xyrc[3, 1]);
    #add code for mouseover locations, basically the annotation info
    #in this case, the name of the node
    jscode <- paste(jscode, paste("circleArray.push({xc:", xyrc[1, 1], ", yc:", xyrc[1, 2],
                                  ", r:", radius, ", lb: \"", lbls[i], "\"})", sep = ""), sep = "\n");
  }

  return(jscode);
}

PlotPathSummary <- function(mSet = NA, imgName, width = NA, height = NA) {
  x <- mSet$analSet$ora.mat[, 8]
  y <- mSet$analSet$ora.mat[, 4]
  names(x) <- names(y) <- rownames(mSet$analSet$ora.mat)
  y = -log(y)
  inx <- order(y, decreasing = T)
  x <- x[inx]
  y <- y[inx]
  sqx <- sqrt(x)
  min.x <- min(sqx, na.rm = TRUE)
  max.x <- max(sqx, na.rm = TRUE)
  if (min.x == max.x) {
    max.x = 1.5 * max.x
    min.x = 0.5 * min.x
  }
  maxR <- (max.x - min.x) / 40
  minR <- (max.x - min.x) / 160

  radi.vec <- minR + (maxR - minR) * (sqx - min.x) / (max.x -
    min.x)
  bg.vec <- heat.colors(length(y))
  path.nms <- names(metpa$path.ids)[match(names(x), metpa$path.ids)]
  bg = "white"
  imgName = imgName
  if (is.na(width)) {
    w <- 7
  }
  else if (width == 0) {
    w <- 7
  }
  else {
    w <- width
  }
  h <- w
  path.overview <- imgName

  width.px <- height.px <- w * 72
  circleInfo <- CalculateCircleInfo(x, y,
                                    radi.vec, width.px, height.px, path.nms)

  plotData <- tibble(x = x, y = y, name = path.nms, r = radi.vec, color = bg.vec) %>%
    arrange(r) %>%
    mutate(pointSize = r * 400) %>%
    filter(is.finite(r))

  print(plotData)

  if (nrow(plotData) == 0) {
    quit(status = 0)
  }

  cols <- plotData %>%
    select(c("name", "color")) %>%
    deframe()
  nameData <- plotData %>%
    filter(y >= -log(0.05))
  p <- ggplot(plotData, aes(x = x, y = y, size = r, fill = name)) +
    theme_classic(base_size = 8.8) +
    theme(axis.text.x = element_text(size = 12, vjust = 0.5),
          axis.text.y = element_text(size = 12), legend.position = 'none',
          axis.title.y = element_text(size = 12), legend.margin = margin(t = 0.3, b = 0.1, unit = 'cm'),
          legend.text = element_text(size = 6), axis.title.x = element_text(size = 12),
          plot.margin = unit(c(1, 1, 1, 1), "cm"), panel.grid.major = element_line(colour = "#BEBEBE",
                                                                                   linetype = 2, size = 0.15,),
          panel.grid.minor = element_blank(),
    ) +
    xlab("Pathway Impact") +
    ylab("-log(p)") +
    geom_point(pch = 21, colour = "black") +
    geom_text_repel(segment.size = 0.2, data = nameData, aes(x = x, y = y, label = name), size = 2, family = "Times",
                    color = "black") +
    scale_fill_manual("", values = cols) +
    scale_radius(range = c(min(plotData$pointSize), max(plotData$pointSize)), name = "")

  ggsave(limitsize = FALSE, imgName, p, width = 7, height = 7)
  return()
}

setRendAttrs = function(g, AllBorder = "transparent",
                        AllFixedsize = FALSE, AllFontsize = 16, AllShape = "rectangle",
                        fillcolor = "lightgreen", ...) {
  nn = KEGGgraph::nodes(g)
  numn = length(nn)
  color = rep(AllBorder, numn)
  names(color) = nn
  fixedsize = rep(AllFixedsize, numn)
  names(fixedsize) = nn
  if (length(fillcolor) == 1) {
    fillcolvec = rep(fillcolor, numn)
    names(fillcolvec) = nn
  } else if (!identical(names(fillcolor), as.vector(KEGGgraph::nodes(g)))) {
    stop("names on vector fillcolor must match nodes(g) exactly")
  } else {
    fillcolvec = fillcolor
  }
  shape = rep(AllShape, numn)
  names(shape) = nn
  fontsize = rep(AllFontsize, numn)
  names(fontsize) = nn;
  list(color = color, fixedsize = fixedsize, fillcolor = fillcolvec, shape = shape,
       fontsize = fontsize)
}

setRendAttrsWithName = function(g, AllBorder = "transparent",
                                AllFixedsize = FALSE, AllFontsize = 16, AllShape = "rectangle",
                                fillcolor = "lightgreen", ...) {
  nn = KEGGgraph::nodes(g)
  numn = length(nn)
  color = rep(AllBorder, numn)
  names(color) = nn
  fixedsize = rep(AllFixedsize, numn)
  names(fixedsize) = nn
  if (length(fillcolor) == 1) {
    fillcolvec = rep(fillcolor, numn)
    names(fillcolvec) = nn
  } else if (!identical(names(fillcolor), as.vector(KEGGgraph::nodes(g)))) {
    stop("names on vector fillcolor must match nodes(g) exactly")
  } else {
    fillcolvec = fillcolor
  }
  shape = rep(AllShape, numn)
  names(shape) = nn
  fontsize = rep(AllFontsize, numn)
  names(fontsize) = nn;
  labels <- names(nn)
  names(labels) <- nn
  list(color = color, fixedsize = fixedsize, fillcolor = fillcolvec, shape = shape,
       fontsize = fontsize, label = labels)
}

GetMetPANodeInfo <- function(pathName, object, tags, histvec, pvec, impvec, width, height, usr = par("usr"), dataDir) {
  nn = sapply(Rgraphviz::AgNode(object), function(x) x@name);
  ## transform user to pixel coordinates
  x.u2p = function(x) { rx = (x - usr[1]) / diff(usr[1:2]); stopifnot(all(rx >= 0 & rx <= 1)); return(rx * width) }
  y.u2p = function(y) { ry = (usr[4] - y) / diff(usr[3:4]); stopifnot(all(ry >= 0 & ry <= 1)); return(ry * height) }

  nxy = Rgraphviz::getNodeXY(object);
  nh = Rgraphviz::getNodeHeight(object) / 2;
  xl = floor(x.u2p(nxy$x - Rgraphviz::getNodeLW(object)));
  xr = ceiling(x.u2p(nxy$x + Rgraphviz::getNodeRW(object)));
  yu = floor(y.u2p(nxy$y - nh));
  yl = ceiling(y.u2p(nxy$y + nh));
  names(xl) = names(xr) = names(yu) = names(yl) = nn;
  # create the javascript code
  jscode <- paste("keggPathLnk=\'<a href=\"javascript:void(0);\" onclick=\"window.open(\\'http://www.genome.jp/kegg-bin/show_pathway?", metpa$path.ids[pathName], "\\',\\'KEGG\\');\">", pathName, "</a>\'", sep = "");
  tag.ids <- names(tags);
  kegg.ids <- names(tags);
  hmdb.ids <- KEGGID2HMDBID(kegg.ids, dataDir = dataDir)
  for (i in 1:length(tag.ids)) {
    nd <- tag.ids[i];
    x1 <- floor(100 * (xl[nd]) / width);
    x2 <- ceiling(100 * (xr[nd]) / width);
    y1 <- floor(100 * (yl[nd]) / height);
    y2 <- ceiling(100 * (yu[nd]) / height);
    #add code for mouseover locations, basically the annotation info
    #in this case, the name of the node
    jscode <- paste(jscode, paste("rectArray.push({x1:", x1, ", y1:", y1, ", x2:", x2, ", y2:", y2,
                                  ", lb: \"", tags[i], "\", kegg: \"", kegg.ids[i], "\", hmdb: \"", hmdb.ids[i],
                                  "\", icon: \"", histvec[i], "\", pvalue: \"", pvec[i], "\", impact: \"", impvec[i], "\"})", sep = ""), sep = "\n");
  }
  return(jscode);
}

PlotMetpaPath <- function(mSetObj = NA, pathName, width = NA, height = NA, dataDir, parent, fcData)
  {
  path.id <- metpa$path.ids[pathName]
  g <- metpa$graph.list[[path.id]]
  print(g)
  tooltip <- names(KEGGgraph::nodes(g))
  nm.vec <- NULL
  fillcolvec <- rep("lightgrey", length(KEGGgraph::nodes(g)))
  pvec <- histvec <- rep("NA", length(KEGGgraph::nodes(g)))
  names(tooltip) <- names(fillcolvec) <- names(histvec) <- names(pvec) <- KEGGgraph::nodes(g)
  if (mSetObj$analSet$type == "pathora") {
    if (!is.null(mSetObj$analSet$ora.hits)) {
      names <- mSetObj$analSet$ora.hits[[path.id]]
      cols <- fcData %>%
        filter(KEGG %in% names) %>%
        .$col
      fillcolvec[mSetObj$analSet$ora.hits[[path.id]]] <- cols
      if (mSetObj$dataSet$use.metabo.filter && !is.null(mSetObj$analSet$ora.filtered.mset)) {
        fillcolvec[!(names(fillcolvec) %in% mSetObj$analSet$ora.filtered.mset[[path.id]])] <- "lightgrey"
      }
    }
  }
  if (is.null(mSetObj$analSet$node.imp) || mSetObj$analSet$node.imp ==
    "rbc") {
    impvec <- metpa$rbc[[path.id]]
  }
  else {
    impvec <- metpa$dgr[[path.id]]
  }
  imgName <- paste(parent, "/", pathName, ".pdf", sep = "")
  pdf(file = imgName, width = width, height = height, bg = "white")
  par(mai = rep(0, 4))
  g.obj <- plot(g, nodeAttrs = setRendAttrs(g, fillcolor = fillcolvec))
  nodeInfo <- GetMetPANodeInfo(pathName, g.obj, tooltip,
                               histvec, pvec, impvec, width, height, dataDir = dataDir)
  dev.off()
  mSetObj$imgSet$current.metpa.graph <- g.obj
  mSetObj$analSet$nodeInfo <- nodeInfo
  if (.on.public.web) {
    .set.mSet(mSetObj)
    return(nodeInfo)
  }
  else {
    return(.set.mSet(mSetObj))
  }
}

PlotMetpaPathByName <- function(mSetObj = NA, pathName, width = NA, height = NA, dataDir, parent, fcData)
  {
  path.id <- metpa$path.ids[pathName]
  g <- metpa$graph.list[[path.id]]
  print(g)
  tooltip <- names(KEGGgraph::nodes(g))
  nm.vec <- NULL
  fillcolvec <- rep("lightgrey", length(KEGGgraph::nodes(g)))
  pvec <- histvec <- rep("NA", length(KEGGgraph::nodes(g)))
  names(tooltip) <- names(fillcolvec) <- names(histvec) <- names(pvec) <- KEGGgraph::nodes(g)
  if (mSetObj$analSet$type == "pathora") {
    if (!is.null(mSetObj$analSet$ora.hits)) {
      names <- mSetObj$analSet$ora.hits[[path.id]]
      cols <- fcData %>%
        filter(KEGG %in% names) %>%
        .$col
      fillcolvec[mSetObj$analSet$ora.hits[[path.id]]] <- cols
      if (mSetObj$dataSet$use.metabo.filter && !is.null(mSetObj$analSet$ora.filtered.mset)) {
        fillcolvec[!(names(fillcolvec) %in% mSetObj$analSet$ora.filtered.mset[[path.id]])] <- "lightgrey"
      }
    }
  }
  if (is.null(mSetObj$analSet$node.imp) || mSetObj$analSet$node.imp ==
    "rbc") {
    impvec <- metpa$rbc[[path.id]]
  }
  else {
    impvec <- metpa$dgr[[path.id]]
  }
  imgName <- paste(parent, "/", pathName, ".pdf", sep = "")
  pdf(file = imgName, width = width, height = height, bg = "white")
  par(mai = rep(0, 4))
  g.obj <- plot(g, nodeAttrs = setRendAttrsWithName(g, fillcolor = fillcolvec))
  nodeInfo <- GetMetPANodeInfo(pathName, g.obj, tooltip,
                               histvec, pvec, impvec, width, height, dataDir = dataDir)
  dev.off()
  mSetObj$imgSet$current.metpa.graph <- g.obj
  mSetObj$analSet$nodeInfo <- nodeInfo
  if (.on.public.web) {
    .set.mSet(mSetObj)
    return(nodeInfo)
  }
  else {
    return(.set.mSet(mSetObj))
  }
}

LoadMsetLib <- function(libname = "kegg_pathway", dataDir)
  {
  if (!exists("current.msetlib") || "current.msetlib$lib.name" !=
    libname) {
    .load.metaboanalyst.lib("msets", libname, dataDir = dataDir)
  }
}

SetCurrentMsetLib <- function(mSetObj = NA, lib.type, excludeNum = 0, dataDir)
  {
  mSetObj <- .get.mSet(mSetObj)
  if (lib.type != "self") {
    LoadMsetLib(lib.type, dataDir = dataDir)
  }
  if (lib.type == "self") {
    ms.list <- mSetObj$dataSet$user.mset
    current.msetlib <- data.frame(name = character(), member = character(),
                                  reference = character(), stringsAsFactors = FALSE)
  }
  else {
    ms.list <- strsplit(current.msetlib[, 3], "; ")
    names(ms.list) <- current.msetlib[, 2]
  }

  if (excludeNum > 0) {
    cmpd.count <- lapply(ms.list, length)
    sel.inx <- cmpd.count >= excludeNum
    ms.list <- ms.list[sel.inx]
    current.msetlib <- current.msetlib[sel.inx,]
  }
  mSetObj$dataSet$uniq.count <- length(unique(unlist(ms.list,
                                                     use.names = FALSE)))
  current.msetlib$member <- ms.list
  current.msetlib <<- current.msetlib
  if (.on.public.web) {
    .set.mSet(mSetObj)
  }
  return(.set.mSet(mSetObj))
}

CalculateHyperScore <- function(mSetObj = NA, dataDir, fileName)
  {
  mSetObj <- .get.mSet(mSetObj)
  nm.map <- GetFinalNameMap(mSetObj)
  
  valid.inx <- !(is.na(nm.map$hmdb) | duplicated(nm.map$hmdb))
  ora.vec <- nm.map$hmdb[valid.inx]
  q.size <- length(ora.vec)
  if (is.na(ora.vec) || q.size == 0) {
    AddErrMsg("No valid HMDB compound names found!")
    return(0)
  }
  current.mset <- current.msetlib$member
  if (mSetObj$dataSet$use.metabo.filter && !is.null(mSetObj$dataSet$metabo.filter.hmdb)) {
    current.mset <- lapply(current.mset, function(x) {
      x[x %in% mSetObj$dataSet$metabo.filter.hmdb]
    })
    mSetObj$dataSet$filtered.mset <- current.mset
  }
  uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)))
  set.size <- length(current.mset)
  if (set.size == 1) {
    AddErrMsg("Cannot perform enrichment analysis on a single metabolite set!")
    return(0)
  }
  hits <- lapply(current.mset, function(x) {
    x[x %in% ora.vec]
  })
  
  hit.num <- unlist(lapply(hits, function(x) length(x)), use.names = FALSE)
  if (sum(hit.num > 0) == 0) {
    AddErrMsg("No match was found to the selected metabolite set library!")
    return(0)
  }
  set.num <- unlist(lapply(current.mset, length), use.names = FALSE)
  res.mat <- matrix(NA, nrow = set.size, ncol = 6)
  rownames(res.mat) <- names(current.mset)
  colnames(res.mat) <- c("Total", "Expected", "Hits", "Raw P",
                         "Holm P", "FDR")
  for (i in 1:set.size) {
    res.mat[i, 1] <- set.num[i]
    res.mat[i, 2] <- q.size * (set.num[i] / uniq.count)
    res.mat[i, 3] <- hit.num[i]
    res.mat[i, 4] <- phyper(hit.num[i] - 1, set.num[i], uniq.count -
      set.num[i], q.size, lower.tail = F)
  }

  enrichNames <- c(1:set.size) %>%
    map(~hits[.x]) %>%
    map_chr(function(x) {
      xVec <- x %>%
        unlist()
      ifelse(length(xVec) > 0, str_c(xVec, collapse = " || "), "")
    })

  res.mat[, 5] <- p.adjust(res.mat[, 4], "holm")
  res.mat[, 6] <- p.adjust(res.mat[, 4], "fdr")

  save.mat <- res.mat %>%
    as.data.frame() %>%
    rownames_to_column("name") %>%
    add_column(Enriched_Compounds = enrichNames) %>%
    filter(Hits > 0) %>%
    arrange(`Raw P`) %>%
    mutate_at(vars(-c("name", "Enriched_Compounds")), ~signif(.x, 3))

  mSetObj$analSet$ora.mat <- save.mat %>%
    select(-c("Enriched_Compounds")) %>%
    column_to_rownames("name") %>%
    as.matrix()
  mSetObj$analSet$ora.hits <- hits

  save.mat <- save.mat %>%
    column_to_rownames("name")

  write.csv(save.mat, file = fileName)
  write.table(save.mat, "MSEA_Result.txt", quote = FALSE, sep = "\t")
  if (.on.public.web) {
    .set.mSet(mSetObj)
    return(1)
  }
  return(.set.mSet(mSetObj))
}

PlotORA <- function(mSetObj = NA, imgName, imgOpt, format = "png",
                    width = NA, parent)
  {
  mSetObj <- .get.mSet(mSetObj)
  folds <- mSetObj$analSet$ora.mat[, 3] / mSetObj$analSet$ora.mat[,
    2]
  names(folds) <- GetShortNames(rownames(mSetObj$analSet$ora.mat))
  pvals <- mSetObj$analSet$ora.mat[, 4]
  imgName = paste(imgName, ".", format, sep = "")
  if (is.na(width)) {
    w <- 9
  }
  else if (width == 0) {
    w <- 7
  }
  else {
    w <- width
  }
  h <- w
  mSetObj$imgSet$ora <- imgName
  mSetObj$imgSet$current.img <- imgName
  Cairo::Cairo(file = imgName, unit = "in", width = w,
               height = h, type = format, bg = "white")
  PlotMSEA.Overview(folds, pvals)
  dev.off()
  if (!.on.public.web) {
    g <- mSetObj$analSet$enrich.net <- PlotEnrichNet.Overview(folds,
                                                              pvals, parent = parent)
    if (is.null(g)) {
      return(.set.mSet(mSetObj))
    }
    pdf(str_c(parent, "/Network.pdf"), 10, 10)
    plot(g, layout = layout.fruchterman.reingold)
    dev.off()
  }
  return(.set.mSet(mSetObj))
}

GetShortNames <- function(nm.vec, max.len = 45) {
  new.nms <- vector(mode = "character", length = length(nm.vec));
  for (i in 1:length(nm.vec)) {
    nm <- nm.vec[i];
    if (nchar(nm) <= max.len) {
      new.nms[i] <- nm;
    }else {
      wrds <- strsplit(nm, "[[:space:]]+")[[1]];
      new.nm <- "";
      if (length(wrds) > 1) {
        for (m in 1:length(wrds)) {
          wrd <- wrds[m];
          if (nchar(new.nm) + 4 + nchar(wrd) <= max.len) {
            new.nm <- paste(new.nm, wrd);
          }else {
            new.nms[i] <- paste(new.nm, "...", sep = "");
            break;
          }
        }
      }else {
        new.nms[i] <- paste(substr(nm, 0, 21), "...", sep = "");
      }
    }
  }
  return(new.nms);
}

PlotMSEA.Overview <- function(folds, pvals)
  {
  title <- "Metabolite Sets Enrichment Overview"
  if (length(folds) > 50) {
    folds <- folds[1:50]
    pvals <- pvals[1:50]
    title <- "Enrichment Overview (top 50)"
  }
  op <- par(mar = c(5, 20, 4, 6), oma = c(0, 0, 0, 4))
  ht.col <- rev(heat.colors(length(folds)))
  barplot(rev(folds), horiz = T, col = ht.col, xlab = "Fold Enrichment",
          las = 1, cex.name = 0.75, space = c(0.5, 0.5), main = title)
  minP <- min(pvals)
  maxP <- max(pvals)
  medP <- (minP + maxP) / 2
  axs.args <- list(at = c(minP, medP, maxP), labels = format(c(maxP,
                                                               medP, minP), scientific = T, digit = 1), tick = F)
  image.plot(legend.only = TRUE, zlim = c(minP, maxP), col = ht.col,
             axis.args = axs.args, legend.shrink = 0.4, legend.lab = "P value")
  par(op)
}

image.plot <- function(..., add = FALSE, nlevel = 64,
                       horizontal = FALSE, legend.shrink = 0.9, legend.width = 1.2,
                       legend.mar = ifelse(horizontal, 3.1, 5.1), legend.lab = NULL,
                       graphics.reset = FALSE, bigplot = NULL, smallplot = NULL,
                       legend.only = FALSE, col = tim.colors(nlevel), lab.breaks = NULL,
                       axis.args = NULL, legend.args = NULL, midpoint = FALSE) {
  old.par <- par(no.readonly = TRUE)
  #  figure out zlim from passed arguments
  info <- image.plot.info(...)
  if (add) {
    big.plot <- old.par$plt
  }
  if (legend.only) {
    graphics.reset <- TRUE
  }
  if (is.null(legend.mar)) {
    legend.mar <- ifelse(horizontal, 3.1, 5.1)
  }
  #
  # figure out how to divide up the plotting real estate.
  #
  temp <- image.plot.plt(add = add, legend.shrink = legend.shrink,
                         legend.width = legend.width, legend.mar = legend.mar,
                         horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
  #
  # bigplot are plotting region coordinates for image
  # smallplot are plotting coordinates for legend
  smallplot <- temp$smallplot
  bigplot <- temp$bigplot
  #
  # draw the image in bigplot, just call the R base function
  # or poly.image for polygonal cells note logical switch
  # for poly.grid parsed out of call from image.plot.info
  if (!legend.only) {
    if (!add) {
      par(plt = bigplot)
    }
    if (!info$poly.grid) {
      image(..., add = add, col = col)
    }
    else {
      poly.image(..., add = add, col = col, midpoint = midpoint)
    }
    big.par <- par(no.readonly = TRUE)
  }
  ##
  ## check dimensions of smallplot
  if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
    par(old.par)
    stop("plot region too small to add legend\n")
  }
  # Following code draws the legend using the image function
  # and a one column image.
  # calculate locations for colors on legend strip
  ix <- 1
  minz <- info$zlim[1]
  maxz <- info$zlim[2]
  binwidth <- (maxz - minz) / nlevel
  midpoints <- seq(minz + binwidth / 2, maxz - binwidth / 2, by = binwidth)
  iy <- midpoints
  iz <- matrix(iy, nrow = 1, ncol = length(iy))
  # extract the breaks from the ... arguments
  # note the breaks delineate intervals of common color
  breaks <- list(...)$breaks
  # draw either horizontal or vertical legends.
  # using either suggested breaks or not -- a total of four cases.
  #
  # next par call sets up a new plotting region just for the legend strip
  # at the smallplot coordinates
  par(new = TRUE, pty = "m", plt = smallplot, err = -1)
  # create the argument list to draw the axis
  #  this avoids 4 separate calls to axis and allows passing extra
  # arguments.
  # then add axis with specified lab.breaks at specified breaks
  if (!is.null(breaks) & !is.null(lab.breaks)) {
    # axis with labels at break points
    axis.args <- c(list(side = ifelse(horizontal, 1, 4),
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2),
                        at = breaks, labels = lab.breaks), axis.args)
  }
  else {
    # If lab.breaks is not specified, with or without breaks, pretty
    # tick mark locations and labels are computed internally,
    # or as specified in axis.args at the function call
    axis.args <- c(list(side = ifelse(horizontal, 1, 4),
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)),
                   axis.args)
  }
  #
  # draw color scales the four cases are horizontal/vertical breaks/no breaks
  # add a label if this is passed.
  if (!horizontal) {
    if (is.null(breaks)) {
      image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
            ylab = "", col = col)
    }
    else {
      image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
            ylab = "", col = col, breaks = breaks)
    }
  }
  else {
    if (is.null(breaks)) {
      image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",
            ylab = "", col = col)
    }
    else {
      image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",
            ylab = "", col = col, breaks = breaks)
    }
  }
  #
  # now add the axis to the legend strip.
  # notice how all the information is in the list axis.args
  #
  do.call("axis", axis.args)
  # add a box around legend strip
  box()
  #
  # add a label to the axis if information has been  supplied
  # using the mtext function. The arguments to mtext are
  # passed as a list like the drill for axis (see above)
  #
  if (!is.null(legend.lab)) {
    legend.args <- list(text = legend.lab, side = ifelse(horizontal,
                                                         1, 3), line = 1)
  }
  #
  # add the label using mtext function
  if (!is.null(legend.args)) {
    do.call(mtext, legend.args)
  }
  #
  #
  # clean up graphics device settings
  # reset to larger plot region with right user coordinates.
  mfg.save <- par()$mfg
  if (graphics.reset | add) {
    par(old.par)
    par(mfg = mfg.save, new = FALSE)
    invisible()
  }
  else {
    par(big.par)
    par(plt = big.par$plt, xpd = FALSE)
    par(mfg = mfg.save, new = FALSE)
    invisible()
  }
}

"image.plot.info" <- function(...) {
  temp <- list(...)
  #
  xlim <- NA
  ylim <- NA
  zlim <- NA
  poly.grid <- FALSE
  #
  # go through various cases of what these can be
  #
  ##### x,y,z list is first argument
  if (is.list(temp[[1]])) {
    xlim <- range(temp[[1]]$x, na.rm = TRUE)
    ylim <- range(temp[[1]]$y, na.rm = TRUE)
    zlim <- range(temp[[1]]$z, na.rm = TRUE)
    if (is.matrix(temp[[1]]$x) &
      is.matrix(temp[[1]]$y) &
      is.matrix(temp[[1]]$z)) {
      poly.grid <- TRUE
    }
  }
  ##### check for polygrid first three arguments should be matrices
  #####
  if (length(temp) >= 3) {
    if (is.matrix(temp[[1]]) &
      is.matrix(temp[[2]]) &
      is.matrix(temp[[3]])) {
      poly.grid <- TRUE
    }
  }
  #####  z is passed without an  x and y  (and not a poly.grid!)
  #####
  if (is.matrix(temp[[1]]) & !poly.grid) {
    xlim <- c(0, 1)
    ylim <- c(0, 1)
    zlim <- range(temp[[1]], na.rm = TRUE)
  }
  ##### if x,y,z have all been passed find their ranges.
  ##### holds if poly.grid or not
  #####
  if (length(temp) >= 3) {
    if (is.matrix(temp[[3]])) {
      xlim <- range(temp[[1]], na.rm = TRUE)
      ylim <- range(temp[[2]], na.rm = TRUE)
      zlim <- range(temp[[3]], na.rm = TRUE)
    }
  }
  #### parse x,y,z if they are  named arguments
  # determine if  this is polygon grid (x and y are matrices)
  if (is.matrix(temp$x) &
    is.matrix(temp$y) &
    is.matrix(temp$z)) {
    poly.grid <- TRUE
  }
  xthere <- match("x", names(temp))
  ythere <- match("y", names(temp))
  zthere <- match("z", names(temp))
  if (!is.na(zthere))
    zlim <- range(temp$z, na.rm = TRUE)
  if (!is.na(xthere))
    xlim <- range(temp$x, na.rm = TRUE)
  if (!is.na(ythere))
    ylim <- range(temp$y, na.rm = TRUE)
  # overwrite zlims with passed values
  if (!is.null(temp$zlim))
    zlim <- temp$zlim
  if (!is.null(temp$xlim))
    xlim <- temp$xlim
  if (!is.null(temp$ylim))
    ylim <- temp$ylim
  list(xlim = xlim, ylim = ylim, zlim = zlim, poly.grid = poly.grid)
}

image.plot.plt <- function(x, add = FALSE, legend.shrink = 0.9,
                           legend.width = 1, horizontal = FALSE, legend.mar = NULL,
                           bigplot = NULL, smallplot = NULL, ...) {
  old.par <- par(no.readonly = TRUE)
  if (is.null(smallplot))
    stick <- TRUE
  else stick <- FALSE
  if (is.null(legend.mar)) {
    legend.mar <- ifelse(horizontal, 3.1, 5.1)
  }
  # compute how big a text character is
  char.size <- ifelse(horizontal, par()$cin[2] / par()$din[2],
                      par()$cin[1] / par()$din[1])
  # This is how much space to work with based on setting the margins in the
  # high level par command to leave between strip and big plot
  offset <- char.size * ifelse(horizontal, par()$mar[1], par()$mar[4])
  # this is the width of the legned strip itself.
  legend.width <- char.size * legend.width
  # this is room for legend axis labels
  legend.mar <- legend.mar * char.size
  # smallplot is the plotting region for the legend.
  if (is.null(smallplot)) {
    smallplot <- old.par$plt
    if (horizontal) {
      smallplot[3] <- legend.mar
      smallplot[4] <- legend.width + smallplot[3]
      pr <- (smallplot[2] - smallplot[1]) * ((1 - legend.shrink) / 2)
      smallplot[1] <- smallplot[1] + pr
      smallplot[2] <- smallplot[2] - pr
    }
    else {
      smallplot[2] <- 1 - legend.mar
      smallplot[1] <- smallplot[2] - legend.width
      pr <- (smallplot[4] - smallplot[3]) * ((1 - legend.shrink) / 2)
      smallplot[4] <- smallplot[4] - pr
      smallplot[3] <- smallplot[3] + pr
    }
  }
  if (is.null(bigplot)) {
    bigplot <- old.par$plt
    if (!horizontal) {
      bigplot[2] <- min(bigplot[2], smallplot[1] - offset)
    }
    else {
      bottom.space <- old.par$mar[1] * char.size
      bigplot[3] <- smallplot[4] + offset
    }
  }
  if (stick & (!horizontal)) {
    dp <- smallplot[2] - smallplot[1]
    smallplot[1] <- min(bigplot[2] + offset, smallplot[1])
    smallplot[2] <- smallplot[1] + dp
  }
  return(list(smallplot = smallplot, bigplot = bigplot))
}

PlotEnrichNet.Overview <- function(folds, pvals, layoutOpt = layout.fruchterman.reingold, parent)
  {
  title <- "Enrichment Network Overview"
  if (length(folds) > 50) {
    folds <- folds[1:50]
    pvals <- pvals[1:50]
    title <- "Enrichment Overview (top 50)"
  }
  if (.on.public.web) {
    load_igraph()
    load_reshape()
  }
  pvalue <- pvals
  id <- names(pvalue)
  geneSets <- current.msetlib$member
  n <- length(pvalue)
  w <- matrix(NA, nrow = n, ncol = n)
  colnames(w) <- rownames(w) <- id
  for (i in 1:n) {
    for (j in i:n) {
      w[i, j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
    }
  }
  wd <- melt(w)
  wd <- wd[wd[, 1] != wd[, 2],]
  wd <- wd[!is.na(wd[, 3]),]
  g <- graph.data.frame(wd[, -3], directed = F)
  E(g)$width <- sqrt(wd[, 3] * 20)
  g <- delete.edges(g, E(g)[wd[, 3] < 0.2])
  idx <- unlist(sapply(V(g)$name, function(x) which(x == id)))
  cols <- color_scale("red", "#E5C494")
  V(g)$color <- cols[sapply(pvalue, getIdx, min(pvalue), max(pvalue))]
  cnt <- folds
  names(cnt) <- id
  cnt2 <- cnt[V(g)$name] + 1
  V(g)$size <- cnt2 / sum(cnt2) * 100
  pos.xy <- layout.fruchterman.reingold(g, niter = 500)
  nodes <- vector(mode = "list")
  node.nms <- V(g)$name
  node.sizes <- V(g)$size
  node.cols <- V(g)$color
  if (length(node.sizes) == 0) {
    return(NULL)
  }

  for (i in 1:length(node.sizes)) {
    nodes[[i]] <- list(id = node.nms[i], label = node.nms[i],
                       size = node.sizes[i], color = node.cols[i], x = pos.xy[i,
        1], y = pos.xy[i, 2])
  }

  edge.mat <- get.edgelist(g)
  edge.mat <- cbind(id = 1:nrow(edge.mat), source = edge.mat[,
    1], target = edge.mat[, 2])
  netData <- list(nodes = nodes, edges = edge.mat)
  sink(str_c(parent, "/MSEA_Network.json"))
  cat(RJSONIO::toJSON(netData))
  sink()
  return(g)
}

overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y)) / length(unique(c(x, y)))
}

color_scale <- function(c1 = "grey", c2 = "red") {
  pal <- colorRampPalette(c(c1, c2))
  colors <- pal(100)
  return(colors)
}

getIdx <- function(v, MIN, MAX) {
  if (MIN == MAX) {
    return(100)
  }
  intervals <- seq(MIN, MAX, length.out = 100)
  max(which(intervals <= v))
}

SetSMPDB.PathLib <- function(mSetObj = NA, smpdb.rda, dataDir)
  {
  mSetObj <- .get.mSet(mSetObj)
  mSetObj$msgSet$lib.msg <- paste("Your selected pathway library code is \\textbf{",
                                  smpdb.rda, "}(KEGG organisms abbreviation).")
  smpdblib <- .load.metaboanalyst.lib("smpdb", smpdb.rda, dataDir = dataDir)
  mSetObj$pathwaylibtype <- "SMPDB"
  if (.on.public.web) {
    .set.mSet(mSetObj)
    return(1)
  }
  return(.set.mSet(mSetObj))
}

SetOrganism <- function(mSetObj = NA, org)
  {
  mSetObj <- .get.mSet(mSetObj)
  pathinteg.org <<- data.org <<- org
  return(.set.mSet(mSetObj))
}

smpdbpw.count <- 0

GeneratePathwayJSON <- function(mSetObj = NA, pathway.nm, dataDir, data.org, parent) {
  mSetObj <- .get.mSet(mSetObj);
  smpdb.path <- paste(dataDir, "/smpdb/", data.org, ".rda", sep = "");
  load(smpdb.path)

  jsons.path <- paste(dataDir, "/smpdb/jsons/", data.org, ".rds", sep = "");
  smpdb.jsons <- readRDS(jsons.path) # no need to be global!

  if (pathway.nm == "top") {
    if (mSetObj$analSet$type == "pathora") {
      pathway.id <- rownames(mSetObj$analSet$ora.mat)[1]
    } else {
      pathway.id <- rownames(mSetObj$analSet$qea.mat)[1]
    }
    pathway.nm <- names(metpa$path.ids)[which(metpa$path.ids == pathway.id)]
  } else {
    pathway.id <- metpa$path.ids[which(names(metpa$path.ids) == pathway.nm)]
  }
  # Get matched metabolites
  if (mSetObj$analSet$type == "pathora") {
    metab.matches <- paste(mSetObj$analSet$ora.hits[[pathway.id]], collapse = ",");
  } else {
    metab.matches <- paste(mSetObj$analSet$qea.hits[[pathway.id]], collapse = ",");
  }

  title <- paste(pathway.id, ";", pathway.nm, sep = "");
  # store json file
  smpdbpw.nm <- paste(parent, "/", pathway.nm, ".json", sep = "");
  smpdbpw.count <<- smpdbpw.count + 1;
  g <- smpdb.jsons[[pathway.id]]
  sink(smpdbpw.nm)
  cat(jsonlite::toJSON(g, pretty = TRUE));
  sink()

  smpdbpw.nm <- paste0(smpdbpw.nm, ";", metab.matches, ";", title)
  return(smpdbpw.nm)
}









