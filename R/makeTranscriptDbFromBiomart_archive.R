### =========================================================================
### makeTranscriptDbFromBiomart_archive()
### -------------------------------------------------------------------------
###
### For people who want to tap BioMart.
### Typical use:
###   txdb <- makeTxDbFromBiomart(biomart="ensembl",
###                                       dataset="hsapiens_gene_ensembl")
### Speed:
###   - for biomart="ensembl" and dataset="hsapiens_gene_ensembl":
###       (1) download takes about 8 min.
###       (2) db creation takes about 60-65 sec.
###


## helper to extract the Genus and species name from the dataset string.
.extractSpeciesFromDatasetDesc <- function(description){
  vals <- unlist(strsplit(description, " "))
  paste(vals[[1]], vals[[2]])
}


.getBiomartDbVersion <- function(mart, host="mar2009.archive.ensembl.org", 
    path="/biomart/martservice", biomart, archive=FALSE)
{
    marts <- listMarts(mart=mart,host=host,path=path,archive=FALSE)
    mart_rowidx <- which(as.character(marts$biomart) == biomart)
    ## This should never happen.
    if (length(mart_rowidx) != 1L)
        stop("found 0 or more than 1 \"", biomart, "\" BioMart database")
    as.character(marts$version)[mart_rowidx]
}

.extractEnsemblReleaseFromDbVersion <- function(db_version)
    # strsplit(db_version," ")[[1]][2]
     substr(db_version,gregexpr('[0-9]',db_version)[[1]][1],gregexpr('[0-9]',db_version)[[1]][2])
   # sub("^ENSEMBL GENES ([^[:space:]]+) \\(SANGER UK\\)", "\\1", db_version)    # for current version

### Groups of BioMart attributes:
###   - A1, A2 and G are required attributes;
###   - B, C and D are optional attributes: C is required for inferring the
###     CDS (they cannot be inferred from D). Therefore, if C is missing,
###     the TxDb object can still be made but won't have any CDS (no
###     row in the cds table). D is only used for sanity check.
.A1_ATTRIBS <- c("ensembl_transcript_id",
                 "chromosome_name",
                 "strand",
                 "transcript_start",
                 "transcript_end")

.A2_ATTRIBS <- c("ensembl_transcript_id",
                 "strand",
                 "rank",
                 "exon_chrom_start",
                 "exon_chrom_end")

.B_ATTRIB <- "ensembl_exon_id"

.C_ATTRIBS <- c("5_utr_start",
                "5_utr_end",
                "3_utr_start",
                "3_utr_end")

.D_ATTRIBS <- c("cds_start",
                "cds_end",
                "cds_length")

.G_ATTRIB <- "ensembl_gene_id"

### 'attribs' can be either a Mart object or a 2-col data frame as returned by
### 'listAttributes()'.
.getDatasetAttrGroups <- function(attribs)
{
    if (is(attribs, "Mart"))
        attribs <- listAttributes(attribs)
    else if (!is.data.frame(attribs) ||
             !identical(colnames(attribs), c("name", "description")))
        stop("invalid 'attribs' object")
    attrgroups <- "none"
    ## Group A: Required attributes.
    attrA <- unique(c(.A1_ATTRIBS, .A2_ATTRIBS))
    if (all(attrA %in% attribs$name))
        attrgroups <- "A"
    ## Groups B, C and D are optional attributes.
    ## C is required for inferring the CDS (they cannot be inferred from D).
    ## Therefore, if C is missing, the TxDb object can still be made
    ## but won't have any CDS (no row in the cds table).
    if (.B_ATTRIB %in% attribs$name)
        attrgroups <- paste(attrgroups, "B", sep="")
    if (all(.C_ATTRIBS %in% attribs$name))
        attrgroups <- paste(attrgroups, "C", sep="")
    if (all(.D_ATTRIBS %in% attribs$name))
        attrgroups <- paste(attrgroups, "D", sep="")
    ## Group G: Required attribute.
    if (.G_ATTRIB %in% attribs$name)
        attrgroups <- paste(attrgroups, "G", sep="")
    attrgroups
}

### 'attrlist' can be a list (as returned by getMartAttribList()), a Mart
### object, or the name of a Mart service (single string).
### Typical use:
###   ensembl_attrgroups <-
###       GenomicFeatures:::.getAllDatasetAttrGroups("ensembl")
.getAllDatasetAttrGroups <- function(attrlist)
{
    if (!is.list(attrlist))
        attrlist <- getMartAttribList(attrlist)
    sapply(attrlist, .getDatasetAttrGroups)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'transcripts' data frame.
###

.makeBiomartTranscripts <- function(filters, values, mart, transcript_ids)
{
    message("Download and preprocess the 'transcripts' data frame ... ",
            appendLF=FALSE)
    bm_table <- getBM(.A1_ATTRIBS, filters=filters, values=values, mart=mart)
    colnames(bm_table) <- .A1_ATTRIBS
    if (!is.null(transcript_ids)) {
        idx <- !(transcript_ids %in% bm_table$ensembl_transcript_id)
        if (any(idx)) {
            bad_ids <- transcript_ids[idx]
            stop("invalid transcript ids: ",
                 paste(bad_ids, collapse=", "), sep="")
        }
    }
    ## Those are the strictly required fields.
    transcripts0 <- data.frame(
        tx_id=integer(0),
        tx_chrom=character(0),
        tx_strand=character(0),
        tx_start=integer(0),
        tx_end=integer(0)
    )
    if (nrow(bm_table) == 0L) {
        message("OK")
        return(transcripts0)
    }
    transcripts_tx_id <- seq_len(nrow(bm_table))
    transcripts_tx_name <- bm_table$ensembl_transcript_id
    ## if (any(duplicated(transcripts_tx_name)))
    ##     stop("the 'ensembl_transcript_id' attribute contains duplicated values")
    if (any(duplicated(bm_table)))
        stop("The 'transcripts' data frame from biomart contains duplicated rows.")
    transcripts <- data.frame(
        tx_id=transcripts_tx_id,
        tx_name=transcripts_tx_name,
        tx_chrom=as.character(bm_table$chromosome_name),
        tx_strand=ifelse(bm_table$strand == 1, "+", "-"),
        tx_start=bm_table$transcript_start,
        tx_end=bm_table$transcript_end
    )
    message("OK")
    transcripts
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'chrominfo' data frame.
###

### Returns NULL if it fails to fetch the chromosome lengths from the
### remote resource.
.makeBiomartChrominfo <- function(mart, extra_seqnames=NULL,
                            circ_seqs=character(0), 
                            host="mar2009.archive.ensembl.org", 
                            path="/biomart/martservice", 
                            archive=FALSE)
{
    biomart <- biomaRt:::martBM(mart)
    dataset <- biomaRt:::martDataset(mart)
    if (biomart == "ENSEMBL_MART_ENSEMBL") {
        message("Download and preprocess the 'chrominfo' data frame ... ",
                appendLF=FALSE)
        db_version <- .getBiomartDbVersion(mart, host=host, path=path, 
                    biomart, archive=FALSE)
        ensembl_release <- .extractEnsemblReleaseFromDbVersion(db_version)
        chromlengths <- try(fetchChromLengthsFromEnsembl(dataset,
                                release=ensembl_release,
                                extra_seqnames=extra_seqnames),
                            silent=TRUE)
        if (is(chromlengths, "try-error")) {
            message("FAILED! (=> skipped)")
            return(NULL)
        }
        chrominfo <- data.frame(
            chrom=chromlengths$name,
            length=chromlengths$length,
            is_circular=matchCircularity(chromlengths$name, circ_seqs)
        )
        message("OK")
        return(chrominfo)
    }
    NULL
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Allow users to discover 'chrominfo' data frame.
###

getChromInfoFromBiomart <- function(biomart="ENSEMBL_MART_ENSEMBL",
                                    dataset="hsapiens_gene_ensembl",host="mar2009.archive.ensembl.org",path="/biomart/martservice",archive=FALSE)
{
    mart <- .parseBMMartParams(biomart=biomart,
                               dataset=dataset,host=host,path=path,archive=FALSE)
    filters <- .parseBMFiltersParams(transcript_ids=NULL)
    values <- .parseBMValuesParams(transcript_ids=NULL)
    transcripts <- .makeBiomartTranscripts(filters, values, mart,
                                           transcript_ids=NULL)
    chrominfo <- .makeBiomartChrominfo(mart,
                                       extra_seqnames=transcripts$tx_chrom)
    chrominfo[,1:2]
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'splicings' data frame.
###

.normUtrCoords <- function(coords)
{
    if (is.numeric(coords))
        return(coords)
    if (is.logical(coords) && all(is.na(coords)))
        return(as.integer(coords))
    stop("BioMart data anomaly: utr coordinates don't ",
         "have a numeric type")
}

.extractCdsRangesFromBiomartTable <- function(bm_table)
{
    if (nrow(bm_table) == 0L)
        return(IRanges())
    strand <- bm_table[["strand"]]
    cds_start <- exon_start <- bm_table[["exon_chrom_start"]]
    cds_end <- exon_end <- bm_table[["exon_chrom_end"]]
    utr5_start <- .normUtrCoords(bm_table[["5_utr_start"]])
    utr5_end <- .normUtrCoords(bm_table[["5_utr_end"]])
    utr3_start <- .normUtrCoords(bm_table[["3_utr_start"]])
    utr3_end <- .normUtrCoords(bm_table[["3_utr_end"]])

    if (!all(strand %in% c(1, -1)))
        stop("BioMart data anomaly: \"strand\" attribute should be 1 or -1")
    if (!is.numeric(exon_start) || !is.numeric(exon_end))
        stop("BioMart data anomaly: exon coordinates don't ",
             "have a numeric type")
    no_utr5 <- is.na(utr5_start)
    if (!identical(no_utr5, is.na(utr5_end)))
        stop("BioMart data anomaly: NAs in \"5_utr_start\" attribute ",
             "don't match NAs in \"5_utr_end\" attribute")
    if (!all(utr5_start <= utr5_end, na.rm=TRUE))
        stop("BioMart data anomaly: some 5' UTR have a start > end")
    if (!all(utr5_start >= exon_start, na.rm=TRUE)
     || !all(utr5_end <= exon_end, na.rm=TRUE))
        stop("BioMart data anomaly: some 5' UTR are not within the exon limits")
    no_utr3 <- is.na(utr3_start)
    if (!identical(no_utr3, is.na(utr3_end)))
        stop("BioMart data anomaly: NAs in \"3_utr_start\" attribute ",
             "don't match NAs in \"3_utr_end\" attribute")
    if (!all(utr3_start <= utr3_end, na.rm=TRUE))
        stop("BioMart data anomaly: some 3' UTR have a start > end")
    if (!all(utr3_start >= exon_start, na.rm=TRUE)
     || !all(utr3_end <= exon_end, na.rm=TRUE))
        stop("BioMart data anomaly: some 3' UTR are not within the exon limits")

    idx <- strand == 1 & !no_utr5
    if (!all(utr5_start[idx] == exon_start[idx]))
        stop("BioMart data anomaly: some 5' UTR on the plus strand ",
             "don't start where the exon starts")
    cds_start[idx] <- utr5_end[idx] + 1L
    idx <- strand == 1 & !no_utr3
    if (!all(utr3_end[idx] == exon_end[idx]))
        stop("BioMart data anomaly: some 3' UTR on the plus strand ",
             "don't end where the exon ends")
    cds_end[idx] <- utr3_start[idx] - 1L
    idx <- strand == -1 & !no_utr3
    if (!all(utr3_start[idx] == exon_start[idx]))
        stop("BioMart data anomaly: some 3' UTR on the minus strand ",
             "don't start where the exon starts")
    cds_start[idx] <- utr3_end[idx] + 1L
    idx <- strand == -1 & !no_utr5
    if (!all(utr5_end[idx] == exon_end[idx]))
        stop("BioMart data anomaly: some 5' UTR on the minus strand ",
             "don't end where the exon ends")
    cds_end[idx] <- utr5_start[idx] - 1L
    ans <- IRanges(start=cds_start, end=cds_end)
    if (length(ans) != 0L) {
        cds_cumlength <-
            sapply(split(width(ans), bm_table$ensembl_transcript_id), sum)
        #if (!all(cds_cumlength[as.vector(bm_table$ensembl_transcript_id)]
        #         == bm_table$cds_length, na.rm=TRUE))
        #    stop("BioMart data anomaly: for some transcripts, the cds ",
        #         "cumulative length inferred from the exon and UTR info ",
        #         "doesn't match the \"cds_length\" attribute from BioMart")
        #if (!all(cds_cumlength %% 3L == 0L))
        #    warning("BioMart data anomaly: for some transcripts, the cds ",
        #            "cumulative length (\"cds_length\" attribute) is not ",
        #            "a multiple of 3")
    }
    ans
}

.makeCdsDataFrameFromRanges <- function(cds_ranges)
{
    nocds_idx <- width(cds_ranges) == 0L
    cds_start <- start(cds_ranges)
    cds_start[nocds_idx] <- NA_integer_
    cds_end <- end(cds_ranges)
    cds_end[nocds_idx] <- NA_integer_
    data.frame(cds_start=cds_start, cds_end=cds_end)
}

### Ironically the cds_start and cds_end attributes that we get from
### BioMart are pretty useless because they are relative to the coding
### mRNA. However, the utr coordinates are relative to the chromosome so
### we use them to infer the cds coordinates. We also retrieve the
### cds_length attribute as a sanity check.
.makeBiomartSplicings <- function(filters, values, mart, transcripts_tx_name)
{
    ## Those are the strictly required fields.
    splicings0 <- data.frame(
        tx_id=integer(0),
        exon_rank=integer(0),
        exon_start=integer(0),
        exon_end=integer(0)
    )
    if (length(transcripts_tx_name) == 0L)
        return(splicings0)
    message("Download and preprocess the 'splicings' data frame ... ",
            appendLF=FALSE)
    allattribs <- listAttributes(mart)$name
    attributes <- .A2_ATTRIBS
    if (.B_ATTRIB %in% allattribs)
        attributes <- c(attributes, .B_ATTRIB)
    if (all(.C_ATTRIBS %in% allattribs))
        attributes <- c(attributes, .C_ATTRIBS)
    if ("cds_length" %in% allattribs)
        attributes <- c(attributes, "cds_length")
    bm_table <- getBM(attributes, filters=filters, values=values, mart=mart)
    colnames(bm_table) <- attributes
    splicings_tx_id <- as.integer(factor(bm_table$ensembl_transcript_id,
                                         levels=transcripts_tx_name))
    splicings <- data.frame(
        tx_id=splicings_tx_id,
        exon_rank=bm_table$rank,
        exon_name=bm_table$ensembl_exon_id,
        exon_start=bm_table$exon_chrom_start,
        exon_end=bm_table$exon_chrom_end
    )
    if (all(.C_ATTRIBS %in% allattribs) && ("cds_length" %in% allattribs)) {
        cds_ranges <- .extractCdsRangesFromBiomartTable(bm_table)
        splicings <- cbind(splicings, .makeCdsDataFrameFromRanges(cds_ranges))
    }
    message("OK")
    splicings
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'genes' data frame.
###

.makeBiomartGenes <- function(filters, values, mart, transcripts_tx_name)
{
    message("Download and preprocess the 'genes' data frame ... ",
            appendLF=FALSE)
    attributes <- c(.G_ATTRIB, "ensembl_transcript_id")
    bm_table <- getBM(attributes, filters=filters, values=values, mart=mart)
    colnames(bm_table) <- attributes
    genes_tx_id <- as.integer(factor(bm_table$ensembl_transcript_id,
                                     levels=transcripts_tx_name))
    message("OK")
    data.frame(
        tx_id=genes_tx_id,
        gene_id=bm_table$ensembl_gene_id
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Prepare the 'metadata' data frame.
###

.prepareBiomartMetadata <- function(mart, is_full_dataset, 
                            host="mar2009.archive.ensembl.org", 
                            path="/biomart/martservice", archive=FALSE)
{
    message("Prepare the 'metadata' data frame ... ",
            appendLF=FALSE)
    biomart <- biomaRt:::martBM(mart)
    dataset <- biomaRt:::martDataset(mart)
    db_version <- .getBiomartDbVersion(mart, host=host, path=path, 
                    biomart, archive=FALSE)
    datasets <- listDatasets(mart)
    dataset_rowidx <- which(as.character(datasets$dataset) == dataset)
    ## This should never happen (the above call to useMart() would have failed
    ## in the first place).
    if (length(dataset_rowidx) != 1L)
        stop("the BioMart database \"", biomaRt:::martBM(mart),
             "\" has no (or more than one) \"", dataset, "\" datasets")
    description <- as.character(datasets$description)[dataset_rowidx]
    dataset_version <- as.character(datasets$version)[dataset_rowidx]
    species <- .extractSpeciesFromDatasetDesc(description)
    message("OK")
    data.frame(
        name=c("Data source",
               "Genus and Species",
               "Resource URL",
               "BioMart database",
               "BioMart database version",
               "BioMart dataset",
               "BioMart dataset description",
               "BioMart dataset version",
               "Full dataset"),
        value=c("BioMart",
                species,
                "http://www.biomart.org/",
                biomart,
                db_version,
                dataset,
                description,
                dataset_version,
                ifelse(is_full_dataset, "yes", "no"))
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTxDbFromBiomart()
###

.parseBMMartParams <- function(biomart="ENSEMBL_MART_ENSEMBL",
                                      dataset="hsapiens_gene_ensembl",host="mar2009.archive.ensembl.org",path="/biomart/martservice",archive=FALSE)
{
    if (is.factor(biomart))
        biomart <- as.character(biomart)
    if (is(dataset, "AsIs"))
        dataset <- as.character(dataset)
    if (!isSingleString(biomart))
        stop("'biomart' must be a single string")
    useMart(biomart=biomart, dataset=dataset,host=host,path=path,archive=FALSE)
}

.parseBMFiltersParams <- function(transcript_ids)
{
    if (is.null(transcript_ids)) {
        filters <- ""
    } else if (is.character(transcript_ids)
            && !any(is.na(transcript_ids))) {
        filters <- "ensembl_transcript_id"
    }
    filters
}

.parseBMValuesParams <- function(transcript_ids)
{
    if (is.null(transcript_ids)) {
        values <- ""
    }else if (is.character(transcript_ids)
            && !any(is.na(transcript_ids))) {
        if (length(transcript_ids) == 0L)
            values <- "____a_very_unlikely_valid_transcript_id____"
        else
            values <- transcript_ids
    } else {
        stop("'transcript_ids' must be a character vector with no NAs")
    }
    values
}


## .testMakeTxDbFromBMParams <- function(biomart="ensembl",
##                                       dataset="hsapiens_gene_ensembl",
##                                       circ_seqs=DEFAULT_CIRC_SEQS,
##                                       transcript_ids=NULL)
## {
    ## if (is.factor(biomart))
    ##     biomart <- as.character(biomart)
    ## if (is(dataset, "AsIs"))
    ##     dataset <- as.character(dataset)
    ## if (!isSingleString(biomart))
    ##     stop("'biomart' must be a single string")
    ## mart <- useMart(biomart=biomart, dataset=dataset)

    ## if (is.null(transcript_ids)) {
    ##     filters <- values <- ""
    ## } else if (is.character(transcript_ids)
    ##         && !any(is.na(transcript_ids))) {
    ##     filters <- "ensembl_transcript_id"
    ##     if (length(transcript_ids) == 0L)
    ##         values <- "____a_very_unlikely_valid_transcript_id____"
    ##     else
    ##         values <- transcript_ids
    ## } else {
    ##     stop("'transcript_ids' must be a character vector with no NAs")
    ## }
## }


### Note that listMarts() and listDatasets() are returning data frames where
### the columns are character factors for the former and "AsIs" character
### vectors for the latter.

makeTranscriptDbFromBiomart_archive <- function(biomart="ENSEMBL_MART_ENSEMBL",
                                        dataset="hsapiens_gene_ensembl",
                                        transcript_ids=NULL,
                                        circ_seqs=DEFAULT_CIRC_SEQS,
                                        host="mar2009.archive.ensembl.org",path="/biomart/martservice",archive=FALSE)
{
    ## Could be that the user got the 'biomart' and/or 'dataset' values
    ## programmatically via calls to listMarts() and/or listDatasets().
    mart <- .parseBMMartParams(biomart=biomart,
                              dataset=dataset,host=host,path=path,archive=FALSE)
    filters <- .parseBMFiltersParams(transcript_ids)
    values <- .parseBMValuesParams(transcript_ids)

    transcripts <- .makeBiomartTranscripts(filters, values, mart,
                                           transcript_ids)
    chrominfo <- .makeBiomartChrominfo(mart,
                                       extra_seqnames=transcripts$tx_chrom,
                                       circ_seqs=circ_seqs,host=host,path=path,archive=FALSE)
    splicings <- .makeBiomartSplicings(filters, values, mart,
                                       transcripts$tx_name)
    genes <- .makeBiomartGenes(filters, values, mart, transcripts$tx_name)
    metadata <- .prepareBiomartMetadata(mart, is.null(transcript_ids),host=host,path=path,archive=FALSE)

    message("Make the TxDb object ... ", appendLF=FALSE)
    txdb <- makeTxDb(transcripts, splicings,
                             genes=genes, chrominfo=chrominfo,
                             metadata=metadata)
    message("OK")
    txdb
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some non-exported tools to help exploring/scanning the BioMart landscape.
###

### 'mart' can be either a Mart object or the name of a Mart service (single
### string). Returns a named list of 2-col data frames with one elt per
### dataset in 'mart'. Each data frame describes the attributes that are
### available for the corresponding dataset.
### Typical use:
###   ensembl_attrlist <- GenomicFeatures:::getMartAttribList("ensembl")
###   sapply(ensembl_attrlist, nrow)
getMartAttribList <- function(mart)
{
    if (!is(mart, "Mart"))
        mart <- useMart(mart)
    datasets <- listDatasets(mart)
    ans_length <- nrow(datasets)
    ans <- vector(mode="list", length=ans_length)
    names(ans) <- as.character(datasets$dataset)
    for (i in seq_len(ans_length)) {
        dataset <- names(ans)[i]
        mart <- useDataset(dataset, mart=mart)
        message("Getting attributes for dataset \"", dataset, "\"... ",
                appendLF=FALSE)
        ans[[i]] <- listAttributes(mart)
        message("OK")
    }
    ans
}

### 'biomart' and 'version' must be single character strings.
scanMart <- function(biomart, version)
{
    cat("Scanning ", biomart, "... ", sep="")
    suppressMessages(attrgroups <- .getAllDatasetAttrGroups(biomart))
    cat("OK\n")
    cat("biomart: ", biomart, "\n", sep="")
    cat("version: ", version, "\n", sep="")
    tmp <- names(attrgroups)
    if (length(tmp) > 3L)
        tmp <- c(tmp[1:3], "...")
    cat("nb of datasets: ", length(attrgroups),
        " (", paste(tmp, collapse=", "), ")\n",
        sep="")
    if (length(attrgroups) != 0L) {
        tbl <- table(attrgroups)
        tbl2 <- as.integer(tbl)
        names(tbl2) <- names(tbl)
        tmp <- paste(names(tbl2), ":", tbl2, sep="", collapse=", ")
        cat("table of attribute groups: ", tmp, "\n", sep="")
    }
    cat("\n")
}

scanMarts <- function(marts=NULL)
{
    if (is.null(marts))
        marts <- listMarts()
    biomarts <- as.character(marts$biomart)
    versions <- as.character(marts$version)
    for (i in seq_len(nrow(marts)))
        scanMart(biomarts[i], versions[i])
}

### scanMarts() output as of 6/28/2010 (only biomarts with at least groups
### A and G are listed):
###
### biomart: ensembl
### version: ENSEMBL GENES 58 (SANGER UK)
### nb of datasets: 51 (hsapiens_gene_ensembl, oanatinus_gene_ensembl,
###                     tguttata_gene_ensembl, cporcellus_gene_ensembl, ...)
### NOTE: the mgallopavo_gene_ensembl dataset seems to be broken!
### table of attribute groups: ABCDG:50
###
### biomart: bacterial_mart_5
### version: ENSEMBL BACTERIA 5 (EBI UK)
### nb of datasets: 183 (str_57_gene, esc_20_gene, myc_25994_gene, ...)
### table of attribute groups: ABG:183
###
### biomart: fungal_mart_5
### version: ENSEMBL FUNGAL 5 (EBI UK)
### nb of datasets: 12 (aniger_eg_gene, aflavus_eg_gene, aterreus_eg_gene, ...)
### table of attribute groups: ABG:12
###
### biomart: metazoa_mart_5
### version: ENSEMBL METAZOA 5 (EBI UK)
### nb of datasets: 23 (dgrimshawi_eg_gene, ppacificus_eg_gene,
###                     dpseudoobscura_eg_gene, ...)
### table of attribute groups: ABG:23
###
### biomart: plant_mart_5
### version: ENSEMBL PLANT 5 (EBI UK)
### nb of datasets: 8 (sbicolor_eg_gene, bdistachyon_eg_gene,
###                    alyrata_eg_gene, ...)
### table of attribute groups: ABG:8
###
### biomart: protist_mart_5
### version: ENSEMBL PROTISTS 5 (EBI UK)
### nb of datasets: 6 (tpseudonana_gene, ptricornutum_gene, pknowlesi_gene, ...)
### table of attribute groups: ABG:6
###
### biomart: ensembl_expressionmart_48
### version: EURATMART (EBI UK)
### nb of datasets: 1 (rnorvegicus_expr_gene_ensembl)
### table of attribute groups: AG:1
###
### biomart: Ensembl56
### version: PANCREATIC EXPRESSION DATABASE (INSTITUTE OF CANCER UK)
### nb of datasets: 1 (hsapiens_gene_pancreas)
### table of attribute groups: ABCDG:1
###
### biomart: ENSEMBL_MART_ENSEMBL
### version: GRAMENE 30 ENSEMBL GENES (CSHL/CORNELL US)
### nb of datasets: 8 (sbicolor_eg_gene, bdistachyon_eg_gene,
###                    alyrata_eg_gene, ...)
### table of attribute groups: ABG:8

