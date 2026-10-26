#' @importFrom utils getFromNamespace
NULL

## ---------------------------------------------------------------------------
## Internal functions borrowed from DOSE / clusterProfiler.
##
## They are deliberately NOT resolved at package-load time any more.
## Doing `foo <- getFromNamespace("foo", "DOSE")` at the top level of an R
## file means the value is computed while the package is being lazily loaded,
## so as soon as DOSE / clusterProfiler renames, relocates or removes one of
## these internal functions, *installing* GseaVis fails with
##
##     Error in get(x, envir = ns, inherits = FALSE) :
##       object 'xxx' not found
##
## Resolving them lazily (at call time) keeps GseaVis installable, and only
## the affected function raises an informative error if it is really missing.
## ---------------------------------------------------------------------------

## Look `fun` up in the internal namespace of the first package that has it.
gsv_internal <- function(fun, pkgs) {
    for (pkg in pkgs) {
        f <- tryCatch(getFromNamespace(fun, pkg), error = function(e) NULL)
        if (is.function(f)) {
            return(f)
        }
    }
    stop(sprintf(
        "GseaVis needs the internal function '%s' from package '%s', but it could not be found.\n  Please update '%s' (and GseaVis) to their latest versions.",
        fun, pkgs[1], pkgs[1]
    ), call. = FALSE)
}

prepare_KEGG <- function(...) {
    gsv_internal("prepare_KEGG", c("clusterProfiler", "DOSE"))(...)
}

get_GO_data <- function(...) {
    gsv_internal("get_GO_data", c("clusterProfiler", "DOSE"))(...)
}

getGeneSet <- function(...) {
    gsv_internal("getGeneSet", c("DOSE", "clusterProfiler"))(...)
}

gson_cpd <- function(...) {
    gsv_internal("gson_cpd", c("clusterProfiler", "DOSE"))(...)
}

organismMapper <- function(...) {
    gsv_internal("organismMapper", c("clusterProfiler", "DOSE"))(...)
}

get_data_from_KEGG_db <- function(...) {
    gsv_internal("get_data_from_KEGG_db", c("clusterProfiler", "DOSE"))(...)
}

build_Anno <- function(...) {
    gsv_internal("build_Anno", c("DOSE", "clusterProfiler"))(...)
}

is.sorted <- function(...) {
    gsv_internal("is.sorted", c("DOSE", "clusterProfiler"))(...)
}

## Turn a character OrgDb name (e.g. "org.Hs.eg.db") into the OrgDb object.
gsv_load_orgdb <- function(OrgDb) {
    if (!is.character(OrgDb)) {
        return(OrgDb)
    }
    nm <- OrgDb[1]
    obj <- tryCatch(get(nm, envir = asNamespace(nm)), error = function(e) NULL)
    if (is.null(obj)) {
        stop(sprintf("Cannot find the OrgDb object '%s', please check that the package is installed.", nm),
             call. = FALSE)
    }
    obj
}

## Species name of an OrgDb object, e.g. "Homo sapiens".
## This used to be DOSE:::get_organism(), which is an *internal* (not exported)
## function and has been dropped in several DOSE releases, so it is
## re-implemented here to keep GseaVis independent of it.
get_organism <- function(OrgDb) {
    OrgDb <- gsv_load_orgdb(OrgDb)

    org <- tryCatch(AnnotationDbi::species(OrgDb), error = function(e) NA_character_)

    if (length(org) == 0 || is.na(org[1])) {
        ## fall back to the ORGANISM field of the database metadata
        md <- tryCatch(as.data.frame(AnnotationDbi::metadata(OrgDb)),
                       error = function(e) NULL)
        if (!is.null(md) && all(c("name", "value") %in% names(md))) {
            v <- md$value[md$name == "ORGANISM"]
            if (length(v) > 0) {
                org <- v[1]
            }
        }
    }

    if (length(org) == 0 || is.na(org[1])) {
        ## last resort: derive it from the sqlite file name (org.Hs.eg.sqlite)
        org <- sub("\\.sqlite$", "",
                   basename(tryCatch(AnnotationDbi::dbfile(OrgDb),
                                     error = function(e) NA_character_)))
    }

    as.character(org[1])
}
