check_blocks <- function(blocks, init = FALSE, n = 2,
                         add_NAlines = FALSE, allow_unnames = TRUE,
                         quiet = FALSE, no_character = FALSE) {
  msg <- ""
  if (is.matrix(blocks)) blocks <- list(blocks)
  if (!is.list(blocks)) stop_rgcca(paste(msg, "is not a list."))
  if (!init && length(blocks) < n) {
    stop_rgcca(paste(msg, "should at least have two blocks."))
  }

  # Completing block names
  if (is.null(names(blocks)) | any(names(blocks) == "")) {
    if (any(names(blocks) == "")) {
      for (x in which(names(blocks) == "")) {
        names(blocks)[x] <- paste0("block", x)
      }
    } else {
      names(blocks) <- paste0("block", seq_along(blocks))
    }
    message("Missing block names are automatically labeled.")
  }

  # Gestion of the case of one variable only
  blocks <- lapply(blocks, as.matrix)
  name_blocks <- names(blocks)

  # Dealing with rownames (if they are all missing)
  if (all(sapply(blocks, function(x) is.null(row.names(x))))) {
    if (sd(sapply(blocks, function(x) NROW(x))) == 0 && allow_unnames) {
      blocks <- lapply(
        blocks,
        function(x) {
          rownames(x) <- paste0("S", 1:(NROW(x)))
          return(x)
        }
      )
      message("Missing rownames are automatically labeled.")
    } else {
      stop_rgcca(paste(msg, "Blocks should have rownames.\n "))
    }
  }

  if (any(sapply(blocks, function(x) is.null(colnames(x))))) {
    message("Missing colnames are automatically labeled.")
    blocks1 <- lapply(
      seq_along(blocks),
      function(x) {
        if (is.null(colnames(blocks[[x]]))) {
          if (NCOL(blocks[[x]]) == 1) {
            colnames(blocks[[x]]) <- names(blocks)[x]
          } else {
            colnames(blocks[[x]]) <-
              paste0("V", x, seq(NCOL(blocks[[x]])))
          }
          return(blocks[[x]])
        } else {
          blocks[[x]] <- blocks[[x]]
          return(blocks[[x]])
        }
      }
    )
    blocks <- blocks1
    names(blocks) <- name_blocks
  }




  # if one colname is identical within or between block
  if (sum(duplicated(unlist(sapply(blocks, colnames)))) != 0) {
    if (!quiet) {
      message("Duplicated colnames are modified to avoid confusion \n")
    }

    blocks_i <- lapply(
      seq_along(blocks),
      function(i) {
        x <- blocks[[i]]
        colnames(x) <- paste(names(blocks)[i],
          colnames(blocks[[i]]),
          sep = "_"
        )
        return(x)
      }
    )
    names(blocks_i) <- names(blocks)
    blocks <- blocks_i
  }

  # If one rownames is missing but the size of blocks is correct
  if (any(sapply(blocks, function(x) is.null(row.names(x))))) {
    matrix_of_rownames <- Reduce(cbind, lapply(blocks, row.names))
    if (sum(!apply(
      matrix_of_rownames, 2,
      function(x) x == matrix_of_rownames[, 1]
    )) == 0) {
      blocks <- lapply(
        blocks,
        function(x) {
          row.names(x) <- matrix_of_rownames[, 1]
          return(x)
        }
      )
    }
  }

  lapply(
    blocks,
    function(x) {
      resdup <- duplicated(rownames(x))
      if (sum(resdup) != 0) {
        if (!quiet) {
          warning(paste0(
            "Duplicated rownames were removed: ",
            rownames(x)[resdup], "\n"
          ))
        }
      }
    }
  )
  inters_rows <- Reduce(intersect, lapply(blocks, row.names))

  if (length(inters_rows) == 0) {
    stop_rgcca(paste(msg, "Elements of the list should have at least
                         one common rowname.\n "))
  }

  equal_rows <- Reduce(identical, lapply(blocks, row.names))

  # If add_NAlines=FALSE, taking the intersection_list
  if (!add_NAlines) {
    if (length(blocks) > 1 && !equal_rows) blocks <- common_rows(blocks)
  }

  if (init) {
    blocks <- remove_null_sd(blocks)$list_m
    for (i in seq(length(blocks))) {
      attributes(blocks[[i]])$nrow <- nrow(blocks[[i]])
    }
  }

  if (no_character) {
    if (any(sapply(blocks, is.character2))) {
      stop(paste(msg, "Blocks contain non-numeric values."))
    }

    for (i in seq(length(blocks))) {
      if (is.character(blocks[[i]])) {
        blocks[[i]] <- to_numeric(blocks[[i]])
      }
    }
  }

  # Add lines if subjects are missing
  if (add_NAlines) {
    union_rows <- Reduce(union, lapply(blocks, row.names))
    blocks2 <- lapply(name_blocks, function(name) {
      # if some subjects are missing (in the rownames)
      if (sum(!union_rows %in% rownames(blocks[[name]])) != 0) {
        message("Some subjects are blockwise missing and NA rows were added.")
        y <- matrix(NA,
          length(union_rows),
          ncol = ifelse(is.null(dim(blocks[[name]])),
            1,
            NCOL(blocks[[name]])
          )
        )
        if (is.null(dim(blocks[[name]]))) {
          colnames(y) <- name
        } else {
          colnames(y) <- colnames(blocks[[name]])
        }
        rownames(y) <- union_rows
        y[rownames(blocks[[name]]), ] <- blocks[[name]]
        return(y)
      } else {
        if (NCOL(blocks[[name]]) == 1) {
          y <- matrix(blocks[[name]][union_rows, ], ncol = 1)
          rownames(y) <- union_rows
          colnames(y) <- name
        } else {
          y <- blocks[[name]][union_rows, ]
        }
        return(y)
      }
    })
    names(blocks2) <- name_blocks
    blocks <- blocks2
  }
  invisible(blocks)
}
