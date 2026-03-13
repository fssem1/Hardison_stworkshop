# helper function for turning sf geometry into lat/lon columns
# Originally written by Josh London here -> https://github.com/r-spatial/sf/issues/231
sfc_as_cols <- function (x, geometry, names = c("x", "y"))
{
  if (missing(geometry)) {
    geometry <- sf::st_geometry(x)
  }
  else {
    geometry <- rlang::eval_tidy(enquo(geometry), x)
  }
  stopifnot(inherits(x, "sf") && inherits(geometry, "sfc_POINT"))
  ret <- sf::st_coordinates(geometry)
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[, !names(x) %in% names]
  ret <- setNames(ret, names)
  dplyr::bind_cols(x, ret)
}

theme_fade <- function(...)
  {
  ggplot2::theme(strip.background = element_blank(), panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), panel.background = element_blank(),
                 panel.border = element_rect(colour = "black", fill = NA,
                                             size = 0.2), legend.key = element_blank(), axis.title = element_text(size = 10))
  }

fct_to_num <-
  function (x)
  {
    as.numeric(levels(x))[x]
  }

theme_minimal2 <-
  function (base_size = 11, base_family = "", header_family = NULL, 
            base_line_size = base_size/22, base_rect_size = base_size/22, 
            ink = "black", paper = alpha(ink, 0), accent = "#3366FF") {
    half_line <- base_size/2
    t <- theme(line = element_blank(), rect = element_rect(fill = paper, 
                                                           colour = NA, linewidth = 0, linetype = 1, inherit.blank = FALSE, 
                                                           linejoin = "round"), polygon = element_blank(), point = element_blank(), 
               text = element_text(family = base_family, face = "plain", 
                                   colour = ink, size = base_size, lineheight = 0.9, 
                                   hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), 
                                   debug = FALSE), title = element_text(family = header_family), 
               spacing = unit(half_line, "pt"), margins = margin_auto(half_line), 
               geom = element_geom(ink = ink, paper = paper, accent = accent, 
                                   linewidth = base_line_size, borderwidth = base_line_size, 
                                   linetype = 1L, bordertype = 1L, family = base_family, 
                                   fontsize = base_size, pointsize = (base_size/11) * 
                                     1.5, pointshape = 19), axis.text = element_blank(), 
               axis.ticks.length = rel(0), 
               axis.ticks.length.x = NULL, axis.ticks.length.x.top = NULL, 
               axis.ticks.length.x.bottom = NULL, axis.ticks.length.y = NULL, 
               axis.ticks.length.y.left = NULL, axis.ticks.length.y.right = NULL, 
               axis.minor.ticks.length = NULL, legend.box = NULL, legend.key.size = unit(1.2, 
                                                                                         "lines"), legend.position = "right", legend.text = element_text(size = rel(0.8)), 
               legend.title = element_text(hjust = 0), legend.key.spacing = rel(1), 
               legend.margin = margin_auto(0), legend.box.margin = margin_auto(0), 
               legend.box.spacing = unit(0.2, "cm"), legend.ticks.length = rel(0.2), 
               legend.background = element_blank(), legend.frame = element_blank(), 
               legend.box.background = element_blank(), strip.clip = "on", 
               strip.text = element_text(size = rel(0.8)), strip.switch.pad.grid = rel(0.5), 
               strip.switch.pad.wrap = rel(0.5), strip.background = element_blank(), 
               panel.ontop = FALSE, panel.spacing = NULL, panel.background = element_blank(), 
               panel.border = element_blank(), plot.margin = margin_auto(0), 
               plot.title = element_text(size = rel(1.2), hjust = 0, 
                                         vjust = 1, margin = margin(t = half_line)), plot.title.position = "panel", 
               plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(t = half_line)), 
               plot.caption = element_text(size = rel(0.8), hjust = 1, 
                                           vjust = 1, margin = margin(t = half_line)), plot.caption.position = "panel", 
               plot.tag = element_text(size = rel(1.2), hjust = 0.5, 
                                       vjust = 0.5), plot.tag.position = "topleft", plot.background = element_rect(), 
               complete = TRUE)
  }


# This is the sanity function for tinyVAST that's currently in the dev branch. Written by James Thorson. 
# I have edited this script to replace g <- object$internal$convergence_diagnostics$final_grads
# with g <- object$sdrep$gradient.fixed that is the sdmTMB equivalent

#' Sanity check of a tinyVAST model
#'
#' @param object Fitted model from [tinyVAST()].
#' @param big_sd_log10 Value to check size of standard errors against. A value
#'   of 2 would indicate that standard errors greater than `10^2` (i.e., 100)
#'   should be flagged.
#' @param gradient_thresh Gradient threshold to issue warning.
#' @param silent Logical: suppress messages? Useful to set to `TRUE` if running
#'   large numbers of models and just interested in returning sanity list
#'   objects.
#'
#' @return An invisible named list of checks.
#' @export
sanity_tv <- 
  function( object, 
            big_sd_log10 = 2, 
            gradient_thresh = 0.001, 
            silent = FALSE) {
    
    # assertClass(object, "tinyVAST")
    hessian_ok <- eigen_values_ok <- gradients_ok <- se_magnitude_ok <- FALSE
    nlminb_ok <- FALSE
    
    simplify_msg <- "Try simplifying the model, adjusting the mesh, or adding priors"
    
    if (identical(object$opt$convergence, 0L)) {
      msg <- "Non-linear minimizer suggests successful convergence"
      if (!silent) cli::cli_alert_success(msg)
      nlminb_ok <- TRUE
    } else {
      msg <- "Non-linear minimizer did not converge: do not trust this model!"
      if (!silent) cli::cli_alert_danger(msg)
      if (!silent) cli::cli_alert_info(simplify_msg)
      cat("\n")
    }
    
    if( is.null(object$sdrep$pdHess) ){
      msg <- "Non-positive-definite Hessian matrix: model may not have converged"
    }else if (isFALSE(object$sdrep$pdHess)) {
      msg <- "Hessian matrix not detected: consider re-running to check if positive definite"
      if (!silent) cli::cli_alert_danger(msg)
      if (!silent) cli::cli_alert_info(simplify_msg)
      cat("\n")
    } else {
      msg <- "Hessian matrix is positive definite"
      if (!silent) cli::cli_alert_success(msg)
      hessian_ok <- TRUE
    }
    
    if (isTRUE(object$convergence_diagnostics$bad_eig)) {
      msg <- "Extreme or very small eigenvalues detected: model may not have converged"
      if (!silent) cli::cli_alert_danger(msg)
      if (!silent) cli::cli_alert_info(simplify_msg)
      cat("\n")
    } else {
      msg <- "No extreme or very small eigenvalues detected"
      if (!silent) cli::cli_alert_success(msg)
      eigen_values_ok <- TRUE
    }

    g <- object$sdrep$gradient.fixed # this is equivalent, 
    # g <- object$internal$convergence_diagnostics$final_grads #not working, prob bc not running dev branch
    np <- names(object$obj$par)
    for (i in seq_along(g)) {
      if (abs(g[i]) > gradient_thresh) {
         if (!silent) {
           cli::cli_alert_danger(c(
             "`", np[i],
             paste0("` gradient > ", gradient_thresh)
           ))
         }
         msg <- "standardize covariates, and/or simplify the model"
         if (!silent) {
           cli::cli_alert_info(msg)
           cat("\n")
         }
       }
     }

    if (all(abs(g) <= gradient_thresh)) {
      msg <- "No gradients with respect to fixed effects are >= "
      if (!silent) cli::cli_alert_success(paste0(msg, gradient_thresh))
      gradients_ok <- TRUE
    }

    obj <- object$tmb_obj
    random <- unique(names(obj$env$par[obj$env$random]))
    
    pl <- as.list(object$sdrep, "Estimate")
    ple <- as.list(object$sdrep, "Std. Error")
    pars <- names(object$obj$par)
    pl <- pl[pars]
    ple <- ple[pars]
    np <- names(ple)
    se_na_ok <- TRUE
    for (i in seq_along(ple)) {
      if (any(is.na(ple[i]))) {
        if (!silent) cli::cli_alert_danger(c("`", np[i], paste0("` standard error is NA")))
        #par_message(np[i])
        if (!silent) {
          cli::cli_alert_info(simplify_msg)
          cat("\n")
        }
        se_na_ok <- FALSE
      }
    }
    if (se_na_ok) {
      msg <- "No fixed-effect standard errors are NA"
      if (!silent) cli::cli_alert_success(msg)
    }
    
    est <- as.list(object$sd_report, "Estimate")
    se <- as.list(object$sd_report, "Std. Error")
    fixed <- !(names(est) %in% random)
    est <- est[fixed]
    se <- se[fixed]
    
    too_big <- function(se) {
      if (any(!is.na(se))) {
        se_max <- max(se, na.rm = TRUE)
        if (any(log10(abs(se_max)) > big_sd_log10)) {
          return(TRUE)
        } else {
          return(NULL)
        }
      } else {
        return(NULL)
      }
    }
    
    se_big <- lapply(se, too_big)
    
    for (i in seq_along(se_big)) {
      if (isTRUE(se_big[[i]])) {
        msg <- paste0("` standard error may be large")
        if (!silent) cli::cli_alert_danger(c("`", names(se_big)[i], msg))
        #par_message(names(se_big)[i])
        if (!silent) {
          cli::cli_alert_info(simplify_msg)
          cat("\n")
        }
      }
    }
    
    if (all(unlist(lapply(se_big, is.null)))) {
      msg <- "No standard errors look unreasonably large"
      if (!silent) cli::cli_alert_success(msg)
      se_magnitude_ok <- TRUE
    }
    
    #  # tidying objects with different names -- fixed won't have model or group_name
    #  b <- tidy(object, conf.int = TRUE, silent = TRUE)
    #  br <- tidy(object, "ran_pars", conf.int = TRUE, silent = TRUE)
    #  b[, names(br)[!names(br) %in% names(b)]] <- NA
    #  b <- rbind(b, br)
    #
    #  if (isTRUE(object$family$delta)) {
    #    b2 <- tidy(object, conf.int = TRUE, model = 2, silent = TRUE)
    #    br2 <- tidy(object, "ran_pars", conf.int = TRUE, model = 2, silent = TRUE)
    #    b2[, names(br2)[!names(br2) %in% names(b2)]] <- NA
    #    b2 <- rbind(b2, br2)
    #
    #    # also check for missing names -- will occur when one model has no random effects
    #    b[, names(b2)[!names(b2) %in% names(b)]] <- NA
    #    b2[, names(b)[!names(b) %in% names(b2)]] <- NA
    #    b <- rbind(b, b2)
    #  }
    #  s <- grep("sigma", b$term)
    #  sigmas_ok <- TRUE
    #  if (length(s)) {
    #    for (i in s) {
    #      if (b$estimate[i] < 0.01) {
    #        msg <- "` is smaller than 0.01"
    #        if (!silent) cli::cli_alert_danger(c("`", b$term[i], msg))
    #        par_message(b$term[i])
    #        msg <- "Consider omitting this part of the model"
    #        if (!silent) {
    #          cli::cli_alert_info(msg)
    #          cat("\n")
    #        }
    #        sigmas_ok <- FALSE
    #      }
    #    }
    #  }
    #  if (sigmas_ok) {
    #    msg <- "No sigma parameters are < 0.01"
    #    if (!silent) cli::cli_alert_success(msg)
    #  }
    #
    #  if (length(s)) {
    #    for (i in s) {
    #      if (b$estimate[i] > 100) {
    #        msg <- "` is larger than 100"
    #        if (!silent) cli::cli_alert_danger(c("`", b$term[i], msg))
    #        par_message(b$term[i])
    #        msg <- "Consider simplifying the model or adding priors"
    #        if (!silent) {
    #          cli::cli_alert_info(msg)
    #          cat("\n")
    #        }
    #        sigmas_ok <- FALSE
    #      }
    #    }
    #  }
    #  if (sigmas_ok) {
    #    msg <- "No sigma parameters are > 100"
    #    if (!silent) cli::cli_alert_success(msg)
    #  }
    
    #  r1 <- diff(range(object$data[[object$spde$xy_cols[1]]]))
    #  r2 <- diff(range(object$data[[object$spde$xy_cols[2]]]))
    #  r <- max(r1, r2)
    #  range_ok <- TRUE
    #  if ("range" %in% b$term) {
    #    if (max(b$estimate[b$term == "range"]) > r * 1.5) {
    #      msg <- "A `range` parameter looks fairly large (> 1.5 the greatest distance in data)"
    #      if (!silent) {
    #        cli::cli_alert_danger(msg)
    #        cli::cli_alert_info(simplify_msg)
    #        cli::cli_alert_info("Also make sure your spatial coordinates are not too big or small (e.g., work in UTM km instead of UTM m)", wrap = TRUE)
    #        cat("\n")
    #      }
    #      range_ok <- FALSE
    #    } else {
    #      nr <- length(grep("range", b$term))
    #      if (nr == 1L) msg <- "Range parameter doesn't look unreasonably large"
    #      if (nr > 1L) msg <- "Range parameters don't look unreasonably large"
    #      if (!silent) cli::cli_alert_success(msg)
    #    }
    #  }
    # 
    # ret <- named_list(
    #   hessian_ok, eigen_values_ok, nlminb_ok, 
    #   gradients_ok, se_magnitude_ok, se_na_ok 
    # )
    
    ret <- sdmTMB:::named_list(
      hessian_ok, 
      eigen_values_ok, 
      nlminb_ok, 
      gradients_ok, 
      se_magnitude_ok, 
      se_na_ok 
    )
    
    all_ok <- all(unlist(ret))
    ret <- c(ret, all_ok = all_ok)
    invisible(ret)
  }
