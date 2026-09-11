#' Copernicus cache directory
#'
#' Returns the path to the package's persistent cache directory for downloaded
#' Copernicus environmental data. Uses [tools::R_user_dir()] so the location
#' survives across sessions and follows platform conventions.
#'
#' @returns Character. Path to cache directory (created if missing).
#' 
#' 
#' @examples
#' copernicus_cache_dir()
#'
#' @export
copernicus_cache_dir <- function() {
  
  path <- file.path(
    tools::R_user_dir(
      "ocean3d",
      which = "cache"
    ),
    "copernicus"
  )
  
  if (!dir.exists(path)) {
    dir.create(
      path,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  path
}

#' Clear the Copernicus cache
#'
#' Remove all files downloaded to the default Copernicus cache.
#'
#' @param confirm Logical. Require interactive confirmation. Default `TRUE`.
#'
#' @returns Invisibly, `TRUE` on success.
#' 
#' @examples
#' \dontrun{
#' copernicus_cache_clear()
#' }
#' 
#' @export
copernicus_cache_clear <- function(confirm = TRUE) {
  path <- copernicus_cache_dir()

  if (confirm && interactive()) {
    ans <- readline(sprintf("Delete all cached Copernicus files in %s? [y/N]: ", path))
    if (!tolower(ans) %in% c("y", "yes")) {
      message("Cancelled.")
      return(invisible(FALSE))
    }
  }

  unlink(path, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

# Internal: one-time interactive consent for writing to persistent cache.
.copernicus_cache_consent <- function(cache_dir) {
  sentinel <- file.path(cache_dir, ".consent")
  if (file.exists(sentinel)) return(invisible(TRUE))

  msg <- paste0(
    "This package will cache downloaded Copernicus files in:\n  ",
    cache_dir, "\n",
    "Copernicus environmental datasets can be large.\n",
    "Pass `output_dir` to use a different location, or call ",
    "copernicus_cache_clear() to reclaim space later."
  )

  if (!interactive()) {
    stop(
      msg,
      "\nNon-interactive session: pass `output_dir` explicitly, or run ",
      "copernicus_load() once interactively to record consent.",
      call. = FALSE
    )
  }

  message(msg)
  ans <- readline("Proceed with caching here? [y/N]: ")
  if (!tolower(ans) %in% c("y", "yes")) stop("Cache consent declined.", call. = FALSE)

  file.create(sentinel)
  invisible(TRUE)
}

# -------------------------------------------------------------------------
# Copernicus Marine Toolbox constants
# -------------------------------------------------------------------------

# Copernicus Marine Toolbox standalone version tested by this package.
.COPERNICUS_TOOLBOX_VERSION <- "2.4.1"

.COPERNICUS_TOOLBOX_INFO_URL <-
  "https://help.marine.copernicus.eu/en/collections/4060068-copernicus-marine-toolbox"

.COPERNICUS_TOOLBOX_REPOSITORY <-
  "https://github.com/mercator-ocean/copernicus-marine-toolbox"

.COPERNICUS_TOOLBOX_CREDENTIALS_URL <-
  "https://help.marine.copernicus.eu/en/articles/8185007-copernicus-marine-toolbox-credentials-configuration"

.COPERNICUS_REGISTER_URL <-
  "https://data.marine.copernicus.eu/register"

# Internal: return the package-managed Copernicus Marine Toolbox directory.
.copernicus_toolbox_dir <- function(
    version = .COPERNICUS_TOOLBOX_VERSION
) {
  .copernicus_validate_string(version, "version")
  
  file.path(
    tools::R_user_dir("ocean3d", which = "data"),
    "copernicusmarine",
    version
  )
}

# Internal: find the Copernicus Marine Toolbox standalone executable.
#
# Search order:
#   1. Explicit path supplied by the caller.
#   2. Package-managed standalone executable.
#   3. Executable available on the system PATH.
#
# Returns NULL when no usable executable can be found.
.copernicus_find_executable <- function(path = NULL) {
  executable_name <- if (.Platform$OS.type == "windows") {
    "copernicusmarine.exe"
  } else {
    "copernicusmarine"
  }
  
  # 1. Explicit path supplied by the user.
  if (!is.null(path)) {
    .copernicus_validate_string(path, "path")
    
    path <- path.expand(path)
    
    if (!file.exists(path)) {
      stop(
        "Copernicus Marine Toolbox executable not found at:\n  ",
        path,
        call. = FALSE
      )
    }
    
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  
  # 2. Package-managed standalone executable.
  managed <- file.path(
    .copernicus_toolbox_dir(),
    executable_name
  )
  
  if (file.exists(managed)) {
    return(normalizePath(managed, winslash = "/", mustWork = TRUE))
  }
  
  # 3. Executable already available on the system PATH.
  on_path <- Sys.which("copernicusmarine")
  
  if (nzchar(on_path)) {
    return(normalizePath(on_path, winslash = "/", mustWork = TRUE))
  }
  
  NULL
}


# Internal: return the version reported by a Copernicus Marine Toolbox
# standalone executable.
.copernicus_executable_version <- function(path) {
  .copernicus_validate_string(path, "path")
  
  path <- path.expand(path)
  
  if (!file.exists(path)) {
    return(NA_character_)
  }
  
  output <- tryCatch(
    system2(
      path,
      args = "--version",
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(output) || !length(output)) {
    return(NA_character_)
  }
  
  status <- attr(output, "status")
  
  if (!is.null(status) && !identical(as.integer(status), 0L)) {
    return(NA_character_)
  }
  
  output <- paste(output, collapse = " ")
  
  match <- regexpr(
    "[0-9]+\\.[0-9]+\\.[0-9]+",
    output,
    perl = TRUE
  )
  
  if (match[1] == -1) {
    return(NA_character_)
  }
  
  regmatches(output, match)
}


# Internal: check whether the Copernicus Marine Toolbox can find valid
# credentials using its own authentication configuration.
#
# Credential discovery is delegated entirely to the official Toolbox. This
# package does not search for or inspect credential files, environment
# variables, usernames, or passwords.
.copernicus_credentials_valid <- function(executable) {
  
  .copernicus_validate_string(executable, "executable")
  
  executable <- path.expand(executable)
  
  if (!file.exists(executable)) {
    stop(
      "The Copernicus Marine Toolbox executable does not exist at:\n  ",
      executable,
      call. = FALSE
    )
  }
  
  result <- tryCatch(
    suppressWarnings(
      system2(
        executable,
        args = c(
          "login",
          "--check-credentials-valid",
          "--log-level",
          "QUIET"
        ),
        stdout = TRUE,
        stderr = TRUE
      )
    ),
    error = function(e) {
      stop(
        "The Copernicus Marine Toolbox could not check the current ",
        "authentication configuration.\n\n",
        "Executable:\n  ",
        executable,
        "\n\n",
        "Original error:\n  ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  status <- attr(result, "status")
  
  # system2() does not attach a status attribute when the command exits
  # successfully.
  if (is.null(status)) {
    status <- 0L
  }
  
  if (identical(as.integer(status), 0L)) {
    return(TRUE)
  }
  
  # The Toolbox uses a non-zero status when no valid credentials are
  # available. This is not an R error: it means interactive login is needed.
  FALSE
}


#' Check Copernicus Marine Toolbox availability
#'
#' Check whether the official Copernicus Marine Toolbox standalone executable
#' is available and report the executable path and version.
#'
#' The executable is searched for in the package-managed Copernicus directory
#' and on the system `PATH`. An explicit executable path can also be supplied.
#'
#' @param path Optional character. Explicit path to a Copernicus Marine Toolbox
#'   executable.
#'
#' @returns A named list containing `available`, `path`, and `version`.
#'
#' @examples
#' copernicus_status()
#'
#' @export
copernicus_status <- function(path = NULL) {
  
  executable <- .copernicus_find_executable(path)
  
  if (is.null(executable)) {
    result <- list(
      available = FALSE,
      path = NULL,
      version = NULL
    )
    
    message(
      "Copernicus Marine Toolbox standalone executable was not found.\n",
      "Run `copernicus_setup()` to install the supported standalone Toolbox."
    )
    
    return(invisible(result))
  }
  
  version <- .copernicus_executable_version(executable)
  
  result <- list(
    available = TRUE,
    path = executable,
    version = version
  )
  
  message(
    "Copernicus Marine Toolbox\n",
    "  Available: yes\n",
    "  Version:   ",
    if (is.na(version)) "unknown" else version,
    "\n",
    "  Path:      ",
    executable
  )
  
  invisible(result)
}


#' Log in to Copernicus Marine
#'
#' Configure or verify authentication for Copernicus Marine using the official
#' standalone Copernicus Marine Toolbox.
#'
#' `copernicus_login()` delegates credential discovery, validation, and
#' storage entirely to the official Copernicus Marine Toolbox. This package
#' does not search for, request, receive, read, or store the user's Copernicus
#' Marine username or password.
#'
#' If valid credentials are already available, the Toolbox reuses them
#' automatically. If authentication still needs to be configured, an
#' interactive Toolbox login is started.
#'
#' Authentication only needs to be configured once. The credentials created
#' by the official Copernicus Marine Toolbox are shared by both the standalone
#' and Python backends, so no separate Python login is required before using
#' `backend = "python"`.
#' 
#' When running inside RStudio, a first-time login is opened in the RStudio
#' Terminal because username and password prompts from external programs may
#' not receive interactive input correctly from the RStudio Console.
#'
#' @param path Optional character string giving the path to a Copernicus Marine
#'   Toolbox standalone executable. If `NULL`, the package-managed executable
#'   and then the system `PATH` are searched automatically.
#'
#' @returns Invisibly, `TRUE` when valid Copernicus Marine credentials were
#'   already available and successfully verified. Invisibly, `FALSE` when an
#'   interactive login has been launched and must be completed by the user.
#'
#' @examples
#' \dontrun{
#' copernicus_login()
#' }
#'
#' @export
copernicus_login <- function(path = NULL) {
  
  if (!interactive()) {
    stop(
      "`copernicus_login()` requires an interactive R session when ",
      "authentication has not already been configured.\n\n",
      "A first-time Copernicus Marine login may require a username and ",
      "password to be entered interactively.\n\n",
      "For automated or non-interactive workflows, configure Copernicus Marine ",
      "authentication separately before making Copernicus Marine requests.\n\n",
      "For more information, see:\n  ",
      .COPERNICUS_TOOLBOX_CREDENTIALS_URL,
      call. = FALSE
    )
  }
  
  executable <- .copernicus_find_executable(path)
  
  if (is.null(executable)) {
    stop(
      "The Copernicus Marine Toolbox standalone executable could not be found.\n\n",
      "`copernicus_login()` relies on the official Toolbox to manage ",
      "authentication, but no usable executable was detected.\n\n",
      "Run `copernicus_setup()` to install the supported standalone Toolbox, ",
      "or provide the path to an existing Copernicus Marine Toolbox ",
      "installation using `path`.\n\n",
      "For more information about the Toolbox, see:\n  ",
      .COPERNICUS_TOOLBOX_INFO_URL,
      call. = FALSE
    )
  }
  
  # Ask the official Toolbox whether it can already find valid credentials.
  credentials_valid <- .copernicus_credentials_valid(executable)
  
  if (credentials_valid) {
    
    message(
      "Valid Copernicus Marine credentials are already available.\n\n",
      "The official Copernicus Marine Toolbox will confirm the authentication ",
      "configuration it is using below. No username or password should need ",
      "to be entered.\n"
    )
    
    status <- tryCatch(
      system2(
        executable,
        args = "login"
      ),
      error = function(e) {
        stop(
          "Valid credentials were detected, but the Copernicus Marine ",
          "Toolbox login command could not be completed.\n\n",
          "Executable:\n  ",
          executable,
          "\n\n",
          "Original error:\n  ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    
    if (!identical(as.integer(status), 0L)) {
      stop(
        "Copernicus Marine authentication could not be confirmed.\n\n",
        "Valid credentials were detected initially, but the official Toolbox ",
        "subsequently returned exit status ",
        status,
        ".\n\n",
        "This may indicate that the authentication configuration changed, ",
        "the network connection was interrupted, or the Copernicus Marine ",
        "authentication service could not be reached.\n\n",
        "Run `copernicus_login()` again. If the problem persists, see:\n  ",
        .COPERNICUS_TOOLBOX_CREDENTIALS_URL,
        call. = FALSE
      )
    }
    
    message(
      "\nCopernicus Marine authentication is ready.\n\n",
      "The official Toolbox has confirmed that valid credentials are ",
      "available for future Copernicus Marine requests.\n\n",
      "Credential discovery, validation, and storage are handled entirely by ",
      "the official Copernicus Marine Toolbox. This R package does not search ",
      "for, read, receive, or store your username or password.\n"
    )
    
    return(invisible(TRUE))
  }
  
  # -----------------------------------------------------------------------
  # No valid credentials were found.
  # -----------------------------------------------------------------------
  
  message(
    "No valid Copernicus Marine credentials are currently available.\n\n",
    
    "Copernicus Marine requires authentication before data can be downloaded. ",
    "The official Copernicus Marine Toolbox will therefore ask for your ",
    "Copernicus Marine username and password.\n\n",
    
    "Your credentials will be entered directly into the official Toolbox and ",
    "managed by the Toolbox for future requests. This R package does not ",
    "receive, read, or store your username or password.\n\n",
    
    "If you do not yet have a Copernicus Marine account, registration is free:\n",
    "  ", .COPERNICUS_REGISTER_URL, "\n\n",
    
    "Authentication documentation:\n",
    "  ", .COPERNICUS_TOOLBOX_CREDENTIALS_URL, "\n"
  )
  
  # -----------------------------------------------------------------------
  # RStudio
  #
  # External interactive password prompts do not reliably receive stdin from
  # the RStudio Console. Open the official Toolbox login in the RStudio
  # Terminal instead.
  # -----------------------------------------------------------------------
  
  rstudio_available <-
    requireNamespace("rstudioapi", quietly = TRUE) &&
    isTRUE(rstudioapi::isAvailable())
  
  if (rstudio_available) {
    
    # Create a small temporary R script that will run inside the RStudio
    # Terminal. Running the Toolbox from a real terminal allows its username
    # and password prompts to receive interactive input correctly.
    login_script <- tempfile(
      pattern = "copernicus_login_",
      fileext = ".R"
    )
    
    executable_for_script <- normalizePath(
      executable,
      winslash = "/",
      mustWork = TRUE
    )
    
    writeLines(
      c(
        sprintf(
          "executable <- %s",
          encodeString(executable_for_script, quote = "\"")
        ),
        
        "",
        
        "cat(",
        "  \"\\n\",",
        "  \"Copernicus Marine login\\n\",",
        "  \"------------------------\\n\",",
        "  \"The official Copernicus Marine Toolbox will now ask for your username and password.\\n\",",
        "  \"Your credentials are handled by the Toolbox and are not read or stored by this R package.\\n\\n\",",
        "  sep = \"\"",
        ")",
        
        "",
        
        "status <- system2(",
        "  executable,",
        "  args = c(",
        "    \"login\",",
        "    \"--log-level\",",
        "    \"ERROR\"",
        "  )",
        ")",
        
        "",
        
        "cat(\"\\n\\n\")",
        
        "",
        
        "if (identical(as.integer(status), 0L)) {",
        
        "  cat(",
        "    \"Copernicus Marine login completed successfully.\\n\\n\",",
        "    \"Your credentials have been configured by the official Copernicus Marine Toolbox.\\n\\n\",",
        "    \"Return to the R Console and run:\\n\\n\",",
        "    \"  copernicus_login()\\n\\n\",",
        "    \"again to verify that authentication is ready.\\n\",",
        "    sep = \"\"",
        "  )",
        
        "} else {",
        
        "  cat(",
        "    \"Copernicus Marine login did not complete successfully.\\n\\n\",",
        "    \"Return to the R Console and run `copernicus_login()` again to retry.\\n\",",
        "    sep = \"\"",
        "  )",
        
        "}",
        
        "",
        
        "quit(",
        "  save = \"no\",",
        "  status = as.integer(status)",
        ")"
      ),
      con = login_script
    )
    
    # Locate the Rscript executable belonging to the current R installation.
    rscript <- file.path(
      R.home("bin"),
      if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
    )
    
    if (!file.exists(rscript)) {
      stop(
        "The Rscript executable required to start the interactive Copernicus ",
        "Marine login could not be found.\n\n",
        "Expected location:\n  ",
        rscript,
        "\n\n",
        "You can still authenticate manually by opening a system terminal and ",
        "running:\n\n  ",
        shQuote(executable),
        " login\n\n",
        "After completing the login, return to R and run ",
        "`copernicus_login()` again.",
        call. = FALSE
      )
    }
    
    rscript_terminal <- normalizePath(
      rscript,
      winslash = "/",
      mustWork = TRUE
    )
    
    login_script_terminal <- normalizePath(
      login_script,
      winslash = "/",
      mustWork = TRUE
    )
    
    terminal_command <- paste(
      shQuote(rscript_terminal),
      shQuote(login_script_terminal)
    )
    
    message(
      "\nThe Copernicus Marine login will now open in the RStudio Terminal.\n\n",
      "Please complete the username and password prompts in the Terminal tab. ",
      "It may take a few seconds for the Toolbox to start, so wait for the ",
      "`Copernicus Marine username:` prompt before typing anything.\n\n",
      "The password may not display any characters while you type it; this is ",
      "normal for secure password entry.\n\n",
      "When the login finishes, the Terminal will tell you exactly how to ",
      "continue."
    )
    
    tryCatch(
      rstudioapi::terminalExecute(
        command = terminal_command,
        workingDir = getwd(),
        show = TRUE
      ),
      error = function(e) {
        stop(
          "The Copernicus Marine login could not be opened in the RStudio ",
          "Terminal.\n\n",
          "You can complete authentication manually by opening a system ",
          "terminal and running:\n\n  ",
          shQuote(executable),
          " login\n\n",
          "After completing the login, return to R and run ",
          "`copernicus_login()` again.\n\n",
          "Original error:\n  ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    
    return(invisible(FALSE))
  } else {
    
    # ---------------------------------------------------------------------
    # Interactive R session outside RStudio.
    #
    # A terminal-based R session can normally pass stdin directly to the
    # Toolbox. If the surrounding R environment cannot do so, provide the
    # exact command that the user can run in a system terminal.
    # ---------------------------------------------------------------------
    
    message(
      "\nThe official Copernicus Marine Toolbox will now start an interactive ",
      "login in the current terminal.\n\n",
      "Wait for the `Copernicus Marine username:` prompt, then enter your ",
      "username and password directly into the Toolbox."
    )
    
    status <- tryCatch(
      system2(
        executable,
        args = "login"
      ),
      error = function(e) {
        stop(
          "The interactive Copernicus Marine login could not be started from ",
          "the current R session.\n\n",
          "This can happen when the R environment does not provide an ",
          "interactive terminal to external programs.\n\n",
          "Open a system terminal and run:\n\n  ",
          shQuote(executable),
          " login\n\n",
          "After the Toolbox reports a successful login, return to R and run ",
          "`copernicus_login()` again.\n\n",
          "Original error:\n  ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    
    if (!identical(as.integer(status), 0L)) {
      stop(
        "Copernicus Marine authentication did not complete successfully.\n\n",
        "The official Toolbox returned exit status ",
        status,
        ".\n\n",
        "If no username or password prompt appeared, the current R environment ",
        "may not support interactive input for external programs.\n\n",
        "In that case, open a system terminal and run:\n\n  ",
        shQuote(executable),
        " login\n\n",
        "Possible causes also include invalid credentials, a cancelled login, ",
        "network connectivity problems, or a temporary authentication service ",
        "problem.\n\n",
        "For more information, see:\n  ",
        .COPERNICUS_TOOLBOX_CREDENTIALS_URL,
        call. = FALSE
      )
    }
    
    message(
      "\nCopernicus Marine authentication completed successfully.\n\n",
      "The official Toolbox has configured credentials that can be reused for ",
      "future Copernicus Marine requests.\n\n",
      "Credential discovery, validation, and storage are handled entirely by ",
      "the official Toolbox. This R package does not read, receive, or store ",
      "your username or password.\n"
    )
    
    invisible(TRUE)
  }
}



# Internal: return information about the standalone binary for this platform.
.copernicus_binary_info <- function(
    version = .COPERNICUS_TOOLBOX_VERSION
) {
  .copernicus_validate_string(version, "version")
  
  os <- Sys.info()[["sysname"]]
  machine <- tolower(Sys.info()[["machine"]])
  
  if (identical(.Platform$OS.type, "windows")) {
    asset <- "copernicusmarine.exe"
    
  } else if (identical(os, "Darwin")) {
    if (machine %in% c("arm64", "aarch64")) {
      asset <- "copernicusmarine_macos-arm64"
    } else if (machine %in% c("x86_64", "amd64")) {
      asset <- "copernicusmarine_macos-x86_64"
    } else {
      stop(
        "Unsupported macOS architecture for the Copernicus Marine Toolbox: ",
        machine,
        ".",
        call. = FALSE
      )
    }
    
  } else if (identical(os, "Linux")) {
    asset <- "copernicusmarine_linux"
    
  } else {
    stop(
      "Automatic Copernicus Marine Toolbox setup is not supported on this ",
      "operating system: ",
      if (is.null(os) || is.na(os)) "unknown" else os,
      ".",
      call. = FALSE
    )
  }
  
  list(
    version = version,
    asset = asset,
    url = paste0(
      .COPERNICUS_TOOLBOX_REPOSITORY,
      "/releases/download/v",
      version,
      "/",
      asset
    )
  )
}

#' Set up the Copernicus Marine Toolbox
#'
#' Download and configure a standalone version of the official Copernicus
#' Marine Toolbox for use by this package. The standalone executable does not
#' require a separate Python installation.
#'
#' The executable is downloaded from the official versioned Copernicus Marine
#' Toolbox release and stored in the user's platform-specific application data
#' directory. Nothing is downloaded when this package itself is installed or
#' loaded.
#'
#' @param version Character. Copernicus Marine Toolbox version to install.
#'   Defaults to the version tested by this package.
#' @param force Logical. If `TRUE`, reinstall the executable even when the
#'   requested version is already available.
#' @param quiet Logical. Suppress setup messages where possible.
#'
#' @returns Invisibly, the result of [copernicus_status()] for the installed
#'   executable.
#'
#' @examples
#' \dontrun{
#' copernicus_setup()
#' }
#'
#' @export
copernicus_setup <- function(
    version = .COPERNICUS_TOOLBOX_VERSION,
    force = FALSE,
    quiet = FALSE
) {
  .copernicus_validate_string(version, "version")
  
  if (!is.logical(force) || length(force) != 1 || is.na(force)) {
    stop("`force` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(quiet) || length(quiet) != 1 || is.na(quiet)) {
    stop("`quiet` must be TRUE or FALSE.", call. = FALSE)
  }
  
  info <- .copernicus_binary_info(version)
  install_dir <- .copernicus_toolbox_dir(version)
  
  executable_name <- if (.Platform$OS.type == "windows") {
    "copernicusmarine.exe"
  } else {
    "copernicusmarine"
  }
  
  dest <- file.path(install_dir, executable_name)
  
  if (file.exists(dest) && !force) {
    installed_version <- .copernicus_executable_version(dest)
    
    if (identical(installed_version, version)) {
      if (!quiet) {
        message(
          "Copernicus Marine Toolbox ", version,
          " is already set up."
        )
      }
      
      return(invisible(copernicus_status(dest)))
    }
    
    stop(
      "A Copernicus Marine Toolbox executable already exists at:\n  ",
      dest,
      "\nbut its version could not be confirmed as ", version,
      ". Use `force = TRUE` to replace it.",
      call. = FALSE
    )
  }
  
  if (!interactive()) {
    stop(
      "`copernicus_setup()` downloads an external executable and therefore ",
      "requires an interactive session. In non-interactive environments, ",
      "install the Copernicus Marine Toolbox separately and make it available ",
      "on the system PATH.",
      call. = FALSE
    )
  }
  
  message(
    "Using `copernicus_load(source = \"marine\", ...)` requires the official ",
    "Copernicus Marine Toolbox.\n\n",
    
    "The Toolbox is the software used to communicate with Copernicus Marine ",
    "and download the spatial, temporal, depth, and variable subsets requested ",
    "by `copernicus_load()`.\n\n",
    
    "To simplify the setup, this package can use the official standalone ",
    "Toolbox executable. This avoids requiring you to install or configure ",
    "Python or the Python `copernicusmarine` package yourself.\n\n",
    
    "The executable will be downloaded directly from the official Copernicus ",
    "Marine Toolbox release and stored in your user application-data directory. ",
    "It is not included in this R package.\n\n",
    
    "More information about the Copernicus Marine Toolbox:\n",
    "  ", .COPERNICUS_TOOLBOX_INFO_URL, "\n\n",
    
    "Installation details:\n",
    "  Version:     ", version, "\n",
    "  Source:      ", .COPERNICUS_TOOLBOX_REPOSITORY, "\n",
    "  Destination: ", dest, "\n\n",
    
    "Nothing will be downloaded until you confirm below."
  )
  
  ans <- readline(
    "Download and configure the Copernicus Marine Toolbox now? [y/N]: "
  )
  
  if (!tolower(trimws(ans)) %in% c("y", "yes")) {
    message("Setup cancelled.")
    return(invisible(FALSE))
  }
  
  dir.create(install_dir, recursive = TRUE, showWarnings = FALSE)
  
  tmp <- tempfile(
    pattern = "copernicusmarine_",
    tmpdir = install_dir
  )
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  
  if (!quiet) {
    message(
      "Downloading Copernicus Marine Toolbox ", version, "..."
    )
  }
  
  tryCatch(
    utils::download.file(
      url = info$url,
      destfile = tmp,
      mode = "wb",
      quiet = quiet
    ),
    error = function(e) {
      stop(
        "Could not download the Copernicus Marine Toolbox.\n",
        "Source: ", info$url, "\n",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  if (!file.exists(tmp) || is.na(file.info(tmp)$size) ||
      file.info(tmp)$size <= 0) {
    stop(
      "The Copernicus Marine Toolbox download did not produce a valid file.",
      call. = FALSE
    )
  }
  
  if (!file.copy(tmp, dest, overwrite = TRUE)) {
    stop(
      "Could not move the downloaded Copernicus Marine Toolbox to:\n  ",
      dest,
      call. = FALSE
    )
  }
  
  if (.Platform$OS.type != "windows") {
    Sys.chmod(dest, mode = "0755")
  }
  
  installed_version <- .copernicus_executable_version(dest)
  
  if (is.na(installed_version)) {
    unlink(dest, force = TRUE)
    
    stop(
      "The downloaded file could not be executed as a Copernicus Marine ",
      "Toolbox executable.",
      call. = FALSE
    )
  }
  
  if (!identical(installed_version, version)) {
    unlink(dest, force = TRUE)
    
    stop(
      "Copernicus Marine Toolbox version verification failed. Expected ",
      version, " but the downloaded executable reported ",
      installed_version, ".",
      call. = FALSE
    )
  }
  
  if (!quiet) {
    message(
      "Copernicus Marine Toolbox ", installed_version,
      " is ready to use."
    )
  }
  
  invisible(copernicus_status(dest))
}

# Internal: validate a single character argument.
.copernicus_validate_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    stop("`", name, "` must be a single non-empty character string.", call. = FALSE)
  }
  invisible(TRUE)
}

# Internal: validate geographical bounds.
.copernicus_validate_bbox <- function(xmin, xmax, ymin, ymax) {
  bbox <- list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  supplied <- !vapply(bbox, is.null, logical(1))

  if (!any(supplied)) return(invisible(TRUE))
  if (!all(supplied)) {
    stop("`xmin`, `xmax`, `ymin`, and `ymax` must be supplied together.", call. = FALSE)
  }

  vals <- unlist(bbox, use.names = FALSE)
  if (!is.numeric(vals) || any(!is.finite(vals)) || any(lengths(bbox) != 1)) {
    stop("Geographical bounds must be single finite numeric values.", call. = FALSE)
  }

  if (xmin >= xmax) {
    stop("Longitude bounds must satisfy `xmin < xmax`.", call. = FALSE)
  }

  if (ymin < -90 || ymax > 90 || ymin >= ymax) {
    stop("Latitude bounds must satisfy -90 <= ymin < ymax <= 90.", call. = FALSE)
  }

  invisible(TRUE)
}

# Internal: validate depth range.
.copernicus_validate_depth <- function(depth_min, depth_max) {
  if (is.null(depth_min) && is.null(depth_max)) return(invisible(TRUE))

  depth <- list(depth_min = depth_min, depth_max = depth_max)
  supplied <- !vapply(depth, is.null, logical(1))
  vals <- unlist(depth[supplied], use.names = FALSE)

  if (!is.numeric(vals) || any(!is.finite(vals)) ||
      any(lengths(depth[supplied]) != 1)) {
    stop("Depth bounds must be single finite numeric values.", call. = FALSE)
  }

  if (any(vals < 0)) {
    stop("Depth bounds must be >= 0 metres (positive downward).", call. = FALSE)
  }

  if (all(supplied) && depth_min > depth_max) {
    stop("`depth_min` must be <= `depth_max`.", call. = FALSE)
  }

  invisible(TRUE)
}

# Internal: convert Date/POSIXt/character input to UTC ISO datetime.
.copernicus_datetime <- function(x, name = "datetime") {
  if (is.null(x)) return(NULL)
  
  if (inherits(x, "Date")) {
    x <- as.POSIXct(x, tz = "UTC")
  } else if (inherits(x, "POSIXt")) {
    x <- as.POSIXct(x, tz = "UTC")
  } else if (is.character(x) && length(x) == 1 && !is.na(x)) {
    if (grepl("^\\d{4}-\\d{2}-\\d{2}$", x)) x <- paste0(x, " 00:00:00")
    
    x <- tryCatch(
      suppressWarnings(as.POSIXct(
        gsub("T", " ", sub("Z$", "", x), fixed = TRUE),
        tz = "UTC",
        tryFormats = c("%Y-%m-%d %H:%M:%OS", "%Y-%m-%d")
      )),
      error = function(e) as.POSIXct(NA, tz = "UTC")
    )
  } else {
    stop(
      "`", name, "` must be a Date, POSIXt, or single character datetime.",
      call. = FALSE
    )
  }
  
  if (length(x) != 1 || is.na(x)) {
    stop("Could not parse `", name, "`.", call. = FALSE)
  }
  
  format(x, "%Y-%m-%dT%H:%M:%S", tz = "UTC")
}

# Internal: prepare and validate temporal bounds.
.copernicus_time_range <- function(start_datetime, end_datetime) {
  start <- .copernicus_datetime(start_datetime, "start_datetime")
  end <- .copernicus_datetime(end_datetime, "end_datetime")

  if (is.null(start) && !is.null(end)) {
    stop("`start_datetime` is required when `end_datetime` is supplied.", call. = FALSE)
  }

  if (!is.null(start) && is.null(end)) end <- start

  if (!is.null(start) &&
      as.POSIXct(end, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC") <
      as.POSIXct(start, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")) {
    stop("`end_datetime` must be equal to or later than `start_datetime`.",
         call. = FALSE)
  }

  list(start = start, end = end)
}

# Internal: extract paths from Copernicus Marine response objects.
.copernicus_marine_paths <- function(result) {
  get_one <- function(x) {
    path <- tryCatch(as.character(x$file_path), error = function(e) character())
    if (!length(path)) return(character())

    file_names <- tryCatch(
      unlist(x$file_names, use.names = FALSE),
      error = function(e) character()
    )

    if (length(path) == 1 && dir.exists(path) && length(file_names)) {
      candidates <- file.path(path, file_names)
      if (all(file.exists(candidates))) return(candidates)
    }

    path
  }

  if (is.list(result) && is.null(result$file_path)) {
    paths <- unlist(lapply(result, get_one), use.names = FALSE)
  } else {
    paths <- get_one(result)
  }

  paths <- unique(as.character(paths))
  if (!length(paths)) {
    stop("Copernicus Marine completed the request but returned no output path.",
         call. = FALSE)
  }

  paths
}

#' Load environmental data from Copernicus services
#'
#' Acquire environmental data from Copernicus Marine or ECMWF-operated
#' Copernicus Data Stores and save the resulting files locally. Supported
#' sources are Copernicus Marine (`"marine"`), Climate Data Store (`"cds"`),
#' and Atmosphere Data Store (`"ads"`).
#'
#' Copernicus Marine requests can use either the official standalone
#' Copernicus Marine Toolbox or the Python `copernicusmarine` package through
#' `reticulate`. The standalone Toolbox is the recommended backend because it
#' does not require users to install or configure Python.
#'
#' With `backend = "auto"`, the standalone Toolbox is used when available.
#' The Python backend is used only when the standalone executable is not
#' available. Backend selection occurs before the request starts; a failed
#' request is never automatically retried using another backend.
#'
#' For Copernicus Marine, each output NetCDF file contains exactly one
#' requested environmental variable. When several variables are supplied,
#' they are automatically separated into different files. `split_by` therefore
#' controls temporal subdivision only.
#'
#' Downloaded Copernicus Marine files can optionally be arranged into a
#' directory hierarchy with `organize_by`. For example,
#'
#' `c("service", "dataset", "variable", "year", "month", "day")`
#'
#' stores each variable independently below the service and dataset, followed
#' by calendar year, month, and day.
#'
#' Temporal organization cannot be finer than `split_by`. For example,
#' `split_by = "month"` cannot be combined with
#' `organize_by = c("year", "month", "day")`, because one monthly output file
#' can contain multiple days. Such incompatible combinations are detected
#' before a download starts.
#'
#' Week splitting and organization use ISO weeks. ISO weeks start on Monday
#' and week 1 is the week containing 4 January. The associated year is the ISO
#' week-year rather than necessarily the ordinary calendar year.
#'
#' Season splitting uses meteorological rather than astronomical seasons.
#' Meteorological seasons are fixed three-month periods beginning on the first
#' day of a month:
#'
#' * DJF: December-February
#' * MAM: March-May
#' * JJA: June-August
#' * SON: September-November
#'
#' These temporal boundaries are identical in both hemispheres.
#' `hemisphere` determines only the human-readable season name. Thus DJF is
#' `"winter_DJF"` in the Northern Hemisphere and `"summer_DJF"` in the
#' Southern Hemisphere. DJF is assigned to the year containing January and
#' February; for example, December 2025 belongs to DJF 2026.
#'
#' The official Copernicus Marine Toolbox natively supports splitting by year,
#' month, day, and hour. Week and meteorological-season splitting are managed
#' by this package by issuing one appropriately bounded Toolbox request per
#' temporal period.
#'
#' CDS and ADS requests use `ecmwfr`. File organization through `organize_by`
#' is currently implemented for Copernicus Marine only.
#'
#' Authentication is managed through the official provider clients. For
#' Copernicus Marine, run [copernicus_login()] once before the first request.
#' The credentials created by the official Copernicus Marine Toolbox are
#' shared by both the standalone and Python backends, so the same login can be
#' used with `backend = "standalone"` or `backend = "python"`.
#'
#' If the Python backend does not find valid credentials, it asks the user to
#' configure authentication with [copernicus_login()] rather than attempting
#' an interactive password prompt through `reticulate`, because such prompts
#' are not reliably supported in R sessions.
#'
#' @param source Character. Copernicus service: `"marine"`, `"cds"`, or
#'   `"ads"`.
#' @param backend Character. Backend used for Copernicus Marine requests:
#'   `"auto"`, `"standalone"`, or `"python"`.
#' @param dataset_id Character. Provider dataset identifier.
#' @param variables Character vector containing one or more environmental
#'   variables to download. For Copernicus Marine this must be supplied
#'   explicitly. Each requested variable is stored in separate NetCDF files.
#' @param start_datetime Optional start datetime.
#' @param end_datetime Optional end datetime.
#' @param xmin,xmax,ymin,ymax Optional geographic bounding box.
#' @param depth_min,depth_max Optional depth range in metres, positive
#'   downward.
#' @param request Named list of provider- or dataset-specific request
#'   arguments.
#' @param output_dir Character. Root destination directory.
#' @param filename Optional base output filename.
#' @param split_by Optional character defining temporal file subdivision for
#'   Copernicus Marine. Supported values are `"year"`, `"season"`, `"month"`,
#'   `"week"`, `"day"`, and `"hour"`. Variables are always split
#'   automatically.
#' @param organize_by Optional character vector defining the directory
#'   hierarchy. Supported components are `"service"`, `"dataset"`,
#'   `"variable"`, `"year"`, `"season"`, `"month"`, `"week"`, `"day"`,
#'   and `"hour"`.
#' @param hemisphere Character. `"north"` or `"south"`. Used only to assign
#'   human-readable names to meteorological seasons.
#' @param concurrent_processes Optional integer >= 1.
#' @param compression Optional integer from 0 to 9.
#' @param dataset_version Optional Copernicus Marine dataset version.
#' @param dataset_part Optional Copernicus Marine dataset part.
#' @param force Logical. Replace existing output files.
#' @param quiet Logical. Suppress progress messages where possible.
#'
#' @returns Character vector containing paths to downloaded files.
#'
#' @examples
#' \dontrun{
#' files <- copernicus_load(
#'   source = "marine",
#'   dataset_id = "cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m",
#'   variables = "thetao",
#'   start_datetime = "2026-07-01",
#'   end_datetime = "2026-07-31",
#'   xmin = -6,
#'   xmax = 10,
#'   ymin = 35,
#'   ymax = 45,
#'   depth_min = 0,
#'   depth_max = 500,
#'   split_by = "day",
#'   organize_by = c(
#'     "service",
#'     "dataset",
#'     "variable",
#'     "year",
#'     "month",
#'     "day"
#'   )
#' )
#' }
#'
#' @export
copernicus_load <- function(
    source = c("marine", "cds", "ads"),
    backend = c("auto", "standalone", "python"),
    dataset_id,
    variables = NULL,
    start_datetime = NULL,
    end_datetime = start_datetime,
    xmin = NULL,
    xmax = NULL,
    ymin = NULL,
    ymax = NULL,
    depth_min = NULL,
    depth_max = NULL,
    request = list(),
    output_dir = NULL,
    filename = NULL,
    split_by = NULL,
    organize_by = NULL,
    hemisphere = c("north", "south"),
    concurrent_processes = NULL,
    compression = NULL,
    dataset_version = NULL,
    dataset_part = NULL,
    force = FALSE,
    quiet = FALSE
) {
  
  source <- match.arg(source)
  hemisphere <- match.arg(hemisphere)
  
  .copernicus_validate_string(
    dataset_id,
    "dataset_id"
  )
  
  if (!is.null(variables) &&
      (!is.character(variables) ||
       !length(variables) ||
       anyNA(variables) ||
       any(!nzchar(variables)))) {
    
    stop(
      "`variables` must be a non-empty character vector.",
      call. = FALSE
    )
  }
  
  if (!is.null(variables) &&
      anyDuplicated(variables)) {
    
    stop(
      "`variables` must not contain duplicated variable names.",
      call. = FALSE
    )
  }
  
  if (!is.list(request) ||
      (length(request) &&
       (is.null(names(request)) ||
        any(!nzchar(names(request)))))) {
    
    stop(
      "`request` must be a named list.",
      call. = FALSE
    )
  }
  
  if (!is.logical(force) ||
      length(force) != 1 ||
      is.na(force)) {
    
    stop(
      "`force` must be TRUE or FALSE.",
      call. = FALSE
    )
  }
  
  if (!is.logical(quiet) ||
      length(quiet) != 1 ||
      is.na(quiet)) {
    
    stop(
      "`quiet` must be TRUE or FALSE.",
      call. = FALSE
    )
  }
  
  if (!is.null(compression) &&
      (!is.numeric(compression) ||
       length(compression) != 1L ||
       !is.finite(compression) ||
       compression < 0 ||
       compression > 9 ||
       compression != as.integer(compression))) {
    stop(
      "`compression` must be a single integer from 0 to 9.",
      call. = FALSE
    )
  }
  .copernicus_validate_bbox(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax
  )
  
  .copernicus_validate_depth(
    depth_min = depth_min,
    depth_max = depth_max
  )
  
  time <- .copernicus_time_range(
    start_datetime,
    end_datetime
  )
  
  if (identical(source, "marine")) {
    
    backend <- match.arg(backend)
    
    if (is.null(variables)) {
      stop(
        "`variables` must be supplied explicitly for Copernicus Marine ",
        "requests.\n\n",
        "This package deliberately stores one environmental variable per ",
        "NetCDF file so that variables with different grids, depths, temporal ",
        "resolution, or metadata are never unintentionally mixed.\n\n",
        "Specify one or more variables, for example:\n\n",
        "  variables = c(\"thetao\", \"so\")",
        call. = FALSE
      )
    }
    
    .copernicus_validate_file_structure(
      split_by = split_by,
      organize_by = organize_by,
      hemisphere = hemisphere,
      concurrent_processes = concurrent_processes,
      variables = variables
    )
    
  } else {
    
    if (!is.null(organize_by)) {
      stop(
        "`organize_by` is currently implemented for ",
        "`source = \"marine\"` only.",
        call. = FALSE
      )
    }
    
    if (!is.null(split_by) &&
        split_by %in% c(
          "season",
          "week"
        )) {
      
      stop(
        "`split_by = \"",
        split_by,
        "\"` is currently implemented for Copernicus Marine only.",
        call. = FALSE
      )
    }
  }
  
  using_default_cache <- is.null(output_dir)
  
  if (using_default_cache) {
    
    cache_root <- copernicus_cache_dir()
    
    .copernicus_cache_consent(
      cache_root
    )
    
    output_dir <- file.path(
      cache_root,
      source
    )
  }
  
  if (!dir.exists(output_dir)) {
    
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  if (identical(source, "marine")) {
    
    marine_backend <- .copernicus_select_marine_backend(
      backend
    )
    
    if (!quiet) {
      message(
        "Using Copernicus Marine backend: ",
        marine_backend,
        "\n"
      )
    }
    
    if (identical(
      marine_backend,
      "standalone"
    )) {
      
      return(
        .copernicus_load_marine_standalone(
          dataset_id = dataset_id,
          variables = variables,
          start_datetime = time$start,
          end_datetime = time$end,
          xmin = xmin,
          xmax = xmax,
          ymin = ymin,
          ymax = ymax,
          depth_min = depth_min,
          depth_max = depth_max,
          request = request,
          output_dir = output_dir,
          filename = filename,
          split_by = split_by,
          organize_by = organize_by,
          hemisphere = hemisphere,
          concurrent_processes = concurrent_processes,
          compression = compression,
          dataset_version = dataset_version,
          dataset_part = dataset_part,
          force = force,
          quiet = quiet
        )
      )
    }
    
    if (!is.null(split_by) &&
        split_by %in% c(
          "week",
          "season"
        )) {
      
      stop(
        "`split_by = \"",
        split_by,
        "\"` requires the standalone Copernicus Marine backend.",
        call. = FALSE
      )
    }
    
    if (!is.null(organize_by)) {
      stop(
        "`organize_by` currently requires the standalone Copernicus Marine ",
        "backend.",
        call. = FALSE
      )
    }
    
    python_paths <- character()
    
    for (i in seq_along(variables)) {
      
      variable <- variables[[i]]
      
      variable_filename <- filename
      
      if (!is.null(filename) &&
          length(variables) > 1L) {
        
        stem <- tools::file_path_sans_ext(
          filename
        )
        
        ext <- tools::file_ext(
          filename
        )
        
        variable_filename <- paste0(
          stem,
          "_",
          .copernicus_safe_component(
            variable
          ),
          if (nzchar(ext)) {
            paste0(
              ".",
              ext
            )
          } else {
            ""
          }
        )
      }
      
      this_path <- .copernicus_load_marine_python(
        dataset_id = dataset_id,
        variables = variable,
        start_datetime = time$start,
        end_datetime = time$end,
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        depth_min = depth_min,
        depth_max = depth_max,
        request = request,
        output_dir = output_dir,
        filename = variable_filename,
        split_by = split_by,
        concurrent_processes = concurrent_processes,
        compression = compression,
        dataset_version = dataset_version,
        dataset_part = dataset_part,
        force = force,
        quiet = quiet
      )
      
      python_paths <- c(
        python_paths,
        this_path
      )
    }
    
    return(
      unique(
        python_paths
      )
    )
  }
  
  .copernicus_load_ecmwf(
    source = source,
    dataset_id = dataset_id,
    variables = variables,
    start_datetime = time$start,
    end_datetime = time$end,
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    depth_min = depth_min,
    depth_max = depth_max,
    request = request,
    output_dir = output_dir,
    filename = filename,
    split_by = split_by,
    concurrent_processes = concurrent_processes,
    compression = compression,
    dataset_version = dataset_version,
    dataset_part = dataset_part,
    force = force,
    quiet = quiet
  )
}

# Internal: select the backend used for Copernicus Marine requests.
#
# The standalone Copernicus Marine Toolbox is preferred because it is the
# officially supported command-line interface and does not require users to
# install or configure Python.
#
# `backend = "auto"` selects the standalone Toolbox when it is available.
# The Python backend is used only when the standalone executable is not
# available and reticulate is installed.
#
# Backend selection happens once before a request starts. A failed request is
# never automatically retried with another backend because doing so could hide
# genuine request, authentication, dataset, or connectivity errors.
.copernicus_select_marine_backend <- function(
    backend = c("auto", "standalone", "python")
) {
  
  backend <- match.arg(backend)
  
  standalone <- .copernicus_find_executable()
  
  standalone_available <-
    !is.null(standalone) &&
    file.exists(standalone)
  
  python_available <-
    requireNamespace("reticulate", quietly = TRUE)
  
  # Explicit standalone request.
  if (identical(backend, "standalone")) {
    
    if (!standalone_available) {
      stop(
        "The standalone Copernicus Marine Toolbox was explicitly requested, ",
        "but no usable executable could be found.\n\n",
        "Run `copernicus_setup()` to install the supported standalone Toolbox, ",
        "or make an existing `copernicusmarine` executable available on the ",
        "system PATH.\n\n",
        "You can inspect the current Toolbox configuration with:\n\n",
        "  copernicus_status()\n\n",
        "For more information, see:\n  ",
        .COPERNICUS_TOOLBOX_INFO_URL,
        call. = FALSE
      )
    }
    
    return("standalone")
  }
  
  # Explicit Python request.
  if (identical(backend, "python")) {
    
    if (!python_available) {
      stop(
        "The Python Copernicus Marine backend was explicitly requested, but ",
        "the R package `reticulate` is not installed.\n\n",
        "The recommended backend is the standalone Copernicus Marine Toolbox, ",
        "which does not require users to install or configure Python. Run:\n\n",
        "  copernicus_setup()\n\n",
        "Alternatively, install `reticulate` if you specifically need the ",
        "Python backend:\n\n",
        "  install.packages(\"reticulate\")",
        call. = FALSE
      )
    }
    
    return("python")
  }
  
  # Automatic selection:
  # prefer the supported standalone Toolbox.
  if (standalone_available) {
    return("standalone")
  }
  
  # Retain the existing Python implementation as a fallback when the
  # standalone Toolbox is unavailable.
  if (python_available) {
    return("python")
  }
  
  stop(
    "No usable Copernicus Marine backend is currently available.\n\n",
    "The recommended option is the official standalone Copernicus Marine ",
    "Toolbox. It can be configured with:\n\n",
    "  copernicus_setup()\n\n",
    "This avoids requiring a separate Python installation or Python ",
    "configuration.\n\n",
    "A legacy Python backend is also available when the R package ",
    "`reticulate` is installed.\n\n",
    "Check the standalone Toolbox configuration with:\n\n",
    "  copernicus_status()\n\n",
    "For more information, see:\n  ",
    .COPERNICUS_TOOLBOX_INFO_URL,
    call. = FALSE
  )
}


# Internal: prepare the Python dependency used by Copernicus Marine.
.copernicus_marine_module <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      "Package 'reticulate' is required for `source = \"marine\"`. ",
      "Install it with install.packages('reticulate').",
      call. = FALSE
    )
  }

  # If Python has not yet been initialized, prefer reticulate's managed
  # environment unless the user has explicitly selected another Python through
  # RETICULATE_PYTHON, RETICULATE_PYTHON_ENV, use_python(), or use_virtualenv().
  # This avoids an unrelated pre-existing `r-reticulate` environment taking
  # precedence over the requirements declared below.
  if (!reticulate::py_available(initialize = FALSE)) {
    Sys.setenv(RETICULATE_USE_MANAGED_VENV = "yes")
  }

  reticulate::py_require(
    packages = "copernicusmarine>=2.2,<3",
    python_version = ">=3.10,<3.14"
  )

  tryCatch(
    reticulate::import("copernicusmarine", convert = TRUE),
    error = function(e) {
      py_initialized <- reticulate::py_available(initialize = FALSE)

      stop(
        "Could not import the Python package 'copernicusmarine'. ",
        if (py_initialized) {
          paste0(
            "Python is already initialized in this R session. Restart R and ",
            "call `copernicus_load()` before initializing another Python ",
            "environment, or explicitly select an environment containing ",
            "copernicusmarine with `reticulate::use_virtualenv()` immediately ",
            "after restarting R. "
          )
        } else {
          paste0(
            "reticulate could not provision a compatible managed Python ",
            "environment. "
          )
        },
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

# Internal: validate Copernicus Marine credentials for the Python backend.
.copernicus_marine_auth <- function(cm) {
  valid <- tryCatch(
    isTRUE(cm$login(check_credentials_valid = TRUE)),
    error = function(e) FALSE
  )
  
  if (valid) {
    return(invisible(TRUE))
  }
  
  stop(
    "Valid Copernicus Marine credentials were not found.\n\n",
    "Configure authentication first with `copernicus_login()`, or run ",
    "`copernicusmarine.login()` directly in a Python terminal, then retry ",
    "`copernicus_load()`.\n\n",
    "Interactive Python login is not attempted through `reticulate` because ",
    "password prompts are not reliably supported in R sessions.",
    call. = FALSE
  )
}

# Internal: Copernicus Marine Python backend.
#
# This backend is retained as a fallback for users who already work with the
# Python `copernicusmarine` package through `reticulate`. The standalone
# Toolbox remains the recommended backend.
#
# `copernicus_load()` calls this function once per requested environmental
# variable, so each output NetCDF file contains exactly one variable.
#
# Week and season splitting, and directory organization through `organize_by`,
# are handled only by the standalone backend and are rejected before this
# function is called.
.copernicus_load_marine_python <- function(
    dataset_id,
    variables,
    start_datetime,
    end_datetime,
    xmin,
    xmax,
    ymin,
    ymax,
    depth_min,
    depth_max,
    request,
    output_dir,
    filename,
    split_by,
    concurrent_processes,
    compression,
    dataset_version,
    dataset_part,
    force,
    quiet
) {
  
  if (!is.character(variables) ||
      length(variables) != 1L ||
      is.na(variables) ||
      !nzchar(variables)) {
    
    stop(
      "Internal error: the Python Copernicus Marine backend must receive ",
      "exactly one environmental variable per request.",
      call. = FALSE
    )
  }
  
  if (!is.null(split_by)) {
    
    allowed_split <- c(
      "year",
      "month",
      "day",
      "hour"
    )
    
    if (!split_by %in% allowed_split) {
      stop(
        "The Python Copernicus Marine backend supports temporal splitting by:\n  ",
        paste(
          allowed_split,
          collapse = ", "
        ),
        "\n\n",
        "Week and season splitting require the standalone backend.",
        call. = FALSE
      )
    }
  }
  
  if (!is.null(concurrent_processes)) {
    
    if (is.null(split_by)) {
      stop(
        "`concurrent_processes` requires `split_by` when using the Python ",
        "Copernicus Marine backend.",
        call. = FALSE
      )
    }
    
    if (!is.numeric(concurrent_processes) ||
        length(concurrent_processes) != 1L ||
        !is.finite(concurrent_processes) ||
        concurrent_processes < 1 ||
        concurrent_processes != as.integer(concurrent_processes)) {
      
      stop(
        "`concurrent_processes` must be a single integer >= 1.",
        call. = FALSE
      )
    }
  }
  
  if (!is.null(compression) &&
      (!is.numeric(compression) ||
       length(compression) != 1L ||
       !is.finite(compression) ||
       compression < 0 ||
       compression > 9 ||
       compression != as.integer(compression))) {
    
    stop(
      "`compression` must be a single integer from 0 to 9.",
      call. = FALSE
    )
  }
  
  cm <- .copernicus_marine_module()
  
  .copernicus_marine_auth(
    cm
  )
  
  common <- list(
    dataset_id = dataset_id,
    output_directory = output_dir,
    overwrite = force,
    skip_existing = !force,
    disable_progress_bar = quiet,
    variables = as.list(
      variables
    )
  )
  
  if (!is.null(start_datetime)) {
    common$start_datetime <- start_datetime
  }
  
  if (!is.null(end_datetime)) {
    common$end_datetime <- end_datetime
  }
  
  if (!is.null(xmin)) {
    common$minimum_longitude <- xmin
  }
  
  if (!is.null(xmax)) {
    common$maximum_longitude <- xmax
  }
  
  if (!is.null(ymin)) {
    common$minimum_latitude <- ymin
  }
  
  if (!is.null(ymax)) {
    common$maximum_latitude <- ymax
  }
  
  if (!is.null(depth_min)) {
    common$minimum_depth <- depth_min
  }
  
  if (!is.null(depth_max)) {
    common$maximum_depth <- depth_max
  }
  
  if (!is.null(filename)) {
    common$output_filename <- filename
  }
  
  if (!is.null(compression)) {
    common$netcdf_compression_level <- as.integer(
      compression
    )
  }
  
  if (!is.null(dataset_version)) {
    common$dataset_version <- dataset_version
  }
  
  if (!is.null(dataset_part)) {
    common$dataset_part <- dataset_part
  }
  
  # Explicit function arguments take precedence over duplicate entries in
  # `request`.
  args <- utils::modifyList(
    request,
    common
  )
  
  if (!quiet) {
    message(
      "Requesting Copernicus Marine dataset ",
      dataset_id,
      " for variable ",
      variables,
      if (is.null(split_by)) {
        "..."
      } else {
        paste0(
          ", split by ",
          split_by,
          "..."
        )
      }
    )
  }
  
  result <- tryCatch(
    {
      if (is.null(split_by)) {
        
        do.call(
          cm$subset,
          args
        )
        
      } else {
        
        split_args <- args
        
        split_args$on_time <- split_by
        
        if (!is.null(concurrent_processes)) {
          split_args$concurrent_processes <- as.integer(
            concurrent_processes
          )
        }
        
        do.call(
          cm$subset_split_on,
          split_args
        )
      }
    },
    error = function(e) {
      stop(
        "Copernicus Marine Python request failed.\n\n",
        "Dataset:\n  ",
        dataset_id,
        "\n\n",
        "Variable:\n  ",
        variables,
        "\n\n",
        "Original error:\n  ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  paths <- .copernicus_marine_paths(
    result
  )
  
  missing <-
    !file.exists(paths) &
    !dir.exists(paths)
  
  if (any(missing)) {
    stop(
      "Copernicus Marine returned output path(s) that do not exist:\n  ",
      paste(
        paths[missing],
        collapse = "\n  "
      ),
      call. = FALSE
    )
  }
  
  normalizePath(
    paths,
    winslash = "/",
    mustWork = TRUE
  )
}

# Internal: validate Copernicus file-splitting and organization options.
#
# Copernicus downloads follow one important package-level rule:
# each output NetCDF file contains one requested environmental variable.
#
# `split_by` therefore controls only temporal splitting. Separation by
# variable is automatic when more than one variable is requested.
#
# Temporal organization must never be finer than the output files created by
# `split_by`. For example, monthly files cannot be organized into daily
# folders because one monthly file may contain several days.
.copernicus_validate_file_structure <- function(
    split_by = NULL,
    organize_by = NULL,
    hemisphere = "north",
    concurrent_processes = NULL,
    variables = NULL
) {
  
  allowed_split <- c(
    "year",
    "season",
    "month",
    "week",
    "day",
    "hour"
  )
  
  allowed_organize <- c(
    "service",
    "dataset",
    "variable",
    "year",
    "season",
    "month",
    "week",
    "day",
    "hour"
  )
  
  # -----------------------------------------------------------------------
  # split_by
  # -----------------------------------------------------------------------
  
  if (!is.null(split_by)) {
    
    .copernicus_validate_string(
      split_by,
      "split_by"
    )
    
    if (!split_by %in% allowed_split) {
      stop(
        "`split_by` must be one of:\n  ",
        paste(
          allowed_split,
          collapse = ", "
        ),
        "\n\n",
        "Variables do not need to be included in `split_by`. ",
        "When multiple variables are requested, this package automatically ",
        "creates separate NetCDF files for each variable.",
        call. = FALSE
      )
    }
  }
  
  # -----------------------------------------------------------------------
  # organize_by
  # -----------------------------------------------------------------------
  
  if (!is.null(organize_by)) {
    
    if (!is.character(organize_by) ||
        !length(organize_by) ||
        anyNA(organize_by) ||
        any(!nzchar(organize_by))) {
      
      stop(
        "`organize_by` must be NULL or a non-empty character vector.",
        call. = FALSE
      )
    }
    
    unknown <- setdiff(
      organize_by,
      allowed_organize
    )
    
    if (length(unknown)) {
      
      stop(
        "Unsupported `organize_by` value",
        if (length(unknown) > 1) "s" else "",
        ":\n  ",
        paste(
          unknown,
          collapse = ", "
        ),
        "\n\nAllowed values are:\n  ",
        paste(
          allowed_organize,
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    
    if (anyDuplicated(organize_by)) {
      stop(
        "`organize_by` must not contain duplicated levels.",
        call. = FALSE
      )
    }
    
    if ("variable" %in% organize_by &&
        is.null(variables)) {
      
      stop(
        "`organize_by` includes \"variable\", but `variables` was not supplied.\n\n",
        "Specify the environmental variable(s) explicitly, for example:\n\n",
        "  variables = c(\"thetao\", \"so\")\n\n",
        "Each requested variable will then be downloaded into separate ",
        "NetCDF files and can be stored in its own folder.",
        call. = FALSE
      )
    }
    
    # ---------------------------------------------------------------------
    # Temporal organization branches
    #
    # Calendar/season hierarchy:
    #   year -> season -> month -> day -> hour
    #
    # ISO-week hierarchy:
    #   year -> week -> day -> hour
    #
    # Month/season and ISO week are intentionally not mixed because an ISO
    # week can cross month and season boundaries.
    # ---------------------------------------------------------------------
    
    has_week <- "week" %in% organize_by
    has_month <- "month" %in% organize_by
    has_season <- "season" %in% organize_by
    
    if (has_week && (has_month || has_season)) {
      
      stop(
        "`organize_by` cannot combine \"week\" with \"month\" or \"season\".\n\n",
        "ISO weeks can cross month and meteorological-season boundaries, so ",
        "these directory structures are not unambiguously nested.\n\n",
        "Use either a calendar/season hierarchy, for example:\n\n",
        "  c(\"year\", \"season\", \"month\", \"day\")\n\n",
        "or an ISO-week hierarchy:\n\n",
        "  c(\"year\", \"week\", \"day\")",
        call. = FALSE
      )
    }
    
    temporal_levels <- organize_by[
      organize_by %in% c(
        "year",
        "season",
        "month",
        "week",
        "day",
        "hour"
      )
    ]
    
    if (length(temporal_levels)) {
      
      expected_order <- if (has_week) {
        c(
          "year",
          "week",
          "day",
          "hour"
        )
      } else {
        c(
          "year",
          "season",
          "month",
          "day",
          "hour"
        )
      }
      
      expected_temporal <- expected_order[
        expected_order %in% temporal_levels
      ]
      
      if (!identical(
        temporal_levels,
        expected_temporal
      )) {
        
        stop(
          "Temporal levels in `organize_by` are not in a valid hierarchical ",
          "order.\n\n",
          if (has_week) {
            paste0(
              "ISO-week organization must follow:\n\n",
              "  year -> week -> day -> hour"
            )
          } else {
            paste0(
              "Calendar/season organization must follow:\n\n",
              "  year -> season -> month -> day -> hour"
            )
          },
          "\n\nNot every level needs to be included, but their relative order ",
          "must be preserved.",
          call. = FALSE
        )
      }
    }
    
    # ---------------------------------------------------------------------
    # Check that organize_by is not finer than split_by.
    # ---------------------------------------------------------------------
    
    if (length(temporal_levels)) {
      
      if (is.null(split_by)) {
        
        stop(
          "`organize_by` contains temporal folders, but `split_by` is NULL.\n\n",
          "Without temporal splitting, one output file may contain multiple ",
          "years, months, weeks, days, or hours.\n\n",
          "Choose an appropriate temporal `split_by`, or remove temporal ",
          "levels from `organize_by`.",
          call. = FALSE
        )
      }
      
      allowed_by_split <- switch(
        split_by,
        
        year = c(
          "year"
        ),
        
        season = c(
          "year",
          "season"
        ),
        
        month = c(
          "year",
          "season",
          "month"
        ),
        
        week = c(
          "year",
          "week"
        ),
        
        day = c(
          "year",
          "season",
          "month",
          "week",
          "day"
        ),
        
        hour = c(
          "year",
          "season",
          "month",
          "week",
          "day",
          "hour"
        )
      )
      
      incompatible <- setdiff(
        temporal_levels,
        allowed_by_split
      )
      
      if (length(incompatible)) {
        
        finest <- incompatible[
          length(incompatible)
        ]
        
        stop(
          "`organize_by = \"",
          finest,
          "\"` is finer than or incompatible with ",
          "`split_by = \"",
          split_by,
          "\"`.\n\n",
          "A ",
          split_by,
          "-level output file can contain data from multiple ",
          finest,
          " units, so it cannot be placed unambiguously into one ",
          finest,
          " folder.\n\n",
          "Use a finer `split_by`, or remove \"",
          finest,
          "\" from `organize_by`.",
          call. = FALSE
        )
      }
    }
  }
  
  # -----------------------------------------------------------------------
  # hemisphere
  # -----------------------------------------------------------------------
  
  if (!is.character(hemisphere) ||
      length(hemisphere) != 1 ||
      is.na(hemisphere)) {
    
    stop(
      "`hemisphere` must be either \"north\" or \"south\".",
      call. = FALSE
    )
  }
  
  hemisphere <- match.arg(
    hemisphere,
    c(
      "north",
      "south"
    )
  )
  
  # -----------------------------------------------------------------------
  # concurrent_processes
  # -----------------------------------------------------------------------
  
  if (!is.null(concurrent_processes)) {
    
    if (!is.numeric(concurrent_processes) ||
        length(concurrent_processes) != 1 ||
        is.na(concurrent_processes) ||
        concurrent_processes < 1 ||
        concurrent_processes != as.integer(concurrent_processes)) {
      
      stop(
        "`concurrent_processes` must be a single integer >= 1.",
        call. = FALSE
      )
    }
    
    multiple_variables <-
      !is.null(variables) &&
      length(variables) > 1
    
    if (is.null(split_by) &&
        !multiple_variables) {
      
      stop(
        "`concurrent_processes` has no effect because the request creates ",
        "only one output file.\n\n",
        "Use it when requesting multiple variables or when `split_by` is ",
        "supplied.",
        call. = FALSE
      )
    }
    
    if (!is.null(split_by) &&
        split_by %in% c(
          "week",
          "season"
        ) &&
        !multiple_variables &&
        concurrent_processes > 1) {
      
      stop(
        "`concurrent_processes > 1` cannot currently be used with ",
        "`split_by = \"",
        split_by,
        "\"` when only one variable is requested.\n\n",
        "Week and season splitting are managed sequentially by this R package ",
        "because the official Toolbox does not provide native week or season ",
        "splitting.",
        call. = FALSE
      )
    }
  }
  
  invisible(
    list(
      split_by = split_by,
      organize_by = organize_by,
      hemisphere = hemisphere
    )
  )
}


# Internal: return ISO week information for a date.
#
# ISO weeks start on Monday. Week 1 is the week containing 4 January.
.copernicus_iso_week <- function(x) {
  
  x <- as.Date(x)
  
  weekday <- as.POSIXlt(x)$wday
  
  days_since_monday <- (
    weekday + 6L
  ) %% 7L
  
  week_start <- x - days_since_monday
  
  thursday <- week_start + 3
  
  iso_year <- as.integer(
    format(
      thursday,
      "%Y"
    )
  )
  
  jan4 <- as.Date(
    paste0(
      iso_year,
      "-01-04"
    )
  )
  
  jan4_weekday <- as.POSIXlt(jan4)$wday
  
  jan4_days_since_monday <- (
    jan4_weekday + 6L
  ) %% 7L
  
  week1_start <- jan4 - jan4_days_since_monday
  
  week <- as.integer(
    (
      as.numeric(
        week_start - week1_start
      ) / 7
    ) + 1
  )
  
  list(
    year = iso_year,
    week = week,
    label = sprintf(
      "week_%02d",
      week
    ),
    start = week_start
  )
}


# Internal: return meteorological season information.
#
# Meteorological seasons use fixed three-month blocks:
#
#   DJF = December-February
#   MAM = March-May
#   JJA = June-August
#   SON = September-November
#
# These divisions are the same in both hemispheres. The hemisphere determines
# only the seasonal name (winter, spring, summer, autumn).
#
# DJF is assigned to the year containing January and February. For example,
# December 2025 belongs to DJF 2026.
.copernicus_season <- function(
    x,
    hemisphere = c("north", "south")
) {
  
  hemisphere <- match.arg(hemisphere)
  
  x <- as.Date(x)
  
  year <- as.integer(
    format(
      x,
      "%Y"
    )
  )
  
  month <- as.integer(
    format(
      x,
      "%m"
    )
  )
  
  block <- if (month %in% c(12, 1, 2)) {
    "DJF"
  } else if (month %in% 3:5) {
    "MAM"
  } else if (month %in% 6:8) {
    "JJA"
  } else {
    "SON"
  }
  
  season_year <- if (
    identical(block, "DJF") &&
    month == 12
  ) {
    year + 1L
  } else {
    year
  }
  
  north_names <- c(
    DJF = "winter",
    MAM = "spring",
    JJA = "summer",
    SON = "autumn"
  )
  
  south_names <- c(
    DJF = "summer",
    MAM = "autumn",
    JJA = "winter",
    SON = "spring"
  )
  
  season_name <- if (
    identical(hemisphere, "north")
  ) {
    north_names[[block]]
  } else {
    south_names[[block]]
  }
  
  list(
    year = season_year,
    block = block,
    season = season_name,
    label = paste0(
      season_name,
      "_",
      block
    )
  )
}


# Internal: convert an ISO datetime character value to UTC POSIXct.
.copernicus_as_posixct <- function(x) {
  
  if (is.null(x)) {
    return(NULL)
  }
  
  as.POSIXct(
    x,
    format = "%Y-%m-%dT%H:%M:%S",
    tz = "UTC"
  )
}


# Internal: create request periods for package-managed week or season splits.
.copernicus_custom_split_periods <- function(
    start_datetime,
    end_datetime,
    split_by,
    hemisphere
) {
  
  if (!split_by %in% c("week", "season")) {
    stop(
      "Internal error: custom split periods are only used for ",
      "`split_by = \"week\"` or `split_by = \"season\"`.",
      call. = FALSE
    )
  }
  
  if (is.null(start_datetime) ||
      is.null(end_datetime)) {
    stop(
      "`start_datetime` and `end_datetime` are required when splitting ",
      "Copernicus Marine downloads by week or season.",
      call. = FALSE
    )
  }
  
  start <- .copernicus_as_posixct(
    start_datetime
  )
  
  end <- .copernicus_as_posixct(
    end_datetime
  )
  
  if (identical(split_by, "week")) {
    
    start_date <- as.Date(
      start,
      tz = "UTC"
    )
    
    week_info <- .copernicus_iso_week(
      start_date
    )
    
    boundary <- as.POSIXct(
      week_info$start,
      tz = "UTC"
    )
    
    boundaries <- boundary
    
    while (tail(boundaries, 1) <= end) {
      boundaries <- c(
        boundaries,
        tail(boundaries, 1) + 7 * 86400
      )
    }
    
    periods <- vector(
      "list",
      length(boundaries) - 1L
    )
    
    for (i in seq_along(periods)) {
      
      period_start <- max(
        start,
        boundaries[i]
      )
      
      period_end <- min(
        end,
        boundaries[i + 1L] - 1
      )
      
      info <- .copernicus_iso_week(
        as.Date(
          boundaries[i],
          tz = "UTC"
        )
      )
      
      periods[[i]] <- list(
        start = format(
          period_start,
          "%Y-%m-%dT%H:%M:%S",
          tz = "UTC"
        ),
        end = format(
          period_end,
          "%Y-%m-%dT%H:%M:%S",
          tz = "UTC"
        ),
        label = paste0(
          info$year,
          "_",
          info$label
        )
      )
    }
    
    return(periods)
  }
  
  # Meteorological seasons:
  # DJF, MAM, JJA, SON.
  start_date <- as.Date(
    start,
    tz = "UTC"
  )
  
  start_year <- as.integer(
    format(
      start_date,
      "%Y"
    )
  )
  
  start_month <- as.integer(
    format(
      start_date,
      "%m"
    )
  )
  
  if (start_month %in% c(12, 1, 2)) {
    
    boundary_year <- if (
      start_month == 12
    ) {
      start_year
    } else {
      start_year - 1L
    }
    
    boundary_month <- 12L
    
  } else if (start_month %in% 3:5) {
    
    boundary_year <- start_year
    boundary_month <- 3L
    
  } else if (start_month %in% 6:8) {
    
    boundary_year <- start_year
    boundary_month <- 6L
    
  } else {
    
    boundary_year <- start_year
    boundary_month <- 9L
  }
  
  boundary <- as.POSIXct(
    sprintf(
      "%04d-%02d-01 00:00:00",
      boundary_year,
      boundary_month
    ),
    tz = "UTC"
  )
  
  boundaries <- boundary
  
  repeat {
    
    next_boundary <- seq(
      tail(boundaries, 1),
      by = "3 months",
      length.out = 2
    )[2]
    
    boundaries <- c(
      boundaries,
      next_boundary
    )
    
    if (next_boundary > end) {
      break
    }
  }
  
  periods <- vector(
    "list",
    length(boundaries) - 1L
  )
  
  for (i in seq_along(periods)) {
    
    period_start <- max(
      start,
      boundaries[i]
    )
    
    period_end <- min(
      end,
      boundaries[i + 1L] - 1
    )
    
    info <- .copernicus_season(
      as.Date(
        boundaries[i],
        tz = "UTC"
      ),
      hemisphere = hemisphere
    )
    
    periods[[i]] <- list(
      start = format(
        period_start,
        "%Y-%m-%dT%H:%M:%S",
        tz = "UTC"
      ),
      end = format(
        period_end,
        "%Y-%m-%dT%H:%M:%S",
        tz = "UTC"
      ),
      label = paste0(
        info$year,
        "_",
        info$label
      )
    )
  }
  
  periods
}


# Internal: identify the temporal values stored in one downloaded netCDF file.
.copernicus_file_times <- function(path) {
  
  if (!requireNamespace(
    "ncdf4",
    quietly = TRUE
  )) {
    stop(
      "Package 'ncdf4' is required to organize downloaded Copernicus files ",
      "by time.",
      call. = FALSE
    )
  }
  
  nc <- ncdf4::nc_open(path)
  
  on.exit(
    ncdf4::nc_close(nc),
    add = TRUE
  )
  
  time_dim <- .copernicus_nc_time_dim(
    nc
  )
  
  .copernicus_nc_time_values(
    nc,
    time_dim
  )
}


# Internal: identify the data variable represented by one downloaded file.
.copernicus_file_variable <- function(
    path,
    requested_variables = NULL
) {
  
  if (!requireNamespace(
    "ncdf4",
    quietly = TRUE
  )) {
    stop(
      "Package 'ncdf4' is required to organize downloaded Copernicus files ",
      "by variable.",
      call. = FALSE
    )
  }
  
  nc <- ncdf4::nc_open(path)
  
  on.exit(
    ncdf4::nc_close(nc),
    add = TRUE
  )
  
  vars <- names(
    nc$var
  )
  
  if (!is.null(requested_variables)) {
    
    candidates <- intersect(
      requested_variables,
      vars
    )
    
  } else {
    
    candidates <- vars
  }
  
  candidates <- candidates[
    !grepl(
      "(bounds?|bnds?|crs|grid_mapping)",
      candidates,
      ignore.case = TRUE
    )
  ]
  
  if (length(candidates) != 1L) {
    stop(
      "Could not identify one unique data variable for file:\n  ",
      path,
      "\n\n",
      "The file contains ",
      length(candidates),
      " possible data variables",
      if (length(candidates)) {
        paste0(
          ":\n  ",
          paste(
            candidates,
            collapse = ", "
          )
        )
      } else {
        "."
      },
      "\n\n",
      "When `organize_by` includes \"variable\", each output file must ",
      "represent exactly one requested environmental variable. ",
      "Multiple requested variables are separated into files automatically.",
      call. = FALSE
    )
  }
  
  candidates[[1]]
}


# Internal: return folder components for one downloaded Copernicus file.
.copernicus_organization_components <- function(
    path,
    source,
    dataset_id,
    requested_variables,
    organize_by,
    hemisphere
) {
  
  if (is.null(organize_by)) {
    return(character())
  }
  
  components <- character()
  
  temporal_levels <- intersect(
    organize_by,
    c(
      "year",
      "season",
      "month",
      "week",
      "day",
      "hour"
    )
  )
  
  times <- NULL
  
  if (length(temporal_levels)) {
    
    times <- .copernicus_file_times(
      path
    )
    
    if (!length(times)) {
      stop(
        "No time values were found in downloaded file:\n  ",
        path,
        call. = FALSE
      )
    }
  }
  
  variable <- NULL
  
  if ("variable" %in% organize_by) {
    variable <- .copernicus_file_variable(
      path,
      requested_variables = requested_variables
    )
  }
  
  for (level in organize_by) {
    
    if (identical(level, "service")) {
      components <- c(
        components,
        source
      )
      
    } else if (identical(level, "dataset")) {
      components <- c(
        components,
        dataset_id
      )
      
    } else if (identical(level, "variable")) {
      components <- c(
        components,
        variable
      )
      
    } else if (identical(level, "season")) {
      
      season_info <- lapply(
        as.Date(
          times,
          tz = "UTC"
        ),
        .copernicus_season,
        hemisphere = hemisphere
      )
      
      labels <- unique(
        vapply(
          season_info,
          `[[`,
          character(1),
          "label"
        )
      )
      
      if (length(labels) != 1L) {
        stop(
          "Cannot organize file by season because it contains data from ",
          "multiple meteorological seasons:\n  ",
          path,
          "\n\n",
          "Use `split_by = \"season\"` or a finer temporal split.",
          call. = FALSE
        )
      }
      
      components <- c(
        components,
        labels
      )
      
    } else if (identical(level, "week")) {
      
      week_info <- lapply(
        as.Date(
          times,
          tz = "UTC"
        ),
        .copernicus_iso_week
      )
      
      week_labels <- unique(
        vapply(
          week_info,
          `[[`,
          character(1),
          "label"
        )
      )
      
      if (length(week_labels) != 1L) {
        stop(
          "Cannot organize file by week because it contains data from ",
          "multiple ISO weeks:\n  ",
          path,
          "\n\n",
          "Use `split_by = \"week\"` or a finer temporal split.",
          call. = FALSE
        )
      }
      
      components <- c(
        components,
        week_labels
      )
      
    } else if (identical(level, "year")) {
      
      dates <- as.Date(
        times,
        tz = "UTC"
      )
      
      if ("season" %in% organize_by) {
        
        years <- unique(
          vapply(
            lapply(
              dates,
              .copernicus_season,
              hemisphere = hemisphere
            ),
            `[[`,
            integer(1),
            "year"
          )
        )
        
      } else if ("week" %in% organize_by) {
        
        years <- unique(
          vapply(
            lapply(
              dates,
              .copernicus_iso_week
            ),
            `[[`,
            integer(1),
            "year"
          )
        )
        
      } else {
        
        years <- unique(
          format(
            dates,
            "%Y"
          )
        )
      }
      
      if (length(years) != 1L) {
        stop(
          "Cannot organize file by year because it contains data from ",
          "multiple years:\n  ",
          path,
          "\n\n",
          "Use `split_by = \"year\"` or a finer temporal split.",
          call. = FALSE
        )
      }
      
      components <- c(
        components,
        as.character(
          years
        )
      )
      
    } else if (identical(level, "month")) {
      
      values <- unique(
        format(
          times,
          "%m",
          tz = "UTC"
        )
      )
      
      if (length(values) != 1L) {
        stop(
          "Cannot organize file by month because it contains data from ",
          "multiple months:\n  ",
          path,
          "\n\n",
          "Use `split_by = \"month\"` or a finer temporal split.",
          call. = FALSE
        )
      }
      
      components <- c(
        components,
        values
      )
      
    } else if (identical(level, "day")) {
      
      values <- unique(
        format(
          times,
          "%d",
          tz = "UTC"
        )
      )
      
      full_dates <- unique(
        format(
          times,
          "%Y-%m-%d",
          tz = "UTC"
        )
      )
      
      if (length(full_dates) != 1L) {
        stop(
          "Cannot organize file by day because it contains data from ",
          "multiple dates:\n  ",
          path,
          "\n\n",
          "Use `split_by = \"day\"` or `split_by = \"hour\"`.",
          call. = FALSE
        )
      }
      
      components <- c(
        components,
        values
      )
      
    } else if (identical(level, "hour")) {
      
      values <- unique(
        format(
          times,
          "%H",
          tz = "UTC"
        )
      )
      
      full_hours <- unique(
        format(
          times,
          "%Y-%m-%dT%H",
          tz = "UTC"
        )
      )
      
      if (length(full_hours) != 1L) {
        stop(
          "Cannot organize file by hour because it contains data from ",
          "multiple hours:\n  ",
          path,
          "\n\n",
          "Use `split_by = \"hour\"`.",
          call. = FALSE
        )
      }
      
      components <- c(
        components,
        values
      )
    }
  }
  
  components
}


# Internal: make a safe file or directory component.
.copernicus_safe_component <- function(x) {
  
  x <- as.character(x)
  
  if (length(x) != 1L ||
      is.na(x) ||
      !nzchar(trimws(x))) {
    stop(
      "Could not create a valid file or directory name.",
      call. = FALSE
    )
  }
  x <- gsub(
    '[<>:"/\\\\|?*]',
    "_",
    x
  )
  
  x <- gsub(
    "[[:space:]]+",
    "_",
    x
  )
  
  x <- gsub(
    "_+",
    "_",
    x
  )
  
  x <- sub(
    "[. ]+$",
    "",
    x
  )
  
  if (!nzchar(x)) {
    stop(
      "Could not create a valid file or directory name.",
      call. = FALSE
    )
  }
  
  x
}


# Internal: create a filename for one package-managed temporal period.
.copernicus_period_filename <- function(
    filename,
    dataset_id,
    variables,
    label
) {
  
  if (!is.null(filename)) {
    
    stem <- tools::file_path_sans_ext(
      basename(
        filename
      )
    )
    
  } else if (length(variables) == 1L) {
    
    stem <- paste0(
      dataset_id,
      "_",
      variables[[1]]
    )
    
  } else {
    
    stem <- dataset_id
  }
  
  stem <- .copernicus_safe_component(
    stem
  )
  
  label <- .copernicus_safe_component(
    label
  )
  
  paste0(
    stem,
    "_",
    label,
    ".nc"
  )
}


# Internal: run one standalone Copernicus Marine Toolbox command.
#
# Output from the official Toolbox is always captured so this package can
# inspect warnings, parse the returned JSON, and validate the resulting files.
.copernicus_run_standalone <- function(
    executable,
    args,
    dataset_id
) {
  
  output <- tryCatch(
    suppressWarnings(
      system2(
        executable,
        args = args,
        stdout = TRUE,
        stderr = TRUE
      )
    ),
    error = function(e) {
      stop(
        "The standalone Copernicus Marine Toolbox could not be started.\n\n",
        "Executable:\n  ",
        executable,
        "\n\n",
        "Dataset:\n  ",
        dataset_id,
        "\n\n",
        "Original error:\n  ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  status <- attr(
    output,
    "status"
  )
  
  if (is.null(status)) {
    status <- 0L
  }
  
  list(
    status = as.integer(status),
    output = as.character(output)
  )
}

# Internal: parse the file information returned by a Copernicus Marine
# standalone dry-run.
#
# The Toolbox may print INFO/WARNING messages before the JSON payload.
# A dry-run returns detailed information about each planned output file,
# including its environmental variable and coordinate extents. These
# metadata are used to determine expected destinations and to validate
# existing files before they are reused.
.copernicus_parse_dry_run <- function(output) {
  
  if (!length(output)) {
    stop(
      "The Copernicus Marine dry-run returned no output.",
      call. = FALSE
    )
  }
  
  json_start <- grep(
    "^\\s*\\[\\s*$",
    output
  )
  
  if (!length(json_start)) {
    stop(
      "Could not find the JSON download plan in the Copernicus Marine ",
      "dry-run output.",
      call. = FALSE
    )
  }
  
  json_text <- paste(
    output[json_start[[1]]:length(output)],
    collapse = "\n"
  )
  
  plan <- tryCatch(
    jsonlite::fromJSON(
      json_text,
      simplifyVector = FALSE
    ),
    error = function(e) {
      stop(
        "Could not parse the Copernicus Marine dry-run download plan.\n\n",
        "Original error:\n  ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  if (!is.list(plan) || !length(plan)) {
    stop(
      "The Copernicus Marine dry-run returned an empty download plan.",
      call. = FALSE
    )
  }
  
  required_fields <- c(
    "filename",
    "variables",
    "coordinates_extent"
  )
  
  for (i in seq_along(plan)) {
    
    missing_fields <- setdiff(
      required_fields,
      names(plan[[i]])
    )
    
    if (length(missing_fields)) {
      stop(
        "Copernicus Marine dry-run item ",
        i,
        " is missing required field",
        if (length(missing_fields) > 1L) "s" else "",
        ":\n  ",
        paste(
          missing_fields,
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    
    if (length(plan[[i]]$variables) != 1L) {
      stop(
        "Copernicus Marine dry-run item ",
        i,
        " does not represent exactly one environmental variable.\n\n",
        "This package requires one environmental variable per output file.",
        call. = FALSE
      )
    }
  }
  
  plan
}


# Internal: extract one coordinate extent from a Copernicus Marine dry-run
# item.
.copernicus_plan_extent <- function(
    plan_item,
    coordinate_id
) {
  
  extents <- plan_item$coordinates_extent
  
  if (!length(extents)) {
    return(NULL)
  }
  
  ids <- vapply(
    extents,
    function(x) {
      if (is.null(x$coordinate_id)) {
        NA_character_
      } else {
        as.character(x$coordinate_id)
      }
    },
    character(1)
  )
  
  match_index <- which(
    ids == coordinate_id
  )
  
  if (!length(match_index)) {
    return(NULL)
  }
  
  extent <- extents[[match_index[[1]]]]
  
  list(
    minimum = extent$minimum,
    maximum = extent$maximum,
    unit = if (is.null(extent$unit)) {
      NULL
    } else {
      extent$unit
    }
  )
}


# Internal: create the final expected path for one planned Copernicus Marine
# file.
#
# This uses the dry-run filename together with the requested directory
# organization. No data file needs to be downloaded in order to construct
# this path.
.copernicus_planned_destination <- function(
    plan_item,
    output_dir,
    source,
    dataset_id,
    organize_by,
    hemisphere
) {
  
  variable <- as.character(
    plan_item$variables[[1]]
  )
  
  time_extent <- .copernicus_plan_extent(
    plan_item,
    "time"
  )
  
  if (is.null(time_extent)) {
    stop(
      "The dry-run plan for `",
      plan_item$filename,
      "` does not contain a time extent.",
      call. = FALSE
    )
  }
  
  representative_time <- .copernicus_as_posixct(
    as.character(
      time_extent$minimum
    )
  )
  
  components <- character()
  
  for (level in organize_by) {
    
    value <- switch(
      level,
      
      service = source,
      
      dataset = dataset_id,
      
      variable = variable,
      
      year = {
        if ("week" %in% organize_by) {
          as.character(
            .copernicus_iso_week(
              representative_time
            )$year
          )
        } else if ("season" %in% organize_by) {
          as.character(
            .copernicus_season(
              representative_time,
              hemisphere = hemisphere
            )$year
          )
        } else {
          format(
            representative_time,
            "%Y",
            tz = "UTC"
          )
        }
      },
      
      season = .copernicus_season(
        representative_time,
        hemisphere = hemisphere
      )$label,
      
      month = format(
        representative_time,
        "%m",
        tz = "UTC"
      ),
      
      week = .copernicus_iso_week(
        representative_time
      )$label,
      
      day = format(
        representative_time,
        "%d",
        tz = "UTC"
      ),
      
      hour = format(
        representative_time,
        "%H",
        tz = "UTC"
      )
    )
    
    components <- c(
      components,
      .copernicus_safe_component(
        value
      )
    )
  }
  
  destination_dir <- if (length(components)) {
    do.call(
      file.path,
      as.list(
        c(
          output_dir,
          components
        )
      )
    )
  } else {
    output_dir
  }
  
  file.path(
    destination_dir,
    basename(
      plan_item$filename
    )
  )
}

# Internal: check whether an existing NetCDF file can safely satisfy one item
# in a Copernicus Marine dry-run plan.
#
# Only NetCDF metadata and coordinate variables are inspected. Environmental
# data values are not loaded into memory, keeping this check inexpensive even
# when files themselves are large.
.copernicus_existing_file_matches <- function(
    path,
    plan_item,
    tolerance = 1e-6
) {
  
  if (!file.exists(path)) {
    return(FALSE)
  }
  
  info <- file.info(
    path
  )
  
  if (is.na(info$size) ||
      info$size <= 0) {
    return(FALSE)
  }
  
  expected_variable <- as.character(
    plan_item$variables[[1]]
  )
  
  nc <- tryCatch(
    ncdf4::nc_open(
      path
    ),
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(nc)) {
    return(FALSE)
  }
  
  on.exit(
    ncdf4::nc_close(
      nc
    ),
    add = TRUE
  )
  
  if (!expected_variable %in% names(nc$var)) {
    return(FALSE)
  }
  
  coordinate_candidates <- list(
    longitude = c(
      "longitude",
      "lon"
    ),
    latitude = c(
      "latitude",
      "lat"
    ),
    depth = c(
      "depth",
      "elevation"
    ),
    time = c(
      "time"
    )
  )
  
  find_coordinate <- function(candidates) {
    
    all_names <- unique(
      c(
        names(nc$dim),
        names(nc$var)
      )
    )
    
    hit <- candidates[
      candidates %in% all_names
    ]
    
    if (!length(hit)) {
      return(NULL)
    }
    
    hit[[1]]
  }
  
  read_coordinate <- function(name) {
    
    if (name %in% names(nc$dim)) {
      return(
        nc$dim[[name]]$vals
      )
    }
    
    if (name %in% names(nc$var)) {
      return(
        ncdf4::ncvar_get(
          nc,
          name
        )
      )
    }
    
    NULL
  }
  
  numeric_match <- function(
    actual,
    expected_min,
    expected_max
  ) {
    
    if (is.null(actual) ||
        !length(actual)) {
      return(FALSE)
    }
    
    actual <- as.numeric(
      actual
    )
    
    actual <- actual[
      is.finite(actual)
    ]
    
    if (!length(actual)) {
      return(FALSE)
    }
    
    expected_min <- as.numeric(
      expected_min
    )
    
    expected_max <- as.numeric(
      expected_max
    )
    
    isTRUE(
      all.equal(
        min(actual),
        expected_min,
        tolerance = tolerance,
        check.attributes = FALSE
      )
    ) &&
      isTRUE(
        all.equal(
          max(actual),
          expected_max,
          tolerance = tolerance,
          check.attributes = FALSE
        )
      )
  }
  
  # -----------------------------------------------------------------------
  # Longitude, latitude, and vertical coordinate
  # -----------------------------------------------------------------------
  
  for (coordinate_id in c(
    "longitude",
    "latitude",
    "depth"
  )) {
    
    expected <- .copernicus_plan_extent(
      plan_item,
      coordinate_id
    )
    
    if (is.null(expected)) {
      next
    }
    
    coordinate_name <- find_coordinate(
      coordinate_candidates[[coordinate_id]]
    )
    
    if (is.null(coordinate_name)) {
      return(FALSE)
    }
    
    actual <- read_coordinate(
      coordinate_name
    )
    
    if (!numeric_match(
      actual = actual,
      expected_min = expected$minimum,
      expected_max = expected$maximum
    )) {
      return(FALSE)
    }
  }
  
  # -----------------------------------------------------------------------
  # Time
  #
  # Reuse the package's existing NetCDF time reader so CF time conversion is
  # handled consistently in one place.
  # -----------------------------------------------------------------------
  
  expected_time <- .copernicus_plan_extent(
    plan_item,
    "time"
  )
  
  if (!is.null(expected_time)) {
    
    actual_time <- tryCatch(
      .copernicus_file_times(
        path
      ),
      error = function(e) {
        NULL
      }
    )
    
    if (is.null(actual_time) ||
        !length(actual_time)) {
      return(FALSE)
    }
    
    expected_min <- .copernicus_as_posixct(
      as.character(
        expected_time$minimum
      )
    )
    
    expected_max <- .copernicus_as_posixct(
      as.character(
        expected_time$maximum
      )
    )
    
    actual_min <- min(
      actual_time,
      na.rm = TRUE
    )
    
    actual_max <- max(
      actual_time,
      na.rm = TRUE
    )
    
    if (!isTRUE(
      all.equal(
        as.numeric(actual_min),
        as.numeric(expected_min),
        tolerance = 1e-6,
        check.attributes = FALSE
      )
    )) {
      return(FALSE)
    }
    
    if (!isTRUE(
      all.equal(
        as.numeric(actual_max),
        as.numeric(expected_max),
        tolerance = 1e-6,
        check.attributes = FALSE
      )
    )) {
      return(FALSE)
    }
  }
  
  TRUE
}

# Internal: move downloaded files into the requested directory structure.
#
# Files are always downloaded to a temporary staging directory first. This
# helper moves them to the final output directory, optionally creating the
# hierarchy requested through `organize_by`.
.copernicus_organize_files <- function(
    paths,
    output_dir,
    source,
    dataset_id,
    requested_variables,
    organize_by,
    hemisphere,
    force,
    quiet
) {
  
  organized <- character(
    length(paths)
  )
  
  for (i in seq_along(paths)) {
    
    components <- .copernicus_organization_components(
      path = paths[i],
      source = source,
      dataset_id = dataset_id,
      requested_variables = requested_variables,
      organize_by = organize_by,
      hemisphere = hemisphere
    )
    
    if (length(components)) {
      components <- vapply(
        components,
        .copernicus_safe_component,
        character(1)
      )
    }
    
    if (length(components)) {
      
      destination_dir <- do.call(
        file.path,
        as.list(
          c(
            output_dir,
            components
          )
        )
      )
      
    } else {
      
      destination_dir <- output_dir
    }
    
    dir.create(
      destination_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    destination <- file.path(
      destination_dir,
      basename(
        paths[i]
      )
    )
    
    source_normalized <- normalizePath(
      paths[i],
      winslash = "/",
      mustWork = FALSE
    )
    
    destination_normalized <- normalizePath(
      destination,
      winslash = "/",
      mustWork = FALSE
    )
    
    same_path <- identical(
      source_normalized,
      destination_normalized
    )
    
    if (!same_path) {
      
      if (file.exists(destination)) {
        
        if (!force) {
          stop(
            "Output file already exists:\n  ",
            destination,
            "\n\n",
            "The new download has been preserved in the temporary staging ",
            "directory rather than overwriting the existing file.\n\n",
            "Use `force = TRUE` if the existing file should be replaced.",
            call. = FALSE
          )
        }
        
        unlink(
          destination,
          force = TRUE
        )
      }
      
      moved <- file.rename(
        paths[i],
        destination
      )
      
      if (!moved) {
        
        copied <- file.copy(
          paths[i],
          destination,
          overwrite = force
        )
        
        if (!copied) {
          stop(
            "Could not move downloaded file to:\n  ",
            destination,
            call. = FALSE
          )
        }
        
        unlink(
          paths[i],
          force = TRUE
        )
      }
    }
    
    organized[i] <- normalizePath(
      destination,
      winslash = "/",
      mustWork = TRUE
    )
    
    if (!quiet) {
      
      destination_relative <- if (length(components)) {
        paste(
          components,
          collapse = "/"
        )
      } else {
        "."
      }
      
      message(
        "  File ",
        i,
        " of ",
        length(paths),
        " -> ",
        destination_relative,
        "/",
        basename(
          destination
        )
      )
    }
  }
  
  organized
}


# -------------------------------------------------------------------------
# Copernicus Marine standalone backend
# -------------------------------------------------------------------------

# Internal: Copernicus Marine standalone backend.
#
# The standalone backend always stores one requested environmental variable per
# output NetCDF file. When several variables are requested, the official
# Toolbox `split-on --on-variables` mechanism is used automatically.
#
# Native Toolbox temporal splitting is used for year/month/day/hour.
# ISO-week and meteorological-season splitting are managed by this package by
# issuing one Toolbox request per temporal period.
.copernicus_load_marine_standalone <- function(
    dataset_id,
    variables,
    start_datetime,
    end_datetime,
    xmin,
    xmax,
    ymin,
    ymax,
    depth_min,
    depth_max,
    request,
    output_dir,
    filename,
    split_by,
    organize_by,
    hemisphere,
    concurrent_processes,
    compression,
    dataset_version,
    dataset_part,
    force,
    quiet
) {
  
  executable <- .copernicus_find_executable()
  
  if (is.null(executable)) {
    stop(
      "The standalone Copernicus Marine Toolbox backend was selected, ",
      "but no usable Toolbox executable could be found.\n\n",
      "Run `copernicus_setup()` to install the supported standalone Toolbox, ",
      "then check the installation with `copernicus_status()`.",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop(
      "Package 'ncdf4' is required to validate Copernicus Marine downloads.\n\n",
      "Install it with:\n\n",
      "  install.packages(\"ncdf4\")",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      "Package 'jsonlite' is required to read Copernicus Marine download ",
      "plans.\n\n",
      "Install it with:\n\n",
      "  install.packages(\"jsonlite\")",
      call. = FALSE
    )
  }
  
  if (!is.null(compression)) {
    
    if (!is.numeric(compression) ||
        length(compression) != 1 ||
        is.na(compression) ||
        compression < 0 ||
        compression > 9 ||
        compression != as.integer(compression)) {
      
      stop(
        "`compression` must be a single integer from 0 to 9.",
        call. = FALSE
      )
    }
  }
  
  if (!is.null(dataset_version)) {
    .copernicus_validate_string(
      dataset_version,
      "dataset_version"
    )
  }
  
  if (!is.null(dataset_part)) {
    .copernicus_validate_string(
      dataset_part,
      "dataset_part"
    )
  }
  
  supported_request_names <- c(
    "service",
    "coordinates_selection_method"
  )
  
  unsupported_request_names <- setdiff(
    names(request),
    supported_request_names
  )
  
  if (length(unsupported_request_names)) {
    
    stop(
      "The standalone Copernicus Marine backend does not currently support ",
      "the following `request` argument",
      if (length(unsupported_request_names) > 1) "s" else "",
      ":\n  ",
      paste(
        unsupported_request_names,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  if (!is.null(request$service)) {
    .copernicus_validate_string(
      request$service,
      "request$service"
    )
  }
  
  if (!is.null(request$coordinates_selection_method)) {
    
    allowed_selection_methods <- c(
      "inside",
      "strict-inside",
      "nearest",
      "outside"
    )
    
    if (!request$coordinates_selection_method %in%
        allowed_selection_methods) {
      
      stop(
        "`request$coordinates_selection_method` must be one of:\n  ",
        paste(
          allowed_selection_methods,
          collapse = ", "
        ),
        call. = FALSE
      )
    }
  }
  
  custom_split <- !is.null(split_by) &&
    split_by %in% c(
      "week",
      "season"
    )
  
  native_time_split <- !is.null(split_by) &&
    split_by %in% c(
      "year",
      "month",
      "day",
      "hour"
    )
  
  multiple_variables <- length(variables) > 1L
  
  staging_dir <- tempfile(
    pattern = ".copernicus_staging_",
    tmpdir = output_dir
  )
  
  dir.create(
    staging_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  completed <- FALSE
  
  on.exit(
    {
      if (completed && dir.exists(staging_dir)) {
        unlink(
          staging_dir,
          recursive = TRUE,
          force = TRUE
        )
      }
    },
    add = TRUE
  )
  
  start_time <- Sys.time()
  
  
  
  # -----------------------------------------------------------------------
  # Decide whether pre-planning is actually necessary.
  #
  # A dry-run is only needed when existing NetCDF files may be reusable.
  # For a completely new download there is nothing to compare, so running
  # the Toolbox twice (dry-run + real download) would only add substantial
  # overhead.
  #
  # `force = TRUE` also skips pre-planning because the user has explicitly
  # requested fresh downloads regardless of existing files.
  #
  # This initial check is deliberately conservative: if any NetCDF file is
  # present anywhere below output_dir, the package performs the rigorous
  # dry-run comparison rather than risking an incorrect reuse decision.
  # -----------------------------------------------------------------------
  
  existing_netcdf_candidates <- list.files(
    output_dir,
    pattern = "\\.(nc|nc4)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  direct_download <-
    isTRUE(force) ||
    length(existing_netcdf_candidates) == 0L
  
  if (!quiet) {
    
    message(
      "\nCopernicus Marine download\n",
      "---------------------------\n",
      "Dataset:   ",
      dataset_id,
      "\n",
      "Variables: ",
      paste(
        variables,
        collapse = ", "
      ),
      "\n",
      if (!is.null(start_datetime)) {
        paste0(
          "Time:      ",
          start_datetime,
          " to ",
          end_datetime,
          "\n"
        )
      } else {
        ""
      },
      if (!is.null(xmin)) {
        paste0(
          "Area:      ",
          xmin,
          " to ",
          xmax,
          " longitude; ",
          ymin,
          " to ",
          ymax,
          " latitude\n"
        )
      } else {
        ""
      },
      if (!is.null(depth_min) || !is.null(depth_max)) {
        paste0(
          "Depth:     ",
          if (is.null(depth_min)) "dataset minimum" else depth_min,
          " to ",
          if (is.null(depth_max)) "dataset maximum" else depth_max,
          " m\n"
        )
      } else {
        ""
      },
      "Variable files: one variable per NetCDF\n",
      "Temporal split: ",
      if (is.null(split_by)) {
        "none"
      } else {
        split_by
      },
      "\n",
      "Organization: ",
      if (is.null(organize_by)) {
        "flat output directory"
      } else {
        paste(
          organize_by,
          collapse = " / "
        )
      },
      "\n",
      "Backend:   standalone\n"
    )
    
    if (direct_download) {
      
      if (isTRUE(force)) {
        
        message(
          "[1/4] Preparing forced download...\n",
          "  `force = TRUE`: pre-planning is skipped because all requested ",
          "files will be downloaded again."
        )
        
      } else {
        
        message(
          "No existing Copernicus Marine NetCDF files found.\n",
          "Starting a new download without pre-planning...\n\n",
          "[1/4] Preparing download..."
        )
      }
      
    } else {
      
      message(
        "Existing Copernicus Marine files were found.\n\n",
        "Checking existing files against the Copernicus Marine request...\n",
        "This may take a little while for large requests because the Toolbox ",
        "needs to identify the expected files and this package checks which ",
        "ones can be safely reused or need to be replaced.\n\n",
        "[1/4] Planning download and checking existing files..."
      )
    }
  }
  
  # -----------------------------------------------------------------------
  # Build the common Toolbox arguments.
  #
  # `request_variables` can contain all requested variables during the
  # dry-run planning stage, or one variable when downloading an individual
  # missing file.
  # -----------------------------------------------------------------------
  
  build_args <- function(
    request_output_dir,
    request_start,
    request_end,
    request_filename = NULL,
    request_variables = variables,
    dry_run = FALSE
  ) {
    
    args <- c(
      "subset",
      "--dataset-id",
      dataset_id
    )
    
    if (!is.null(dataset_version)) {
      args <- c(
        args,
        "--dataset-version",
        dataset_version
      )
    }
    
    if (!is.null(dataset_part)) {
      args <- c(
        args,
        "--dataset-part",
        dataset_part
      )
    }
    
    for (variable in request_variables) {
      args <- c(
        args,
        "--variable",
        variable
      )
    }
    
    if (!is.null(xmin)) {
      args <- c(
        args,
        "--minimum-longitude",
        format(
          xmin,
          scientific = FALSE,
          trim = TRUE
        ),
        "--maximum-longitude",
        format(
          xmax,
          scientific = FALSE,
          trim = TRUE
        ),
        "--minimum-latitude",
        format(
          ymin,
          scientific = FALSE,
          trim = TRUE
        ),
        "--maximum-latitude",
        format(
          ymax,
          scientific = FALSE,
          trim = TRUE
        )
      )
    }
    
    if (!is.null(depth_min)) {
      args <- c(
        args,
        "--minimum-depth",
        format(
          depth_min,
          scientific = FALSE,
          trim = TRUE
        )
      )
    }
    
    if (!is.null(depth_max)) {
      args <- c(
        args,
        "--maximum-depth",
        format(
          depth_max,
          scientific = FALSE,
          trim = TRUE
        )
      )
    }
    
    if (!is.null(request_start)) {
      args <- c(
        args,
        "--start-datetime",
        request_start
      )
    }
    
    if (!is.null(request_end)) {
      args <- c(
        args,
        "--end-datetime",
        request_end
      )
    }
    
    args <- c(
      args,
      "--output-directory",
      request_output_dir,
      "--file-format",
      "netcdf"
    )
    
    if (!is.null(request_filename)) {
      args <- c(
        args,
        "--output-filename",
        request_filename
      )
    }
    
    if (!is.null(compression)) {
      args <- c(
        args,
        "--netcdf-compression-level",
        as.character(
          as.integer(
            compression
          )
        )
      )
    }
    
    if (!is.null(request$service)) {
      args <- c(
        args,
        "--service",
        request$service
      )
    }
    
    if (!is.null(request$coordinates_selection_method)) {
      args <- c(
        args,
        "--coordinates-selection-method",
        request$coordinates_selection_method
      )
    }
    
    if (isTRUE(force) && !dry_run) {
      args <- c(
        args,
        "--overwrite"
      )
    }
    
    if (dry_run) {
      args <- c(
        args,
        "--dry-run"
      )
    }
    if (isTRUE(quiet)) {
      args <- c(
        args,
        "--log-level",
        "ERROR"
      )
    }
    args
  }
  
  # -----------------------------------------------------------------------
  # Direct-download path.
  #
  # When there are no existing NetCDF files to reuse, or when `force = TRUE`,
  # the real Toolbox request is executed immediately without a preliminary
  # dry-run.
  #
  # The JSON returned by the completed Toolbox request becomes the
  # authoritative description of the downloaded files. These files are then
  # validated, organized and validated again at their final destinations.
  #
  # Native year/month/day/hour splitting is handled by a single Toolbox
  # request. Week and season remain package-managed because the Toolbox does
  # not provide those temporal split units.
  # -----------------------------------------------------------------------
  
  if (direct_download) {
    
    if (!quiet) {
      message(
        "[2/4] Downloading data from Copernicus Marine..."
      )
    }
    
    download_items <- list()
    
    # ---------------------------------------------------------------------
    # Native or unsplit request: one Toolbox call.
    # ---------------------------------------------------------------------
    
    if (!custom_split) {
      
      download_args <- build_args(
        request_output_dir = staging_dir,
        request_start = start_datetime,
        request_end = end_datetime,
        request_filename = filename,
        request_variables = variables,
        dry_run = FALSE
      )
      
      needs_split_on <-
        multiple_variables ||
        native_time_split
      
      if (needs_split_on) {
        
        download_args <- c(
          download_args,
          "split-on"
        )
        
        if (multiple_variables) {
          download_args <- c(
            download_args,
            "--on-variables"
          )
        }
        
        if (native_time_split) {
          download_args <- c(
            download_args,
            "--on-time",
            split_by
          )
        }
        
        if (!is.null(concurrent_processes)) {
          download_args <- c(
            download_args,
            "--concurrent-processes",
            as.character(
              as.integer(
                concurrent_processes
              )
            )
          )
        }
      }
      
      download_result <- .copernicus_run_standalone(
        executable = executable,
        args = download_args,
        dataset_id = dataset_id
      )
      
      if (!identical(download_result$status, 0L)) {
        
        toolbox_output <- paste(
          download_result$output,
          collapse = "\n"
        )
        
        stop(
          "The Copernicus Marine download failed.\n\n",
          "Dataset:\n  ",
          dataset_id,
          "\n\n",
          "The official Toolbox returned exit status ",
          download_result$status,
          ".\n\n",
          if (nzchar(toolbox_output)) {
            paste0(
              "Toolbox output:\n",
              toolbox_output
            )
          } else {
            ""
          },
          call. = FALSE
        )
      }
      
      if (!quiet) {
        
        toolbox_warnings <- grep(
          "WARNING",
          download_result$output,
          value = TRUE
        )
        
        if (length(toolbox_warnings)) {
          
          toolbox_warnings <- unique(
            trimws(
              sub(
                "^WARNING\\s*-\\s*\\d{4}-\\d{2}-\\d{2}T[^ ]+\\s*-\\s*",
                "",
                toolbox_warnings
              )
            )
          )
          
          for (warning_message in toolbox_warnings) {
            message(
              "        WARNING: ",
              warning_message
            )
          }
        }
      }
      
      # The completed Toolbox request does not report output filenames.
      # Once the process has completed successfully, use the NetCDF files
      # actually written to the staging directory as the authoritative result.
      staging_candidates <- list.files(
        staging_dir,
        pattern = "\\.(nc|nc4)$",
        full.names = TRUE,
        recursive = TRUE,
        ignore.case = TRUE
      )
      if (!length(staging_candidates)) {
        stop(
          "The Copernicus Marine Toolbox completed the custom temporal split ",
          "request successfully, but no NetCDF files were found in the ",
          "temporary staging directory.",
          call. = FALSE
        )
      }
      download_items <- lapply(
        staging_candidates,
        function(path) {
          list(
            filename = basename(path)
          )
        }
      )
      
    } else {
      
      # -------------------------------------------------------------------
      # Week / season request.
      #
      # These periods are generated locally and downloaded directly, without
      # performing one dry-run per period.
      # -------------------------------------------------------------------
      
      periods <- .copernicus_custom_split_periods(
        start_datetime = start_datetime,
        end_datetime = end_datetime,
        split_by = split_by,
        hemisphere = hemisphere
      )
      
      for (period_i in seq_along(periods)) {
        
        period <- periods[[period_i]]
        
        if (!quiet) {
          message(
            "  [",
            period_i,
            "/",
            length(periods),
            "] ",
            period$label
          )
        }
        
        period_filename <- .copernicus_period_filename(
          filename = filename,
          dataset_id = dataset_id,
          variables = variables,
          label = period$label
        )
        
        period_args <- build_args(
          request_output_dir = staging_dir,
          request_start = period$start,
          request_end = period$end,
          request_filename = period_filename,
          request_variables = variables,
          dry_run = FALSE
        )
        
        if (multiple_variables) {
          
          period_args <- c(
            period_args,
            "split-on",
            "--on-variables"
          )
          
          if (!is.null(concurrent_processes)) {
            period_args <- c(
              period_args,
              "--concurrent-processes",
              as.character(
                as.integer(
                  concurrent_processes
                )
              )
            )
          }
        }
        
        period_result <- .copernicus_run_standalone(
          executable = executable,
          args = period_args,
          dataset_id = dataset_id
        )
        
        if (!identical(period_result$status, 0L)) {
          
          toolbox_output <- paste(
            period_result$output,
            collapse = "\n"
          )
          
          stop(
            "The Copernicus Marine download failed while processing ",
            split_by,
            " period ",
            period_i,
            " of ",
            length(periods),
            ".\n\n",
            "Period:\n  ",
            period$label,
            "\n\n",
            if (nzchar(toolbox_output)) {
              paste0(
                "Toolbox output:\n",
                toolbox_output
              )
            } else {
              ""
            },
            call. = FALSE
          )
        }
        
        if (!quiet) {
          
          toolbox_warnings <- grep(
            "WARNING",
            period_result$output,
            value = TRUE
          )
          
          if (length(toolbox_warnings)) {
            
            toolbox_warnings <- unique(
              trimws(
                sub(
                  "^WARNING\\s*-\\s*\\d{4}-\\d{2}-\\d{2}T[^ ]+\\s*-\\s*",
                  "",
                  toolbox_warnings
                )
              )
            )
            
            for (warning_message in toolbox_warnings) {
              message(
                "        WARNING: ",
                warning_message
              )
            }
          }
        }
        
      }
    }
    
    # Completed Toolbox requests do not reliably report output filenames.
    # If no download items are available yet, use the NetCDF files
    # actually written to the staging directory as the authoritative result.
    if (!length(download_items)) {
      staging_candidates <- list.files(
        staging_dir,
        pattern = "\\.(nc|nc4)$",
        full.names = TRUE,
        recursive = TRUE,
        ignore.case = TRUE
      )
      if (length(staging_candidates)) {
        download_items <- lapply(
          staging_candidates,
          function(path) {
            list(
              filename = basename(path)
            )
          }
        )
      } else {
        stop(
          "The Copernicus Marine Toolbox completed the request successfully, ",
          "but no NetCDF files were found in the temporary staging directory.",
          call. = FALSE
        )
      }
    }
    
    # ---------------------------------------------------------------------
    # Match the expected download items to the NetCDF files actually
    # written to the staging directory.
    # ---------------------------------------------------------------------
    
    staging_candidates <- list.files(
      staging_dir,
      pattern = "\\.(nc|nc4)$",
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )
    
    if (length(staging_candidates) != length(download_items)) {
      stop(
        "The Copernicus Marine Toolbox completed the download, but the ",
        "number of NetCDF files found does not match the reported result.\n\n",
        "Expected download items:\n  ",
        length(download_items),
        "\n\n",
        "NetCDF files found:\n  ",
        length(staging_candidates),
        "\n\n",
        "The downloaded files have been left in the temporary staging ",
        "directory for inspection.",
        call. = FALSE
      )
    }
    
    downloaded_staging_paths <- character(
      length(download_items)
    )
    
    unmatched_candidates <- staging_candidates
    
    for (i in seq_along(download_items)) {
      
      expected_name <- basename(
        download_items[[i]]$filename
      )
      
      filename_matches <- unmatched_candidates[
        basename(unmatched_candidates) == expected_name
      ]
      
      if (length(filename_matches) != 1L) {
        stop(
          "A downloaded Copernicus Marine file could not be matched ",
          "uniquely to the Toolbox result.\n\n",
          "Expected filename:\n  ",
          expected_name,
          "\n\n",
          "Matching files found:\n  ",
          length(filename_matches),
          "\n\n",
          "The downloaded files have been left in the temporary staging ",
          "directory for inspection.",
          call. = FALSE
        )
      }
      
      downloaded_staging_paths[[i]] <- filename_matches[[1]]
      
      unmatched_candidates <- setdiff(
        unmatched_candidates,
        filename_matches[[1]]
      )
    }
    
    # ---------------------------------------------------------------------
    # Validate the downloaded NetCDF files before moving them.
    #
    # A real Toolbox download does not return the detailed coordinate
    # metadata provided by dry-run. Instead, we inspect the files that were
    # actually created and verify that they are readable and contain exactly
    # one of the requested environmental variables.
    # ---------------------------------------------------------------------
    
    downloaded_variables <- character(
      length(downloaded_staging_paths)
    )
    
    for (i in seq_along(downloaded_staging_paths)) {
      
      path <- downloaded_staging_paths[[i]]
      
      if (!file.exists(path)) {
        stop(
          "A Copernicus Marine output file expected from the download does not ",
          "exist:\n  ",
          path,
          call. = FALSE
        )
      }
      
      info <- file.info(
        path
      )
      
      if (is.na(info$size) || info$size <= 0) {
        stop(
          "A downloaded Copernicus Marine NetCDF file is empty:\n  ",
          path,
          call. = FALSE
        )
      }
      
      downloaded_variables[[i]] <- .copernicus_file_variable(
        path = path,
        requested_variables = variables
      )
      
      if (!downloaded_variables[[i]] %in% variables) {
        stop(
          "A downloaded Copernicus Marine file contains an unexpected ",
          "environmental variable.\n\n",
          "File:\n  ",
          path,
          "\n\n",
          "Variable found:\n  ",
          downloaded_variables[[i]],
          "\n\n",
          "Requested variables:\n  ",
          paste(
            variables,
            collapse = ", "
          ),
          call. = FALSE
        )
      }
    }
    
    downloaded_variable_counts <- table(
      factor(
        downloaded_variables,
        levels = variables
      )
    )
    
    if (any(downloaded_variable_counts == 0L)) {
      
      missing_variables <- names(
        downloaded_variable_counts[
          downloaded_variable_counts == 0L
        ]
      )
      
      stop(
        "The Copernicus Marine download completed, but no output file was ",
        "found for requested variable",
        if (length(missing_variables) > 1L) "s" else "",
        ":\n  ",
        paste(
          missing_variables,
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    
    # ---------------------------------------------------------------------
    # Organize downloaded files.
    # ---------------------------------------------------------------------
    
    organized_files <- character(
      length(downloaded_staging_paths)
    )
    
    if (!quiet) {
      message(
        "[3/4] Organizing ",
        length(downloaded_staging_paths),
        if (length(downloaded_staging_paths) == 1L) {
          " downloaded file..."
        } else {
          " downloaded files..."
        }
      )
    }
    
    for (i in seq_along(downloaded_staging_paths)) {
      
      organized <- .copernicus_organize_files(
        paths = downloaded_staging_paths[[i]],
        output_dir = output_dir,
        source = "marine",
        dataset_id = dataset_id,
        requested_variables = variables,
        organize_by = organize_by,
        hemisphere = hemisphere,
        force = TRUE,
        quiet = TRUE
      )
      
      organized_files[[i]] <- organized
      
      if (!quiet) {
        
        organized_normalized <- normalizePath(
          organized,
          winslash = "/",
          mustWork = TRUE
        )
        
        relative_path <- sub(
          paste0(
            "^",
            gsub(
              "([][{}()+*^$|\\\\?.])",
              "\\\\\\1",
              normalizePath(
                output_dir,
                winslash = "/",
                mustWork = TRUE
              )
            ),
            "/?"
          ),
          "",
          organized_normalized
        )
        
        message(
          "  [",
          i,
          "/",
          length(downloaded_staging_paths),
          "] ",
          relative_path
        )
      }
    }
    
    # ---------------------------------------------------------------------
    # Final validation.
    # ---------------------------------------------------------------------
    
    if (!quiet) {
      message(
        "[4/4] Validating ",
        length(download_items),
        if (length(download_items) == 1L) {
          " output file..."
        } else {
          " output files..."
        }
      )
    }
    
    output_files <- normalizePath(
      organized_files,
      winslash = "/",
      mustWork = TRUE
    )
    
    for (i in seq_along(output_files)) {
      
      if (!file.exists(output_files[[i]])) {
        stop(
          "A final Copernicus Marine output file is missing:\n  ",
          output_files[[i]],
          call. = FALSE
        )
      }
      
      final_variable <- .copernicus_file_variable(
        path = output_files[[i]],
        requested_variables = variables
      )
      
      if (!identical(
        final_variable,
        downloaded_variables[[i]]
      )) {
        stop(
          "A Copernicus Marine file changed unexpectedly during ",
          "organization.\n\n",
          "File:\n  ",
          output_files[[i]],
          "\n\n",
          "Variable before organization:\n  ",
          downloaded_variables[[i]],
          "\n\n",
          "Variable after organization:\n  ",
          final_variable,
          call. = FALSE
        )
      }
    }
    
    completed <- TRUE
    
    if (!quiet) {
      
      elapsed <- as.numeric(
        difftime(
          Sys.time(),
          start_time,
          units = "secs"
        )
      )
      
      message(
        "\nDone in ",
        format(
          round(
            elapsed,
            1
          ),
          nsmall = 1
        ),
        " s.\n",
        "Requested:  ",
        length(download_items),
        "\n",
        "Reused:     0\n",
        "Downloaded: ",
        length(download_items)
      )
    }
    
    return(
      sort(
        unique(
          output_files
        )
      )
    )
  }
  
  # -----------------------------------------------------------------------
  # Create the authoritative download plan using the official Toolbox.
  #
  # Native year/month/day/hour splitting is planned by the Toolbox itself.
  # Week and season splitting first creates package-managed temporal periods,
  # and each period is then planned by the Toolbox.
  # -----------------------------------------------------------------------
  
  plan <- list()
  
  if (!custom_split) {
    
    plan_args <- build_args(
      request_output_dir = staging_dir,
      request_start = start_datetime,
      request_end = end_datetime,
      request_filename = filename,
      request_variables = variables,
      dry_run = TRUE
    )
    
    needs_split_on <-
      multiple_variables ||
      native_time_split
    
    if (needs_split_on) {
      
      plan_args <- c(
        plan_args,
        "split-on"
      )
      
      if (multiple_variables) {
        plan_args <- c(
          plan_args,
          "--on-variables"
        )
      }
      
      if (native_time_split) {
        plan_args <- c(
          plan_args,
          "--on-time",
          split_by
        )
      }
      
      if (!is.null(concurrent_processes)) {
        plan_args <- c(
          plan_args,
          "--concurrent-processes",
          as.character(
            as.integer(
              concurrent_processes
            )
          )
        )
      }
    }
    
    plan_result <- .copernicus_run_standalone(
      executable = executable,
      args = plan_args,
      dataset_id = dataset_id
    )
    
    if (!identical(plan_result$status, 0L)) {
      
      toolbox_output <- paste(
        plan_result$output,
        collapse = "\n"
      )
      
      stop(
        "The Copernicus Marine download plan could not be created.\n\n",
        "Dataset:\n  ",
        dataset_id,
        "\n\n",
        "The official Toolbox returned exit status ",
        plan_result$status,
        ".\n\n",
        if (nzchar(toolbox_output)) {
          paste0(
            "Toolbox output:\n",
            toolbox_output
          )
        } else {
          ""
        },
        call. = FALSE
      )
    }
    
    plan <- .copernicus_parse_dry_run(
      plan_result$output
    )
    
  } else {
    
    periods <- .copernicus_custom_split_periods(
      start_datetime = start_datetime,
      end_datetime = end_datetime,
      split_by = split_by,
      hemisphere = hemisphere
    )
    
    for (i in seq_along(periods)) {
      
      period <- periods[[i]]
      
      period_filename <- .copernicus_period_filename(
        filename = filename,
        dataset_id = dataset_id,
        variables = variables,
        label = period$label
      )
      
      period_plan_args <- build_args(
        request_output_dir = staging_dir,
        request_start = period$start,
        request_end = period$end,
        request_filename = period_filename,
        request_variables = variables,
        dry_run = TRUE
      )
      
      if (multiple_variables) {
        
        period_plan_args <- c(
          period_plan_args,
          "split-on",
          "--on-variables"
        )
        
        if (!is.null(concurrent_processes)) {
          period_plan_args <- c(
            period_plan_args,
            "--concurrent-processes",
            as.character(
              as.integer(
                concurrent_processes
              )
            )
          )
        }
      }
      
      period_plan_result <- .copernicus_run_standalone(
        executable = executable,
        args = period_plan_args,
        dataset_id = dataset_id
      )
      
      if (!identical(period_plan_result$status, 0L)) {
        
        toolbox_output <- paste(
          period_plan_result$output,
          collapse = "\n"
        )
        
        stop(
          "The Copernicus Marine download plan failed while processing ",
          split_by,
          " period ",
          i,
          " of ",
          length(periods),
          ".\n\n",
          "Period:\n  ",
          period$label,
          "\n\n",
          if (nzchar(toolbox_output)) {
            paste0(
              "Toolbox output:\n",
              toolbox_output
            )
          } else {
            ""
          },
          call. = FALSE
        )
      }
      
      this_plan <- .copernicus_parse_dry_run(
        period_plan_result$output
      )
      
      plan <- c(
        plan,
        this_plan
      )
    }
  }
  
  if (!length(plan)) {
    stop(
      "The Copernicus Marine Toolbox returned an empty download plan.",
      call. = FALSE
    )
  }
  
  # -----------------------------------------------------------------------
  # Determine the final destination of every planned file and check whether
  # a valid copy is already available there.
  # -----------------------------------------------------------------------
  
  planned_destinations <- vapply(
    plan,
    function(item) {
      .copernicus_planned_destination(
        plan_item = item,
        output_dir = output_dir,
        source = "marine",
        dataset_id = dataset_id,
        organize_by = organize_by,
        hemisphere = hemisphere
      )
    },
    character(1)
  )
  
  destination_exists <- file.exists(
    planned_destinations
  )
  
  reusable <- rep(
    FALSE,
    length(plan)
  )
  
  if (!force) {
    
    for (i in seq_along(plan)) {
      
      if (!destination_exists[[i]]) {
        next
      }
      
      reusable[[i]] <- .copernicus_existing_file_matches(
        path = planned_destinations[[i]],
        plan_item = plan[[i]]
      )
    }
  }
  
  invalid_existing <-
    !force &
    destination_exists &
    !reusable
  
  pending <- if (force) {
    rep(
      TRUE,
      length(plan)
    )
  } else {
    !reusable
  }
  
  n_requested <- length(plan)
  n_reused <- sum(reusable)
  n_pending <- sum(pending)
  n_invalid <- sum(invalid_existing)
  
  # Estimate the amount of data that still needs to be transferred.
  #
  # The official Toolbox provides `data_transfer_size` during the dry-run.
  # This is only used for user information; it does not affect the request.
  pending_transfer_size <- 0
  
  if (n_pending > 0L) {
    
    pending_indices <- which(pending)
    
    pending_transfer_values <- vapply(
      pending_indices,
      function(i) {
        
        value <- plan[[i]]$data_transfer_size
        
        if (is.null(value) || !length(value)) {
          return(NA_character_)
        }
        
        as.character(value[[1]])
      },
      character(1)
    )
    
    parse_transfer_mb <- function(x) {
      
      if (is.na(x) || !nzchar(x)) {
        return(NA_real_)
      }
      
      number <- suppressWarnings(
        as.numeric(
          sub(
            "^\\s*([0-9.]+).*$",
            "\\1",
            x
          )
        )
      )
      
      if (!is.finite(number)) {
        return(NA_real_)
      }
      
      unit <- toupper(
        trimws(
          sub(
            "^\\s*[0-9.]+\\s*",
            "",
            x
          )
        )
      )
      
      multiplier <- switch(
        unit,
        "B" = 1 / (1024^2),
        "KB" = 1 / 1024,
        "MB" = 1,
        "GB" = 1024,
        "TB" = 1024^2,
        NA_real_
      )
      
      number * multiplier
    }
    
    pending_transfer_mb <- vapply(
      pending_transfer_values,
      parse_transfer_mb,
      numeric(1)
    )
    
    if (all(is.finite(pending_transfer_mb))) {
      pending_transfer_size <- sum(pending_transfer_mb)
    } else {
      pending_transfer_size <- NA_real_
    }
  }
  
  format_transfer_size <- function(size_mb) {
    
    if (!is.finite(size_mb)) {
      return(NULL)
    }
    
    if (size_mb >= 1024^2) {
      return(
        paste0(
          format(round(size_mb / 1024^2, 2), nsmall = 2),
          " TB"
        )
      )
    }
    
    if (size_mb >= 1024) {
      return(
        paste0(
          format(round(size_mb / 1024, 2), nsmall = 2),
          " GB"
        )
      )
    }
    
    paste0(
      format(round(size_mb, 2), nsmall = 2),
      " MB"
    )
  }
  
  if (!quiet) {
    
    message(
      "  Planned:    ",
      n_requested,
      if (n_requested == 1L) " file" else " files"
    )
    
    if (!force) {
      message(
        "  Reusable:   ",
        n_reused
      )
      
      if (n_invalid > 0L) {
        message(
          "  Invalid:    ",
          n_invalid,
          " existing ",
          if (n_invalid == 1L) "file" else "files",
          " will be replaced"
        )
      }
      
    } else {
      message(
        "  Force:      all planned files will be downloaded again"
      )
    }
    
    if (n_pending > 0L) {
      
      message(
        "  To download: ",
        n_pending
      )
      
      transfer_label <- format_transfer_size(
        pending_transfer_size
      )
      
      if (!is.null(transfer_label)) {
        message(
          "  Transfer:   ~",
          transfer_label
        )
      }
    }
  }
  
  # -----------------------------------------------------------------------
  # If everything is already available and valid, return immediately without
  # downloading any data.
  # -----------------------------------------------------------------------
  
  if (!any(pending)) {
    
    output_files <- normalizePath(
      planned_destinations,
      winslash = "/",
      mustWork = TRUE
    )
    
    completed <- TRUE
    
    if (!quiet) {
      
      elapsed <- as.numeric(
        difftime(
          Sys.time(),
          start_time,
          units = "secs"
        )
      )
      
      message(
        "\n[2/4] Downloading data from Copernicus Marine..."
      )
      
      message(
        "  Nothing to download; all planned files are already available ",
        "and valid."
      )
      
      message(
        "[3/4] Organizing files..."
      )
      
      message(
        "  Nothing to organize."
      )
      
      message(
        "[4/4] Validating output..."
      )
      
      message(
        "\nDone in ",
        format(
          round(
            elapsed,
            1
          ),
          nsmall = 1
        ),
        " s.\n",
        "Requested:  ",
        n_requested,
        "\n",
        "Reused:     ",
        n_reused,
        "\n",
        "Downloaded: 0"
      )
    }
    
    return(
      sort(
        unique(
          output_files
        )
      )
    )
  }
  
  if (!quiet) {
    message(
      "[2/4] Downloading ",
      n_pending,
      if (n_pending == 1L) {
        " missing file..."
      } else {
        " missing files..."
      }
    )
  }
  
  # -----------------------------------------------------------------------
  # Download only the missing or invalid files.
  #
  # Each dry-run plan item already describes exactly one environmental
  # variable and one temporal output file. We therefore submit one small
  # Toolbox request for each pending plan item.
  #
  # The original spatial/depth request is retained. The time range comes from
  # the authoritative dry-run plan, and the planned filename is supplied
  # explicitly so the final filename remains identical to the Toolbox plan.
  # -----------------------------------------------------------------------
  
  downloaded_indices <- which(
    pending
  )
  
  # -----------------------------------------------------------------------
  # Bulk path when every planned file requires download.
  #
  # This branch is reached after pre-planning when existing NetCDF files were
  # found below `output_dir`, but none of the files required by the current
  # request can be safely reused.
  #
  # When every planned file is pending, the official Toolbox can perform the
  # complete request in a single process and split the result internally by
  # variable and/or time. This is substantially faster than launching one
  # Toolbox process for every output file.
  #
  # Partial resumes continue to use the individual-file path below so valid
  # existing files are never downloaded again.
  #
  # Custom week/season splitting remains package-managed because these split
  # units are not provided natively by the Toolbox.
  # -----------------------------------------------------------------------
  
  use_bulk_download <-
    !custom_split &&
    length(downloaded_indices) == length(plan)
  
  if (use_bulk_download) {
    
    if (!quiet) {
      message(
        "\n  All planned files require download.",
        "\n  Using one bulk Toolbox request for faster transfer..."
      )
    }
    
    bulk_args <- build_args(
      request_output_dir = staging_dir,
      request_start = start_datetime,
      request_end = end_datetime,
      request_filename = filename,
      request_variables = variables,
      dry_run = FALSE
    )
    
    needs_split_on <-
      multiple_variables ||
      native_time_split
    
    if (needs_split_on) {
      
      bulk_args <- c(
        bulk_args,
        "split-on"
      )
      
      if (multiple_variables) {
        bulk_args <- c(
          bulk_args,
          "--on-variables"
        )
      }
      
      if (native_time_split) {
        bulk_args <- c(
          bulk_args,
          "--on-time",
          split_by
        )
      }
      
      if (!is.null(concurrent_processes)) {
        bulk_args <- c(
          bulk_args,
          "--concurrent-processes",
          as.character(
            as.integer(
              concurrent_processes
            )
          )
        )
      }
    }
    
    bulk_result <- .copernicus_run_standalone(
      executable = executable,
      args = bulk_args,
      dataset_id = dataset_id
    )
    
    if (!identical(bulk_result$status, 0L)) {
      
      toolbox_output <- paste(
        bulk_result$output,
        collapse = "\n"
      )
      
      stop(
        "The Copernicus Marine bulk download failed.\n\n",
        "Dataset:\n  ",
        dataset_id,
        "\n\n",
        "The official Toolbox returned exit status ",
        bulk_result$status,
        ".\n\n",
        if (nzchar(toolbox_output)) {
          paste0(
            "Toolbox output:\n",
            toolbox_output
          )
        } else {
          ""
        },
        call. = FALSE
      )
    }
    
    if (!quiet) {
      
      toolbox_warnings <- grep(
        "WARNING",
        bulk_result$output,
        value = TRUE
      )
      
      if (length(toolbox_warnings)) {
        
        toolbox_warnings <- unique(
          trimws(
            sub(
              "^WARNING\\s*-\\s*\\d{4}-\\d{2}-\\d{2}T[^ ]+\\s*-\\s*",
              "",
              toolbox_warnings
            )
          )
        )
        
        for (warning_message in toolbox_warnings) {
          message(
            "        WARNING: ",
            warning_message
          )
        }
      }
    }
    
    downloaded_staging_paths <- character(
      length(downloaded_indices)
    )
    
    staging_candidates <- list.files(
      staging_dir,
      pattern = "\\.(nc|nc4)$",
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )
    
    if (length(staging_candidates) != length(downloaded_indices)) {
      stop(
        "The Copernicus Marine Toolbox completed the bulk request, but the ",
        "number of downloaded NetCDF files does not match the authoritative ",
        "download plan.\n\n",
        "Expected:\n  ",
        length(downloaded_indices),
        " files\n\n",
        "Found:\n  ",
        length(staging_candidates),
        " files\n\n",
        "The files have been left in the temporary staging directory for ",
        "inspection.",
        call. = FALSE
      )
    }
    
    # Match every downloaded file back to its authoritative dry-run plan.
    #
    # We deliberately validate metadata rather than relying only on filenames.
    # This makes sure each downloaded file contains the expected variable,
    # spatial extent, depth range and time range.
    
    unmatched_candidates <- staging_candidates
    
    for (j in seq_along(downloaded_indices)) {
      
      i <- downloaded_indices[[j]]
      item <- plan[[i]]
      
      matching_candidate <- NULL
      
      # The Toolbox normally uses exactly the filename reported by dry-run.
      # Try that first because it is cheap and deterministic.
      
      expected_name <- basename(
        item$filename
      )
      
      filename_matches <- unmatched_candidates[
        basename(unmatched_candidates) == expected_name
      ]
      
      if (length(filename_matches) == 1L) {
        
        candidate_valid <- .copernicus_existing_file_matches(
          path = filename_matches[[1]],
          plan_item = item
        )
        
        if (candidate_valid) {
          matching_candidate <- filename_matches[[1]]
        }
      }
      
      # If the filename cannot identify the file uniquely, fall back to
      # metadata validation against the remaining downloaded files.
      
      if (is.null(matching_candidate)) {
        
        candidate_matches <- vapply(
          unmatched_candidates,
          function(candidate) {
            .copernicus_existing_file_matches(
              path = candidate,
              plan_item = item
            )
          },
          logical(1)
        )
        
        if (sum(candidate_matches) != 1L) {
          stop(
            "A file from the Copernicus Marine bulk download could not be ",
            "matched uniquely to the authoritative dry-run plan.\n\n",
            "Planned file:\n  ",
            item$filename,
            "\n\n",
            "Matching downloaded files:\n  ",
            sum(candidate_matches),
            "\n\n",
            "The bulk download has been left in the temporary staging ",
            "directory for inspection.",
            call. = FALSE
          )
        }
        
        matching_candidate <- unmatched_candidates[
          candidate_matches
        ][[1]]
      }
      
      downloaded_staging_paths[[j]] <- matching_candidate
      
      unmatched_candidates <- setdiff(
        unmatched_candidates,
        matching_candidate
      )
    }
    
  } else {
    
    downloaded_staging_paths <- character(
      length(downloaded_indices)
    )
  
  for (j in seq_along(downloaded_indices)) {
    
    i <- downloaded_indices[[j]]
    item <- plan[[i]]
    
    item_variable <- as.character(
      item$variables[[1]]
    )
    
    item_time <- .copernicus_plan_extent(
      item,
      "time"
    )
    
    if (is.null(item_time)) {
      stop(
        "The Copernicus Marine download plan does not contain a time extent ",
        "for:\n  ",
        item$filename,
        call. = FALSE
      )
    }
    
    item_start <- .copernicus_datetime(
      .copernicus_as_posixct(
        as.character(
          item_time$minimum
        )
      ),
      "planned start datetime"
    )
    
    item_end <- .copernicus_datetime(
      .copernicus_as_posixct(
        as.character(
          item_time$maximum
        )
      ),
      "planned end datetime"
    )
    
    item_dir <- file.path(
      staging_dir,
      sprintf(
        "file_%06d",
        i
      )
    )
    
    dir.create(
      item_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    if (!quiet) {
      
      item_transfer <- plan[[i]]$data_transfer_size
      
      if (is.null(item_transfer) || !length(item_transfer)) {
        item_transfer <- NULL
      } else {
        item_transfer <- as.character(item_transfer[[1]])
      }
      
      message(
        "\n  [",
        j,
        "/",
        length(downloaded_indices),
        "] ",
        item_variable,
        "\n",
        "        Time:     ",
        item_start,
        if (!identical(item_start, item_end)) {
          paste0(
            " to ",
            item_end
          )
        } else {
          ""
        },
        if (!is.null(item_transfer)) {
          paste0(
            "\n        Transfer: ~",
            item_transfer
          )
        } else {
          ""
        }
      )
    }
    
    item_args <- build_args(
      request_output_dir = item_dir,
      request_start = item_start,
      request_end = item_end,
      request_filename = basename(
        item$filename
      ),
      request_variables = item_variable,
      dry_run = FALSE
    )
    
    item_result <- .copernicus_run_standalone(
      executable = executable,
      args = item_args,
      dataset_id = dataset_id
    )
    
    if (!identical(item_result$status, 0L)) {
      
      toolbox_output <- paste(
        item_result$output,
        collapse = "\n"
      )
      
      stop(
        "The Copernicus Marine request failed while downloading file ",
        j,
        " of ",
        length(downloaded_indices),
        ".\n\n",
        "Variable:\n  ",
        item_variable,
        "\n\n",
        "Planned file:\n  ",
        item$filename,
        "\n\n",
        "The official Toolbox returned exit status ",
        item_result$status,
        ".\n\n",
        if (nzchar(toolbox_output)) {
          paste0(
            "Toolbox output:\n",
            toolbox_output
          )
        } else {
          ""
        },
        call. = FALSE
      )
    }
    
    if (!quiet) {
      
      toolbox_warnings <- grep(
        "WARNING",
        item_result$output,
        value = TRUE
      )
      
      if (length(toolbox_warnings)) {
        
        toolbox_warnings <- unique(
          trimws(
            sub(
              "^WARNING\\s*-\\s*\\d{4}-\\d{2}-\\d{2}T[^ ]+\\s*-\\s*",
              "",
              toolbox_warnings
            )
          )
        )
        
        for (warning_message in toolbox_warnings) {
          message(
            "        WARNING: ",
            warning_message
          )
        }
      }
    }
    
    expected_staging_path <- file.path(
      item_dir,
      basename(
        item$filename
      )
    )
    
    if (!file.exists(expected_staging_path)) {
      
      candidates <- list.files(
        item_dir,
        pattern = "\\.(nc|nc4)$",
        full.names = TRUE,
        recursive = TRUE,
        ignore.case = TRUE
      )
      
      if (length(candidates) != 1L) {
        stop(
          "The Copernicus Marine Toolbox reported a successful download, ",
          "but the expected NetCDF file could not be identified.\n\n",
          "Expected:\n  ",
          expected_staging_path,
          "\n\n",
          "NetCDF files found:\n  ",
          if (length(candidates)) {
            paste(
              candidates,
              collapse = "\n  "
            )
          } else {
            "none"
          },
          call. = FALSE
        )
      }
      
      expected_staging_path <- candidates[[1]]
    }
    
    staging_valid <- .copernicus_existing_file_matches(
      path = expected_staging_path,
      plan_item = item
    )
    
    if (!staging_valid) {
      stop(
        "A Copernicus Marine file was downloaded successfully but does not ",
        "match the authoritative dry-run plan.\n\n",
        "File:\n  ",
        expected_staging_path,
        "\n\n",
        "Variable:\n  ",
        item_variable,
        "\n\n",
        "The file has been left in the temporary staging directory for ",
        "inspection rather than being accepted as valid output.",
        call. = FALSE
      )
    }
    
    downloaded_staging_paths[[j]] <- expected_staging_path
  }
    
  }
  
  # -----------------------------------------------------------------------
  # Organize newly downloaded files.
  #
  # Only pending files reach this stage. Therefore an existing destination
  # encountered here is either invalid or is intentionally being replaced
  # because `force = TRUE`.
  # -----------------------------------------------------------------------
  
  if (!quiet) {
    message(
      "[3/4] Organizing ",
      length(downloaded_staging_paths),
      if (length(downloaded_staging_paths) == 1L) {
        " downloaded file..."
      } else {
        " downloaded files..."
      }
    )
  }
  
  organized_downloads <- character(
    length(downloaded_staging_paths)
  )
  
  for (j in seq_along(downloaded_staging_paths)) {
    
    i <- downloaded_indices[[j]]
    
    organized <- .copernicus_organize_files(
      paths = downloaded_staging_paths[[j]],
      output_dir = output_dir,
      source = "marine",
      dataset_id = dataset_id,
      requested_variables = variables,
      organize_by = organize_by,
      hemisphere = hemisphere,
      force = TRUE,
      quiet = TRUE
    )
    
    expected_destination <- planned_destinations[[i]]
    
    organized_normalized <- normalizePath(
      organized,
      winslash = "/",
      mustWork = TRUE
    )
    
    expected_normalized <- normalizePath(
      expected_destination,
      winslash = "/",
      mustWork = TRUE
    )
    
    if (!identical(
      organized_normalized,
      expected_normalized
    )) {
      stop(
        "Internal Copernicus file-organization mismatch.\n\n",
        "The dry-run plan expected:\n  ",
        expected_normalized,
        "\n\n",
        "but the downloaded file was organized as:\n  ",
        organized_normalized,
        call. = FALSE
      )
    }
    
    organized_downloads[[j]] <- organized_normalized
    
    if (!quiet) {
      
      relative_path <- sub(
        paste0(
          "^",
          gsub(
            "([][{}()+*^$|\\\\?.])",
            "\\\\\\1",
            normalizePath(
              output_dir,
              winslash = "/",
              mustWork = TRUE
            )
          ),
          "/?"
        ),
        "",
        organized_normalized
      )
      
      message(
        "  [",
        j,
        "/",
        length(downloaded_staging_paths),
        "] ",
        relative_path
      )
    }
  }
  
  # -----------------------------------------------------------------------
  # Final validation.
  #
  # Every planned output is checked again using metadata only. This catches
  # accidental file moves, invalid replacements, or other inconsistencies
  # without reading the environmental data arrays themselves.
  # -----------------------------------------------------------------------
  
  if (!quiet) {
    message(
      "[4/4] Validating ",
      n_requested,
      if (n_requested == 1L) {
        " output file..."
      } else {
        " output files..."
      }
    )
  }
  
  final_valid <- logical(
    length(plan)
  )
  
  for (i in seq_along(plan)) {
    
    final_valid[[i]] <- .copernicus_existing_file_matches(
      path = planned_destinations[[i]],
      plan_item = plan[[i]]
    )
  }
  
  if (!all(final_valid)) {
    
    failed <- planned_destinations[
      !final_valid
    ]
    
    stop(
      "One or more final Copernicus Marine files failed validation:\n  ",
      paste(
        failed,
        collapse = "\n  "
      ),
      call. = FALSE
    )
  }
  
  output_files <- normalizePath(
    planned_destinations,
    winslash = "/",
    mustWork = TRUE
  )
  
  completed <- TRUE
  
  if (!quiet) {
    
    elapsed <- as.numeric(
      difftime(
        Sys.time(),
        start_time,
        units = "secs"
      )
    )
    
    message(
      "\nDone in ",
      format(
        round(
          elapsed,
          1
        ),
        nsmall = 1
      ),
      " s.\n",
      "Requested:  ",
      n_requested,
      "\n",
      "Reused:     ",
      n_reused,
      "\n",
      "Downloaded: ",
      n_pending,
      if (n_invalid > 0L && !force) {
        paste0(
          "\nReplaced:   ",
          n_invalid,
          " invalid existing ",
          if (n_invalid == 1L) "file" else "files"
        )
      } else {
        ""
      }
    )
  }
  
  sort(
    unique(
      output_files
    )
  )
}

# Internal: identify CERRA datasets requiring local geographic cropping.
.copernicus_is_cerra <- function(dataset_id) {
  grepl("^reanalysis-cerra-", dataset_id, ignore.case = TRUE)
}

# Internal: crop CERRA locally while preserving its native Lambert grid.
.copernicus_crop_cerra <- function(input_file, output_file,
                                   xmin, xmax, ymin, ymax,
                                   overwrite = FALSE) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required to crop CERRA downloads locally.", call. = FALSE)
  }

  r <- tryCatch(
    terra::rast(input_file),
    error = function(e) stop(
      "Could not open the downloaded CERRA file with terra: ",
      conditionMessage(e), call. = FALSE
    )
  )

  bbox <- terra::vect(
    matrix(c(xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin),
           ncol = 2, byrow = TRUE),
    type = "polygons",
    crs = "EPSG:4326"
  )

  bbox <- tryCatch(
    terra::project(bbox, terra::crs(r)),
    error = function(e) stop(
      "Could not project the requested bounding box to the CERRA CRS: ",
      conditionMessage(e), call. = FALSE
    )
  )

  cropped <- tryCatch(
    terra::crop(r, bbox),
    error = function(e) stop(
      "Could not crop the CERRA raster to the requested bounding box: ",
      conditionMessage(e), call. = FALSE
    )
  )

  if (terra::ncell(cropped) == 0) {
    stop("The requested bounding box does not overlap the CERRA domain.", call. = FALSE)
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  ext <- tolower(tools::file_ext(output_file))
  filetype <- switch(
    ext,
    grib = "GRIB",
    grb = "GRIB",
    grib2 = "GRIB",
    grb2 = "GRIB",
    nc = "netCDF",
    nc4 = "netCDF",
    tif = "GTiff",
    tiff = "GTiff",
    NULL
  )

  if (is.null(filetype)) {
    stop(
      "Cannot determine an output raster format from `output_file`: ",
      output_file,
      ". Use a supported extension such as .grib, .nc, or .tif.",
      call. = FALSE
    )
  }

  terra::writeRaster(
    cropped,
    output_file,
    overwrite = overwrite,
    filetype = filetype
  )

  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

# Internal: validate/configure ECMWF Data Store authentication.
.copernicus_ecmwf_auth <- function() {
  key <- tryCatch(ecmwfr::wf_get_key(), error = function(e) NULL)
  
  if (is.character(key) && length(key) == 1 && !is.na(key) && nzchar(key)) {
    return(invisible(TRUE))
  }
  
  if (!interactive()) {
    stop(
      "ECMWF Data Store authentication is required. `ecmwfr >= 2.0.0` uses ",
      "a Personal Access Token (PAT), not the legacy username/password scheme. ",
      "Run `ecmwfr::wf_set_key()` once before using CDS or ADS ",
      "non-interactively.",
      call. = FALSE
    )
  }
  
  message(
    "ECMWF Data Store authentication is required.\n",
    "CDS and ADS use the same Personal Access Token (PAT).\n",
    "Opening the official CDS API setup page. Log in, copy your PAT, then paste ",
    "it into the authentication dialog.\n",
    "Legacy ERA5/CERRA username/password files are not used."
  )
  
  try(
    utils::browseURL("https://cds.climate.copernicus.eu/how-to-api"),
    silent = TRUE
  )
  
  tryCatch(
    ecmwfr::wf_set_key(),
    error = function(e) {
      stop(
        "Could not configure the ECMWF Data Store token: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  key <- tryCatch(ecmwfr::wf_get_key(), error = function(e) NULL)
  
  if (!is.character(key) || length(key) != 1 || is.na(key) || !nzchar(key)) {
    stop(
      "ECMWF Data Store authentication was not configured successfully. ",
      "Obtain a Personal Access Token from the official CDS API setup page ",
      "and run `ecmwfr::wf_set_key()`.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

# Internal: CDS/ADS backend using ecmwfr.
.copernicus_load_ecmwf <- function(source, dataset_id, variables,
                                   start_datetime, end_datetime,
                                   xmin, xmax, ymin, ymax,
                                   depth_min, depth_max,
                                   request, output_dir, filename,
                                   split_by, concurrent_processes,
                                   compression, dataset_version,
                                   dataset_part, force, quiet) {

  if (!requireNamespace("ecmwfr", quietly = TRUE)) {
    stop(
      "Package 'ecmwfr' is required for `source = \"", source, "\"`. ",
      "Install it with install.packages('ecmwfr').",
      call. = FALSE
    )
  }

  .copernicus_ecmwf_auth()

  unsupported <- c(
    split_by = !is.null(split_by),
    concurrent_processes = !is.null(concurrent_processes),
    compression = !is.null(compression),
    dataset_version = !is.null(dataset_version),
    dataset_part = !is.null(dataset_part)
  )

  if (any(unsupported)) {
    stop(
      "The following arguments are specific to Copernicus Marine and cannot ",
      "be used with `source = \"", source, "\"`: ",
      paste(names(unsupported)[unsupported], collapse = ", "),
      ". Supply dataset-specific ECMWF parameters through `request` instead.",
      call. = FALSE
    )
  }

  if (!is.null(depth_min) || !is.null(depth_max)) {
    stop(
      "`depth_min` and `depth_max` are currently supported directly only for ",
      "Copernicus Marine. For ", toupper(source),
      " datasets, specify vertical dimensions such as `pressure_level`, ",
      "`model_level`, or equivalent through `request`.",
      call. = FALSE
    )
  }

  req <- request
  req$dataset_short_name <- dataset_id

  # Explicit common arguments override duplicate request entries.
  if (!is.null(variables)) req$variable <- variables

  if (!is.null(start_datetime)) {
    start_date <- substr(start_datetime, 1, 10)
    end_date <- substr(end_datetime, 1, 10)
    start_time <- substr(start_datetime, 12, 16)
    end_time <- substr(end_datetime, 12, 16)

    req$date <- if (identical(start_date, end_date)) {
      start_date
    } else {
      paste(start_date, end_date, sep = "/")
    }

    if (!identical(start_time, end_time)) {
      stop(
        "`start_datetime` and `end_datetime` contain different clock times. ",
        "For CDS/ADS requests requiring multiple times, supply the time ",
        "selection explicitly through `request = list(time = ...)` and omit ",
        "`start_datetime`/`end_datetime`.",
        call. = FALSE
      )
    }

    req$time <- start_time
  }

  is_cerra <- .copernicus_is_cerra(dataset_id)
  local_crop <- is_cerra && !is.null(xmin)

  if (!is.null(xmin) && !local_crop) req$area <- c(ymax, xmin, ymin, xmax)
  if (!local_crop && !is.null(filename)) req$target <- filename

  if (local_crop && !quiet) {
    message(
      "CERRA dataset detected: downloading the native full domain temporarily ",
      "and cropping locally to the requested bounding box."
    )
  }

  if (!is.null(req$target)) {
    .copernicus_validate_string(req$target, "request$target")
    dest <- file.path(output_dir, req$target)

    if (file.exists(dest) && !force) {
      if (!quiet) message("Cached: ", req$target)
      return(normalizePath(dest, winslash = "/", mustWork = TRUE))
    }

    if (file.exists(dest) && force) unlink(dest, recursive = TRUE, force = TRUE)
  }

  if (!quiet) message("Requesting ", toupper(source), " dataset ", dataset_id, "...")

  if (local_crop) {
    if (is.null(filename)) {
      ext <- if (!is.null(req$data_format) &&
                 tolower(as.character(req$data_format)[1]) %in% c("netcdf", "netcdf4")) ".nc" else ".grib"
      var_tag <- if (!is.null(variables) && length(variables)) {
        gsub("[^A-Za-z0-9._-]+", "_", variables[1])
      } else "cerra"
      filename <- paste0(var_tag, "_cropped", ext)
    }

    .copernicus_validate_string(filename, "filename")
    dest <- file.path(output_dir, filename)

    if (file.exists(dest) && !force) {
      if (!quiet) message("Cached: ", filename)
      return(normalizePath(dest, winslash = "/", mustWork = TRUE))
    }
    if (file.exists(dest) && force) unlink(dest, recursive = TRUE, force = TRUE)

    tmp_dir <- tempfile(pattern = "copernicus_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

    tmp_ext <- if (!is.null(req$data_format) &&
                  tolower(as.character(req$data_format)[1]) %in% c("netcdf", "netcdf4")) ".nc" else ".grib"
    req$target <- paste0("cerra_full_", format(Sys.time(), "%Y%m%d%H%M%S"), tmp_ext)

    full_path <- tryCatch(
      ecmwfr::wf_request(
        request = req,
        transfer = TRUE,
        path = tmp_dir,
        verbose = !quiet
      ),
      error = function(e) stop(
        toupper(source), " request failed: ", conditionMessage(e), call. = FALSE
      )
    )

    if (!is.character(full_path) || length(full_path) != 1 || !file.exists(full_path)) {
      stop("CERRA request completed but did not return one readable local file.", call. = FALSE)
    }

    cropped <- .copernicus_crop_cerra(
      input_file = full_path,
      output_file = dest,
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      overwrite = force
    )

    if (!quiet) {
      message(sprintf(
        "CERRA local crop complete: %.2f MB full domain -> %.2f MB cropped.",
        file.info(full_path)$size / 1024^2,
        file.info(cropped)$size / 1024^2
      ))
    }

    return(cropped)
  }

  path <- tryCatch(
    ecmwfr::wf_request(
      request = req,
      transfer = TRUE,
      path = output_dir,
      verbose = !quiet
    ),
    error = function(e) {
      stop(
        toupper(source), " request failed: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (!is.character(path) || !length(path)) {
    stop(
      toupper(source),
      " request was submitted but no completed local file path was returned.",
      call. = FALSE
    )
  }

  missing <- !file.exists(path)
  if (any(missing)) {
    stop(
      toupper(source), " returned output path(s) that do not exist:\n  ",
      paste(path[missing], collapse = "\n  "),
      call. = FALSE
    )
  }

  normalizePath(path, winslash = "/", mustWork = TRUE)
}


# Internal: validate netCDF files supplied to copernicus_summarise().
.copernicus_summary_files <- function(file_paths) {
  if (!is.character(file_paths) || !length(file_paths) || anyNA(file_paths) ||
      any(!nzchar(file_paths))) {
    stop("`file_paths` must be a non-empty character vector.", call. = FALSE)
  }

  missing <- !file.exists(file_paths)
  if (any(missing)) {
    stop(
      "Copernicus netCDF file(s) not found:\n  ",
      paste(file_paths[missing], collapse = "\n  "),
      call. = FALSE
    )
  }

  ext <- tolower(tools::file_ext(file_paths))
  if (any(!ext %in% c("nc", "nc4"))) {
    stop("`copernicus_summarise()` currently supports netCDF (.nc/.nc4) files only.",
         call. = FALSE)
  }

  normalizePath(file_paths, winslash = "/", mustWork = TRUE)
}

# Internal: identify the temporal dimension in a netCDF file.
.copernicus_nc_time_dim <- function(nc) {
  dims <- names(nc$dim)
  lower <- tolower(dims)

  exact <- which(lower %in% c("time", "valid_time"))
  if (length(exact) == 1) return(dims[exact])

  by_units <- which(vapply(nc$dim, function(x) {
    units <- if (is.null(x$units)) "" else tolower(x$units)
    grepl("^(seconds|minutes|hours|days|months|years) since ", units)
  }, logical(1)))

  candidates <- unique(c(exact, by_units))
  if (length(candidates) == 1) return(dims[candidates])

  if (!length(candidates)) {
    stop("Could not identify a time dimension in the netCDF file.", call. = FALSE)
  }

  stop(
    "Multiple possible time dimensions were found: ",
    paste(dims[candidates], collapse = ", "),
    ".",
    call. = FALSE
  )
}

# Internal: convert a netCDF time dimension to UTC POSIXct values.
.copernicus_nc_time_values <- function(nc, time_dim) {
  d <- nc$dim[[time_dim]]
  vals <- as.numeric(d$vals)
  units <- if (is.null(d$units)) "" else as.character(d$units)

  m <- regexec(
    "^(seconds?|minutes?|hours?|days?) since (.+)$",
    units,
    ignore.case = TRUE
  )
  parts <- regmatches(units, m)[[1]]

  if (length(parts) != 3) {
    stop(
      "Unsupported or missing time units for dimension '", time_dim, "': ",
      units, ". Expected units such as 'hours since YYYY-MM-DD HH:MM:SS'.",
      call. = FALSE
    )
  }

  origin_txt <- sub("Z$", "", parts[3])
  origin <- tryCatch(
    as.POSIXct(
      origin_txt,
      tz = "UTC",
      tryFormats = c(
        "%Y-%m-%d %H:%M:%OS",
        "%Y-%m-%dT%H:%M:%OS",
        "%Y-%m-%d"
      )
    ),
    error = function(e) as.POSIXct(NA, tz = "UTC")
  )

  if (is.na(origin)) {
    stop("Could not parse netCDF time origin: ", parts[3], call. = FALSE)
  }

  mult <- switch(
    tolower(parts[2]),
    second = 1, seconds = 1,
    minute = 60, minutes = 60,
    hour = 3600, hours = 3600,
    day = 86400, days = 86400
  )

  origin + vals * mult
}

# Internal: index temporal slices across one or more netCDF files.
.copernicus_summary_time_index <- function(file_paths, start_datetime, end_datetime) {
  rows <- vector("list", length(file_paths))

  for (i in seq_along(file_paths)) {
    nc <- ncdf4::nc_open(file_paths[i])
    rows[[i]] <- tryCatch(
      {
        time_dim <- .copernicus_nc_time_dim(nc)
        time <- .copernicus_nc_time_values(nc, time_dim)

        data.frame(
          file = file_paths[i],
          file_index = i,
          time_index = seq_along(time),
          datetime = as.POSIXct(time, origin = "1970-01-01", tz = "UTC"),
          stringsAsFactors = FALSE
        )
      },
      finally = ncdf4::nc_close(nc)
    )
  }

  idx <- do.call(rbind, rows)
  idx <- idx[order(idx$datetime, idx$file_index, idx$time_index), , drop = FALSE]
  rownames(idx) <- NULL

  start <- if (is.null(start_datetime)) NULL else
    as.POSIXct(.copernicus_datetime(start_datetime, "start_datetime"),
               format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
  end <- if (is.null(end_datetime)) NULL else
    as.POSIXct(.copernicus_datetime(end_datetime, "end_datetime"),
               format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")

  if (!is.null(start) && is.null(end)) end <- start
  if (is.null(start) && !is.null(end)) {
    stop("`start_datetime` is required when `end_datetime` is supplied.",
         call. = FALSE)
  }
  if (!is.null(start) && end < start) {
    stop("`end_datetime` must be equal to or later than `start_datetime`.",
         call. = FALSE)
  }

  selected <- rep(TRUE, nrow(idx))
  if (!is.null(start)) selected <- selected & idx$datetime >= start
  if (!is.null(end)) selected <- selected & idx$datetime <= end

  idx$selected <- selected
  idx
}

# Internal: print available and selected time information.
.copernicus_summary_time_message <- function(index) {
  selected <- index[index$selected, , drop = FALSE]

  message(
    "Available time steps: ", nrow(index),
    " (", format(min(index$datetime), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    " to ", format(max(index$datetime), "%Y-%m-%d %H:%M:%S", tz = "UTC"), ")"
  )

  if (!nrow(selected)) {
    message("Selected time steps: 0")
    return(invisible(NULL))
  }

  message(
    "Selected time steps: ", nrow(selected),
    " (", format(min(selected$datetime), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    " to ", format(max(selected$datetime), "%Y-%m-%d %H:%M:%S", tz = "UTC"), ")"
  )
  message("Dates/times to summarise:")
  for (x in format(selected$datetime, "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")) {
    message("  - ", x)
  }

  invisible(NULL)
}

# Internal: return the names of variables that should be summarised over time.
.copernicus_nc_summary_vars <- function(nc, time_dim) {
  vars <- names(nc$var)

  has_time <- vapply(vars, function(v) {
    time_dim %in% vapply(nc$var[[v]]$dim, `[[`, character(1), "name")
  }, logical(1))

  vars <- vars[has_time]
  vars <- vars[!grepl("(^|_)(time_?bounds?|bounds?|bnds?)$", vars, ignore.case = TRUE)]

  vars <- vars[vapply(vars, function(v) {
    dims <- vapply(nc$var[[v]]$dim, `[[`, character(1), "name")
    length(setdiff(dims, time_dim)) > 0
  }, logical(1))]

  if (!length(vars)) {
    stop("No data variables with a time dimension were found.", call. = FALSE)
  }

  vars
}

# Internal: create output paths for requested summaries.
.copernicus_summary_paths <- function(file_paths, fun, output_dir, filename) {
  if (is.null(output_dir)) output_dir <- dirname(file_paths[1])
  .copernicus_validate_string(output_dir, "output_dir")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(filename)) {
    stem <- tools::file_path_sans_ext(basename(file_paths[1]))
    out <- file.path(output_dir, paste0(stem, "_", fun, ".nc"))
  } else {
    .copernicus_validate_string(filename, "filename")
    ext <- tools::file_ext(filename)
    if (!tolower(ext) %in% c("nc", "nc4")) {
      stop("`filename` must end in .nc or .nc4.", call. = FALSE)
    }

    stem <- tools::file_path_sans_ext(filename)
    suffix <- paste0(".", ext)

    if (length(fun) == 1) {
      out <- file.path(output_dir, filename)
    } else {
      out <- file.path(output_dir, paste0(stem, "_", fun, suffix))
    }
  }

  stats::setNames(out, fun)
}

# Internal: copy netCDF attributes, excluding attributes handled separately.
.copernicus_nc_copy_atts <- function(nc_in, nc_out, var_in = 0, var_out = 0,
                                     exclude = character()) {
  atts <- tryCatch(ncdf4::ncatt_get(nc_in, var_in), error = function(e) list())
  if (!length(atts)) return(invisible(TRUE))

  keep <- setdiff(names(atts), exclude)
  for (nm in keep) {
    value <- atts[[nm]]
    if (is.null(value) || !length(value)) next
    try(ncdf4::ncatt_put(nc_out, var_out, nm, value), silent = TRUE)
  }

  invisible(TRUE)
}

# Internal: build dimension definitions for a summarised netCDF.
.copernicus_nc_dim_defs <- function(nc, dim_names) {
  defs <- lapply(dim_names, function(nm) {
    d <- nc$dim[[nm]]
    calendar <- if (is.null(d$calendar) || !length(d$calendar)) NA else d$calendar
    longname <- if (is.null(d$longname) || !length(d$longname)) nm else d$longname

    ncdf4::ncdim_def(
      name = nm,
      units = if (is.null(d$units)) "" else d$units,
      vals = d$vals,
      unlim = FALSE,
      create_dimvar = TRUE,
      calendar = calendar,
      longname = longname
    )
  })

  stats::setNames(defs, dim_names)
}

# Internal: initialize streaming accumulators for one variable.
.copernicus_summary_accumulator <- function(dim_lengths, fun) {
  need_mean <- any(fun %in% c("mean", "sd"))

  list(
    dim = dim_lengths,
    count = array(0, dim = dim_lengths),
    mean = if (need_mean) array(0, dim = dim_lengths) else NULL,
    m2 = if ("sd" %in% fun) array(0, dim = dim_lengths) else NULL,
    min = if ("min" %in% fun) array(Inf, dim = dim_lengths) else NULL,
    max = if ("max" %in% fun) array(-Inf, dim = dim_lengths) else NULL,
    any_na = array(FALSE, dim = dim_lengths)
  )
}

# Internal: update one variable accumulator from one time slice.
.copernicus_summary_update <- function(acc, x, fun, na.rm) {
  x <- array(as.numeric(x), dim = acc$dim)
  valid <- !is.na(x)

  if (!na.rm) acc$any_na <- acc$any_na | !valid

  if (any(valid)) {
    old_count <- acc$count[valid]
    new_count <- old_count + 1

    if (any(fun %in% c("mean", "sd"))) {
      delta <- x[valid] - acc$mean[valid]
      acc$mean[valid] <- acc$mean[valid] + delta / new_count

      if ("sd" %in% fun) {
        delta2 <- x[valid] - acc$mean[valid]
        acc$m2[valid] <- acc$m2[valid] + delta * delta2
      }
    }

    if ("min" %in% fun) acc$min[valid] <- pmin(acc$min[valid], x[valid])
    if ("max" %in% fun) acc$max[valid] <- pmax(acc$max[valid], x[valid])
    acc$count[valid] <- new_count
  }

  acc
}

# Internal: finalize one summary statistic from an accumulator.
.copernicus_summary_finalize <- function(acc, fun, na.rm) {
  out <- switch(
    fun,
    mean = acc$mean,
    min = acc$min,
    max = acc$max,
    sd = {
      z <- array(NA_real_, dim = acc$dim)
      ok <- acc$count > 1
      z[ok] <- sqrt(acc$m2[ok] / (acc$count[ok] - 1))
      z
    }
  )

  if (fun %in% c("mean", "min", "max")) out[acc$count == 0] <- NA_real_
  if (!na.rm) out[acc$any_na] <- NA_real_
  out
}

#' Summarise Copernicus netCDF data across time
#'
#' Summarise one or more Copernicus netCDF files across their temporal
#' dimension while preserving all non-temporal dimensions (for example
#' longitude, latitude, depth, or projected x/y dimensions). Processing is
#' performed directly on the native netCDF data using `ncdf4`; conversion to a
#' `SpatRaster` is not required.
#'
#' All variables containing the detected time dimension are summarised. Variables
#' without a time dimension are copied unchanged from the first input file so
#' that auxiliary coordinates, grid mappings, and other static metadata remain
#' available in the output.
#'
#' The calculation is performed one time step at a time. This avoids loading the
#' complete multidimensional time series into memory and is suitable for large
#' Copernicus files. When multiple input files are supplied they must have
#' compatible non-temporal dimensions and variables.
#'
#' One netCDF file is written per requested summary statistic. The time dimension
#' is removed from the summarised variables; spatial and depth dimensions are
#' retained.
#'
#' @param file_paths Character vector. Paths to one or more Copernicus netCDF
#'   (`.nc` or `.nc4`) files.
#' @param fun Character vector of temporal summary statistics. Supported values
#'   are `"mean"`, `"min"`, `"max"`, and `"sd"`.
#' @param start_datetime Optional start datetime used to filter the available
#'   netCDF time steps before summarising. A `Date`, `POSIXt`, or character
#'   value interpreted in UTC.
#' @param end_datetime Optional end datetime. Defaults to `start_datetime` when
#'   only a start is supplied. If both are `NULL`, all available time steps are
#'   summarised.
#' @param output_dir Optional character. Directory in which summarised netCDF
#'   files are written. Defaults to the directory containing the first input
#'   file.
#' @param filename Optional character. Output filename. When multiple statistics
#'   are requested, the statistic is appended before the extension
#'   (for example `summary_mean.nc`, `summary_max.nc`).
#' @param na.rm Logical. If `TRUE` (default), missing values are ignored within
#'   each grid cell/depth combination. If `FALSE`, any missing temporal value
#'   produces a missing summary value at that location.
#' @param force Logical. Overwrite existing summary files. Default `FALSE`.
#' @param quiet Logical. Suppress progress messages. Default `FALSE`.
#'
#' @returns Named character vector containing paths to the summarised netCDF
#'   files, with names corresponding to `fun`.
#'
#' @examples
#' \dontrun{
#' summaries <- copernicus_summarise(
#'   file_paths = c("thetao_2020_01.nc", "thetao_2020_02.nc"),
#'   fun = c("mean", "min", "max", "sd"),
#'   output_dir = "summary"
#' )
#' }
#'
#' @export
copernicus_summarise <- function(file_paths,
                                 fun = c("mean", "min", "max", "sd"),
                                 start_datetime = NULL,
                                 end_datetime = NULL,
                                 output_dir = NULL,
                                 filename = NULL,
                                 na.rm = TRUE,
                                 force = FALSE,
                                 quiet = FALSE) {

  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop(
      "Package 'ncdf4' is required for `copernicus_summarise()`. ",
      "Install it with install.packages('ncdf4').",
      call. = FALSE
    )
  }

  file_paths <- .copernicus_summary_files(file_paths)

  if (!is.character(fun) || !length(fun) || anyNA(fun)) {
    stop("`fun` must be a non-empty character vector.", call. = FALSE)
  }

  allowed <- c("mean", "min", "max", "sd")
  unknown <- setdiff(fun, allowed)
  if (length(unknown)) {
    stop(
      "Unsupported summary function(s): ", paste(unknown, collapse = ", "),
      ". Supported values are: ", paste(allowed, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  fun <- unique(fun)

  if (!is.logical(na.rm) || length(na.rm) != 1 || is.na(na.rm)) {
    stop("`na.rm` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(force) || length(force) != 1 || is.na(force)) {
    stop("`force` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(quiet) || length(quiet) != 1 || is.na(quiet)) {
    stop("`quiet` must be TRUE or FALSE.", call. = FALSE)
  }

  paths <- .copernicus_summary_paths(file_paths, fun, output_dir, filename)

  time_index <- .copernicus_summary_time_index(
    file_paths = file_paths,
    start_datetime = start_datetime,
    end_datetime = end_datetime
  )

  if (!quiet) .copernicus_summary_time_message(time_index)

  selected_index <- time_index[time_index$selected, , drop = FALSE]
  if (!nrow(selected_index)) {
    stop("No netCDF time steps fall within the requested datetime range.",
         call. = FALSE)
  }

  if (all(file.exists(paths)) && !force) {
    if (!quiet) message("Cached summaries: ", paste(basename(paths), collapse = ", "))

    out <- normalizePath(paths, winslash = "/", mustWork = TRUE)
    names(out) <- names(paths)
    return(out)
  }

  nc0 <- ncdf4::nc_open(file_paths[1])
  on.exit(ncdf4::nc_close(nc0), add = TRUE)

  time_dim <- .copernicus_nc_time_dim(nc0)
  summary_vars <- .copernicus_nc_summary_vars(nc0, time_dim)

  specs <- lapply(summary_vars, function(v) {
    var <- nc0$var[[v]]
    dim_names <- vapply(var$dim, `[[`, character(1), "name")
    time_idx <- match(time_dim, dim_names)
    out_dims <- dim_names[-time_idx]
    out_lengths <- vapply(var$dim[-time_idx], `[[`, numeric(1), "len")

    list(
      name = v,
      dim_names = dim_names,
      time_idx = time_idx,
      out_dims = out_dims,
      out_lengths = out_lengths,
      units = if (is.null(var$units)) "" else var$units,
      longname = if (is.null(var$longname)) v else var$longname,
      prec = if (identical(var$prec, "double")) "double" else "float",
      missval = if (is.null(var$missval) || !is.finite(var$missval)) 1e20 else var$missval
    )
  })
  names(specs) <- summary_vars

  accumulators <- lapply(specs, function(s) {
    .copernicus_summary_accumulator(s$out_lengths, fun)
  })

  # Static variables are copied unchanged from the first file. Dimension
  # coordinate variables are created automatically by ncdim_def().
  first_dim_names <- names(nc0$dim)
  static_vars <- names(nc0$var)[vapply(nc0$var, function(v) {
    dims <- vapply(v$dim, `[[`, character(1), "name")
    !time_dim %in% dims
  }, logical(1))]
  static_vars <- setdiff(static_vars, first_dim_names)

  if (!quiet) {
    message(
      "Summarising ", length(summary_vars), " variable(s) across ",
      length(file_paths), " netCDF file(s): ",
      paste(summary_vars, collapse = ", ")
    )
  }

  # Stream through only the selected files/time steps, updating accumulators.
  selected_files <- unique(selected_index$file_index)

  for (f in selected_files) {
    nc <- if (f == 1) nc0 else ncdf4::nc_open(file_paths[f])
    close_after <- f > 1

    tryCatch(
      {
        this_time <- .copernicus_nc_time_dim(nc)
        if (!identical(this_time, time_dim)) {
          stop(
            "Input files use different time dimensions ('", time_dim, "' vs '",
            this_time, "').",
            call. = FALSE
          )
        }

        file_times <- selected_index$time_index[selected_index$file_index == f]

        for (v in summary_vars) {
          if (!v %in% names(nc$var)) {
            stop("Variable '", v, "' is missing from file: ", file_paths[f],
                 call. = FALSE)
          }

          var <- nc$var[[v]]
          dims <- vapply(var$dim, `[[`, character(1), "name")
          time_idx <- match(time_dim, dims)

          if (is.na(time_idx)) {
            stop("Variable '", v, "' has no time dimension in file: ",
                 file_paths[f], call. = FALSE)
          }

          out_dims <- dims[-time_idx]
          out_lengths <- vapply(var$dim[-time_idx], `[[`, numeric(1), "len")

          if (!identical(out_dims, specs[[v]]$out_dims) ||
              !identical(as.numeric(out_lengths),
                         as.numeric(specs[[v]]$out_lengths))) {
            stop(
              "Non-temporal dimensions for variable '", v,
              "' are not compatible across input files.",
              call. = FALSE
            )
          }

          if (f > 1) {
            for (dn in out_dims) {
              ref <- nc0$dim[[dn]]$vals
              cur <- nc$dim[[dn]]$vals
              if (length(ref) != length(cur) ||
                  !isTRUE(all.equal(ref, cur, tolerance = 1e-10,
                                   check.attributes = FALSE))) {
                stop(
                  "Coordinate values for dimension '", dn,
                  "' differ across input files.",
                  call. = FALSE
                )
              }
            }
          }

          start <- rep(1, length(dims))
          count <- vapply(var$dim, `[[`, numeric(1), "len")
          count[time_idx] <- 1

          for (tt in file_times) {
            start[time_idx] <- tt

            x <- ncdf4::ncvar_get(
              nc,
              v,
              start = start,
              count = count,
              collapse_degen = FALSE
            )

            slice <- array(as.numeric(x), dim = specs[[v]]$out_lengths)
            accumulators[[v]] <- .copernicus_summary_update(
              accumulators[[v]], slice, fun, na.rm
            )
          }
        }
      },
      finally = {
        if (close_after) ncdf4::nc_close(nc)
      }
    )

    if (!quiet) {
      message(
        "Processed ", sum(selected_index$file_index == f),
        " selected time step(s) from: ", basename(file_paths[f])
      )
    }
  }

  # Determine all non-time dimensions needed by summary and static variables.
  used_dims <- unique(unlist(lapply(specs, `[[`, "out_dims"), use.names = FALSE))
  if (length(static_vars)) {
    static_dims <- unlist(lapply(static_vars, function(v) {
      vapply(nc0$var[[v]]$dim, `[[`, character(1), "name")
    }), use.names = FALSE)
    used_dims <- unique(c(used_dims, static_dims))
  }
  used_dims <- setdiff(used_dims, time_dim)

  dim_defs <- .copernicus_nc_dim_defs(nc0, used_dims)

  # Build definitions for summarised variables.
  summary_defs <- lapply(specs, function(s) {
    ncdf4::ncvar_def(
      name = s$name,
      units = s$units,
      dim = unname(dim_defs[s$out_dims]),
      missval = s$missval,
      longname = s$longname,
      prec = s$prec
    )
  })
  names(summary_defs) <- summary_vars

  # Build definitions for static auxiliary variables.
  static_defs <- lapply(static_vars, function(v) {
    var <- nc0$var[[v]]
    dims <- vapply(var$dim, `[[`, character(1), "name")
    missval <- if (is.null(var$missval) || !is.finite(var$missval)) 1e20 else var$missval

    ncdf4::ncvar_def(
      name = v,
      units = if (is.null(var$units)) "" else var$units,
      dim = unname(dim_defs[dims]),
      missval = missval,
      longname = if (is.null(var$longname)) v else var$longname,
      prec = var$prec
    )
  })
  names(static_defs) <- static_vars

  # Write one native netCDF per requested summary statistic.
  for (stat in fun) {
    dest <- paths[[stat]]

    if (file.exists(dest) && !force) {
      if (!quiet) message("Cached: ", basename(dest))
      next
    }
    if (file.exists(dest) && force) unlink(dest, force = TRUE)

    nc_out <- ncdf4::nc_create(
      dest,
      vars = c(unname(summary_defs), unname(static_defs)),
      force_v4 = TRUE
    )
    

    # Copy global metadata before adding this package-specific provenance.
    .copernicus_nc_copy_atts(
      nc0, nc_out, 0, 0,
      exclude = c("history")
    )

    history_in <- tryCatch(ncdf4::ncatt_get(nc0, 0, "history")$value,
                           error = function(e) NULL)
    history_new <- paste0(
      format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      " copernicus_summarise(): temporal ", stat,
      " across ", nrow(selected_index), " selected time step(s) from ",
      length(unique(selected_index$file_index)), " file(s)."
    )
    history <- if (is.null(history_in) || !length(history_in) || is.na(history_in)) {
      history_new
    } else {
      paste(history_in, history_new, sep = "\n")
    }

    ncdf4::ncatt_put(nc_out, 0, "history", history)
    ncdf4::ncatt_put(nc_out, 0, "temporal_summary", stat)
    ncdf4::ncatt_put(
      nc_out, 0, "source_files",
      paste(basename(file_paths), collapse = ", ")
    )
    ncdf4::ncatt_put(
      nc_out, 0, "summary_start_datetime",
      format(min(selected_index$datetime), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    )
    ncdf4::ncatt_put(
      nc_out, 0, "summary_end_datetime",
      format(max(selected_index$datetime), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    )
    ncdf4::ncatt_put(nc_out, 0, "summary_time_steps", nrow(selected_index))

    for (v in summary_vars) {
      values <- .copernicus_summary_finalize(accumulators[[v]], stat, na.rm)
      ncdf4::ncvar_put(nc_out, v, values)

      .copernicus_nc_copy_atts(
        nc0, nc_out, v, v,
        exclude = c(
          "_FillValue", "missing_value", "scale_factor", "add_offset",
          "cell_methods"
        )
      )

      original_cell_methods <- tryCatch(
        ncdf4::ncatt_get(nc0, v, "cell_methods")$value,
        error = function(e) NULL
      )
      summary_method <- paste0("time: ", stat)
      if (!is.null(original_cell_methods) && length(original_cell_methods) &&
          !is.na(original_cell_methods) && nzchar(original_cell_methods)) {
        summary_method <- paste(original_cell_methods, summary_method)
      }

      ncdf4::ncatt_put(nc_out, v, "cell_methods", summary_method)
      ncdf4::ncatt_put(nc_out, v, "temporal_summary", stat)
    }

    for (v in static_vars) {
      ncdf4::ncvar_put(nc_out, v, ncdf4::ncvar_get(nc0, v, collapse_degen = FALSE))
      .copernicus_nc_copy_atts(
        nc0, nc_out, v, v,
        exclude = c("_FillValue", "missing_value", "scale_factor", "add_offset")
      )
    }

    ncdf4::nc_close(nc_out)

    if (!file.exists(dest) || file.info(dest)$size == 0) {
      stop("Failed to create summarised netCDF file: ", dest, call. = FALSE)
    }

    if (!quiet) message("Written: ", basename(dest))
  }

  out <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  names(out) <- names(paths)
  out
}
