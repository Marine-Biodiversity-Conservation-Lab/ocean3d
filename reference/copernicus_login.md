# Log in to Copernicus Marine

Configure or verify authentication for Copernicus Marine using the
official standalone Copernicus Marine Toolbox.

## Usage

``` r
copernicus_login(path = NULL)
```

## Arguments

- path:

  Optional character string giving the path to a Copernicus Marine
  Toolbox standalone executable. If `NULL`, the package-managed
  executable and then the system `PATH` are searched automatically.

## Value

Invisibly, `TRUE` when valid Copernicus Marine credentials were already
available and successfully verified. Invisibly, `FALSE` when an
interactive login has been launched and must be completed by the user.

## Details

`copernicus_login()` delegates credential discovery, validation, and
storage entirely to the official Copernicus Marine Toolbox. This package
does not search for, request, receive, read, or store the user's
Copernicus Marine username or password.

If valid credentials are already available, the Toolbox reuses them
automatically. If authentication still needs to be configured, an
interactive Toolbox login is started.

Authentication only needs to be configured once. The credentials created
by the official Copernicus Marine Toolbox are shared by both the
standalone and Python backends, so no separate Python login is required
before using `backend = "python"`.

When running inside RStudio, a first-time login is opened in the RStudio
Terminal because username and password prompts from external programs
may not receive interactive input correctly from the RStudio Console.

## Examples

``` r
if (FALSE) { # \dontrun{
copernicus_login()
} # }
```
