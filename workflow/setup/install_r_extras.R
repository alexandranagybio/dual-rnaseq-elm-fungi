# Additional R packages not installed through environment.yml

required <- c(
  ltc = "0.4.0"
)

for (pkg in names(required)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }

  installed <- as.character(packageVersion(pkg))

  if (installed != required[[pkg]]) {
    warning(
      sprintf(
        "%s version %s is installed; manuscript environment used %s",
        pkg,
        installed,
        required[[pkg]]
      )
    )
  }
}
