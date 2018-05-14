.onAttach <- function(libname, pkgname) {
  packageStartupMessage(rep("-", 10), "\nsinglearm: Design and Analysis of Single-Arm Clinical Trials\n",
                        rep("-", 10), "\n\nv.1.0: For an overview of the package's functionality enter: ?singlearm\n\n",
                        "For news on the latest updates enter: news(package = \"singlearm\")")
}

.onUnload <- function (libpath) {
  library.dynam.unload("singlearm", libpath)
}
