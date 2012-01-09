
# this function loads the compiled C-code for spacodi.calc() 
## hopefully, this alleviates issues with .so or .dll in cross-platform installation
.First.lib <- function(lib, pkg) {
  library.dynam("spacodi",pkg,lib)
}
