# 
# setGeneric("GWLD", function(geno012, method, cores, ...){standardGeneric("GWLD")})
# 
# 
# 
# 
# setMethod("GWLD", "data.frame",
#           function(geno012, ...){
#             return(as.matrix(geno012))
#           })
# 
# setMethod("GWLD", "matrix",
#           function(geno012, ...){
#             return(as.matrix(geno012))
#           })








