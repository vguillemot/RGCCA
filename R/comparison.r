#'Compares 2 RGCCA
#'
#' @param rgcca1 A result of a RGCCA function 
#' @param rgcca2 Another result of RGCCA function
#' @param nAxe=1 Number of axes taken into account during the comparison
#' @param selec Number of biomarkers to be selected
#' @param selectPatient A vector allowing to select only some patients for the RV calculation
#' @return \item{A} A list containing: a: the correlation between axes, rv: the rv coefficient, bm the biomarkers
#' @return \item{crit} Convergence criterion : abs(1-obj_k/obj_{k-1})
#' @return \item{obj} Vector containing the mean square error between the predict values and the original non missing values at each iteration
#' @title comparison of two RGCCA results
#' @examples 
#'  data();...
#'  

comparison=function(rgcca1,rgcca2,nAxe=1,selec=10,selectPatient=NULL,indNA=NULL)
{
  diffNorm2=function(vec1,vec2)
  {
    if(!is.null(vec1)&!is.null(vec2))
    {
      c1= cor(vec1,vec2)
      c2= cor(vec1,-vec2)
      return(max(c1,c2))	
    }
    else{ return(NA)}
  }
  selectAllPatient=intersect(rownames(rgcca1[["Y"]][[1]]),rownames(rgcca2[["Y"]][[1]]))
  J=length(rgcca1$A)
  com=rv=pctBm=rvComplete=rmse=rep(NA,J)
  for(i in 1:J)
  {
    refBm=biomarker(resRGCCA=rgcca1,block=i,axes=1,selec=selec)
    com[i]=diffNorm2(rgcca1[["astar"]][[i]][,nAxe],rgcca2[["astar"]][[i]][,nAxe])
    if(dim(rgcca1[["Y"]][[i]])[2]>1)
    {
      rvComplete[i]=coeffRV(rgcca1[["Y"]][[i]][selectPatient,1:2],rgcca2[["Y"]][[i]][selectPatient,1:2])$rv
      rv[i]=coeffRV(rgcca1[["Y"]][[i]][selectAllPatient,1:2],rgcca2[["Y"]][[i]][selectAllPatient,1:2])$rv
    }
   if(dim(rgcca1[["Y"]][[i]])[2]==1)
   {
     rvComplete[i]=diffNorm2(as.vector(rgcca1[["Y"]][[i]][selectPatient,]),as.vector(rgcca2[["Y"]][[i]][selectPatient,]))
      rv[i]=diffNorm2(as.vector(rgcca1[["Y"]][[i]][selectAllPatient,]),as.vector(rgcca2[["Y"]][[i]][selectAllPatient,]))
    }
    
    testBm=biomarker(resRGCCA=rgcca2,block=i,axes=1,selec=selec)
    pctBm[i]=sum(names(testBm)%in%names(refBm))/length(refBm)	
    
    if(!is.null(indNA[[i]]))
    {
      if(dim(rgcca1$Y[[i]])[1]==dim(rgcca2$Y[[i]])[1])
      {
        stdev=apply(rgcca1$A[[i]],2,sd)
        denom=matrix(rep(stdev,dim(rgcca1$Y[[i]])[1]),nrow=dim(rgcca1$Y[[i]])[1],ncol=dim(rgcca1$A[[i]])[2],byrow = TRUE)
        difRel=(rgcca1$A[[i]]-rgcca2$A[[i]]) /denom  
        rmse[i]=sqrt(sum(difRel[indNA[[i]]]^2)/prod(dim(indNA[[i]])))
      }
    }
  }
  return(list(a=com,rv=rv,bm=pctBm,rvComplete=rvComplete,rmse=rmse))
}
