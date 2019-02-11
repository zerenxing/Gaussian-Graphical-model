rm(list=ls())

code.generate = function(par.vals, par.vals2,par.vals3, par.vals4,
                         pat.nams, pat.nams2,pat.nams3, pat.nams4, path="", path.opt=F,
                         templatename, out_name){ 
  
  for(i in 1:length(par.vals)){
     for(j in 1:length(par.vals2)){
       for(k in 1:length(par.vals3)){
         SourceCode=readLines(templatename) #read desired code
         SourceCode=gsub(pattern=pat.nams[i],replacement=par.vals[i],x=SourceCode) 
         SourceCode=gsub(pattern=pat.nams2[j],replacement=par.vals2[j],x=SourceCode) 
         SourceCode=gsub(pattern=pat.nams3[k],replacement=par.vals3[k],x=SourceCode)
         SourceCode=gsub(pattern=pat.nams4,replacement=par.vals4,x=SourceCode)
         
         #output changed code
         if(i==4){
           print(SourceCode)
         }
         path=ifelse(path==""|path.opt==F,"",paste(path,"/",sep=""))
         write(SourceCode,file=paste(path, out_name,
                                   k+(j-1)*length(par.vals3)+(i-1)*length(par.vals3)*length(par.vals2), 
                                   '.R',sep=""),sep='\n')  
       }
     }
  }
}  

######################################################################

#ntree
numvals = c(50,100,200)
numvals2 = c(200,400,800,1600,3200)
numvals3 = 0.5
numvals4 = "\"random\""
#numvals=4

nams = rep("PP",times=length(numvals))
nams2 = rep("NN",times=length(numvals2))
nams3 = rep("ZZ",times=length(numvals3))
nams4 = "TRUEL"

code.generate(numvals,numvals2,numvals3,numvals4,nams,nams2,nams3,nams4,templatename='Template-0mean.R', out_name='random_')

numvals4 = "\"AR3\""
code.generate(numvals,numvals2,numvals3,numvals4,nams,nams2,nams3,nams4,templatename='Template-0mean.R', out_name='AR3_')
