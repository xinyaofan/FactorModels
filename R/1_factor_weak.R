library(FactorModels)
library(VineCopula)
library(CopulaModel)
#'nsize simulation size
#'n sample size
#'d data dimension
#'thlow lower bound for theta
#'thupp upper bound for theta

gl=gausslegendre(25)
xl=gl$nodes
wl=gl$weights

simu1fctweak<-function(nsize,n,d,seed=1234){
	set.seed(seed)
	result=matrix(NA,nsize,15)
	res_list=list()
	vlat_list=list()
	thv_mat=matrix(NA,nsize,2*d-1)
	for(k in 1:nsize){
		data=simWeakRes(n=n,d=d,thlow1=1.67,thupp1=5.0,
										thlow2=0.95,thupp2=2.95,fam1=4,fam2=5)
		udata=data$simdata[,1:d]
		vlat=data$simdata[,d+1]
    thv=data$thv
    thv_mat[k,]=thv
		#exact approach
		param0=c(rbind(c(rep(2,d),rep(2,d-1)),rep(0,2*d-1))) #initial value
		ifixed=rep(c(FALSE,TRUE),2*d-1)
		dstruct=list(data=udata,edg1=c(1:(d-1)),edg2=c(2:d),
								 xl=xl,wl=wl,fam=c(rep(4,d),rep(5,d-1)))
		LB=c(rbind(c(rep(1,d),rep(0.0,d-1)),rep(-0.999,2*d-1)))
		UB=c(rbind(c(rep(20,d),rep(20,d-1)),rep(0.001,2*d-1)))

		mlex1=pdhessminb(param=param0,objfn=f901f1tnllk,dstruct=dstruct,
									ifixed=ifixed,eps=1e-05,
		 							LB=LB,UB=UB,iprint = F)$parmin[!ifixed]

    fam=c(rep(4,d),rep(5,d-1))
		out2=proxy1fctweak(udata=udata,fam=fam,start=param0,
									LB=LB,UB=UB,xl=xl,wl=wl,
									iprint = T,ifixed=ifixed)

		proxyMean=out2$proxyMean
		proxyNew=out2$proxyNew
		mlpx2=out2$mlpx1
		mlpx3=out2$mlpx2
		mae=c(mean(abs(proxyMean-vlat)),
					mean(abs(proxyNew-vlat)),
					mean(abs(proxyMean-proxyNew)))
		rmse=c(sqrt(mean((proxyMean-vlat)^2)),
					 sqrt(mean((proxyNew-vlat)^2)),
					 sqrt(mean((proxyNew-proxyMean)^2)))
		res=list("v1"=proxyMean,"v2"=proxyNew,
						 "dif1"=proxyMean-vlat,"dif2"=proxyNew-vlat,
						 "theta1"=mlex1,"theta2"=mlpx2,"theta3"=mlpx3,
						 "rmse1"=sqrt(mean((mlex1-thv)^2)),
						 "rmse2"=sqrt(mean((mlpx2-thv)^2)),
						 "rmse3"=sqrt(mean((mlpx3-thv)^2)),
						 "mae1"=mean(abs(mlex1-thv)),
						 "mae2"=mean(abs(mlpx2-thv)),
						 "mae3"=mean(abs(mlpx3-thv)),
						 "mae_proxy"=mae,"rmse_proxy"=rmse,
						 "dif_para"=mean(abs(mlpx2-mlex1)),
						 "dif_para2"=mean(abs(mlpx3-mlex1)),
						 "dif_para3"=mean(abs(mlpx3-mlpx2)))
		res_list[[k]]=res
		vlat_list[[k]]=vlat
		result[k,]=c(res$rmse1,res$rmse2,res$rmse3,res$mae1,res$mae2,
								 res$mae3,res$mae_proxy,res$rmse_proxy,res$dif_para,
								 res$dif_para2,res$dif_para3)
		print(k)
	}
	colnames(result)=c("rmse1","rmse2","rmse3","mae1","mae2","mae3",
										 "mae_proxy1","mae_proxy2","mae_proxy12",
										 "rmse_proxy1","rmse_proxy2","rmse_proxy12",
										 "dif_para","dif_para2","dif_para3")
  return(list(result=result,res=res_list,vlat_list=vlat_list,
  						thv_mat=thv_mat))
	}



thv10=simu1fctweak(nsize=1000,n=500,d=10,seed=1234)
thv20=simu1fctweak(nsize=1000,n=500,d=20,seed=1234)
thv30=simu1fctweak(nsize=1000,n=500,d=30,seed=1234)
thv40=simu1fctweak(nsize=1000,n=500,d=40,seed=1234)
thv50=simu1fctweak(nsize=1000,n=500,d=50,seed=1234)
save(file="thv_list.RData",thv10,thv20,thv30,thv40,thv50)

res10=simu1fctweak(nsize=1000,n=500,d=10,seed=1234)

save(file="res10_1factweak.RData",res10)

res20=simu1fctweak(nsize=1000,n=500,d=20,seed=1234)
save(file="res20_1factweak.RData",res20)

res30=simu1fctweak(nsize=1000,n=500,d=30,seed=1234)
save(file="res30_1factweak.RData",res30)

res40=simu1fctweak(nsize=1000,n=500,d=40,seed=1234)
save(file="res40_1factweak.RData",res40)


res50=simu1fctweak(nsize=1000,n=500,d=50,seed=1234)
save(file="res50_1factweak.RData",res50)
