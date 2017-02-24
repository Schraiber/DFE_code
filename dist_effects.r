require(gtools)
require(DAAG)
require(hypergeo)
require(msm)

read.simsfs = function(path,is.zero=1) {
    files = list.files(path,pattern="selmathum",full.names=T)
    files_sorted = mixedsort(files) #note that files_sorted[is.zero] corresponds to s = 0!!
    test_scan = scan(files[1],quiet=T)
    print(files[1])
    print(test_scan)
    num_dat = length(test_scan)
    print(paste("Reading simulated SFS from ",length(files), " files, each for ", num_dat, " individuals",sep=""))
    sfs_mat = matrix(nrow=length(files_sorted),ncol=num_dat)
    for (i in 1:length(files_sorted)) {
        sfs_mat[i,] = scan(files_sorted[i],quiet=T)
    }
    s = c()
    for (i in 1:length(files_sorted)) {
	splitName = unlist(strsplit(files_sorted[i],"_"))
        cur_s = as.numeric(splitName[length(splitName)-2])
        s = c(s,cur_s)
    }
    good_s = c(s[is.zero],rev(10^(-s[-is.zero])))
    good_sfs = rbind(sfs_mat[1,],sfs_mat[length(files_sorted):2,])
    return(list(sfs=good_sfs,s=good_s))
}

riemann_integral = function(f,x) {
	sum(filter(f,rep(1/2,2))[-length(f)]*diff(x))
}

#NOTE: ONLY FOLD THE UNNORMALIED SFS
fold_sfs = function(sfs) {
	n = length(sfs)+1
	nHalf = floor(n/2) 
	eta = sapply(1:nHalf,function(i){(sfs[i]+sfs[n-i])/(1+as.numeric(i==(n-i)))})
	return(eta)
}

#evalulate the integral for a given frequency
#assumes that freq_counts and s are sorted in the same way
#and that s[1] = 0, and 0 < s[1] < s[2] < s[3] < ...
#*_extra is in case the simulations don't have enough. Provide some more s!
#ONLY FOR FREQ 1 ALLELES
int_freq_gamma = function(freq_counts, s, alpha, beta,freq_counts_extra = c(), s_extra = c()) { 
    #compute up to s[1] using a linear interpolation and pgamma
    m = (freq_counts[2]-freq_counts[1])/(s[2]-s[1])
    b = freq_counts[1]
    #i1 = freq_counts[1]*pgamma(s[2],alpha,beta)
    i1 = m*alpha/beta*pgamma(s[2],alpha+1,beta)+b*pgamma(s[2],alpha,beta)
    #print(i1)
    #compute the rest using trapezoid rule
    i2 = riemann_integral(freq_counts[2:length(freq_counts)]*dgamma(s[2:length(s)],alpha,beta),s[2:length(s)])
    #plot(s[2:length(s)],freq_counts[2:length(freq_counts)]*dgamma(s[2:length(s)],alpha,beta))
    if (length(freq_counts_extra)&&length(s_extra)) {
	i3 = riemann_integral(freq_counts_extra,s_extra)
	i2 = i2 + i3	
    }
    #pause()
    #print(c(i1,i2))
    #add together
    return(i1+i2)
}

#truncated normal distribution of fitness effects
int_freq_normal = function(freq_counts,s,mu,sigma) {
	riemann_integral(freq_counts*dtnorm(s,mu,sigma,0,Inf),s)
}


#p1 and p2 are the parameters of the distributions
compute_sfs = function(sfs_mat,s,p1,p2,norm=TRUE,gamma=TRUE,fold=FALSE) {
    #print(c(alpha, beta))
    sfs = vector(length=ncol(sfs_mat))
    for (i in 1:ncol(sfs_mat)) {
	#print(i)
	if (gamma) {
        	cur_freq = int_freq_gamma(sfs_mat[,i],s,p1,p2)
	} else {
		cur_freq = int_freq_normal(sfs_mat[,i],s,p1,p2)
	}
	sfs[i] = cur_freq
	#print(c(i,sfs[i]))
    }
    if (fold) {
        sfs = fold_sfs(sfs)
    }
    if (norm) {
        sfs = sfs/sum(sfs)
    }
    #print(head(s))
    #print(tail(s))
    #print(head(sfs_mat[,1]))
    #print(tail(sfs_mat[,1]))
    #print(c(p1,p2))
    #print(sfs)
    return(sfs)
}

fit_single_coefficient = function(data,sims,norm = TRUE,fold=FALSE,num_sam = max(data[,1])+1,min_c=0,max_c=40,info=FALSE) {
    LL_mat = matrix(nrow=(max_c-min_c+1),ncol=length(sims$s))
    for (i in min_c:max_c) {
        #print(i)
        cur_sfs = get_c_sfs(i,data)
	if (fold) {
 	    cur_sfs = fold_sfs(cur_sfs)
	}
        cur_LL = c()
        for (j in 1:nrow(sims$sfs)) {
            cur_sim = sims$sfs[j,]
	    if (fold) {
		cur_sim = fold_sfs(cur_sim)
            }
            if (norm) {
                cur_sim = cur_sim/sum(cur_sim)
            }
            #print(c(length(cur_sfs),length(cur_sim)))
            #pause()
            cur_LL = c(cur_LL, -sum(cur_sfs*log(cur_sim)))
        }
        LL_mat[i+as.numeric(min_c==0),] = cur_LL
    }
    best_per_c = apply(LL_mat,1,which.min)
    best_s = sims$s[best_per_c]
    best_LL = apply(LL_mat,1,function(x){x[which.min(x)]})
    return(list(LL_mat=LL_mat,s=best_s,LL=best_LL,c=min_c:max_c))
}

#Fit a single coefficient to the CORRECTED spectrum
fit_single_coefficient_corrected = function(data,sims,post.p,norm = TRUE,fold=FALSE,num_sam = max(data[,1])+1,min_c=1,max_c=max(data[,2])) {
    LL_mat = matrix(nrow=(max_c-min_c+1),ncol=length(sims$s))
    for (i in min_c:max_c) {
        #print(i)
        cur_sfs = c()
        for (k in 1:(num_sam-1)) {
            cur_freq = data[data[,2]==i&data[,1]==k,3]
            if (length(cur_freq)) {
                cur_sfs = c(cur_sfs,cur_freq)
            } else {
                cur_sfs = c(cur_sfs,0)
            }
        }	
	cur_sfs = rev(post.p$mis*cur_sfs)+post.p$notmis*cur_sfs
	if (fold) {
 	    cur_sfs = fold_sfs(cur_sfs)
	}
        cur_LL = c()
        for (j in 1:nrow(sims$sfs)) {
            cur_sim = sims$sfs[j,]
	    if (fold) {
		cur_sim = fold_sfs(cur_sim)
            }
            if (norm) {
                cur_sim = cur_sim/sum(cur_sim)
            }
            #print(c(length(cur_sfs),length(cur_sim)))
            #pause()
            cur_LL = c(cur_LL, -sum(cur_sfs*log(cur_sim)))
        }
        LL_mat[i+as.numeric(min_c==0),] = cur_LL
    }
    best_per_c = apply(LL_mat,1,which.min)
    best_s = sims$s[best_per_c]
    best_LL = apply(LL_mat,1,function(x){x[which.min(x)]})
    return(list(LL_mat=LL_mat,s=best_s,LL=best_LL,c=min_c:max_c))
}

#fits the mean and variance parameters
fit_gamma_dist = function(data,sims,min_c = 0, max_c = 40,fold=FALSE) {
    optims = list()
    for (i in min_c:max_c) {
        print(i)
        cur_sfs = get_c_sfs(i,data)
	if (fold) {
		cur_sfs = fold_sfs(cur_sfs)
	}
        optims[[as.character(i)]] = optim(par=c(1e-5,1e-5),function(pars) {-sum(cur_sfs*log(compute_sfs(sims$sfs,sims$s,pars[1]^2/pars[2]^2,pars[1]/pars[2]^2,fold=fold)))}, method="L-BFGS-B",lower=c(1e-7,1e-6),upper=c(.01,.01),control=list(factr=1e-2,ndeps=c(1e-6,1e-6)))
	print(optims[[as.character(i)]]$par)
    }
    return(list(opt=optims,c=min_c:max_c))
}

#fits the mean and variance parameters, takes account of sfs correction!!!
#NB: assumes a FIXED correction, based on demography
fit_gamma_dist_corrected = function(data,sims,post.p,min_c = 1, max_c = max(data[,2]),fold=FALSE) {
    optims = list()
    for (i in min_c:max_c) {
        print(i)
        cur_sfs = c()
        for (k in 1:max(data[,1])) {
            cur_freq = data[data[,2]==i&data[,1]==k,3]
            if (length(cur_freq)) {
                cur_sfs = c(cur_sfs,cur_freq)
            } else {
                cur_sfs = c(cur_sfs,0)
            }
        }
	#CORRECT THE SFS
	cur_sfs = rev(post.p$mis*cur_sfs)+post.p$notmis*cur_sfs
	if (fold) {
		cur_sfs = fold_sfs(cur_sfs)
	}
        optims[[i+as.numeric(min_c==0)]] = optim(par=c(1e-5,1e-5),function(pars) {-sum(cur_sfs*log(compute_sfs(sims$sfs,sims$s,pars[1]^2/pars[2]^2,pars[1]/pars[2]^2,fold=fold)))}, method="L-BFGS-B",lower=c(1e-7,1e-6),upper=c(.005,.005),control=list(factr=1e-2,ndeps=c(1e-6,1e-6)))
	print(optims[[i+as.numeric(min_c==0)]]$par)
    }
    return(list(opt=optims,c=min_c:max_c))
}

mixture_of_gammas = function(fit_dist,countData = c(), tsvData = c(),min_c=0,max_c=40) {
	if (length(countData) && !length(tsvData)) {
		frac_per_c = sapply(min_c:max_c, function(c) {sum(countData[countData[,2]==c,3])})
		frac_per_c = frac_per_c/sum(frac_per_c)
	} else if (length(tsvData) && !length(countData)) {
		frac_per_c = hist(round(tsvData[,7]),breaks=0:100,plot=FALSE,right=FALSE)$density[min_c:max_c+1]
	} else {
		print("Give me either countData or tsvData!")
		return(0)
	}
	best_alpha = sapply(fit_dist,function(x){x$par[1]^2/x$par[2]^2})[min_c:max_c] 
	best_beta = sapply(fit_dist,function(x){x$par[1]/x$par[2]^2})[min_c:max_c]
	return( list(dist=function(sel_coef) {sapply(sel_coef,function(s){sum(frac_per_c*dgamma(s,best_alpha,best_beta))})} ,frac=frac_per_c,best_alpha=best_alpha,best_beta=best_beta))
}

mixture_of_gammas_em = function(fit_dist,countData = c(), tsvData = c(),min_c=0,max_c=40) {
	if (length(countData) && !length(tsvData)) {
		frac_per_c = sapply(min_c:max_c, function(c) {sum(countData[countData[,2]==c,3])})
		frac_per_c = frac_per_c/sum(frac_per_c)
	} else if (length(tsvData) && !length(countData)) {
		frac_per_c = hist(round(tsvData[,7]),breaks=0:100,plot=FALSE,right=FALSE)$density[(min_c+as.numeric(min_c==0)):(max_c+as.numeric(min_c==0))]
		frac_per_c = frac_per_c/sum(frac_per_c)
	} else {
		print("Give me either countData or tsvData!")
		return(0)
	}
	best_alpha = sapply(fit_dist,function(x){x$mean^2/x$sd^2})[(min_c+as.numeric(min_c==0)):(max_c+as.numeric(min_c==0))] 
	best_beta = sapply(fit_dist,function(x){x$mean/x$sd^2})[(min_c+as.numeric(min_c==0)):(max_c+as.numeric(min_c==0))]
	return( list(dist=function(sel_coef) {sapply(sel_coef,function(s){sum(frac_per_c*dgamma(s,best_alpha,best_beta))})} ,frac=frac_per_c,best_alpha=best_alpha,best_beta=best_beta))
}

fit_tnorm_dist = function(data,sims,num_sam = max(data[,1])-1) {
    optims = list()
    for (i in min(data[,2]):max(data[,2])) {
        print(i)
        cur_sfs = c()
        for (k in 1:(num_sam-1)) {
            cur_freq = data[data[,2]==i&data[,1]==k,3]
            if (length(cur_freq)) {
                cur_sfs = c(cur_sfs,cur_freq)
            } else {
                cur_sfs = c(cur_sfs,0)
            }
        }
	print(cur_sfs)
        optims[[i]] = optim(par=c(.001,.001),function(pars) {-sum(cur_sfs*log(compute_sfs(sims$sfs,sims$s,pars[1],pars[2],gamma=FALSE)))}, method="L-BFGS-B",lower=1e-7,upper=c(.01,.005),control=list(ndeps=c(1e-6,1e-6)))
	print(optims[[i]]$par)
    }
    return(optims)
}

get_c_sfs = function(c,data) {
    cur_sfs = c()
    for (k in 1:max(data[,1])) {
        cur_freq = data[data[,2]==c&data[,1]==k,3]
        if (length(cur_freq)) {
            cur_sfs = c(cur_sfs,cur_freq)
        } else {
            cur_sfs = c(cur_sfs,0)
        }
    }
    return(cur_sfs)
}

pop_sfs = function(x,alpha) {
	if (alpha) {
		if (alpha < -350) {
			1/(x*(1-x))*exp(2*alpha*x)
		} else {
			1/(x*(1-x))*(1-exp(-2*alpha*(1-x)))/(1-exp(-2*alpha))
		}
	} else {
		1/x
	}
}

sfs_integrand = function(x,n,k,alpha) {
	choose(n,k)*x^k*(1-x)^(n-k)*pop_sfs(x,alpha)
}

exp_spec_numerical = function(alpha,n,norm=TRUE) { 
	sample = sapply(1:(n-1),function(k){integrate(sfs_integrand,0,1,n=n,k=k,alpha=alpha)$value})
	if (norm) {
		return(sample/sum(sample))
	} else {
		return(sample)
	}
}

exp_spec = function(s,n,norm=TRUE,maxiter=2000,series=TRUE) {
	#generates the expected spectrum with constant demography
	if (s == 0) {
		sample = 1/1:(n-1)
	} else {
		hypergeos = sapply(1:(n-1),function(k){genhypergeo(k,n,2*s,maxiter=maxiter,series=series)})
		k = 1:(n-1)
		sample = n/(2*k*(n-k))*(1/tanh(s)-1)*(exp(2*s)-hypergeos)
	}
	if (any(sample<0)) {
		print(s)
		print(sample)
		return(NA)
	}
	if (norm) {
		return(sample/sum(sample))
	} else {
		return(sample)
	}		
}

const_spec_like = function(s,dat) {
	n = length(dat)+1
	mySpec = exp_spec(s,n)
	LL = sum(dat*log(mySpec))
	return(-LL)	
}

rfreq = function(s,n) {
	spec = exp_spec(s,n)
	which(rmultinom(1,1,spec)==1)	
}

simulate_data = function(n,all_c_scores,c_to_s,...,num_sites=length(all_c_scores)) {
	#n is the sample size
	#all_c_scores is all the c_scores from the data
	#c_to_s is a fuenction that takes c scores and maps them to selection coefficients. SHOULD BE VECTORIZED
	#... is parameters of c_to_s
	#num_sites is how many sites to simulate
	new_c = sample(all_c_scores,num_sites,replace=TRUE)
	true_s = c_to_s(new_c,...)
	freqs = sapply(true_s,function(s){rfreq(s,n)})	
	return(list(freq=freqs,c=new_c))
}

fit_dist_sim = function(sim_data,n,lower=0,upper=40.5,by=1) {
	#create the bins
	bin_ends = c(lower,seq(lower+by,upper,by))
	bin_sfs = matrix(nrow=(length(bin_ends)-1),ncol=n-1) #MAYBE length(bin_ends)-1???
	for (i in 1:(length(bin_ends)-1)) {
		bin_sfs[i,] = hist(sim_data$freq[sim_data$c>=bin_ends[i]&sim_data$c<bin_ends[i+1]],breaks=0:(n-1),plot=F)$counts		
	}
	#return(bin_sfs)
	#fit to each bin
	best_s = c()
	for (i in 1:nrow(bin_sfs)) {
		print(i)
		if (sum(bin_sfs[i,]) <= 5) {
			print(i)
			print(bin_sfs[i,])
			best_s = c(best_s,NA)
		} else {
			best_s = c(best_s, optim(0,const_spec_like,gr=NULL,dat=bin_sfs[i,],method="L-BFGS-B",lower=-3,upper=3)$par)
		}	
		print(best_s)
	}
	return(list(sfs=bin_sfs,s=best_s))	
}

read.demographic.sfs = function(path) {
    files = list.files(path,full.names=T,pattern="exp")
    test_scan = scan(files[1],quiet=T)
    num_dat = length(test_scan)
    sfs_mat = matrix(nrow=length(files),ncol=num_dat)
    for (i in 1:length(files)) {
        sfs_mat[i,] = scan(files[i],quiet=T)
    }
    time = c()
    rate = c()
    for (i in 1:length(files)) {
	subName = gsub("exp","",files[i])
	subName = gsub(".txt","_.txt",subName)
	subName = gsub("_t","_",subName)
	splitName = unlist(strsplit(subName,"_"))
        cur_rate = splitName[length(splitName)-1]
	if (cur_rate == "constant") {
		cur_rate = 0
	} else {
		cur_rate = as.numeric(cur_rate)
	}
	cur_time = as.numeric(splitName[length(splitName)-2])
	if (is.na(cur_time) || is.na(cur_rate)) {
		print(files[i])
		pause()
	}
        time = c(time,cur_time)
	rate = c(rate,cur_rate)
    }  
    return(list(sfs=sfs_mat,time=time,rate=rate))

}

correct_sfs = function(data,post_probs) {
	post_probs$notmis*data + rev(post_probs$mis*data)
}

#returns the prob for the sfs supplied, ASSUMED NORMALIZED
#a vector, where each entry corresponds to an allele frequency
post_probs = function(sfs,cur_probs) {
	#print(sfs)
	sfs = sfs/sum(sfs)
	n = length(sfs)+1
	i = 1:(n-1)
	mis = sfs[n-i]*cur_probs$mis
	notmis = sfs[i]*cur_probs$notmis
	#print(sfs)
	#print(mis)
	#print(notmis)
	#print(mis/(mis+notmis))
	#print(notmis/(mis+notmis))
	#pause()
	return(list(mis=mis/(mis+notmis),notmis=notmis/(mis+notmis)))
}

#data is an SFS, UNNORMALIZED!!!
new_probs = function(data,post_probs) {
	N = sum(data)
	#print(data)
	#print(post_probs$mis)
	#print(post_probs$notmis)
	#print(data*post_probs$mis)
	#print(data*post_probs$notmis)
	#pause()
	new_mis = sum(data*post_probs$mis)/N
	new_not_mis = sum(data*post_probs$notmis)/N
	return(list(mis=new_mis,notmis=new_not_mis))
}

#data is an SFS, UNNORMALIZED
#theory is a THEORETICAL sfs
weighted_sfs_LL = function(data,theory,post_probs) {
	theory = theory/sum(theory)
	LLmis = sum(data*post_probs$mis*log(rev(theory)))
	LLnotmis = sum(data*post_probs$notmis*log(theory))
	return(LLmis+LLnotmis)
}

weighted_sfs_LL_int = function(data,sims,post_probs,p1,p2) {
	theory = compute_sfs(sims$sfs,sims$s,p1,p2)
	LLmis = -sum(data*post_probs$mis*log(rev(theory)))
	LLnotmis = -sum(data*post_probs$notmis*log(theory))
	#print(c(p1,p2))
	#print(data)
	#print(theory)
	#print(c(post_probs$mis,post_probs$notmis))
	#print(c(LLmis,LLnotmis))
	#pause()
	return(LLmis+LLnotmis)
}

new_sfs_no_int = function(data,sims,post_probs) {
	LL = c()
	for (i in 1:nrow(sims$sfs)) {
		LL = c(LL,weighted_sfs_LL(data,sims$sfs[i,],post_probs))
	}
	return(list(max_ind=which.max(LL),LL=LL))
}

new_sfs_int = function(data,sims,post_probs,p1,p2,lower=c(1e-7,1e-5),upper=c(.02,.02)) {
	opt = optim(c(p1,p2),function(pars){weighted_sfs_LL_int(data,sims,post_probs,pars[1]^2/pars[2]^2,pars[1]/pars[2]^2)},gr=NULL,method="L-BFGS-B",lower=lower,upper=upper,control=list(factr=1e-2,ndeps=c(1e-6,1e-6)))
	cur_mean = opt$par[1]
	cur_sd = opt$par[2]
	cur_alpha = cur_mean^2/cur_sd^2
	cur_beta = cur_mean/cur_sd^2
	cur_sfs = compute_sfs(sims$sfs,sims$s,p1^2/p2^2,p1/p2^2)
	return(list(sfs=cur_sfs,mean=cur_mean,sd=cur_sd))	
}

EM_sfs_demo = function(data,sims,num_iter = 50, cur_sfs = 1,cur_probs=list(mis=.5,notmis=.5),info=FALSE) {
	for (i in 1:num_iter) {
		if (info) {
			print(cur_probs)
			print(c(sims$time[cur_sfs],sims$rate[cur_sfs]))
		}
		#print(sims$sfs[cur_sfs,])
		#pause()
		sfs = sims$sfs[cur_sfs,]
		post.p = post_probs(sfs,cur_probs)
		cur_probs = new_probs(data,post.p)
		cur_sfs = new_sfs_no_int(data,sims,post.p)$max_ind 	
	}
	sfs = sims$sfs[cur_sfs,]
	post.p = post_probs(sfs,cur_probs)
	cor_sfs = correct_sfs(data,post.p)
	return(list(theory=sfs,input=data,corrected=cor_sfs,probs=cur_probs,best_sfs_ind = cur_sfs,rate=sims$rate[cur_sfs],time=sims$time[cur_sfs]))
}

EM_sfs_no_int = function(data,sims,num_iter = 50, cur_sfs = 1,cur_probs=list(mis=.5,notmis=.5),info=FALSE) {
	for (i in 1:num_iter) {
		if (info) {
			print(cur_probs)
			print(sims$s[cur_sfs])
		}
		#print(sims$sfs[cur_sfs,])
		#pause()
		sfs = sims$sfs[cur_sfs,]
		post.p = post_probs(sfs,cur_probs)
		cur_probs = new_probs(data,post.p)
		cur_sfs = new_sfs_no_int(data,sims,post.p)$max_ind 	
	}
	sfs = sims$sfs[cur_sfs,]
	post.p = post_probs(sfs,cur_probs)
	cor_sfs = correct_sfs(data,post.p)
	return(list(theory=sfs,input=data,corrected=cor_sfs,probs=cur_probs,best_sfs_ind = cur_sfs,s=sims$s[cur_sfs]))
}

 
EM_sfs_int = function(data,sims,num_iter = 50, cur_mean = .00001, cur_sd = .00001,cur_probs=list(mis=.5,notmis=.5),info=FALSE,lower=c(1e-7,5e-6),upper=c(.02,.02)) {
	sfs = compute_sfs(sims$sfs,sims$s,cur_mean^2/cur_sd^2,cur_mean/cur_sd^2)
	for (i in 1:num_iter) {
		if (info) {
			print(paste("Iteration = ", i, sep=""))
			print(cur_probs)
			print(c(cur_mean,cur_sd))
			#pause()
		}
		post.p = post_probs(sfs,cur_probs)
		cur_probs = new_probs(data,post.p)
		sfs_opt = new_sfs_int(data,sims,post.p,cur_mean,cur_sd,lower=lower,upper=upper)	
		sfs = sfs_opt$sfs
		cur_mean = sfs_opt$mean
		cur_sd = sfs_opt$sd
	}
	post.p = post_probs(sfs,cur_probs)
	cor_sfs = correct_sfs(data,post.p)
	return(list(theory=sfs,input=data,corrected=cor_sfs,probs=cur_probs,mean=cur_mean,sd=cur_sd))
}

#NB: this returns a list and the keys are STRINGS!!!!
fit_em_int = function(data,sims,min_c=0,max_c=40,num_iter=50,cur_mean=.00001,cur_sd=.00001,cur_probs=list(mis=.5,notmis=.5),info=FALSE,lower=c(1e-7,5e-6),upper=c(.02,.02)) {
	fit_dist = list()
	for (c in min_c:max_c) {
		if (info) { 
			print(paste("C = ",c,sep=""))
		}
		cur_sfs = get_c_sfs(c=c,data=data)
		fit_dist[[as.character(c)]] = EM_sfs_int(cur_sfs,sims,num_iter=num_iter,info=info,lower=lower,upper=upper)
		if (info) {
			print("")
		}
	}	
	return(fit_dist)
}

fit_em_no_int = function(data,sims,min_c = 0, max_c = 40, num_iter = 50, cur_sfs = 1,cur_probs=list(mis=.5,notmis=.5),info=FALSE) {
	fit_coef = list()
	for (c in min_c:max_c) { 
		if (info) {
			print(paste("C = ",c,sep=""))
		}
		this_sfs = get_c_sfs(c=c,data=data)
		fit_coef[[as.character(c)]] = EM_sfs_no_int(data=this_sfs,sims=sims,num_iter = num_iter, cur_sfs = cur_sfs,cur_probs=cur_probs,info=info)
		if (info) {
			print("")
		} 
	}
	return(fit_coef)
}

fit_em_no_int_unique = function(data,sims,min_c = 0, max_c = 40, num_iter = 50, cur_sfs = 1,cur_probs=list(mis=.5,notmis=.5),info=FALSE) {
	fit_coef = list()
	for (c in unique(data[data[,2]<=max_c&data[,2]>=min_c,2])) { 
		if (info) {
			print(paste("C = ",c,sep=""))
		}
		this_sfs = get_c_sfs(c=c,data=data)
		fit_coef[[as.character(c)]] = EM_sfs_no_int(data=this_sfs,sims=sims,num_iter = num_iter, cur_sfs = cur_sfs,cur_probs=cur_probs,info=info)
		if (info) {
			print("")
		} 
	}
	return(fit_coef)
}

#computes the pseudo-likelihood of the corrected observations
LRT_em = function(fit_dist,fit_coef) {
	LRT = c()
	for (i in 1:length(fit_dist)) {
		dist_theory = fit_dist[[i]]$theory
		dist_theory = dist_theory/sum(dist_theory)
		dist_corrected = fit_dist[[i]]$corrected
		dist_likelihood = sum(dist_corrected*log(dist_theory))
		coef_theory = fit_coef[[i]]$theory
		coef_theory = coef_theory/sum(coef_theory)
		coef_corrected = fit_coef[[i]]$corrected
		coef_likelihood = sum(coef_corrected*log(coef_theory))
		LRT = c(LRT,2*(dist_likelihood-coef_likelihood))
	}
	return(LRT)
}

#this computes the true likelihood of the data under the inferred parameters
likelihood_em = function(fit_dist) {
	mis = fit_dist$probs$mis
	notmis = fit_dist$probs$notmis
	data = fit_dist$input
	theory = fit_dist$theory
	theory = theory/sum(theory)
	sum(data*log(notmis*theory + mis*rev(theory)))
}

plot_sfs_barplot = function(sfs_to_plot,log="",xlab="Number of derived alleles",ylab="Fraction of sites",legend=c()) {
	#normalize
	for (i in 1:nrow(sfs_to_plot)) {
		sfs_to_plot[i,] = sfs_to_plot[i,]/sum(sfs_to_plot[i,])
	}
	barplot(sfs_to_plot,beside=T,log=log,names.arg=1:ncol(sfs_to_plot),xlab=xlab,ylab=ylab)
	if (length(legend)) {
		legend("topright",fill=grey.colors(nrow(sfs_to_plot)),legend=legend)
	}
}
