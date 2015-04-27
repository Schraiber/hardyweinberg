require(EMT) #this has the exact multinomial test

multinomial = function(n_vec) {
	#n_vec is vector of desired guys
	n = sum(n_vec)
	to_return = lgamma(n+1)
	for (n_k in n_vec) {
		to_return = to_return - lgamma(n_k+1)
	}	
	return(exp(to_return))
}

hw_exact = function(n12,n,n1) {
	#n12 is number of hets, n is number of samples and n1 is number of copies of allele 1
	n11 = (n1 - n12)/2
	n22 = n-n11-n12
	if (n11 %% 1 != 0 || n22 %% 1 != 0) {
		return(0)
	}
	if (n11 < 0 || n22 < 0) {
		return(0)
	}
	return(multinomial(c(n11,n12,n22))/choose(2*n,n1)*2^(n12))
}

hw_exact_f = function(f,n12,n,n1) {
	#computes the exact probability with inbreeding coefficient f (derived by Monty)
	if (n1 %% 2 == 0) {
		np = n
	} else {
		np = n-1
	}
	
	#my n1 is Monty's i
	if (n1 %% 2 == 0) {
		n1p = n1
	} else {
		n1p = n1-1
	}
	
	#j is number of individuals who are ibd at that site
	
	#k is the number of individuals who are ibd for A (j-k is number who are ibd for a)
	
	to_return = 0
	for (j in 0:np) {
		for (k in 0:j) {
			#print(c(n-j, n1-2*k))
			to_return = to_return + hw_exact(n12,n-j,n1-2*k)*dbinom(j,np,f)*dhyper(k,n1p/2,np-n1p/2,j)
                        #print(c(np,n1p))
                        #print(c(k,n1p/2,np-n1p/2,j))
                        #print(dhyper(k,n1p/2,np-n1p/2,j))
		}
	}
	return(to_return)
}

#a function to generate the exact distribution assuming a beta distribution of inbreeding coefficients
hw_exact_f_beta = function(a,b,n12,n,n1) {
  if (n1 %% 2 == 0) {
    np = n
  } else {
    np = n - 1
  }
  if (n1 %% 2 == 0) {
    n1p = n1
  } else {
    n1p = n1-1
  }

  to_return = 0
  for (j in 0:np) {
    for (k in 0:j) {
      to_return = to_return + hw_exact(n12,n-j,n1-2*k)*dhyper(k,n1p/2,np-n1p/2,j)*choose(np,j)*exp(lbeta(j+a,np-j+b)-lbeta(a,b))
    }
  }
  return(to_return)
}

LL_HW = function(f,het_obs,n1) {
	#note that this is reindexing to go to 0!!
	ll = 0.0
	if (n1 %% 2 == 0) {
		het_exist = 2*(1:((length(het_obs))/2))
	} else {
		het_exist = 2*(1:((length(het_obs))/2)) + 1
	}
        #print(n1)
        #print(het_exist)
        #print(het_exist-1)
        #print(het_exist-2)
	for (i in het_exist) {
          #note that the number of samples is length(het_obs)-1
		to_add = log(hw_exact_f(f,i-2,length(het_obs)-1,n1))
		#print(c(i,to_add))
		if (to_add != -Inf) {
			#print(c(as.numeric(het_obs[i-1]),to_add,as.numeric(het_obs[i-1])*to_add))
			ll = ll + as.numeric(het_obs[i-1])*to_add
		}
	}
	###### NOTE THE HACK WITH i-2!!!!! ############################
	return(-ll)
}

#goodness of fit for distribution
GOF_HW = function(F,het_obs,n1) {
  GOF = 0.0
  exp = gen_hw(F,length(het_obs)-1,n1)*sum(het_obs)
  non_zero = which(exp!=0)
  return(sum((het_obs[non_zero]-exp[non_zero])^2/exp[non_zero]))
}



#copy of above function using the beta distributed F
LL_HW_beta = function(pars,het_obs,n1) {
  a = pars[1]
  b = pars[2]
  ll = 0.0
  if (n1 %% 2 == 0) {
    het_exist = 2*(1:((length(het_obs))/2))
  } else {
    het_exist = 2*(1:((length(het_obs))/2)) + 1
  }
  for (i in het_exist) {
    to_add = log(hw_exact_f_beta(a,b,i-2,length(het_obs)-1,n1))
    if (to_add != -Inf) {
      ll = ll + as.numeric(het_obs[i-1])*to_add
    }
  }
  return(-ll)
}

#this function computes the likelihood lumping all individuals together
#het_obs is a matrix with nrow=2*number of indivduals and ncol= number of individuals
#het_obs[i,j] is how many alleles with count i-1 show up in j-1 heterozygotes
LL_combined = function(f,het_obs) {
	num_freq = nrow(het_obs)
	ll = 0.0
	total = sum(het_obs)
	for (i in 1:num_freq) {
          #print(het_obs[i,])
		cur_count = sum(het_obs[i,])
		#NOTE THE i-1 !!!!!
		ll = ll + LL_HW(f,het_obs[i,],i-1)
	}
	return(ll)
}

LL_combined_2 = function(f,het_obs) {
  expected = colSums(gen_hw_mat(f,ncol(het_obs)-1)*rowSums(het_obs))
  expected = expected/sum(expected)
  obs = colSums(het_obs)
  return(-sum(obs*log(expected)))
}

LL_combined_3 = function(f,het_obs) {
  num_freq = nrow(het_obs)
  ll = 0.0
  total = sum(het_obs)
  for (i in 1:num_freq) {
    #print(het_obs[i,])
    cur_count = sum(het_obs[i,])
    #NOTE THE i-1 !!!!!
    ll = ll + cur_count*LL_HW(f,het_obs[i,],i-1)
  }
  return(ll)
}

GOF_combined = function(f,het_obs) {
  num_freq = nrow(het_obs)
  GOF = 0.0
  for (i in 1:num_freq) {
    GOF = GOF + GOF_HW(f,het_obs[i,],i-1)
  }
  return(GOF)
}

GOF_combined_2 = function(f,het_obs) {
  expected = colSums(gen_hw_mat(f,ncol(het_obs)-1)*rowSums(het_obs))
  obs = colSums(het_obs)
  return(sum((expected-obs)^2/expected))
}

#copy of above function using the beta distributed F
LL_combined_beta = function(pars,het_obs) {
  a = pars[1]
  b = pars[2]
  num_freq = nrow(het_obs)
  ll = 0.0
  total = sum(het_obs)
  for (i in 1:num_freq) {
    cur_count = sum(het_obs[i,])
    ll = ll + LL_HW_beta(c(a,b),het_obs[i,],i-1)
  }
  return(ll)
}

#a function to simulate a sample in HW equilibrium by random sampling for just one allele frequency
#n is number of samples, n1 is number of alleles, N is number of sites
sim_random_single = function(n,n1,N) {
  num_alleles = 2*n
  hets = rep(0,n+1)
  for (i in 1:N) {
    config = sample(c(rep(1,n1),rep(0,num_alleles-n1)),num_alleles)
    num_hets = sum(config[2*1:n]==1&config[2*1:n-1]==0) + sum(config[2*1:n]==0&config[2*1:n-1]==1)
    hets[num_hets+1] = hets[num_hets+1]+1
  }
  return(hets)
}

#A function to simulate a sample in HW equilibrium by random sampling
#note that this allows for 0!
sim_random = function(n,freqs) {
  num_alleles = 2*n
  hets = matrix(0,nrow=num_alleles, ncol=n+1)
  for (i in 1:length(freqs)) {
    #print(c(i, freqs[i]))
    for (j in 1:freqs[i]) {
      config = sample(c(rep(1,i),rep(0,num_alleles-i)),num_alleles)
      #determine how many hets
      num_hets = sum(config[2*1:n]==1&config[2*1:n-1]==0) + sum(config[2*1:n]==0&config[2*1:n-1]==1)
      hets[i,num_hets+1] = hets[i,num_hets+1] + 1
    }
  }
  return(hets)
}

#A function to simulate a sample in HW equilibrium by sampling from the multinomial
sim_multi = function(n, freqs,f=0) {
  num_alleles = 2*n
  hets = matrix(0,nrow=num_alleles,ncol=n+1)
  for (i in 1:length(freqs)) {
    theor = gen_hw(f,n,i)
    hets[i,] = rmultinom(1,freqs[i],theor)
  }
  return(hets)
}


#a function to make a theoretical vector for any f and any allele frequency
gen_hw = function(f,n,n1) {
  theor = vector();
  for (i in 1:(n+1)) {
    if (f != 0) {
      theor[i] = hw_exact_f(f,i-1,n,n1)
    } else {
      theor[i] = hw_exact(i-1,n,n1)
    }
  }
  return(theor)
}

#copy of above function for beta distributed f
gen_hw_beta = function(a,b,n,n1) {
  theor = vector()
  for (i in 1:(n+1)) {
    theor[i] = hw_exact_f_beta(a,b,i-1,n,n1)
  }
  return(theor)
}

#generates the entire matrix
gen_hw_mat = function(f,n) {
  theor = matrix(nrow=(2*n+1),ncol=n+1)
  for (i in 1:(2*n+1)) {
    theor[i,] = gen_hw(f,n,i-1)
  }
  return(theor)
}

#copy of above function for beta distributed f
gen_hw_mat_beta = function(a,b,n) {
  theor = matrix(nrow=(2*n+1),ncol=n+1)
  for (i in 1:(2*n+1)) {
    theor[i,] = gen_hw_beta(a,b,n,i-1)
  }
  return(theor)
}

#generates an hw mat where each frequency has a different F
gen_hw_mat_af = function(f_vec,n) {
  theor = matrix(nrow=(2*n+1),ncol=n+1)
  for (i in 1:(2*n+1)) {
    theor[i,] = gen_hw(f_vec[i],n,i-1)
  }
  return(theor)
}

sim_like_ratio = function(obs,n1,n,ntest=10000) {
  num_obs = sum(obs)
  prob_obs = as.vector(obs/num_obs)
  prob_vec = gen_hw(0,n,n1)
  good_ones = which(obs>0)
  LL_ratio = -2*sum(obs[good_ones]*log(prob_vec[good_ones]/(obs[good_ones]/num_obs)))
  test_stats = vector()
  for (i in 1:ntest) {
    cur_draw = rmultinom(1,num_obs,prob=prob_vec)
    good_ones = which(cur_draw>0)
    test_stats[i] = -2*sum(cur_draw[good_ones]*log(prob_vec[good_ones]/(cur_draw[good_ones]/num_obs)))
  }
  return(list(LL_ratio=LL_ratio,stats=test_stats,prob_exact=prob_vec,prob_obs=prob_obs))
}

#this one makes use of each allele frequency separately
sim_like_ratio_combined = function(obs,exp,ntest=10000) {
  n = ncol(obs)
  num_obs = rowSums(obs)
  obs_freq = obs/rowSums(obs)
  exp_freq = exp/rowSums(exp)
  LL_ratio = -2*sum(rowSums(obs*log(exp_freq/obs_freq),na.rm=T),na.rm=T)
  test_stats = vector()
  for (i in 1:ntest) {
    cur_draw = t(apply(cbind(exp_freq,num_obs),1,function(x){rmultinom(1,x[n+1],x[1:n])}))
    cur_draw_freq = cur_draw/rowSums(cur_draw)
    test_stats[i] = -2*sum(rowSums(cur_draw*log(exp_freq/cur_draw_freq),na.rm=T),na.rm=T)
  }
  return(list(LL.ratio = LL_ratio, stats=test_stats))
}

#this one is if you pool all the data first
sim_like_ratio_combined_2 = function(obs,n,ntest=10000) {
  freqs = rowSums(obs)
  th.0 = gen_hw_mat(0,n)
  prob_vec = freqs%*%th.0
  prob_vec = prob_vec/sum(prob_vec)
  obs_combined = colSums(obs)
  num_obs = sum(obs_combined)
  prob_obs = obs_combined/num_obs
  LL_ratio = -2*sum(obs_combined*log(prob_vec/(obs_combined/num_obs)))
  test_stats = vector()
  for (i in 1:ntest) {
    cur_draw = t(rmultinom(1,num_obs,prob=prob_vec))
    test_stats[i] = -2*sum(cur_draw*log(prob_vec/(cur_draw/num_obs)))
  }
  return(list(LL_ratio=LL_ratio,stats=test_stats,prob_exact=prob_vec,prob_obs=prob_obs))
}

site_test = function(site,theor) {
  site = as.numeric(site)
  return(sum(theor[site[4]+1,theor[site[4]+1,]<=theor[site[4]+1,site[3]+1]]))
}

exact_test_site = function(sites,n) {
  theor = matrix(nrow=(2*n+1),ncol=n+1)
  for (i in 1:(2*n+1)) {
    theor[i,] = gen_hw(f=0,n=n,n1=i-1)
  }
  print("made theor")
  p.vals = apply(sites,1,site_test,theor)
  return(p.vals)
}

F_per_ind_old = function(F,total,individual) {
  #note the indexing stupidity. Here is where 1-offset is really stupid.
  #NOTE THAT THIS MEANS YOU NEED TO HAVE THE COUNT OF 0 FREQUENCY ALLELES!!! 
  n = (length(total)-1)/2
  LL = -sum(individual[2:(2*n)]*log((1-F)*(1:(2*n-1)*(2*n-1:(2*n-1)))/(n*(2*n-1)))+(total[2:(2*n)]-individual[2:(2*n)])*log(1-(1-F)*(1:(2*n-1)*(2*n-1:(2*n-1)))/(n*(2*n-1))))
  return(LL)
}

F_per_ind = function(F,total,individual) {
   #note the indexing stupidity. Here is where 1-offset is really stupid.
  #NOTE THAT THIS MEANS YOU NEED TO HAVE THE COUNT OF 0 FREQUENCY ALLELES!!! 
  n = (length(total)-1)/2
  LL = -sum(individual[2:(2*n)]*log(pmax(exp(-300),(1-F)*(1:(2*n-1)*(2*n-1:(2*n-1)))/(n*(2*n-1))))+(total[2:(2*n)]-individual[2:(2*n)])*log(pmax(exp(-300),1-(1-F)*(1:(2*n-1)*(2*n-1:(2*n-1)))/(n*(2*n-1)))))
  return(LL)
}

F_per_ind_detailed = function(F,total,individual) {
  n = (length(total)-1)/2
  print(paste("n =",n))
  term1 = individual[2:(2*n)]
  print(paste("term1 =",term1))
  term2 = (1-F)*(1:(2*n-1)*(2*n-1:(2*n-1)))/(n*(2*n-1))
  print(paste("term2 =",term2))
  term3 = total[2:(2*n)]-individual[2:(2*n)]
  print(paste("term3 =",term3))
  term4 = pmax(10^-80,1-(1-F)*(1:(2*n-1)*(2*n-1:(2*n-1)))/(n*(2*n-1)))
  print(paste("term4 =",term4))
}

F_per_ind_err_combined = function(pars,data) {
  #data is a matrix, data[1,] is total, data[2:nrow(data),] is all the individuals
  npar = length(pars)
  LL = sum(sapply(1:(nrow(data)-1), function(x) {F_per_ind_err(c(pars[x],pars[npar-1],pars[npar]),data[1,],data[x+1,])}))
  return(LL)
}

#the X^2 for the model with error per individual
GOF_err_ind = function(pars,total,individiual) {
  n = (length(total)-1)/2
  exp = total[2:(2*n)]*prob_per_ind_err(pars,n)
  x2 = sum((individual[2:(2*n)]-exp)^2/exp)
  return(x2)
}

F_per_ind_err = function(pars,total,individual) {
  n = (length(total)-1)/2
  p = prob_per_ind_err(pars,n)
  LL = -sum(individual[2:(2*n)]*log(p)+(total[2:(2*n)]-individual[2:(2*n)])*log(1-p))
  return(LL)
}

prob_per_ind = function(F,n) {
  (1-F)*(0:(2*n))*(2*n-0:(2*n))/(n*(2*n-1))
}

prob_per_ind_err_test = function(pars,n,i) {
  #one frequency at a time
  (1-pars[3])*((1-pars[1])*i*(2*n-i)/(n*(2*n-1))) + pars[2]*(1/2*(1-(1-pars[1])*(i+1)*(2*n-(i+1))/(n*(2*n-1)))+1/2*(1-(1-pars[1])*(i-1)*(2*n-(i-1))/(n*(2*n-1))))
}

prob_per_ind_err = function(pars,n) {
  #pars[1] = F, pars[2] = eps_1 (hom called het), pars[3] = eps_2 (het called hom)
  inbetween = (1-pars[3])*((1-pars[1])*1:(2*n-1)*(2*n-1:(2*n-1))/(n*(2*n-1))) + pars[2]*(1/2*(1-(1-pars[1])*(1:(2*n-1)+1)*(2*n-(1:(2*n-1)+1))/(n*(2*n-1)))+1/2*(1-(1-pars[1])*(1:(2*n-1)-1)*(2*n-(1:(2*n-1)-1))/(n*(2*n-1))))
  return(inbetween)
}

#goodness of fit for individual
GOF_HW_ind = function(F,total,individual) {
  n = (length(total)-1)/2
  sum((as.numeric(total)*prob_per_ind(F,n)-individual)^2/(as.numeric(total)*prob_per_ind(F,n)),na.rm=T)
}

#swaps the chromosomes randomly to create pseudoindividuals and estimates stuff
swap_chr = function(F.mat,nrep=10000) {
  num_chr = nrow(F.mat)
  cur_ind = vector()
  means = vector()
  vars = vector()
  for (i in 1:nrep) {
    for (j in 1:num_chr) {
      cur_ind[j] = sample(F.mat[j,],1)
    }
    #change so it's weighted by # of polymorphisms 
    means[i] = mean(cur_ind)
    vars[i] = var(cur_ind)
  }
  return(list(mean=means,var=vars))
}

#this function makes the matrices needed by swap_site_ind
make_ind_matrix = function(dat,ind) {
  ind.all = lapply(dat,function(x){x[ind+1,]})
  num_chr = length(ind.all)
  num_freq = length(ind.all[[1]])
  ind.all = matrix(data=unlist(ind.all),byrow=T,nrow=num_chr,ncol=num_freq)
  return(ind.all)
}

#makes the SFS from the chromosome/individual list
make_sfs = function(dat) {
  sfs.all = lapply(dat,function(x){x[1,]})
  num_chr = length(sfs.all)
  num_freq = length(sfs.all[[1]])
  sfs.all = matrix(data=unlist(sfs.all),byrow=T,nrow=num_chr,ncol=num_freq)
  return(sfs.all)
}

#this function draws fake chromosomes from the genomic backround of an individual
#basically a bootstrap from the genomic background frequencies
swap_site_ind = function(ind,sfs,chr,nreps=10000,numsites=sum(sfs[chr,])) {
  F.est = vector()
  for (i in 1:nreps) {
    sfs.sample = rmultinom(1,numsites,colSums(sfs)/sum(sfs))
    prob = colSums(ind)/colSums(sfs)
    dat = rbinom(length(prob),sfs.sample,prob)
    #print(dat)
    F.est[i] = optim(0,F_per_ind,gr=NULL,sfs.sample,dat,method="L-BFGS-B",lower=-.5,upper=.5)$par
  }
  return(F.est)
}

#A function to take the reduced VCF file and make counts so I can estimate F per individual
#used to jackknife
#note that I uselessly have to append a single site of 0 frequency because I am an idiot
make_ind_counts = function(VCF,ind) {
  #first make total counts
  total = c(1,sapply(1:max(VCF[,4]),function(x){sum(VCF[,4]==x)}),1)
  #now make the individual counts
  #I think that the 1:17 is hard coded for the particular data we analyzed
  ind = c(0,sapply(1:17,function(x){sum(VCF[VCF[,4]==x,ind+4]==1)}),0)
  return(list(total=total,ind=ind))
}

#estimates the F per individual for each chromosome from the reduced VCF file
F_per_chrom = function(VCF) {
  num_ind = ncol(VCF)-4
  res = matrix(ncol=22,nrow=num_ind)
  for (j in 1:22) {
    cur_VCF = VCF[VCF[,1]==j,]
    for (i in 1:num_ind) {
      cur_ind = make_ind_counts(cur_VCF,i)
      res[i,j] = optim(0,F_per_ind,gr=NULL,cur_ind$total,cur_ind$ind,method="L-BFGS-B",lower=-.5,upper=.5)$par
    }
  }
  return(res)
}

jackknife_F = function(VCF,ind,block_size=2000000) {
  cur_start = 0
  cur_end = block_size
  chr_size = max(VCF[,2])
  est = vector()
  i = 1
  while (cur_end < chr_size) {
    if (i%%100 == 0) {
      print(i)
    }
    cur_knife = make_ind_counts(VCF[VCF[,2]<cur_start|VCF[,2]>cur_end,],ind)
    est[i] = optim(0,F_per_ind,gr=NULL,cur_knife$total,cur_knife$ind,method="L-BFGS-B",lower=-.5,upper=.5)$par
    cur_start = cur_end
    cur_end = cur_end + block_size
    i = i + 1
  }
  return(est)
}

#jackknifes but across the whole genome. Thus, VCF should have EVERY chromosome
jackknife_F_genome = function(VCF,ind,block_size=2000000,print_freq=100) {
  est = vector()
  i = 1
  for (j in 1:max(VCF[,1])) {
    #loop over chromosomes
    print(j)
    cur_start = 0
    cur_end = block_size
    chr_size = max(VCF[VCF[,1]==j,2])
    while (cur_end < chr_size) {
      if (i%%print_freq == 0) {
        print(i)
      }
      cur_knife = make_ind_counts(VCF[VCF[,2]<cur_start|VCF[,2]>cur_end|VCF[,1]!=j,],ind)
      est[i] = optim(0,F_per_ind,gr=NULL,cur_knife$total,cur_knife$ind,method="L-BFGS-B",lower=-.5,upper=.5)$par
      cur_start = cur_end
      cur_end = cur_end + block_size
      i = i + 1
    }
  }
  return(est)
}

sites_from_VCF_file = function(file,n=9) {
  dat = read.table(file)
  dat = dat[dat[,4]%in%1:(2*n-1),]
  dat = dat[dat[,1]%in%1:22, ]
  dat.counts = counts_from_VCF(dat)
  return(list(dat=dat,counts=dat.counts))
}

#this function can make the input for LL_combined from the reduced VCF file
counts_from_VCF = function(VCF) {
  num_freq = max(VCF[,4])
  num_ind = max(VCF[,3])
  mat = matrix(0,nrow=num_freq+2,ncol=num_ind+1)
  for (i in 0:(num_freq+1)) {
    for (j in 0:(num_ind)) {
      mat[i+1,j+1] = nrow(VCF[VCF[,4]==i&VCF[,3]==j,])
    }
  }
  return(mat)
}

require(DAAG)

#outputs the good AND bad regions, contrary to the name. The good regions are "outside genes" and the bad regions are "inside genes"
good_regions = function(VCF,gene_list,dist,print=0) {
  #make things to store the output
  good_regions = data.frame()
  bad_regions = data.frame()
  #loop over every site in VCF while keeping track of where you are in the gene_list
  gene_ind = 1
  cur_chr = gene_list[gene_ind,2]
  cur_start = gene_list[gene_ind,4]
  cur_end = gene_list[gene_ind,5]
  i = 1
  while (i <= nrow(VCF)) {
    if (print) {
      if (i %% print == 0) {
        print(i)
      }
    }
    if (VCF[i,1] == cur_chr) {
      if (VCF[i,2] < cur_start - dist) {
        #it's good
        #print("good")
        #print(VCF[i,])
        #print(gene_list[gene_ind,])
        #pause()
        good_regions = rbind(good_regions,VCF[i,])
        i = i+1
      } else if (VCF[i,2] > cur_end + dist) {
        #need to change the gene_index and check if it's good
        #print("Change")
        #print(VCF[i,])
        #print(gene_list[gene_ind,])
        #pause()
        gene_ind = gene_ind + 1
        cur_chr = gene_list[gene_ind,2]
        cur_start = gene_list[gene_ind,4]
        cur_end = gene_list[gene_ind,5]
      } else {
        #it's bad
        #print("bad")
        #print(VCF[i,])
        #print(gene_list[gene_ind,])
        #pause()
        bad_regions = rbind(bad_regions,VCF[i,])
        i = i+1
      }
    } else if (VCF[i,1] == cur_chr + 1) {
      #need to change the gene_index and check if it's good
      #print("Change chr")
      #print(VCF[i,])
      #print(gene_list[gene_ind,])
      #pause()
      gene_ind = min(which(gene_list[,2]==(cur_chr+1)))
      cur_chr = gene_list[gene_ind,2]
      cur_start = gene_list[gene_ind,4]
      cur_end = gene_list[gene_ind,5]
    }
  }
  return(list(good = good_regions,bad=bad_regions))
}

#do this for all inds at once so you don't have to worry about making the good regions so much
F_gene_distance = function(VCF,gene_list,dist=5000) {
  VCF_genes = 5
}



#computes F in a bunch of blocks along a chromosome, sort of like the inverse of the jackknife
F_per_block = function(VCF,ind,block_size=2000000,slide=2000000) {
  cur_start = 0
  cur_end = block_size
  chr_size = max(VCF[,2])
  est = vector()
  snp = vector()
  i = 1
  while (cur_end < chr_size) {
    if (i%%100 == 0) {
      print(i)
    }
    #print(c(i,cur_start,cur_end))
    cur_block_VCF = VCF[VCF[,2]>cur_start&VCF[,2]<cur_end,]
    if (nrow(cur_block_VCF)==0) {
      est[i] = NA
      cur_start = cur_end
      cur_end = cur_end + block_size
      i = i + 1
      next
    }
    nsnp = nrow(cur_block_VCF)
    snp[i] = nsnp
    cur_block = make_ind_counts(cur_block_VCF,ind)
    est[i] = optim(0,F_per_ind,gr=NULL,cur_block$total,cur_block$ind,method="L-BFGS-B",lower=-.6,upper=.6)$par
    cur_start = cur_start + slide
    cur_end = cur_end + slide
    i = i+1
  }
  return(list(est=est,snp=snp))
}

#computes F using the equation toward the bottom of page 241, left column of Keller, Visscher and Goddard 2011, Genetics
canonical_F_per_ind = function(VCF,ind) {
  n = ncol(VCF)-4
  #compute expected homozygosity
  exp_hom = sum(1-VCF[,4]/n*(2*n-VCF[,4])/(2*n-1))
  #compute observed homozygosity
  obs_hom = sum(1-VCF[,4+ind])
  #number of sites
  m = nrow(VCF)
  #compute F
  return((obs_hom-exp_hom)/(m-exp_hom))
}

#this one doesn't correct for small sample size
canonical_F_per_ind_p = function(VCF,ind) {
  n = ncol(VCF)-4
  p = VCF[,4]/(2*n)
  exp_hom = sum(2*p*(1-p))
  obs_hom = sum(1-VCF[,4+ind])
  m = nrow(VCF)
  return((obs_hom-exp_hom)/(m-exp_hom))
}
