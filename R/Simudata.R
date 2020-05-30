Simulate<-function(scenario=1,setting=NULL,overlap=NULL,seed=1){
	  # scenario: We have three scenarios as decrisbed in the paper. Scenario 1 has 5 settings, and scenarios 2 and 3 have both two choices: overlap and no overlap
	  # setting: Setiings for secnario 1 only. It's a total of 4 settings
	  # overlap: If yes, then there are some overlap between components, otherwise, there is no overlap. This input only works wen cenario is either 2 or 3
	  # seed: Seed to generate random numbers 
	  #
	  if ((scenario==1) && is.null(setting)){
		      stop("You must specify the setting for scenario 1")
  }
  
  if (((scenario==2) || (scenario==3)) && (is.null(overlap))){
	      stop("You must specify overlap for scenario 2 or 3. It's yes or no for overlap")
    }
#  library(Matrix);
  P=500;
  Q=500;
  r=4; #%4 latent components;
  p=P/5;
  q=Q/5;
  Group1=matrix(0,P, 11);
  Group1[1:10,1]=1;
  Group1[11:20,2]=1;
  Group1[21:30,3]=1;
  Group1[31:40,4]=1;
  Group1[41:50,5]=1;
  Group1[51:60,6]=1;
  Group1[61:70,7]=1;
  Group1[71:80,8]=1;
  Group1[81:90,9]=1;
  Group1[91:100,10]=1;
  Group1[101:P,11]=1;
  Group2=Group1;

  #In general, you can generate N random numbers in the interval (a,b)
  #with the formula r = a + (b-a).*rand(N,1).
  a1=-.5;
  b1=-.3;
  a2=.5;
  b2=.3;
  #latent component
  V1a=V1b=matrix(0,p/2,r);V2a=V2b=matrix(0,q/2,r);
  for (j in 1:r){
	    V1a[,j]=runif(p/2,a1,b1);
    V1b[,j]=runif(p/2,b2,a2);
      V11=rbind(V1a,V1b);
      # for second dataset;
      V2a[,j]=runif(q/2,a1,b1);
        V2b[,j]=runif(q/2,b2,a2);
        V22=rbind(V2a,V2b);
  }

  #multiply effect size for main networks by 2;
  for (j in c(1,11, 21, 31, 41, 51, 61, 71, 81, 91)){
	    V11[j,]=2*V11[j,];
    V22[j,]=2*V22[j,];
  }

  Sx=matrix(1,10,10);
  for (j in 2:10){
	    Sx[1,j]=.7;
    Sx[j,1]=.7;
      for (ii in 3:10){
	          Sx[ii,j]=.7^2;
        Sx[j,ii]=.7^2;
	    Sx[ii,ii]=1;
	  }
  }
  Sigma1=Matrix::bdiag(Sx,Sx,Sx,Sx,Sx,Sx,Sx,Sx,Sx,Sx);
  s=nrow(Sigma1);
  Sigma1=Matrix::bdiag(Sigma1,diag(P-s));
  Sigma2=Sigma1;

  if (scenario==1) {
	    
	    if (setting==1){
		      #                                       %there are 11 groups; groups 1-10 are signal variables, group 11 are noise
		      #                               %variables;
		      #                               %generate a P by 11 indicator matrix where there is 1 in group k if
		      #                               %variable j is in that group;


		    #orthonormalize;
		    #A1=[orth(V11) ;zeros(P-p,4)];
		    #A2=[orth(V22) ;zeros(Q-q,4)];

		    #not orthonomalize
		    A1=rbind(V11 ,matrix(0,P-p,4));
		    A2=rbind(V22 ,matrix(0,Q-q,4));


		    A1[abs(A1)<=10^-8]=0;
		    A2[abs(A2)<=10^-8]=0;

		    TrueVariables1=rep(0,P);
		    TrueVariables1[abs(A1[,1])!=0]=1;

		    TrueVariables2=rep(0,Q);
		    TrueVariables2[abs(A2[,1])!=0]=1;
	    } else if (setting==2) {
		      #%there are 11 groups; first 3 groups corresponding to the first 30
		     # %variables contribute to associatin between X1 and X2, and to the outcome.
		      #%The remaining groups are noise. 
		      #%orthonormalize;
		      #%A1=[orth(V11(1:30,:)) ;zeros(P-30,4)];
		      #%A2=[orth(V22(1:30,:)) ;zeros(Q-30,4)];
		      
		      #%not orthonomalize
		      A1=rbind(V11[1:30,] ,matrix(0,P-30,4));
	      A2=rbind(V22[1:30,] ,matrix(0,Q-30,4));


	        A1[abs(A1)<=10^-8]=0;
	        A2[abs(A2)<=10^-8]=0;
		  
		  TrueVariables1=rep(0,P);
		  TrueVariables1[abs(A1[,1])!=0]=1;
		    
		    TrueVariables2=rep(0,Q);
		    TrueVariables2[abs(A2[,1])!=0]=1;
		      
		      
	    } else if (setting==3){
		    randvec=NULL;
		    #at most 5 random varibles in each group do not contribute to A;
		    for (j in 1:10){
			    randvec=c(randvec,sample((10*(j-1)+1):(10*j),5,1)); 
			    #randvec=[randvec;randi([10*(j-1)+1,10*j],5,1)]; 
		    }
		    V11[randvec,]=0;
		    V22[randvec,]=0;


		    A1=rbind(V11 ,matrix(0,P-nrow(V11),4));
		    A2=rbind(V22 ,matrix(0,Q-nrow(V22),4));


		    A1[abs(A1)<=10^-8]=0;
		    A2[abs(A2)<=10^-8]=0;

		    TrueVariables1=rep(0,P);
		    TrueVariables1[abs(A1[,1])!=0]=1;

		    TrueVariables2=rep(0,Q);
		    TrueVariables2[abs(A2[,1])!=0]=1;


	    } else if (setting==4){
		    randvec=NULL;
		    for (j in 1:3){
			    randvec=c(randvec,sample((10*(j-1)+1):(10*j),5,1));       
		    }

		    A1=rbind(V11[1:30,] ,matrix(0,P-30,4));
		    A2=rbind(V22[1:30,] ,matrix(0,Q-30,4));

		    V11=V11[1:30,] ;
		    V22=V22[1:30,];

		    V11[randvec,]=0;
		    V22[randvec,]=0;


		    A1=rbind(V11 ,matrix(0,P-nrow(V11),4));
		    A2=rbind(V22 ,matrix(0,Q-nrow(V22),4));

		    A1[abs(A1)<=10^-8]=0;
		    A2[abs(A2)<=10^-8]=0;

		    TrueVariables1=rep(0,P);
		    TrueVariables1[abs(A1[,1])!=0]=1;

		    TrueVariables2=rep(0,Q);
		    TrueVariables2[abs(A2[,1])!=0]=1;


	    } else if (setting==5) {
		      #%100 important variables but 25 each loaded on first to fourth factors
		     


		      A1=rbind(V11 ,matrix(0,P-p,4));
	      A2=rbind(V22 ,matrix(0,Q-q,4));
	        
	      whic=matrix(0,P,4);
	      whic[1:25,1]=1;
	      whic[26:50,2]=1;
	      whic[56:75,3]=1; 
	      whic[76:100,4]=1;
	             A1[whic[,1]!=1,1]=0;
	             A1[whic[,2]!=1,2]=0;
		            A1[whic[,3]!=1,3]=0;
		            A1[whic[,4]!=1,4]=0;
			           
			           A2[whic[,1]!=1,1]=0;
			           A2[whic[,2]!=1,2]=0;
				          A2[whic[,3]!=1,3]=0;
				          A2[whic[,4]!=1,4]=0;
					         
					         A1[abs(A1)<=10^-8]=0;
					         A2[abs(A2)<=10^-8]=0;
						        
						        
						        TrueVariables1=rep(0,P);
						        TrueVariables1[apply(A1,1,function(x) sum(abs(x)))!=0]=1;
							       #TrueVariables1(sum(abs(A1),2)~=0)=1;
							       
							       TrueVariables2=rep(0,Q);
							       TrueVariables2[apply(A2,1,function(x) sum(abs(x)))!=0]=1;
	    } 
  } else if (scenario==2){

	    A1=rbind(V11 ,matrix(0,P-p,4));
    A2=rbind(V22 ,matrix(0,Q-q,4)); 
           if (overlap=="yes"){
		           # %overlap
		            whic=matrix(0,P,4);
             whic[1:100,1]=1;
	              whic[1:100,2]=1;
	              whic[1:100,3]=1;
		               whic[1:100,4]=1;
		               A1[whic[,1]!=1,1]=0;
			                A1[whic[,2]!=1,2]=0;
			                A1[,3:4]=0;
					         
					         
					         A2[,1:2]=0;
					         A2[whic[,3]!=1,3]=0;
						          A2[whic[,4]!=1,4]=0;
						          
						          A1[abs(A1)<=10^-8]=0;
							           A2[abs(A2)<=10^-8]=0;
							           TrueVariables1=rep(0,P);
								            TrueVariables1[apply(A1,1,function(x) sum(abs(x)))!=0]=1;
								            TrueVariables2=rep(0,Q);
									             TrueVariables2[apply(A2,1,function(x) sum(abs(x)))!=0]=1;
									           } else if (overlap=="no"){
											            #%No overlap    
											             whic=matrix(0,P,4);
										             whic[1:50,1]=1;
											               whic[51:100,2]=1;
											               whic[1:50,3]=1;
												                 whic[51:100,4]=1;
												              
												                 A1[whic[,1]!=1,1]=0;
														           A1[whic[,2]!=1,2]=0;
														           A1[,3:4]=0;
															             
															             A2[,1:2]=0;
															             A2[whic[,3]!=1,3]=0;
																               A2[whic[,4]!=1,4]=0;
																            
																               A1[abs(A1)<=10^-8]=0;
																	                 A2[abs(A2)<=10^-8]=0;
																	                 TrueVariables1=rep(0,P);
																			           TrueVariables1[apply(A1,1,function(x) sum(abs(x)))!=0]=1;
																			           TrueVariables2=rep(0,Q);
																				             TrueVariables2[apply(A2,1,function(x) sum(abs(x)))!=0]=1;
																				          }
      
  } else if (scenario==3){
	   
	    A1=rbind(V11 ,matrix(0,P-p,4));
    A2=rbind(V22 ,matrix(0,Q-q,4)); 
      if (overlap=="no"){
	          whic=matrix(0,P,4);
        whic[1:25,1]=1;
	    whic[26:50,2]=1;
	    whic[51:100,3]=1;
	        whic[51:100,4]=1;
	        A1[whic[,1]!=1,1]=0;
		    A1[whic[,2]!=1,2]=0;
		    A1[whic[,3]!=1,3]=0;
		        A1[,4]=0;
		        
		        A2[whic[,1]!=1,1]=0;
			    A2[whic[,2]!=1,2]=0;
			    A2[whic[,3]!=1,4]=0;
			        A2[,3]=0;
			        
			        A1[abs(A1)<=10^-8]=0;
				    A2[abs(A2)<=10^-8]=0;
				  } else if (overlap=="yes"){
					      whic=matrix(0,P,4);
				      whic[1:50,1]=1;
				          whic[1:50,2]=1;
				          whic[51:100,3]=1;
					      whic[51:100,4]=1;
					      A1[whic[,1]!=1,1]=0;
					          A1[whic[,2]!=1,2]=0;
					          A1[whic[,3]!=1,3]=0;
						      A1[,4]=0;
						      A2[whic[,1]!=1,1]=0;
						          A2[whic[,2]!=1,2]=0;
						          A2[whic[,3]!=1,4]=0;
							      A2[,3]=0;
							      A1[abs(A1)<=10^-8]=0;
							          A2[abs(A2)<=10^-8]=0;
							          
							        }
      
      TrueVariables1=rep(0,P);
        TrueVariables1[apply(A1,1,function(x) sum(abs(x)))!=0]=1;
        TrueVariables2=rep(0,Q);
	  TrueVariables2[apply(A2,1,function(x) sum(abs(x)))!=0]=1;
  }
  ### Generate data
  set.seed(seed);
 # library(MASS)
  n=200;
  sigma2u=1;
  U=mvrnorm(n,mu=rep(0,r),sigma2u*diag(r))
  A3=c(1, 0, 1, 1)#'; %components 1, 3, and 4 affect response
  E1=MASS::mvrnorm(n, mu=rep(0,P), Sigma1);
  E2=MASS::mvrnorm(n, mu=rep(0,Q), Sigma2);
  E3=MASS::mvrnorm(n, mu=1, sigma2u);
  X1=U%*%t(A1)+E1;X2=U%*%t(A2)+E2;Y=U%*%A3+E3
  list(X1=X1,X2=X2,Y=Y,TrueVar1=TrueVariables1,TrueVar2=TrueVariables2,U=U,LoadA1=A1,LoadA2=A2,Group1=Group1,Group2=Group2)
}


  
