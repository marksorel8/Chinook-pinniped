#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  //Data
    DATA_INTEGER(fit_HMM);                    //"flag" for whether to fit the HMM
    DATA_INTEGER(fit_pop_Ast);                //"flag" for whether or not to fit the population- and year-specific Astoria departure timing
    DATA_INTEGER(penalized_complexity);       //"flag" for whether or not to use penalized complexity priors on coefficients
    DATA_INTEGER(Na);                         //total number of possible Astoria departure days
    DATA_INTEGER(start_day);
    DATA_INTEGER(Nb);                         //total number of possible Bonneville arrival days
    DATA_INTEGER(Nyears);                     //total number of years of data
    DATA_IVECTOR(relDOY);                     //release day-of-yeaR for each fish
    DATA_IVECTOR(rel_Year);                   //index of release year for each fish
    DATA_MATRIX(pinCounts);                   //matrix of Z-scored Astoria pinniped counts for survival model
    DATA_ARRAY(phi_x);                        //array of other covariate values (besides sea lions) for survival model
    DATA_ARRAY(X);                            //array of covariate values for travel time model from Astoria to Bonneville
    DATA_INTEGER(dayCov);                     //indivator variable for whether log-release-day is a covariate in the travel time model?
    DATA_VECTOR(dayScales);                   //scaling values (mean and sd) for log release day covariate
    DATA_INTEGER(Npops);                      //number of populations in the data of population-specific Bonneville Detection dates
    DATA_IVECTOR(bonDOY_unk_pop);             //day of bonneville detection for fish tagged in Astoria (unknown population)
    DATA_IVECTOR(relOrigin);                  //veccot of origin (hatchery/natural) of fish released at Astoria, based on presence or abesnce of adipose fin.
    DATA_IVECTOR(bonDOY_known_pop);           //vector of detection DOYs at Bonneville dam for fish tagged as juveniles (known population) Dataset #2
    DATA_IVECTOR(pop);                        //vector of populations of fish (tagged as juveniles) detected at Bonenville dam in Dataset #2
    DATA_IVECTOR(bon_Years);                  //vector of detection Years  at Bonneville dam of fish tagged as juveniles datset #2
    DATA_VECTOR(exp_rate);                    //rate paramater(s) for the exponential prior(s) on the standard deviations of the normal (penalized-complexity) priors (centered at zero) on coefficients. Can have seperate rate paramaters for survival and transition model covariates
  
  //paramaters
    PARAMETER_VECTOR(log_sigma);              //log-standard deviations of penalized-complexity priors on coefficients
    PARAMETER(logit_lambda);                  //logit-scale paramater in Weibull baseline hazard model for transition probability
    PARAMETER(log_alpha);                     //log- shape parameter  in Weibull baseline hazard model for transition probability
    PARAMETER_VECTOR(B);                      //coefficients in hazard model of daily transition probability (of passing Bonneville dam)
    PARAMETER(phi_intercept);                 //intercept in logit-linear model of daily survival probability
    PARAMETER_VECTOR(phi_betas);              //coefficients in logit-linear model of daily survival 
    PARAMETER_VECTOR(A_mu);                   //log-means of population-specific distributions of Astoria departure day
    PARAMETER_VECTOR(A_sigma);                //log_standard deviations for distribution of Astoria departure day
    PARAMETER_VECTOR(A_year_beta);            //years-effect (annual shifts) on log-means of Astoria departure day. Act on all populations the same in a given year.
    //PARAMETER_VECTOR(A_hypers);              //years-effect (annual shifts) on log-means of Astoria departure day. Act on all populations the same in a given year.
    

  // Variable transofrmation
     
    //transform paramaters for hazard (travel-time) model.
    Type lambda = invlogit( logit_lambda);     //scale parameter in hazard model
    Type alpha = exp(log_alpha);               //shape parameter in hazard model
   
    vector<Type> sigma = exp(log_sigma);       //transform the standard deviations of the normal (penalized-complexity) priors on the coefficients
  
  
  
  //Start of penalized-complexity prior on covariate coefficients
  vector<Type> sigma_B = sigma.head(B.size()) ;               //subset of standard deviations of priors on coefficients in transition-probabilty model.
  vector<Type> sigma_phi_beta = sigma.tail(phi_betas.size()); //standard deviations of priors on coefficients in survival model.
  
  Type l3=0;                                                  // initialize likelihood for beta paramater complexity likelihood
  Type l4=0;                                                  // initialize likelihood for penalized complexity prior SD prior
  if(penalized_complexity==1){
       
      l3= -sum(dnorm(B,Type(0.0),sigma_B, true))        
       -sum(dnorm(phi_betas ,Type(0.0),sigma_phi_beta, true));      //normal (centered on zero) priors on coefficients 
     
       l4 = -sum(dexp(sigma_B,exp_rate(0),true))  
       -sum(dexp(sigma_phi_beta ,exp_rate(1),true))  
       -sum(log_sigma);                                        //exponential priors on sigma. Using change in variables formula to calculate probabilty of exp-transformed parameters.
   }
//end of penalized complexity prior     
     
     
    //starting matrix                         //Initial state probability matrix for HMM multistate model [1,0,0] 
    matrix<Type> delta(1,3);
    delta.setZero();                          //initialize with zeros
    delta(0,0)=1.0;                           // and 1 for "alive and downstream on Bonneville" state
  
    //transition probability matrix... first row filled in in likelihood evaluation loop below
    matrix<Type> transMat(3,3);
    transMat.setZero();                      //initialize with zeroes
    transMat(1,1) = Type(1.0);               // if you are passed Bonneville you stay there
    transMat(2,2) = Type(1.0);               // If you are dead you remain dead
  
    //observation matrix
    matrix<Type> obs(3,2);                  // 3*2 matrix first column is not detected second is detected
    obs(0,0)=Type(1.0);                     // if alive and downstream then 100% not detected
    obs(0,1)=Type(0.0);
    obs(1,0)=Type(0.0);
    obs(1,1)=Type(1.0);                     // if passing Bonneville Dam, then 100% detected (~perfect detection probability in fish ladders)
    obs(2,0)=Type(1.0);                     // if dead, definitely not detected
    obs(2,1)=Type(0.0);
  
    //P matrix                              //using Zucchini et al 2016 Book notation the P matrix is a diagonal matrix 
                                            //with the column of "observation-matrix" that corresponds to the observation
                                            //along the diagonal.
    matrix<Type> P(3,3);
    P.setZero();
    P.diagonal()=obs.col(0);               //We will evaluate the foreward probability for all possible capture histories in our
                                           //data set, which we will need eventually for calculating derived survival estimates. 
                                           //So, here I set the P matrix to be what it should be for a 0 (non-observation) in the capture history. 
                                           //I will use this below in calculating the probability of being detected at Bonneville on each
                                           //day for fish departing Astoria on each day, and then at the end the probability of never being detected
                                           //will also be stored for fish departing on each day.
    
    matrix<Type> l_t(1,3);                 //declaration of the matrix of foreward probabilties of a fish being in each state on a given day in the HMM model 

    Type XB;                               //declaration of covariate effects (for all covariates combined) on the log hazard model of transition probability
    int i ;                                //delaration of counter for transition model covariate effects calculations in likelihood
    int j ;                                //delaration of counter for transition model covariate effects calculations in likelihood
    Type  H_t;                             //declaration of piecewise integrate of cumulative Hazard function for a given day.
    Type psi;                              //declaration of daily probability of remaining below Bonneville
  
    vector<matrix<Type> > all_zeros_prob(2); //construct 3-dimensional array (i.e. vector of matrices) to hold all-zero probabilities for hatchery and natural origin (these are the two levels)
    matrix<Type> all_zeros_prob_or(Na,(Nyears+2));//declare matrix to hold the conditional (all-zero) probability of never being detected (died at some point or never passed during the whole study) given that you were released in Astoria on a given day.
                                              //this will be filled in each interation of the loop over (hatchery/natural) origin and then stored in the above 3-d array.
   
    vector<vector<matrix<Type> > > BcondA(2); //begin declaration of 4 dimensional array (i.e. vector of vectors of matrices) to hold [B|A,Y,P,O], where O is hatchery/natural origin
    vector<matrix<Type> > BcondA_Or((Nyears+2));  //begin declartion of 3-dimensional array (i.e. vector of matrices) that will hold the conditional probability of arriving at Bonneville given that you departed Astoria on a given day and year i.e. [B|A,Y] 
    matrix<Type> BcondA_Or_yr(Na,Nb);         //declare matrix that will hold the conditional probability of arriving at Bonneville given that a fish departed Astoria on a given day and year i.e. [B|A,Y] 
    BcondA_Or_yr.setZero();                   //initialize matrix with all zeros because we wont fill every cell (i.e. B < A) but we will sum over those cells and so we want them to be 0.0
    
    vector<matrix<Type> > A_dens_array((Nyears+2));//begin declare 3-dimensional array (i.e. vector of matrices) that will hold the probability of departing Astoria given that you are part of a particular population in a given year [A|P,Y] 
    matrix<Type> A_dens(Na,Npops);             //declare 3matrix that will hold the probability of departing Astoria on each day for each population in each year
    
    Type A_year_effect=0;                  //initialize year-effect on Astoria departure time at 0 for the first (reference year).
    
    vector<matrix<Type> > BP(Nyears+2);      //begin declaration of 3-dimensional array (i.e. vector of 2-d matrices) that will hold the probability of detection at bonneville given that you are in a certain population and in a certain year [B|P,Y]
   
    matrix<Type> BP_yr(Nb,Npops);          //declare matrix that will hold the conditional probability of arriving at Bonneville given that you are from a particular population in a given year i.e. [B|P,Y] 
   
    vector<Type> colSums (Nb);             //declare vector that will hold the scaling (i.e. normalizing) factors for the [B|P,Y], probability of passing bonneville on each date given population membership 
    
    matrix<Type> pop_yr_mort(Npops,Nyears+2);     //matrix of population- and year-specific mortality
   
   Type phi=0;
   
   
 //Begin foreward probability calculation
 
   //We will first calculate and store the likelihood for all possible capture histories,
   //which we will then use to fit the three modules: 1) daily survival (phi),
   // 2)transition probability (psi), and 3) Population-specific Astoria-departure timing.
   
   for (int M=1; M>-1; M--){               //loop through 1) natural- and 2) hatchery-origin.
     int year = Nyears;
     if(M==0) year+=2;
   for ( int K=0; K< year; K++){         //loop through years
     
     if (K>0 && K<Nyears && M==0) A_year_effect =     // if in "natural-origin" loop and after first (reference)year,
       A_year_beta(K-1);                  //update year-effect on Astoria departure time to parameter value. Year 1 is the reference year.
  
  if ( K>=Nyears && M==0) A_year_effect = A_year_beta.sum()/6;   //average year 

   for (int I=0; I<Na; I++){              //loop through possible Astoria departure days
     
     l_t = delta;                         //start HMM multi-state likelihood vector at [1,0,0]
     
     
     for(int J=0; J<(Nb-I); J++){          // loop through days of study for each departure date through end of study
       
       //calculate daily survival probability:
       
          phi = phi_intercept+                                      //intercept
         pinCounts(I+J,K)*phi_betas(0)+                            //sea lion effect
         phi_betas(1)*Type(M);                                     //hatchery (factor) effect
    
         
         i=0;                 
         while (i< (phi_betas.size()-2)){                          //begin loop over other (non-sea lion related) covariates on survival
           phi+=phi_betas(2+i)*phi_x(I+J,K,i);
            i++;
          }
         
         
         phi = invlogit(phi);              //transfor survival from logit to linear (0-1) space.
         

       //Begin calculating covariate effects for hazard model of transiton probability
       XB = Type(1.0);                    //start muliplicative covariate effect at 1.0 (no effect)
       i=0;                               //start counter for covariate number at 0
       if (dayCov==1){                    //if log tag day is a covariate, indicated by dayCov = 1
         XB*=Type(exp((((log(I+start_day)-dayScales(0))/
           dayScales(1))*B(0))));         //multiply the coefficient by the z-scored log calendar tag day; exponentiate; and multiply by XB (which is equal to 1)
         i ++;                            //advance counter
       }
       
       
       j = 0;                             //start counter for non-tag-day covariates (e.g. temperature) 
       while (i < B.size()){              //loop through non-tag-day covariates
         XB *= Type(exp(B(i)*X(I+J,K,j))); //multiply the coefficient by the covariate; exponentiate; and multiply times XB
         i++;                             // advance counter
         j++;                             // advance counter
       }
       
       //calculate piecewise integral of cumulative hazard
       H_t=  (pow(lambda*(J+1),alpha)
                 -pow(lambda*(J),alpha))*
                   XB;                    
       
       psi = Type(1.0)-exp(-H_t);         //calculate probability of passing Bonneville dam given that you haven't yet passed,
                                          //which is equal to 1 minus the survivorship function (-exp(H_t)) of the piecewise integral
       
       
      
       //fill first row of transition matrix
       transMat(0,0) = phi*(Type(1.0)-psi);//probability of surviving and remaining below bonneville
       transMat(0,1) = phi*(psi);         //probability of surviving and passing bonneville
       transMat(0,2) = Type(1.0)-phi;     //probability of dying
       
       
       l_t *= transMat;                   // (matrix) multiply the probabilities that a fish was in each state the 
                                          //previous day by the transition matrix (gamma in HMM notation)
       
       BcondA_Or_yr(I,I+J) = l_t(1)+1e-80;//store the likelihood that a fish arrived at Bonneville on a given day conditional on departing Astoria on a given day [B|A], adding a trivial amoun to prevent numerical issues when taking log(0)        
       
       l_t *=  P;                         //multiply the vector of state probabilities by the P (observation) matrix for when a fish wasn't observed passing Bonneville. This way we move to the next day and estimate the probabilities that it survives and doesn't pass, passes, or dies the next day

     }                                    //end of "J" loop over days of study starting at release/departure day 
     
     all_zeros_prob_or(I,K)= sum(l_t)+1e-80;    //After looping through all the days, the likelihood of the all-zero capture history is the probabiliity that a fish is still in state 1 (alive below) or 3 (dead).
     
     if (M==0 && fit_pop_Ast==1){                           //for natural-origin fish,

     //Loop through the populations and calculate and store the probability that a fish from a given population departed astoria in the current year
     for (int L =0; L<Npops; L++){
         A_dens(I,L)=Type(pbeta(((Type(I)+Type(1.0))/Na),exp(A_mu(L)+A_year_effect),exp(A_sigma(L)))-
           pbeta(((Type(I))/Na),exp(A_mu(L)+A_year_effect),exp(A_sigma(L))))
                                  ;       //subtract the CDF of the lognorm Astoria-departure distribution through the start of the day
                                          //from the CDF through the end of the day to calculate and store the probability that a fish from a given population arrived on a given day in a given year  [A|P,Y]

       } //end of "L" loop over populations
     }
     
     
   } //end of "I" loop over possible Astoria Departure days
    
    BcondA_Or(K)=BcondA_Or_yr;                  //store [B|A,Y] matrix from the current year in the loop in the 3-dimensional array 

   if (M==0 && fit_pop_Ast==1){                //for natural origin fish,
    A_dens_array(K)=A_dens;                   //store matrix of population a year specific Astoria arrival probabilities [A|P,Y]
   // 
     BP_yr=BcondA_Or_yr.transpose() *A_dens; //calculate the probability of arriving at Bonneville on each day in a given year if you are in a particular population [B|P,Y],
   //                                           //which is equal to the [B|A,Y]*[A|P] marginilzed over A, done here with matrix multiplication
     colSums =BP_yr.colwise().sum() ;        //in order to turn [B|P,Y] into a "proper" probability distribution, we will divide each daily value by the sum accross all daily values.
   //                                           //here I calculate the sum accross all daily values for each population
    for (int L =0; L<Npops; L++){          //here I loop through populations and normalize so that [B|P,Y] sums to 1,
        BP_yr.col(L)=BP_yr.col(L)/colSums(L);//making [B|P,Y] a proper probability density function.
     }                                      //end of "L" loop over populations

     BP(K)=BP_yr;                           //store [B|P,Y] for the current year in the loop in the 3-dimensional array for use in likelihod below

     pop_yr_mort.col(K)=
     (  vector<Type> (A_dens.transpose()*
       all_zeros_prob_or.col(K)));           //calculate logit population and year-specific mortality rate

  }
   
       }//end of "K" loop over years
   BcondA(M)=BcondA_Or;                     //store [B|A,O,Y], probability of passing bonneville given departure day, origin, and year, for use in likelihood below
     all_zeros_prob(M)=all_zeros_prob_or;   //store probability of death (all-zero CH) given departure day, year, and origin for likelihood below
}
   
   //Begin likelihood  
  //likelihood for data on fish tagged in Astoria (travel-time (hazard) and survival modules)
  Type l1=0;                              //initialize likelihood for first dataset at 0.0
  vector<int> sim(relDOY.size());
    if(fit_HMM ==1){
   for( int I =0; I<relDOY.size() ; I++){ //loop through individual capture histories
     if(bonDOY_unk_pop(I)==Nb)            //if the fish was never seen again after if was released
       l1 -= log(all_zeros_prob(relOrigin(I)) //what was the log-probability of the all-zero capture history
                   (relDOY(I),rel_Year(I)));//subtract from the negative log-likelihood the the probability f the all-zero capture history, ()the fish either died at some point during the study, or lived and never passed Bonneville) 
     else                                 //otherwise
       l1 -= log(BcondA(relOrigin(I))
                   (rel_Year(I))(relDOY(I),bonDOY_unk_pop(I)))
                                        ; //subtract the log of the probability of the fish passing Bonneville on the day that it did
     

     
   } //end of loop over fish in dataset #1

  }
    
    SIMULATE{//rmultinom
      for( int I =0; I<relDOY.size() ; I++){ 
      int j = relDOY(I);
      Type s = 0;
      Type p_cum =0;
      while ( s==0 && j<Nb ){
        Type p= BcondA(relOrigin(I)) (rel_Year(I))(relDOY(I),j);
        s=rbinom(Type(1.0),p/(1-p_cum));
        p_cum+=p;
        j++;
      }
      if (s==1)sim(I)=j-1;
      else sim(I)=j;
      }
      }
   REPORT(sim); 
   
 
   //likelihood for dataset #2 of fish tagged as juveniles and detected at Bonneville (for population-specific Astoria-departure time module)
   Type l2 = Type(0.0);                   //initialize likelihood for second dataset at 0.0
  if(fit_pop_Ast==1){
   for( int I = 0;
        I<bonDOY_known_pop.size() ; I++){ //loop through population-specific observations (days) of fish passing Bonneville Dam
     l2-=  log(BP(bon_Years(I))(bonDOY_known_pop(I),pop(I))+1e-80)
     ;                                    //subtract the log of [B|P,Y] for each observation from the negative log likelihood component for this dataset

   }//end of loop over fish in dataset #2

  }

   Type Obj_fun = l1+l2+l3+l4;          //add the first through fourth negative log likelihood components to get the full negative log likelihood

    


  
  //REPORTS (These specify what variables are reported back to R)
  
   //ADREPORT means standard deviations are calculated when SDREPORT() is run
    ADREPORT(log_sigma);              //log-standard deviations of penalized-complexity priors on coefficients
    ADREPORT(logit_lambda);                  //logit-scale paramater in Weibull baseline hazard model for transition probability
    ADREPORT(log_alpha);                     //log- shape parameter  in Weibull baseline hazard model for transition probability
    ADREPORT(B);                      //coefficients in hazard model of daily transition probability (of passing Bonneville dam)
    ADREPORT(phi_intercept);                 //intercept in logit-linear model of daily survival probability
    ADREPORT(phi_betas);              //coefficients in logit-linear model of daily survival 
    ADREPORT(A_mu);                   //log-means of population-specific distributions of Astoria departure day
    ADREPORT(A_sigma);                //log_standard deviations for distribution of Astoria departure day
    ADREPORT(A_year_beta);
    
  //REPORT but dont necessariy calculate SDs when SDREPORT() is run
    REPORT(phi);
    REPORT(A_dens);
    REPORT(A_dens_array);
    REPORT(delta);
    REPORT(all_zeros_prob);
    REPORT(BP_yr);
    REPORT(BP);
    REPORT(transMat);
    REPORT(l_t);
    REPORT(l2);
    REPORT(l3);
    REPORT(l1);
    REPORT(l4);
    REPORT(P);
    REPORT(obs);
    REPORT(Obj_fun);
    REPORT(BcondA);
    REPORT(XB);
    REPORT(pop_yr_mort);
      return(Obj_fun);    

    
    }