data {
   int<lower=1> nSubjects;
   int<lower=1> nConditions;
   int<lower=1> nTotalTrials;
   int subject[nTotalTrials];
   int condition[nSubjects];
   int<lower=0,upper=1> correct[nTotalTrials];
   vector[nTotalTrials] trialOrder;
   real<lower=0,upper=1> chanceCorrect;
}
transformed data{
   vector[nTotalTrials] logTrialOrder;
   logTrialOrder = log(trialOrder);

}
parameters {
   vector<lower=0>[nConditions] midpoint_group;
   vector<lower=0>[nConditions] slope_group;

   vector<lower=0>[nSubjects] midpoint_subject;
   vector<lower=0>[nSubjects] slope_subject;

   //real<lower=0> kappaMinusTwo;
   //real<lower=0> omega;

   //vector<lower=0, upper=.1>[nConditions] initialLevel;
   //vector<lower=0, upper=.01>[nSubjects] initial_level;

}
model {
   real trialEffect;
   real pCorrect;
   real maxLevel;
   real initialLevel_tr;
   real groupMean;
   real groupSD;
   real pIER;
   real slope_tr;
   real midpoint_tr;

   maxLevel = 1;
   groupMean = 1;
   groupSD = .5;


   midpoint_group ~ normal(groupMean, groupSD);
   slope_group ~ normal(groupMean, groupSD);

   for (sj in 1:nSubjects) { 
      midpoint_subject[sj] ~ normal(midpoint_group[condition[sj]], groupSD);
      slope_subject[sj] ~ normal(slope_group[condition[sj]], groupSD);
   }


   for (tr in 1:nTotalTrials) { 

      //midpoint_tr = midpoint_group[condition[subject[tr]]];
      //slope_tr = slope_group[condition[subject[tr]]];

      midpoint_tr = midpoint_subject[subject[tr]];
      slope_tr = slope_subject[subject[tr]];

      initialLevel_tr = 0; //initial_level[subject[tr]];

      trialEffect = 1-exp(-exp(2 * slope_tr * midpoint_tr / log2()
                               * (logTrialOrder[tr] - log(midpoint_tr)) 
                               + log(log2())));
      pIER = initialLevel_tr + (maxLevel-initialLevel_tr) * trialEffect; 
      //pCorrect = (1-trialEffect) + trialEffect * chanceCorrect;
      pCorrect = (1 - pIER) + pIER * chanceCorrect;

      correct[tr] ~ bernoulli(pCorrect);
   }
   
}
