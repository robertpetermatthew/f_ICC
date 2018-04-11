function ICC = f_ICC(data,alpha)
% Computes the Intra Class Correllation Coefficients ICC1, ICC2, ICC3,
% ICC1k, ICC2k, and ICC3k.
% Based on the development by Shrout1979, presentation by McGraw including
% errata corrections
% Data is returned as presented in the form returned by the ICC function in the 
% R package 'DescTools' 

% Shrout, Patrick E. and Fleiss, Joseph L. Intraclass correlations: uses in 
% assessing rater reliability. Psychological Bulletin, 1979, 86, 420-3428.
% McGraw, Kenneth O. and Wong, S. P. (1996), Forming inferences about some 
% intraclass correlation coefficients. Psychological Methods, 1, 30-46. 
% and errata in Psychological Methods, 4, page 390.

% Syntax: 
% Input: 
%   data is an k by m matrix of k targets by m raters
%   alpha is the alpha level for significance using the confidence intevals
%   rho0 is the hypothesised value of ICC (set to zero if unsure)
rho0    = 0;
% Output: 
%   ICC
%   F-value 
%   degrees of freedom df1, df2 
%   p-value corresponding to the F-Value 
%   confidence limits at the significance given in alpha

% Example: (data from McGraw Table 6, 
% data = [103,119;% 82,65;116,106;102,102;99,105;98,100;104,107;
%           62,85;97,101;107,110];
% alpha     = 0.05;
% ICC = f_ICC(data,alpha);

% Verification from R:
% Intraclass correlation coefficients 
%                          type   est F-val df1 df2   p-val lwr.ci upr.ci
% Single_raters_absolute   ICC1 0.722  6.18   9  10 0.00437  0.241  0.922
% Single_random_raters     ICC2 0.720  6.00   9   9 0.00677  0.228  0.922
% Single_fixed_raters      ICC3 0.714  6.00   9   9 0.00677  0.197  0.920
% Average_raters_absolute ICC1k 0.838  6.18   9  10 0.00437  0.389  0.959
% Average_random_raters   ICC2k 0.837  6.00   9   9 0.00677  0.371  0.959
% Average_fixed_raters    ICC3k 0.833  6.00   9   9 0.00677  0.329  0.959

% RPMatthew 20180410



[numTargets,numJudges]=size(data);
meanTarget      = mean(data,2);
meanJudge       = mean(data);
meanTotal       = mean(data(:));

% Within Targets Mean Square (MSW in McGraw1996)
tmp = (data - repmat(meanTarget,1,numJudges)).^2;
WSS = sum(tmp(:));
WMS = WSS / (numTargets*(numJudges - 1));

% Within Judges (Raters) Mean Square (MSC in McGraw1996)
RSS = sum((meanJudge - meanTotal).^2) * numTargets;
CMS = RSS / (numJudges - 1);

% Between Targets Mean Square (MSR in McGraw1996)
BSS = sum((meanTarget - meanTotal).^2) * numJudges;
RMS = BSS / (numTargets - 1);

% Residual Mean Square (MSE in McGraw1996)
ESS = WSS - RSS;
EMS = ESS / ((numJudges - 1) * (numTargets - 1));



ICC_1   = (RMS-WMS)/(RMS+(numJudges-1)*WMS);
ICC_2   = (RMS-EMS)/(RMS+(CMS-EMS)*numJudges/numTargets+(numJudges-1)*EMS);
ICC_3   = (RMS-EMS)/(RMS+(numJudges-1)*EMS);

ICC_1k  = (RMS-WMS)/(RMS);
ICC_2k  = (RMS-EMS)/(RMS+(CMS-EMS)/numTargets);
ICC_3k  = (RMS-EMS)/(RMS);

FObs_1      = RMS/WMS;%anova1(data,[],'off');
FTabled_1L  = finv((1-0.5*alpha),numTargets-1,numTargets*(numJudges-1));
FTabled_1U  = finv((1-0.5*alpha),numTargets*(numJudges-1),numTargets-1);
FL_1        = FObs_1/FTabled_1L;
FU_1        = FObs_1*FTabled_1U;
ICC_1CI     = [(FL_1-1)/(FL_1+(numJudges-1)),(FU_1-1)/(FU_1+(numJudges-1))];
ICC_1kCI    = [1-1/FL_1,1-1/FU_1];

ICC_1F     = RMS/WMS*((1-rho0)/(1+(numJudges-1)*rho0));
ICC_1pVal  = 1-fcdf(ICC_1F,numTargets-1,numTargets*(numJudges-1));
sprintf('%2.5f',ICC_1pVal);
ICC_1df1   = numTargets-1;
ICC_1df2   = numTargets*(numJudges-1);

ICC_1kF     = RMS/WMS*(1-rho0);
ICC_1kpVal  = 1-fcdf(ICC_1kF,numTargets-1,numTargets*(numJudges-1));
sprintf('%2.5f',ICC_1kpVal);
ICC_1kdf1   = numTargets-1;
ICC_1kdf2   = numTargets*(numJudges-1);

FObs_2      = RMS/EMS;
a_2         = (numJudges*ICC_2)/(numTargets*(1-ICC_2));
b_2         = 1+(numJudges*ICC_2*(numTargets-1))/(numTargets*(1-ICC_2));
c_2         = ICC_2/(numTargets*(1-ICC_2));
d_2         = 1+(ICC_2*(numTargets-1))/(numTargets*(1-ICC_2));
v_2         = ((a_2*CMS+b_2*EMS)^2)/((a_2*CMS)^2/(numJudges-1)+((b_2*EMS)^2/((numTargets-1)*(numJudges-1))));
FLstar_2    = finv((1-0.5*alpha),numTargets-1,v_2);
FUstar_2    = finv((1-0.5*alpha),v_2,numTargets-1);
ICC_2CI     = [(numTargets*(RMS-FLstar_2*EMS))/(FLstar_2*(numJudges*CMS+(numJudges*numTargets-numJudges-numTargets)*EMS)+numTargets*RMS),...
                (numTargets*(FUstar_2*RMS-EMS))/(numJudges*CMS+(numJudges*numTargets-numJudges-numTargets)*EMS+numTargets*FUstar_2*RMS)];
           
FObs_2k     = RMS/EMS;
c_2k        = ICC_2k/(numTargets*(1-ICC_2k));
d_2k        = 1+(ICC_2k*(numTargets-1))/(numTargets*(1-ICC_2k));
v_2k        = (c_2k*CMS+d_2k*EMS)^2/((c_2k*CMS)^2/(numJudges-1)+((d_2k*EMS)^2/((numTargets-1)*(numJudges-1))));
FLstar_2k   = finv((1-0.5*alpha),numTargets-1,v_2k);
FUstar_2k   = finv((1-0.5*alpha),v_2k,numTargets-1);
ICC_2kCI    = [(numTargets*(RMS-FLstar_2k*EMS))/(FLstar_2k*(CMS-EMS)+numTargets*RMS),...
                (numTargets*(FUstar_2k*RMS-EMS))/(CMS-EMS+numTargets*FUstar_2k*RMS)];

a   = (numJudges*rho0)/(numTargets*(1-rho0));
b   = 1+(numJudges*rho0*(numTargets-1))/(numTargets*(1-rho0));
c   = (rho0)/(numTargets*(1-rho0));
d   = 1+(rho0*(numTargets-1))/(numTargets*(1-rho0));

ICC_2F     = (RMS)/(a*CMS+b*EMS);
ICC_2pVal  = 1-fcdf(ICC_2F,numTargets-1,(numTargets-1)*(numJudges-1));
sprintf('%2.5f',ICC_2pVal);
ICC_2df1   = numTargets-1;
ICC_2df2   = (numTargets-1)*(numJudges-1);

ICC_2kF     = (RMS)/(c*CMS+d*EMS);
ICC_2kpVal  = 1-fcdf(ICC_2kF,numTargets-1,(numTargets-1)*(numJudges-1));
sprintf('%2.5f',ICC_2kpVal);
ICC_2kdf1   = numTargets-1;
ICC_2kdf2   = (numTargets-1)*(numJudges-1);
            

FObs_3      = RMS/EMS;%anova1(data,[],'off');
FTabled_3L  = finv((1-0.5*alpha),numTargets-1,(numTargets-1)*(numJudges-1));
FTabled_3U  = finv((1-0.5*alpha),(numTargets-1)*(numJudges-1),numTargets-1);
FL_3        = FObs_3/FTabled_3L;
FU_3        = FObs_3*FTabled_3U;
ICC_3CI     = [(FL_3-1)/(FL_3+(numJudges-1)),(FU_3-1)/(FU_3+(numJudges-1))];
ICC_3kCI    = [1-1/FL_3,1-1/FU_3];


ICC_3F     = RMS/EMS*((1-rho0)/(1+(numJudges-1)*rho0));
ICC_3pVal  = 1-fcdf(ICC_3F,numTargets-1,(numTargets-1)*(numJudges-1));
sprintf('%2.5f',ICC_3pVal);
ICC_3df1   = numTargets-1;
ICC_3df2   = (numTargets-1)*(numJudges-1);

ICC_3kF     = RMS/EMS*(1-rho0);
ICC_3kpVal  = 1-fcdf(ICC_3kF,numTargets-1,(numTargets-1)*(numJudges-1));
sprintf('%2.5f',ICC_3kpVal);
ICC_3kdf1   = numTargets-1;
ICC_3kdf2   = (numTargets-1)*(numJudges-1);


clear ICC
ICC{1}.nameMcGraw   = 'Single Raters Absolute';
ICC{1}.nameShrout   = 'ICC(1,1)';
ICC{1}.est          = ICC_1;
ICC{1}.FValue       = ICC_1F;
ICC{1}.df1          = ICC_1df1;
ICC{1}.df2          = ICC_1df2;
ICC{1}.pVal         = ICC_1pVal;
ICC{1}.confInterval = ICC_1CI;

ICC{2}.nameMcGraw   = 'Single Random Raters';
ICC{2}.nameShrout   = 'ICC(2,1)';
ICC{2}.est          = ICC_2;
ICC{2}.FValue       = ICC_2F;
ICC{2}.df1          = ICC_2df1;
ICC{2}.df2          = ICC_2df2;
ICC{2}.pVal         = ICC_2pVal;
ICC{2}.confInterval = ICC_2CI;

ICC{3}.nameMcGraw   = 'Single Fixed Raters';
ICC{3}.nameShrout   = 'ICC(3,1)';
ICC{3}.est          = ICC_3;
ICC{3}.FValue       = ICC_3F;
ICC{3}.df1          = ICC_3df1;
ICC{3}.df2          = ICC_3df2;
ICC{3}.pVal         = ICC_3pVal;
ICC{3}.confInterval = ICC_3CI;

ICC{4}.nameMcGraw   = 'Average Raters Absolute';
ICC{4}.nameShrout   = 'ICC(1,k)';
ICC{4}.est          = ICC_1k;
ICC{4}.FValue       = ICC_1kF;
ICC{4}.df1          = ICC_1kdf1;
ICC{4}.df2          = ICC_1kdf2;
ICC{4}.pVal         = ICC_1kpVal;
ICC{4}.confInterval = ICC_1kCI;

ICC{5}.nameMcGraw   = 'Average Random Raters';
ICC{5}.nameShrout   = 'ICC(2,k)';
ICC{5}.est          = ICC_2k;
ICC{5}.FValue       = ICC_2kF;
ICC{5}.df1          = ICC_2kdf1;
ICC{5}.df2          = ICC_2kdf2;
ICC{5}.pVal         = ICC_2kpVal;
ICC{5}.confInterval = ICC_2kCI;

ICC{6}.nameMcGraw   = 'Average Fixed Raters';
ICC{6}.nameShrout   = 'ICC(3,k)';
ICC{6}.est          = ICC_3k;
ICC{6}.FValue       = ICC_3kF;
ICC{6}.df1          = ICC_3kdf1;
ICC{6}.df2          = ICC_3kdf2;
ICC{6}.pVal         = ICC_3kpVal;
ICC{6}.confInterval = ICC_3kCI;


end

