function [accuracy,spec,sens,PPV,NPV]=getClassificationResult(predict,groundTrues)
% % % % precision=length(intersect(find(predict==1),find(groundTrues==1)))/length(find(predict==1));
% % % % accuracy=length(find(predict==groundTrues))/length(predict);
% % % % recall=length(intersect(find(predict==1),find(groundTrues==1)))/length(find(groundTrues==1));
% % % % specific=length(intersect(find(predict==-1),find(groundTrues==-1)))/length(find(groundTrues==-1));

TP=length(intersect(find(predict==1),find(groundTrues==1)));
TN=length(intersect(find(predict==-1),find(groundTrues==-1)));
FN=length(intersect(find(predict==-1),find(groundTrues==1)));
FP=length(intersect(find(predict==1),find(groundTrues==-1)));


spec=TN/(TN + FP);
sens = TP/(TP + FN);
accuracy= (TP + TN) / (TP + FP + TN + FN);

PPV=TP/(TP+FP); 
NPV = TN/(TN+FN);


% specific=TN/(TN+FP);
% F=(2*precision*recall)/(precision+recall);


% % % % % % % % % result.precision=precision;
% % % % % % % % % result.accuracy=accuracy;
% % % % % % % % % result.recall=recall;
% % % % % % % % % result.specific=specific;
% % % % % % % % % result.F=(2*precision*recall)/(precision+recall);
% % % temp=(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN);
% % % result.MCC=((TP*TN)-(FP*FN))/temp^(0.5);
% % % sum1=0;
% % % sum2=0;
% % % sum3=0;
% % % 
% % % num1=length(predict);
% % % num2=length(find(groundTrues==1));
% % % num3=length(find(groundTrues==-1));
% % % 
% % % for i=1:length(predict)
% % %   if predict(i)==
% % % 
% % % end