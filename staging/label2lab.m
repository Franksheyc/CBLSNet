function y1=label2lab(test_result)
% for i=1:length(test_result)
%     for j=1:max(test_result)
%     switch test_result(i)
%         case j
%             y1(i,j)=1;
% %         case 2
% %             y1(i,2)=1;
% %         case 3
% %             y1(i,3)=1;
% %         case 4
% %             y1(i,4)=1;
% %         case 5
% %             y1(i,5)=1;
%     end
%     end
% end
% end

m=length(test_result);
n=max(test_result);
y1=zeros(m,n);
for i=1:m
    y1(i,test_result(i))=1;
end