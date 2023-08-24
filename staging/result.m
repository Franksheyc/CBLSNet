function y=result(x)
y=[];
% for i=1:size(x,1)
% [~,y(i)]=max(x(i,:));
% end
[~,y]=max(x');
if ~isempty(y)
y=y';
end
