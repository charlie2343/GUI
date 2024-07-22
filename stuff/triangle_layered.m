% num = input("How many: ")
% rad = input("rad: ")
% 
% 
% for n = 1:num
% 
%     p = nsidedpoly(3, 'Center', [1 ,1], 'SideLength', rad-n);
% 
%     plot(p);
% end
for n = 0:0.2:0.5
patch([0+n,2-n,1-n],[0+n,4-n,0+n],'r')
end

