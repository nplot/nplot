function list = subset(datalist,mat,values)


index = true(size(datalist.coordlist,1),1);
for i=1:size(datalist.coordlist,1)
    index(i) = all(abs((mat * datalist.coordlist(i,:)') - values(:))<1E-6);
    %index = index & (datalist.coordlist(:,dims(i)) == values(i));
end

list= sublist(datalist,index);


