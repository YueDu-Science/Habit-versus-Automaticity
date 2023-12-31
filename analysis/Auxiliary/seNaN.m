function s = seNaN(x)
for i=1:size(x,2)
    igood = find(~isnan(x(:,i,:)));
    if(isempty(igood))
        s(i) = NaN;
    else
        s(i) = std(x(:,i,igood))/sqrt(length(igood));
    end
end