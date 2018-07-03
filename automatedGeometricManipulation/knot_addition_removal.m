function sp_out=knot_addition_removal(sp,data,d,C)
%%-----------------------------------------------------------------------
% Function to Perform knot addition
% ref:
% Zhao&Shen 2001 J AM Stat Assoc
% Sangalli et al. 2009 Appl Statist
% Stein 1981 Ann Statist
% sp=[spx spy spz];
% data=[x y z]; dimension n*3
% d dimension n
% C scalar

% sp=[lsx lsy lsz];
%knot addition
initial_tau=sp(1).knots;
idata=[fnval(sp(1),d) fnval(sp(2),d) fnval(sp(3),d)]; % interpolated data
sse=trace((data-idata)'*(data-idata)); % sum of squared errors

% % sigma=var(sqrt(sum((data-idata)'.^2)));
Sure=sse+C*(sp(1).order+sp(1).number); % Stein unbiased risk estimate, Stein (1981) Appl Statist - see Sangalli et al 2009 Appl Statist
Sure_ka=Sure;
sp_ka=sp;
idata_ka=idata;

% Knots are added in the segment between two knots where the error is max
while Sure>=Sure_ka
    sp=sp_ka;
    Sure=Sure_ka;
    for i=sp(1).order:sp(1).number
        ind=find(d>=sp(1).knots(i) & d<=sp(1).knots(i+1)); % knots are the same for the 3 coordinates
        err(i-sp(1).order+1)=trace((data(ind,:)-idata_ka(ind,:))'*(data(ind,:)-idata_ka(ind,:)));
    end
    tau_ka=sort([sp(1).knots mean([sp(1).knots(sp(1).order+find(err==max(err))-1) sp(1).knots(sp(1).order+find(err==max(err)))])]);
    sp_ka=[spap2(tau_ka,sp(1).order,d,data(:,1)); spap2(tau_ka,sp(2).order,d,data(:,2)); spap2(tau_ka,sp(3).order,d,data(:,3))];
    idata_ka=[fnval(sp_ka(1),d) fnval(sp_ka(2),d) fnval(sp_ka(3),d)];
    sse_ka=trace((data-idata_ka)'*(data-idata_ka));
    Sure_ka=(sse_ka)+C*(sp_ka(1).order+sp_ka(1).number);
end
[~,addedknti]=setdiff(sp(1).knots,initial_tau);

for i=1:length(addedknti)
    for k=-2:2
        tau_kr=sp(1).knots;
        tau_kr(max(sp(1).order+1,min(sp(1).number,addedknti(i)+k)))=[];
        sp_kr=[spap2(tau_kr,sp(1).order,d,data(:,1)); spap2(tau_kr,sp(2).order,d,data(:,2)); spap2(tau_kr,sp(3).order,d,data(:,3))];
        idata_kr=[fnval(sp_kr(1),d) fnval(sp_kr(2),d) fnval(sp_kr(3),d)];
        sse_kr=trace((data-idata_kr)'*(data-idata_kr));
        Sure_kr=(sse_kr)+C*(sp_kr(1).order+sp_kr(1).number);
        if Sure>=Sure_kr
           sp=sp_kr;
           Sure_kr;
           Sure=Sure_kr;
        end
    end
end
sp_out=sp;
