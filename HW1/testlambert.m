clear all

close all

case2 = 0
global mu

mu = 132712.0 *1e6;

AU = 1.496e+8;
day = 86400.0;


%#1I/%Oumouamoua
r1 = [3.515868886595499 *1e-2, -3.162046390773074, 4.493983111703389]*AU;
v1 = [-2.317577766980901*1e-3,9.843360903693031*1e-3,-1.541856855538041*1e-2]*AU/day;


%#2I/Borisov
r2 = [7.249472033259724, 14.61063037906177, 14.24274452216359]*AU;
v2 = [-8.241709369476881 *1e-3,-1.156219024581502 *1e-2,-1.317135977481448 * 1e-2]*AU/day;

%#Earth
rE = [-1.796136509111975 *1e-1, 9.667949206859814 *1e-1,-3.668681017942158 *1e-5]*AU;
vE = [-1.720038360888334 *1e-2,-3.211186197806460 *1e-3, 7.927736735960840 *1e-7]*AU/day;

dt = 86400;

%#Assume year is 365 days 2020 is a leap year
year = 365;
secInYear = year*day;
jan = 31;
feb1 = 28;
feb2 = 29;
mar = 31;
apr = 30;
may = 31;
jun = 30;
jul = 31;
aug = 31;
sept = 30;
oct = 31;
nov = 30;
dec = 31;
numberdays1 = year + (aug+sept +oct+nov+dec+jan);
numberdays2 = year + (jun+jul+aug+sept+oct+nov+dec +jan+366);
numberdays3 = year + year + year + (jan + feb2 + mar + apr+ may + jun +jul);

%#arr1 aug 2017 to jan 2019 5 + 12 + 1
arr1 = linspace(0,(secInYear + (aug+sept +oct+nov+dec+jan)*day),numberdays1+1);
arr1day = linspace(0,(year + (aug+sept +oct+nov+dec+jan)),numberdays1+1);

%#arr2 june 2019 - jan 2022 leap year 2020
arr2 = linspace(0,(secInYear + (jun+jul+aug+sept+oct+nov+dec +jan+366)*day),numberdays2+1);
arr2day = linspace(0,(year + (jun+jul+aug+sept+oct+nov+dec +jan+366)),numberdays2+1);

%#dep1 jan 2017 to dec 2017
dep1 = linspace(0,(secInYear),year+1);
dep1day = linspace(0,year,year+1);
XDates = [datetime(2017,1,1:31) datetime(2017,2,1:28) datetime(2017,3,1:31) datetime(2017,4,1:30) datetime(2017,5,1:31) datetime(2017,6,1:30) datetime(2017,7,1:31) datetime(2017,8,1:31) datetime(2017,9,1:30) datetime(2017,10,1:31) datetime(2017,11,1:30) datetime(2017,12,1:31) datetime(2018,1,1)];
Xdates = string(XDates);

%#dep2 jan 2017 - jul 2020
dep2 = linspace(0,(secInYear*3 +(jan + feb2 + mar + apr+ may + jun +jul)*day ),numberdays3+1);
dep2day = linspace(0,(year*3 +(jan + feb2 + mar + apr+ may + jun +jul) ),numberdays3+1);
XDates2 = [datetime(2017,1,1:31) datetime(2017,2,1:28) datetime(2017,3,1:31) datetime(2017,4,1:30) datetime(2017,5,1:31) datetime(2017,6,1:30) datetime(2017,7,1:31) datetime(2017,8,1:31) datetime(2017,9,1:30) datetime(2017,10,1:31) datetime(2017,11,1:30) datetime(2017,12,1:31) ...
    datetime(2018,1,1:31) datetime(2018,2,1:28) datetime(2018,3,1:31) datetime(2018,4,1:30) datetime(2018,5,1:31) datetime(2018,6,1:30) datetime(2018,7,1:31) datetime(2018,8,1:31) datetime(2018,9,1:30) datetime(2018,10,1:31) datetime(2018,11,1:30) datetime(2018,12,1:31)...
    datetime(2019,1,1:31) datetime(2019,2,1:28) datetime(2019,3,1:31) datetime(2019,4,1:30) datetime(2019,5,1:31) datetime(2019,6,1:30) datetime(2019,7,1:31) datetime(2019,8,1:31) datetime(2019,9,1:30) datetime(2019,10,1:31) datetime(2019,11,1:30) datetime(2019,12,1:31)...
    datetime(2020,1,1:31) datetime(2020,2,1:29) datetime(2020,3,1:31) datetime(2020,4,1:30) datetime(2020,5,1:31) datetime(2020,6,1:30) datetime(2020,7,1:31) datetime(2020,8,1)];
Xdates2 = string(XDates2);

%#arrRV1 = uniVari(r1,v1,8000,mu)
 

arrRV1 = zeros(6,length(arr1));
depRVE = zeros(6,length(dep1));

dvi = zeros(length(dep1),length(arr1));
dvr = zeros(length(dep1),length(arr1));

arrRV2 = zeros(6,length(arr2));
depRVE2 = zeros(6,length(dep2));

dvi1 = zeros(length(dep2),length(arr2));
dvr1 = zeros(length(dep2),length(arr2));



output1 = rv2oe(r1,v1,mu)

output2 = rv2oe(r2,v2,mu)



if case2 ==1

for i =1:length(dep2)
    for j =1:length(arr2)

        arrRV2(:,j) = uniVari(r2.',v2.',arr2(j),mu);

        depRVE2(:,i) = uniVari(rE.',vE.',dep2(i),mu);


        dt = arr2(j) - dep2(i);
        if dt<=30*day
            dvi1(i,j) = NaN;
            dvr1(i,j) = NaN ;   
        else
            vel = lambertCurtis(depRVE2(1:3,i),arrRV2(1:3,j),dt,mu,1);
            
            
            dvitest = norm(vel(:,1)-depRVE2(4:6,i));
            

            dvrtest = norm(vel(:,1)-depRVE2(4:6,i)) + norm(arrRV2(4:6,j)-vel(:,2));

            %#intercept case
            if dvitest <= 20
                dvi1(i,j) = dvitest;
                %#print(dvi1[i,j])
            else
                dvi1(i,j) = NaN;
            end
            
            %#rendezvous case
            if dvrtest <= 60
                dvr1(i,j) =  dvrtest;
                
            else
                dvr1(i,j) =  NaN;
            end
        end

        
        
        %print(arrRV2(1:3,j))
    end
    
end


figure(1)
contourf(dep2day,arr2day,dvi1.',10)
colorbar
grid on 
grid minor 
xticklabels([Xdates2(1) Xdates2(201) Xdates2(401) Xdates2(601) Xdates2(801) Xdates2(1001) Xdates2(1201)])
xlabel('Departure Dates')
ylabel('Time of Flight (days)')
ylim([900 976])
xlim([0 700])
title('Borisov \Delta{v} for Flyby Case')

figure(2)
contourf(dep2day,arr2day,dvr1.',10)
colorbar
grid on 
grid minor 
xticklabels([Xdates2(1) Xdates2(201) Xdates2(401) Xdates2(601) Xdates2(801) Xdates2(1001) Xdates2(1201)])
xlabel('Departure Dates')
ylabel('Time of Flight (days)')
title('Borisov \Delta{v} for Rendezvous Case')

ylim([850 976])
xlim([0 750])





else

for i =1:length(dep1)
    for j =1:length(arr1)

        %arrRV2[:,j] = uniVari(r2,v2,arr2[j],mu)
        
        RV = uniVari(r1.',v1.',arr1(j),mu);
%         R1 = RV(1:3);
        V2 = RV(4:6);
        arrRV1(:,j) = [RV];
        %depRVE2[:,i] = uniVari(rE,vE,dep2[i],mu)

        
        RV = uniVari(rE.',vE.',dep1(i),mu);
        depRVE(:,i) = RV;
        dt = arr1(j) - dep1(i);
        if dt<=10*day
            dvi(i,j) = NaN;
            dvr(i,j) = NaN ;   
        else
%             vel = lambertCurtis(depRVE2[1:3,i],arrRV2[1:3,j],dt,mu,1)
            
            output =lambertCurtis(depRVE(1:3,i),arrRV1(1:3,j),dt,mu,1);
            V1 = output(:,1);
            V2 = output(:,2);
            
            dvitest = norm(V1-depRVE(4:6,i));
            

            dvrtest = norm(V1-depRVE(4:6,i)) + norm(arrRV1(4:6,j)-V2);

            %#intercept case
            if dvitest < 20
                dvi(i,j) = dvitest;
                %#print(dvi1[i,j])
            else
                dvi(i,j) = NaN;
            end
            
            %#rendezvous case
            if dvrtest <= 50
                dvr(i,j) =  dvrtest;
                
            else
                dvr(i,j) =  NaN;
            end
        end

        
        
        %print(arrRV2(1:3,j))
    end
    
end

figure(1)
contourf(dep1day,arr1day,dvi.',10)
colorbar
grid on 
grid minor 
xticklabels([Xdates(1) Xdates(51) Xdates(101) Xdates(151) Xdates(201) Xdates(251) Xdates(301) Xdates(351)])
xlabel('Departure Dates')
ylabel('Time of Flight (days)')
title('Oumouamoua \Delta{v} for Flyby Case')

figure(2)
contourf(dep1day,arr1day,dvr.',10)
colorbar
grid on 
grid minor 
xticklabels([Xdates(1) Xdates(51) Xdates(101) Xdates(151) Xdates(201) Xdates(251) Xdates(301) Xdates(351)])
xlabel('Departure Dates')
ylabel('Time of Flight (days)')
title('Oumouamoua \Delta{v} for Rendezvous Case')

end


























