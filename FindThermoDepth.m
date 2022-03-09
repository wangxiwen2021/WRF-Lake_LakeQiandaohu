function [ thermoD,thermoInd,dt_dz,SthermoD,SthermoInd ] = ...
    FindThermoDepth( tVar,depths,Smin )
%----Author: Jordan S Read 2009 ----
% updated 3 march 2011

% removed signal processing toolbox function 'findpeaks.m' replaced with
% 'locPeaks.m' which was written to provide the same functionality

seasonal = false;
if nargout > 3
    seasonal = 1;
end
if nargin < 3
    Smin = -0.2;
end
dtPerc = 0.15; %min percentage max for unique thermocline step
numDepths = length(depths);
dt_dz = NaN(1,numDepths-1);

for i = 1:numDepths-1
    dt_dz(i) = (tVar(i+1)-tVar(i))/(depths(i+1)-depths(i));
end

if seasonal
    %look for two distinct maximum slopes, lower one assumed to be seasonal
    [mDtZ,thermoInd] = min(dt_dz);          %find max slope
    % disp(thermoInd)
        thermoD = mean([depths(thermoInd)...
        depths(thermoInd+1)]);                  %depth of max slope
    if thermoInd > 1 && thermoInd < numDepths-1 %if within range, 
        Sdn = -(depths(thermoInd+1)-depths(thermoInd))/...
            (dt_dz(thermoInd+1)-dt_dz(thermoInd));
        Sup = (depths(thermoInd)-depths(thermoInd-1))/...
            (dt_dz(thermoInd)-dt_dz(thermoInd-1));
        upD  = depths(thermoInd);
        dnD  = depths(thermoInd+1);
        if ~any([isinf(Sup) isinf(Sdn)])
            thermoD = dnD*(Sdn/(Sdn+Sup))+upD*(Sup/(Sdn+Sup));
        end
    end
    dtCut = max([dtPerc*mDtZ Smin]);
    [pks,locs] = locPeaks(dt_dz,dtCut);
    if isempty(pks)
        SthermoD = thermoD;
        SthermoInd = thermoInd;
    else
        mDtZ = pks(length(pks));
        SthermoInd = locs(length(pks));
        if SthermoInd > thermoInd+1
            SthermoD = mean([depths(SthermoInd)...
                depths(SthermoInd+1)]);
            if SthermoInd > 1 && SthermoInd < numDepths-1
                Sdn = -(depths(SthermoInd+1)-depths(SthermoInd))/...
                    (dt_dz(SthermoInd+1)-dt_dz(SthermoInd));
                Sup = (depths(SthermoInd)-depths(SthermoInd-1))/...
                    (dt_dz(SthermoInd)-dt_dz(SthermoInd-1));
                upD  = depths(SthermoInd);
                dnD  = depths(SthermoInd+1);
                if ~any([isinf(Sup) isinf(Sdn)])
                    SthermoD = dnD*(Sdn/(Sdn+Sup))+upD*(Sup/(Sdn+Sup));
                end
            end
        else
            SthermoD = thermoD;
            SthermoInd = thermoInd;
        end
    end
    if SthermoD < thermoD;
        SthermoD = thermoD;
        SthermoInd = thermoInd;
    end
    
else
    [mDtZ,thermoInd] = min(dt_dz);          %find max slope
    % disp(mDtZ)
    % disp(thermoInd)
        thermoD = mean([depths(thermoInd)...
        depths(thermoInd+1)]);                  %depth of max slope
    if thermoInd > 1 && thermoInd < numDepths-1 %if within range, 
        Sdn = -(depths(thermoInd+1)-depths(thermoInd))/...
            (dt_dz(thermoInd+1)-dt_dz(thermoInd));
        Sup = (depths(thermoInd)-depths(thermoInd-1))/...
            (dt_dz(thermoInd)-dt_dz(thermoInd-1));
        upD  = depths(thermoInd);
        dnD  = depths(thermoInd+1);
        if ~any([isinf(Sup) isinf(Sdn)])
            thermoD = dnD*(Sdn/(Sdn+Sup))+upD*(Sup/(Sdn+Sup));
        end
    end
end

end
