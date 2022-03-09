function [ metaTop_depth ] = FindMetaTop(dt_dz,thermoD,depths,slope)
%----Author: Jordan S Read 2009 ----
%updated 12/16/2009 with thermoD pass

numDepths = length(depths);
metaTop_depth = mean([depths(2) depths(1)]);
Tdepth = NaN(1,numDepths-1);
for i = 1:numDepths-1
    Tdepth(i) = mean([depths(i+1) depths(i)]);
end
[sortDepth,sortInd] = sort([Tdepth thermoD+1e-6]);
dt_dz = interp1(Tdepth,dt_dz,sortDepth);

% ------------
sumtmp = cumsum(ones(size(dt_dz)));
dt_dz_e = sumtmp.*dt_dz*eps;
for i = 1:numDepths-1 
    if dt_dz(i) == 0
        dt_dz_e(i) = sumtmp(i)*eps;
    end
end
dt_dz_local = dt_dz + dt_dz_e;

[tmp, tmpind] = unique(dt_dz_local);
for i = 1:numDepths-1
    if ~ismember(i, tmpind) 
        dt_dz_local(i) = dt_dz(i) + sumtmp(i)*dt_dz(i)*eps*10;
    end
end
% ------------

thermo_index = 1;
thermoId = numDepths;
for i = 1:numDepths
    if eq(thermoId,sortInd(i))
        thermo_index = i;
        break
    end
end
disp(dt_dz_local)
for i = thermo_index:-1:1 %moving up from thermocline index
    % if dt_dz_local(i) > slope %top of metalimnion
    %     metaTop_depth = sortDepth(i);
    %     break
    % end
    if i-8 < 1 && dt_dz_local(i) > slope    % 连续4m的温度梯度都小于阈值
        metaTop_depth = sortDepth(i);
        break
    end
    if dt_dz_local(i) > slope && dt_dz_local(i-1) > slope ...,
        && dt_dz_local(i-2) > slope && dt_dz_local(i-3) > slope ...,
        && dt_dz_local(i-4) > slope && dt_dz_local(i-5) > slope ...,
        && dt_dz_local(i-6) > slope && dt_dz_local(i-7) > slope ...,        % bottom of metalimnion
        && dt_dz_local(i-8) > slope
        metaTop_depth = sortDepth(i);
        disp(metaTop_depth)
        break
    end
end
%             % if interp can happen
% if thermo_index-i > 1 && dt_dz_local(thermo_index) > slope  
%     metaTop_depth = interp1(dt_dz_local(i:thermo_index),sortDepth(i:thermo_index),slope);
% end


if isnan(metaTop_depth)
    metaTop_depth = min(depths);
end
if metaTop_depth > thermoD
    metaTop_depth = thermoD;
end

end

