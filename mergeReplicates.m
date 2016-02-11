function [ chdirs ] = mergeReplicates( platesRes,sigIdStructs)
%Merge chdir replicates and computes the significance.
%   Detailed explanation goes here
disp('start merge replicates');

chdirStructsAllPlates = [platesRes{:}];

% chdirLm = cellfun(@(x)x.chdirLm, chdirStructsAllPlates,'UniformOutput',false);
chdirLm = cellfun(@(x)x.chdirLm, chdirStructsAllPlates,'UniformOutput',false);
chdirLm = [chdirLm{:}];

chdirFull = cellfun(@(x)x.chdir, chdirStructsAllPlates,'UniformOutput',false);
chdirFull = [chdirFull{:}];

repCount = cellfun(@(x)x.replicateCount,sigIdStructs,'UniformOutput',false);
repCount = [repCount{:}];

% caculate cosine distance null distribution.
% about 90 seconds execution time for 360 experiments
cosDistLm = permuByRepCount(chdirLm,repCount);
cosDistFull = permuByRepCountFull(chdirFull,repCount);

sigLevel = 0.01;

rpKeys = fieldnames(cosDistFull);
for i = 1:numel(rpKeys)
    eval(sprintf('rankNull = sort(cosDistFull.%s);',rpKeys{i}));
    sampleCount = numel(rankNull);
    cutoff = rankNull(sampleCount*sigLevel);
    eval(sprintf('cosDistFull.%scutoff = cutoff;',rpKeys{i}));
end


chdirs = cell(numel(sigIdStructs),1);
% meanDists = NaN(numel(sigIdStructs),1);

for i = 1:numel(sigIdStructs)
    sigIdStruct = sigIdStructs{i};
    chdir = sigIdStruct.x0x5F_id;
    chdir.replicateCount = sigIdStruct.replicateCount;
    if sigIdStruct.replicateCount==1
        chdirReplicate = dict(chdirStructsAllPlates,@(x)strcmp(x.SM_Center_Compound_ID,sigIdStruct.x0x5F_id.pert_id)&&(x.SM_Dose==sigIdStruct.x0x5F_id.pert_dose));
        chdir.chdirLm = chdirReplicate.chdirLm';
        chdir.chdirFull = chdirReplicate.chdir';
        chdir.geo_id = {chdirReplicate.cid};
        chdirs{i} = addFields(chdir,chdirReplicate);
    else
        chdirReplicates = dict(chdirStructsAllPlates,@(x)strcmp(x.SM_Center_Compound_ID,sigIdStruct.x0x5F_id.pert_id)&&(x.SM_Dose==sigIdStruct.x0x5F_id.pert_dose));
        assert(numel(chdirReplicates)==sigIdStruct.replicateCount);
        chdir.geo_id = cellfun(@(x)x.cid,chdirReplicates,'UniformOutput',false);
        chdirVectorsLm = cellfun(@(x)x.chdirLm,chdirReplicates,'UniformOutput',false);
        chdirVectorsLm = [chdirVectorsLm{:}];
        meanChdirVectorLm = mean(chdirVectorsLm,2);
        chdir.chdirLm = (meanChdirVectorLm/norm(meanChdirVectorLm))';
        chdirMeanDistLm = mean(pdist(chdirVectorsLm','cosine'));
        eval(sprintf('cosDistByRepCount = cosDistLm.repCount%d;',chdir.replicateCount));
        chdir.pvalue = distribution2pval(cosDistByRepCount,chdirMeanDistLm);
        chdir.chdirMeanDistLm = chdirMeanDistLm;
%         meanDists(i) = chdirMeanDistLm;
        
         % meanChdirFull criterion
        chdirVectorsFull = cellfun(@(x)x.chdir,chdirReplicates,'UniformOutput',false);
        chdirVectorsFull = [chdirVectorsFull{:}];
        meanChdirVectorFull = mean(chdirVectorsFull,2);
        chdir.chdirFull = (meanChdirVectorFull/norm(meanChdirVectorFull))';
        [~,meanSortIdx] = sort(meanChdirVectorFull.^2,'descend');
        
        % rp criterion
        [~,sortIdx] = sort(chdirVectorsFull.^2,1,'descend');
        [~,rank] = sort(sortIdx,1);
        rp = prod(rank,2).^(1/size(rank,2));
        [~,rpSortIdx] = sort(rp);
        eval(sprintf('diffIdx = rp<=cosDistFull.rp%dcutoff;',sigIdStruct.replicateCount));
        diffCount = sum(diffIdx);
        
        % equal sign criterion
        allPosIdx = sum(chdirVectorsFull>=0,2)== sigIdStruct.replicateCount;
        allNegIdx = sum(chdirVectorsFull<=0,2) == sigIdStruct.replicateCount;
        sameSignIdx = allPosIdx | allNegIdx;
        sameSignIdx = find(sameSignIdx);
        
        % intersection that meets all three criterions.
        intersection = intersect(meanSortIdx(1:diffCount),rpSortIdx(1:diffCount));
        intersection = intersect(intersection,sameSignIdx);
        
        % now intersection is sorted by meanSortIdx order.
        intersectionMemberIdx =  ismember(meanSortIdx(1:diffCount), intersection);
        chdir.sigIdx = meanSortIdx(intersectionMemberIdx)';
        
        chdirs{i} = addFields(chdir,chdirReplicates{1});
    end
end


% adjustIdx = (meanDists>=1);
% adjustCount = mean(sigCount(adjustIdx));
% 
% for i = 1:numel(sigIdStructs)    
%     currentSigCount = floor(sigCount(i)-adjustCount);
%     if currentSigCount<=0 
%         chdirs{i} = rmfield(chdirs{i},'sigIdx');
%     elseif ~isnan(currentSigCount)
%         chdirs{i}.sigIdx=chdirs{i}.sigIdx(1:currentSigCount)';
%     end
% end

disp('replicates Merged');
end




function [elems] = dict(arr,matchFunc)

   matchIdx = cellfun(matchFunc,arr,'UniformOutput',false);
   matchIdx = [matchIdx{:}];
   elems = arr(matchIdx);
   if numel(elems)==1
       elems = elems{1};
   end

end