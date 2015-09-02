import com.mongodb.BasicDBObject;
import java.util.regex.Pattern;
import com.mongodb.util.JSON;
[readColl,db,client] = getdbcoll('LINCS_L1000_LJP2015','data');
writeColl = db.getCollection('CD');

% filter = BasicDBObject();
% filter.append('det_plate',Pattern.compile('CPC'));
batches = readColl.distinct('batch');
batches = j2m(batches);

load('D:\Qiaonan backup\LINCS data\newData\id2gene');
load('D:\Qiaonan backup\LINCS data\newData\2015\LJP59Rid');
geneSymbols = cell(22268,1);
lmIdx = false(22268,1);
for i = 1:numel(rid)
    geneSymbols{i} = dict(rid{i}).gene;
    lmIdx(i) = dict(rid{i}).islm;
end
%%

for i = 53:-1:33
    batch = batches{i};
    fprintf('%s %d\n',batch,i);
    filter = BasicDBObject();
    filter.append('det_plate',Pattern.compile(batch));
    plates = readColl.distinct('det_plate',filter);
    plates = j2m(plates);
    
%     % replcate data on this plate has been removed from db
%     idx = strcmp('LJP005_MCF10A_3H_X2_B17',plates);
%     plates(idx) = [];
    platesRes = cell(1,numel(plates));
    
    tic
    % compute chdir for replicates on each plate
    for j = 1:numel(plates)
        plate = plates{j};
        query = BasicDBObject();
        query.append('det_plate',plate);
        cursor = readColl.find(query);
        arr = cell(cursor.count,1);
        for k = 1:cursor.count
            arr{k} = j2m(cursor.next());
        end
        plateRes = getChdir_2(arr,lmIdx);
        platesRes{j} = plateRes';
    end
    toc
    
    % get unique experiments' sig_ids. Make sure ctl_vehicle is the only
    % denomination of ctrl replicates. Refer to Main-linux.m
    jsonQuery = sprintf('[{$match:{batch:"%s",SM_Pert_Type:{$ne:"ctl_vehicle"}}},{$group:{_id:{"batch":"$batch","pert_id":"$SM_Center_Compound_ID","pert_dose":"$SM_Dose"},replicateCount:{$sum:1}}}]',batch);
    aggregateOutput = readColl.aggregate(JSON.parse(jsonQuery));
    sigIdStructs = j2m(aggregateOutput.results());
    
    % merge chdir replicates and computeSignficance
    chdirArr = mergeReplicates(platesRes,sigIdStructs);
    
    % save chdir arr to db.
    tic
    for j = 1:numel(chdirArr)
        chdirStruct = chdirArr{j};
        writeColl.save(json2dbobj(savejson('',chdirStruct)));
    end
    toc
end
