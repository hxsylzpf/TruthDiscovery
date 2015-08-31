clear;
THRESHOLD = 0.5;

fp = fopen('mturklabels.txt', 'rt');
data = textscan(fp, '%d %s %d');
image_ids = data{1};
labeler_ids = data{2};
labels = data{3};
fclose(fp);
truth = importdata('groundtruth.txt');

% We restrict ourselves to only those images for which we have ground truth. However, it's actually conceivable that the other data might
% help us to learn who the good labelers are, just by their agremeent with each other.
idxs = find(image_ids >= 2000);  
image_ids = image_ids(idxs);
labeler_ids = labeler_ids(idxs);
labels = labels(idxs);

[ imageStats, labelerStats ] = em(image_ids, labeler_ids, labels, 0.5, ones(length(unique(labeler_ids)), 1), ones(length(unique(image_ids)), 1));
[dummy,idxsEM,idxsTruth] = intersect(imageStats{1}, truth(:,1));                                                                                 
sum(((imageStats{2}(idxsEM)) >= THRESHOLD) == truth(idxsTruth,2))/length(idxsTruth)

% Estimate labels by majority vote
count = 1;
for i = 1:size(truth, 1)
	idxs = find(ismember(image_ids, truth(i,1)));
	if isempty(idxs) % Not all images were actually labeled
		continue;
	end
	labelsMajorityVote(count) = sum(labels(idxs))/length(idxs) >= THRESHOLD;
	count = count + 1;
end
sum(labelsMajorityVote' == truth(idxsTruth,2))/length(idxsTruth)
