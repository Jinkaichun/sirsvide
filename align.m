function immune=align(n,sequnces,escape_score)
sub=sequnces-sequnces(n,:);
immune=sum(escape_score.*(sub~=0),2)/length(sequnces(n,:));
immune(find(immune>1))=1;
end