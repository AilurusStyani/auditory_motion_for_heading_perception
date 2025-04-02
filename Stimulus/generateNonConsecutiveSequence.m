function randomSequence = generateNonConsecutiveSequence(N,min,max,target,times)
if nargin<=3
    target=NaN;
    times = 0;
end

if N <= 0
    error('N 必须为正整数');
end

if min >= max
    error('Invalid number of min and max number.');
elseif target < min || target > max
    error('Invalid number of target.');
elseif N < 2*times
    error('Invalid number of times');
end

randomSequence = nan(1, N); % 初始化序列
randomSequence(1) = randi([min,max]); % 生成第一个随机数字

while randomSequence(1) == target
    randomSequence(1) = randi([min,max]);
end
 
numberCandidate = ones(1,times)*target;
seq = mod(1:N-times+ceil(N/10),10);
del = seq ==target;
seq(del) = [];
numberCandidate = cat(2,numberCandidate,seq(1:N-times-1));
seqlength = length(numberCandidate)+1;
for i = 2:seqlength
    while true
        pushTarget = sum(numberCandidate == target) >= sum(numberCandidate ~= target);
        if pushTarget && randomSequence(i-1)~=target
            ranIndex = 1;
        else
            ranIndex = randi(length(numberCandidate));
        end
        nextNum = numberCandidate(ranIndex);
        if nextNum ~= randomSequence(i-1)
            randomSequence(i) = nextNum;
            numberCandidate(ranIndex) = [];
            break
        end
    end
end
end