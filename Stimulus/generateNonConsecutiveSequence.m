function random_sequence = generateNonConsecutiveSequence(N,min,max)
    if N <= 0
        error('N 必须为正整数');
    end
    
    random_sequence = zeros(1, N); % 初始化序列
    random_sequence(1) = randi([min,max]); % 生成第一个随机数字
    
    for i = 2:N
        while true
            next_num = randi([min,max]); % 生成下一个随机数字
            if next_num ~= random_sequence(i-1) % 确保与前一个数字不同
                random_sequence(i) = next_num;
                break;
            end
        end
    end
end