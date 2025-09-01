function logError(ME)
    logFile = 'matlab_errors.log';
    
    % 如果文件不存在则创建，存在则追加
    fid = fopen(logFile, 'a+');
    
    fprintf(fid, '[%s] 错误发生\n', datestr(now));
    fprintf(fid, '标识符: %s\n', ME.identifier);
    fprintf(fid, '消息: %s\n', ME.message);
    fprintf(fid, '堆栈跟踪:\n');
    
    for k = 1:length(ME.stack)
        fprintf(fid, '\t文件: %s\n\t函数: %s\n\t行号: %d\n\n',...
                ME.stack(k).file,...
                ME.stack(k).name,...
                ME.stack(k).line);
    end
    fprintf(fid, '----------------------------------------\n');
    fclose(fid);
end