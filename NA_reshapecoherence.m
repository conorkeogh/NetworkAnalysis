function coherencematrix = RTT_reshapecoherence(data)

    testgrid = zeros(8,8);
    %testlabels = 1:28;
    ind = 1;

    for i = 1:7
        for j = i:7
            testgrid(i,j+1) = data(ind);
            testgrid(j+1,i) = data(ind);
            if i == j
                testgrid(i,j) = 1;
            end
            ind = ind+1;
        end
    end
    testgrid(8,8) = 1;
    
    coherencematrix = testgrid;
    
end